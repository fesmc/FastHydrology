module fast_hydrology_k24
    ! K24 effective-pressure / water-flux model (Kazmierczak et al 2024,
    ! https://doi.org/10.5194/tc-18-5887-2024).
    !
    ! Diagnostic only: reads ice geometry, bed, melt, sliding speed and Glen
    ! rate factor; returns the distributed water flux magnitude q (and its
    ! components q_x, q_y), the effective pressure N, the water pressure
    ! p_w = Po - N, and the reportable water-layer thickness W. Does not touch
    ! the till storage W_til, which is owned by the bucket model.
    !
    ! ============================================================
    ! Reference implementation
    ! ============================================================
    ! This module is a line-by-line port of the `kazmierczak2024` model in
    ! TakisAngelides/FastHydrology.jl (src/models/kazmierczak2024/{model,
    ! water_flux,effective_pressure,sliding_law}.jl). Where the two could
    ! differ, the Julia side is the source of truth. The call sequence mirrors
    ! `update_steady_state!` in run.jl:
    !
    !     update_q  ->  update_N  ->  update_W
    !
    ! Parameter defaults below reproduce the `KazmierczakHydroModel`
    ! constructor's keyword defaults exactly, including the four KORI-ULB
    ! clamps it deliberately leaves off (W_min = 0, W_max = Inf,
    ! min_pressure_fraction = 0, q_max = Inf). Pass KORI-ULB's own values
    ! (1e-8 / 0.015 / 0.02 / 1e5 m2/yr) through the namelist if you want them.
    !
    ! ============================================================
    ! Units convention (differs from the Julia reference on purpose)
    ! ============================================================
    ! FastHydrology.jl carries the melt rate as a MASS rate [kg/m2/s] and
    ! divides by rho_w when seeding psi_out. This library carries it as a
    ! water-equivalent VOLUME rate `mdot` [m/s], matching bucket.f90 and the
    ! rest of the Fortran API, so the seed is simply mdot*dx*dy and every
    ! melt-like source term picked up from the Julia side acquires an extra
    ! 1/rho_w:
    !
    !     Julia   mdot_total += tau_b*v_b/L_w            [kg/m2/s]
    !     Fortran mdot_total += tau_b*v_b/(L_w*rho_w)    [m/s]
    !
    !     Julia   mdot_total += |q*grad(phi0)|/L_w       [kg/m2/s]
    !     Fortran mdot_total += |q*grad(phi0)|/(L_w*rho_w)
    !
    ! Everything else (q [m2/s], Q [m3/s], N/Po/phi0 [Pa], W [m]) is SI and
    ! identical to the Julia side.
    !
    ! ============================================================
    ! Deliberate deviations from FastHydrology.jl
    ! ============================================================
    ! Each is marked with a "DEVIATION:" comment at its site. In summary:
    !   1. Degenerate 0/0 and x/0 forms that Julia leaves to IEEE arithmetic
    !      (S_inf with Q == 0, N_inf with S_inf == 0, N with N_inf == 0) are
    !      resolved by an explicit branch here, taking the same limit Julia's
    !      own `overwrite_where!` cleanups take. This library is built with
    !      -Ofast (-fp-model fast=2), under which Inf/NaN propagation is not
    !      reliable, so the branch is both safer and closer to intent.
    !   2. `masked_mean` over an empty mask returns 0 here rather than NaN.
    !   3. psi_out is zeroed over the whole grid at the start of every sweep.
    !      Julia's persists across calls, so its non-grounded cells carry stale
    !      values; nothing downstream reads them on either side.

    use nml

    implicit none

    integer, parameter :: dp = kind(1.d0)

    real(dp), parameter :: K24_PI = 3.14159265358979323846_dp

    ! Seconds per year used by the parameter defaults FastHydrology.jl writes
    ! as perYear2perSecond(...) -- eta_w and the sliding laws' u0. Julia's
    ! SECONDS_PER_YEAR is 60^2*24*365.25 (Julian year). Deliberately NOT the
    ! same constant as fast_hydrology::SEC_PER_YEAR (3.1556926e7, a tropical
    ! year), which converts the public API's time argument: that one is
    ! calendar bookkeeping, this one reproduces a Julia default. They differ
    ! by 0.002%.
    real(dp), parameter :: K24_SEC_PER_YEAR = 60.0_dp * 60.0_dp * 24.0_dp * 365.25_dp

    ! ---------- Substrate-type enum (par%substrate_type) ----------
    integer, parameter, public :: K24_SUBSTRATE_HARD  = 0
    integer, parameter, public :: K24_SUBSTRATE_SOFT  = 1
    integer, parameter, public :: K24_SUBSTRATE_MIXED = 2

    ! ---------- Flow-routing enum (par%flux_solver) ----------
    ! Mirrors AbstractPsiOutAlgorithm in FastHydrology.jl/model.jl.
    ! NOTE: the numbering changed in this release -- TOPOSORT moved from 1 to
    ! 2 to make room for ITERATIVE, matching Julia's own ordering.
    integer, parameter, public :: K24_FLUX_RECURSIVE = 0   ! RecursivePsiOut
    integer, parameter, public :: K24_FLUX_ITERATIVE = 1   ! IterativePsiOut
    integer, parameter, public :: K24_FLUX_TOPOSORT  = 2   ! TopologicalPsiOut

    ! ---------- Drainage-mode enum (par%drainage_mode) ----------
    ! Mirrors AbstractDrainageMode. Note the Q_c limits below are the
    ! *corrected* version of the paper's Sect. 3.2 description (as published it
    ! swaps them; the authors confirmed the typo by email).
    integer, parameter, public :: K24_DRAINAGE_BOTH        = 0  ! BothDrainage
    integer, parameter, public :: K24_DRAINAGE_EFFICIENT   = 1  ! EfficientOnly   (Q_c -> 0)
    integer, parameter, public :: K24_DRAINAGE_INEFFICIENT = 2  ! InefficientOnly (Q_c -> Inf)

    ! ---------- Water-thickness closure enum (par%water_thickness_algorithm) ----------
    ! Mirrors AbstractWaterThicknessAlgorithm.
    integer, parameter, public :: K24_WTHICK_DARCY_WEISBACH = 0  ! DarcyWeisbachThickness (default)
    integer, parameter, public :: K24_WTHICK_LAMINAR        = 1  ! LaminarThickness
    integer, parameter, public :: K24_WTHICK_AREAL_CONDUIT  = 2  ! ArealConduitThickness

    ! ---------- Gradient-convention enum (par%gradient_convention) ----------
    ! Mirrors AbstractGradientConvention; applies to the two sheet-flow
    ! closures (Darcy-Weisbach, laminar) only.
    integer, parameter, public :: K24_GRAD_MEAN  = 0  ! MeanGradient (default)
    integer, parameter, public :: K24_GRAD_LOCAL = 1  ! LocalGradient

    ! ---------- Sliding-law enum (par%sliding_law) ----------
    ! Mirrors AbstractSlidingLaw. NONE/WEERTMAN do not depend on N; the other
    ! two do, and switch resolve_q to the joint (q, N) Picard loop.
    integer, parameter, public :: K24_SLIDING_NONE          = 0  ! NoSlidingLaw (default)
    integer, parameter, public :: K24_SLIDING_WEERTMAN      = 1  ! WeertmanSlidingLaw
    integer, parameter, public :: K24_SLIDING_POWER_PLASTIC = 2  ! PowerPlasticSlidingLaw
    integer, parameter, public :: K24_SLIDING_REG_COULOMB   = 3  ! RegularizedCoulombSlidingLaw

    ! ============================================================
    ! Parameters (runtime, namelist-overridable)
    ! ============================================================
    ! The Julia field each entry corresponds to is named in the trailing
    ! comment. Namelist keys keep their historical Fortran spellings so
    ! existing configuration files stay valid, even where Julia's name is the
    ! clearer one (k24_manning_exponent is Glen's n; k24_bed_thickness is the
    ! bed obstacle height h_b).
    type k24_param_class

        ! -- mode selectors --
        integer  :: substrate_type              ! (Fortran-only; builds kappa)
        integer  :: flux_solver                 ! psi_out_algorithm
        integer  :: drainage_mode               ! drainage_mode
        integer  :: water_thickness_algorithm   ! water_thickness_algorithm
        integer  :: gradient_convention         ! DarcyWeisbach/Laminar gradient_convention
        integer  :: sliding_law                 ! sliding_law
        logical  :: toposort_allow_cycles       ! TopologicalPsiOut.allow_cycles

        ! -- physical constants (water/ice/gravity come from the top level) --
        real(dp) :: water_density               ! rho_w  [kg/m3]
        real(dp) :: ice_density                 ! rho_i  [kg/m3]
        real(dp) :: gravity                     ! g      [m/s2]
        real(dp) :: manning_exponent            ! n      Glen's flow-law exponent
        real(dp) :: latent_heat_water           ! L_w    [J/kg]
        real(dp) :: bed_thickness               ! h_b    bed obstacle height [m]
        real(dp) :: manning_coefficient_exponent! alpha
        real(dp) :: bed_friction_exponent       ! beta
        real(dp) :: friction_factor             ! f      Darcy-Weisbach friction factor
        real(dp) :: till_factor                 ! F_till
        real(dp) :: critical_discharge          ! Q_c    [m3/s]
        real(dp) :: initial_cavity_height       ! H_0    [m]
        real(dp) :: coupling_length             ! l_c    conduit spacing [m]
        real(dp) :: long_coupling_water         ! longcoupwater
        real(dp) :: min_pressure_fraction       ! sigmat
        real(dp) :: eta_w                       ! eta_w  [Pa s]

        ! -- clamps --
        real(dp) :: W_min                       ! Wmin  [m]
        real(dp) :: W_max                       ! Wmax  [m]
        real(dp) :: q_min                       ! q_min [m2/s]
        real(dp) :: q_max                       ! q_max [m2/s]

        ! -- solver configuration --
        integer  :: fill_iters                  ! fill_iters
        integer  :: max_psi_out_calls           ! max_psi_out_calls
        integer  :: max_dissipation_iters       ! max_dissipation_iters
        real(dp) :: dissipation_rtol            ! dissipation_rtol
        logical  :: dissipation_melt            ! dissipation_melt
        logical  :: dissipation_verbose         ! dissipation_verbose
        integer  :: max_coupling_iters          ! max_coupling_iters
        real(dp) :: coupling_rtol               ! coupling_rtol
        logical  :: coupling_verbose            ! coupling_verbose

        ! -- sliding-law parameters, one set per law so each keeps Julia's own
        !    per-law default (the velocity exponent q differs between them) --
        real(dp) :: weertman_C                  ! WeertmanSlidingLaw.C      [Pa (s/m)^q]
        real(dp) :: weertman_q                  ! WeertmanSlidingLaw.q
        real(dp) :: power_plastic_c_till        ! PowerPlasticSlidingLaw.c_till
        real(dp) :: power_plastic_q             ! PowerPlasticSlidingLaw.q
        real(dp) :: power_plastic_u0            ! PowerPlasticSlidingLaw.u0 [m/s]
        real(dp) :: reg_coulomb_c_till          ! RegularizedCoulombSlidingLaw.c_till
        real(dp) :: reg_coulomb_q               ! RegularizedCoulombSlidingLaw.q
        real(dp) :: reg_coulomb_u0              ! RegularizedCoulombSlidingLaw.u0 [m/s]

        ! -- derived, filled by k24_par_load --
        real(dp) :: K                           ! K = (2/pi)^(1/4)*sqrt((pi+2)/(rho_w*f))

    end type

    ! ============================================================
    ! Scratch workspace
    ! ============================================================
    ! Mirrors KazmierczakWorkspace. Allocated at the top of calc_k24 and
    ! released at the end, the same lifetime the previous implementation gave
    ! its local arrays -- the cost is negligible next to the FFT convolution
    ! and the flow routing, and keeping it call-scoped keeps calc_k24 free of
    ! hidden state.
    type k24_work_class
        integer :: nx = 0
        integer :: ny = 0
        ! geometric potential
        real(dp), allocatable :: phi0(:,:), phi0_tmp(:,:), h(:,:)
        real(dp), allocatable :: gx(:,:),  gy(:,:),  abs_g(:,:)      ! unsmoothed
        real(dp), allocatable :: gsx(:,:), gsy(:,:), abs_gs(:,:)     ! smoothed
        ! water flux
        real(dp), allocatable :: mdot_total(:,:), psi_out(:,:), corfac(:,:)
        real(dp), allocatable :: q_prev(:,:), N_prev(:,:), tau_b(:,:)
        integer,  allocatable :: visited(:,:)
        ! effective pressure
        real(dp), allocatable :: Q(:,:), S_inf(:,:)
        real(dp), allocatable :: H_hard(:,:), H_soft(:,:), H_cond(:,:)
        real(dp), allocatable :: N_inf(:,:), Po(:,:)
        ! routing scratch (iterative / topological solvers)
        integer,  allocatable :: in_degree(:,:)
        integer,  allocatable :: stack_i(:), stack_j(:), stack_k(:)
    end type

    ! Direction offsets, in the order accumulate_psi_out! iterates them.
    ! Accumulation is order-dependent in floating point, so this order is part
    ! of the numerics, not an implementation detail.
    integer, parameter :: K24_DIRS(2,4) = reshape([ -1, 0,  1, 0,  0, -1,  0, 1 ], [2,4])

    private
    public :: k24_param_class
    public :: k24_par_load
    public :: k24_finalize_par
    public :: initialize_kappa
    public :: calc_k24
    public :: update_psi_out                 ! dispatcher
    public :: update_psi_out_recursive       ! exposed for testing/comparison
    public :: update_psi_out_iterative       ! exposed for testing/comparison
    public :: update_psi_out_toposort        ! exposed for testing/comparison

contains

    ! ============================================================
    ! Namelist load
    ! ============================================================
    subroutine k24_par_load(par, filename, group, init)

        implicit none

        type(k24_param_class), intent(INOUT) :: par
        character(len=*),      intent(IN)    :: filename
        character(len=*),      intent(IN)    :: group
        logical, optional,     intent(IN)    :: init

        logical :: init_pars

        character(len=*), parameter :: def_file  = "input/yelmo_defaults.nml"
        character(len=*), parameter :: def_group = "yhyd"

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        ! Defaults reproduce KazmierczakHydroModel's keyword defaults
        ! (FastHydrology.jl/src/models/kazmierczak2024/model.jl).
        par%substrate_type                = K24_SUBSTRATE_HARD
        par%flux_solver                   = K24_FLUX_RECURSIVE
        par%drainage_mode                 = K24_DRAINAGE_BOTH
        par%water_thickness_algorithm     = K24_WTHICK_DARCY_WEISBACH
        par%gradient_convention           = K24_GRAD_MEAN
        par%sliding_law                   = K24_SLIDING_NONE
        par%toposort_allow_cycles         = .FALSE.

        par%water_density                 = 1000.0_dp
        par%ice_density                   =  917.0_dp
        par%gravity                       =    9.81_dp
        par%manning_exponent              =    3.0_dp
        ! L_w = 3.34e5 (KAZMIERCZAK_DEFAULT_L_W). Was 3.35e5 before this port.
        par%latent_heat_water             =    3.34e5_dp
        par%bed_thickness                 =    0.1_dp
        par%manning_coefficient_exponent  =    1.25_dp      ! alpha = 5/4
        par%bed_friction_exponent         =    1.5_dp       ! beta  = 3/2
        par%friction_factor               =    0.1_dp
        par%till_factor                   =    1.1_dp
        par%critical_discharge            =    1.0_dp
        par%initial_cavity_height         =    0.1_dp
        par%coupling_length               =    1.0e4_dp
        par%long_coupling_water           =    5.0_dp
        ! sigmat = 0 -- no N_inf floor. Pass 0.02 for KORI-ULB's own value.
        par%min_pressure_fraction         =    0.0_dp
        ! eta_w = perYear2perSecond(1.8e-3): KORI-ULB's par.waterviscosity is a
        ! per-year quantity, so it needs the same conversion as every other
        ! per-year input. Was the bare 1.8e-3 before this port -- a factor of
        ! ~3.16e7 too large.
        par%eta_w                         = 1.8e-3_dp / K24_SEC_PER_YEAR

        ! The four KORI-ULB clamps, off by default (Julia's Wmin/Wmax/q_min/
        ! q_max defaults). huge() stands in for Inf: min(huge, x) == x for any
        ! finite x, without relying on Inf surviving -Ofast.
        par%W_min                         =    0.0_dp
        par%W_max                         = huge(1.0_dp)
        par%q_min                         =    0.0_dp
        par%q_max                         = huge(1.0_dp)

        par%fill_iters                    =   10
        par%max_psi_out_calls             = 50000
        par%max_dissipation_iters         =   20
        par%dissipation_rtol              = 1.0e-12_dp
        par%dissipation_melt              = .TRUE.
        par%dissipation_verbose           = .TRUE.
        par%max_coupling_iters            =   20
        par%coupling_rtol                 = 1.0e-8_dp
        par%coupling_verbose              = .TRUE.

        ! Sliding-law parameters. C and c_till have no Julia default (they are
        ! mandatory keywords there); 0 here makes an unconfigured law a no-op
        ! rather than a silent wrong answer.
        par%weertman_C                    = 0.0_dp
        par%weertman_q                    = 1.0_dp / 3.0_dp
        par%power_plastic_c_till          = 0.0_dp
        par%power_plastic_q               = 1.0_dp
        par%power_plastic_u0              = 100.0_dp / K24_SEC_PER_YEAR
        par%reg_coulomb_c_till            = 0.0_dp
        par%reg_coulomb_q                 = 1.0_dp / 3.0_dp
        par%reg_coulomb_u0                = 100.0_dp / K24_SEC_PER_YEAR

        call nml_read(filename,group,"k24_substrate_type",                par%substrate_type,                init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_flux_solver",                   par%flux_solver,                   init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_drainage_mode",                 par%drainage_mode,                 init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_water_thickness_algorithm",     par%water_thickness_algorithm,     init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_gradient_convention",           par%gradient_convention,           init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_sliding_law",                   par%sliding_law,                   init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_toposort_allow_cycles",         par%toposort_allow_cycles,         init=init_pars,defaults_file=def_file,defaults_group=def_group)
        ! water_density, ice_density, gravity are set from top-level
        ! rho_w / rho_ice / g in hydro_par_load (single source of truth).
        call nml_read(filename,group,"k24_manning_exponent",              par%manning_exponent,              init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_latent_heat_water",             par%latent_heat_water,             init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_bed_thickness",                 par%bed_thickness,                 init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_manning_coefficient_exponent",  par%manning_coefficient_exponent,  init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_bed_friction_exponent",         par%bed_friction_exponent,         init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_friction_factor",               par%friction_factor,               init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_till_factor",                   par%till_factor,                   init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_critical_discharge",            par%critical_discharge,            init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_initial_cavity_height",         par%initial_cavity_height,         init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_coupling_length",               par%coupling_length,               init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_long_coupling_water",           par%long_coupling_water,           init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_min_pressure_fraction",         par%min_pressure_fraction,         init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_eta_w",                         par%eta_w,                         init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_W_min",                         par%W_min,                         init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_W_max",                         par%W_max,                         init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_q_min",                         par%q_min,                         init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_q_max",                         par%q_max,                         init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_fill_iters",                    par%fill_iters,                    init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_max_psi_out_calls",             par%max_psi_out_calls,             init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_max_dissipation_iters",         par%max_dissipation_iters,         init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_dissipation_rtol",              par%dissipation_rtol,              init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_dissipation_melt",              par%dissipation_melt,              init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_dissipation_verbose",           par%dissipation_verbose,           init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_max_coupling_iters",            par%max_coupling_iters,            init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_coupling_rtol",                 par%coupling_rtol,                 init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_coupling_verbose",              par%coupling_verbose,              init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_weertman_C",                    par%weertman_C,                    init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_weertman_q",                    par%weertman_q,                    init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_power_plastic_c_till",          par%power_plastic_c_till,          init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_power_plastic_q",               par%power_plastic_q,               init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_power_plastic_u0",              par%power_plastic_u0,              init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_reg_coulomb_c_till",            par%reg_coulomb_c_till,            init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_reg_coulomb_q",                 par%reg_coulomb_q,                 init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_reg_coulomb_u0",                par%reg_coulomb_u0,                init=init_pars,defaults_file=def_file,defaults_group=def_group)

        ! Derived. hydro_par_load overwrites water_density afterwards with the
        ! top-level RHO_W and calls k24_finalize_par to refresh K.
        call k24_finalize_par(par)

        return

    end subroutine k24_par_load

    subroutine k24_finalize_par(par)
        ! Recompute every parameter derived from another one. Called at the end
        ! of k24_par_load and again by hydro_par_load once it has stamped the
        ! top-level rho_w / rho_ice / g into the sub-struct, so K always
        ! reflects the densities actually in use.
        implicit none
        type(k24_param_class), intent(INOUT) :: par

        ! Manning-Strickler / Darcy-Weisbach conductivity coefficient,
        ! K = (2/pi)^(1/4) * sqrt((pi+2)/(rho_w*f)).
        par%K = (2.0_dp / K24_PI)**0.25_dp * &
                sqrt((K24_PI + 2.0_dp) / (par%water_density * par%friction_factor))

    end subroutine k24_finalize_par

    ! ============================================================
    ! Substrate indicator (one-time, called from hydro_init_state)
    ! ============================================================
    subroutine initialize_kappa(kappa, b, substrate_type)

        implicit none

        real(dp),   intent(OUT) :: kappa(:,:)
        real(dp),   intent(IN)  :: b(:,:)
        integer,    intent(IN)  :: substrate_type

        integer :: i, j, nx, ny

        nx = size(kappa,1)
        ny = size(kappa,2)

        select case (substrate_type)
            case (K24_SUBSTRATE_HARD)
                kappa = 0.0_dp
            case (K24_SUBSTRATE_SOFT)
                kappa = 1.0_dp
            case (K24_SUBSTRATE_MIXED)
                !$omp parallel do default(shared) private(i,j) schedule(static)
                do j = 1, ny
                    do i = 1, nx
                        if (b(i,j) < -1000.0_dp) then
                            kappa(i,j) = 1.0_dp
                        else
                            kappa(i,j) = 0.0_dp
                        end if
                    end do
                end do
                !$omp end parallel do
            case default
                write(*,*) "initialize_kappa:: error: substrate_type must be one of [0,1,2]."
                write(*,*) "substrate_type = ", substrate_type
                stop
        end select

        return

    end subroutine initialize_kappa

    ! ============================================================
    ! Top-level K24 driver
    ! ============================================================
    subroutine calc_k24(q_x, q_y, N, p_w, W, q, &
                        H_ice, z_bed, mask, mdot, uxy_b, A_glen, kappa, &
                        dx, dy, par, gsx_out, gsy_out, absgs_out, absg_out, phi0_out)
        ! Mirrors update_steady_state! (FastHydrology.jl/.../run.jl):
        ! update_q, then update_N, then update_W -- W last because its
        ! ArealConduitThickness closure reads S_inf, which update_N is what
        ! keeps current.
        !
        ! `q` and `N` are INOUT, not OUT: FastHydrology.jl's model.q and
        ! state.N persist between solves, and both Picard loops below warm-start
        ! from them. Pass the caller's own persistent fields (hyd%now%q,
        ! hyd%now%N), zero-initialised on the first call.
        !
        ! The source `mdot` is the water source rate fed into the transport
        ! solver [m/s, water-equivalent]. In sequential bucket->K24 coupling
        ! this is the bucket overflow (till-saturation spill); in TIL_NONE mode
        ! it can be the raw basal-melt water equivalent.
        !
        ! The trailing optional arguments expose internals that are otherwise
        ! scoped to this call. They exist so tests/k24_synth.f90 can compare
        ! the intermediate fields against FastHydrology.jl stage by stage
        ! rather than only the end result -- which is how the drainage-mode
        ! opening-coefficient bug was found. Ordinary callers omit them and
        ! pay nothing.

        implicit none

        real(dp), intent(OUT)   :: q_x(:,:), q_y(:,:)   ! water flux components [m2/s]
        real(dp), intent(INOUT) :: N(:,:)               ! effective pressure [Pa]
        real(dp), intent(OUT)   :: p_w(:,:)             ! water pressure (Po - N) [Pa]
        real(dp), intent(OUT)   :: W(:,:)               ! water layer thickness [m]
        real(dp), intent(INOUT) :: q(:,:)               ! distributed flux magnitude [m2/s]
        real(dp), intent(IN)    :: H_ice(:,:)           ! [m]
        real(dp), intent(IN)    :: z_bed(:,:)           ! [m]
        real(dp), intent(IN)    :: mask(:,:)            ! 1 on grounded ice, else 0
        real(dp), intent(IN)    :: mdot(:,:)            ! source rate [m/s, water-equiv.]
        real(dp), intent(IN)    :: uxy_b(:,:)           ! basal sliding speed magnitude [m/s]
        real(dp), intent(IN)    :: A_glen(:,:)          ! Glen's A [Pa^-n s^-1]
        real(dp), intent(IN)    :: kappa(:,:)           ! bed type indicator (0 hard, 1 soft)
        real(dp), intent(IN)    :: dx, dy               ! [m] grid spacing
        type(k24_param_class), intent(IN) :: par
        ! Optional diagnostics (see above): the smoothed gradient components
        ! and magnitude, the unsmoothed magnitude, and the filled potential.
        real(dp), intent(OUT), optional :: gsx_out(:,:), gsy_out(:,:), absgs_out(:,:)
        real(dp), intent(OUT), optional :: absg_out(:,:), phi0_out(:,:)

        type(k24_work_class) :: wk
        integer  :: nx, ny, i, j
        real(dp) :: gmag

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        call k24_work_alloc(wk, nx, ny)

        ! ================= update_q! =================

        ! Geometric potential from the RAW ice thickness.
        call update_phi0(wk%phi0, H_ice, z_bed, par)

        ! Fill local minima so water does not get stuck.
        call potential_filling(wk%phi0, wk%phi0_tmp, par%fill_iters)

        ! Ice thickness consistent with the filled potential. Stored separately
        ! so it does not leak into the effective-pressure calculation, which
        ! keeps using the raw H_ice (see update_Po below).
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
            do i = 1, nx
                wk%h(i,j) = (wk%phi0(i,j) - par%water_density * par%gravity * z_bed(i,j)) &
                          / (par%ice_density * par%gravity)
            end do
        end do
        !$omp end parallel do

        call update_potential_gradients(wk, dx, dy)
        call update_smoothed_potential_gradients(wk, dx, dy, mask, par)

        ! Correction factor from psi_out to q. Depends only on the (already
        ! updated) smoothed gradients, so it is fixed for the whole Picard
        ! loop below. Anisotropy-aware: the per-cell outflow width is
        ! |gsx|*dy + |gsy|*dx, exact for dx /= dy (this replaces the earlier
        ! sqrt(dx*dy) placeholder, which was only correct for square cells).
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
            do i = 1, nx
                wk%corfac(i,j) = (abs(wk%gsx(i,j)) * dy + abs(wk%gsy(i,j)) * dx) &
                    / (sqrt(wk%gsx(i,j)*wk%gsx(i,j) + wk%gsy(i,j)*wk%gsy(i,j)) + 1.0e-15_dp)
            end do
        end do
        !$omp end parallel do

        call resolve_q(q, N, wk, mask, mdot, uxy_b, A_glen, kappa, H_ice, dx, dy, par)

        ! ================= update_N! =================
        call update_N(N, q, wk, uxy_b, A_glen, kappa, H_ice, par)

        ! ================= update_W! =================
        call update_W(W, q, wk, mask, par)

        ! ================= diagnostics =================
        ! q_x / q_y have no counterpart in FastHydrology.jl, which carries only
        ! the scalar q. They are a Fortran-side diagnostic: the magnitude q
        ! resolved along the UNSMOOTHED potential gradient, unchanged from the
        ! previous implementation. (The routing itself uses the smoothed
        ! gradient; if you need the components to point along the routing
        ! direction, use gsx/gsy here instead.)
        !$omp parallel do default(shared) private(i,j,gmag) schedule(static)
        do j = 1, ny
            do i = 1, nx
                gmag = wk%abs_g(i,j)
                if (gmag > 1.0e-12_dp) then
                    q_x(i,j) = q(i,j) * wk%gx(i,j) / gmag
                    q_y(i,j) = q(i,j) * wk%gy(i,j) / gmag
                else
                    q_x(i,j) = 0.0_dp
                    q_y(i,j) = 0.0_dp
                end if
                p_w(i,j) = wk%Po(i,j) - N(i,j)
            end do
        end do
        !$omp end parallel do

        if (present(gsx_out))   gsx_out   = wk%gsx
        if (present(gsy_out))   gsy_out   = wk%gsy
        if (present(absgs_out)) absgs_out = wk%abs_gs
        if (present(absg_out))  absg_out  = wk%abs_g
        if (present(phi0_out))  phi0_out  = wk%phi0

        call k24_work_free(wk)

        return

    end subroutine calc_k24

    ! ============================================================
    ! Workspace management
    ! ============================================================
    subroutine k24_work_alloc(wk, nx, ny)
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        integer,              intent(IN)    :: nx, ny

        call k24_work_free(wk)

        wk%nx = nx
        wk%ny = ny

        allocate(wk%phi0(nx,ny), wk%phi0_tmp(nx,ny), wk%h(nx,ny))
        allocate(wk%gx(nx,ny),  wk%gy(nx,ny),  wk%abs_g(nx,ny))
        allocate(wk%gsx(nx,ny), wk%gsy(nx,ny), wk%abs_gs(nx,ny))
        allocate(wk%mdot_total(nx,ny), wk%psi_out(nx,ny), wk%corfac(nx,ny))
        allocate(wk%q_prev(nx,ny), wk%N_prev(nx,ny), wk%tau_b(nx,ny))
        allocate(wk%visited(nx,ny))
        allocate(wk%Q(nx,ny), wk%S_inf(nx,ny))
        allocate(wk%H_hard(nx,ny), wk%H_soft(nx,ny), wk%H_cond(nx,ny))
        allocate(wk%N_inf(nx,ny), wk%Po(nx,ny))

        wk%psi_out = 0.0_dp
        wk%tau_b   = 0.0_dp

    end subroutine k24_work_alloc

    subroutine k24_work_free(wk)
        implicit none
        type(k24_work_class), intent(INOUT) :: wk

        if (allocated(wk%phi0))       deallocate(wk%phi0)
        if (allocated(wk%phi0_tmp))   deallocate(wk%phi0_tmp)
        if (allocated(wk%h))          deallocate(wk%h)
        if (allocated(wk%gx))         deallocate(wk%gx)
        if (allocated(wk%gy))         deallocate(wk%gy)
        if (allocated(wk%abs_g))      deallocate(wk%abs_g)
        if (allocated(wk%gsx))        deallocate(wk%gsx)
        if (allocated(wk%gsy))        deallocate(wk%gsy)
        if (allocated(wk%abs_gs))     deallocate(wk%abs_gs)
        if (allocated(wk%mdot_total)) deallocate(wk%mdot_total)
        if (allocated(wk%psi_out))    deallocate(wk%psi_out)
        if (allocated(wk%corfac))     deallocate(wk%corfac)
        if (allocated(wk%q_prev))     deallocate(wk%q_prev)
        if (allocated(wk%N_prev))     deallocate(wk%N_prev)
        if (allocated(wk%tau_b))      deallocate(wk%tau_b)
        if (allocated(wk%visited))    deallocate(wk%visited)
        if (allocated(wk%Q))          deallocate(wk%Q)
        if (allocated(wk%S_inf))      deallocate(wk%S_inf)
        if (allocated(wk%H_hard))     deallocate(wk%H_hard)
        if (allocated(wk%H_soft))     deallocate(wk%H_soft)
        if (allocated(wk%H_cond))     deallocate(wk%H_cond)
        if (allocated(wk%N_inf))      deallocate(wk%N_inf)
        if (allocated(wk%Po))         deallocate(wk%Po)
        if (allocated(wk%in_degree))  deallocate(wk%in_degree)
        if (allocated(wk%stack_i))    deallocate(wk%stack_i)
        if (allocated(wk%stack_j))    deallocate(wk%stack_j)
        if (allocated(wk%stack_k))    deallocate(wk%stack_k)

        wk%nx = 0
        wk%ny = 0

    end subroutine k24_work_free

    ! ============================================================
    ! Masked reductions (mirroring grid.jl's masked_* helpers)
    ! ============================================================
    function masked_mean(field, mask) result(m)
        implicit none
        real(dp), intent(IN) :: field(:,:), mask(:,:)
        real(dp) :: m
        integer  :: i, j, nx, ny, cnt
        real(dp) :: s

        nx = size(field,1); ny = size(field,2)
        s = 0.0_dp; cnt = 0
        !$omp parallel do default(shared) private(i,j) reduction(+:s,cnt) schedule(static)
        do j = 1, ny
            do i = 1, nx
                if (mask(i,j) == 1.0_dp) then
                    s   = s + field(i,j)
                    cnt = cnt + 1
                end if
            end do
        end do
        !$omp end parallel do

        ! DEVIATION: Julia divides by cnt unconditionally and returns NaN for
        ! an empty mask. Returning 0 keeps a fully-floating domain finite.
        if (cnt > 0) then
            m = s / real(cnt, dp)
        else
            m = 0.0_dp
        end if
    end function masked_mean

    function masked_max_abs(field, mask) result(m)
        implicit none
        real(dp), intent(IN) :: field(:,:), mask(:,:)
        real(dp) :: m
        integer  :: i, j, nx, ny

        nx = size(field,1); ny = size(field,2)
        m = 0.0_dp
        !$omp parallel do default(shared) private(i,j) reduction(max:m) schedule(static)
        do j = 1, ny
            do i = 1, nx
                if (mask(i,j) == 1.0_dp) m = max(m, abs(field(i,j)))
            end do
        end do
        !$omp end parallel do
    end function masked_max_abs

    function masked_max_abs_diff(a, b, mask) result(m)
        implicit none
        real(dp), intent(IN) :: a(:,:), b(:,:), mask(:,:)
        real(dp) :: m
        integer  :: i, j, nx, ny

        nx = size(a,1); ny = size(a,2)
        m = 0.0_dp
        !$omp parallel do default(shared) private(i,j) reduction(max:m) schedule(static)
        do j = 1, ny
            do i = 1, nx
                if (mask(i,j) == 1.0_dp) m = max(m, abs(a(i,j) - b(i,j)))
            end do
        end do
        !$omp end parallel do
    end function masked_max_abs_diff

    ! ============================================================
    ! Hydraulic potential
    ! ============================================================
    subroutine update_phi0(phi0, h, b, par)
        implicit none
        real(dp), intent(OUT) :: phi0(:,:)
        real(dp), intent(IN)  :: h(:,:), b(:,:)
        type(k24_param_class), intent(IN) :: par
        integer :: i, j, nx, ny

        nx = size(phi0,1); ny = size(phi0,2)
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
            do i = 1, nx
                phi0(i,j) = par%ice_density   * par%gravity * h(i,j) &
                          + par%water_density * par%gravity * b(i,j)
            end do
        end do
        !$omp end parallel do
    end subroutine update_phi0

    ! ============================================================
    ! Iterative hollow-filling for spurious sinks
    ! ============================================================
    subroutine potential_filling(phi0, phi0_tmp, iterations)
        ! Mirrors potential_filling! (water_flux.jl). Every cell is visited,
        ! including the domain edges, whose out-of-range neighbours are
        ! edge-replicated (index clamped) -- the same convention
        ! minus_gradient_x!/minus_gradient_y! use. The previous Fortran
        ! implementation skipped the border ring entirely.
        implicit none
        real(dp), intent(INOUT) :: phi0(:,:)
        real(dp), intent(INOUT) :: phi0_tmp(:,:)
        integer,  intent(IN)    :: iterations

        integer  :: iter, i, j, nx, ny, im1, ip1, jm1, jp1
        real(dp) :: p, p1, p2, p3, p4

        nx = size(phi0,1); ny = size(phi0,2)

        phi0_tmp = phi0

        do iter = 1, iterations
            !$omp parallel do default(shared) private(i,j,im1,ip1,jm1,jp1,p,p1,p2,p3,p4) schedule(static)
            do j = 1, ny
                do i = 1, nx
                    p   = phi0(i,j)
                    im1 = max(i-1, 1);  ip1 = min(i+1, nx)
                    jm1 = max(j-1, 1);  jp1 = min(j+1, ny)
                    p1  = phi0(ip1,j);  p2 = phi0(im1,j)
                    p3  = phi0(i,jp1);  p4 = phi0(i,jm1)
                    if (p < p1 .and. p < p2 .and. p < p3 .and. p < p4) then
                        phi0_tmp(i,j) = (p1 + p2 + p3 + p4) / 4.0_dp
                    end if
                end do
            end do
            !$omp end parallel do
            phi0 = phi0_tmp
        end do

    end subroutine potential_filling

    ! ============================================================
    ! Potential gradients (unsmoothed)
    ! ============================================================
    subroutine update_potential_gradients(wk, dx, dy)
        ! Central difference with the neighbour INDEX clamped at the domain
        ! edge (minus_gradient_x!/minus_gradient_y! in grid.jl), i.e. an edge
        ! cell differences itself against its single interior neighbour over
        ! the full 2*dx. The previous Fortran implementation instead copied the
        ! adjacent interior cell's gradient into the border ring, which is a
        ! different value.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: dx, dy

        integer :: i, j, nx, ny, im1, ip1, jm1, jp1

        nx = wk%nx; ny = wk%ny

        !$omp parallel do default(shared) private(i,j,im1,ip1) schedule(static)
        do j = 1, ny
            do i = 1, nx
                im1 = max(i-1, 1);  ip1 = min(i+1, nx)
                wk%gx(i,j) = -(wk%phi0(ip1,j) - wk%phi0(im1,j)) / (2.0_dp * dx)
            end do
        end do
        !$omp end parallel do

        !$omp parallel do default(shared) private(i,j,jm1,jp1) schedule(static)
        do j = 1, ny
            jm1 = max(j-1, 1);  jp1 = min(j+1, ny)
            do i = 1, nx
                wk%gy(i,j) = -(wk%phi0(i,jp1) - wk%phi0(i,jm1)) / (2.0_dp * dy)
            end do
        end do
        !$omp end parallel do

        ! Euclidean magnitude for the unsmoothed field (the smoothed one below
        ! uses the L1 sum instead -- that asymmetry is Julia's, not a typo).
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
            do i = 1, nx
                wk%abs_g(i,j) = sqrt(wk%gx(i,j)*wk%gx(i,j) + wk%gy(i,j)*wk%gy(i,j))
            end do
        end do
        !$omp end parallel do

    end subroutine update_potential_gradients

    ! ============================================================
    ! Stress-gradient-coupling smoothing of the potential gradients
    ! ============================================================
    subroutine update_smoothed_potential_gradients(wk, dx, dy, mask, par)
        ! Mirrors update_smoothed_potential_gradients! (water_flux.jl).
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: dx, dy, mask(:,:)
        type(k24_param_class),intent(IN)    :: par

        real(dp), allocatable :: kernel(:,:)
        real(dp) :: h_avg, scale, width, delta_min, dist, kernel_sum
        integer  :: maxlevel_x, maxlevel_y, frb_x, frb_y, i, j, ni, nj, nx, ny

        nx = wk%nx; ny = wk%ny

        ! longcoupwater == 0 disables the smoothing entirely. Use this when the
        ! coupling length is smaller than one grid cell (e.g. 16-32 km grids).
        if (par%long_coupling_water == 0.0_dp) then
            wk%gsx = wk%gx
            wk%gsy = wk%gy
            !$omp parallel do default(shared) private(i,j) schedule(static)
            do j = 1, ny
                do i = 1, nx
                    wk%abs_gs(i,j) = abs(wk%gx(i,j)) + abs(wk%gy(i,j))
                end do
            end do
            !$omp end parallel do
            return
        end if

        ! Mean grounded-ice thickness, from the potential-filled h.
        h_avg = max(masked_mean(wk%h, mask), 10.0_dp)

        scale = h_avg * par%long_coupling_water * 2.0_dp

        ! Radius of the cone base. The effective coupling length is the
        ! kernel's weighted mean distance from the centre, width/3 =
        ! (4/3)*h_avg*longcoupwater -- about 6.7x ice thickness at
        ! longcoupwater = 5, consistent with Kamb & Echelmeyer (1986).
        width = 2.0_dp * scale

        delta_min = min(dx, dy)
        if (width <= delta_min) then
            ! Bump the cone scale so the kernel spans at least ~1 cell in the
            ! tighter direction instead of degenerating to a no-op.
            ! NOTE: `width` is deliberately NOT recomputed here, matching
            ! FastHydrology.jl exactly. The previous Fortran implementation did
            ! recompute it, which made the kernel a different size in this
            ! (rare, very-fine-grid) branch.
            scale = delta_min / 2.0_dp + 1.0_dp
        end if

        ! Kernel size, sized independently per axis so a dx /= dy grid does not
        ! pay for the finer axis's resolution in both directions.
        maxlevel_x = 2 * round_half_even(width / dx - 0.5_dp) + 1
        maxlevel_y = 2 * round_half_even(width / dy - 0.5_dp) + 1
        frb_x = (maxlevel_x - 1) / 2
        frb_y = (maxlevel_y - 1) / 2

        allocate(kernel(maxlevel_x, maxlevel_y))
        do nj = 1, maxlevel_y
            do ni = 1, maxlevel_x
                ! True physical (Euclidean) distance from the kernel centre,
                ! using dx and dy separately, so the support is a circle in
                ! physical space for any cell aspect ratio.
                dist = sqrt( (dx * real(ni - frb_x - 1, dp))**2 + &
                             (dy * real(nj - frb_y - 1, dp))**2 ) / scale
                kernel(ni,nj) = max(0.0_dp, 1.0_dp - dist / 2.0_dp)
            end do
        end do
        kernel_sum = sum(kernel)
        if (kernel_sum > 0.0_dp) kernel = kernel / kernel_sum

        call imfilter_replicate_fftw(wk%gx, nx, ny, kernel, frb_x, frb_y, wk%gsx)
        call imfilter_replicate_fftw(wk%gy, nx, ny, kernel, frb_x, frb_y, wk%gsy)

        ! L1 magnitude, matching abs_grad_phi0_s in water_flux.jl.
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
            do i = 1, nx
                wk%abs_gs(i,j) = abs(wk%gsx(i,j)) + abs(wk%gsy(i,j))
            end do
        end do
        !$omp end parallel do

        deallocate(kernel)

    end subroutine update_smoothed_potential_gradients

    integer function round_half_even(x) result(r)
        ! Round half to even, matching Julia's round(Int, x). Fortran's nint
        ! rounds half away from zero, which would disagree on the exact .5
        ! boundary that decides the kernel size.
        implicit none
        real(dp), intent(IN) :: x
        real(dp) :: f, diff

        f    = floor(x)
        diff = x - f
        if (diff > 0.5_dp) then
            r = int(f) + 1
        else if (diff < 0.5_dp) then
            r = int(f)
        else
            ! Exactly .5: pick the even neighbour.
            if (modulo(int(f), 2) == 0) then
                r = int(f)
            else
                r = int(f) + 1
            end if
        end if
    end function round_half_even

    ! ============================================================
    ! 2-D convolution via FFTW3, "replicate" border.
    ! ============================================================
    subroutine imfilter_replicate_fftw(input, nx, ny, kernel, frb_x, frb_y, output)
        ! Cached_fft_convolve! (fft_convolution.jl) in Fortran: embed the input
        ! in an (nx+2*frb_x, ny+2*frb_y) array with REPLICATE padding (nearest
        ! edge extended), wrap the centred kernel around the array's origin,
        ! multiply in the frequency domain, and crop back.
        !
        ! Two fixes relative to the previous implementation:
        !   * the padding was "reflect"; ImageFiltering's default -- which
        !     FastHydrology.jl reproduces -- is "replicate".
        !   * the crop offset was (i + 2*frb, j + 2*frb). The circular
        !     convolution places the centred result at (i + frb, j + frb), so
        !     the smoothed field came out translated by frb cells.
        !
        ! The FFT array is zero-padded from (nx+2*frb_x, ny+2*frb_y) up to the
        ! next 5-smooth size purely for FFTW speed. That cannot change the
        ! result: for an output cell p in [frb+1, frb+n] the circular sum only
        ! reads indices p-d for |d| <= frb, i.e. [1, n+2*frb] -- always real
        ! data, never the zeros or a wrapped-around neighbour.
        !
        ! Link with -lfftw3.
        use, intrinsic :: iso_c_binding
        implicit none
        include 'fftw3.f03'

        integer,  intent(IN)  :: nx, ny, frb_x, frb_y
        real(dp), intent(IN)  :: input(nx, ny)
        real(dp), intent(IN)  :: kernel(2*frb_x+1, 2*frb_y+1)
        real(dp), intent(OUT) :: output(nx, ny)

        integer :: Npx, Npy, Mfft, Nfft, i, j, ii, jj, ni, nj
        real(C_DOUBLE),            allocatable :: work_a(:,:), work_b(:,:), work_out(:,:)
        complex(C_DOUBLE_COMPLEX), allocatable :: A_hat(:,:), B_hat(:,:)
        type(C_PTR) :: plan_fwd_a, plan_fwd_b, plan_bwd

        Npx = nx + 2*frb_x
        Npy = ny + 2*frb_y

        Mfft = next_smooth_size(Npx)
        Nfft = next_smooth_size(Npy)

        allocate(work_a(Mfft, Nfft), work_b(Mfft, Nfft), work_out(Mfft, Nfft))
        allocate(A_hat(Mfft/2+1, Nfft), B_hat(Mfft/2+1, Nfft))
        work_a = 0.0_dp
        work_b = 0.0_dp

        ! Replicate padding: clamp the source index to the domain.
        !$omp parallel do default(shared) private(i,j,ii,jj) schedule(static)
        do j = 1, Npy
            do i = 1, Npx
                ii = min(max(i - frb_x, 1), nx)
                jj = min(max(j - frb_y, 1), ny)
                work_a(i, j) = input(ii, jj)
            end do
        end do
        !$omp end parallel do

        ! Kernel offset (di,dj) goes to wrapped index mod(di,M)+1.
        do nj = 1, 2*frb_y+1
            do ni = 1, 2*frb_x+1
                ii = modulo(ni - frb_x - 1, Mfft) + 1
                jj = modulo(nj - frb_y - 1, Nfft) + 1
                work_b(ii, jj) = kernel(ni, nj)
            end do
        end do

        plan_fwd_a = fftw_plan_dft_r2c_2d(Nfft, Mfft, work_a, A_hat, FFTW_ESTIMATE)
        call fftw_execute_dft_r2c(plan_fwd_a, work_a, A_hat)
        call fftw_destroy_plan(plan_fwd_a)

        plan_fwd_b = fftw_plan_dft_r2c_2d(Nfft, Mfft, work_b, B_hat, FFTW_ESTIMATE)
        call fftw_execute_dft_r2c(plan_fwd_b, work_b, B_hat)
        call fftw_destroy_plan(plan_fwd_b)

        A_hat = A_hat * B_hat

        plan_bwd = fftw_plan_dft_c2r_2d(Nfft, Mfft, A_hat, work_out, FFTW_ESTIMATE)
        call fftw_execute_dft_c2r(plan_bwd, A_hat, work_out)
        call fftw_destroy_plan(plan_bwd)

        work_out = work_out / real(Mfft, dp) / real(Nfft, dp)

        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
            do i = 1, nx
                output(i, j) = work_out(i + frb_x, j + frb_y)
            end do
        end do
        !$omp end parallel do

        deallocate(work_a, work_b, work_out, A_hat, B_hat)
    end subroutine imfilter_replicate_fftw

    integer function next_smooth_size(n) result(m)
        ! Smallest m >= n whose only prime factors are 2, 3 and 5 -- the sizes
        ! FFTW has dedicated codelets for. Much less padding than rounding up
        ! to a power of two (never more than ~20% over n, versus up to 100%).
        implicit none
        integer, intent(IN) :: n
        integer :: k

        m = max(n, 1)
        do
            k = m
            do while (modulo(k, 2) == 0); k = k / 2; end do
            do while (modulo(k, 3) == 0); k = k / 3; end do
            do while (modulo(k, 5) == 0); k = k / 5; end do
            if (k == 1) return
            m = m + 1
        end do
    end function next_smooth_size

    ! ============================================================
    ! Sliding law -> basal shear stress
    ! ============================================================
    subroutine update_tau_b(tau_b, N, uxy_b, par)
        ! Mirrors update_tau_b! (sliding_law.jl). tau_b feeds the frictional
        ! heating term tau_b*v_b/(L_w*rho_w) of the melt rate (Eq. 3, Sec.
        ! 2.2.1 of Kazmierczak et al 2024).
        implicit none
        real(dp), intent(OUT) :: tau_b(:,:)
        real(dp), intent(IN)  :: N(:,:), uxy_b(:,:)
        type(k24_param_class), intent(IN) :: par

        integer  :: i, j, nx, ny
        real(dp) :: C, qe, c_till, u0

        nx = size(tau_b,1); ny = size(tau_b,2)

        select case (par%sliding_law)

            case (K24_SLIDING_NONE)
                tau_b = 0.0_dp

            case (K24_SLIDING_WEERTMAN)
                C = par%weertman_C; qe = par%weertman_q
                !$omp parallel do default(shared) private(i,j) schedule(static)
                do j = 1, ny
                    do i = 1, nx
                        tau_b(i,j) = C * uxy_b(i,j)**qe
                    end do
                end do
                !$omp end parallel do

            case (K24_SLIDING_POWER_PLASTIC)
                c_till = par%power_plastic_c_till
                qe     = par%power_plastic_q
                u0     = par%power_plastic_u0
                !$omp parallel do default(shared) private(i,j) schedule(static)
                do j = 1, ny
                    do i = 1, nx
                        tau_b(i,j) = c_till * N(i,j) * (uxy_b(i,j) / u0)**qe
                    end do
                end do
                !$omp end parallel do

            case (K24_SLIDING_REG_COULOMB)
                c_till = par%reg_coulomb_c_till
                qe     = par%reg_coulomb_q
                u0     = par%reg_coulomb_u0
                !$omp parallel do default(shared) private(i,j) schedule(static)
                do j = 1, ny
                    do i = 1, nx
                        tau_b(i,j) = c_till * N(i,j) * &
                                     (uxy_b(i,j) / (uxy_b(i,j) + u0))**qe
                    end do
                end do
                !$omp end parallel do

            case default
                write(*,*) "update_tau_b:: error: k24_sliding_law must be one of [0,1,2,3]."
                write(*,*) "sliding_law = ", par%sliding_law
                stop

        end select

    end subroutine update_tau_b

    ! ============================================================
    ! Water-flux fixed point
    ! ============================================================
    subroutine resolve_q(q, N, wk, mask, mdot, uxy_b, A_glen, kappa, H_ice, dx, dy, par)
        ! Mirrors the three resolve_q! methods (water_flux.jl):
        !
        !   * N-independent law (NONE/WEERTMAN), dissipation off -- one pass.
        !   * N-independent law, dissipation on -- Picard on q alone.
        !   * N-dependent law (POWER_PLASTIC/REG_COULOMB) -- joint (q, N)
        !     Picard, regardless of the dissipation setting.
        implicit none
        real(dp),             intent(INOUT) :: q(:,:), N(:,:)
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), mdot(:,:), uxy_b(:,:)
        real(dp),             intent(IN)    :: A_glen(:,:), kappa(:,:), H_ice(:,:)
        real(dp),             intent(IN)    :: dx, dy
        type(k24_param_class),intent(IN)    :: par

        logical :: pressure_dependent

        pressure_dependent = (par%sliding_law == K24_SLIDING_POWER_PLASTIC) .or. &
                             (par%sliding_law == K24_SLIDING_REG_COULOMB)

        if (pressure_dependent) then
            call resolve_q_coupled(q, N, wk, mask, mdot, uxy_b, A_glen, kappa, H_ice, dx, dy, par)
        else if (par%dissipation_melt) then
            call resolve_q_dissipation(q, N, wk, mask, mdot, uxy_b, dx, dy, par)
        else
            call resolve_q_single(q, N, wk, mask, mdot, uxy_b, dx, dy, par)
        end if

    end subroutine resolve_q

    subroutine set_q_from_psi_out(q, wk, par)
        implicit none
        real(dp),             intent(OUT) :: q(:,:)
        type(k24_work_class), intent(IN)  :: wk
        type(k24_param_class),intent(IN)  :: par
        integer :: i, j

        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, wk%ny
            do i = 1, wk%nx
                ! DEVIATION: corfac is exactly zero only where the smoothed
                ! gradient vanishes in both directions, i.e. where there is no
                ! flow direction at all. Julia divides regardless and lets the
                ! cell go to Inf or NaN; taking q = 0 there is the physical
                ! limit and keeps -Ofast from producing something arbitrary.
                if (wk%corfac(i,j) > 0.0_dp) then
                    q(i,j) = min(max(wk%psi_out(i,j) / wk%corfac(i,j), par%q_min), par%q_max)
                else
                    q(i,j) = min(max(0.0_dp, par%q_min), par%q_max)
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine set_q_from_psi_out

    subroutine resolve_q_single(q, N, wk, mask, mdot, uxy_b, dx, dy, par)
        ! Dissipation off and an N-independent sliding law: the water source
        ! depends on neither q nor N, so a single routing pass is exact.
        implicit none
        real(dp),             intent(INOUT) :: q(:,:)
        real(dp),             intent(IN)    :: N(:,:)
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), mdot(:,:), uxy_b(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par

        integer  :: i, j
        real(dp) :: rL

        call update_tau_b(wk%tau_b, N, uxy_b, par)

        rL = par%latent_heat_water * par%water_density
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, wk%ny
            do i = 1, wk%nx
                wk%mdot_total(i,j) = mdot(i,j) + wk%tau_b(i,j) * uxy_b(i,j) / rL
            end do
        end do
        !$omp end parallel do

        call update_psi_out(wk, mask, dx, dy, par)
        call set_q_from_psi_out(q, wk, par)

    end subroutine resolve_q_single

    subroutine resolve_q_dissipation(q, N, wk, mask, mdot, uxy_b, dx, dy, par)
        ! Dissipation on and an N-independent sliding law: mdot_total depends
        ! on q through |q*grad(phi0)|/(L_w*rho_w) only, so Picard-iterate on q.
        implicit none
        real(dp),             intent(INOUT) :: q(:,:)
        real(dp),             intent(IN)    :: N(:,:)
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), mdot(:,:), uxy_b(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par

        integer  :: i, j, iter, n_iters
        real(dp) :: rL, q_scale
        logical  :: converged

        call update_tau_b(wk%tau_b, N, uxy_b, par)

        rL        = par%latent_heat_water * par%water_density
        converged = .FALSE.
        n_iters   = par%max_dissipation_iters

        do iter = 1, par%max_dissipation_iters

            wk%q_prev = q

            ! Total source: basal melt, the (fixed) frictional-heating term,
            ! and the dissipation melt from the current q. Zero on the first
            ! sweep of a cold start, since q begins at zero.
            !$omp parallel do default(shared) private(i,j) schedule(static)
            do j = 1, wk%ny
                do i = 1, wk%nx
                    wk%mdot_total(i,j) = mdot(i,j) &
                        + wk%tau_b(i,j) * uxy_b(i,j) / rL &
                        + abs(q(i,j) * wk%abs_g(i,j)) / rL
                end do
            end do
            !$omp end parallel do

            call update_psi_out(wk, mask, dx, dy, par)
            call set_q_from_psi_out(q, wk, par)

            q_scale = max(masked_max_abs(q, mask), 1.0e-15_dp)
            if (masked_max_abs_diff(q, wk%q_prev, mask) <= par%dissipation_rtol * q_scale) then
                converged = .TRUE.
                n_iters   = iter
                exit
            end if

        end do

        if (par%dissipation_verbose) then
            if (converged) then
                write(*,'(a,i0,a)') " k24: water flux Picard loop converged after ", n_iters, " iteration(s)."
            else
                write(*,'(a,i0,a)') " k24: water flux Picard loop did NOT converge (hit max_dissipation_iters = ", n_iters, ")."
            end if
        end if

    end subroutine resolve_q_dissipation

    subroutine resolve_q_coupled(q, N, wk, mask, mdot, uxy_b, A_glen, kappa, H_ice, dx, dy, par)
        ! N-dependent sliding law: tau_b depends on N, which is downstream of
        ! q, so q and N form a joint fixed point. Each sweep recomputes tau_b
        ! from the current N, routes q, then refreshes N from the new q.
        !
        ! N starts from whatever the caller passed in (zero on a cold start),
        ! so the first sweep's tau_b is zero and ramps up -- ordinary Picard
        ! behaviour, not a bug.
        implicit none
        real(dp),             intent(INOUT) :: q(:,:), N(:,:)
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), mdot(:,:), uxy_b(:,:)
        real(dp),             intent(IN)    :: A_glen(:,:), kappa(:,:), H_ice(:,:)
        real(dp),             intent(IN)    :: dx, dy
        type(k24_param_class),intent(IN)    :: par

        integer  :: i, j, iter, n_iters
        real(dp) :: rL, q_scale, N_scale
        logical  :: converged, q_ok, N_ok

        rL        = par%latent_heat_water * par%water_density
        converged = .FALSE.
        n_iters   = par%max_coupling_iters

        do iter = 1, par%max_coupling_iters

            wk%q_prev = q
            wk%N_prev = N

            call update_tau_b(wk%tau_b, N, uxy_b, par)

            !$omp parallel do default(shared) private(i,j) schedule(static)
            do j = 1, wk%ny
                do i = 1, wk%nx
                    wk%mdot_total(i,j) = mdot(i,j) + wk%tau_b(i,j) * uxy_b(i,j) / rL
                end do
            end do
            !$omp end parallel do

            if (par%dissipation_melt) then
                !$omp parallel do default(shared) private(i,j) schedule(static)
                do j = 1, wk%ny
                    do i = 1, wk%nx
                        wk%mdot_total(i,j) = wk%mdot_total(i,j) &
                            + abs(q(i,j) * wk%abs_g(i,j)) / rL
                    end do
                end do
                !$omp end parallel do
            end if

            call update_psi_out(wk, mask, dx, dy, par)
            call set_q_from_psi_out(q, wk, par)

            call update_N(N, q, wk, uxy_b, A_glen, kappa, H_ice, par)

            q_scale = max(masked_max_abs(q, mask), 1.0e-15_dp)
            N_scale = max(masked_max_abs(N, mask), 1.0e-15_dp)
            q_ok = masked_max_abs_diff(q, wk%q_prev, mask) <= par%coupling_rtol * q_scale
            N_ok = masked_max_abs_diff(N, wk%N_prev, mask) <= par%coupling_rtol * N_scale

            if (q_ok .and. N_ok) then
                converged = .TRUE.
                n_iters   = iter
                exit
            end if

        end do

        if (par%coupling_verbose) then
            if (converged) then
                write(*,'(a,i0,a)') " k24: (q,N) coupling Picard loop converged after ", n_iters, " iteration(s)."
            else
                write(*,'(a,i0,a)') " k24: (q,N) coupling Picard loop did NOT converge (hit max_coupling_iters = ", n_iters, ")."
            end if
        end if

    end subroutine resolve_q_coupled

    ! ============================================================
    ! Flow routing: dispatcher + three implementations
    ! ============================================================
    ! All three compute the same field -- psi_out, the accumulated upstream
    ! water-potential outflow per cell [m3/s]. Mirrors route_psi_out!
    ! (water_flux.jl) dispatching on AbstractPsiOutAlgorithm.
    !
    !   * RECURSIVE (default) -- depth-first with memoization. Fastest, but
    !     recurses as deep as the longest flow chain, which on a real ice-sheet
    !     grid can exhaust the process stack.
    !   * ITERATIVE -- the same traversal with an explicit stack; no recursion
    !     depth limit, reproduces RECURSIVE cell for cell (max_psi_out_calls
    !     cutoff included).
    !   * TOPOSORT -- Kahn's algorithm over the flow-direction graph, a genuine
    !     single pass. Exact only if that graph is acyclic, which real
    !     topography usually is not (confirmed on Thwaites-2km at every
    !     longcoupwater). Errors on a detected cycle unless
    !     k24_toposort_allow_cycles is set.
    subroutine update_psi_out(wk, mask, dx, dy, par)
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par

        select case (par%flux_solver)
            case (K24_FLUX_RECURSIVE)
                call update_psi_out_recursive(wk, mask, dx, dy, par)
            case (K24_FLUX_ITERATIVE)
                call update_psi_out_iterative(wk, mask, dx, dy, par)
            case (K24_FLUX_TOPOSORT)
                call update_psi_out_toposort(wk, mask, dx, dy, par)
            case default
                write(*,*) "update_psi_out:: error: k24_flux_solver must be one of [0,1,2]."
                write(*,*) "flux_solver = ", par%flux_solver
                stop
        end select
    end subroutine update_psi_out

    ! --- Recursive (DFS + memoization) --------------------------
    subroutine update_psi_out_recursive(wk, mask, dx, dy, par)
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par

        integer  :: i, j, call_count
        logical  :: warned
        real(dp) :: dummy

        ! DEVIATION: psi_out is zeroed everywhere. FastHydrology.jl's persists
        ! between solves, so its non-grounded cells keep stale values; nothing
        ! downstream reads them on either side.
        wk%psi_out = 0.0_dp
        wk%visited = 0
        call_count = 0
        warned     = .FALSE.

        do j = 1, wk%ny
            do i = 1, wk%nx
                if (mask(i,j) == 1.0_dp) then
                    dummy = accumulate_psi_out(wk, i, j, mask, dx, dy, par, call_count, warned)
                end if
            end do
        end do
    end subroutine update_psi_out_recursive

    recursive function accumulate_psi_out(wk, i, j, mask, dx, dy, par, call_count, warned) result(psi_value)
        ! Mirrors accumulate_psi_out! (water_flux.jl). `call_count` caps the
        ! number of cells one sweep may visit (KORI-ULB's funcnt <= 5e4 in
        ! DpareaWarGds.m); once tripped the current cell is treated as a
        ! terminal source and returned WITHOUT the trailing max(0, .) clamp,
        ! exactly as the reference does.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        integer,              intent(IN)    :: i, j
        real(dp),             intent(IN)    :: mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par
        integer,              intent(INOUT) :: call_count
        logical,              intent(INOUT) :: warned
        real(dp) :: psi_value

        integer  :: ni, nj, d
        real(dp) :: w

        if (mask(i,j) /= 1.0_dp) then
            psi_value = 0.0_dp
            return
        end if

        if (wk%visited(i,j) == 1) then
            psi_value = wk%psi_out(i,j)
            return
        end if

        wk%visited(i,j) = 1
        wk%psi_out(i,j) = wk%mdot_total(i,j) * dx * dy

        call_count = call_count + 1
        if (call_count > par%max_psi_out_calls) then
            if (.not. warned) then
                write(*,'(a,i0,a)') " k24: WARNING accumulate_psi_out hit k24_max_psi_out_calls = ", &
                    par%max_psi_out_calls, " cells in one sweep; cutting the flow routing off early." // &
                    " Raise k24_max_psi_out_calls if this grid genuinely has more grounded cells."
                warned = .TRUE.
            end if
            psi_value = wk%psi_out(i,j)
            return
        end if

        do d = 1, 4
            ni = i + K24_DIRS(1,d)
            nj = j + K24_DIRS(2,d)
            if (ni < 1 .or. ni > wk%nx .or. nj < 1 .or. nj > wk%ny) cycle

            w = -(wk%gsx(ni,nj) * real(K24_DIRS(1,d), dp) + &
                  wk%gsy(ni,nj) * real(K24_DIRS(2,d), dp)) / (wk%abs_gs(ni,nj) + 1.0e-15_dp)

            if (w > 0.0_dp) then
                wk%psi_out(i,j) = wk%psi_out(i,j) + &
                    accumulate_psi_out(wk, ni, nj, mask, dx, dy, par, call_count, warned) * w
            end if
        end do

        ! If mdot is negative enough that all the flux refreezes, floor at zero.
        wk%psi_out(i,j) = max(0.0_dp, wk%psi_out(i,j))
        psi_value = wk%psi_out(i,j)

    end function accumulate_psi_out

    ! --- Iterative (explicit stack, same traversal) --------------
    subroutine update_psi_out_iterative(wk, mask, dx, dy, par)
        ! Mirrors update_psi_out_iterative! (water_flux.jl). Stack entry
        ! (i, j, k): k == 0 means unvisited; 1 <= k <= 4 means neighbours
        ! 1..k-1 are folded in and k is next; k == 5 means finalize. A
        ! not-yet-visited neighbour is pushed WITHOUT advancing k, so the
        ! parent frame is re-entered and takes the "already visited" branch --
        ! the explicit form of a return value flowing back to a paused caller.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par

        integer :: i, j, si, sj, sk, ni, nj, d, top, capacity, call_count
        real(dp) :: w
        logical  :: warned

        wk%psi_out = 0.0_dp
        wk%visited = 0
        call_count = 0
        warned     = .FALSE.

        ! Every cell is pushed at most once (a push is guarded by
        ! visited /= 1, and the cell is marked visited the moment it is first
        ! examined at the top of the loop), so the grounded-cell count plus one
        ! is a hard bound on the stack depth.
        capacity = count(mask == 1.0_dp) + 1
        allocate(wk%stack_i(capacity), wk%stack_j(capacity), wk%stack_k(capacity))
        top = 0

        do j = 1, wk%ny
            do i = 1, wk%nx

                if (mask(i,j) /= 1.0_dp) cycle
                if (wk%visited(i,j) == 1) cycle

                top = 1
                wk%stack_i(1) = i; wk%stack_j(1) = j; wk%stack_k(1) = 0

                do while (top > 0)

                    si = wk%stack_i(top); sj = wk%stack_j(top); sk = wk%stack_k(top)

                    if (sk == 0) then

                        wk%visited(si,sj) = 1
                        wk%psi_out(si,sj) = wk%mdot_total(si,sj) * dx * dy

                        call_count = call_count + 1
                        if (call_count > par%max_psi_out_calls) then
                            if (.not. warned) then
                                write(*,'(a,i0,a)') " k24: WARNING update_psi_out_iterative hit k24_max_psi_out_calls = ", &
                                    par%max_psi_out_calls, " cells in one sweep; cutting the flow routing off early."
                                warned = .TRUE.
                            end if
                            ! Pop WITHOUT clamping, matching the recursive
                            ! form's cap-trip branch, which returns before its
                            ! trailing max(0, .). A cut-off cell with a
                            ! negative local source therefore stays negative.
                            top = top - 1
                        else
                            wk%stack_k(top) = 1
                        end if

                    else if (sk <= 4) then

                        d  = sk
                        ni = si + K24_DIRS(1,d)
                        nj = sj + K24_DIRS(2,d)

                        if (ni < 1 .or. ni > wk%nx .or. nj < 1 .or. nj > wk%ny) then
                            wk%stack_k(top) = sk + 1
                            cycle
                        end if

                        w = -(wk%gsx(ni,nj) * real(K24_DIRS(1,d), dp) + &
                              wk%gsy(ni,nj) * real(K24_DIRS(2,d), dp)) / (wk%abs_gs(ni,nj) + 1.0e-15_dp)

                        if (w <= 0.0_dp) then
                            wk%stack_k(top) = sk + 1
                            cycle
                        end if

                        ! A non-grounded neighbour contributes 0, matching the
                        ! recursive form's mask check (which only runs once
                        ! w > 0 has already sent us into the call).
                        if (mask(ni,nj) /= 1.0_dp) then
                            wk%stack_k(top) = sk + 1
                            cycle
                        end if

                        if (wk%visited(ni,nj) == 1) then
                            wk%psi_out(si,sj) = wk%psi_out(si,sj) + wk%psi_out(ni,nj) * w
                            wk%stack_k(top) = sk + 1
                        else
                            if (top >= capacity) then
                                write(*,*) "update_psi_out_iterative:: error: routing stack overflow."
                                stop
                            end if
                            top = top + 1
                            wk%stack_i(top) = ni; wk%stack_j(top) = nj; wk%stack_k(top) = 0
                        end if

                    else

                        wk%psi_out(si,sj) = max(0.0_dp, wk%psi_out(si,sj))
                        top = top - 1

                    end if

                end do

            end do
        end do

        deallocate(wk%stack_i, wk%stack_j, wk%stack_k)

    end subroutine update_psi_out_iterative

    ! --- Topological sort (Kahn's algorithm, single pass) --------
    subroutine update_psi_out_toposort(wk, mask, dx, dy, par)
        ! Mirrors update_psi_out_topological! (water_flux.jl). The edge test is
        ! evaluated from the SOURCE cell's own gradient -- w = (gs(A).d)/|gs(A)|
        ! for the edge A -> A+d -- which is algebraically the same test the
        ! recursive form applies from the receiving side, just rearranged so a
        ! cell's outgoing edges can be built from its own data in one pass.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par

        integer :: i, j, ni, nj, d, head, tail, capacity
        integer :: total_masked, processed, n_stuck
        real(dp) :: w

        wk%psi_out = 0.0_dp   ! accumulated via += below, so must start clean

        allocate(wk%in_degree(wk%nx, wk%ny))
        wk%in_degree = 0

        total_masked = count(mask == 1.0_dp)

        ! In-degree: how many grounded neighbours flow into each grounded cell.
        !$omp parallel do default(shared) private(i,j,d,ni,nj,w) schedule(static)
        do j = 1, wk%ny
            do i = 1, wk%nx
                if (mask(i,j) /= 1.0_dp) cycle
                do d = 1, 4
                    ni = i + K24_DIRS(1,d)
                    nj = j + K24_DIRS(2,d)
                    if (ni < 1 .or. ni > wk%nx .or. nj < 1 .or. nj > wk%ny) cycle
                    if (mask(ni,nj) /= 1.0_dp) cycle

                    ! Edge (ni,nj) -> (i,j) exists iff (ni,nj) flows toward us.
                    w = (wk%gsx(ni,nj) * real(-K24_DIRS(1,d), dp) + &
                         wk%gsy(ni,nj) * real(-K24_DIRS(2,d), dp)) / (wk%abs_gs(ni,nj) + 1.0e-15_dp)

                    if (w > 0.0_dp) wk%in_degree(i,j) = wk%in_degree(i,j) + 1
                end do
            end do
        end do
        !$omp end parallel do

        capacity = max(total_masked, 1)
        allocate(wk%stack_i(capacity), wk%stack_j(capacity))
        head = 1
        tail = 0

        do j = 1, wk%ny
            do i = 1, wk%nx
                if (mask(i,j) == 1.0_dp .and. wk%in_degree(i,j) == 0) then
                    tail = tail + 1
                    wk%stack_i(tail) = i
                    wk%stack_j(tail) = j
                end if
            end do
        end do

        processed = 0
        do while (head <= tail)

            i = wk%stack_i(head)
            j = wk%stack_j(head)
            head = head + 1
            processed = processed + 1

            ! Every upstream contribution is already folded in by construction
            ! (this cell only reached the queue once all of them were final),
            ! so all that is left is its own source term and the clamp.
            wk%psi_out(i,j) = max(0.0_dp, wk%psi_out(i,j) + wk%mdot_total(i,j) * dx * dy)

            do d = 1, 4
                ni = i + K24_DIRS(1,d)
                nj = j + K24_DIRS(2,d)
                if (ni < 1 .or. ni > wk%nx .or. nj < 1 .or. nj > wk%ny) cycle
                if (mask(ni,nj) /= 1.0_dp) cycle

                w = (wk%gsx(i,j) * real(K24_DIRS(1,d), dp) + &
                     wk%gsy(i,j) * real(K24_DIRS(2,d), dp)) / (wk%abs_gs(i,j) + 1.0e-15_dp)

                if (w > 0.0_dp) then
                    wk%psi_out(ni,nj) = wk%psi_out(ni,nj) + wk%psi_out(i,j) * w
                    wk%in_degree(ni,nj) = wk%in_degree(ni,nj) - 1
                    if (wk%in_degree(ni,nj) == 0) then
                        tail = tail + 1
                        wk%stack_i(tail) = ni
                        wk%stack_j(tail) = nj
                    end if
                end if
            end do

        end do

        if (processed < total_masked) then
            n_stuck = total_masked - processed
            if (par%toposort_allow_cycles) then
                write(*,'(a,i0,a)') " k24: WARNING update_psi_out_toposort left ", n_stuck, &
                    " grounded cell(s) unprocessed (a cycle in the flow-direction graph);" // &
                    " they keep only their partial upstream contributions."
            else
                write(*,'(a,i0,a,i0,a)') " update_psi_out_toposort:: error: cycle in the flow-direction graph -- ", &
                    n_stuck, " of ", total_masked, " grounded cell(s) never reached in-degree zero."
                write(*,*) "Expect this on real topography at any k24_long_coupling_water."
                write(*,*) "Use k24_flux_solver = 0 (recursive) or 1 (iterative), or set k24_toposort_allow_cycles = .TRUE."
                stop
            end if
        end if

        deallocate(wk%in_degree, wk%stack_i, wk%stack_j)

    end subroutine update_psi_out_toposort

    ! ============================================================
    ! Effective pressure
    ! ============================================================
    subroutine update_N(N, q, wk, uxy_b, A_glen, kappa, H_ice, par)
        ! Mirrors update_N! (effective_pressure.jl): Q, S_inf, H, Po, N_inf,
        ! then the complementary-error-function transition to N.
        implicit none
        real(dp),             intent(INOUT) :: N(:,:)
        real(dp),             intent(IN)    :: q(:,:)
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: uxy_b(:,:), A_glen(:,:), kappa(:,:), H_ice(:,:)
        type(k24_param_class),intent(IN)    :: par

        integer  :: i, j, nx, ny
        real(dp) :: alpha, beta, n_glen, K, sqrt_pi
        real(dp) :: Q_c, expo, ratio, numer, denom_const, arg
        real(dp) :: sliding_coeff, melt_coeff

        nx = wk%nx; ny = wk%ny
        alpha  = par%manning_coefficient_exponent
        beta   = par%bed_friction_exponent
        n_glen = par%manning_exponent
        K      = par%K

        ! ---- Q: volumetric flux per conduit ----
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
            do i = 1, nx
                wk%Q(i,j) = q(i,j) * par%coupling_length
            end do
        end do
        !$omp end parallel do

        ! ---- S_inf: far-field conduit cross-section ----
        ! DEVIATION: Julia evaluates the formula everywhere and then overwrites
        ! Q == 0 cells with 0 (its own comment notes that those cells otherwise
        ! give 0^negative * 0^positive = Inf*0 = NaN). The branch here takes
        ! the same limit without relying on Inf/NaN surviving -Ofast.
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
            do i = 1, nx
                if (wk%Q(i,j) == 0.0_dp) then
                    wk%S_inf(i,j) = 0.0_dp
                else
                    wk%S_inf(i,j) = K**(-1.0_dp/alpha) &
                                  * wk%abs_g(i,j)**((1.0_dp - beta)/alpha) &
                                  * wk%Q(i,j)**(1.0_dp/alpha)
                end if
            end do
        end do
        !$omp end parallel do

        ! ---- H: conduit thickness ----
        ! Q_c is chosen by the drainage mode (AbstractDrainageMode): the
        ! configured value for BOTH, the exp(-Q/Q_c) -> 0 limit for EFFICIENT,
        ! the exp(-Q/Q_c) -> 1 limit for INEFFICIENT. Taking the two limits
        ! analytically avoids dividing by zero / by huge().
        Q_c = par%critical_discharge

        !$omp parallel do default(shared) private(i,j,expo) schedule(static)
        do j = 1, ny
            do i = 1, nx

                wk%H_hard(i,j) = sqrt(wk%S_inf(i,j))

                select case (par%drainage_mode)

                    case (K24_DRAINAGE_INEFFICIENT)
                        expo = 1.0_dp

                    case (K24_DRAINAGE_EFFICIENT)
                        ! Q_c -> 0: exp(-Q/Q_c) -> 0 for Q > 0. At Q == 0 the
                        ! ratio is 0/0; Julia special-cases it to H_soft = 0,
                        ! which is also the sqrt(S_inf)/F_till limit there
                        ! since S_inf == 0 whenever Q == 0.
                        if (wk%Q(i,j) == 0.0_dp) then
                            wk%H_soft(i,j)  = 0.0_dp
                            wk%H_cond(i,j)  = (1.0_dp - kappa(i,j)) * wk%H_hard(i,j)
                            cycle
                        end if
                        expo = 0.0_dp

                    case default   ! K24_DRAINAGE_BOTH
                        if (Q_c == 0.0_dp) then
                            if (wk%Q(i,j) == 0.0_dp) then
                                wk%H_soft(i,j) = 0.0_dp
                                wk%H_cond(i,j) = (1.0_dp - kappa(i,j)) * wk%H_hard(i,j)
                                cycle
                            end if
                            expo = 0.0_dp
                        else
                            expo = exp(-wk%Q(i,j) / Q_c)
                        end if

                end select

                wk%H_soft(i,j) = max(0.0_dp, par%initial_cavity_height &
                    + (wk%H_hard(i,j) / par%till_factor - par%initial_cavity_height) * expo)

                wk%H_cond(i,j) = (1.0_dp - kappa(i,j)) * wk%H_hard(i,j) &
                               + kappa(i,j) * wk%H_soft(i,j)

            end do
        end do
        !$omp end parallel do

        ! ---- Po: ice overburden pressure ----
        ! From the RAW ice thickness, not the potential-filled h: update_Po!
        ! reads state.h, which potential_filling! deliberately leaves alone
        ! (it writes its correction into model.h instead).
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
            do i = 1, nx
                wk%Po(i,j) = max(par%ice_density * par%gravity * H_ice(i,j), 1.0e5_dp)
            end do
        end do
        !$omp end parallel do

        ! ---- N_inf: far-field effective pressure ----
        denom_const = 2.0_dp * n_glen**(-n_glen) * par%ice_density * par%latent_heat_water

        ! The drainage mode gates Eq. (5b)'s two opening terms as well as
        ! update_H's Q_c (opening_coefficients in effective_pressure.jl):
        ! EFFICIENT drops the sliding-over-obstacles term, INEFFICIENT drops
        ! the melt-driven one, BOTH keeps them both.
        select case (par%drainage_mode)
            case (K24_DRAINAGE_EFFICIENT)
                sliding_coeff = 0.0_dp
                melt_coeff    = 1.0_dp
            case (K24_DRAINAGE_INEFFICIENT)
                sliding_coeff = 1.0_dp
                melt_coeff    = 0.0_dp
            case default   ! K24_DRAINAGE_BOTH
                sliding_coeff = 1.0_dp
                melt_coeff    = 1.0_dp
        end select

        !$omp parallel do default(shared) private(i,j,ratio,numer) schedule(static)
        do j = 1, ny
            do i = 1, nx
                ! DEVIATION: branch instead of computing H^2/S_inf^2 with
                ! S_inf == 0 and overwriting afterwards, as Julia does.
                if (wk%S_inf(i,j) == 0.0_dp) then
                    wk%N_inf(i,j) = wk%Po(i,j)
                    cycle
                end if

                ratio = wk%H_cond(i,j) / wk%S_inf(i,j)
                numer = sliding_coeff * par%ice_density * par%latent_heat_water &
                                      * uxy_b(i,j) * par%bed_thickness &
                      + melt_coeff * wk%Q(i,j) * wk%abs_g(i,j)

                wk%N_inf(i,j) = min(max( &
                    (ratio * ratio * numer / (denom_const * A_glen(i,j)))**(1.0_dp/n_glen), &
                    par%min_pressure_fraction * wk%Po(i,j)), wk%Po(i,j))
            end do
        end do
        !$omp end parallel do

        ! ---- N ----
        sqrt_pi = sqrt(K24_PI)
        !$omp parallel do default(shared) private(i,j,arg) schedule(static)
        do j = 1, ny
            do i = 1, nx
                ! DEVIATION: N_inf == 0 (possible once min_pressure_fraction is
                ! 0, the new default) makes the erf argument diverge. Julia
                ! evaluates erf(+-Inf)*0 = 0 there, and 0/0 -> NaN in the
                ! doubly-degenerate phi0 == 0 case; both give 0 here.
                if (wk%N_inf(i,j) > 0.0_dp) then
                    arg    = sqrt_pi * wk%phi0(i,j) / (2.0_dp * wk%N_inf(i,j))
                    N(i,j) = max(0.0_dp, erf(arg) * wk%N_inf(i,j))
                else
                    N(i,j) = 0.0_dp
                end if
            end do
        end do
        !$omp end parallel do

    end subroutine update_N

    ! ============================================================
    ! Water-layer thickness closures
    ! ============================================================
    subroutine update_W(W, q, wk, mask, par)
        ! Mirrors update_W! (water_flux.jl), dispatching on
        ! water_thickness_algorithm. Must run after update_N: the areal-conduit
        ! closure reads S_inf.
        implicit none
        real(dp),             intent(OUT) :: W(:,:)
        real(dp),             intent(IN)  :: q(:,:), mask(:,:)
        type(k24_work_class), intent(IN)  :: wk
        type(k24_param_class),intent(IN)  :: par

        integer  :: i, j, nx, ny
        real(dp) :: g_mean, num, third

        nx = wk%nx; ny = wk%ny
        third = 1.0_dp / 3.0_dp

        select case (par%water_thickness_algorithm)

            case (K24_WTHICK_AREAL_CONDUIT)
                ! S_inf / l_c: the conduit cross-section smeared over the
                ! inter-conduit spacing. Deliberately unclamped -- the
                ! [W_min, W_max] bounds were chosen for a thin sheet.
                !$omp parallel do default(shared) private(i,j) schedule(static)
                do j = 1, ny
                    do i = 1, nx
                        W(i,j) = wk%S_inf(i,j) / par%coupling_length
                    end do
                end do
                !$omp end parallel do

            case (K24_WTHICK_DARCY_WEISBACH)
                ! Turbulent parallel-plate inversion for a wide slot,
                ! d = (f*rho_w*q^2 / (4*|grad phi0|))^(1/3), using the
                ! UNSMOOTHED gradient. f is the same Darcy-Weisbach friction
                ! factor that K folds into S_inf -- one parameter, not two.
                if (par%gradient_convention == K24_GRAD_MEAN) then
                    g_mean = masked_mean(wk%abs_g, mask)
                    !$omp parallel do default(shared) private(i,j,num) schedule(static)
                    do j = 1, ny
                        do i = 1, nx
                            num = par%friction_factor * par%water_density * q(i,j) * q(i,j) &
                                / (4.0_dp * g_mean + 1.0e-15_dp)
                            W(i,j) = min(par%W_max, max(par%W_min, num**third))
                        end do
                    end do
                    !$omp end parallel do
                else
                    !$omp parallel do default(shared) private(i,j,num) schedule(static)
                    do j = 1, ny
                        do i = 1, nx
                            num = par%friction_factor * par%water_density * q(i,j) * q(i,j) &
                                / (4.0_dp * wk%abs_g(i,j) + 1.0e-15_dp)
                            W(i,j) = min(par%W_max, max(par%W_min, num**third))
                        end do
                    end do
                    !$omp end parallel do
                end if

            case (K24_WTHICK_LAMINAR)
                ! Le Brocq / Weertman laminar inversion (Eq. 8, Kazmierczak et
                ! al 2022), d = (12*eta_w*q / |grad phi0_s|)^(1/3), using the
                ! SMOOTHED gradient.
                if (par%gradient_convention == K24_GRAD_MEAN) then
                    ! Matches Kori-ULB's mean(gdsmag(...)) in SubWaterFlux.m.
                    ! Julia adds no 1e-15 guard in this one branch (unlike the
                    ! local variant and both Darcy-Weisbach ones); the max()
                    ! below keeps a zero-mean domain from producing a negative
                    ! or NaN W rather than reproducing an infinity.
                    g_mean = masked_mean(wk%abs_gs, mask)
                    if (g_mean <= 0.0_dp) then
                        W = min(par%W_max, max(par%W_min, 0.0_dp))
                    else
                        !$omp parallel do default(shared) private(i,j,num) schedule(static)
                        do j = 1, ny
                            do i = 1, nx
                                num = 12.0_dp * par%eta_w * q(i,j) / g_mean
                                W(i,j) = min(par%W_max, max(par%W_min, num**third))
                            end do
                        end do
                        !$omp end parallel do
                    end if
                else
                    !$omp parallel do default(shared) private(i,j,num) schedule(static)
                    do j = 1, ny
                        do i = 1, nx
                            num = 12.0_dp * par%eta_w * q(i,j) / (wk%abs_gs(i,j) + 1.0e-15_dp)
                            W(i,j) = min(par%W_max, max(par%W_min, num**third))
                        end do
                    end do
                    !$omp end parallel do
                end if

            case default
                write(*,*) "update_W:: error: k24_water_thickness_algorithm must be one of [0,1,2]."
                write(*,*) "water_thickness_algorithm = ", par%water_thickness_algorithm
                stop

        end select

    end subroutine update_W

end module fast_hydrology_k24
