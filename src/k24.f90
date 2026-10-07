module fast_hydrology_k24
    ! K24 effective-pressure / water-flux model (Kazmierczak et al 2024,
    ! https://doi.org/10.5194/tc-18-5887-2024).
    !
    ! Diagnostic only: reads ice geometry, bed, the terms of the basal melt
    ! rate, sliding speed and Glen rate factor; returns the distributed water
    ! flux magnitude q (and its components q_x, q_y), the effective pressure N,
    ! the water pressure p_w = Po - N, the reportable water-layer thickness W
    ! and the frictional and dissipation heat it used. Does not touch the till
    ! storage W_til, which is owned by the bucket model.
    !
    ! ============================================================
    ! Reference implementation
    ! ============================================================
    ! This module is a line-by-line port of the `kazmierczak2024` model in
    ! TakisAngelides/FastHydrology.jl (src/models/kazmierczak2024/{model,
    ! water_flux,routing,effective_pressure,sliding_law}.jl). Where the two
    ! could differ, the Julia side is the source of truth. The call sequence
    ! mirrors `update_steady_state!` in run.jl:
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
    ! Water source: built from its terms
    ! ============================================================
    ! The melt rate is never supplied whole. calc_k24 takes the geothermal
    ! heat G and the conductive heat into the ice q_T [W/m2], and the water
    ! reaching the bed from above i_eb (drained englacial water, surface
    ! input; Sommers et al 2018's i_eb) [kg/m2/s], and computes the frictional
    ! heat Q_b and the dissipation heat Q_diss [W/m2] itself:
    !
    !     mdot_fixed = (G - q_T) / L_w
    !     mdot_total = mdot_fixed + (Q_b + Q_diss) / L_w + i_eb     [kg/m2/s]
    !
    ! mdot_fixed + (Q_b + Q_diss)/L_w is the basal melt rate; i_eb is routed
    ! with it but is not melt. Q_b and Q_diss depend on the water flux (through
    ! N and q), so they are recomputed every Picard sweep. The routed volume
    ! flux is seeded with mdot_total*dx*dy/rho_w [m3/s], as in Julia.
    !
    ! ============================================================
    ! Deliberate deviations from FastHydrology.jl
    ! ============================================================
    ! Each is marked with a "DEVIATION:" comment at its site. In summary:
    !   1. Degenerate 0/0 and x/0 forms that Julia leaves to IEEE arithmetic
    !      (S_inf with Q == 0, N_inf with S_inf == 0, N with N_inf == 0) are
    !      resolved by an explicit branch here, taking the same limit Julia's
    !      own cleanups take. This library is built with -Ofast
    !      (-fp-model fast=2), under which Inf/NaN propagation is not reliable.
    !   2. `masked_mean` over an empty mask returns 0 here rather than NaN.
    !   3. psi_out is zeroed over the whole grid at the start of every sweep.
    !      Julia's persists across calls, so its non-grounded cells carry stale
    !      values; nothing downstream reads them on either side.
    !   4. Periodic domains (par%periodic_x/periodic_y, set by hydro_init) have
    !      no Julia counterpart: every neighbour stencil (potential filling,
    !      gradients, smoothing padding, routing weights, flow routing, face
    !      fluxes, staggered friction) wraps in a periodic direction, and the
    !      priority-flood fill has no outlet on a periodic edge. With both
    !      .FALSE. (the default) nothing changes.
    !   5. The workspace is call-scoped, so the face-assembled dissipation that
    !      Julia keeps in model.routing_tape between solves is passed in and out
    !      through the optional `diss_io` argument instead (hydro_update keeps
    !      it in hyd%now%diss_face).

    use nml
    use phys_constants, only : sec_year_julian

    implicit none

    integer, parameter :: dp = kind(1.d0)

    real(dp), parameter :: K24_PI = 3.14159265358979323846_dp

    ! Seconds per year used by the parameter defaults FastHydrology.jl writes
    ! as perYear2perSecond(...) -- eta_w, the sliding laws' u0 and the
    ! staggered-friction velocity floor. Julia's SECONDS_PER_YEAR is the Julian
    ! year, 60^2*24*365.25. Deliberately NOT par%sec_year, the host's calendar
    ! year, which converts the public API's time argument.
    real(dp), parameter :: K24_SEC_PER_YEAR = sec_year_julian

    ! ---------- Substrate-type enum (par%substrate_type) ----------
    integer, parameter, public :: K24_SUBSTRATE_HARD  = 0
    integer, parameter, public :: K24_SUBSTRATE_SOFT  = 1
    integer, parameter, public :: K24_SUBSTRATE_MIXED = 2

    ! ---------- Flow-routing algorithm enum (par%flux_solver) ----------
    ! Mirrors AbstractPsiOutAlgorithm in FastHydrology.jl/model.jl.
    ! RECURSIVE/ITERATIVE/TOPOSORT implement the GDS_WARNER routing only;
    ! every other routing scheme needs TAPED (the default).
    integer, parameter, public :: K24_FLUX_RECURSIVE = 0   ! RecursivePsiOut
    integer, parameter, public :: K24_FLUX_ITERATIVE = 1   ! IterativePsiOut
    integer, parameter, public :: K24_FLUX_TOPOSORT  = 2   ! TopologicalPsiOut
    integer, parameter, public :: K24_FLUX_TAPED     = 3   ! TapedPsiOut (default)

    ! ---------- Routing-scheme enum (par%routing_scheme) ----------
    ! Mirrors AbstractRoutingScheme: the schemes compared by Le Brocq, Payne &
    ! Siegert (2006). GDS_WARNER is the original K24/KORI routing.
    integer, parameter, public :: K24_ROUTE_WARNER            = 0  ! Warner (default)
    integer, parameter, public :: K24_ROUTE_GDS_WARNER        = 1  ! GDSWarner
    integer, parameter, public :: K24_ROUTE_QUINN             = 2  ! Quinn (k24_quinn_original)
    integer, parameter, public :: K24_ROUTE_TARBOTON          = 3  ! Tarboton
    integer, parameter, public :: K24_ROUTE_MODIFIED_TARBOTON = 4  ! ModifiedTarboton
    integer, parameter, public :: K24_ROUTE_GDS_TARBOTON      = 5  ! GDSTarboton

    ! ---------- Fill-algorithm enum (par%fill_algorithm) ----------
    ! Mirrors AbstractFillAlgorithm. AUTO follows the routing scheme: JACOBI
    ! for the GDS schemes, PRIORITY_FLOOD otherwise.
    integer, parameter, public :: K24_FILL_AUTO             = -1
    integer, parameter, public :: K24_FILL_JACOBI           =  0  ! JacobiFill
    integer, parameter, public :: K24_FILL_LOWEST_NEIGHBOUR =  1  ! LowestNeighbourFill
    integer, parameter, public :: K24_FILL_PRIORITY_FLOOD   =  2  ! PriorityFloodFill (k24_priority_flood_epsilon)

    ! ---------- q-conversion enum (par%q_conversion) ----------
    ! Mirrors AbstractQConversion. AUTO: FACE_AVERAGE for WARNER, else OUTFLOW.
    integer, parameter, public :: K24_QCONV_AUTO         = -1
    integer, parameter, public :: K24_QCONV_OUTFLOW      =  0  ! QFromOutflow
    integer, parameter, public :: K24_QCONV_FACE_AVERAGE =  1  ! QFromFaceAverage

    ! ---------- Dissipation-discretization enum (par%dissipation_discretization) ----------
    ! Mirrors AbstractDissipationDiscretization. AUTO: FACE for WARNER, else CELL.
    integer, parameter, public :: K24_DISS_AUTO = -1
    integer, parameter, public :: K24_DISS_CELL =  0  ! CellCentredDissipation
    integer, parameter, public :: K24_DISS_FACE =  1  ! FaceDissipation

    ! ---------- Friction-discretization enum (par%friction_discretization) ----------
    ! Mirrors AbstractFrictionDiscretization. STAGGERED needs the C-grid
    ! velocities ux_b/uy_b (acx/acy) passed to calc_k24.
    integer, parameter, public :: K24_FRICTION_CELL      = 0  ! CellCentredFriction (default)
    integer, parameter, public :: K24_FRICTION_STAGGERED = 1  ! StaggeredFriction (k24_friction_quadrature)

    ! ---------- Drainage-mode enum (par%drainage_mode) ----------
    ! Mirrors AbstractDrainageMode. Note the Q_c limits below are the
    ! *corrected* version of the paper's Sect. 3.2 description (as published it
    ! swaps them; the authors confirmed the typo by email).
    integer, parameter, public :: K24_DRAINAGE_BOTH        = 0  ! BothDrainage
    integer, parameter, public :: K24_DRAINAGE_EFFICIENT   = 1  ! EfficientOnly   (Q_c -> 0)
    integer, parameter, public :: K24_DRAINAGE_INEFFICIENT = 2  ! InefficientOnly (Q_c -> Inf)

    ! ---------- Water-thickness closure enum (par%water_thickness_algorithm) ----------
    integer, parameter, public :: K24_WTHICK_DARCY_WEISBACH = 0  ! DarcyWeisbachThickness (default)
    integer, parameter, public :: K24_WTHICK_LAMINAR        = 1  ! LaminarThickness
    integer, parameter, public :: K24_WTHICK_AREAL_CONDUIT  = 2  ! ArealConduitThickness

    ! ---------- Gradient-convention enum (par%gradient_convention) ----------
    integer, parameter, public :: K24_GRAD_MEAN  = 0  ! MeanGradient (default)
    integer, parameter, public :: K24_GRAD_LOCAL = 1  ! LocalGradient

    ! ---------- Sliding-law enum (par%sliding_law) ----------
    ! Mirrors AbstractSlidingLaw. NO_FRICTION/WEERTMAN/PRESCRIBED_FIELD do not
    ! depend on N; the others do, and switch resolve_q to the joint (q, N)
    ! Picard loop. PRESCRIBED_FIELD needs tau_b_in and REG_COULOMB_FIELD needs
    ! c_till_in passed to calc_k24; SHAKTI_REG_COULOMB uses
    ! lambda = k24_shakti_lambda_coeff * A_glen.
    integer, parameter, public :: K24_SLIDING_NO_FRICTION       = 0  ! NoFrictionSlidingLaw (default, tau_b = 0)
    integer, parameter, public :: K24_SLIDING_WEERTMAN          = 1  ! WeertmanSlidingLaw
    integer, parameter, public :: K24_SLIDING_POWER_PLASTIC     = 2  ! PowerPlasticSlidingLaw
    integer, parameter, public :: K24_SLIDING_REG_COULOMB       = 3  ! RegularizedCoulombSlidingLaw
    integer, parameter, public :: K24_SLIDING_PRESCRIBED_FIELD  = 4  ! PrescribedFieldSlidingLaw
    integer, parameter, public :: K24_SLIDING_REG_COULOMB_FIELD = 5  ! RegularizedCoulombFieldSlidingLaw
    integer, parameter, public :: K24_SLIDING_SHAKTI_REG_COULOMB = 6 ! ShaktiRegularizedCoulombSlidingLaw

    ! ============================================================
    ! Parameters (runtime, namelist-overridable)
    ! ============================================================
    ! The Julia field each entry corresponds to is named in the trailing
    ! comment.
    type k24_param_class

        ! -- mode selectors --
        integer  :: substrate_type              ! (Fortran-only; builds kappa)
        integer  :: flux_solver                 ! psi_out_algorithm
        integer  :: routing_scheme              ! routing_scheme
        logical  :: quinn_original              ! Quinn(original = ...)
        integer  :: fill_algorithm              ! fill_algorithm (AUTO follows routing_scheme)
        real(dp) :: priority_flood_epsilon      ! PriorityFloodFill(epsilon = ...) [Pa]
        integer  :: q_conversion                ! q_conversion (AUTO follows routing_scheme)
        integer  :: dissipation_discretization  ! dissipation_discretization (AUTO follows routing_scheme)
        integer  :: friction_discretization     ! friction_discretization
        logical  :: friction_quadrature         ! StaggeredFriction(quadrature = ...)
        real(dp) :: friction_u_floor            ! StaggeredFriction(u_floor = ...) [m/s]
        integer  :: drainage_mode               ! drainage_mode
        integer  :: water_thickness_algorithm   ! water_thickness_algorithm
        integer  :: gradient_convention         ! DarcyWeisbach/Laminar gradient_convention
        integer  :: sliding_law                 ! sliding_law
        logical  :: toposort_allow_cycles       ! TopologicalPsiOut.allow_cycles
        logical  :: ub_hook                     ! (Fortran-only) host re-evaluates N from its sliding speed inside its velocity iteration (hydro_N_from_ub)

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
        real(dp) :: coupling_length_kamb86      ! coupling_length_kamb86: Kamb & Echelmeyer (1986) stress-gradient-
                                                 ! coupling length as a multiple of mean ice thickness (>= 0; 0 disables)
        real(dp) :: min_pressure_fraction       ! sigmat
        real(dp) :: eta_w                       ! eta_w  [Pa s]

        ! -- clamps --
        real(dp) :: W_min                       ! Wmin  [m]
        real(dp) :: W_max                       ! Wmax  [m]
        real(dp) :: q_min                       ! q_min [m2/s]
        real(dp) :: q_max                       ! q_max [m2/s]

        ! -- solver configuration --
        integer  :: fill_iters                  ! fill_iters (JACOBI / LOWEST_NEIGHBOUR)
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
        real(dp) :: reg_coulomb_q               ! RegularizedCoulomb(Field)SlidingLaw.q
        real(dp) :: reg_coulomb_u0              ! RegularizedCoulomb(Field)SlidingLaw.u0 [m/s]
        real(dp) :: shakti_C                    ! ShaktiRegularizedCoulombSlidingLaw.C
        real(dp) :: shakti_n                    ! ShaktiRegularizedCoulombSlidingLaw.n
        real(dp) :: shakti_lambda_coeff         ! lambda = shakti_lambda_coeff * A_glen

        ! -- derived, filled by k24_finalize_par --
        real(dp) :: K                           ! K = (2/pi)^(1/4)*sqrt((pi+2)/(rho_w*f))
        real(dp) :: long_coupling_water         ! longcoupwater = coupling_length_kamb86/2 (KORI-ULB's own parameter)
        integer  :: fill_alg                    ! fill_algorithm with AUTO resolved
        integer  :: q_conv                      ! q_conversion with AUTO resolved
        integer  :: diss_disc                   ! dissipation_discretization with AUTO resolved
        integer  :: n_dirs                      ! 4 or 8 neighbour directions of routing_scheme

        ! -- domain topology (Fortran-only; set by hydro_init, not the namelist) --
        logical  :: periodic_x                  ! x wraps with period nx (no halo); else edge-clamped
        logical  :: periodic_y                  ! y wraps with period ny (no halo); else edge-clamped

    end type

    ! ============================================================
    ! Scratch workspace
    ! ============================================================
    ! Mirrors KazmierczakWorkspace (and its RoutingTape). Allocated at the
    ! top of calc_k24 and released at the end.
    type k24_work_class
        integer :: nx = 0
        integer :: ny = 0
        ! geometric potential
        real(dp), allocatable :: phi0(:,:), phi0_filled(:,:), phi0_tmp(:,:), h(:,:)
        real(dp), allocatable :: gx(:,:),  gy(:,:),  abs_g(:,:)      ! gx/gy: filled potential; abs_g: true potential
        real(dp), allocatable :: gsx(:,:), gsy(:,:), abs_gs(:,:)     ! routing directions (smoothed for the GDS schemes)
        ! water flux
        real(dp), allocatable :: mdot_fixed(:,:), mdot_total(:,:), psi_out(:,:), corfac(:,:)
        real(dp), allocatable :: q_prev(:,:), N_prev(:,:), tau_b(:,:), lambda(:,:)
        real(dp), allocatable :: Q_b(:,:), Q_diss(:,:)
        integer,  allocatable :: visited(:,:)
        ! routing weights / face fluxes / face dissipation (routing.jl)
        real(dp), allocatable :: w8(:,:,:)        ! w8(d,i,j): fraction of (i,j)'s outflow sent in direction d
        real(dp), allocatable :: Fx(:,:)          ! (nx+1,ny) net x-face volume fluxes [m3/s]
        real(dp), allocatable :: Fy(:,:)          ! (nx,ny+1) net y-face volume fluxes [m3/s]
        real(dp), allocatable :: diss(:,:)        ! face-assembled dissipation melt [kg/m2/s]
        ! routing tape (TapedPsiOut)
        logical :: tape_valid = .FALSE.
        integer :: tape_n = 0
        integer,  allocatable :: t_dst_i(:), t_dst_j(:), t_src_i(:), t_src_j(:)
        real(dp), allocatable :: t_w(:)
        ! effective pressure
        real(dp), allocatable :: Q(:,:), S_inf(:,:)
        real(dp), allocatable :: H_hard(:,:), H_soft(:,:), H_cond(:,:)
        real(dp), allocatable :: N_inf(:,:), Po(:,:)
        ! routing scratch (iterative / topological solvers, tape recorder, priority flood)
        integer,  allocatable :: in_degree(:,:)
        integer,  allocatable :: stack_i(:), stack_j(:), stack_k(:)
    end type

    ! Direction offsets, ROUTE_OFFSETS in routing.jl. The first four are the
    ! order accumulate_psi_out! iterates them. Accumulation is order-dependent
    ! in floating point, so this order is part of the numerics.
    integer, parameter :: K24_DIRS(2,8) = reshape([ -1, 0,  1, 0,  0, -1,  0, 1, &
                                                     -1,-1,  1,-1, -1,  1,  1, 1 ], [2,8])
    integer, parameter :: K24_OPPOSITE(8) = [2, 1, 4, 3, 8, 7, 6, 5]
    ! Tarboton (1997) triangular facets as (cardinal, diagonal) direction pairs.
    integer, parameter :: K24_FACETS(2,8) = reshape([2,6, 2,8, 4,8, 4,7, 1,7, 1,5, 3,5, 3,6], [2,8])

    private
    public :: k24_param_class
    public :: k24_par_load
    public :: k24_finalize_par
    public :: initialize_kappa
    public :: calc_k24
    public :: k24_N_from_ub

contains

    ! ============================================================
    ! Namelist load
    ! ============================================================
    subroutine k24_par_load(par, filename, group, init, skip_phys_const)
        ! skip_phys_const: the host supplied the physical constants, so
        ! k24_latent_heat_water is not read from &yhyd -- hydro_par_load sets
        ! latent_heat_water from the host's L_ice instead, keeping the melt and
        ! opening terms consistent with the host's energy budget.

        implicit none

        type(k24_param_class), intent(INOUT) :: par
        character(len=*),      intent(IN)    :: filename
        character(len=*),      intent(IN)    :: group
        logical, optional,     intent(IN)    :: init
        logical, optional,     intent(IN)    :: skip_phys_const

        logical :: init_pars
        logical :: read_phys_const

        character(len=*), parameter :: def_file  = "input/yelmo_defaults.nml"
        character(len=*), parameter :: def_group = "yhyd"

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        read_phys_const = .TRUE.
        if (present(skip_phys_const)) read_phys_const = .not. skip_phys_const

        par%substrate_type                = K24_SUBSTRATE_HARD
        par%flux_solver                   = K24_FLUX_TAPED
        par%routing_scheme                = K24_ROUTE_WARNER
        par%quinn_original                = .FALSE.
        par%fill_algorithm                = K24_FILL_AUTO
        par%priority_flood_epsilon        = 1.0_dp
        par%q_conversion                  = K24_QCONV_AUTO
        par%dissipation_discretization    = K24_DISS_AUTO
        par%friction_discretization       = K24_FRICTION_CELL
        par%friction_quadrature           = .FALSE.
        ! u_floor = perYear2perSecond(1e-3), Yelmo's ub_sq_min
        par%friction_u_floor              = 1.0e-3_dp / K24_SEC_PER_YEAR
        par%drainage_mode                 = K24_DRAINAGE_BOTH
        par%water_thickness_algorithm     = K24_WTHICK_DARCY_WEISBACH
        par%gradient_convention           = K24_GRAD_MEAN
        par%sliding_law                   = K24_SLIDING_NO_FRICTION
        par%toposort_allow_cycles         = .FALSE.
        par%ub_hook                       = .TRUE.

        par%water_density                 = 1000.0_dp
        par%ice_density                   =  917.0_dp
        par%gravity                       =    9.81_dp
        par%manning_exponent              =    3.0_dp
        ! L_w = 3.34e5 (KAZMIERCZAK_DEFAULT_L_W).
        par%latent_heat_water             =    3.34e5_dp
        par%bed_thickness                 =    0.1_dp
        par%manning_coefficient_exponent  =    1.25_dp      ! alpha = 5/4
        par%bed_friction_exponent         =    1.5_dp       ! beta  = 3/2
        par%friction_factor               =    0.1_dp
        par%till_factor                   =    1.1_dp
        par%critical_discharge            =    1.0_dp
        par%initial_cavity_height         =    0.1_dp
        par%coupling_length               =    1.0e4_dp
        ! Upper edge of Kamb & Echelmeyer's 4-10x ice-thickness range for ice
        ! sheets (~1-3x valley glaciers, ~12x surging). 0 at grids coarser
        ! than the coupling length.
        par%coupling_length_kamb86        =   10.0_dp
        ! sigmat = 0 -- no N_inf floor. Pass 0.02 for KORI-ULB's own value.
        par%min_pressure_fraction         =    0.0_dp
        ! eta_w = perYear2perSecond(1.8e-3): KORI-ULB's par.waterviscosity is a
        ! per-year quantity.
        par%eta_w                         = 1.8e-3_dp / K24_SEC_PER_YEAR

        ! The four KORI-ULB clamps, off by default. huge() stands in for Inf.
        par%W_min                         =    0.0_dp
        par%W_max                         = huge(1.0_dp)
        par%q_min                         =    0.0_dp
        par%q_max                         = huge(1.0_dp)

        par%fill_iters                    =   10
        par%max_psi_out_calls             = 100000
        par%max_dissipation_iters         =   20
        par%dissipation_rtol              = 1.0e-12_dp
        par%dissipation_melt              = .TRUE.
        par%dissipation_verbose           = .TRUE.
        par%max_coupling_iters            =   20
        par%coupling_rtol                 = 1.0e-8_dp
        par%coupling_verbose              = .TRUE.

        ! Sliding-law parameters. C and c_till have no Julia default (they are
        ! mandatory there); 0 here makes an unconfigured law a no-op.
        par%weertman_C                    = 0.0_dp
        par%weertman_q                    = 1.0_dp / 3.0_dp
        par%power_plastic_c_till          = 0.0_dp
        par%power_plastic_q               = 1.0_dp
        par%power_plastic_u0              = 100.0_dp / K24_SEC_PER_YEAR
        par%reg_coulomb_c_till            = 0.0_dp
        par%reg_coulomb_q                 = 1.0_dp / 3.0_dp
        par%reg_coulomb_u0                = 100.0_dp / K24_SEC_PER_YEAR
        par%shakti_C                      = 0.0_dp
        par%shakti_n                      = 3.0_dp
        par%shakti_lambda_coeff           = 1.5_dp

        par%periodic_x                    = .FALSE.
        par%periodic_y                    = .FALSE.

        call nml_read(filename,group,"k24_substrate_type",                par%substrate_type,                init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_flux_solver",                   par%flux_solver,                   init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_routing_scheme",                par%routing_scheme,                init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_quinn_original",                par%quinn_original,                init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_fill_algorithm",                par%fill_algorithm,                init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_priority_flood_epsilon",        par%priority_flood_epsilon,        init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_q_conversion",                  par%q_conversion,                  init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_dissipation_discretization",    par%dissipation_discretization,    init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_friction_discretization",       par%friction_discretization,       init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_friction_quadrature",           par%friction_quadrature,           init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_friction_u_floor",              par%friction_u_floor,              init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_drainage_mode",                 par%drainage_mode,                 init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_water_thickness_algorithm",     par%water_thickness_algorithm,     init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_gradient_convention",           par%gradient_convention,           init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_sliding_law",                   par%sliding_law,                   init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_toposort_allow_cycles",         par%toposort_allow_cycles,         init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_ub_hook",                       par%ub_hook,                       init=init_pars,defaults_file=def_file,defaults_group=def_group)
        ! water_density, ice_density, gravity are set from top-level
        ! rho_w / rho_ice / g in hydro_par_load (single source of truth).
        call nml_read(filename,group,"k24_manning_exponent",              par%manning_exponent,              init=init_pars,defaults_file=def_file,defaults_group=def_group)
        if (read_phys_const) then
            call nml_read(filename,group,"k24_latent_heat_water",             par%latent_heat_water,             init=init_pars,defaults_file=def_file,defaults_group=def_group)
        end if
        call nml_read(filename,group,"k24_bed_thickness",                 par%bed_thickness,                 init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_manning_coefficient_exponent",  par%manning_coefficient_exponent,  init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_bed_friction_exponent",         par%bed_friction_exponent,         init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_friction_factor",               par%friction_factor,               init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_till_factor",                   par%till_factor,                   init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_critical_discharge",            par%critical_discharge,            init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_initial_cavity_height",         par%initial_cavity_height,         init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_coupling_length",               par%coupling_length,               init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_coupling_length_kamb86",        par%coupling_length_kamb86,        init=init_pars,defaults_file=def_file,defaults_group=def_group)
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
        call nml_read(filename,group,"k24_shakti_C",                      par%shakti_C,                      init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_shakti_n",                      par%shakti_n,                      init=init_pars,defaults_file=def_file,defaults_group=def_group)
        call nml_read(filename,group,"k24_shakti_lambda_coeff",           par%shakti_lambda_coeff,           init=init_pars,defaults_file=def_file,defaults_group=def_group)

        ! Derived. hydro_par_load overwrites water_density afterwards with the
        ! top-level par%rho_w and calls k24_finalize_par again.
        call k24_finalize_par(par)

        return

    end subroutine k24_par_load

    subroutine k24_finalize_par(par)
        ! Recompute every parameter derived from another one and check the
        ! option combinations the KazmierczakHydroModel constructor checks.
        implicit none
        type(k24_param_class), intent(INOUT) :: par

        logical :: gds

        ! Manning-Strickler / Darcy-Weisbach conductivity coefficient,
        ! K = (2/pi)^(1/4) * sqrt((pi+2)/(rho_w*f)).
        par%K = (2.0_dp / K24_PI)**0.25_dp * &
                sqrt((K24_PI + 2.0_dp) / (par%water_density * par%friction_factor))

        ! The smoothing kernel's 2D area-weighted effective coupling length is
        ! 2*longcoupwater*h_avg, so longcoupwater = coupling_length_kamb86/2
        ! makes the effective length coupling_length_kamb86*h_avg exactly.
        if (par%coupling_length_kamb86 < 0.0_dp) then
            write(*,*) "k24_finalize_par:: error: k24_coupling_length_kamb86 must be >= 0 (0 disables the smoothing)."
            write(*,*) "k24_coupling_length_kamb86 = ", par%coupling_length_kamb86
            stop
        end if
        par%long_coupling_water = par%coupling_length_kamb86 / 2.0_dp

        if (par%routing_scheme < K24_ROUTE_WARNER .or. par%routing_scheme > K24_ROUTE_GDS_TARBOTON) then
            write(*,*) "k24_finalize_par:: error: k24_routing_scheme must be one of [0,1,2,3,4,5]."
            write(*,*) "routing_scheme = ", par%routing_scheme
            stop
        end if

        gds = (par%routing_scheme == K24_ROUTE_GDS_WARNER .or. par%routing_scheme == K24_ROUTE_GDS_TARBOTON)

        if (par%routing_scheme == K24_ROUTE_GDS_WARNER .or. par%routing_scheme == K24_ROUTE_WARNER) then
            par%n_dirs = 4
        else
            par%n_dirs = 8
        end if

        ! Unset options follow the routing scheme (AbstractRoutingScheme):
        ! Warner gets its face-based partners, the GDS schemes the original
        ! K24 choices.
        par%fill_alg = par%fill_algorithm
        if (par%fill_alg == K24_FILL_AUTO) then
            if (gds) then
                par%fill_alg = K24_FILL_JACOBI
            else
                par%fill_alg = K24_FILL_PRIORITY_FLOOD
            end if
        end if
        par%q_conv = par%q_conversion
        if (par%q_conv == K24_QCONV_AUTO) then
            if (par%routing_scheme == K24_ROUTE_WARNER) then
                par%q_conv = K24_QCONV_FACE_AVERAGE
            else
                par%q_conv = K24_QCONV_OUTFLOW
            end if
        end if
        par%diss_disc = par%dissipation_discretization
        if (par%diss_disc == K24_DISS_AUTO) then
            if (par%routing_scheme == K24_ROUTE_WARNER) then
                par%diss_disc = K24_DISS_FACE
            else
                par%diss_disc = K24_DISS_CELL
            end if
        end if

        if (par%routing_scheme /= K24_ROUTE_GDS_WARNER .and. par%flux_solver /= K24_FLUX_TAPED) then
            write(*,*) "k24_finalize_par:: error: k24_routing_scheme = ", par%routing_scheme, &
                       " requires k24_flux_solver = 3 (taped); the other flux solvers only implement GDS_WARNER (1)."
            stop
        end if
        if ((par%q_conv == K24_QCONV_FACE_AVERAGE .or. par%diss_disc == K24_DISS_FACE) .and. par%n_dirs /= 4) then
            write(*,*) "k24_finalize_par:: error: face-average q (k24_q_conversion = 1) and face dissipation"
            write(*,*) "(k24_dissipation_discretization = 1) need a 4-neighbour routing scheme (0 WARNER or 1 GDS_WARNER)."
            stop
        end if

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
    ! Main entry point
    ! ============================================================
    subroutine calc_k24(q_x, q_y, N, p_w, W, q, Q_b, Q_diss, &
                        H_ice, z_bed, mask, G, q_T, i_eb, uxy_b, A_glen, kappa, &
                        dx, dy, par, ux_b, uy_b, tau_b_in, c_till_in, diss_io, C_frz, &
                        gsx_out, gsy_out, absgs_out, absg_out, phi0_out)
        ! Mirrors update_steady_state! (FastHydrology.jl/.../run.jl):
        ! update_q, then update_N, then update_W.
        !
        ! `q` and `N` are INOUT: FastHydrology.jl's model.q and state.N persist
        ! between solves, and both Picard loops warm-start from them. Pass the
        ! caller's own persistent fields, zero-initialised on the first call.
        ! `diss_io` (optional, [kg/m2/s]) plays the same role for the
        ! face-assembled dissipation (deviation 5).
        !
        ! Water source terms (see the module header): G, q_T [W/m2] and i_eb
        ! [kg/m2/s]. Q_b and Q_diss [W/m2] are returned as used in the last
        ! Picard sweep.
        !
        ! `C_frz` (optional, [m/s ice equivalent]) is the freeze-on capacity
        ! for the host's capacity basal boundary condition: the water that
        ! reaches each grounded cell, routed in from upstream (Psi_in, from the
        ! final routing) plus the water from above i_eb,
        ! (rho_w*Psi_in/(dx*dy) + i_eb)/rho_i. The cell's own melt and
        ! dissipation are not in it: they are heat, which the host counts in
        ! its freezing demand. i_eb is water and in no heat balance.
        ! Fortran-only for now.
        !
        ! Optional inputs needed by some options only:
        !   ux_b, uy_b : C-grid basal velocities [m/s] (acx: face between i and
        !                i+1; acy: face between j and j+1) -- STAGGERED friction
        !   tau_b_in   : basal shear stress magnitude [Pa] -- PRESCRIBED_FIELD
        !   c_till_in  : per-cell Coulomb coefficient -- REG_COULOMB_FIELD
        !
        ! The trailing optional outputs expose internals so tests/k24_synth.f90
        ! can compare the intermediate fields against FastHydrology.jl stage by
        ! stage. Ordinary callers omit them.

        implicit none

        real(dp), intent(OUT)   :: q_x(:,:), q_y(:,:)   ! water flux components [m2/s]
        real(dp), intent(INOUT) :: N(:,:)               ! effective pressure [Pa]
        real(dp), intent(OUT)   :: p_w(:,:)             ! water pressure (Po - N) [Pa]
        real(dp), intent(OUT)   :: W(:,:)               ! water layer thickness [m]
        real(dp), intent(INOUT) :: q(:,:)               ! distributed flux magnitude [m2/s]
        real(dp), intent(OUT)   :: Q_b(:,:)             ! frictional heat [W/m2]
        real(dp), intent(OUT)   :: Q_diss(:,:)          ! dissipation heat [W/m2]
        real(dp), intent(IN)    :: H_ice(:,:)           ! [m]
        real(dp), intent(IN)    :: z_bed(:,:)           ! [m]
        real(dp), intent(IN)    :: mask(:,:)            ! 1 on grounded ice, else 0
        real(dp), intent(IN)    :: G(:,:)               ! geothermal heat flux into the bed [W/m2]
        real(dp), intent(IN)    :: q_T(:,:)             ! conductive heat flux from the bed into the ice [W/m2]
        real(dp), intent(IN)    :: i_eb(:,:)            ! water reaching the bed from above [kg/m2/s]
        real(dp), intent(IN)    :: uxy_b(:,:)           ! basal sliding speed magnitude [m/s]
        real(dp), intent(IN)    :: A_glen(:,:)          ! Glen's A [Pa^-n s^-1]
        real(dp), intent(IN)    :: kappa(:,:)           ! bed type indicator (0 hard, 1 soft)
        real(dp), intent(IN)    :: dx, dy               ! [m] grid spacing
        type(k24_param_class), intent(IN) :: par
        real(dp), intent(IN),    optional :: ux_b(:,:), uy_b(:,:), tau_b_in(:,:), c_till_in(:,:)
        real(dp), intent(INOUT), optional :: diss_io(:,:)
        real(dp), intent(OUT),   optional :: C_frz(:,:)
        real(dp), intent(OUT),   optional :: gsx_out(:,:), gsy_out(:,:), absgs_out(:,:)
        real(dp), intent(OUT),   optional :: absg_out(:,:), phi0_out(:,:)

        type(k24_work_class) :: wk
        integer  :: nx, ny, i, j
        real(dp) :: gmag, epsT

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        if (par%friction_discretization == K24_FRICTION_STAGGERED) then
            if (.not. (present(ux_b) .and. present(uy_b))) then
                write(*,*) "calc_k24:: error: k24_friction_discretization = 1 (staggered) needs ux_b and uy_b."
                stop
            end if
        end if
        if (par%sliding_law == K24_SLIDING_PRESCRIBED_FIELD .and. .not. present(tau_b_in)) then
            write(*,*) "calc_k24:: error: k24_sliding_law = 4 (prescribed field) needs tau_b_in."
            stop
        end if
        if (par%sliding_law == K24_SLIDING_REG_COULOMB_FIELD .and. .not. present(c_till_in)) then
            write(*,*) "calc_k24:: error: k24_sliding_law = 5 (regularized Coulomb field) needs c_till_in."
            stop
        end if

        call k24_work_alloc(wk, nx, ny)

        ! Fixed part of the melt rate, (G - q_T)/L_w [kg/m2/s].
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
            do i = 1, nx
                wk%mdot_fixed(i,j) = (G(i,j) - q_T(i,j)) / par%latent_heat_water
            end do
        end do
        !$omp end parallel do

        if (par%sliding_law == K24_SLIDING_SHAKTI_REG_COULOMB) then
            wk%lambda = par%shakti_lambda_coeff * A_glen
        end if

        if (present(diss_io)) wk%diss = diss_io

        ! ================= update_q! =================

        ! True geometric potential from the RAW ice thickness.
        call update_phi0(wk%phi0, H_ice, z_bed, par)

        ! Routing potential and flow directions (prepare_routing!, routing.jl).
        call prepare_routing(wk, H_ice, z_bed, mask, dx, dy, par)

        ! Correction factor from psi_out to q. Depends only on the routing
        ! directions, so it is fixed for the whole Picard loop below.
        epsT = epsilon(1.0_dp)
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
            do i = 1, nx
                wk%corfac(i,j) = (abs(wk%gsx(i,j)) * dy + abs(wk%gsy(i,j)) * dx) &
                    / (sqrt(wk%gsx(i,j)*wk%gsx(i,j) + wk%gsy(i,j)*wk%gsy(i,j)) + epsT)
            end do
        end do
        !$omp end parallel do

        ! The routing graph just changed: any recorded tape is stale.
        wk%tape_valid = .FALSE.

        call resolve_q(q, N, wk, mask, i_eb, uxy_b, A_glen, kappa, H_ice, dx, dy, par, &
                       ux_b, uy_b, tau_b_in, c_till_in)

        ! ================= update_N! =================
        call update_N(N, q, wk, mask, uxy_b, A_glen, kappa, H_ice, par)

        ! ================= update_W! =================
        call update_W(W, q, wk, mask, par)

        ! ================= outputs / diagnostics =================
        Q_b    = wk%Q_b
        Q_diss = wk%Q_diss
        if (present(diss_io)) diss_io = wk%diss

        ! q_x / q_y have no counterpart in FastHydrology.jl, which carries only
        ! the scalar q. They are a Fortran-side diagnostic: q resolved along
        ! the (unsmoothed) gradient of the routing potential.
        !$omp parallel do default(shared) private(i,j,gmag) schedule(static)
        do j = 1, ny
            do i = 1, nx
                gmag = sqrt(wk%gx(i,j)*wk%gx(i,j) + wk%gy(i,j)*wk%gy(i,j))
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

        if (present(C_frz)) then
            ! GDS-Warner routes with weights computed on the fly, so w8 is
            ! only filled when face fluxes need it.
            if (par%routing_scheme == K24_ROUTE_GDS_WARNER .and. .not. needs_face_fluxes(par)) &
                call compute_routing_weights(wk, mask, dx, dy, par)
            call calc_capacity(C_frz, wk, mask, i_eb, dx, dy, par)
        end if

        if (present(gsx_out))   gsx_out   = wk%gsx
        if (present(gsy_out))   gsy_out   = wk%gsy
        if (present(absgs_out)) absgs_out = wk%abs_gs
        if (present(absg_out))  absg_out  = wk%abs_g
        if (present(phi0_out))  phi0_out  = wk%phi0

        call k24_work_free(wk)

        return

    end subroutine calc_k24

    subroutine k24_N_from_ub(N, q, abs_g, phi0, mask, uxy_b, A_glen, kappa, H_ice, par)
        ! K24's effective pressure for a new sliding speed uxy_b, with the routing of the last
        ! calc_k24 held: q (distributed flux), abs_g (routing-potential gradient) and phi0 (true
        ! geometric potential). This is update_N alone. A host calls it inside its velocity
        ! iteration, so that N and u_b are solved together rather than lagged by a step (a lag
        ! makes them alternate between two states every step: N rises with u_b through the
        ! cavity opening in N_inf, and u_b falls steeply with N through the friction law).
        !
        ! Holding the routing is exact when q does not depend on u_b, i.e. under
        ! K24_SLIDING_NO_FRICTION (no frictional heat in the water source); with a friction law
        ! the routing lags the velocity iteration by one call of calc_k24.
        implicit none
        real(dp),              intent(INOUT) :: N(:,:)
        real(dp),              intent(IN)    :: q(:,:), abs_g(:,:), phi0(:,:), mask(:,:)
        real(dp),              intent(IN)    :: uxy_b(:,:), A_glen(:,:), kappa(:,:), H_ice(:,:)
        type(k24_param_class), intent(IN)    :: par
        type(k24_work_class) :: wk

        call k24_work_alloc(wk, size(N,1), size(N,2))
        wk%abs_g = abs_g
        wk%phi0  = phi0
        call update_N(N, q, wk, mask, uxy_b, A_glen, kappa, H_ice, par)
        call k24_work_free(wk)

        return

    end subroutine k24_N_from_ub

    subroutine calc_capacity(C_frz, wk, mask, i_eb, dx, dy, par)
        ! Freeze-on capacity [m/s ice equivalent] of each grounded cell: the
        ! water that reaches it, i.e. the inflow Psi_in = sum over neighbours n
        ! of psi_out(n) times the fraction of its outflow sent toward the cell,
        ! plus the water from above i_eb [kg/m2/s],
        ! C = (rho_w*Psi_in/(dx*dy) + i_eb)/rho_i.
        ! The cell's own melt and dissipation are heat, already in the host's
        ! freezing demand, so they are left out; i_eb is water, not heat.
        implicit none
        real(dp),             intent(OUT) :: C_frz(:,:)
        type(k24_work_class), intent(IN)  :: wk
        real(dp),             intent(IN)  :: mask(:,:), i_eb(:,:), dx, dy
        type(k24_param_class),intent(IN)  :: par
        integer  :: i, j, d, ni, nj
        real(dp) :: psi_in, supply

        !$omp parallel do default(shared) private(i,j,d,ni,nj,psi_in,supply) schedule(static)
        do j = 1, wk%ny
            do i = 1, wk%nx
                psi_in = 0.0_dp
                supply = 0.0_dp
                if (mask(i,j) == 1.0_dp) then
                    do d = 1, par%n_dirs
                        ni = wrap_index(i + K24_DIRS(1,d), wk%nx, par%periodic_x)
                        nj = wrap_index(j + K24_DIRS(2,d), wk%ny, par%periodic_y)
                        if (ni < 1 .or. ni > wk%nx .or. nj < 1 .or. nj > wk%ny) cycle
                        if (mask(ni,nj) /= 1.0_dp) cycle
                        ! The neighbour sends toward (i,j) in the direction opposite to d.
                        psi_in = psi_in + wk%psi_out(ni,nj) * wk%w8(K24_OPPOSITE(d),ni,nj)
                    end do
                    supply = par%water_density * psi_in / (dx * dy) + i_eb(i,j)
                end if
                C_frz(i,j) = max(supply, 0.0_dp) / par%ice_density
            end do
        end do
        !$omp end parallel do

    end subroutine calc_capacity

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

        allocate(wk%phi0(nx,ny), wk%phi0_filled(nx,ny), wk%phi0_tmp(nx,ny), wk%h(nx,ny))
        allocate(wk%gx(nx,ny),  wk%gy(nx,ny),  wk%abs_g(nx,ny))
        allocate(wk%gsx(nx,ny), wk%gsy(nx,ny), wk%abs_gs(nx,ny))
        allocate(wk%mdot_fixed(nx,ny), wk%mdot_total(nx,ny), wk%psi_out(nx,ny), wk%corfac(nx,ny))
        allocate(wk%q_prev(nx,ny), wk%N_prev(nx,ny), wk%tau_b(nx,ny), wk%lambda(nx,ny))
        allocate(wk%Q_b(nx,ny), wk%Q_diss(nx,ny))
        allocate(wk%visited(nx,ny))
        allocate(wk%w8(8,nx,ny), wk%Fx(nx+1,ny), wk%Fy(nx,ny+1), wk%diss(nx,ny))
        allocate(wk%Q(nx,ny), wk%S_inf(nx,ny))
        allocate(wk%H_hard(nx,ny), wk%H_soft(nx,ny), wk%H_cond(nx,ny))
        allocate(wk%N_inf(nx,ny), wk%Po(nx,ny))

        wk%psi_out = 0.0_dp
        wk%tau_b   = 0.0_dp
        wk%lambda  = 0.0_dp
        wk%Q_b     = 0.0_dp
        wk%Q_diss  = 0.0_dp
        wk%w8      = 0.0_dp
        wk%Fx      = 0.0_dp
        wk%Fy      = 0.0_dp
        wk%diss    = 0.0_dp
        wk%tape_valid = .FALSE.
        wk%tape_n     = 0

    end subroutine k24_work_alloc

    subroutine k24_work_free(wk)
        implicit none
        type(k24_work_class), intent(INOUT) :: wk

        if (allocated(wk%phi0))        deallocate(wk%phi0)
        if (allocated(wk%phi0_filled)) deallocate(wk%phi0_filled)
        if (allocated(wk%phi0_tmp))    deallocate(wk%phi0_tmp)
        if (allocated(wk%h))           deallocate(wk%h)
        if (allocated(wk%gx))          deallocate(wk%gx)
        if (allocated(wk%gy))          deallocate(wk%gy)
        if (allocated(wk%abs_g))       deallocate(wk%abs_g)
        if (allocated(wk%gsx))         deallocate(wk%gsx)
        if (allocated(wk%gsy))         deallocate(wk%gsy)
        if (allocated(wk%abs_gs))      deallocate(wk%abs_gs)
        if (allocated(wk%mdot_fixed))  deallocate(wk%mdot_fixed)
        if (allocated(wk%mdot_total))  deallocate(wk%mdot_total)
        if (allocated(wk%psi_out))     deallocate(wk%psi_out)
        if (allocated(wk%corfac))      deallocate(wk%corfac)
        if (allocated(wk%q_prev))      deallocate(wk%q_prev)
        if (allocated(wk%N_prev))      deallocate(wk%N_prev)
        if (allocated(wk%tau_b))       deallocate(wk%tau_b)
        if (allocated(wk%lambda))      deallocate(wk%lambda)
        if (allocated(wk%Q_b))         deallocate(wk%Q_b)
        if (allocated(wk%Q_diss))      deallocate(wk%Q_diss)
        if (allocated(wk%visited))     deallocate(wk%visited)
        if (allocated(wk%w8))          deallocate(wk%w8)
        if (allocated(wk%Fx))          deallocate(wk%Fx)
        if (allocated(wk%Fy))          deallocate(wk%Fy)
        if (allocated(wk%diss))        deallocate(wk%diss)
        if (allocated(wk%t_dst_i))     deallocate(wk%t_dst_i)
        if (allocated(wk%t_dst_j))     deallocate(wk%t_dst_j)
        if (allocated(wk%t_src_i))     deallocate(wk%t_src_i)
        if (allocated(wk%t_src_j))     deallocate(wk%t_src_j)
        if (allocated(wk%t_w))         deallocate(wk%t_w)
        if (allocated(wk%Q))           deallocate(wk%Q)
        if (allocated(wk%S_inf))       deallocate(wk%S_inf)
        if (allocated(wk%H_hard))      deallocate(wk%H_hard)
        if (allocated(wk%H_soft))      deallocate(wk%H_soft)
        if (allocated(wk%H_cond))      deallocate(wk%H_cond)
        if (allocated(wk%N_inf))       deallocate(wk%N_inf)
        if (allocated(wk%Po))          deallocate(wk%Po)
        if (allocated(wk%in_degree))   deallocate(wk%in_degree)
        if (allocated(wk%stack_i))     deallocate(wk%stack_i)
        if (allocated(wk%stack_j))     deallocate(wk%stack_j)
        if (allocated(wk%stack_k))     deallocate(wk%stack_k)

        wk%nx = 0
        wk%ny = 0
        wk%tape_valid = .FALSE.
        wk%tape_n     = 0

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

        ! DEVIATION: Julia returns NaN for an empty mask. 0 keeps a
        ! fully-floating domain finite.
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
    ! Routing potential and directions (prepare_routing!, routing.jl)
    ! ============================================================
    subroutine prepare_routing(wk, H_ice, z_bed, mask, dx, dy, par)
        ! GDS schemes (GDS_WARNER, GDS_TARBOTON): fill phi0, take its
        ! gradient, smooth the gradient components with the Kamb & Echelmeyer
        ! (1986) kernel. Potential schemes (WARNER, QUINN, TARBOTON,
        ! MODIFIED_TARBOTON): smooth the potential itself with the same
        ! kernel, then fill it; gsx/gsy then hold the unsmoothed gradient of
        ! that surface. abs_g (N, dissipation) is always the gradient of the
        ! TRUE potential.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: H_ice(:,:), z_bed(:,:), mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par

        real(dp), allocatable :: kernel(:,:)
        integer :: i, j, frb_x, frb_y

        select case (par%routing_scheme)

            case (K24_ROUTE_GDS_WARNER, K24_ROUTE_GDS_TARBOTON)

                wk%phi0_filled = wk%phi0
                call potential_filling(wk, z_bed, mask, par)
                call update_potential_gradients(wk, dx, dy, par%periodic_x, par%periodic_y)
                call update_smoothed_potential_gradients(wk, dx, dy, mask, par)
                if (par%routing_scheme /= K24_ROUTE_GDS_WARNER .or. needs_face_fluxes(par)) then
                    call compute_routing_weights(wk, mask, dx, dy, par)
                end if

            case default

                if (par%long_coupling_water == 0.0_dp) then
                    wk%phi0_filled = wk%phi0
                else
                    call coupling_kernel(kernel, frb_x, frb_y, max(masked_mean(H_ice, mask), 10.0_dp), dx, dy, par)
                    call imfilter_replicate_fftw(wk%phi0, wk%nx, wk%ny, kernel, frb_x, frb_y, wk%phi0_filled, &
                                                 par%periodic_x, par%periodic_y)
                    deallocate(kernel)
                end if
                call potential_filling(wk, z_bed, mask, par)
                call update_potential_gradients(wk, dx, dy, par%periodic_x, par%periodic_y)
                wk%gsx = wk%gx
                wk%gsy = wk%gy
                !$omp parallel do default(shared) private(i,j) schedule(static)
                do j = 1, wk%ny
                    do i = 1, wk%nx
                        wk%abs_gs(i,j) = abs(wk%gsx(i,j)) + abs(wk%gsy(i,j))
                    end do
                end do
                !$omp end parallel do
                call compute_routing_weights(wk, mask, dx, dy, par)

        end select

    end subroutine prepare_routing

    logical function needs_face_fluxes(par)
        ! Face fluxes are needed each sweep for face-average q, or for
        ! face-assembled dissipation with the dissipation melt on.
        implicit none
        type(k24_param_class), intent(IN) :: par
        needs_face_fluxes = (par%q_conv == K24_QCONV_FACE_AVERAGE) .or. &
                            (par%dissipation_melt .and. par%diss_disc == K24_DISS_FACE)
    end function needs_face_fluxes

    ! ============================================================
    ! Pit filling of the routing potential (potential_filling!)
    ! ============================================================
    subroutine potential_filling(wk, z_bed, mask, par)
        ! Removes local minima of wk%phi0_filled (which holds the surface to
        ! fill on entry) per par%fill_alg, then sets wk%h to the ice thickness
        ! consistent with the filled potential (only used for the mean
        ! thickness in the GDS smoothing kernel). wk%phi0 is untouched.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: z_bed(:,:), mask(:,:)
        type(k24_param_class),intent(IN)    :: par
        integer :: i, j

        select case (par%fill_alg)
            case (K24_FILL_JACOBI)
                call fill_jacobi(wk%phi0_filled, wk%phi0_tmp, par%fill_iters, .FALSE., par%periodic_x, par%periodic_y)
            case (K24_FILL_LOWEST_NEIGHBOUR)
                call fill_jacobi(wk%phi0_filled, wk%phi0_tmp, par%fill_iters, .TRUE., par%periodic_x, par%periodic_y)
            case (K24_FILL_PRIORITY_FLOOD)
                call fill_priority_flood(wk, mask, par)
            case default
                write(*,*) "potential_filling:: error: k24_fill_algorithm must be one of [-1,0,1,2]."
                write(*,*) "fill_algorithm = ", par%fill_algorithm
                stop
        end select

        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, wk%ny
            do i = 1, wk%nx
                wk%h(i,j) = (wk%phi0_filled(i,j) - par%water_density * par%gravity * z_bed(i,j)) &
                          / (par%ice_density * par%gravity)
            end do
        end do
        !$omp end parallel do

    end subroutine potential_filling

    subroutine fill_jacobi(phi0, phi0_tmp, iterations, lowest, periodic_x, periodic_y)
        ! JacobiFill / LowestNeighbourFill: `iterations` passes of raising every
        ! strict local minimum to the mean of its 4 neighbours (or, `lowest`,
        ! to its lowest neighbour). Out-of-range neighbours are edge-replicated
        ! (index clamped), so an edge cell is never a strict minimum; in a
        ! periodic direction they wrap instead (Fortran-only). Julia only
        ! re-checks the cells next to the last pass's fills, which gives the
        ! same result as this full scan.
        implicit none
        real(dp), intent(INOUT) :: phi0(:,:)
        real(dp), intent(INOUT) :: phi0_tmp(:,:)
        integer,  intent(IN)    :: iterations
        logical,  intent(IN)    :: lowest, periodic_x, periodic_y

        integer  :: iter, i, j, nx, ny, im1, ip1, jm1, jp1
        real(dp) :: p, p1, p2, p3, p4

        nx = size(phi0,1); ny = size(phi0,2)

        phi0_tmp = phi0

        do iter = 1, iterations
            !$omp parallel do default(shared) private(i,j,im1,ip1,jm1,jp1,p,p1,p2,p3,p4) schedule(static)
            do j = 1, ny
                do i = 1, nx
                    p   = phi0(i,j)
                    im1 = max(wrap_index(i-1, nx, periodic_x), 1);  ip1 = min(wrap_index(i+1, nx, periodic_x), nx)
                    jm1 = max(wrap_index(j-1, ny, periodic_y), 1);  jp1 = min(wrap_index(j+1, ny, periodic_y), ny)
                    p1  = phi0(ip1,j);  p2 = phi0(im1,j)
                    p3  = phi0(i,jp1);  p4 = phi0(i,jm1)
                    if (p < p1 .and. p < p2 .and. p < p3 .and. p < p4) then
                        if (lowest) then
                            phi0_tmp(i,j) = min(p1, p2, p3, p4)
                        else
                            phi0_tmp(i,j) = (p1 + p2 + p3 + p4) / 4.0_dp
                        end if
                    end if
                end do
            end do
            !$omp end parallel do
            phi0 = phi0_tmp
        end do

    end subroutine fill_jacobi

    subroutine fill_priority_flood(wk, mask, par)
        ! PriorityFloodFill: Priority-Flood+epsilon (Barnes, Lehman & Mulla
        ! 2014) restricted to grounded cells. Outlets (seeds): every
        ! non-grounded cell next to a grounded one and every grounded cell on
        ! a (non-periodic) domain edge. Cells are visited in increasing order
        ! of filled potential; each newly reached grounded neighbour not
        ! higher than the cell it was reached from is raised to that value
        ! plus epsilon. The heap is the same binary min-heap as Julia's
        ! heap_push!/heap_pop!, so ties pop in the same order.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:)
        type(k24_param_class),intent(IN)    :: par

        integer  :: nx, ny, i, j, ni, nj, d, k, n
        integer, allocatable  :: heap_k(:)
        real(dp), allocatable :: heap_p(:)
        real(dp) :: pc, eps_pf
        logical  :: seed

        nx = wk%nx; ny = wk%ny
        eps_pf = par%priority_flood_epsilon

        allocate(heap_k(nx*ny), heap_p(nx*ny))
        allocate(wk%in_degree(nx,ny))   ! reused as the "already queued" flag
        wk%in_degree = 0
        n = 0

        do j = 1, ny
            do i = 1, nx
                if (mask(i,j) == 1.0_dp) then
                    seed = ((.not. par%periodic_x) .and. (i == 1 .or. i == nx)) .or. &
                           ((.not. par%periodic_y) .and. (j == 1 .or. j == ny))
                else
                    seed = .FALSE.
                    do d = 1, 4
                        ni = wrap_index(i + K24_DIRS(1,d), nx, par%periodic_x)
                        nj = wrap_index(j + K24_DIRS(2,d), ny, par%periodic_y)
                        if (ni < 1 .or. ni > nx .or. nj < 1 .or. nj > ny) cycle
                        if (mask(ni,nj) == 1.0_dp) then
                            seed = .TRUE.
                            exit
                        end if
                    end do
                end if
                if (seed) then
                    wk%in_degree(i,j) = 1
                    call heap_push(heap_p, heap_k, n, wk%phi0_filled(i,j), i + (j-1)*nx)
                end if
            end do
        end do

        do while (n > 0)
            call heap_pop(heap_p, heap_k, n, pc, k)
            i = modulo(k-1, nx) + 1
            j = (k-1) / nx + 1
            do d = 1, 4
                ni = wrap_index(i + K24_DIRS(1,d), nx, par%periodic_x)
                nj = wrap_index(j + K24_DIRS(2,d), ny, par%periodic_y)
                if (ni < 1 .or. ni > nx .or. nj < 1 .or. nj > ny) cycle
                if (wk%in_degree(ni,nj) /= 0 .or. mask(ni,nj) /= 1.0_dp) cycle
                wk%in_degree(ni,nj) = 1
                if (wk%phi0_filled(ni,nj) <= pc) wk%phi0_filled(ni,nj) = pc + eps_pf
                call heap_push(heap_p, heap_k, n, wk%phi0_filled(ni,nj), ni + (nj-1)*nx)
            end do
        end do

        deallocate(heap_k, heap_p, wk%in_degree)

    end subroutine fill_priority_flood

    subroutine heap_push(p, k, n, pv, kv)
        ! heap_push! (water_flux.jl): append, then sift up.
        implicit none
        real(dp), intent(INOUT) :: p(:)
        integer,  intent(INOUT) :: k(:), n
        real(dp), intent(IN)    :: pv
        integer,  intent(IN)    :: kv
        integer  :: c, par_i, ktmp
        real(dp) :: ptmp

        n = n + 1
        p(n) = pv; k(n) = kv
        c = n
        do while (c > 1)
            par_i = c / 2
            if (p(par_i) <= p(c)) exit
            ptmp = p(par_i); p(par_i) = p(c); p(c) = ptmp
            ktmp = k(par_i); k(par_i) = k(c); k(c) = ktmp
            c = par_i
        end do
    end subroutine heap_push

    subroutine heap_pop(p, k, n, top_p, top_k)
        ! heap_pop! (water_flux.jl): take the root, move the last element to
        ! the root and sift down.
        implicit none
        real(dp), intent(INOUT) :: p(:)
        integer,  intent(INOUT) :: k(:), n
        real(dp), intent(OUT)   :: top_p
        integer,  intent(OUT)   :: top_k
        integer  :: c, l, r, m, last_k
        real(dp) :: last_p

        top_p = p(1); top_k = k(1)
        last_p = p(n); last_k = k(n)
        n = n - 1
        if (n > 0) then
            c = 1
            do
                l = 2*c; r = l + 1
                if (l > n) exit
                if (r <= n) then
                    if (p(r) < p(l)) then
                        m = r
                    else
                        m = l
                    end if
                else
                    m = l
                end if
                if (.not. (p(m) < last_p)) exit
                p(c) = p(m); k(c) = k(m)
                c = m
            end do
            p(c) = last_p; k(c) = last_k
        end if
    end subroutine heap_pop

    ! ============================================================
    ! Potential gradients
    ! ============================================================
    subroutine minus_gradients(f, gx, gy, dx, dy, periodic_x, periodic_y)
        ! minus_gradient_x_kernel!/minus_gradient_y_kernel! (grid.jl): centred
        ! differences in the interior, one-sided (full-gradient) differences
        ! at the domain edges, 0 on a single-cell axis. In a periodic
        ! direction the centred difference wraps instead (Fortran-only).
        implicit none
        real(dp), intent(IN)  :: f(:,:)
        real(dp), intent(OUT) :: gx(:,:), gy(:,:)
        real(dp), intent(IN)  :: dx, dy
        logical,  intent(IN)  :: periodic_x, periodic_y
        integer :: i, j, nx, ny

        nx = size(f,1); ny = size(f,2)

        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
            do i = 1, nx
                if (nx == 1) then
                    gx(i,j) = 0.0_dp
                else if (periodic_x) then
                    gx(i,j) = -(f(wrap_index(i+1,nx,.TRUE.),j) - f(wrap_index(i-1,nx,.TRUE.),j)) / (2.0_dp*dx)
                else if (i == 1) then
                    gx(i,j) = -(f(2,j) - f(1,j)) / dx
                else if (i == nx) then
                    gx(i,j) = -(f(nx,j) - f(nx-1,j)) / dx
                else
                    gx(i,j) = -(f(i+1,j) - f(i-1,j)) / (2.0_dp*dx)
                end if

                if (ny == 1) then
                    gy(i,j) = 0.0_dp
                else if (periodic_y) then
                    gy(i,j) = -(f(i,wrap_index(j+1,ny,.TRUE.)) - f(i,wrap_index(j-1,ny,.TRUE.))) / (2.0_dp*dy)
                else if (j == 1) then
                    gy(i,j) = -(f(i,2) - f(i,1)) / dy
                else if (j == ny) then
                    gy(i,j) = -(f(i,ny) - f(i,ny-1)) / dy
                else
                    gy(i,j) = -(f(i,j+1) - f(i,j-1)) / (2.0_dp*dy)
                end if
            end do
        end do
        !$omp end parallel do

    end subroutine minus_gradients

    subroutine update_potential_gradients(wk, dx, dy, periodic_x, periodic_y)
        ! update_potential_gradients! (water_flux.jl): abs_g is the magnitude
        ! of the gradient of the TRUE potential (S_inf, N_inf, dissipation);
        ! gx/gy are then the components of the FILLED potential's gradient,
        ! which only feed the routing.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: dx, dy
        logical,              intent(IN)    :: periodic_x, periodic_y
        integer :: i, j

        call minus_gradients(wk%phi0, wk%gx, wk%gy, dx, dy, periodic_x, periodic_y)

        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, wk%ny
            do i = 1, wk%nx
                wk%abs_g(i,j) = sqrt(wk%gx(i,j)*wk%gx(i,j) + wk%gy(i,j)*wk%gy(i,j))
            end do
        end do
        !$omp end parallel do

        call minus_gradients(wk%phi0_filled, wk%gx, wk%gy, dx, dy, periodic_x, periodic_y)

    end subroutine update_potential_gradients

    ! ============================================================
    ! Stress-gradient-coupling kernel and smoothing
    ! ============================================================
    subroutine coupling_kernel(kernel, frb_x, frb_y, h_avg, dx, dy, par)
        ! coupling_kernel (water_flux.jl): the normalised Kamb & Echelmeyer
        ! (1986) cone for mean grounded-ice thickness h_avg.
        implicit none
        real(dp), allocatable, intent(OUT) :: kernel(:,:)
        integer,               intent(OUT) :: frb_x, frb_y
        real(dp),              intent(IN)  :: h_avg, dx, dy
        type(k24_param_class), intent(IN)  :: par

        real(dp) :: scale, width, delta_min, dist, kernel_sum
        integer  :: maxlevel_x, maxlevel_y, ni, nj

        scale = h_avg * par%long_coupling_water * 2.0_dp

        ! Radius of the cone base (4*h_avg*longcoupwater). The effective
        ! coupling length is the kernel's 2D area-weighted mean distance from
        ! the centre, width/2 = 2*h_avg*longcoupwater = coupling_length_kamb86
        ! * h_avg (10x ice thickness at the default). A 1D average over r,
        ! ignoring the r dr area element, would give width/3 instead.
        width = 2.0_dp * scale

        delta_min = min(dx, dy)
        if (width <= delta_min) then
            ! Bump the cone scale so the kernel spans at least ~1 cell in the
            ! tighter direction. `width` is deliberately NOT recomputed,
            ! matching FastHydrology.jl.
            scale = delta_min / 2.0_dp + 1.0_dp
        end if

        maxlevel_x = 2 * round_half_even(width / dx - 0.5_dp) + 1
        maxlevel_y = 2 * round_half_even(width / dy - 0.5_dp) + 1
        frb_x = (maxlevel_x - 1) / 2
        frb_y = (maxlevel_y - 1) / 2

        allocate(kernel(maxlevel_x, maxlevel_y))
        do nj = 1, maxlevel_y
            do ni = 1, maxlevel_x
                dist = sqrt( (dx * real(ni - frb_x - 1, dp))**2 + &
                             (dy * real(nj - frb_y - 1, dp))**2 ) / scale
                kernel(ni,nj) = max(0.0_dp, 1.0_dp - dist / 2.0_dp)
            end do
        end do
        kernel_sum = sum(kernel)
        if (kernel_sum > 0.0_dp) kernel = kernel / kernel_sum

    end subroutine coupling_kernel

    subroutine update_smoothed_potential_gradients(wk, dx, dy, mask, par)
        ! update_smoothed_potential_gradients! (water_flux.jl).
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: dx, dy, mask(:,:)
        type(k24_param_class),intent(IN)    :: par

        real(dp), allocatable :: kernel(:,:)
        integer  :: frb_x, frb_y, i, j, nx, ny

        nx = wk%nx; ny = wk%ny

        if (par%long_coupling_water == 0.0_dp) then
            wk%gsx = wk%gx
            wk%gsy = wk%gy
        else
            ! Mean grounded-ice thickness, from the potential-filled h.
            call coupling_kernel(kernel, frb_x, frb_y, max(masked_mean(wk%h, mask), 10.0_dp), dx, dy, par)
            call imfilter_replicate_fftw(wk%gx, nx, ny, kernel, frb_x, frb_y, wk%gsx, par%periodic_x, par%periodic_y)
            call imfilter_replicate_fftw(wk%gy, nx, ny, kernel, frb_x, frb_y, wk%gsy, par%periodic_x, par%periodic_y)
            deallocate(kernel)
        end if

        ! L1 magnitude, matching abs_grad_phi0_s in water_flux.jl.
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
            do i = 1, nx
                wk%abs_gs(i,j) = abs(wk%gsx(i,j)) + abs(wk%gsy(i,j))
            end do
        end do
        !$omp end parallel do

    end subroutine update_smoothed_potential_gradients

    integer function round_half_even(x) result(r)
        ! Round half to even, matching Julia's round(Int, x).
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
    subroutine imfilter_replicate_fftw(input, nx, ny, kernel, frb_x, frb_y, output, periodic_x, periodic_y)
        ! cached_fft_convolve! (fft_convolution.jl) in Fortran: embed the input
        ! in an (nx+2*frb_x, ny+2*frb_y) array with REPLICATE padding, wrap the
        ! centred kernel around the array's origin, multiply in the frequency
        ! domain, and crop back. In a periodic direction the padding wraps
        ! instead (Fortran-only).
        !
        ! The FFT array is zero-padded up to the next 7-smooth size (only
        ! prime factors 2, 3, 5, 7), Julia's fft_size, for FFTW speed. That
        ! cannot change the result: an output cell only reads inputs within
        ! frb of it, never the zeros or a wrapped-around neighbour.
        !
        ! Link with -lfftw3.
        use, intrinsic :: iso_c_binding
        implicit none
        include 'fftw3.f03'

        integer,  intent(IN)  :: nx, ny, frb_x, frb_y
        real(dp), intent(IN)  :: input(nx, ny)
        real(dp), intent(IN)  :: kernel(2*frb_x+1, 2*frb_y+1)
        real(dp), intent(OUT) :: output(nx, ny)
        logical,  intent(IN)  :: periodic_x, periodic_y

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

        !$omp parallel do default(shared) private(i,j,ii,jj) schedule(static)
        do j = 1, Npy
            do i = 1, Npx
                ii = min(max(wrap_index(i - frb_x, nx, periodic_x), 1), nx)
                jj = min(max(wrap_index(j - frb_y, ny, periodic_y), 1), ny)
                work_a(i, j) = input(ii, jj)
            end do
        end do
        !$omp end parallel do

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
        ! Smallest m >= n whose only prime factors are 2, 3, 5 and 7 --
        ! fft_size in fft_convolution.jl (nextprod((2, 3, 5, 7), n)).
        implicit none
        integer, intent(IN) :: n
        integer :: k

        m = max(n, 1)
        do
            k = m
            do while (modulo(k, 2) == 0); k = k / 2; end do
            do while (modulo(k, 3) == 0); k = k / 3; end do
            do while (modulo(k, 5) == 0); k = k / 5; end do
            do while (modulo(k, 7) == 0); k = k / 7; end do
            if (k == 1) return
            m = m + 1
        end do
    end function next_smooth_size

    pure integer function wrap_index(k, n, periodic) result(kw)
        ! Neighbour index k of a 1..n axis: wrapped with period n in a
        ! periodic direction, otherwise returned as is (the caller clamps or
        ! skips an out-of-range index).
        implicit none
        integer, intent(IN) :: k, n
        logical, intent(IN) :: periodic

        if (periodic) then
            kw = modulo(k-1, n) + 1
        else
            kw = k
        end if
    end function wrap_index

    ! ============================================================
    ! Routing weights (routing.jl)
    ! ============================================================
    pure real(dp) function routing_weight(sx, sy, di, dj, dx, dy) result(w)
        ! Fraction of a cell's psi_out leaving through the face toward the
        ! neighbour in direction (di, dj), given the cell's routing direction
        ! (sx, sy): (sx*di*dy + sy*dj*dx) / (|sx|*dy + |sy|*dx). Positive only
        ! when the flow leaves toward that neighbour; the denominator is
        ! corfac*|s|, consistent with q = psi_out/corfac.
        implicit none
        real(dp), intent(IN) :: sx, sy, dx, dy
        integer,  intent(IN) :: di, dj
        w = (sx * real(di, dp) * dy + sy * real(dj, dp) * dx) &
            / (abs(sx) * dy + abs(sy) * dx + epsilon(1.0_dp))
    end function routing_weight

    subroutine compute_routing_weights(wk, mask, dx, dy, par)
        ! compute_routing_weights!: w8(d,i,j), the fraction of grounded cell
        ! (i,j)'s outflow sent to its neighbour in direction d. Fractions sum
        ! to 1, or 0 for a sink. Water sent to a non-grounded neighbour or off
        ! the domain edge leaves the system.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par
        integer :: i, j

        wk%w8 = 0.0_dp
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, wk%ny
            do i = 1, wk%nx
                if (mask(i,j) /= 1.0_dp) cycle
                select case (par%routing_scheme)
                    case (K24_ROUTE_GDS_WARNER)
                        call weights_gds_warner(wk, i, j, dx, dy)
                    case (K24_ROUTE_WARNER)
                        call weights_warner(wk, i, j, dx, dy, par)
                    case (K24_ROUTE_QUINN)
                        call weights_quinn(wk, i, j, dx, dy, par)
                    case (K24_ROUTE_TARBOTON)
                        call weights_tarboton(wk, i, j, dx, dy, par)
                    case default   ! MODIFIED_TARBOTON, GDS_TARBOTON
                        call angle_split(wk, i, j, wk%gsx(i,j), wk%gsy(i,j), dx, dy)
                end select
            end do
        end do
        !$omp end parallel do

    end subroutine compute_routing_weights

    subroutine neighbour(i, j, d, nx, ny, par, ni, nj, inside)
        implicit none
        integer, intent(IN)  :: i, j, d, nx, ny
        type(k24_param_class), intent(IN) :: par
        integer, intent(OUT) :: ni, nj
        logical, intent(OUT) :: inside
        ni = wrap_index(i + K24_DIRS(1,d), nx, par%periodic_x)
        nj = wrap_index(j + K24_DIRS(2,d), ny, par%periodic_y)
        inside = (ni >= 1 .and. ni <= nx .and. nj >= 1 .and. nj <= ny)
    end subroutine neighbour

    subroutine weights_gds_warner(wk, i, j, dx, dy)
        ! GDS-Warner: the routing direction's component toward each of the 4
        ! neighbours (identical to routing_weight as the routing uses it).
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        integer,  intent(IN) :: i, j
        real(dp), intent(IN) :: dx, dy
        integer :: d
        do d = 1, 4
            wk%w8(d,i,j) = max(0.0_dp, routing_weight(wk%gsx(i,j), wk%gsy(i,j), K24_DIRS(1,d), K24_DIRS(2,d), dx, dy))
        end do
    end subroutine weights_gds_warner

    subroutine weights_warner(wk, i, j, dx, dy, par)
        ! Warner (Budd & Warner 1996; Le Brocq Eq. 8): shared among the
        ! downhill 4-neighbours in proportion to the potential drop, weighted
        ! by face length / distance.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        integer,  intent(IN) :: i, j
        real(dp), intent(IN) :: dx, dy
        type(k24_param_class), intent(IN) :: par
        integer  :: d, ni, nj
        logical  :: inside
        real(dp) :: p, tot, drop, v

        p = wk%phi0_filled(i,j)
        tot = 0.0_dp
        do d = 1, 4
            call neighbour(i, j, d, wk%nx, wk%ny, par, ni, nj, inside)
            if (.not. inside) cycle
            drop = p - wk%phi0_filled(ni,nj)
            if (drop > 0.0_dp) then
                if (d <= 2) then
                    v = drop * (dy / dx)
                else
                    v = drop * (dx / dy)
                end if
                wk%w8(d,i,j) = v
                tot = tot + v
            end if
        end do
        if (tot > 0.0_dp) then
            do d = 1, 4
                wk%w8(d,i,j) = wk%w8(d,i,j) / tot
            end do
        end if
    end subroutine weights_warner

    subroutine weights_quinn(wk, i, j, dx, dy, par)
        ! Quinn et al. (1991): as Warner over all 8 neighbours. Not
        ! original: shares proportional to the drop (Le Brocq Eq. 8).
        ! Original: slope times effective contour length.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        integer,  intent(IN) :: i, j
        real(dp), intent(IN) :: dx, dy
        type(k24_param_class), intent(IN) :: par
        integer  :: d, ni, nj
        logical  :: inside
        real(dp) :: p, tot, drop, v, ddiag, Ldiag

        p = wk%phi0_filled(i,j)
        tot = 0.0_dp
        ddiag = hypot(dx, dy)
        Ldiag = sqrt(2.0_dp) / 4.0_dp * sqrt(dx * dy)
        do d = 1, 8
            call neighbour(i, j, d, wk%nx, wk%ny, par, ni, nj, inside)
            if (.not. inside) cycle
            drop = p - wk%phi0_filled(ni,nj)
            if (drop > 0.0_dp) then
                if (.not. par%quinn_original) then
                    v = drop
                else if (d <= 2) then
                    v = drop / dx * (dy / 2.0_dp)
                else if (d <= 4) then
                    v = drop / dy * (dx / 2.0_dp)
                else
                    v = drop / ddiag * Ldiag
                end if
                wk%w8(d,i,j) = v
                tot = tot + v
            end if
        end do
        if (tot > 0.0_dp) then
            do d = 1, 8
                wk%w8(d,i,j) = wk%w8(d,i,j) / tot
            end do
        end if
    end subroutine weights_quinn

    subroutine weights_tarboton(wk, i, j, dx, dy, par)
        ! Tarboton (1997) D-infinity: steepest downhill direction over the 8
        ! triangular facets; the outflow goes to the facet's cardinal and
        ! diagonal neighbours in proportion to the flow angle.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        integer,  intent(IN) :: i, j
        real(dp), intent(IN) :: dx, dy
        type(k24_param_class), intent(IN) :: par
        integer  :: f, c1, c2, best, i1, j1, i2, j2
        logical  :: in1, in2
        real(dp) :: e0, e1, e2, d1, d2, s1, s2, amax, r, s
        real(dp) :: best_s, best_r, best_amax

        e0 = wk%phi0_filled(i,j)
        best_s = 0.0_dp; best = 0; best_r = 0.0_dp; best_amax = 1.0_dp
        do f = 1, 8
            c1 = K24_FACETS(1,f); c2 = K24_FACETS(2,f)
            call neighbour(i, j, c1, wk%nx, wk%ny, par, i1, j1, in1)
            call neighbour(i, j, c2, wk%nx, wk%ny, par, i2, j2, in2)
            if (.not. (in1 .and. in2)) cycle
            e1 = wk%phi0_filled(i1,j1); e2 = wk%phi0_filled(i2,j2)
            if (c1 <= 2) then
                d1 = dx; d2 = dy
            else
                d1 = dy; d2 = dx
            end if
            s1 = (e0 - e1) / d1
            s2 = (e1 - e2) / d2
            amax = atan2(d2, d1)
            r = atan2(s2, s1)
            s = hypot(s1, s2)
            if (r < 0.0_dp) then
                r = 0.0_dp; s = s1
            else if (r > amax) then
                r = amax; s = (e0 - e2) / hypot(d1, d2)
            end if
            if (s > best_s) then
                best_s = s; best = f; best_r = r; best_amax = amax
            end if
        end do
        if (best == 0) return
        c1 = K24_FACETS(1,best); c2 = K24_FACETS(2,best)
        wk%w8(c2,i,j) = best_r / best_amax
        wk%w8(c1,i,j) = 1.0_dp - best_r / best_amax
    end subroutine weights_tarboton

    subroutine angle_split(wk, i, j, fx, fy, dx, dy)
        ! Modified Tarboton / GDS-Tarboton (Le Brocq Sec. 3): one flow
        ! direction (fx, fy) split between the two of the 8 neighbours whose
        ! directions bracket it, in proportion to the angles.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        integer,  intent(IN) :: i, j
        real(dp), intent(IN) :: fx, fy, dx, dy
        integer, parameter :: dirs(8) = [2, 8, 4, 7, 1, 5, 3, 6]
        real(dp) :: a, angs(8), theta, lo, hi, frac
        integer  :: k, m, k2

        if (fx == 0.0_dp .and. fy == 0.0_dp) return
        a = atan2(dy, dx)
        angs = [0.0_dp, a, K24_PI / 2.0_dp, K24_PI - a, K24_PI * 1.0_dp, K24_PI + a, &
                3.0_dp * K24_PI / 2.0_dp, 2.0_dp * K24_PI - a]
        ! mod(atan(fy, fx), 2pi): atan2 is in (-pi, pi], so this is exact.
        theta = atan2(fy, fx)
        if (theta < 0.0_dp) theta = theta + 2.0_dp * K24_PI
        k = 8
        do m = 1, 7
            if (angs(m) <= theta .and. theta < angs(m+1)) then
                k = m
                exit
            end if
        end do
        lo = angs(k)
        if (k == 8) then
            hi = 2.0_dp * K24_PI
            k2 = 1
        else
            hi = angs(k+1)
            k2 = k + 1
        end if
        frac = (theta - lo) / (hi - lo)
        wk%w8(dirs(k),i,j)  = wk%w8(dirs(k),i,j)  + (1.0_dp - frac)
        wk%w8(dirs(k2),i,j) = wk%w8(dirs(k2),i,j) + frac
    end subroutine angle_split

    ! ============================================================
    ! Sliding law -> basal shear stress, and the frictional heat
    ! ============================================================
    subroutine update_tau_b(tau_b, N, uxy_b, lambda, par, tau_b_in, c_till_in)
        ! Mirrors update_tau_b! (sliding_law.jl).
        implicit none
        real(dp), intent(OUT) :: tau_b(:,:)
        real(dp), intent(IN)  :: N(:,:), uxy_b(:,:), lambda(:,:)
        type(k24_param_class), intent(IN) :: par
        real(dp), intent(IN), optional :: tau_b_in(:,:), c_till_in(:,:)

        integer  :: i, j, nx, ny
        real(dp) :: C, qe, c_till, u0, n_s, inv_n

        nx = size(tau_b,1); ny = size(tau_b,2)

        select case (par%sliding_law)

            case (K24_SLIDING_NO_FRICTION)
                tau_b = 0.0_dp

            case (K24_SLIDING_PRESCRIBED_FIELD)
                tau_b = tau_b_in

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

            case (K24_SLIDING_REG_COULOMB_FIELD)
                qe = par%reg_coulomb_q
                u0 = par%reg_coulomb_u0
                !$omp parallel do default(shared) private(i,j) schedule(static)
                do j = 1, ny
                    do i = 1, nx
                        tau_b(i,j) = c_till_in(i,j) * N(i,j) * &
                                     (uxy_b(i,j) / (uxy_b(i,j) + u0))**qe
                    end do
                end do
                !$omp end parallel do

            case (K24_SLIDING_SHAKTI_REG_COULOMB)
                C     = par%shakti_C
                n_s   = par%shakti_n
                inv_n = 1.0_dp / n_s
                !$omp parallel do default(shared) private(i,j) schedule(static)
                do j = 1, ny
                    do i = 1, nx
                        tau_b(i,j) = C * N(i,j) * &
                            (uxy_b(i,j) / (uxy_b(i,j) + abs(N(i,j))**n_s * lambda(i,j)))**inv_n
                    end do
                end do
                !$omp end parallel do

            case default
                write(*,*) "update_tau_b:: error: k24_sliding_law must be one of [0,1,2,3,4,5,6]."
                write(*,*) "sliding_law = ", par%sliding_law
                stop

        end select

    end subroutine update_tau_b

    logical function pressure_dependent_law(par)
        implicit none
        type(k24_param_class), intent(IN) :: par
        pressure_dependent_law = (par%sliding_law == K24_SLIDING_POWER_PLASTIC)     .or. &
                                 (par%sliding_law == K24_SLIDING_REG_COULOMB)       .or. &
                                 (par%sliding_law == K24_SLIDING_REG_COULOMB_FIELD) .or. &
                                 (par%sliding_law == K24_SLIDING_SHAKTI_REG_COULOMB)
    end function pressure_dependent_law

    subroutine add_friction_term(wk, uxy_b, par, ux_b, uy_b)
        ! add_friction_term!: Q_b [W/m2] from the current tau_b, per
        ! friction_discretization, and Q_b/L_w added to mdot_total.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: uxy_b(:,:)
        type(k24_param_class),intent(IN)    :: par
        real(dp), intent(IN), optional      :: ux_b(:,:), uy_b(:,:)
        integer :: i, j

        if (par%friction_discretization == K24_FRICTION_STAGGERED) then
            call staggered_friction(wk%Q_b, wk%tau_b, uxy_b, ux_b, uy_b, par)
        else
            !$omp parallel do default(shared) private(i,j) schedule(static)
            do j = 1, wk%ny
                do i = 1, wk%nx
                    wk%Q_b(i,j) = wk%tau_b(i,j) * uxy_b(i,j)
                end do
            end do
            !$omp end parallel do
        end if

        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, wk%ny
            do i = 1, wk%nx
                wk%mdot_total(i,j) = wk%mdot_total(i,j) + wk%Q_b(i,j) / par%latent_heat_water
            end do
        end do
        !$omp end parallel do

    end subroutine add_friction_term

    subroutine staggered_friction(Q_b, tau_b, uxy_b, ux, uy, par)
        ! staggered_friction_kernel! (sliding_law.jl): C-grid frictional heat.
        ! beta = tau_b/max(|u_b|, u_floor) at centres is averaged to each face
        ! to form the face tractions; the heat is formed on the faces
        ! (quadrature off) or at 2x2 Gauss points (Yelmo's qb_method = 2).
        ! Domain-edge faces use the one cell they have; periodic directions
        ! wrap (Fortran-only).
        implicit none
        real(dp), intent(OUT) :: Q_b(:,:)
        real(dp), intent(IN)  :: tau_b(:,:), uxy_b(:,:), ux(:,:), uy(:,:)
        type(k24_param_class), intent(IN) :: par

        real(dp), allocatable :: tx(:,:), ty(:,:), beta(:,:)
        real(dp) :: cux(4), cuy(4), ctx(4), cty(4), acc, s3, pts(2,4)
        integer  :: i, j, nx, ny, im1, ip1, jm1, jp1, k

        nx = size(Q_b,1); ny = size(Q_b,2)
        allocate(tx(nx,ny), ty(nx,ny), beta(nx,ny))

        do j = 1, ny
            do i = 1, nx
                beta(i,j) = tau_b(i,j) / max(uxy_b(i,j), par%friction_u_floor)
            end do
        end do

        do j = 1, ny
            do i = 1, nx
                ip1 = min(wrap_index(i+1, nx, par%periodic_x), nx)
                jp1 = min(wrap_index(j+1, ny, par%periodic_y), ny)
                tx(i,j) = (beta(i,j) + beta(ip1,j)) / 2.0_dp * ux(i,j)
                ty(i,j) = (beta(i,j) + beta(i,jp1)) / 2.0_dp * uy(i,j)
            end do
        end do

        s3 = 1.0_dp / sqrt(3.0_dp)
        pts = reshape([-s3,-s3,  s3,-s3,  s3,s3,  -s3,s3], [2,4])

        do j = 1, ny
            do i = 1, nx
                im1 = max(wrap_index(i-1, nx, par%periodic_x), 1); ip1 = min(wrap_index(i+1, nx, par%periodic_x), nx)
                jm1 = max(wrap_index(j-1, ny, par%periodic_y), 1); jp1 = min(wrap_index(j+1, ny, par%periodic_y), ny)
                if (par%friction_quadrature) then
                    call acx_corners(ux, i, j, im1, jm1, jp1, cux)
                    call acy_corners(uy, i, j, im1, ip1, jm1, cuy)
                    call acx_corners(tx, i, j, im1, jm1, jp1, ctx)
                    call acy_corners(ty, i, j, im1, ip1, jm1, cty)
                    acc = 0.0_dp
                    do k = 1, 4
                        acc = acc + hypot(gq_interp(cux, pts(:,k)), gq_interp(cuy, pts(:,k))) &
                                  * hypot(gq_interp(ctx, pts(:,k)), gq_interp(cty, pts(:,k)))
                    end do
                    Q_b(i,j) = acc / 4.0_dp
                else
                    Q_b(i,j) = (tx(im1,j) * ux(im1,j) + tx(i,j) * ux(i,j)) / 2.0_dp &
                             + (ty(i,jm1) * uy(i,jm1) + ty(i,j) * uy(i,j)) / 2.0_dp
                end if
            end do
        end do

        deallocate(tx, ty, beta)

    end subroutine staggered_friction

    pure real(dp) function gq_interp(c, p) result(v)
        ! Bilinear shape functions of the SW, SE, NE, NW corners at point p of
        ! the reference cell [-1, 1]^2 (Yelmo's gq2D).
        implicit none
        real(dp), intent(IN) :: c(4), p(2)
        v = ((1.0_dp - p(1)) * (1.0_dp - p(2)) * c(1) + (1.0_dp + p(1)) * (1.0_dp - p(2)) * c(2) + &
             (1.0_dp + p(1)) * (1.0_dp + p(2)) * c(3) + (1.0_dp - p(1)) * (1.0_dp + p(2)) * c(4)) / 4.0_dp
    end function gq_interp

    pure subroutine acx_corners(F, i, j, im1, jm1, jp1, c)
        ! Corner (ab-node) values of an acx field around cell (i, j): SW, SE, NE, NW.
        implicit none
        real(dp), intent(IN)  :: F(:,:)
        integer,  intent(IN)  :: i, j, im1, jm1, jp1
        real(dp), intent(OUT) :: c(4)
        c(1) = (F(im1,jm1) + F(im1,j)) / 2.0_dp
        c(2) = (F(i,jm1)   + F(i,j))   / 2.0_dp
        c(3) = (F(i,j)     + F(i,jp1)) / 2.0_dp
        c(4) = (F(im1,j)   + F(im1,jp1)) / 2.0_dp
    end subroutine acx_corners

    pure subroutine acy_corners(F, i, j, im1, ip1, jm1, c)
        ! Corner (ab-node) values of an acy field around cell (i, j): SW, SE, NE, NW.
        implicit none
        real(dp), intent(IN)  :: F(:,:)
        integer,  intent(IN)  :: i, j, im1, ip1, jm1
        real(dp), intent(OUT) :: c(4)
        c(1) = (F(im1,jm1) + F(i,jm1))   / 2.0_dp
        c(2) = (F(i,jm1)   + F(ip1,jm1)) / 2.0_dp
        c(3) = (F(i,j)     + F(ip1,j))   / 2.0_dp
        c(4) = (F(im1,j)   + F(i,j))     / 2.0_dp
    end subroutine acy_corners

    ! ============================================================
    ! Water-flux fixed point
    ! ============================================================
    subroutine resolve_q(q, N, wk, mask, i_eb, uxy_b, A_glen, kappa, H_ice, dx, dy, par, &
                         ux_b, uy_b, tau_b_in, c_till_in)
        ! Mirrors the three resolve_q! methods (water_flux.jl):
        !   * N-independent law (NO_FRICTION/WEERTMAN/PRESCRIBED_FIELD),
        !     dissipation off -- one pass.
        !   * N-independent law, dissipation on -- Picard on q alone.
        !   * N-dependent law -- joint (q, N) Picard, regardless of the
        !     dissipation setting.
        ! Every sweep rebuilds the source from its terms:
        !   mdot_total = mdot_fixed + i_eb, += Q_b/L_w, += Q_diss/L_w.
        implicit none
        real(dp),             intent(INOUT) :: q(:,:), N(:,:)
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), i_eb(:,:), uxy_b(:,:)
        real(dp),             intent(IN)    :: A_glen(:,:), kappa(:,:), H_ice(:,:)
        real(dp),             intent(IN)    :: dx, dy
        type(k24_param_class),intent(IN)    :: par
        real(dp), intent(IN), optional      :: ux_b(:,:), uy_b(:,:), tau_b_in(:,:), c_till_in(:,:)

        integer  :: iter, n_iters
        real(dp) :: q_scale, N_scale
        logical  :: converged, q_ok, N_ok

        if (pressure_dependent_law(par)) then

            converged = .FALSE.
            n_iters   = par%max_coupling_iters

            do iter = 1, par%max_coupling_iters

                wk%q_prev = q
                wk%N_prev = N

                call update_tau_b(wk%tau_b, N, uxy_b, wk%lambda, par, tau_b_in, c_till_in)
                call reset_source(wk, i_eb)
                call add_friction_term(wk, uxy_b, par, ux_b, uy_b)
                call add_dissipation_term(wk, q, par)

                call update_psi_out(wk, mask, dx, dy, par)
                call update_q_from_psi_out(q, wk, mask, dx, dy, par)

                call update_N(N, q, wk, mask, uxy_b, A_glen, kappa, H_ice, par)

                q_scale = max(masked_max_abs(q, mask), epsilon(1.0_dp))
                N_scale = max(masked_max_abs(N, mask), epsilon(1.0_dp))
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

        else if (par%dissipation_melt) then

            converged = .FALSE.
            n_iters   = par%max_dissipation_iters

            call update_tau_b(wk%tau_b, N, uxy_b, wk%lambda, par, tau_b_in, c_till_in)

            do iter = 1, par%max_dissipation_iters

                wk%q_prev = q

                ! Dissipation from the current q (zero on the first sweep of a
                ! cold start, since q begins at zero).
                call reset_source(wk, i_eb)
                call add_friction_term(wk, uxy_b, par, ux_b, uy_b)
                call add_dissipation_term(wk, q, par)

                call update_psi_out(wk, mask, dx, dy, par)
                call update_q_from_psi_out(q, wk, mask, dx, dy, par)

                q_scale = max(masked_max_abs(q, mask), epsilon(1.0_dp))
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

        else

            ! The source depends on neither q nor N: one routing pass is exact.
            call update_tau_b(wk%tau_b, N, uxy_b, wk%lambda, par, tau_b_in, c_till_in)
            call reset_source(wk, i_eb)
            call add_friction_term(wk, uxy_b, par, ux_b, uy_b)
            wk%Q_diss = 0.0_dp

            call update_psi_out(wk, mask, dx, dy, par)
            call update_q_from_psi_out(q, wk, mask, dx, dy, par)

        end if

    end subroutine resolve_q

    subroutine reset_source(wk, i_eb)
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: i_eb(:,:)
        integer :: i, j
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, wk%ny
            do i = 1, wk%nx
                wk%mdot_total(i,j) = wk%mdot_fixed(i,j) + i_eb(i,j)
            end do
        end do
        !$omp end parallel do
    end subroutine reset_source

    subroutine add_dissipation_term(wk, q, par)
        ! add_dissipation_term!: Q_diss [W/m2] and Q_diss/L_w added to
        ! mdot_total, per dissipation_discretization; Q_diss = 0 when off.
        ! FACE uses the face-assembled dissipation of the previous q
        ! conversion (wk%diss, [kg/m2/s]).
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: q(:,:)
        type(k24_param_class),intent(IN)    :: par
        integer :: i, j

        if (.not. par%dissipation_melt) then
            wk%Q_diss = 0.0_dp
            return
        end if

        if (par%diss_disc == K24_DISS_FACE) then
            !$omp parallel do default(shared) private(i,j) schedule(static)
            do j = 1, wk%ny
                do i = 1, wk%nx
                    wk%Q_diss(i,j)     = wk%diss(i,j) * par%latent_heat_water
                    wk%mdot_total(i,j) = wk%mdot_total(i,j) + wk%diss(i,j)
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(i,j) schedule(static)
            do j = 1, wk%ny
                do i = 1, wk%nx
                    wk%Q_diss(i,j)     = abs(q(i,j) * wk%abs_g(i,j))
                    wk%mdot_total(i,j) = wk%mdot_total(i,j) + wk%Q_diss(i,j) / par%latent_heat_water
                end do
            end do
            !$omp end parallel do
        end if

    end subroutine add_dissipation_term

    ! ============================================================
    ! Routed flux -> q (routing.jl)
    ! ============================================================
    subroutine update_q_from_psi_out(q, wk, mask, dx, dy, par)
        ! update_q_from_psi_out!: face fluxes first if any option needs them,
        ! then q per q_conversion, clamped to [q_min, q_max].
        implicit none
        real(dp),             intent(OUT)   :: q(:,:)
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par
        integer  :: i, j
        real(dp) :: v, qx, qy, epsT

        if (needs_face_fluxes(par)) call update_face_fluxes(wk, mask, dx, dy, par)

        if (par%q_conv == K24_QCONV_FACE_AVERAGE) then
            ! Face fluxes -> face-normal q (Fx/dy, Fy/dx) -> averaged to the
            ! centre per component -> |q|.
            !$omp parallel do default(shared) private(i,j,v,qx,qy) schedule(static)
            do j = 1, wk%ny
                do i = 1, wk%nx
                    v = 0.0_dp
                    if (mask(i,j) == 1.0_dp) then
                        qx = (wk%Fx(i,j) + wk%Fx(i+1,j)) / (2.0_dp * dy)
                        qy = (wk%Fy(i,j) + wk%Fy(i,j+1)) / (2.0_dp * dx)
                        v  = hypot(qx, qy)
                    end if
                    q(i,j) = min(max(v, par%q_min), par%q_max)
                end do
            end do
            !$omp end parallel do
        else
            ! Le Brocq Eq. 9 / K24: q = psi_out / corfac. eps keeps a cell with
            ! no flow direction finite (0 for psi_out == 0).
            epsT = epsilon(1.0_dp)
            !$omp parallel do default(shared) private(i,j) schedule(static)
            do j = 1, wk%ny
                do i = 1, wk%nx
                    q(i,j) = min(max(wk%psi_out(i,j) / (wk%corfac(i,j) + epsT), par%q_min), par%q_max)
                end do
            end do
            !$omp end parallel do
        end if


    end subroutine update_q_from_psi_out

    subroutine update_face_fluxes(wk, mask, dx, dy, par)
        ! face_fluxes_kernel!: net volume flux through every cell face [m3/s]
        ! from psi_out and the outflow fractions (4-neighbour schemes).
        ! Fx(a,j) is the flux through the face between (a-1,j) and (a,j),
        ! positive in +x; Fy likewise in +y. With face dissipation, each cell
        ! gets half the energy F*(phi_up - phi_down) of each of its faces
        ! (true potential), divided by its area and L_w, into wk%diss. In a
        ! periodic direction the two edge faces are the same face
        ! (Fortran-only).
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par
        integer  :: i, j, a, b, nx, ny, il, ir, jl, jr, im1, ip1, jm1, jp1
        real(dp) :: f, P

        nx = wk%nx; ny = wk%ny

        !$omp parallel do default(shared) private(a,j,f,il,ir) schedule(static)
        do j = 1, ny
            do a = 1, nx+1
                f  = 0.0_dp
                il = a - 1
                ir = a
                if (par%periodic_x) then
                    il = wrap_index(il, nx, .TRUE.)
                    ir = wrap_index(ir, nx, .TRUE.)
                end if
                if (il >= 1) then
                    if (mask(il,j) == 1.0_dp) f = f + wk%psi_out(il,j) * wk%w8(2,il,j)
                end if
                if (ir <= nx) then
                    if (mask(ir,j) == 1.0_dp) f = f - wk%psi_out(ir,j) * wk%w8(1,ir,j)
                end if
                wk%Fx(a,j) = f
            end do
        end do
        !$omp end parallel do

        !$omp parallel do default(shared) private(i,b,f,jl,jr) schedule(static)
        do b = 1, ny+1
            do i = 1, nx
                f  = 0.0_dp
                jl = b - 1
                jr = b
                if (par%periodic_y) then
                    jl = wrap_index(jl, ny, .TRUE.)
                    jr = wrap_index(jr, ny, .TRUE.)
                end if
                if (jl >= 1) then
                    if (mask(i,jl) == 1.0_dp) f = f + wk%psi_out(i,jl) * wk%w8(4,i,jl)
                end if
                if (jr <= ny) then
                    if (mask(i,jr) == 1.0_dp) f = f - wk%psi_out(i,jr) * wk%w8(3,i,jr)
                end if
                wk%Fy(i,b) = f
            end do
        end do
        !$omp end parallel do

        if (.not. (par%dissipation_melt .and. par%diss_disc == K24_DISS_FACE)) return

        !$omp parallel do default(shared) private(i,j,P,im1,ip1,jm1,jp1) schedule(static)
        do j = 1, ny
            do i = 1, nx
                if (mask(i,j) /= 1.0_dp) then
                    wk%diss(i,j) = 0.0_dp
                    cycle
                end if
                im1 = wrap_index(i-1, nx, par%periodic_x); ip1 = wrap_index(i+1, nx, par%periodic_x)
                jm1 = wrap_index(j-1, ny, par%periodic_y); jp1 = wrap_index(j+1, ny, par%periodic_y)
                P = 0.0_dp
                if (im1 >= 1)  P = P + wk%Fx(i,j)   * (wk%phi0(im1,j) - wk%phi0(i,j))
                if (ip1 <= nx) P = P + wk%Fx(i+1,j) * (wk%phi0(i,j)   - wk%phi0(ip1,j))
                if (jm1 >= 1)  P = P + wk%Fy(i,j)   * (wk%phi0(i,jm1) - wk%phi0(i,j))
                if (jp1 <= ny) P = P + wk%Fy(i,j+1) * (wk%phi0(i,j)   - wk%phi0(i,jp1))
                wk%diss(i,j) = P / 2.0_dp / (dx * dy * par%latent_heat_water)
            end do
        end do
        !$omp end parallel do

    end subroutine update_face_fluxes

    ! ============================================================
    ! Flow routing: dispatcher + four implementations
    ! ============================================================
    ! All compute psi_out, the accumulated upstream outflow per cell [m3/s],
    ! seeded with mdot_total*dx*dy/rho_w. Mirrors route_psi_out!.
    !
    !   * TAPED (default) -- the recursion's traversal recorded once per
    !     calc_k24 call (it depends only on the mask and the routing
    !     weights) and replayed every sweep. Bit-identical to RECURSIVE for
    !     GDS_WARNER; the only algorithm for the other routing schemes.
    !   * RECURSIVE -- depth-first with memoization (GDS_WARNER).
    !   * ITERATIVE -- the same traversal with an explicit stack.
    !   * TOPOSORT -- Kahn's algorithm; exact only on an acyclic graph.
    !
    ! At a domain edge the out-of-range neighbour is skipped; in a periodic
    ! direction the neighbour index wraps instead.
    subroutine update_psi_out(wk, mask, dx, dy, par)
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par

        select case (par%flux_solver)
            case (K24_FLUX_TAPED)
                if (.not. wk%tape_valid) call record_routing_tape(wk, mask, dx, dy, par)
                call replay_routing_tape(wk, mask, dx, dy, par)
            case (K24_FLUX_RECURSIVE)
                call update_psi_out_recursive(wk, mask, dx, dy, par)
            case (K24_FLUX_ITERATIVE)
                call update_psi_out_iterative(wk, mask, dx, dy, par)
            case (K24_FLUX_TOPOSORT)
                call update_psi_out_toposort(wk, mask, dx, dy, par)
            case default
                write(*,*) "update_psi_out:: error: k24_flux_solver must be one of [0,1,2,3]."
                write(*,*) "flux_solver = ", par%flux_solver
                stop
        end select
    end subroutine update_psi_out

    ! --- Taped (record once, replay every sweep) ------------------
    subroutine record_routing_tape(wk, mask, dx, dy, par)
        ! record_routing_tape_kernel! / record_routing_tape_weights_kernel!:
        ! the explicit-stack traversal of update_psi_out_iterative, recording
        ! only the operations. Op k is psi(dst) += psi(src)*w, or, when
        ! src_i == 0, the clamp psi(dst) = max(0, psi(dst)). GDS_WARNER reads
        ! the routing directions directly (routing_weight from the
        ! neighbour's side); the other schemes read w8(opposite(d), n).
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par

        integer  :: i, j, si, sj, sk, ni, nj, top, capacity, call_count, ndirs, done, n_ground, cap_ops
        real(dp) :: w
        logical  :: gds, hit_cap, inside

        gds   = (par%routing_scheme == K24_ROUTE_GDS_WARNER)
        ndirs = par%n_dirs
        done  = ndirs + 1

        n_ground = count(mask == 1.0_dp)
        ! Each grounded cell folds in at most ndirs neighbours and is clamped once.
        cap_ops = (ndirs + 1) * n_ground + 1
        if (allocated(wk%t_dst_i)) deallocate(wk%t_dst_i, wk%t_dst_j, wk%t_src_i, wk%t_src_j, wk%t_w)
        allocate(wk%t_dst_i(cap_ops), wk%t_dst_j(cap_ops), wk%t_src_i(cap_ops), wk%t_src_j(cap_ops), wk%t_w(cap_ops))
        wk%tape_n = 0

        wk%visited = 0
        call_count = 0
        hit_cap    = .FALSE.
        capacity   = n_ground + 1
        allocate(wk%stack_i(capacity), wk%stack_j(capacity), wk%stack_k(capacity))

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
                        call_count = call_count + 1
                        if (call_count > par%max_psi_out_calls) then
                            hit_cap = .TRUE.
                            call push_tape_op(wk, si, sj, 0, 0, 0.0_dp)
                            top = top - 1
                        else
                            wk%stack_k(top) = 1
                        end if

                    else if (sk < done) then

                        ! Overwritten below if the neighbour must be resolved first.
                        wk%stack_k(top) = sk + 1

                        call neighbour(si, sj, sk, wk%nx, wk%ny, par, ni, nj, inside)
                        if (.not. inside) cycle

                        if (gds) then
                            w = routing_weight(wk%gsx(ni,nj), wk%gsy(ni,nj), -K24_DIRS(1,sk), -K24_DIRS(2,sk), dx, dy)
                        else
                            w = wk%w8(K24_OPPOSITE(sk), ni, nj)
                        end if
                        if (.not. (w > 0.0_dp .and. mask(ni,nj) == 1.0_dp)) cycle

                        if (wk%visited(ni,nj) == 1) then
                            call push_tape_op(wk, si, sj, ni, nj, w)
                        else
                            ! Resolve the neighbour first, then come back to this
                            ! same neighbour (sk unchanged).
                            wk%stack_k(top) = sk
                            if (top >= capacity) then
                                write(*,*) "record_routing_tape:: error: routing stack overflow."
                                stop
                            end if
                            top = top + 1
                            wk%stack_i(top) = ni; wk%stack_j(top) = nj; wk%stack_k(top) = 0
                        end if

                    else

                        call push_tape_op(wk, si, sj, 0, 0, 0.0_dp)
                        top = top - 1

                    end if

                end do

            end do
        end do

        deallocate(wk%stack_i, wk%stack_j, wk%stack_k)

        if (hit_cap) then
            write(*,'(a,i0,a)') " k24: WARNING the taped routing hit k24_max_psi_out_calls = ", &
                par%max_psi_out_calls, " cells in one sweep; cutting the flow routing off early." // &
                " Raise k24_max_psi_out_calls if this grid genuinely has more grounded cells."
        end if

        wk%tape_valid = .TRUE.

    end subroutine record_routing_tape

    subroutine push_tape_op(wk, ci, cj, ni, nj, w)
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        integer,  intent(IN) :: ci, cj, ni, nj
        real(dp), intent(IN) :: w
        wk%tape_n = wk%tape_n + 1
        wk%t_dst_i(wk%tape_n) = ci; wk%t_dst_j(wk%tape_n) = cj
        wk%t_src_i(wk%tape_n) = ni; wk%t_src_j(wk%tape_n) = nj
        wk%t_w(wk%tape_n)     = w
    end subroutine push_tape_op

    subroutine replay_routing_tape(wk, mask, dx, dy, par)
        ! replay_routing_tape!: every grounded cell's psi_out set to its own
        ! source term, then the recorded operations applied in order.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par
        integer :: i, j, k, ci, cj, ni

        ! DEVIATION 3: zeroed everywhere first.
        wk%psi_out = 0.0_dp
        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, wk%ny
            do i = 1, wk%nx
                if (mask(i,j) == 1.0_dp) wk%psi_out(i,j) = wk%mdot_total(i,j) * dx * dy / par%water_density
            end do
        end do
        !$omp end parallel do

        do k = 1, wk%tape_n
            ci = wk%t_dst_i(k); cj = wk%t_dst_j(k); ni = wk%t_src_i(k)
            if (ni == 0) then
                wk%psi_out(ci,cj) = max(0.0_dp, wk%psi_out(ci,cj))
            else
                wk%psi_out(ci,cj) = wk%psi_out(ci,cj) + wk%psi_out(ni, wk%t_src_j(k)) * wk%t_w(k)
            end if
        end do

    end subroutine replay_routing_tape

    ! --- Recursive (DFS + memoization) --------------------------
    subroutine update_psi_out_recursive(wk, mask, dx, dy, par)
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par

        integer  :: i, j, call_count
        logical  :: warned
        real(dp) :: dummy

        ! DEVIATION 3: psi_out is zeroed everywhere.
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
        ! number of cells one sweep may visit (KORI-ULB's funcnt <= 5e4); once
        ! tripped the current cell is treated as a terminal source, clamped at
        ! zero like the normal exit.
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
        logical  :: inside

        if (mask(i,j) /= 1.0_dp) then
            psi_value = 0.0_dp
            return
        end if

        if (wk%visited(i,j) == 1) then
            psi_value = wk%psi_out(i,j)
            return
        end if

        wk%visited(i,j) = 1
        wk%psi_out(i,j) = wk%mdot_total(i,j) * dx * dy / par%water_density

        call_count = call_count + 1
        if (call_count > par%max_psi_out_calls) then
            if (.not. warned) then
                write(*,'(a,i0,a)') " k24: WARNING accumulate_psi_out hit k24_max_psi_out_calls = ", &
                    par%max_psi_out_calls, " cells in one sweep; cutting the flow routing off early." // &
                    " Raise k24_max_psi_out_calls if this grid genuinely has more grounded cells."
                warned = .TRUE.
            end if
            wk%psi_out(i,j) = max(0.0_dp, wk%psi_out(i,j))
            psi_value = wk%psi_out(i,j)
            return
        end if

        do d = 1, 4
            call neighbour(i, j, d, wk%nx, wk%ny, par, ni, nj, inside)
            if (.not. inside) cycle

            w = routing_weight(wk%gsx(ni,nj), wk%gsy(ni,nj), -K24_DIRS(1,d), -K24_DIRS(2,d), dx, dy)

            if (w > 0.0_dp) then
                wk%psi_out(i,j) = wk%psi_out(i,j) + &
                    accumulate_psi_out(wk, ni, nj, mask, dx, dy, par, call_count, warned) * w
            end if
        end do

        ! If the source is negative enough that all the flux refreezes, floor at zero.
        wk%psi_out(i,j) = max(0.0_dp, wk%psi_out(i,j))
        psi_value = wk%psi_out(i,j)

    end function accumulate_psi_out

    ! --- Iterative (explicit stack, same traversal) --------------
    subroutine update_psi_out_iterative(wk, mask, dx, dy, par)
        ! Mirrors update_psi_out_iterative! (water_flux.jl). Stack entry
        ! (i, j, k): k == 0 unvisited; 1 <= k <= 4 neighbours 1..k-1 folded in,
        ! k next; k == 5 finalize. A not-yet-visited neighbour is pushed
        ! WITHOUT advancing k, so the parent frame is re-entered and takes the
        ! "already visited" branch.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par

        integer :: i, j, si, sj, sk, ni, nj, d, top, capacity, call_count
        real(dp) :: w
        logical  :: warned, inside

        wk%psi_out = 0.0_dp
        wk%visited = 0
        call_count = 0
        warned     = .FALSE.

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
                        wk%psi_out(si,sj) = wk%mdot_total(si,sj) * dx * dy / par%water_density

                        call_count = call_count + 1
                        if (call_count > par%max_psi_out_calls) then
                            if (.not. warned) then
                                write(*,'(a,i0,a)') " k24: WARNING update_psi_out_iterative hit k24_max_psi_out_calls = ", &
                                    par%max_psi_out_calls, " cells in one sweep; cutting the flow routing off early."
                                warned = .TRUE.
                            end if
                            ! Clamp before popping, matching the recursive cap-trip.
                            wk%psi_out(si,sj) = max(0.0_dp, wk%psi_out(si,sj))
                            top = top - 1
                        else
                            wk%stack_k(top) = 1
                        end if

                    else if (sk <= 4) then

                        d = sk
                        call neighbour(si, sj, d, wk%nx, wk%ny, par, ni, nj, inside)

                        if (.not. inside) then
                            wk%stack_k(top) = sk + 1
                            cycle
                        end if

                        w = routing_weight(wk%gsx(ni,nj), wk%gsy(ni,nj), -K24_DIRS(1,d), -K24_DIRS(2,d), dx, dy)

                        if (w <= 0.0_dp) then
                            wk%stack_k(top) = sk + 1
                            cycle
                        end if

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
        ! evaluated from the SOURCE cell's own direction, routing_weight(gs(A),
        ! d) for the edge A -> A+d.
        implicit none
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: mask(:,:), dx, dy
        type(k24_param_class),intent(IN)    :: par

        integer :: i, j, ni, nj, d, head, tail, capacity
        integer :: total_masked, processed, n_stuck
        real(dp) :: w
        logical  :: inside

        wk%psi_out = 0.0_dp

        allocate(wk%in_degree(wk%nx, wk%ny))
        wk%in_degree = 0

        total_masked = count(mask == 1.0_dp)

        !$omp parallel do default(shared) private(i,j,d,ni,nj,w,inside) schedule(static)
        do j = 1, wk%ny
            do i = 1, wk%nx
                if (mask(i,j) /= 1.0_dp) cycle
                do d = 1, 4
                    call neighbour(i, j, d, wk%nx, wk%ny, par, ni, nj, inside)
                    if (.not. inside) cycle
                    if (mask(ni,nj) /= 1.0_dp) cycle
                    ! Edge (ni,nj) -> (i,j) exists iff (ni,nj) flows toward us.
                    w = routing_weight(wk%gsx(ni,nj), wk%gsy(ni,nj), -K24_DIRS(1,d), -K24_DIRS(2,d), dx, dy)
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

            wk%psi_out(i,j) = max(0.0_dp, wk%psi_out(i,j) + wk%mdot_total(i,j) * dx * dy / par%water_density)

            do d = 1, 4
                call neighbour(i, j, d, wk%nx, wk%ny, par, ni, nj, inside)
                if (.not. inside) cycle
                if (mask(ni,nj) /= 1.0_dp) cycle

                w = routing_weight(wk%gsx(i,j), wk%gsy(i,j), K24_DIRS(1,d), K24_DIRS(2,d), dx, dy)

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
                write(*,*) "Expect this on real topography at any k24_coupling_length_kamb86."
                write(*,*) "Use k24_flux_solver = 3 (taped), 0 (recursive) or 1 (iterative), or set k24_toposort_allow_cycles = .TRUE."
                stop
            end if
        end if

        deallocate(wk%in_degree, wk%stack_i, wk%stack_j)

    end subroutine update_psi_out_toposort

    ! ============================================================
    ! Effective pressure
    ! ============================================================
    subroutine update_N(N, q, wk, mask, uxy_b, A_glen, kappa, H_ice, par)
        ! Mirrors update_N! (effective_pressure.jl), the fused grounded-only
        ! pass: Po = rho_i*g*H everywhere; on grounded cells Q, S_inf, H,
        ! N_inf and N, each evaluated in the same order as Julia; every other
        ! cell gets N = 0 and zero conduit fields. abs_g and phi0 are those of
        ! the TRUE potential.
        implicit none
        real(dp),             intent(INOUT) :: N(:,:)
        real(dp),             intent(IN)    :: q(:,:), mask(:,:)
        type(k24_work_class), intent(INOUT) :: wk
        real(dp),             intent(IN)    :: uxy_b(:,:), A_glen(:,:), kappa(:,:), H_ice(:,:)
        type(k24_param_class),intent(IN)    :: par

        integer  :: i, j, nx, ny
        real(dp) :: K_fac, grad_exp, Q_exp, denom_const, inv_n, sqrt_pi
        real(dp) :: sliding_coeff, melt_coeff, Q_c
        real(dp) :: Po_ij, Q_ij, S_ij, Hh, Hs, H_ij, Ninf_ij, expo, k

        nx = wk%nx; ny = wk%ny

        K_fac       = par%K**(-1.0_dp / par%manning_coefficient_exponent)
        grad_exp    = (1.0_dp - par%bed_friction_exponent) / par%manning_coefficient_exponent
        Q_exp       = 1.0_dp / par%manning_coefficient_exponent
        denom_const = 2.0_dp * par%manning_exponent**(-par%manning_exponent) * par%ice_density * par%latent_heat_water
        inv_n       = 1.0_dp / par%manning_exponent
        sqrt_pi     = sqrt(K24_PI)
        Q_c         = par%critical_discharge

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

        !$omp parallel do default(shared) private(i,j,Po_ij,Q_ij,S_ij,Hh,Hs,H_ij,Ninf_ij,expo,k) schedule(static)
        do j = 1, ny
            do i = 1, nx

                ! From the RAW ice thickness, not the potential-filled h.
                Po_ij = par%ice_density * par%gravity * H_ice(i,j)
                wk%Po(i,j) = Po_ij

                if (mask(i,j) /= 1.0_dp) then
                    wk%Q(i,j) = 0.0_dp; wk%S_inf(i,j) = 0.0_dp
                    wk%H_hard(i,j) = 0.0_dp; wk%H_soft(i,j) = 0.0_dp; wk%H_cond(i,j) = 0.0_dp
                    wk%N_inf(i,j) = 0.0_dp; N(i,j) = 0.0_dp
                    cycle
                end if

                Q_ij = q(i,j) * par%coupling_length

                ! DEVIATION 1: zero flux means zero conduit cross-section.
                if (Q_ij == 0.0_dp) then
                    S_ij = 0.0_dp
                else
                    S_ij = K_fac * wk%abs_g(i,j)**grad_exp * Q_ij**Q_exp
                end if

                ! Q_c per drainage mode: the configured value (BOTH), the
                ! exp(-Q/Q_c) -> 0 limit (EFFICIENT, Q_c -> 0) or -> 1 limit
                ! (INEFFICIENT, Q_c -> Inf), taken analytically. The
                ! Q_c == 0, Q == 0 case is the 0/0 limit, H_soft = 0.
                Hh = sqrt(S_ij)
                select case (par%drainage_mode)
                    case (K24_DRAINAGE_INEFFICIENT)
                        expo = 1.0_dp
                    case (K24_DRAINAGE_EFFICIENT)
                        expo = 0.0_dp
                    case default
                        if (Q_c == 0.0_dp) then
                            expo = 0.0_dp
                        else
                            expo = exp(-Q_ij / Q_c)
                        end if
                end select
                Hs = max(0.0_dp, par%initial_cavity_height &
                     + (sqrt(S_ij) / par%till_factor - par%initial_cavity_height) * expo)
                if (Q_ij == 0.0_dp .and. (par%drainage_mode == K24_DRAINAGE_EFFICIENT .or. &
                    (par%drainage_mode == K24_DRAINAGE_BOTH .and. Q_c == 0.0_dp))) then
                    Hs = 0.0_dp
                end if
                k = kappa(i,j)
                H_ij = (1.0_dp - k) * Hh + k * Hs

                ! DEVIATION 1: S_inf == 0 gives N_inf = Po.
                if (S_ij == 0.0_dp) then
                    Ninf_ij = Po_ij
                else
                    Ninf_ij = min(max( &
                        ((H_ij * H_ij) / (S_ij * S_ij) * (sliding_coeff * par%ice_density * par%latent_heat_water &
                            * uxy_b(i,j) * par%bed_thickness + melt_coeff * Q_ij * wk%abs_g(i,j)) &
                         / (denom_const * A_glen(i,j)))**inv_n, &
                        par%min_pressure_fraction * Po_ij), Po_ij)
                end if

                ! DEVIATION 1: N_inf == 0 gives N = 0.
                if (Ninf_ij == 0.0_dp) then
                    N(i,j) = 0.0_dp
                else
                    N(i,j) = max(0.0_dp, erf(sqrt_pi * wk%phi0(i,j) / (2.0_dp * Ninf_ij)) * Ninf_ij)
                end if

                wk%Q(i,j) = Q_ij; wk%S_inf(i,j) = S_ij; wk%H_hard(i,j) = Hh; wk%H_soft(i,j) = Hs
                wk%H_cond(i,j) = H_ij; wk%N_inf(i,j) = Ninf_ij

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
                !$omp parallel do default(shared) private(i,j) schedule(static)
                do j = 1, ny
                    do i = 1, nx
                        W(i,j) = wk%S_inf(i,j) / par%coupling_length
                    end do
                end do
                !$omp end parallel do

            case (K24_WTHICK_DARCY_WEISBACH)
                ! d = (f*rho_w*q^2 / (4*|grad phi0|))^(1/3), true potential.
                if (par%gradient_convention == K24_GRAD_MEAN) then
                    g_mean = masked_mean(wk%abs_g, mask)
                    !$omp parallel do default(shared) private(i,j,num) schedule(static)
                    do j = 1, ny
                        do i = 1, nx
                            num = par%friction_factor * par%water_density * q(i,j) * q(i,j) &
                                / (4.0_dp * g_mean + epsilon(1.0_dp))
                            W(i,j) = min(par%W_max, max(par%W_min, num**third))
                        end do
                    end do
                    !$omp end parallel do
                else
                    !$omp parallel do default(shared) private(i,j,num) schedule(static)
                    do j = 1, ny
                        do i = 1, nx
                            num = par%friction_factor * par%water_density * q(i,j) * q(i,j) &
                                / (4.0_dp * wk%abs_g(i,j) + epsilon(1.0_dp))
                            W(i,j) = min(par%W_max, max(par%W_min, num**third))
                        end do
                    end do
                    !$omp end parallel do
                end if

            case (K24_WTHICK_LAMINAR)
                ! d = (12*eta_w*q / |grad phi0_s|)^(1/3), routing directions.
                if (par%gradient_convention == K24_GRAD_MEAN) then
                    ! Julia adds no eps guard in this branch; the max() below
                    ! keeps a zero-mean domain finite.
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
                            num = 12.0_dp * par%eta_w * q(i,j) / (wk%abs_gs(i,j) + epsilon(1.0_dp))
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
