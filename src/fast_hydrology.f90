module fast_hydrology
    ! Top-level wrapper for basal-hydrology models. Splits the model into
    ! two orthogonal components run sequentially each step:
    !
    !   method_til       : how the till storage W_til is evolved
    !       TIL_BUCKET   : local mass-balance bucket (BvP15-style till storage)
    !                      W_attempt = W_til + dt*(mdot - till_rate)
    !                      W_til     = clamp(W_attempt, 0, W_til_max)
    !                      overflow  = max(0, (W_attempt - W_til_max)/dt)
    !                      With W_til_max(i,j) = 0 the bucket holds nothing
    !                      and everything passes through as overflow.
    !       TIL_EXTERNAL : host owns W_til; bucket does not touch it.
    !                      Overflow = mdot (raw pass-through to transport).
    !
    !   method_transport : how the distributed sheet W is computed
    !       TRANSPORT_NONE : W = 0, q_x = q_y = 0. N comes from the bucket's
    !                        N_closure (overburden / marine / till / ...).
    !       TRANSPORT_K24  : Kazmierczak 2024 distributed model with
    !                        FFT-smoothed gradients. Writes W, q_x, q_y, N,
    !                        p_w. The K24 source is the bucket overflow
    !                        (TIL_BUCKET) or raw mdot (TIL_EXTERNAL).
    !
    ! Notation follows van Pelt & Bueler 2015 (BvP15):
    !   W_til : till water storage thickness   [m]
    !   W     : distributed sheet thickness    [m]
    !   mdot  : source rate from ice base      [m/s, water-equivalent]
    !
    ! Time / unit convention on the public API:
    !   hydro_update(... time, ...) : time in YEARS (matches host model)
    !   mdot, uxy_b                 : SI [m/s]
    !   par%dx, par%dy              : [m]
    ! Internally everything is SI. dt_year (time - hyd%now%time) is computed
    ! in years for state bookkeeping; dt_sec = dt_year * SEC_PER_YEAR is the
    ! denominator used by the bucket update and dW_til_dt / overflow rates.
    !
    ! In all cases, hyd%now%dW_til_dt is computed as the natural-sign change
    ! in W_til over the step: (W_til_new - W_til_old) / dt_sec  [m/s].

    use nml
    use fast_hydrology_k24
    use fast_hydrology_bucket
    use fast_hydrology_closures

    implicit none

    integer, parameter :: dp = kind(1.d0)
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: wp = sp

    ! Seconds per year used to convert the API's time argument (years) to
    ! the internal dt_sec. Match the value in bucket.f90.
    real(wp), parameter, public :: SEC_PER_YEAR = 3.1556926e7_wp

    ! ---------- method_til enum (par%method_til) ----------
    integer, parameter, public :: TIL_BUCKET   = 0    ! default
    integer, parameter, public :: TIL_EXTERNAL = 1    ! host owns W_til

    ! ---------- method_transport enum (par%method_transport) ----------
    integer, parameter, public :: TRANSPORT_NONE = 0
    integer, parameter, public :: TRANSPORT_K24  = 1

    ! ---------- Init-method enum (par%init_method) ----------
    integer, parameter, public :: HYDRO_INIT_ZERO     = 0
    integer, parameter, public :: HYDRO_INIT_EXTERNAL = 1

    type hydro_param_class
        integer  :: method_til
        integer  :: method_transport
        integer  :: init_method
        real(wp) :: dx
        real(wp) :: dy
        real(wp) :: W_til_max         ! [m]  scalar default for hyd%now%W_til_max
        integer  :: mask_bc           ! see bucket::MASK_BC_*
        real(wp) :: W_til_bc          ! [m]  imposed W_til at the domain border (MASK_BC_IMPOSED)
        type(bucket_param_class)  :: bucket
        type(k24_param_class)     :: k24
        type(closure_param_class) :: closures
    end type

    type hydro_state_class
        real(wp) :: time
        real(wp) :: dt
        logical  :: initialized
        real(wp), allocatable :: W_til(:,:)        ! [m]      till storage (bucket)
        real(wp), allocatable :: W_til_max(:,:)    ! [m]      per-cell cap (defaults to par%W_til_max)
        real(wp), allocatable :: dW_til_dt(:,:)    ! [m/s]    natural-sign change in W_til over the step
        real(wp), allocatable :: overflow(:,:)     ! [m/s]    till-saturation spill -> K24 source
        real(wp), allocatable :: W(:,:)            ! [m]      K24 distributed sheet thickness; 0 when TRANSPORT_NONE
        real(wp), allocatable :: p_w(:,:)          ! [Pa]
        real(wp), allocatable :: q_x(:,:)          ! [m2/s]
        real(wp), allocatable :: q_y(:,:)          ! [m2/s]
        real(wp), allocatable :: N(:,:)            ! [Pa]
        real(wp), allocatable :: kappa(:,:)
    end type

    type hydro_class
        type(hydro_param_class) :: par
        type(hydro_state_class) :: now
    end type

    private
    public :: hydro_class
    public :: hydro_init
    public :: hydro_init_state
    public :: hydro_update
    public :: wp

contains

    subroutine hydro_init(hyd, filename, nx, ny, group)

        implicit none

        type(hydro_class), intent(INOUT)        :: hyd
        character(len=*),  intent(IN)           :: filename
        integer,           intent(IN)           :: nx, ny
        character(len=*),  intent(IN), optional :: group

        character(len=32) :: nml_group

        ! Resolve the namelist group name (default = "fhyd").
        if (present(group)) then
            nml_group = trim(group)
        else
            nml_group = "fhyd"
        end if

        call hydro_par_load(hyd%par, filename, nml_group)

        call hydro_allocate(hyd%now, nx, ny)

        hyd%now%time        = 0.0_wp
        hyd%now%dt          = 0.0_wp
        hyd%now%initialized = .false.

        return

    end subroutine hydro_init

    subroutine hydro_init_state(hyd, z_bed, f_ice, f_grnd, time)
        ! Initialize state for the first update. Kappa is filled only when
        ! method_transport == K24. W_til is set according to par%init_method
        ! on grounded cells; the floating + adjacent-to-floating override is
        ! applied for TIL_BUCKET. Per-cell W_til_max is filled from
        ! par%W_til_max (host can overwrite hyd%now%W_til_max(:,:) between
        ! init and the first update if it wants a spatial cap). Diagnostic
        ! fields are zeroed.

        implicit none

        type(hydro_class), intent(INOUT) :: hyd
        real(wp),          intent(IN)    :: z_bed(:,:)
        real(wp),          intent(IN)    :: f_ice(:,:)
        real(wp),          intent(IN)    :: f_grnd(:,:)
        real(wp),          intent(IN)    :: time

        real(dp), allocatable :: z_bed_dp(:,:), kappa_dp(:,:)
        integer :: nx, ny

        nx = size(z_bed,1)
        ny = size(z_bed,2)

        if (hyd%par%method_transport == TRANSPORT_K24) then
            allocate(z_bed_dp(nx,ny), kappa_dp(nx,ny))
            z_bed_dp = real(z_bed, dp)
            call initialize_kappa(kappa_dp, z_bed_dp, hyd%par%k24%substrate_type)
            hyd%now%kappa = real(kappa_dp, wp)
            deallocate(z_bed_dp, kappa_dp)
        else
            hyd%now%kappa = 0.0_wp
        end if

        ! Per-cell W_til_max defaults to the scalar parameter. Host can
        ! overwrite the array between init and the first update if a
        ! spatial cap is needed.
        hyd%now%W_til_max = hyd%par%W_til_max

        select case (hyd%par%init_method)
            case (HYDRO_INIT_ZERO)
                hyd%now%W_til = 0.0_wp
            case (HYDRO_INIT_EXTERNAL)
                ! Leave W_til untouched (host filled it)
                continue
            case default
                write(*,*) "hydro_init_state:: error: init_method must be 0 or 1."
                write(*,*) "init_method = ", hyd%par%init_method
                stop
        end select

        ! Zero W_til outside the active (grounded ice) set and apply the
        ! configurable domain-border BC. Skipped for TIL_EXTERNAL (host
        ! manages W_til including overrides).
        if (hyd%par%method_til == TIL_BUCKET) then
            call apply_floating_override(hyd%now%W_til, f_ice, f_grnd)
            call apply_mask_bc(hyd%now%W_til, hyd%par%mask_bc, hyd%par%W_til_bc)
        end if

        hyd%now%dW_til_dt = 0.0_wp
        hyd%now%overflow  = 0.0_wp
        hyd%now%W         = 0.0_wp
        hyd%now%p_w       = 0.0_wp
        hyd%now%q_x       = 0.0_wp
        hyd%now%q_y       = 0.0_wp
        hyd%now%N         = 0.0_wp

        hyd%now%time        = time
        hyd%now%dt          = 0.0_wp
        hyd%now%initialized = .true.

        return

    end subroutine hydro_init_state

    subroutine hydro_update(hyd, H_ice, z_bed, z_sl, f_ice, f_grnd, mask, &
                            mdot, uxy_b, A_glen, time)
        ! Advance one step. Flow:
        !   1. Till step (method_til):
        !        BUCKET   -> updates W_til; produces hyd%now%overflow
        !        EXTERNAL -> leaves W_til alone; overflow = mdot
        !   2. Transport step (method_transport):
        !        NONE -> W, q_x, q_y zeroed
        !        K24  -> calc_k24 with source = hyd%now%overflow,
        !                writes W, q_x, q_y, N, p_w
        !   3. N closure:
        !        K24 active -> N already set by K24 (skip)
        !        otherwise  -> apply par%bucket%N_closure on the current W_til
        !   4. W_til overrides (BUCKET only): floating + mask-BC.
        !   5. dW_til_dt = (W_til_new - W_til_old) / dt.

        implicit none

        type(hydro_class), intent(INOUT) :: hyd
        real(wp),          intent(IN)    :: H_ice(:,:)
        real(wp),          intent(IN)    :: z_bed(:,:)
        real(wp),          intent(IN)    :: z_sl(:,:)
        real(wp),          intent(IN)    :: f_ice(:,:)
        real(wp),          intent(IN)    :: f_grnd(:,:)
        real(wp),          intent(IN)    :: mask(:,:)
        real(wp),          intent(IN)    :: mdot(:,:)
        real(wp),          intent(IN)    :: uxy_b(:,:)
        real(wp),          intent(IN)    :: A_glen(:,:)
        real(wp),          intent(IN)    :: time

        ! dp scratch arrays for the K24 boundary
        real(dp), allocatable :: H_ice_dp(:,:), z_bed_dp(:,:), mask_dp(:,:)
        real(dp), allocatable :: src_dp(:,:), uxy_b_dp(:,:), A_glen_dp(:,:)
        real(dp), allocatable :: kappa_dp(:,:)
        real(dp), allocatable :: q_x_dp(:,:), q_y_dp(:,:), N_dp(:,:), p_w_dp(:,:), W_dp(:,:)

        real(wp), allocatable :: W_til_old(:,:)
        real(wp) :: dt_year, dt_sec
        integer  :: nx, ny

        if (.not. hyd%now%initialized) then
            write(*,*) "hydro_update:: error: hydro_init_state must be called before hydro_update."
            stop
        end if

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(W_til_old(nx,ny))
        W_til_old = hyd%now%W_til

        ! API time is in years; convert step to seconds for internal SI.
        dt_year      = max(time - hyd%now%time, 0.0_wp)
        dt_sec       = dt_year * SEC_PER_YEAR
        hyd%now%time = time
        hyd%now%dt   = dt_sec

        ! ---- Step 1: till update (W_til + overflow) ----
        if (dt_sec > 0.0_wp) then
            select case (hyd%par%method_til)

                case (TIL_BUCKET)
                    call calc_bucket(hyd%now%W_til, hyd%now%overflow, &
                                     f_ice, f_grnd, mask, mdot, &
                                     dt_sec, hyd%par%bucket, hyd%now%W_til_max)

                case (TIL_EXTERNAL)
                    ! Host owns W_til; do not touch it. K24 sees raw source.
                    hyd%now%overflow = mdot

                case default
                    write(*,*) "hydro_update:: error: method_til must be 0 (BUCKET) or 1 (EXTERNAL)."
                    write(*,*) "method_til = ", hyd%par%method_til
                    stop

            end select
        else
            ! dt = 0: no step. Overflow is zero by convention.
            hyd%now%overflow = 0.0_wp
        end if

        ! ---- Step 2: transport (W, q_x, q_y) ----
        select case (hyd%par%method_transport)

            case (TRANSPORT_NONE)
                hyd%now%W   = 0.0_wp
                hyd%now%q_x = 0.0_wp
                hyd%now%q_y = 0.0_wp

            case (TRANSPORT_K24)
                if (dt_sec > 0.0_wp) then
                    allocate(H_ice_dp(nx,ny), z_bed_dp(nx,ny), mask_dp(nx,ny))
                    allocate(src_dp(nx,ny), uxy_b_dp(nx,ny), A_glen_dp(nx,ny))
                    allocate(kappa_dp(nx,ny))
                    allocate(q_x_dp(nx,ny), q_y_dp(nx,ny), N_dp(nx,ny), p_w_dp(nx,ny), W_dp(nx,ny))

                    H_ice_dp  = real(H_ice,           dp)
                    z_bed_dp  = real(z_bed,           dp)
                    mask_dp   = real(mask,            dp)
                    src_dp    = real(hyd%now%overflow,dp)
                    uxy_b_dp  = real(uxy_b,           dp)
                    A_glen_dp = real(A_glen,          dp)
                    kappa_dp  = real(hyd%now%kappa,   dp)

                    call calc_k24(q_x_dp, q_y_dp, N_dp, p_w_dp, W_dp, &
                                  H_ice_dp, z_bed_dp, mask_dp, src_dp, uxy_b_dp, A_glen_dp, &
                                  kappa_dp, &
                                  real(hyd%par%dx, dp), real(hyd%par%dy, dp), hyd%par%k24)

                    hyd%now%q_x = real(q_x_dp, wp)
                    hyd%now%q_y = real(q_y_dp, wp)
                    hyd%now%N   = real(N_dp,   wp)
                    hyd%now%p_w = real(p_w_dp, wp)
                    hyd%now%W   = real(W_dp,   wp)

                    deallocate(H_ice_dp, z_bed_dp, mask_dp)
                    deallocate(src_dp, uxy_b_dp, A_glen_dp)
                    deallocate(kappa_dp)
                    deallocate(q_x_dp, q_y_dp, N_dp, p_w_dp, W_dp)
                end if

            case default
                write(*,*) "hydro_update:: error: method_transport must be 0 (NONE) or 1 (K24)."
                write(*,*) "method_transport = ", hyd%par%method_transport
                stop

        end select

        ! ---- Step 3: N closure ----
        ! K24 produces N on its own; in all other cases run the configured
        ! N-closure on the current W_til (and geometry).
        if (hyd%par%method_transport /= TRANSPORT_K24) then
            call apply_N_closure(hyd, H_ice, z_bed, z_sl, f_ice, f_grnd)
        end if

        ! ---- Step 4: W_til overrides ----
        if (hyd%par%method_til == TIL_BUCKET .and. dt_sec > 0.0_wp) then
            call apply_floating_override(hyd%now%W_til, f_ice, f_grnd)
            call apply_mask_bc(hyd%now%W_til, hyd%par%mask_bc, hyd%par%W_til_bc)
        end if

        ! ---- Step 5: dW_til_dt ----
        if (dt_sec > 0.0_wp) then
            hyd%now%dW_til_dt = (hyd%now%W_til - W_til_old) / dt_sec
        else
            hyd%now%dW_til_dt = 0.0_wp
        end if

        deallocate(W_til_old)

        return

    end subroutine hydro_update

    ! ------------------------------------------------------------
    ! N-closure post-step. Writes hyd%now%N and derives p_w = Po - N.
    ! Called when method_transport /= K24. No-op when N_closure == NONE.
    ! ------------------------------------------------------------
    subroutine apply_N_closure(hyd, H_ice, z_bed, z_sl, f_ice, f_grnd)

        implicit none

        type(hydro_class), intent(INOUT) :: hyd
        real(wp),          intent(IN)    :: H_ice(:,:), z_bed(:,:), z_sl(:,:)
        real(wp),          intent(IN)    :: f_ice(:,:), f_grnd(:,:)

        integer  :: i, j, nx, ny
        real(wp) :: Po, H_eff

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        select case (hyd%par%bucket%N_closure)

            case (N_CLOSURE_NONE)
                return

            case (N_CLOSURE_OVERBURDEN)
                !$omp parallel do default(shared) private(i,j) schedule(static)
                do j = 1, ny
                do i = 1, nx
                    call hydro_calc_N_overburden(hyd%now%N(i,j), H_ice(i,j), f_ice(i,j), f_grnd(i,j), &
                                                 hyd%par%closures%rho_ice, hyd%par%closures%g)
                end do
                end do
                !$omp end parallel do

            case (N_CLOSURE_MARINE)
                !$omp parallel do default(shared) private(i,j) schedule(static)
                do j = 1, ny
                do i = 1, nx
                    call hydro_calc_N_marine(hyd%now%N(i,j), H_ice(i,j), f_ice(i,j), z_bed(i,j), z_sl(i,j), &
                                             hyd%par%closures%marine%p, hyd%par%closures%rho_ice, &
                                             hyd%par%closures%marine%rho_sw, hyd%par%closures%g)
                end do
                end do
                !$omp end parallel do

            case (N_CLOSURE_TILL)
                !$omp parallel do default(shared) private(i,j) schedule(static)
                do j = 1, ny
                do i = 1, nx
                    call hydro_calc_N_till(hyd%now%N(i,j), hyd%now%W_til(i,j), H_ice(i,j), &
                                           f_ice(i,j), f_grnd(i,j), hyd%now%W_til_max(i,j), &
                                           hyd%par%closures%till%N0, hyd%par%closures%till%delta, &
                                           hyd%par%closures%till%e0, hyd%par%closures%till%Cc, &
                                           hyd%par%closures%rho_ice, hyd%par%closures%g)
                end do
                end do
                !$omp end parallel do

            case default
                write(*,*) "apply_N_closure:: error: unsupported post-step closure ", &
                           hyd%par%bucket%N_closure
                write(*,*) "(N_CLOSURE_TWO_VALUE is standalone-only.)"
                stop

        end select

        ! Derive p_w = Po - N on grounded cells
        !$omp parallel do default(shared) private(i,j,Po,H_eff) schedule(static)
        do j = 1, ny
        do i = 1, nx
            if (f_grnd(i,j) > 0.0_wp .and. f_ice(i,j) > 0.0_wp) then
                if (f_ice(i,j) > 0.0_wp) then
                    H_eff = H_ice(i,j) / f_ice(i,j)
                else
                    H_eff = H_ice(i,j)
                end if
                if (f_ice(i,j) < 1.0_wp) H_eff = 0.0_wp
                Po           = hyd%par%closures%rho_ice * hyd%par%closures%g * H_eff
                hyd%now%p_w(i,j) = Po - hyd%now%N(i,j)
            else
                hyd%now%p_w(i,j) = 0.0_wp
            end if
        end do
        end do
        !$omp end parallel do

        return

    end subroutine apply_N_closure

    subroutine hydro_par_load(par, filename, group, init)

        implicit none

        type(hydro_param_class), intent(INOUT) :: par
        character(len=*),        intent(IN)    :: filename
        character(len=*),        intent(IN)    :: group
        logical, optional,       intent(IN)    :: init

        logical :: init_pars

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        par%method_til       = TIL_BUCKET
        par%method_transport = TRANSPORT_NONE
        par%init_method      = HYDRO_INIT_ZERO
        par%dx               = 0.0_wp
        par%dy               = 0.0_wp
        par%W_til_max        = 2.0_wp
        par%mask_bc          = MASK_BC_ZERO
        par%W_til_bc         = 0.0_wp

        call nml_read(filename,group,"method_til",        par%method_til,        init=init_pars)
        call nml_read(filename,group,"method_transport",  par%method_transport,  init=init_pars)
        call nml_read(filename,group,"init_method",       par%init_method,       init=init_pars)
        call nml_read(filename,group,"dx",                par%dx,                init=init_pars)
        call nml_read(filename,group,"dy",                par%dy,                init=init_pars)
        call nml_read(filename,group,"W_til_max",         par%W_til_max,         init=init_pars)
        call nml_read(filename,group,"mask_bc",           par%mask_bc,           init=init_pars)
        call nml_read(filename,group,"W_til_bc",          par%W_til_bc,          init=init_pars)

        call bucket_par_load (par%bucket,   filename, group, init=init_pars)
        call k24_par_load    (par%k24,      filename, group, init=init_pars)
        call closure_par_load(par%closures, filename, group, init=init_pars)

        return

    end subroutine hydro_par_load

    subroutine hydro_allocate(now, nx, ny)

        implicit none

        type(hydro_state_class), intent(INOUT) :: now
        integer,                 intent(IN)    :: nx, ny

        call hydro_deallocate(now)

        allocate(now%W_til(nx,ny))
        allocate(now%W_til_max(nx,ny))
        allocate(now%dW_til_dt(nx,ny))
        allocate(now%overflow(nx,ny))
        allocate(now%W(nx,ny))
        allocate(now%p_w(nx,ny))
        allocate(now%q_x(nx,ny))
        allocate(now%q_y(nx,ny))
        allocate(now%N(nx,ny))
        allocate(now%kappa(nx,ny))

        now%W_til     = 0.0_wp
        now%W_til_max = 0.0_wp
        now%dW_til_dt = 0.0_wp
        now%overflow  = 0.0_wp
        now%W         = 0.0_wp
        now%p_w       = 0.0_wp
        now%q_x       = 0.0_wp
        now%q_y       = 0.0_wp
        now%N         = 0.0_wp
        now%kappa     = 0.0_wp

        return

    end subroutine hydro_allocate

    subroutine hydro_deallocate(now)

        implicit none

        type(hydro_state_class), intent(INOUT) :: now

        if (allocated(now%W_til))     deallocate(now%W_til)
        if (allocated(now%W_til_max)) deallocate(now%W_til_max)
        if (allocated(now%dW_til_dt)) deallocate(now%dW_til_dt)
        if (allocated(now%overflow))  deallocate(now%overflow)
        if (allocated(now%W))         deallocate(now%W)
        if (allocated(now%p_w))       deallocate(now%p_w)
        if (allocated(now%q_x))       deallocate(now%q_x)
        if (allocated(now%q_y))       deallocate(now%q_y)
        if (allocated(now%N))         deallocate(now%N)
        if (allocated(now%kappa))     deallocate(now%kappa)

        return

    end subroutine hydro_deallocate

end module fast_hydrology
