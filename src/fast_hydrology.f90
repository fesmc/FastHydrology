module fast_hydrology
    ! Top-level wrapper for basal-hydrology models. Dispatches to one of:
    !   NONE   - no-op (state frozen)
    !   BUCKET - local mass-balance bucket (yelmo-equivalent)
    !   K24    - steady-state distributed model with FFT-smoothed gradients
    !
    ! Method-specific contracts:
    !   NONE   : writes nothing; H_w externally managed.
    !   BUCKET : writes H_w (with floating/grounded-ice-free overrides);
    !            optionally writes N, p_w via par%bucket%N_closure.
    !   K24    : writes q_x, q_y, N, p_w always; writes H_w only if
    !            par%k24%update_H_w (instant if tau_H=0, else relaxed).
    !
    ! In all cases, hyd%now%dHwdt is computed as the natural-sign change
    ! over the step: (H_w_new - H_w_old) / dt.

    use nml
    use fast_hydrology_k24
    use fast_hydrology_bucket
    use fast_hydrology_closures

    implicit none

    integer, parameter :: dp = kind(1.d0)
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: wp = sp

    ! ---------- Method enum (par%method) ----------
    integer, parameter, public :: HYDRO_METHOD_NONE   = 0
    integer, parameter, public :: HYDRO_METHOD_BUCKET = 1
    integer, parameter, public :: HYDRO_METHOD_K24    = 2

    ! ---------- Init-method enum (par%init_method) ----------
    integer, parameter, public :: HYDRO_INIT_ZERO     = 0
    integer, parameter, public :: HYDRO_INIT_EXTERNAL = 1

    type hydro_param_class
        integer  :: method
        integer  :: init_method
        real(wp) :: dx
        real(wp) :: dy
        real(wp) :: H_w_max
        type(bucket_param_class)  :: bucket
        type(k24_param_class)     :: k24
        type(closure_param_class) :: closures
        real(wp) :: tau_H              ! [a] K24 H_w relaxation timescale (0 = instant)
        logical  :: k24_update_H_w     ! K24 writes H_w if true
    end type

    type hydro_state_class
        real(wp) :: time
        real(wp) :: dt
        logical  :: initialized
        real(wp), allocatable :: H_w(:,:)
        real(wp), allocatable :: dHwdt(:,:)
        real(wp), allocatable :: p_w(:,:)
        real(wp), allocatable :: q_x(:,:)
        real(wp), allocatable :: q_y(:,:)
        real(wp), allocatable :: N(:,:)
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

    subroutine hydro_init(hyd, filename, nx, ny)

        implicit none

        type(hydro_class), intent(INOUT) :: hyd
        character(len=*),  intent(IN)    :: filename
        integer,           intent(IN)    :: nx, ny

        call hydro_par_load(hyd%par, filename)

        if (hyd%par%dx /= hyd%par%dy) then
            write(*,*) "hydro_init:: error: K24 requires dx == dy."
            write(*,*) "dx = ", hyd%par%dx, "  dy = ", hyd%par%dy
            stop
        end if

        call hydro_allocate(hyd%now, nx, ny)

        hyd%now%time        = 0.0_wp
        hyd%now%dt          = 0.0_wp
        hyd%now%initialized = .false.

        return

    end subroutine hydro_init

    subroutine hydro_init_state(hyd, z_bed, f_ice, f_grnd, time)
        ! Initialize state for the first update. Kappa is filled only when
        ! method = K24. H_w is set according to par%init_method on grounded
        ! cells; the floating + adjacent-to-floating override is always
        ! applied. Diagnostic fields are zeroed.

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

        if (hyd%par%method == HYDRO_METHOD_K24) then
            allocate(z_bed_dp(nx,ny), kappa_dp(nx,ny))
            z_bed_dp = real(z_bed, dp)
            call initialize_kappa(kappa_dp, z_bed_dp, hyd%par%k24%substrate_type)
            hyd%now%kappa = real(kappa_dp, wp)
            deallocate(z_bed_dp, kappa_dp)
        else
            hyd%now%kappa = 0.0_wp
        end if

        select case (hyd%par%init_method)
            case (HYDRO_INIT_ZERO)
                hyd%now%H_w = 0.0_wp
            case (HYDRO_INIT_EXTERNAL)
                ! Leave H_w untouched on grounded cells (host filled it)
                continue
            case default
                write(*,*) "hydro_init_state:: error: init_method must be 0 or 1."
                write(*,*) "init_method = ", hyd%par%init_method
                stop
        end select

        ! Always-applied floating + adjacent override
        call apply_floating_override(hyd%now%H_w, f_ice, f_grnd, hyd%par%H_w_max)

        hyd%now%dHwdt = 0.0_wp
        hyd%now%p_w   = 0.0_wp
        hyd%now%q_x   = 0.0_wp
        hyd%now%q_y   = 0.0_wp
        hyd%now%N     = 0.0_wp

        hyd%now%time        = time
        hyd%now%dt          = 0.0_wp
        hyd%now%initialized = .true.

        return

    end subroutine hydro_init_state

    subroutine hydro_update(hyd, H_ice, z_bed, z_sl, f_ice, f_grnd, mask, &
                            bmb_w, uxy_b, A_glen, time)

        implicit none

        type(hydro_class), intent(INOUT) :: hyd
        real(wp),          intent(IN)    :: H_ice(:,:)
        real(wp),          intent(IN)    :: z_bed(:,:)
        real(wp),          intent(IN)    :: z_sl(:,:)
        real(wp),          intent(IN)    :: f_ice(:,:)
        real(wp),          intent(IN)    :: f_grnd(:,:)
        real(wp),          intent(IN)    :: mask(:,:)
        real(wp),          intent(IN)    :: bmb_w(:,:)
        real(wp),          intent(IN)    :: uxy_b(:,:)
        real(wp),          intent(IN)    :: A_glen(:,:)
        real(wp),          intent(IN)    :: time

        ! dp scratch arrays for the K24 boundary
        real(dp), allocatable :: H_ice_dp(:,:), z_bed_dp(:,:), mask_dp(:,:)
        real(dp), allocatable :: bmb_w_dp(:,:), uxy_b_dp(:,:), A_glen_dp(:,:)
        real(dp), allocatable :: kappa_dp(:,:)
        real(dp), allocatable :: q_x_dp(:,:), q_y_dp(:,:), N_dp(:,:), p_w_dp(:,:), H_w_eq_dp(:,:)

        real(wp), allocatable :: H_w_old(:,:), H_w_eq(:,:)
        real(wp) :: dt_step, relax
        integer  :: nx, ny

        if (.not. hyd%now%initialized) then
            write(*,*) "hydro_update:: error: hydro_init_state must be called before hydro_update."
            stop
        end if

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        allocate(H_w_old(nx,ny))
        H_w_old = hyd%now%H_w

        dt_step      = max(time - hyd%now%time, 0.0_wp)
        hyd%now%time = time
        hyd%now%dt   = dt_step

        select case (hyd%par%method)

            case (HYDRO_METHOD_NONE)
                continue

            case (HYDRO_METHOD_BUCKET)
                if (dt_step > 0.0_wp) then
                    call calc_bucket(hyd%now%H_w, f_ice, f_grnd, mask, bmb_w, &
                                     dt_step, hyd%par%bucket%till_rate, hyd%par%H_w_max)
                end if
                call apply_N_closure(hyd, H_ice, z_bed, z_sl, f_ice, f_grnd)

            case (HYDRO_METHOD_K24)
                if (dt_step > 0.0_wp) then
                    allocate(H_ice_dp(nx,ny), z_bed_dp(nx,ny), mask_dp(nx,ny))
                    allocate(bmb_w_dp(nx,ny), uxy_b_dp(nx,ny), A_glen_dp(nx,ny))
                    allocate(kappa_dp(nx,ny))
                    allocate(q_x_dp(nx,ny), q_y_dp(nx,ny), N_dp(nx,ny), p_w_dp(nx,ny), H_w_eq_dp(nx,ny))

                    H_ice_dp  = real(H_ice,         dp)
                    z_bed_dp  = real(z_bed,         dp)
                    mask_dp   = real(mask,          dp)
                    bmb_w_dp  = real(bmb_w,         dp)
                    uxy_b_dp  = real(uxy_b,         dp)
                    A_glen_dp = real(A_glen,        dp)
                    kappa_dp  = real(hyd%now%kappa, dp)

                    call calc_k24(q_x_dp, q_y_dp, N_dp, p_w_dp, H_w_eq_dp, &
                                  H_ice_dp, z_bed_dp, mask_dp, bmb_w_dp, uxy_b_dp, A_glen_dp, &
                                  kappa_dp, &
                                  real(hyd%par%dx, dp), hyd%par%k24)

                    hyd%now%q_x = real(q_x_dp, wp)
                    hyd%now%q_y = real(q_y_dp, wp)
                    hyd%now%N   = real(N_dp,   wp)
                    hyd%now%p_w = real(p_w_dp, wp)

                    if (hyd%par%k24_update_H_w) then
                        allocate(H_w_eq(nx,ny))
                        H_w_eq = real(H_w_eq_dp, wp)
                        if (hyd%par%tau_H <= 0.0_wp) then
                            hyd%now%H_w = H_w_eq
                        else
                            relax       = 1.0_wp - exp(-dt_step / hyd%par%tau_H)
                            hyd%now%H_w = hyd%now%H_w + (H_w_eq - hyd%now%H_w) * relax
                        end if
                        deallocate(H_w_eq)
                    end if

                    deallocate(H_ice_dp, z_bed_dp, mask_dp)
                    deallocate(bmb_w_dp, uxy_b_dp, A_glen_dp)
                    deallocate(kappa_dp)
                    deallocate(q_x_dp, q_y_dp, N_dp, p_w_dp, H_w_eq_dp)
                end if

            case default
                write(*,*) "hydro_update:: error: method must be one of [0,1,2]."
                write(*,*) "method = ", hyd%par%method
                stop

        end select

        if (dt_step > 0.0_wp) then
            hyd%now%dHwdt = (hyd%now%H_w - H_w_old) / dt_step
        else
            hyd%now%dHwdt = 0.0_wp
        end if

        deallocate(H_w_old)

        return

    end subroutine hydro_update

    ! ------------------------------------------------------------
    ! N-closure post-step for BUCKET mode. Writes hyd%now%N and derives
    ! p_w = Po - N. No-op when par%bucket%N_closure == N_CLOSURE_NONE.
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
                do j = 1, ny
                do i = 1, nx
                    call hydro_calc_N_overburden(hyd%now%N(i,j), H_ice(i,j), f_ice(i,j), f_grnd(i,j), &
                                                 hyd%par%closures%rho_ice, hyd%par%closures%g)
                end do
                end do

            case (N_CLOSURE_MARINE)
                do j = 1, ny
                do i = 1, nx
                    call hydro_calc_N_marine(hyd%now%N(i,j), H_ice(i,j), f_ice(i,j), z_bed(i,j), z_sl(i,j), &
                                             hyd%par%closures%marine%p, hyd%par%closures%rho_ice, &
                                             hyd%par%closures%marine%rho_sw, hyd%par%closures%g)
                end do
                end do

            case (N_CLOSURE_TILL)
                do j = 1, ny
                do i = 1, nx
                    call hydro_calc_N_till(hyd%now%N(i,j), hyd%now%H_w(i,j), H_ice(i,j), &
                                           f_ice(i,j), f_grnd(i,j), hyd%par%H_w_max, &
                                           hyd%par%closures%till%N0, hyd%par%closures%till%delta, &
                                           hyd%par%closures%till%e0, hyd%par%closures%till%Cc, &
                                           hyd%par%closures%rho_ice, hyd%par%closures%g)
                end do
                end do

            case default
                write(*,*) "apply_N_closure:: error: unsupported post-step closure ", &
                           hyd%par%bucket%N_closure
                write(*,*) "(N_CLOSURE_TWO_VALUE is standalone-only.)"
                stop

        end select

        ! Derive p_w = Po - N on grounded cells
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

        return

    end subroutine apply_N_closure

    subroutine hydro_par_load(par, filename, init)

        implicit none

        type(hydro_param_class), intent(INOUT) :: par
        character(len=*),        intent(IN)    :: filename
        logical, optional,       intent(IN)    :: init

        logical :: init_pars

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        par%method         = HYDRO_METHOD_NONE
        par%init_method    = HYDRO_INIT_ZERO
        par%dx             = 0.0_wp
        par%dy             = 0.0_wp
        par%H_w_max        = 2.0_wp
        par%tau_H          = 0.0_wp
        par%k24_update_H_w = .true.

        call nml_read(filename,"fast_hydrology","method",         par%method,         init=init_pars)
        call nml_read(filename,"fast_hydrology","init_method",    par%init_method,    init=init_pars)
        call nml_read(filename,"fast_hydrology","dx",             par%dx,             init=init_pars)
        call nml_read(filename,"fast_hydrology","dy",             par%dy,             init=init_pars)
        call nml_read(filename,"fast_hydrology","H_w_max",        par%H_w_max,        init=init_pars)
        call nml_read(filename,"fast_hydrology","tau_H",          par%tau_H,          init=init_pars)
        call nml_read(filename,"fast_hydrology","k24_update_H_w", par%k24_update_H_w, init=init_pars)

        call bucket_par_load (par%bucket,   filename, init=init_pars)
        call k24_par_load    (par%k24,      filename, init=init_pars)
        call closure_par_load(par%closures, filename, init=init_pars)

        return

    end subroutine hydro_par_load

    subroutine hydro_allocate(now, nx, ny)

        implicit none

        type(hydro_state_class), intent(INOUT) :: now
        integer,                 intent(IN)    :: nx, ny

        call hydro_deallocate(now)

        allocate(now%H_w(nx,ny))
        allocate(now%dHwdt(nx,ny))
        allocate(now%p_w(nx,ny))
        allocate(now%q_x(nx,ny))
        allocate(now%q_y(nx,ny))
        allocate(now%N(nx,ny))
        allocate(now%kappa(nx,ny))

        now%H_w   = 0.0_wp
        now%dHwdt = 0.0_wp
        now%p_w   = 0.0_wp
        now%q_x   = 0.0_wp
        now%q_y   = 0.0_wp
        now%N     = 0.0_wp
        now%kappa = 0.0_wp

        return

    end subroutine hydro_allocate

    subroutine hydro_deallocate(now)

        implicit none

        type(hydro_state_class), intent(INOUT) :: now

        if (allocated(now%H_w))   deallocate(now%H_w)
        if (allocated(now%dHwdt)) deallocate(now%dHwdt)
        if (allocated(now%p_w))   deallocate(now%p_w)
        if (allocated(now%q_x))   deallocate(now%q_x)
        if (allocated(now%q_y))   deallocate(now%q_y)
        if (allocated(now%N))     deallocate(now%N)
        if (allocated(now%kappa)) deallocate(now%kappa)

        return

    end subroutine hydro_deallocate

end module fast_hydrology
