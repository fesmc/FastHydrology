module fast_hydrology_closures
    ! Effective-pressure (N) closures parameterized as N = f(H_w, geometry).
    ! All ported from yelmo's basal_dragging.f90 with parameters re-routed
    ! into per-closure parameter sub-structs.
    !
    ! The standalone subs are public and callable regardless of the host's
    ! choice of par%method; they take their inputs explicitly. The wrapper
    ! (fast_hydrology) selects one of them as a post-step when
    ! par%method = BUCKET and par%bucket%N_closure /= NONE.

    use nml

    implicit none

    integer, parameter :: dp = kind(1.d0)
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: wp = sp

    ! ---------- Closure enum (par%bucket%N_closure) ----------
    integer, parameter, public :: N_CLOSURE_NONE       = 0
    integer, parameter, public :: N_CLOSURE_OVERBURDEN = 1
    integer, parameter, public :: N_CLOSURE_MARINE     = 2  ! Leguy 2014
    integer, parameter, public :: N_CLOSURE_TILL       = 3  ! van Pelt & Bueler 2015
    integer, parameter, public :: N_CLOSURE_TWO_VALUE  = 4  ! standalone-only

    type closure_marine_param_class
        real(wp) :: p           ! [0:1] ocean connectivity exponent
        real(wp) :: rho_sw      ! [kg/m3]
    end type

    type closure_till_param_class
        real(wp) :: N0          ! [Pa]  reference effective pressure
        real(wp) :: delta       ! [-]   fraction of overburden for saturated till
        real(wp) :: e0          ! [-]   reference void ratio at N0
        real(wp) :: Cc          ! [-]   till compressibility
    end type

    type closure_two_value_param_class
        real(wp) :: delta       ! [-]   fraction reduction at PMP
    end type

    type closure_param_class
        ! constants shared across closures
        real(wp) :: rho_ice
        real(wp) :: g
        type(closure_marine_param_class)    :: marine
        type(closure_till_param_class)      :: till
        type(closure_two_value_param_class) :: two_value
    end type

    private
    public :: closure_param_class
    public :: closure_marine_param_class
    public :: closure_till_param_class
    public :: closure_two_value_param_class
    public :: closure_par_load
    public :: hydro_calc_N_overburden
    public :: hydro_calc_N_marine
    public :: hydro_calc_N_till
    public :: hydro_calc_N_two_value

contains

    subroutine closure_par_load(par, filename, group, init)

        implicit none

        type(closure_param_class), intent(INOUT) :: par
        character(len=*),          intent(IN)    :: filename
        character(len=*),          intent(IN)    :: group
        logical, optional,         intent(IN)    :: init

        logical :: init_pars

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        par%rho_ice         = 917.0_wp
        par%g               =   9.81_wp
        par%marine%p        =   1.0_wp
        par%marine%rho_sw   = 1028.0_wp
        par%till%N0         = 1000.0_wp
        par%till%delta      =    0.02_wp
        par%till%e0         =    0.69_wp
        par%till%Cc         =    0.12_wp
        par%two_value%delta =    0.02_wp

        call nml_read(filename,group,"rho_ice",        par%rho_ice,         init=init_pars)
        call nml_read(filename,group,"g",              par%g,               init=init_pars)
        call nml_read(filename,group,"marine_p",       par%marine%p,        init=init_pars)
        call nml_read(filename,group,"marine_rho_sw",  par%marine%rho_sw,   init=init_pars)
        call nml_read(filename,group,"till_N0",        par%till%N0,         init=init_pars)
        call nml_read(filename,group,"till_delta",     par%till%delta,      init=init_pars)
        call nml_read(filename,group,"till_e0",        par%till%e0,         init=init_pars)
        call nml_read(filename,group,"till_Cc",        par%till%Cc,         init=init_pars)
        call nml_read(filename,group,"two_value_delta",par%two_value%delta, init=init_pars)

        return

    end subroutine closure_par_load

    ! ============================================================
    ! Overburden: N = rho_ice * g * H_eff
    ! ============================================================
    elemental subroutine hydro_calc_N_overburden(N, H_ice, f_ice, f_grnd, rho_ice, g)

        real(wp), intent(OUT) :: N
        real(wp), intent(IN)  :: H_ice, f_ice, f_grnd, rho_ice, g

        real(wp) :: H_eff

        call calc_H_eff(H_eff, H_ice, f_ice)

        if (f_grnd > 0.0_wp) then
            N = rho_ice * g * H_eff
        else
            N = 0.0_wp
        end if

    end subroutine hydro_calc_N_overburden

    ! ============================================================
    ! Marine: Leguy 2014, Eq. 14 (geometry only; no H_w dependence)
    ! ============================================================
    elemental subroutine hydro_calc_N_marine(N, H_ice, f_ice, z_bed, z_sl, p, rho_ice, rho_sw, g)

        real(wp), intent(OUT) :: N
        real(wp), intent(IN)  :: H_ice, f_ice, z_bed, z_sl, p, rho_ice, rho_sw, g

        real(wp) :: H_eff, H_float, p_w, x, rho_sw_ice

        rho_sw_ice = rho_sw / rho_ice
        H_float    = max(0.0_wp, rho_sw_ice * (z_sl - z_bed))

        call calc_H_eff(H_eff, H_ice, f_ice)

        if (H_eff == 0.0_wp) then
            p_w = 0.0_wp
        else if (H_eff < H_float) then
            p_w = rho_ice * g * H_eff
        else
            x   = min(1.0_wp, H_float / H_eff)
            p_w = rho_ice * g * H_eff * (1.0_wp - (1.0_wp - x)**p)
        end if

        N = rho_ice * g * H_eff - p_w

    end subroutine hydro_calc_N_marine

    ! ============================================================
    ! Till: van Pelt & Bueler 2015, Eq. 23
    ! ============================================================
    elemental subroutine hydro_calc_N_till(N, H_w, H_ice, f_ice, f_grnd, H_w_max, &
                                          N0, delta, e0, Cc, rho_ice, g)

        real(wp), intent(OUT) :: N
        real(wp), intent(IN)  :: H_w, H_ice, f_ice, f_grnd
        real(wp), intent(IN)  :: H_w_max, N0, delta, e0, Cc, rho_ice, g

        real(wp) :: H_eff, P0, s, q1

        if (f_grnd == 0.0_wp) then
            N = 0.0_wp
            return
        end if

        call calc_H_eff(H_eff, H_ice, f_ice)

        P0 = rho_ice * g * H_eff
        s  = min(H_w / H_w_max, 1.0_wp)
        q1 = (e0 / Cc) * (1.0_wp - s)
        q1 = min(q1, 10.0_wp)

        N = min( N0 * (delta * P0 / N0)**s * 10.0_wp**q1, P0 )

    end subroutine hydro_calc_N_till

    ! ============================================================
    ! Two-value: weighted average of overburden and delta*overburden by f_pmp
    ! ============================================================
    elemental subroutine hydro_calc_N_two_value(N, f_pmp, H_ice, f_ice, f_grnd, delta, rho_ice, g)

        real(wp), intent(OUT) :: N
        real(wp), intent(IN)  :: f_pmp, H_ice, f_ice, f_grnd, delta, rho_ice, g

        real(wp) :: H_eff, P0, P1

        if (f_grnd == 0.0_wp) then
            N = 0.0_wp
            return
        end if

        call calc_H_eff(H_eff, H_ice, f_ice)

        P0 = rho_ice * g * H_eff
        P1 = P0 * delta
        N  = P0 * (1.0_wp - f_pmp) + P1 * f_pmp

    end subroutine hydro_calc_N_two_value

    ! ------------------------------------------------------------
    ! Helper: H_eff with margin treatment (set_frac_zero=true semantics)
    ! ------------------------------------------------------------
    elemental subroutine calc_H_eff(H_eff, H_ice, f_ice)

        real(wp), intent(OUT) :: H_eff
        real(wp), intent(IN)  :: H_ice, f_ice

        if (f_ice > 0.0_wp) then
            H_eff = H_ice / f_ice
        else
            H_eff = H_ice
        end if
        if (f_ice < 1.0_wp) H_eff = 0.0_wp

    end subroutine calc_H_eff

end module fast_hydrology_closures
