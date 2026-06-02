module fast_hydrology_bucket
    ! Local "bucket" basal-hydrology model: per-cell water mass balance with
    ! linear till drainage and a hard cap. Ported verbatim from yelmo's
    ! calc_basal_water_local (yelmo/src/physics/thermodynamics.f90).
    !
    ! Cell logic:
    !   floating  (f_grnd == 0)                                  : H_w = H_w_max
    !   grounded fully-ice-covered & has floating neighbor       : H_w = H_w_max
    !   grounded ice-free (f_ice < 1)                            : H_w = 0
    !   grounded fully-ice-covered & no floating neighbor (else) :
    !       H_w := H_w + dt * (bmb_w - till_rate), clamped to [0, H_w_max]
    !
    ! Host-supplied `mask` is intersected with these cases: where mask /= 1
    ! the cell is left untouched (host has decided it is outside the active
    ! domain).

    use nml

    implicit none

    integer, parameter :: dp = kind(1.d0)
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: wp = sp

    type bucket_param_class
        real(wp) :: till_rate          ! [m/a] drainage rate (water equiv.)
        integer  :: N_closure          ! see fast_hydrology_closures::N_CLOSURE_*
    end type

    private
    public :: bucket_param_class
    public :: bucket_par_load
    public :: calc_bucket
    public :: apply_floating_override

contains

    subroutine bucket_par_load(par, filename, init)

        implicit none

        type(bucket_param_class), intent(INOUT) :: par
        character(len=*),         intent(IN)    :: filename
        logical, optional,        intent(IN)    :: init

        logical :: init_pars

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        par%till_rate = 1.0e-3_wp
        par%N_closure = 0

        call nml_read(filename,"fast_hydrology_bucket","till_rate", par%till_rate, init=init_pars)
        call nml_read(filename,"fast_hydrology_bucket","N_closure", par%N_closure, init=init_pars)

        return

    end subroutine bucket_par_load

    subroutine calc_bucket(H_w, f_ice, f_grnd, mask, bmb_w, dt, till_rate, H_w_max)

        implicit none

        real(wp), intent(INOUT) :: H_w(:,:)
        real(wp), intent(IN)    :: f_ice(:,:)
        real(wp), intent(IN)    :: f_grnd(:,:)
        real(wp), intent(IN)    :: mask(:,:)
        real(wp), intent(IN)    :: bmb_w(:,:)
        real(wp), intent(IN)    :: dt
        real(wp), intent(IN)    :: till_rate
        real(wp), intent(IN)    :: H_w_max

        integer :: i, j, nx, ny, im1, ip1, jm1, jp1

        nx = size(f_ice,1)
        ny = size(f_ice,2)

        do j = 1, ny
        do i = 1, nx

            if (mask(i,j) /= 1.0_wp) cycle

            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)

            if (f_grnd(i,j) == 0.0_wp) then
                H_w(i,j) = H_w_max
            else if (f_grnd(i,j) > 0.0_wp .and. f_ice(i,j) == 1.0_wp .and. &
                    (f_grnd(im1,j) == 0.0_wp .or. f_grnd(ip1,j) == 0.0_wp .or. &
                     f_grnd(i,jm1) == 0.0_wp .or. f_grnd(i,jp1) == 0.0_wp) ) then
                H_w(i,j) = H_w_max
            else if (f_grnd(i,j) > 0.0_wp .and. f_ice(i,j) < 1.0_wp) then
                H_w(i,j) = 0.0_wp
            else
                H_w(i,j) = H_w(i,j) + dt * (bmb_w(i,j) - till_rate)
                H_w(i,j) = max(H_w(i,j), 0.0_wp)
                H_w(i,j) = min(H_w(i,j), H_w_max)
            end if

        end do
        end do

        return

    end subroutine calc_bucket

    ! ------------------------------------------------------------
    ! Init-time floating override: independent of the time-stepping bucket
    ! call, this routine sets H_w = H_w_max on floating + adjacent-to-
    ! floating cells, leaving all other cells untouched. Grounded ice-free
    ! cells are also forced to H_w = 0 (matches yelmo's bucket convention).
    ! ------------------------------------------------------------
    subroutine apply_floating_override(H_w, f_ice, f_grnd, H_w_max)

        implicit none

        real(wp), intent(INOUT) :: H_w(:,:)
        real(wp), intent(IN)    :: f_ice(:,:)
        real(wp), intent(IN)    :: f_grnd(:,:)
        real(wp), intent(IN)    :: H_w_max

        integer :: i, j, nx, ny, im1, ip1, jm1, jp1

        nx = size(f_ice,1)
        ny = size(f_ice,2)

        do j = 1, ny
        do i = 1, nx

            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)

            if (f_grnd(i,j) == 0.0_wp) then
                H_w(i,j) = H_w_max
            else if (f_grnd(i,j) > 0.0_wp .and. f_ice(i,j) == 1.0_wp .and. &
                    (f_grnd(im1,j) == 0.0_wp .or. f_grnd(ip1,j) == 0.0_wp .or. &
                     f_grnd(i,jm1) == 0.0_wp .or. f_grnd(i,jp1) == 0.0_wp) ) then
                H_w(i,j) = H_w_max
            else if (f_grnd(i,j) > 0.0_wp .and. f_ice(i,j) < 1.0_wp) then
                H_w(i,j) = 0.0_wp
            end if

        end do
        end do

        return

    end subroutine apply_floating_override

end module fast_hydrology_bucket
