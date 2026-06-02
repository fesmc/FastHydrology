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

    ! ---------- Boundary-condition enum (par%bucket%mask_ice) ----------
    ! Controls how H_w is set at floating + adjacent-to-floating cells.
    ! Grounded-ice-free cells are always set to 0 (this is not a BC choice
    ! but a physical statement that no water sits there).
    integer, parameter, public :: MASK_ICE_ZERO    = 0  ! H_w = 0
    integer, parameter, public :: MASK_ICE_IMPOSED = 1  ! H_w = par%bucket%H_w_bc (default H_w_max)
    integer, parameter, public :: MASK_ICE_MIRROR  = 2  ! H_w mirrored from interior neighbor

    type bucket_param_class
        real(wp) :: till_rate          ! [m/a] drainage rate (water equiv.)
        integer  :: N_closure          ! see fast_hydrology_closures::N_CLOSURE_*
        integer  :: mask_ice           ! see MASK_ICE_* above
        real(wp) :: H_w_bc             ! imposed H_w at floating cells (MASK_ICE_IMPOSED)
    end type

    private
    public :: bucket_param_class
    public :: bucket_par_load
    public :: calc_bucket
    public :: apply_floating_override
    public :: apply_mask_ice_bc

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
        par%mask_ice  = MASK_ICE_IMPOSED
        par%H_w_bc    = 2.0_wp

        call nml_read(filename,"fast_hydrology_bucket","till_rate", par%till_rate, init=init_pars)
        call nml_read(filename,"fast_hydrology_bucket","N_closure", par%N_closure, init=init_pars)
        call nml_read(filename,"fast_hydrology_bucket","mask_ice",  par%mask_ice,  init=init_pars)
        call nml_read(filename,"fast_hydrology_bucket","H_w_bc",    par%H_w_bc,    init=init_pars)

        return

    end subroutine bucket_par_load

    subroutine calc_bucket(H_w, f_ice, f_grnd, mask, bmb_w, dt, par, H_w_max)

        implicit none

        real(wp),                 intent(INOUT) :: H_w(:,:)
        real(wp),                 intent(IN)    :: f_ice(:,:)
        real(wp),                 intent(IN)    :: f_grnd(:,:)
        real(wp),                 intent(IN)    :: mask(:,:)
        real(wp),                 intent(IN)    :: bmb_w(:,:)
        real(wp),                 intent(IN)    :: dt
        type(bucket_param_class), intent(IN)    :: par
        real(wp),                 intent(IN)    :: H_w_max

        integer :: i, j, nx, ny

        nx = size(f_ice,1)
        ny = size(f_ice,2)

        do j = 1, ny
        do i = 1, nx

            if (mask(i,j) /= 1.0_wp) cycle
            if (.not. is_grounded_ice_interior(i, j, f_ice, f_grnd)) cycle
            if (f_grnd(i,j) > 0.0_wp .and. f_ice(i,j) < 1.0_wp)       cycle

            H_w(i,j) = H_w(i,j) + dt * (bmb_w(i,j) - par%till_rate)
            H_w(i,j) = max(H_w(i,j), 0.0_wp)
            H_w(i,j) = min(H_w(i,j), H_w_max)

        end do
        end do

        ! Re-apply BCs (floating + grounded-ice-free) after the evolution
        ! step so the interior bucket update can't leak into BC cells.
        call apply_mask_ice_bc(H_w, f_ice, f_grnd, par)

        return

    end subroutine calc_bucket

    ! ------------------------------------------------------------
    ! Init-time floating override: independent of the time-stepping bucket
    ! call, this routine sets H_w = H_w_max on floating + adjacent-to-
    ! floating cells, leaving all other cells untouched. Grounded ice-free
    ! cells are also forced to H_w = 0 (matches yelmo's bucket convention).
    ! ------------------------------------------------------------
    subroutine apply_floating_override(H_w, f_ice, f_grnd, par)
        ! Wrapper kept for the init path. Delegates to apply_mask_ice_bc.

        implicit none

        real(wp),                 intent(INOUT) :: H_w(:,:)
        real(wp),                 intent(IN)    :: f_ice(:,:)
        real(wp),                 intent(IN)    :: f_grnd(:,:)
        type(bucket_param_class), intent(IN)    :: par

        call apply_mask_ice_bc(H_w, f_ice, f_grnd, par)

    end subroutine apply_floating_override

    subroutine apply_mask_ice_bc(H_w, f_ice, f_grnd, par)
        ! Apply boundary condition at non-grounded-ice cells:
        !   floating (f_grnd == 0):                      treated per par%mask_ice
        !   grounded ice-covered & adjacent to floating: treated per par%mask_ice
        !   grounded ice-free (f_grnd > 0, f_ice < 1):   H_w = 0 (fixed)
        ! Interior grounded-ice cells are left untouched.

        implicit none

        real(wp),                 intent(INOUT) :: H_w(:,:)
        real(wp),                 intent(IN)    :: f_ice(:,:)
        real(wp),                 intent(IN)    :: f_grnd(:,:)
        type(bucket_param_class), intent(IN)    :: par

        integer :: i, j, nx, ny, im1, ip1, jm1, jp1
        logical :: floating_bc

        nx = size(f_ice,1)
        ny = size(f_ice,2)

        do j = 1, ny
        do i = 1, nx

            im1 = max(i-1,1)
            ip1 = min(i+1,nx)
            jm1 = max(j-1,1)
            jp1 = min(j+1,ny)

            floating_bc = .false.
            if (f_grnd(i,j) == 0.0_wp) then
                floating_bc = .true.
            else if (f_grnd(i,j) > 0.0_wp .and. f_ice(i,j) == 1.0_wp .and. &
                    (f_grnd(im1,j) == 0.0_wp .or. f_grnd(ip1,j) == 0.0_wp .or. &
                     f_grnd(i,jm1) == 0.0_wp .or. f_grnd(i,jp1) == 0.0_wp) ) then
                floating_bc = .true.
            end if

            if (floating_bc) then
                select case (par%mask_ice)
                    case (MASK_ICE_ZERO)
                        H_w(i,j) = 0.0_wp
                    case (MASK_ICE_IMPOSED)
                        H_w(i,j) = par%H_w_bc
                    case (MASK_ICE_MIRROR)
                        H_w(i,j) = mirror_from_interior(i, j, H_w, f_ice, f_grnd, nx, ny)
                    case default
                        write(*,*) "apply_mask_ice_bc:: error: unknown mask_ice = ", par%mask_ice
                        stop
                end select
            else if (f_grnd(i,j) > 0.0_wp .and. f_ice(i,j) < 1.0_wp) then
                H_w(i,j) = 0.0_wp
            end if

        end do
        end do

        return

    end subroutine apply_mask_ice_bc

    pure function is_grounded_ice_interior(i, j, f_ice, f_grnd) result(is_int)
        integer,  intent(IN) :: i, j
        real(wp), intent(IN) :: f_ice(:,:), f_grnd(:,:)
        logical              :: is_int

        integer :: im1, ip1, jm1, jp1, nx, ny

        nx = size(f_ice,1)
        ny = size(f_ice,2)
        im1 = max(i-1,1); ip1 = min(i+1,nx)
        jm1 = max(j-1,1); jp1 = min(j+1,ny)

        is_int = f_grnd(i,j) > 0.0_wp .and. f_ice(i,j) == 1.0_wp .and. &
                 f_grnd(im1,j) > 0.0_wp .and. f_grnd(ip1,j) > 0.0_wp .and. &
                 f_grnd(i,jm1) > 0.0_wp .and. f_grnd(i,jp1) > 0.0_wp

    end function is_grounded_ice_interior

    pure function mirror_from_interior(i, j, H_w, f_ice, f_grnd, nx, ny) result(val)
        ! Average over neighbors that are grounded-ice-interior; fall back
        ! to H_w(i,j) itself if no eligible neighbor exists.
        integer,  intent(IN) :: i, j, nx, ny
        real(wp), intent(IN) :: H_w(:,:), f_ice(:,:), f_grnd(:,:)
        real(wp)             :: val

        integer  :: nn
        real(wp) :: acc

        acc = 0.0_wp
        nn  = 0
        if (i > 1)  then
            if (is_grounded_ice_interior(i-1, j, f_ice, f_grnd)) then
                acc = acc + H_w(i-1,j); nn = nn + 1
            end if
        end if
        if (i < nx) then
            if (is_grounded_ice_interior(i+1, j, f_ice, f_grnd)) then
                acc = acc + H_w(i+1,j); nn = nn + 1
            end if
        end if
        if (j > 1)  then
            if (is_grounded_ice_interior(i, j-1, f_ice, f_grnd)) then
                acc = acc + H_w(i,j-1); nn = nn + 1
            end if
        end if
        if (j < ny) then
            if (is_grounded_ice_interior(i, j+1, f_ice, f_grnd)) then
                acc = acc + H_w(i,j+1); nn = nn + 1
            end if
        end if

        if (nn > 0) then
            val = acc / real(nn, wp)
        else
            val = H_w(i,j)
        end if

    end function mirror_from_interior

end module fast_hydrology_bucket
