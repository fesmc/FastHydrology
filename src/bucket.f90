module fast_hydrology_bucket
    ! Local "bucket" basal-hydrology model: per-cell water mass balance with
    ! linear till drainage and a hard cap. Adapted from yelmo's
    ! calc_basal_water_local (yelmo/src/physics/thermodynamics.f90), with the
    ! H_w_max overrides on floating / adjacent-to-floating cells removed so
    ! the bucket and K24 methods behave comparably.
    !
    ! Cell logic (within host mask == 1):
    !   f_grnd > 0 .and. f_ice > 0 :
    !       H_w := H_w + dt * (bmb_w - till_rate), clamped to [0, H_w_max]
    !   otherwise (ocean / floating / grounded ice-free) :
    !       H_w = 0
    !
    ! Host-supplied `mask` is intersected with these cases: where mask /= 1
    ! the cell is left untouched (host has decided it is outside the active
    ! domain).

    use nml

    implicit none

    integer, parameter :: dp = kind(1.d0)
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: wp = sp

    ! ---------- Domain-border BC enum (par%mask_bc) ----------
    ! Controls how H_w is treated on the outer halo of the domain (the
    ! i=1, i=nx, j=1, j=ny rim of cells). Yelmo's floating-cell logic
    ! (H_w = H_w_max on floating + adjacent-to-floating) is independent
    ! and always applied; grounded-ice-free cells are always set to 0.
    integer, parameter, public :: MASK_BC_ZERO    = 0  ! H_w = 0 at the rim
    integer, parameter, public :: MASK_BC_IMPOSED = 1  ! H_w = par%H_w_bc at the rim
    integer, parameter, public :: MASK_BC_MIRROR  = 2  ! H_w mirrored from the inward neighbor (Neumann)

    type bucket_param_class
        real(wp) :: till_rate          ! [m/a] drainage rate (water equiv.)
        integer  :: N_closure          ! see fast_hydrology_closures::N_CLOSURE_*
    end type

    private
    public :: bucket_param_class
    public :: bucket_par_load
    public :: calc_bucket
    public :: apply_floating_override
    public :: apply_mask_bc

contains

    subroutine bucket_par_load(par, filename, group, init)

        implicit none

        type(bucket_param_class), intent(INOUT) :: par
        character(len=*),         intent(IN)    :: filename
        character(len=*),         intent(IN)    :: group
        logical, optional,        intent(IN)    :: init

        logical :: init_pars

        init_pars = .FALSE.
        if (present(init)) init_pars = init

        par%till_rate = 1.0e-3_wp
        par%N_closure = 0

        call nml_read(filename,group,"bkt_till_rate", par%till_rate, init=init_pars)
        call nml_read(filename,group,"bkt_N_closure", par%N_closure, init=init_pars)

        return

    end subroutine bucket_par_load

    subroutine calc_bucket(H_w, f_ice, f_grnd, mask, bmb_w, dt, par, H_w_max)
        ! Local bucket update: evolve H_w on any partially-or-fully grounded
        ! ice cell (f_grnd > 0 .and. f_ice > 0); zero elsewhere. No H_w_max
        ! fill on floating / adjacent-to-floating cells. Domain-border BC is
        ! applied separately by apply_mask_bc.

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

        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
        do i = 1, nx

            if (mask(i,j) /= 1.0_wp) cycle

            if (f_grnd(i,j) > 0.0_wp .and. f_ice(i,j) > 0.0_wp) then
                H_w(i,j) = H_w(i,j) + dt * (bmb_w(i,j) - par%till_rate)
                H_w(i,j) = max(H_w(i,j), 0.0_wp)
                H_w(i,j) = min(H_w(i,j), H_w_max)
            else
                H_w(i,j) = 0.0_wp
            end if

        end do
        end do
        !$omp end parallel do

        return

    end subroutine calc_bucket

    ! ------------------------------------------------------------
    ! Zero H_w on cells that are not partially-or-fully grounded ice (i.e.
    ! ocean, floating, and grounded ice-free). The active set for both
    ! bucket and K24 is the same: f_grnd > 0 .and. f_ice > 0. Active cells
    ! are left untouched. No H_w_max fill is applied.
    ! ------------------------------------------------------------
    subroutine apply_floating_override(H_w, f_ice, f_grnd)

        implicit none

        real(wp), intent(INOUT) :: H_w(:,:)
        real(wp), intent(IN)    :: f_ice(:,:)
        real(wp), intent(IN)    :: f_grnd(:,:)

        integer :: i, j, nx, ny

        nx = size(f_ice,1)
        ny = size(f_ice,2)

        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
        do i = 1, nx
            if (.not. (f_grnd(i,j) > 0.0_wp .and. f_ice(i,j) > 0.0_wp)) then
                H_w(i,j) = 0.0_wp
            end if
        end do
        end do
        !$omp end parallel do

        return

    end subroutine apply_floating_override

    subroutine apply_mask_bc(H_w, mask_bc, H_w_bc)
        ! Apply the domain-border BC at the i=1, i=nx, j=1, j=ny rim.
        ! Acts after the main update and after the floating-cell logic.

        implicit none

        real(wp), intent(INOUT) :: H_w(:,:)
        integer,  intent(IN)    :: mask_bc
        real(wp), intent(IN)    :: H_w_bc

        integer :: nx, ny

        nx = size(H_w,1)
        ny = size(H_w,2)

        select case (mask_bc)

            case (MASK_BC_ZERO)
                H_w(1,:)  = 0.0_wp
                H_w(nx,:) = 0.0_wp
                H_w(:,1)  = 0.0_wp
                H_w(:,ny) = 0.0_wp

            case (MASK_BC_IMPOSED)
                H_w(1,:)  = H_w_bc
                H_w(nx,:) = H_w_bc
                H_w(:,1)  = H_w_bc
                H_w(:,ny) = H_w_bc

            case (MASK_BC_MIRROR)
                if (nx >= 2) then
                    H_w(1,:)  = H_w(2,:)
                    H_w(nx,:) = H_w(nx-1,:)
                end if
                if (ny >= 2) then
                    H_w(:,1)  = H_w(:,2)
                    H_w(:,ny) = H_w(:,ny-1)
                end if

            case default
                write(*,*) "apply_mask_bc:: error: unknown mask_bc = ", mask_bc
                stop

        end select

        return

    end subroutine apply_mask_bc

end module fast_hydrology_bucket
