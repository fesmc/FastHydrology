module fast_hydrology_bucket
    ! Local "bucket" basal-hydrology model for till water storage. Per-cell
    ! mass balance with background drainage to a deeper substrate and a hard
    ! cap on storage. Adapted from yelmo's calc_basal_water_local
    ! (yelmo/src/physics/thermodynamics.f90), with the W_til_max overrides on
    ! floating / adjacent-to-floating cells removed so the bucket and K24
    ! transport behave comparably.
    !
    ! Notation follows van Pelt & Bueler 2015 (BvP15): W_til is the till
    ! water storage layer thickness (m); bkt_till_rate is a constant
    ! background drainage to a deeper unspecified substrate (m/a in the
    ! namelist; converted to m/s internally -- see fast_hydrology). The
    ! distributed sheet thickness W (between till and ice) is produced by
    ! the transport model (K24), not this bucket.
    !
    ! Cell logic (within host mask == 1):
    !   f_grnd > 0 .and. f_ice > 0 :
    !       W_attempt = W_til + dt * (mdot - till_rate)
    !       W_til     = clamp(W_attempt, 0, W_til_max(i,j))
    !       overflow  = max(0, (W_attempt - W_til_max(i,j)) / dt)
    !   otherwise (ocean / floating / grounded ice-free) :
    !       W_til    = 0
    !       overflow = 0
    !
    ! `overflow` is the till-saturation overflow rate (units: same as mdot;
    ! see fast_hydrology for the m/s convention). It is the source term the
    ! transport model (K24) consumes.
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
    ! Controls how W_til is treated on the outer halo of the domain (the
    ! i=1, i=nx, j=1, j=ny rim of cells). Floating-cell logic is independent
    ! and always applied; grounded-ice-free cells are always set to 0.
    integer, parameter, public :: MASK_BC_ZERO    = 0  ! W_til = 0 at the rim
    integer, parameter, public :: MASK_BC_IMPOSED = 1  ! W_til = par%W_til_bc at the rim
    integer, parameter, public :: MASK_BC_MIRROR  = 2  ! W_til mirrored from the inward neighbor (Neumann)

    type bucket_param_class
        real(wp) :: till_rate          ! background till drainage [units: same as mdot]
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

    subroutine calc_bucket(W_til, overflow, f_ice, f_grnd, mask, mdot, dt, par, W_til_max)
        ! Sequential bucket step. Updates W_til in place on any partially-or-
        ! fully grounded ice cell (f_grnd > 0 .and. f_ice > 0); zero elsewhere.
        ! Returns the per-cell overflow rate (same units as mdot): the share
        ! of the source that did NOT fit in the till and is available to feed
        ! the transport model downstream. No W_til_max fill is applied on
        ! floating / adjacent-to-floating cells; that is the caller's job
        ! via apply_floating_override + apply_mask_bc.

        implicit none

        real(wp),                 intent(INOUT) :: W_til(:,:)
        real(wp),                 intent(OUT)   :: overflow(:,:)
        real(wp),                 intent(IN)    :: f_ice(:,:)
        real(wp),                 intent(IN)    :: f_grnd(:,:)
        real(wp),                 intent(IN)    :: mask(:,:)
        real(wp),                 intent(IN)    :: mdot(:,:)
        real(wp),                 intent(IN)    :: dt
        type(bucket_param_class), intent(IN)    :: par
        real(wp),                 intent(IN)    :: W_til_max(:,:)

        integer  :: i, j, nx, ny
        real(wp) :: W_attempt, cap

        nx = size(f_ice,1)
        ny = size(f_ice,2)

        overflow = 0.0_wp

        !$omp parallel do default(shared) private(i,j,W_attempt,cap) schedule(static)
        do j = 1, ny
        do i = 1, nx

            if (mask(i,j) /= 1.0_wp) cycle

            if (f_grnd(i,j) > 0.0_wp .and. f_ice(i,j) > 0.0_wp) then

                cap       = W_til_max(i,j)
                W_attempt = W_til(i,j) + dt * (mdot(i,j) - par%till_rate)

                if (W_attempt > cap) then
                    overflow(i,j) = (W_attempt - cap) / dt
                    W_til(i,j)    = cap
                else if (W_attempt < 0.0_wp) then
                    overflow(i,j) = 0.0_wp
                    W_til(i,j)    = 0.0_wp
                else
                    overflow(i,j) = 0.0_wp
                    W_til(i,j)    = W_attempt
                end if

            else
                W_til(i,j)    = 0.0_wp
                overflow(i,j) = 0.0_wp
            end if

        end do
        end do
        !$omp end parallel do

        return

    end subroutine calc_bucket

    ! ------------------------------------------------------------
    ! Zero W_til on cells that are not partially-or-fully grounded ice (i.e.
    ! ocean, floating, and grounded ice-free). The active set for both
    ! bucket and K24 is the same: f_grnd > 0 .and. f_ice > 0. Active cells
    ! are left untouched. No W_til_max fill is applied.
    ! ------------------------------------------------------------
    subroutine apply_floating_override(W_til, f_ice, f_grnd)

        implicit none

        real(wp), intent(INOUT) :: W_til(:,:)
        real(wp), intent(IN)    :: f_ice(:,:)
        real(wp), intent(IN)    :: f_grnd(:,:)

        integer :: i, j, nx, ny

        nx = size(f_ice,1)
        ny = size(f_ice,2)

        !$omp parallel do default(shared) private(i,j) schedule(static)
        do j = 1, ny
        do i = 1, nx
            if (.not. (f_grnd(i,j) > 0.0_wp .and. f_ice(i,j) > 0.0_wp)) then
                W_til(i,j) = 0.0_wp
            end if
        end do
        end do
        !$omp end parallel do

        return

    end subroutine apply_floating_override

    subroutine apply_mask_bc(W_til, mask_bc, W_til_bc)
        ! Apply the domain-border BC at the i=1, i=nx, j=1, j=ny rim.
        ! Acts after the main update and after the floating-cell logic.

        implicit none

        real(wp), intent(INOUT) :: W_til(:,:)
        integer,  intent(IN)    :: mask_bc
        real(wp), intent(IN)    :: W_til_bc

        integer :: nx, ny

        nx = size(W_til,1)
        ny = size(W_til,2)

        select case (mask_bc)

            case (MASK_BC_ZERO)
                W_til(1,:)  = 0.0_wp
                W_til(nx,:) = 0.0_wp
                W_til(:,1)  = 0.0_wp
                W_til(:,ny) = 0.0_wp

            case (MASK_BC_IMPOSED)
                W_til(1,:)  = W_til_bc
                W_til(nx,:) = W_til_bc
                W_til(:,1)  = W_til_bc
                W_til(:,ny) = W_til_bc

            case (MASK_BC_MIRROR)
                if (nx >= 2) then
                    W_til(1,:)  = W_til(2,:)
                    W_til(nx,:) = W_til(nx-1,:)
                end if
                if (ny >= 2) then
                    W_til(:,1)  = W_til(:,2)
                    W_til(:,ny) = W_til(:,ny-1)
                end if

            case default
                write(*,*) "apply_mask_bc:: error: unknown mask_bc = ", mask_bc
                stop

        end select

        return

    end subroutine apply_mask_bc

end module fast_hydrology_bucket
