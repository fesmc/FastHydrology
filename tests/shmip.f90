program shmip
    ! Stand-alone SHMIP-style driver for the FastHydrology library.
    !
    ! Currently implements SHMIP case A1 only (steady square slab with uniform
    ! distributed melt). Other SHMIP cases (A2..A6, B*, C*, D*) plug in by
    ! adding new `case` branches that fill (z_bed, H_ice, bmb_w, ...) and
    ! set the relevant time-stepping.
    !
    ! Configuration is via namelist file passed on the command line:
    !     ./shmip.x par/shmip.nml
    ! All FastHydrology parameters live in the same file (read by hydro_init).

    use nml
    use fast_hydrology

    implicit none

    integer, parameter :: wp_local = kind(1.0)
    integer, parameter :: dp = kind(1.d0)

    character(len=128) :: nml_file
    character(len=16)  :: shmip_case
    integer            :: nx, ny
    real(wp_local)     :: dx, dy
    real(wp_local)     :: t_start, t_end, dt_out, dt_step
    real(wp_local)     :: A_glen_const, uxy_b_const, melt_const

    real(wp_local), allocatable :: H_ice(:,:), z_bed(:,:), z_sl(:,:)
    real(wp_local), allocatable :: f_ice(:,:), f_grnd(:,:), mask(:,:)
    real(wp_local), allocatable :: bmb_w(:,:), uxy_b(:,:), A_glen(:,:)

    type(hydro_class) :: hyd
    real(wp_local)    :: time
    integer           :: nargs, n_out, i_out

    ! ---- Command line ----
    nargs = command_argument_count()
    if (nargs < 1) then
        write(*,*) "usage: shmip.x <namelist>"
        stop
    end if
    call get_command_argument(1, nml_file)

    ! ---- Driver config ----
    call nml_read(nml_file, "shmip", "case",         shmip_case)
    call nml_read(nml_file, "shmip", "nx",           nx)
    call nml_read(nml_file, "shmip", "ny",           ny)
    call nml_read(nml_file, "shmip", "t_start",      t_start)
    call nml_read(nml_file, "shmip", "t_end",        t_end)
    call nml_read(nml_file, "shmip", "dt_step",      dt_step)
    call nml_read(nml_file, "shmip", "dt_out",       dt_out)
    call nml_read(nml_file, "shmip", "A_glen_const", A_glen_const)
    call nml_read(nml_file, "shmip", "uxy_b_const",  uxy_b_const)
    call nml_read(nml_file, "shmip", "melt_const",   melt_const)

    write(*,'(a)') "FastHydrology SHMIP driver"
    write(*,'(a,a)')      "  case    = ", trim(shmip_case)
    write(*,'(a,i6,a,i6)')"  grid    = ", nx, " x ", ny

    ! ---- Build geometry per case ----
    allocate(H_ice(nx,ny), z_bed(nx,ny), z_sl(nx,ny))
    allocate(f_ice(nx,ny), f_grnd(nx,ny), mask(nx,ny))
    allocate(bmb_w(nx,ny), uxy_b(nx,ny), A_glen(nx,ny))

    select case (trim(shmip_case))
        case ("A1")
            ! SHMIP A1: square slab on a slope. 100km x 20km domain, 1500m
            ! upper-ice thickness, surface follows sqrt-profile per SHMIP
            ! protocol. Bed is flat at z=0; surface elevation s(x) = 6*(sqrt(x+5000) - sqrt(5000)) + 1.
            ! (Standard SHMIP profile.)
            call setup_shmip_A(nx, ny, H_ice, z_bed, dx)
            dy = dx
            melt_const = max(melt_const, 7.93e-11_wp_local * 31556926.0_wp_local) ! 7.93e-11 m/s -> m/a
        case default
            write(*,*) "shmip:: case '"//trim(shmip_case)//"' not implemented in this driver yet."
            write(*,*) "         Supported: A1"
            stop
    end select

    ! Uniform host state
    z_sl   = 0.0_wp_local
    f_ice  = 1.0_wp_local
    f_grnd = 1.0_wp_local
    mask   = 1.0_wp_local
    bmb_w  = melt_const
    uxy_b  = uxy_b_const
    A_glen = A_glen_const

    ! ---- Initialize library ----
    call hydro_init(hyd, nml_file, nx, ny)
    hyd%par%dx = dx
    hyd%par%dy = dy
    call hydro_init_state(hyd, z_bed, f_ice, f_grnd, t_start)

    write(*,'(a,i0,a,i0,a)') "  method  = ", hyd%par%method,    " (0=NONE 1=BUCKET 2=K24)"
    write(*,'(a,i0)')        "  N_clos  = ", hyd%par%bucket%N_closure

    ! ---- Time-step loop ----
    time  = t_start
    i_out = 0
    n_out = max(1, nint(dt_out / dt_step))

    write(*,'(a)') ""
    write(*,'(a)') "       time      H_w_max     H_w_mean       N_max      N_mean"

    do while (time < t_end)
        time = time + dt_step
        call hydro_update(hyd, H_ice, z_bed, z_sl, f_ice, f_grnd, mask, &
                          bmb_w, uxy_b, A_glen, time)

        i_out = i_out + 1
        if (mod(i_out, n_out) == 0 .or. time >= t_end) then
            write(*,'(f12.4, 4es13.4)') time, &
                maxval(hyd%now%H_w), sum(hyd%now%H_w)/(nx*ny), &
                maxval(hyd%now%N),   sum(hyd%now%N)/(nx*ny)
        end if
    end do

    write(*,'(a)') ""
    write(*,'(a)') "shmip: done."

contains

    subroutine setup_shmip_A(nx, ny, H_ice, z_bed, dx)
        ! SHMIP A geometry: 100km x 20km, sqrt-profile surface, flat bed.
        ! z_b = 0 ; z_s = 6*(sqrt(x+5000) - sqrt(5000)) + 1 ; H_ice = z_s - z_b.

        implicit none

        integer,        intent(IN)  :: nx, ny
        real(wp_local), intent(OUT) :: H_ice(nx,ny), z_bed(nx,ny)
        real(wp_local), intent(OUT) :: dx

        real(wp_local) :: x, z_s, L
        integer        :: i, j

        L  = 100.0e3_wp_local                          ! [m] x-length
        dx = L / real(nx - 1, wp_local)

        do j = 1, ny
        do i = 1, nx
            x         = real(i-1, wp_local) * dx
            z_s       = 6.0_wp_local * (sqrt(x + 5.0e3_wp_local) - sqrt(5.0e3_wp_local)) + 1.0_wp_local
            z_bed(i,j) = 0.0_wp_local
            H_ice(i,j) = max(z_s - z_bed(i,j), 0.0_wp_local)
        end do
        end do

    end subroutine setup_shmip_A

end program shmip
