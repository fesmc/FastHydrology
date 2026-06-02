program shmip
    ! Stand-alone SHMIP driver for the FastHydrology library.
    !
    ! Cases:
    !   A1..A6  square slab + uniform distributed melt at SHMIP-spec rate
    !   B1..B5  square slab + moulin point sources (NOT YET IMPLEMENTED)
    !   C1..C4  A5 setup + diurnal moulin (NOT YET IMPLEMENTED)
    !   D1..D5  A5 setup + seasonal melt (NOT YET IMPLEMENTED)
    !
    ! Forcing & geometry are set from the case label; namelist controls
    ! library parameters, grid, and time stepping. Output written to a
    ! NetCDF file with H_w, dHwdt, N, p_w, q_x, q_y per output step.

    use nml
    use ncio
    use fast_hydrology

    implicit none

    integer, parameter :: wp_local = kind(1.0)
    integer, parameter :: dp = kind(1.d0)

    real(wp_local), parameter :: SEC_PER_YEAR = 31536000.0_wp_local   ! 365 day year (SHMIP convention)

    character(len=128) :: nml_file
    character(len=128) :: out_file
    character(len=16)  :: shmip_case
    integer            :: nx, ny
    real(wp_local)     :: dx, dy
    real(wp_local)     :: t_start, t_end, dt_out, dt_step
    real(wp_local)     :: A_glen_const, uxy_b_const
    real(wp_local)     :: melt_rate                  ! [m/a] water-equivalent

    real(wp_local), allocatable :: H_ice(:,:), z_bed(:,:), z_sl(:,:)
    real(wp_local), allocatable :: f_ice(:,:), f_grnd(:,:), mask(:,:)
    real(wp_local), allocatable :: bmb_w(:,:), uxy_b(:,:), A_glen(:,:)
    real(wp_local), allocatable :: xc(:), yc(:)

    type(hydro_class) :: hyd
    real(wp_local)    :: time
    integer           :: nargs, i_out, n_step

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
    call nml_read(nml_file, "shmip", "out_file",     out_file)

    write(*,'(a)') "FastHydrology SHMIP driver"
    write(*,'(a,a)')      "  case    = ", trim(shmip_case)
    write(*,'(a,i6,a,i6)')"  grid    = ", nx, " x ", ny
    write(*,'(a,a)')      "  output  = ", trim(out_file)

    ! ---- Allocate ----
    allocate(H_ice(nx,ny), z_bed(nx,ny), z_sl(nx,ny))
    allocate(f_ice(nx,ny), f_grnd(nx,ny), mask(nx,ny))
    allocate(bmb_w(nx,ny), uxy_b(nx,ny), A_glen(nx,ny))
    allocate(xc(nx), yc(ny))

    ! ---- Build geometry + forcing per case ----
    call setup_case(shmip_case, nx, ny, H_ice, z_bed, dx, dy, melt_rate, xc, yc)

    z_sl   = 0.0_wp_local
    f_ice  = 1.0_wp_local
    f_grnd = 1.0_wp_local
    mask   = 1.0_wp_local
    bmb_w  = melt_rate
    uxy_b  = uxy_b_const
    A_glen = A_glen_const

    ! ---- Initialize library ----
    call hydro_init(hyd, nml_file, nx, ny)
    hyd%par%dx = dx
    hyd%par%dy = dy
    call hydro_init_state(hyd, z_bed, f_ice, f_grnd, t_start)

    write(*,'(a,i0,a)') "  method  = ", hyd%par%method,    " (0=NONE 1=BUCKET 2=K24)"
    write(*,'(a,i0)')   "  N_clos  = ", hyd%par%bucket%N_closure
    write(*,'(a,es12.3,a)') "  melt    = ", melt_rate, " m/a (water equiv.)"

    ! ---- Initialize output file ----
    call output_init(out_file, xc, yc)
    call output_step(out_file, hyd, t_start, 1)

    ! ---- Time-step loop ----
    time   = t_start
    i_out  = 1
    n_step = 0

    write(*,'(a)') ""
    write(*,'(a)') "       time      H_w_max     H_w_mean       N_max      N_mean"

    do while (time < t_end - 0.5_wp_local*dt_step)
        time = time + dt_step
        n_step = n_step + 1
        call hydro_update(hyd, H_ice, z_bed, z_sl, f_ice, f_grnd, mask, &
                          bmb_w, uxy_b, A_glen, time)

        if (mod(real(n_step,wp_local)*dt_step, dt_out) < 0.5_wp_local*dt_step .or. time >= t_end) then
            i_out = i_out + 1
            call output_step(out_file, hyd, time, i_out)
            write(*,'(f12.4, 4es13.4)') time, &
                maxval(hyd%now%H_w), sum(hyd%now%H_w)/(nx*ny), &
                maxval(hyd%now%N),   sum(hyd%now%N)/(nx*ny)
        end if
    end do

    write(*,'(a)') ""
    write(*,'(a)') "shmip: done."

contains

    subroutine setup_case(case_id, nx, ny, H_ice, z_bed, dx, dy, melt_rate, xc, yc)
        ! Build geometry, baseline forcing, and coordinate axes per SHMIP case.

        implicit none

        character(len=*), intent(IN)  :: case_id
        integer,          intent(IN)  :: nx, ny
        real(wp_local),   intent(OUT) :: H_ice(nx,ny), z_bed(nx,ny)
        real(wp_local),   intent(OUT) :: dx, dy
        real(wp_local),   intent(OUT) :: melt_rate
        real(wp_local),   intent(OUT) :: xc(nx), yc(ny)

        real(wp_local) :: m_s
        integer        :: i

        select case (trim(case_id))
            case ("A1") ; m_s = 7.93e-11_wp_local
            case ("A2") ; m_s = 1.59e-9_wp_local
            case ("A3") ; m_s = 5.79e-9_wp_local
            case ("A4") ; m_s = 2.5e-8_wp_local
            case ("A5") ; m_s = 4.5e-8_wp_local
            case ("A6") ; m_s = 5.79e-7_wp_local
            case default
                write(*,*) "shmip:: case '"//trim(case_id)//"' not implemented in this driver yet."
                write(*,*) "         Implemented: A1..A6"
                write(*,*) "         Stubs for B/C/D forthcoming."
                stop
        end select

        melt_rate = m_s * SEC_PER_YEAR              ! m/a water-equivalent

        ! A geometry: square slab, 100km x 20km, surface profile s(x) =
        ! 6*(sqrt(x+5000) - sqrt(5000)) + 1, flat bed at z=0.
        call setup_shmip_A(nx, ny, H_ice, z_bed, dx, dy, xc, yc)

    end subroutine setup_case

    subroutine setup_shmip_A(nx, ny, H_ice, z_bed, dx, dy, xc, yc)

        implicit none

        integer,        intent(IN)  :: nx, ny
        real(wp_local), intent(OUT) :: H_ice(nx,ny), z_bed(nx,ny)
        real(wp_local), intent(OUT) :: dx, dy
        real(wp_local), intent(OUT) :: xc(nx), yc(ny)

        real(wp_local) :: x, z_s, Lx, Ly
        integer        :: i, j

        Lx = 100.0e3_wp_local
        Ly =  20.0e3_wp_local
        dx = Lx / real(nx - 1, wp_local)
        dy = Ly / real(ny - 1, wp_local)

        do i = 1, nx
            xc(i) = real(i-1, wp_local) * dx
        end do
        do j = 1, ny
            yc(j) = real(j-1, wp_local) * dy
        end do

        do j = 1, ny
        do i = 1, nx
            x          = xc(i)
            z_s        = 6.0_wp_local * (sqrt(x + 5.0e3_wp_local) - sqrt(5.0e3_wp_local)) + 1.0_wp_local
            z_bed(i,j) = 0.0_wp_local
            H_ice(i,j) = max(z_s - z_bed(i,j), 0.0_wp_local)
        end do
        end do

    end subroutine setup_shmip_A

    subroutine output_init(filename, xc, yc)

        implicit none

        character(len=*), intent(IN) :: filename
        real(wp_local),   intent(IN) :: xc(:), yc(:)

        call nc_create(filename)
        call nc_write_dim(filename, "xc", x=xc, units="m", long_name="x-coordinate")
        call nc_write_dim(filename, "yc", x=yc, units="m", long_name="y-coordinate")
        call nc_write_dim(filename, "time", x=0.0_wp_local, dx=1.0_wp_local, nx=1, &
                          units="years", unlimited=.true.)

    end subroutine output_init

    subroutine output_step(filename, hyd, time, n)

        implicit none

        character(len=*), intent(IN) :: filename
        type(hydro_class), intent(IN) :: hyd
        real(wp_local),   intent(IN) :: time
        integer,          intent(IN) :: n

        call nc_write(filename, "time",  time,             dim1="time", start=[n], count=[1])
        call nc_write(filename, "H_w",   hyd%now%H_w,      dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], units="m",     long_name="Basal water layer thickness")
        call nc_write(filename, "dHwdt", hyd%now%dHwdt,    dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], units="m a-1", long_name="Rate of change of H_w")
        call nc_write(filename, "N",     hyd%now%N,        dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], units="Pa",    long_name="Effective pressure")
        call nc_write(filename, "p_w",   hyd%now%p_w,      dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], units="Pa",    long_name="Water pressure")
        call nc_write(filename, "q_x",   hyd%now%q_x,      dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], units="m2 s-1",long_name="Water flux x-component")
        call nc_write(filename, "q_y",   hyd%now%q_y,      dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], units="m2 s-1",long_name="Water flux y-component")

    end subroutine output_step

end program shmip
