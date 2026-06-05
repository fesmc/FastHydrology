program shmip
    ! Stand-alone SHMIP driver for the FastHydrology library.
    !
    ! Cases:
    !   A1..A6  square slab + uniform distributed melt at SHMIP-spec rate
    !   B1..B5  square slab + moulin point sources (background = A1 rate)
    !   D1..D5  A5 setup + seasonal melt forcing
    ! (SHMIP C is not currently supported: real C resolution needs dt ~ minutes,
    !  which is out of scope for this driver.)
    !
    ! Output written to a NetCDF file with W_til, dW_til_dt, overflow, W, N, p_w, q_x, q_y
    ! per output step.

    use nml
    use ncio
    use fast_hydrology

    implicit none

    integer, parameter :: wp_local = kind(1.0)
    integer, parameter :: dp = kind(1.d0)

    real(wp_local), parameter :: SEC_PER_YEAR = 31536000.0_wp_local
    real(wp_local), parameter :: PI           = 4.0_wp_local * atan(1.0_wp_local)
    integer, parameter        :: N_MOULIN     = 10

    type case_state_t
        character(len=16) :: case_id
        real(wp_local)    :: m_s_background           ! [m/s] distributed melt
        real(wp_local)    :: Q_moulin_total           ! [m3/s] total moulin input
        real(wp_local)    :: seasonal_amp             ! [-] for D-cases
        integer           :: ix_moulin(N_MOULIN)
        integer           :: jy_moulin(N_MOULIN)
        logical           :: has_moulins
        logical           :: is_seasonal
        real(wp_local)    :: dx, dy
    end type

    character(len=128) :: nml_file
    character(len=128) :: out_file
    character(len=16)  :: shmip_case
    integer            :: nx, ny
    real(wp_local)     :: t_start, t_end, dt_out, dt_step
    real(wp_local)     :: A_glen_const, uxy_b_const

    real(wp_local), allocatable :: H_ice(:,:), z_bed(:,:), z_sl(:,:)
    real(wp_local), allocatable :: f_ice(:,:), f_grnd(:,:), mask(:,:)
    real(wp_local), allocatable :: mdot(:,:), uxy_b(:,:), A_glen(:,:)
    real(wp_local), allocatable :: xc(:), yc(:)

    type(hydro_class) :: hyd
    type(case_state_t) :: cs
    real(wp_local)    :: time
    integer           :: nargs, i_out, n_step

    nargs = command_argument_count()
    if (nargs < 1) then
        write(*,*) "usage: shmip.x <namelist>"
        stop
    end if
    call get_command_argument(1, nml_file)

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
    write(*,'(a,a)')       "  case    = ", trim(shmip_case)
    write(*,'(a,i6,a,i6)') "  grid    = ", nx, " x ", ny
    write(*,'(a,a)')       "  output  = ", trim(out_file)

    allocate(H_ice(nx,ny), z_bed(nx,ny), z_sl(nx,ny))
    allocate(f_ice(nx,ny), f_grnd(nx,ny), mask(nx,ny))
    allocate(mdot(nx,ny), uxy_b(nx,ny), A_glen(nx,ny))
    allocate(xc(nx), yc(ny))

    call setup_case(shmip_case, nx, ny, H_ice, z_bed, xc, yc, cs)

    z_sl   = 0.0_wp_local
    f_ice  = 1.0_wp_local
    f_grnd = 1.0_wp_local
    mask   = 1.0_wp_local
    uxy_b  = uxy_b_const
    A_glen = A_glen_const

    call update_forcing(cs, t_start, mdot)

    call hydro_init(hyd, nml_file, nx, ny)
    hyd%par%dx = cs%dx
    hyd%par%dy = cs%dy
    call hydro_init_state(hyd, z_bed, f_ice, f_grnd, t_start)

    write(*,'(a,i0,a)')      "  method_til       = ", hyd%par%method_til, " (0=BUCKET 1=EXTERNAL)"
    write(*,'(a,i0,a)')      "  method_transport = ", hyd%par%method_transport, " (0=NONE 1=K24)"
    write(*,'(a,i0)')        "  N_clos           = ", hyd%par%bucket%N_closure
    write(*,'(a,es12.3,a)')  "  bg melt          = ", cs%m_s_background, " m/s"
    if (cs%has_moulins) then
        write(*,'(a,es12.3,a,i0,a)') "  moulins = ", cs%Q_moulin_total, " m3/s total over ", N_MOULIN, " points"
    end if

    call output_init(out_file, xc, yc)
    call output_step(out_file, hyd, t_start, 1)

    time   = t_start
    i_out  = 1
    n_step = 0

    write(*,'(a)') ""
    write(*,'(a)') "       time    W_til_max   W_til_mean       N_max      N_mean"

    do while (time < t_end - 0.5_wp_local*dt_step)
        time = time + dt_step
        n_step = n_step + 1

        call update_forcing(cs, time, mdot)
        call hydro_update(hyd, H_ice, z_bed, z_sl, f_ice, f_grnd, mask, &
                          mdot, uxy_b, A_glen, time)

        if (mod(real(n_step,wp_local)*dt_step, dt_out) < 0.5_wp_local*dt_step .or. time >= t_end) then
            i_out = i_out + 1
            call output_step(out_file, hyd, time, i_out)
            write(*,'(f12.4, 4es13.4)') time, &
                maxval(hyd%now%W_til), sum(hyd%now%W_til)/(nx*ny), &
                maxval(hyd%now%N),     sum(hyd%now%N)/(nx*ny)
        end if
    end do

    write(*,'(a)') ""
    write(*,'(a)') "shmip: done."

contains

    subroutine setup_case(case_id, nx, ny, H_ice, z_bed, xc, yc, cs)
        ! Build geometry, coordinate axes, and the persistent case_state
        ! struct (moulin positions, melt rates, time-variation flags).

        implicit none

        character(len=*),    intent(IN)  :: case_id
        integer,             intent(IN)  :: nx, ny
        real(wp_local),      intent(OUT) :: H_ice(nx,ny), z_bed(nx,ny)
        real(wp_local),      intent(OUT) :: xc(nx), yc(ny)
        type(case_state_t),  intent(OUT) :: cs

        cs%case_id        = case_id
        cs%m_s_background = 0.0_wp_local
        cs%Q_moulin_total = 0.0_wp_local
        cs%seasonal_amp   = 0.0_wp_local
        cs%has_moulins    = .false.
        cs%is_seasonal    = .false.

        select case (trim(case_id))
            case ("A1") ; cs%m_s_background = 7.93e-11_wp_local
            case ("A2") ; cs%m_s_background = 1.59e-9_wp_local
            case ("A3") ; cs%m_s_background = 5.79e-9_wp_local
            case ("A4") ; cs%m_s_background = 2.5e-8_wp_local
            case ("A5") ; cs%m_s_background = 4.5e-8_wp_local
            case ("A6") ; cs%m_s_background = 5.79e-7_wp_local
            case ("B1") ; cs%m_s_background = 7.93e-11_wp_local ; cs%Q_moulin_total = 0.0463_wp_local ; cs%has_moulins = .true.
            case ("B2") ; cs%m_s_background = 7.93e-11_wp_local ; cs%Q_moulin_total = 0.926_wp_local  ; cs%has_moulins = .true.
            case ("B3") ; cs%m_s_background = 7.93e-11_wp_local ; cs%Q_moulin_total = 4.63_wp_local   ; cs%has_moulins = .true.
            case ("B4") ; cs%m_s_background = 7.93e-11_wp_local ; cs%Q_moulin_total = 23.15_wp_local  ; cs%has_moulins = .true.
            case ("B5") ; cs%m_s_background = 7.93e-11_wp_local ; cs%Q_moulin_total = 46.3_wp_local   ; cs%has_moulins = .true.
            ! D-cases: A5 setup with seasonal melt scaling.
            case ("D1") ; cs%m_s_background = 4.5e-8_wp_local ; cs%is_seasonal = .true. ; cs%seasonal_amp = 0.25_wp_local
            case ("D2") ; cs%m_s_background = 4.5e-8_wp_local ; cs%is_seasonal = .true. ; cs%seasonal_amp = 0.5_wp_local
            case ("D3") ; cs%m_s_background = 4.5e-8_wp_local ; cs%is_seasonal = .true. ; cs%seasonal_amp = 1.0_wp_local
            case ("D4") ; cs%m_s_background = 4.5e-8_wp_local ; cs%is_seasonal = .true. ; cs%seasonal_amp = 2.0_wp_local
            case ("D5") ; cs%m_s_background = 4.5e-8_wp_local ; cs%is_seasonal = .true. ; cs%seasonal_amp = 4.0_wp_local
            case default
                write(*,*) "shmip:: case '"//trim(case_id)//"' not recognized."
                write(*,*) "         Supported: A1..A6, B1..B5, D1..D5"
                stop
        end select

        ! All SHMIP A/B/C/D cases share the square-slab geometry.
        call setup_shmip_A(nx, ny, H_ice, z_bed, cs%dx, cs%dy, xc, yc)

        if (cs%has_moulins) call init_moulin_positions(nx, ny, cs)

    end subroutine setup_case

    subroutine update_forcing(cs, time, mdot)
        ! Fill mdot (m/a, water-equivalent) from the case state at the
        ! given time. Distributed background scaled per case; moulin
        ! contributions added as point sources.

        implicit none

        type(case_state_t), intent(IN)  :: cs
        real(wp_local),     intent(IN)  :: time
        real(wp_local),     intent(OUT) :: mdot(:,:)

        real(wp_local) :: m_s, Q, q_per_cell, scale
        integer        :: k, i, j, nx, ny

        nx = size(mdot,1)
        ny = size(mdot,2)

        ! Distributed background, optionally seasonal
        if (cs%is_seasonal) then
            ! Smooth seasonal cycle peaking at mid-summer (t mod 1 = 0.5).
            scale = max(0.0_wp_local, &
                        1.0_wp_local + cs%seasonal_amp * sin(2.0_wp_local*PI*(time - 0.25_wp_local)))
            m_s = cs%m_s_background * scale
        else
            m_s = cs%m_s_background
        end if

        mdot = m_s * SEC_PER_YEAR

        ! Moulin point sources
        if (cs%has_moulins) then
            Q          = cs%Q_moulin_total
            q_per_cell = Q / real(N_MOULIN, wp_local) / (cs%dx * cs%dy) * SEC_PER_YEAR
            do k = 1, N_MOULIN
                i = cs%ix_moulin(k)
                j = cs%jy_moulin(k)
                mdot(i,j) = mdot(i,j) + q_per_cell
            end do
        end if

    end subroutine update_forcing

    subroutine init_moulin_positions(nx, ny, cs)
        ! Pseudo-random moulin positions inside the SHMIP-B active domain:
        !   5 km < x < 95 km    (excludes the 5 km halo at each end of the slab)
        !   1 km < y < 19 km
        ! Generated with a Park-Miller LCG seeded at SEED for reproducibility.
        ! This follows the SHMIP convention of fixed-seed random placement
        ! but is NOT bit-equivalent to the protocol's reference positions --
        ! patch those in if a specific intercomparison is needed.

        implicit none

        integer,            intent(IN)    :: nx, ny
        type(case_state_t), intent(INOUT) :: cs

        integer(8), parameter :: SEED   = 4242_8
        integer(8), parameter :: PM_A   = 48271_8
        integer(8), parameter :: PM_M   = 2147483647_8

        integer(8)     :: state
        integer        :: k
        real(wp_local) :: u, x_frac, y_frac, x, y, x_lo, x_hi, y_lo, y_hi
        real(wp_local) :: Lx, Ly

        Lx   = (nx - 1) * cs%dx
        Ly   = (ny - 1) * cs%dy
        x_lo =  5.0e3_wp_local ; x_hi = max(x_lo + 1.0_wp_local, Lx - 5.0e3_wp_local)
        y_lo =  1.0e3_wp_local ; y_hi = max(y_lo + 1.0_wp_local, Ly - 1.0e3_wp_local)

        state = SEED
        do k = 1, N_MOULIN
            state = mod(state * PM_A, PM_M) ; u = real(state, wp_local) / real(PM_M, wp_local)
            x_frac = u
            state = mod(state * PM_A, PM_M) ; u = real(state, wp_local) / real(PM_M, wp_local)
            y_frac = u

            x = x_lo + x_frac * (x_hi - x_lo)
            y = y_lo + y_frac * (y_hi - y_lo)

            cs%ix_moulin(k) = max(2, min(nx-1, nint(x / cs%dx) + 1))
            cs%jy_moulin(k) = max(2, min(ny-1, nint(y / cs%dy) + 1))
        end do

    end subroutine init_moulin_positions

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

        call nc_write(filename, "time",     time,               dim1="time", start=[n], count=[1])
        call nc_write(filename, "W_til",    hyd%now%W_til,      dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], units="m",     long_name="Till water storage thickness")
        call nc_write(filename, "dW_til_dt",hyd%now%dW_til_dt,  dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], units="m a-1", long_name="Rate of change of W_til")
        call nc_write(filename, "overflow", hyd%now%overflow,   dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], units="m a-1", long_name="Till-saturation overflow rate")
        call nc_write(filename, "W",        hyd%now%W,          dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], units="m",     long_name="Distributed sheet water-layer thickness (K24)")
        call nc_write(filename, "N",        hyd%now%N,          dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], units="Pa",    long_name="Effective pressure")
        call nc_write(filename, "p_w",      hyd%now%p_w,        dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], units="Pa",    long_name="Water pressure")
        call nc_write(filename, "q_x",      hyd%now%q_x,        dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], units="m2 s-1",long_name="Water flux x-component")
        call nc_write(filename, "q_y",      hyd%now%q_y,        dim1="xc", dim2="yc", dim3="time", &
                      start=[1,1,n], units="m2 s-1",long_name="Water flux y-component")

    end subroutine output_step

end program shmip
