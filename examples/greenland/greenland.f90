program greenland
    ! Load a Yelmo Greenland-16km restart file, run the FastHydrology
    ! library to a steady-state H_w field with either BUCKET or K24, and
    ! write a NetCDF output. Configured via a namelist that points at the
    ! restart and selects par%method.

    use nml
    use ncio
    use fast_hydrology

    implicit none

    integer, parameter :: wp_local = kind(1.0)

    character(len=256) :: nml_file
    character(len=256) :: restart_file, out_file
    real(wp_local)     :: t_end, dt_step, dt_out

    real(wp_local), allocatable :: xc(:), yc(:)
    real(wp_local), allocatable :: H_ice(:,:), z_bed(:,:), z_srf(:,:), z_sl(:,:)
    real(wp_local), allocatable :: f_ice(:,:), f_grnd(:,:), mask(:,:)
    real(wp_local), allocatable :: bmb_grnd(:,:), uxy_b(:,:), A_glen(:,:)
    real(wp_local), allocatable :: bmb_w(:,:)

    type(hydro_class) :: hyd
    real(wp_local)    :: time, dx_km, dy_km
    integer           :: nx, ny, nargs, n_step, i_out, n_out

    nargs = command_argument_count()
    if (nargs < 1) then
        write(*,*) "usage: greenland.x <namelist>"
        stop
    end if
    call get_command_argument(1, nml_file)

    call nml_read(nml_file, "greenland", "restart_file", restart_file)
    call nml_read(nml_file, "greenland", "out_file",     out_file)
    call nml_read(nml_file, "greenland", "t_end",        t_end)
    call nml_read(nml_file, "greenland", "dt_step",      dt_step)
    call nml_read(nml_file, "greenland", "dt_out",       dt_out)

    write(*,'(a)') "FastHydrology Greenland example"
    write(*,'(a,a)') "  restart = ", trim(restart_file)
    write(*,'(a,a)') "  output  = ", trim(out_file)

    ! ---- Read grid dimensions ----
    nx = nc_size(restart_file, "xc")
    ny = nc_size(restart_file, "yc")
    write(*,'(a,i0,a,i0)') "  grid    = ", nx, " x ", ny

    allocate(xc(nx), yc(ny))
    allocate(H_ice(nx,ny), z_bed(nx,ny), z_srf(nx,ny), z_sl(nx,ny))
    allocate(f_ice(nx,ny), f_grnd(nx,ny), mask(nx,ny))
    allocate(bmb_grnd(nx,ny), uxy_b(nx,ny), A_glen(nx,ny))
    allocate(bmb_w(nx,ny))

    call nc_read(restart_file, "xc",       xc)
    call nc_read(restart_file, "yc",       yc)
    ! 2D fields from a single-time restart (start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "H_ice",    H_ice,    start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "z_bed",    z_bed,    start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "z_srf",    z_srf,    start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "z_sl",     z_sl,     start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "f_ice",    f_ice,    start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "f_grnd",   f_grnd,   start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "uxy_b",    uxy_b,    start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "ATT_bar",  A_glen,   start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "bmb_grnd", bmb_grnd, start=[1,1,1], count=[nx,ny,1])

    ! Convert bmb_grnd (ice-equivalent, sign convention: positive = melt) to
    ! water-equivalent (m/a) for the library: bmb_w = -bmb_grnd * (rho_ice/rho_w).
    ! Yelmo's sign convention has bmb_grnd negative when melting; flip it.
    bmb_w = -bmb_grnd * (917.0_wp_local / 1000.0_wp_local)

    ! Build mask: K24 / BUCKET active where grounded ice exists.
    where (f_grnd > 0.0_wp_local .and. f_ice > 0.0_wp_local)
        mask = 1.0_wp_local
    elsewhere
        mask = 0.0_wp_local
    end where

    ! Grid spacing from the xc axis. yelmo restart stores xc in km.
    dx_km = xc(2) - xc(1)
    dy_km = yc(2) - yc(1)
    write(*,'(a,f7.3,a,f7.3,a)') "  dx,dy   = ", dx_km, " km, ", dy_km, " km"

    ! ---- Initialize library ----
    call hydro_init(hyd, nml_file, nx, ny)
    hyd%par%dx = dx_km * 1000.0_wp_local
    hyd%par%dy = dy_km * 1000.0_wp_local
    call hydro_init_state(hyd, z_bed, f_ice, f_grnd, 0.0_wp_local)

    write(*,'(a,i0,a)') "  method  = ", hyd%par%method, " (0=NONE 1=BUCKET 2=K24)"
    write(*,'(a,i0)')   "  N_clos  = ", hyd%par%bucket%N_closure
    write(*,'(a,i0)')   "  mask_bc = ", hyd%par%mask_bc
    write(*,'(a,es12.3,a)') "  bmb_max = ", maxval(bmb_w), " m/a (water equiv.)"

    ! ---- Output ----
    call output_init(out_file, xc, yc)
    call output_step(out_file, hyd, 0.0_wp_local, 1)

    ! ---- Time-step loop ----
    time   = 0.0_wp_local
    n_step = 0
    i_out  = 1
    n_out  = max(1, nint(dt_out / dt_step))

    write(*,'(a)') ""
    write(*,'(a)') "       time      H_w_max     H_w_mean       N_max      N_mean"

    do while (time < t_end - 0.5_wp_local*dt_step)
        time   = time + dt_step
        n_step = n_step + 1
        call hydro_update(hyd, H_ice, z_bed, z_sl, f_ice, f_grnd, mask, &
                          bmb_w, uxy_b, A_glen, time)

        if (mod(n_step, n_out) == 0 .or. time >= t_end) then
            i_out = i_out + 1
            call output_step(out_file, hyd, time, i_out)
            write(*,'(f12.4, 4es13.4)') time, &
                maxval(hyd%now%H_w), sum(hyd%now%H_w * mask) / max(1.0_wp_local, sum(mask)), &
                maxval(hyd%now%N),   sum(hyd%now%N   * mask) / max(1.0_wp_local, sum(mask))
        end if
    end do

    write(*,'(a)') ""
    write(*,'(a)') "greenland: done."

contains

    subroutine output_init(filename, xc, yc)

        implicit none

        character(len=*), intent(IN) :: filename
        real(wp_local),   intent(IN) :: xc(:), yc(:)

        call nc_create(filename)
        call nc_write_dim(filename, "xc", x=xc, units="km", long_name="x-coordinate")
        call nc_write_dim(filename, "yc", x=yc, units="km", long_name="y-coordinate")
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

end program greenland
