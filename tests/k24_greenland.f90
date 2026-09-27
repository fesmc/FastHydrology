program k24_compare
    ! Cross-validation driver: run the Fortran K24 on the Greenland-16km
    ! restart with EXACTLY the inputs examples/greenland/greenland.jl feeds
    ! FastHydrology.jl, so the two can be compared field by field.
    !
    ! Alignment with greenland.jl:
    !   * mdot = -bmb_grnd. greenland.jl sets mdot_mass = -bmb_grnd*1000
    !     [kg/m2/s]; the Julia routing seeds psi_out with mdot_mass*dx*dy/rho_w,
    !     i.e. a volume rate of -bmb_grnd. This library's mdot IS that volume
    !     rate, so no conversion.
    !   * method_til = 0 (TIL_NONE) so the K24 source is mdot itself, not the
    !     bucket overflow. greenland.jl has no bucket.
    !   * one step, matching run!(SteadyStateSimulation) = one
    !     update_steady_state!.
    !   * mask = f_ice*f_grnd > 0, h = H_ice, b = z_bed, as greenland.jl.

    use nml
    use ncio
    use fast_hydrology

    implicit none

    integer, parameter :: wp_local = kind(1.0)

    character(len=256) :: nml_file, restart_file, out_file, scale_arg
    real(wp_local)     :: mdot_scale

    real(wp_local), allocatable :: xc(:), yc(:)
    real(wp_local), allocatable :: H_ice(:,:), z_bed(:,:), z_sl(:,:)
    real(wp_local), allocatable :: f_ice(:,:), f_grnd(:,:), mask(:,:)
    real(wp_local), allocatable :: bmb_grnd(:,:), uxy_b(:,:), A_glen(:,:)
    real(wp_local), allocatable :: mdot(:,:)

    type(hydro_class) :: hyd
    real(wp_local)    :: dx_km, dy_km
    integer           :: nx, ny

    call get_command_argument(1, nml_file)
    call get_command_argument(2, out_file)
    if (len_trim(out_file) == 0) out_file = "k24_compare_fortran.nc"

    ! Optional 3rd argument: multiplier applied to mdot, so the same driver can
    ! be run with the raw restart units or with a physically-scaled forcing.
    call get_command_argument(3, scale_arg)
    if (len_trim(scale_arg) == 0) then
        mdot_scale = 1.0_wp_local
    else
        read(scale_arg, *) mdot_scale
    end if

    restart_file = "input/GRL-16KM_yelmo_restart.nc"

    nx = nc_size(restart_file, "xc")
    ny = nc_size(restart_file, "yc")

    allocate(xc(nx), yc(ny))
    allocate(H_ice(nx,ny), z_bed(nx,ny), z_sl(nx,ny))
    allocate(f_ice(nx,ny), f_grnd(nx,ny), mask(nx,ny))
    allocate(bmb_grnd(nx,ny), uxy_b(nx,ny), A_glen(nx,ny), mdot(nx,ny))

    call nc_read(restart_file, "xc",       xc)
    call nc_read(restart_file, "yc",       yc)
    call nc_read(restart_file, "H_ice",    H_ice,    start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "z_bed",    z_bed,    start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "z_sl",     z_sl,     start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "f_ice",    f_ice,    start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "f_grnd",   f_grnd,   start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "uxy_b",    uxy_b,    start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "ATT_bar",  A_glen,   start=[1,1,1], count=[nx,ny,1])
    call nc_read(restart_file, "bmb_grnd", bmb_grnd, start=[1,1,1], count=[nx,ny,1])

    ! Match greenland.jl's volume source rate exactly (see header).
    mdot = -bmb_grnd * mdot_scale

    where (f_grnd > 0.0_wp_local .and. f_ice > 0.0_wp_local)
        mask = 1.0_wp_local
    elsewhere
        mask = 0.0_wp_local
    end where

    dx_km = xc(2) - xc(1)
    dy_km = yc(2) - yc(1)

    call hydro_init(hyd, nml_file, nx, ny, &
                    dx_km * 1000.0_wp_local, dy_km * 1000.0_wp_local)
    call hydro_init_state(hyd, H_ice, z_bed, f_ice, f_grnd, 0.0_wp_local)

    write(*,'(a,i0,a,i0)') "k24_compare: grid ", nx, " x ", ny
    write(*,'(a,es14.6)') "  mdot_scale = ", mdot_scale
    write(*,'(a,es14.6)') "  mdot_max   = ", maxval(mdot)
    write(*,'(a,i0,a,i0)') "  method_til = ", hyd%par%method_til, &
                           "  method_transport = ", hyd%par%method_transport

    ! One step. Any dt > 0 works: with TIL_NONE the source is mdot and K24 is
    ! purely diagnostic, so the result does not depend on dt.
    call hydro_update(hyd, H_ice, z_bed, z_sl, f_ice, f_grnd, mask, &
                      mdot, uxy_b, A_glen, 1.0_wp_local)

    call nc_create(out_file)
    call nc_write_dim(out_file, "xc", x=xc, units="km")
    call nc_write_dim(out_file, "yc", x=yc, units="km")

    call nc_write(out_file, "W",    hyd%now%W,    dim1="xc", dim2="yc")
    call nc_write(out_file, "N",    hyd%now%N,    dim1="xc", dim2="yc")
    call nc_write(out_file, "q",    hyd%now%q,    dim1="xc", dim2="yc")
    call nc_write(out_file, "p_w",  hyd%now%p_w,  dim1="xc", dim2="yc")
    call nc_write(out_file, "q_x",  hyd%now%q_x,  dim1="xc", dim2="yc")
    call nc_write(out_file, "q_y",  hyd%now%q_y,  dim1="xc", dim2="yc")
    call nc_write(out_file, "mask", mask,         dim1="xc", dim2="yc")

    write(*,'(a,a)') "k24_compare: wrote ", trim(out_file)

end program k24_compare
