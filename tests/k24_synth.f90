program k24_synth
    ! Synthetic cross-validation of the Fortran K24 against FastHydrology.jl.
    !
    ! Complements k24_compare.f90 (real Greenland data), which cannot reach
    ! several code paths: on that dataset N_inf saturates at Po on 100% of
    ! grounded cells, so the kappa blend, the drainage modes and the interior
    ! N_inf branch are all invisible there. This driver picks A_glen and the
    ! forcing so N_inf lands strictly inside (sigmat*Po, Po), and uses
    ! dx /= dy so the anisotropic corfac and the rectangular smoothing kernel
    ! are exercised too.
    !
    ! Calls calc_k24 directly rather than going through hydro_update, so
    ! everything stays in double precision and the comparison is not limited by
    ! the public API's real(sp) state arrays.
    !
    ! usage: k24_synth.x <namelist> <out.nc> [A_glen]

    use nml
    use ncio
    use fast_hydrology_k24

    implicit none

    integer,  parameter :: dp = kind(1.d0)
    real(dp), parameter :: PI = 3.14159265358979323846_dp

    ! Grid size and spacing come from the input file, so the same driver runs
    ! the small synthetic case and a full ice-sheet dataset (e.g. Thwaites 2km).
    integer  :: nx, ny

    character(len=256) :: nml_file, out_file, in_file

    real(dp) :: dx, dy
    real(dp), allocatable :: h(:,:), b(:,:), mask(:,:), mdot(:,:)
    real(dp), allocatable :: uxy_b(:,:), A_glen(:,:), kappa(:,:)
    real(dp), allocatable :: q_x(:,:), q_y(:,:), N(:,:), p_w(:,:), W(:,:), q(:,:)
    real(dp), allocatable :: xc(:), yc(:)
    real(dp), allocatable :: gsx(:,:), gsy(:,:), absgs(:,:), absg(:,:), phi0(:,:)
    integer  :: i, j

    type(k24_param_class) :: par

    call get_command_argument(1, nml_file)
    call get_command_argument(2, out_file)
    call get_command_argument(3, in_file)
    if (len_trim(in_file) == 0) in_file = "tests/k24_synth_input.nc"

    call k24_par_load(par, nml_file, "yhyd")
    ! k24_par_load leaves rho/g at its own defaults; stamp the same values the
    ! Julia constructor uses and refresh K, exactly as hydro_par_load does.
    par%water_density = 1000.0_dp
    par%ice_density   =  917.0_dp
    par%gravity       =    9.81_dp
    call k24_finalize_par(par)

    ! ---- fields read from the shared input file so both implementations see
    !      bit-identical data (see k24_synth_gen.jl / thw_export.jl) ----
    nx = nc_size(in_file, "xc")
    ny = nc_size(in_file, "yc")

    allocate(xc(nx), yc(ny))
    allocate(h(nx,ny), b(nx,ny), mask(nx,ny), mdot(nx,ny))
    allocate(uxy_b(nx,ny), A_glen(nx,ny), kappa(nx,ny))
    allocate(q_x(nx,ny), q_y(nx,ny), N(nx,ny), p_w(nx,ny), W(nx,ny), q(nx,ny))
    allocate(gsx(nx,ny), gsy(nx,ny), absgs(nx,ny), absg(nx,ny), phi0(nx,ny))

    call nc_read(in_file, "xc",   xc)
    call nc_read(in_file, "yc",   yc)
    dx = xc(2) - xc(1)
    dy = yc(2) - yc(1)

    call nc_read(in_file, "h",    h)
    call nc_read(in_file, "b",    b)
    call nc_read(in_file, "mask", mask)
    call nc_read(in_file, "vb",   uxy_b)
    call nc_read(in_file, "A",    A_glen)
    call nc_read(in_file, "mdot", mdot)

    call initialize_kappa(kappa, b, par%substrate_type)

    q = 0.0_dp
    N = 0.0_dp

    call calc_k24(q_x, q_y, N, p_w, W, q, &
                  h, b, mask, mdot, uxy_b, A_glen, kappa, dx, dy, par, &
                  gsx, gsy, absgs, absg, phi0)

    call nc_create(out_file)
    call nc_write_dim(out_file, "xc", x=xc, units="m")
    call nc_write_dim(out_file, "yc", x=yc, units="m")
    call nc_write(out_file, "W",    W,    dim1="xc", dim2="yc")
    call nc_write(out_file, "N",    N,    dim1="xc", dim2="yc")
    call nc_write(out_file, "q",    q,    dim1="xc", dim2="yc")
    call nc_write(out_file, "p_w",  p_w,  dim1="xc", dim2="yc")
    call nc_write(out_file, "mask", mask, dim1="xc", dim2="yc")
    call nc_write(out_file, "gsx",   gsx,   dim1="xc", dim2="yc")
    call nc_write(out_file, "gsy",   gsy,   dim1="xc", dim2="yc")
    call nc_write(out_file, "absgs", absgs, dim1="xc", dim2="yc")
    call nc_write(out_file, "absg",  absg,  dim1="xc", dim2="yc")
    call nc_write(out_file, "phi0",  phi0,  dim1="xc", dim2="yc")

    write(*,'(a,i0,a,i0,a,f9.2,a,f9.2)') "k24_synth: ", nx, " x ", ny, "  dx=", dx, " dy=", dy
    write(*,'(a,a)') "k24_synth: wrote ", trim(out_file)

end program k24_synth
