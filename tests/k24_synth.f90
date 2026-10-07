program k24_synth
    ! Synthetic cross-validation of the Fortran K24 against FastHydrology.jl.
    !
    ! Complements k24_greenland.f90 (real Greenland data), which cannot reach
    ! several code paths: on that dataset N_inf saturates at Po on 100% of
    ! grounded cells, so the kappa blend, the drainage modes and the interior
    ! N_inf branch are all invisible there. This driver picks A_glen and the
    ! forcing so N_inf lands strictly inside (sigmat*Po, Po), and uses
    ! dx /= dy so the anisotropic corfac, routing weights and the rectangular
    ! smoothing kernel are exercised too.
    !
    ! Calls calc_k24 directly rather than going through hydro_update, so
    ! everything stays in double precision.
    !
    ! The water source is built from terms, from the input file's volume rate
    ! `mdot` [m/s]: G = mdot*rho_w*L_w [W/m2], q_T = qtfrac*G, i_eb =
    ! iebfrac*mdot*rho_w [kg/m2/s]. The fields some options need are derived
    ! from vb the same way on both sides: ux_b = vb, uy_b = 0.5*vb (C-grid),
    ! tau_b = 2e4 + 1e10*vb, c_till = reg_coulomb_c_till*(1 + 1e5*vb).
    !
    ! usage: k24_synth.x <namelist> <out.nc> [input.nc] [key=value ...]
    !   keys read here: qtfrac=<float> iebfrac=<float> (default 0);
    !   ubfac=<float>: also write N_ub, N for the sliding speed ubfac*vb with the routing of
    !   the solve above held (k24_N_from_ub); every other key=value (the Julia twin's options)
    !   is ignored.

    use nml
    use ncio
    use fast_hydrology_k24

    implicit none

    integer,  parameter :: dp = kind(1.d0)

    integer  :: nx, ny, nargs, ia, ieq

    character(len=256) :: nml_file, out_file, in_file, arg

    real(dp) :: dx, dy, qtfrac, iebfrac, ubfac
    real(dp), allocatable :: N_ub(:,:), uxy_b_new(:,:)
    real(dp), allocatable :: h(:,:), b(:,:), mask(:,:), mdot(:,:)
    real(dp), allocatable :: uxy_b(:,:), A_glen(:,:), kappa(:,:)
    real(dp), allocatable :: G(:,:), q_T(:,:), i_eb(:,:), ux_b(:,:), uy_b(:,:), taub(:,:), c_till(:,:)
    real(dp), allocatable :: q_x(:,:), q_y(:,:), N(:,:), p_w(:,:), W(:,:), q(:,:), Q_b(:,:), Q_diss(:,:), C_frz(:,:)
    real(dp), allocatable :: xc(:), yc(:)
    real(dp), allocatable :: gsx(:,:), gsy(:,:), absgs(:,:), absg(:,:), phi0(:,:)

    type(k24_param_class) :: par

    call get_command_argument(1, nml_file)
    call get_command_argument(2, out_file)
    in_file = "tests/k24_synth_input.nc"
    qtfrac  = 0.0_dp
    iebfrac = 0.0_dp
    ubfac   = -1.0_dp
    nargs = command_argument_count()
    do ia = 3, nargs
        call get_command_argument(ia, arg)
        ieq = index(arg, "=")
        if (ieq == 0) then
            if (ia == 3) in_file = arg
            cycle
        end if
        if (arg(1:ieq-1) == "qtfrac")  read(arg(ieq+1:), *) qtfrac
        if (arg(1:ieq-1) == "iebfrac") read(arg(ieq+1:), *) iebfrac
        if (arg(1:ieq-1) == "ubfac")   read(arg(ieq+1:), *) ubfac
    end do

    call k24_par_load(par, nml_file, "yhyd")
    ! k24_par_load leaves rho/g at its own defaults; stamp the same values the
    ! Julia constructor uses and refresh the derived ones, as hydro_par_load does.
    par%water_density = 1000.0_dp
    par%ice_density   =  917.0_dp
    par%gravity       =    9.81_dp
    call k24_finalize_par(par)

    nx = nc_size(in_file, "xc")
    ny = nc_size(in_file, "yc")

    allocate(xc(nx), yc(ny))
    allocate(h(nx,ny), b(nx,ny), mask(nx,ny), mdot(nx,ny))
    allocate(uxy_b(nx,ny), A_glen(nx,ny), kappa(nx,ny))
    allocate(G(nx,ny), q_T(nx,ny), i_eb(nx,ny), ux_b(nx,ny), uy_b(nx,ny), taub(nx,ny), c_till(nx,ny))
    allocate(q_x(nx,ny), q_y(nx,ny), N(nx,ny), p_w(nx,ny), W(nx,ny), q(nx,ny), Q_b(nx,ny), Q_diss(nx,ny), C_frz(nx,ny))
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

    ! Same operations, in the same order, as tests/k24_synth.jl.
    G      = mdot * 1000.0_dp * 3.34e5_dp
    q_T    = qtfrac * G
    i_eb   = iebfrac * mdot * 1000.0_dp
    ux_b   = uxy_b
    uy_b   = 0.5_dp * uxy_b
    taub   = 2.0e4_dp + uxy_b * 1.0e10_dp
    c_till = par%reg_coulomb_c_till * (1.0_dp + uxy_b * 1.0e5_dp)

    call initialize_kappa(kappa, b, par%substrate_type)

    q = 0.0_dp
    N = 0.0_dp

    call calc_k24(q_x, q_y, N, p_w, W, q, Q_b, Q_diss, &
                  h, b, mask, G, q_T, i_eb, uxy_b, A_glen, kappa, dx, dy, par, &
                  ux_b=ux_b, uy_b=uy_b, tau_b_in=taub, c_till_in=c_till, C_frz=C_frz, &
                  gsx_out=gsx, gsy_out=gsy, absgs_out=absgs, absg_out=absg, phi0_out=phi0)

    ! N for a new sliding speed with the routing (q, absg, phi0) held, as a host's velocity
    ! iteration calls it
    if (ubfac >= 0.0_dp) then
        allocate(N_ub(nx,ny), uxy_b_new(nx,ny))
        N_ub = N
        uxy_b_new = ubfac * uxy_b
        call k24_N_from_ub(N_ub, q, absg, phi0, mask, uxy_b_new, A_glen, kappa, h, par)
    end if

    call nc_create(out_file)
    call nc_write_dim(out_file, "xc", x=xc, units="m")
    call nc_write_dim(out_file, "yc", x=yc, units="m")
    call nc_write(out_file, "W",      W,      dim1="xc", dim2="yc")
    call nc_write(out_file, "N",      N,      dim1="xc", dim2="yc")
    call nc_write(out_file, "q",      q,      dim1="xc", dim2="yc")
    call nc_write(out_file, "p_w",    p_w,    dim1="xc", dim2="yc")
    call nc_write(out_file, "Q_b",    Q_b,    dim1="xc", dim2="yc")
    call nc_write(out_file, "Q_diss", Q_diss, dim1="xc", dim2="yc")
    call nc_write(out_file, "C_frz",  C_frz,  dim1="xc", dim2="yc")
    call nc_write(out_file, "mask",   mask,   dim1="xc", dim2="yc")
    call nc_write(out_file, "gsx",    gsx,    dim1="xc", dim2="yc")
    call nc_write(out_file, "gsy",    gsy,    dim1="xc", dim2="yc")
    call nc_write(out_file, "absgs",  absgs,  dim1="xc", dim2="yc")
    call nc_write(out_file, "absg",   absg,   dim1="xc", dim2="yc")
    call nc_write(out_file, "phi0",   phi0,   dim1="xc", dim2="yc")
    if (allocated(N_ub)) call nc_write(out_file, "N_ub", N_ub, dim1="xc", dim2="yc")

    write(*,'(a,i0,a,i0,a,f9.2,a,f9.2)') "k24_synth: ", nx, " x ", ny, "  dx=", dx, " dy=", dy
    write(*,'(a,a)') "k24_synth: wrote ", trim(out_file)

end program k24_synth
