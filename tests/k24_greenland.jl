# Cross-validation reference: run FastHydrology.jl's kazmierczak2024 model on
# the Greenland-16km restart with the same inputs k24_compare.f90 feeds the
# Fortran library, and dump the result for a field-by-field comparison.
#
# Uses ArrayHydroGrid (plain arrays) rather than OGRectHydroGrid: the
# grid-interface DEFAULT implementations -- edge-clamped central differences,
# replicate-padded cached FFT convolution -- are exactly what the Fortran port
# mirrors, so this is the direct reference.
#
# Usage: julia k24_compare.jl <out.nc> [key=value ...]
#   wthick=darcy|laminar|areal   grad=mean|local   longcoup=<float>
#   psi=recursive|iterative|topological            dissip=true|false
#   sliding=none|weertman|powerplastic|regcoulomb
#   substrate=hard|soft|mixed

using NCDatasets
using FastHydrology

const RESTART = "input/GRL-16KM_yelmo_restart.nc"
out_path = length(ARGS) >= 1 ? ARGS[1] : "k24_compare_julia.nc"

cfg = Dict{String,String}()
for a in ARGS[2:end]
    k, v = split(a, "="; limit = 2)
    cfg[k] = v
end
get_cfg(k, d) = get(cfg, k, d)

ds = NCDataset(RESTART)

x = ds["xc"][:] .* 1000.0     # km -> m
y = ds["yc"][:] .* 1000.0
Nx, Ny = length(x), length(y)

t = 1
f_ice  = ds["f_ice"][:, :, t]
f_grnd = ds["f_grnd"][:, :, t]
mask   = Float64.((f_ice .> 0.0) .& (f_grnd .> 0.0))

h        = Float64.(ds["H_ice"][:, :, t])
b        = Float64.(ds["z_bed"][:, :, t])
abs_v_b  = Float64.(ds["uxy_b"][:, :, t])
A_visc   = Float64.(ds["ATT_bar"][:, :, t])
bmb_grnd = Float64.(ds["bmb_grnd"][:, :, t])

close(ds)

# greenland.jl's source: mass rate [kg/m2/s]. The routing divides by rho_w, so
# the volume rate is -bmb_grnd -- exactly what k24_compare.f90 passes as mdot.
# `mdotscale` lets the same script run with the raw restart units or a
# physically-scaled forcing; it must match k24_compare.f90's 3rd argument.
mdot_scale = parse(Float64, get_cfg("mdotscale", "1.0"))
mdot = -bmb_grnd .* mdot_scale .* 1000.0

# Bed-type indicator, mirroring initialize_kappa in k24.f90.
substrate = get_cfg("substrate", "hard")
kappa = zeros(Nx, Ny)
if substrate == "soft"
    kappa .= 1.0
elseif substrate == "mixed"
    kappa[b .< -1000.0] .= 1.0
elseif substrate != "hard"
    error("unknown substrate $substrate")
end

grad_conv = get_cfg("grad", "mean") == "local" ? LocalGradient() : MeanGradient()

wthick = get_cfg("wthick", "darcy")
water_thickness_algorithm =
    wthick == "laminar" ? LaminarThickness(gradient_convention = grad_conv) :
    wthick == "areal"   ? ArealConduitThickness() :
                          DarcyWeisbachThickness(gradient_convention = grad_conv)

psi = get_cfg("psi", "recursive")
psi_out_algorithm =
    psi == "iterative"   ? IterativePsiOut() :
    psi == "topological" ? TopologicalPsiOut(allow_cycles = true) :
                           RecursivePsiOut()

dm = get_cfg("drainage", "both")
drainage_mode =
    dm == "efficient"   ? EfficientOnly() :
    dm == "inefficient" ? InefficientOnly() :
                          BothDrainage()

ctill = parse(Float64, get_cfg("ctill", "0.2"))

sl = get_cfg("sliding", "none")
sliding_law =
    sl == "weertman"     ? WeertmanSlidingLaw(C = 1.0e5) :
    sl == "powerplastic" ? PowerPlasticSlidingLaw(c_till = ctill) :
    sl == "regcoulomb"   ? RegularizedCoulombSlidingLaw(c_till = ctill) :
                           NoSlidingLaw()

xlims, ylims = FastHydrology.compute_lims(x, y)

grid  = ArrayHydroGrid(Nx, Ny, xlims, ylims; T = Float64)
@info "config" Nx Ny dx=grid.dx substrate wthick grad_conv=get_cfg("grad","mean") psi sl

model = KazmierczakHydroModel(grid, kappa, abs_v_b, A_visc, mdot;
    longcoupwater             = parse(Float64, get_cfg("longcoup", "5.0")),
    water_thickness_algorithm = water_thickness_algorithm,
    psi_out_algorithm         = psi_out_algorithm,
    drainage_mode             = drainage_mode,
    sliding_law               = sliding_law,
    dissipation_melt          = parse(Bool, get_cfg("dissip", "true")))

state = HydroState(grid, mask, h, b)
sim   = SteadyStateSimulation(model, grid, state)

FastHydrology.run!(sim)

# update_Po!'s own definition (rho_i*g*h, no floor as of the Po-floor removal --
# see effective_pressure.jl), which is also what the Fortran uses for p_w and
# what greenland.jl recomputes independently; the two no longer diverge at thin
# ice, unlike when Po carried a 1e5 floor here that greenland.jl's rho_i*g*h
# recomputation didn't.
p_w = Array(model.Po) .- Array(state.N)

NCDataset(out_path, "c") do out
    defDim(out, "xc", Nx)
    defDim(out, "yc", Ny)
    defVar(out, "xc",   x ./ 1000.0, ("xc",))
    defVar(out, "yc",   y ./ 1000.0, ("yc",))
    defVar(out, "W",    Array(state.W), ("xc", "yc"))
    defVar(out, "N",    Array(state.N), ("xc", "yc"))
    defVar(out, "q",    Array(model.q), ("xc", "yc"))
    defVar(out, "p_w",  p_w,            ("xc", "yc"))
    defVar(out, "mask", mask,           ("xc", "yc"))
end

println("wrote $out_path")
