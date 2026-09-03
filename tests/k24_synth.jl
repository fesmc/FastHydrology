# Julia twin of k24_synth.f90 -- builds bit-identical synthetic fields and runs
# FastHydrology.jl's kazmierczak2024 model on them.
#
# usage: julia k24_synth.jl <out.nc> [key=value ...]
#   aglen=<float>  substrate=hard|soft|mixed  drainage=both|efficient|inefficient
#   wthick=darcy|laminar|areal  grad=mean|local  longcoup=<float>
#   psi=recursive|iterative|topological  dissip=true|false
#   sliding=none|weertman|powerplastic|regcoulomb  ctill=<float>  report=0|1
#   maxpsi=<int>   filliters=<int>   input=<path>
#
# Every option must be matched by the corresponding k24_* key in the namelist
# given to k24_synth.x -- an option exposed on one side but not the other
# silently compares two different configurations.

using NCDatasets
using FastHydrology
using Printf

out_path = ARGS[1]
cfg = Dict{String,String}()
for a in ARGS[2:end]
    k, v = split(a, "="; limit = 2)
    cfg[k] = v
end
get_cfg(k, d) = get(cfg, k, d)

const IN = get_cfg("input", "tests/k24_synth_input.nc")

ds = NCDataset(IN)
# Grid size and spacing come from the file, so the same script runs the small
# synthetic case and a full ice-sheet dataset (e.g. Thwaites 2km).
xc = Array{Float64}(ds["xc"][:]); yc = Array{Float64}(ds["yc"][:])
dx = xc[2] - xc[1];               dy = yc[2] - yc[1]
h    = Array{Float64}(ds["h"][:, :]);    b    = Array{Float64}(ds["b"][:, :])
mask = Array{Float64}(ds["mask"][:, :]); vb   = Array{Float64}(ds["vb"][:, :])
A    = Array{Float64}(ds["A"][:, :]);    mdot = Array{Float64}(ds["mdot"][:, :])
close(ds)
const Nx, Ny = size(h)

# initialize_kappa's three cases
substrate = get_cfg("substrate", "hard")
kappa = zeros(Nx, Ny)
substrate == "soft"  && (kappa .= 1.0)
substrate == "mixed" && (kappa[b .< -1000.0] .= 1.0)

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

ctill = parse(Float64, get_cfg("ctill", "0.02"))
sl = get_cfg("sliding", "none")
sliding_law =
    sl == "weertman"     ? WeertmanSlidingLaw(C = 1.0e5) :
    sl == "powerplastic" ? PowerPlasticSlidingLaw(c_till = ctill) :
    sl == "regcoulomb"   ? RegularizedCoulombSlidingLaw(c_till = ctill) :
                           PrescribedFrictionSlidingLaw()

# Half-cell-padded limits so grid.dx/dy come out exactly dx/dy.
xlims = (xc[1] - dx/2, xc[Nx] + dx/2)
ylims = (yc[1] - dy/2, yc[Ny] + dy/2)
grid  = ArrayHydroGrid(Nx, Ny, xlims, ylims; T = Float64)
@assert isapprox(grid.dx, dx; rtol = 1e-12) && isapprox(grid.dy, dy; rtol = 1e-12)

model = KazmierczakHydroModel(grid, kappa, vb, A, mdot .* 1000.0;   # mass rate = volume rate * rho_w
    longcoupwater             = parse(Float64, get_cfg("longcoup", "5.0")),
    water_thickness_algorithm = water_thickness_algorithm,
    psi_out_algorithm         = psi_out_algorithm,
    drainage_mode             = drainage_mode,
    sliding_law               = sliding_law,
    max_psi_out_calls         = parse(Int, get_cfg("maxpsi", "50000")),
    fill_iters                = parse(Int, get_cfg("filliters", "10")),
    dissipation_melt          = parse(Bool, get_cfg("dissip", "true")))

state = HydroState(grid, mask, h, b)
FastHydrology.run!(SteadyStateSimulation(model, grid, state))

if get_cfg("report", "0") == "1"
    mk   = Bool.(mask .== 1)
    Ninf = Array(model.N_inf)[mk]; Po = Array(model.Po)[mk]; S = Array(model.S_inf)[mk]
    @printf("  grounded=%d  N_inf==Po: %d (%.1f%%)  S_inf==0: %d  mean(N_inf/Po)=%.5f\n",
            length(Po), count(Ninf .>= Po .- 1e-9), 100*count(Ninf .>= Po .- 1e-9)/length(Po),
            count(S .== 0.0), sum(Ninf ./ Po)/length(Po))
end

p_w = Array(model.Po) .- Array(state.N)

NCDataset(out_path, "c") do out
    defDim(out, "xc", Nx); defDim(out, "yc", Ny)
    defVar(out, "xc", xc, ("xc",))
    defVar(out, "yc", yc, ("yc",))
    defVar(out, "W",    Array(state.W), ("xc", "yc"))
    defVar(out, "N",    Array(state.N), ("xc", "yc"))
    defVar(out, "q",    Array(model.q), ("xc", "yc"))
    defVar(out, "p_w",  p_w,            ("xc", "yc"))
    defVar(out, "mask", mask,           ("xc", "yc"))
    defVar(out, "gsx",   Array(model.minus_grad_phi0_sx), ("xc", "yc"))
    defVar(out, "gsy",   Array(model.minus_grad_phi0_sy), ("xc", "yc"))
    defVar(out, "absgs", Array(model.abs_grad_phi0_s),    ("xc", "yc"))
    defVar(out, "absg",  Array(model.abs_grad_phi0),      ("xc", "yc"))
    defVar(out, "phi0",  Array(model.phi0),               ("xc", "yc"))
end
println("wrote $out_path")
