# Julia twin of k24_synth.f90 -- reads the shared synthetic fields and runs
# FastHydrology.jl's kazmierczak2024 model on them.
#
# usage: julia k24_synth.jl <out.nc> [key=value ...]
#   aglen=<float>  substrate=hard|soft|mixed  drainage=both|efficient|inefficient
#   wthick=darcy|laminar|areal  grad=mean|local  kamb86=<float>
#   routing=warner|gdswarner|quinn|quinnorig|tarboton|modtarboton|gdstarboton
#   fill=auto|jacobi|lowest|flood  qconv=auto|outflow|face  dissdisc=auto|cell|face
#   friction=cell|staggered|staggeredquad
#   psi=taped|recursive|iterative|topological  dissip=true|false
#   sliding=none|weertman|powerplastic|regcoulomb|field|regcoulombfield|shakti  ctill=<float>
#   qtfrac=<float>  iebfrac=<float>  report=0|1  maxpsi=<int>  filliters=<int>  input=<path>
#
# Every option must be matched by the corresponding k24_* key in the namelist
# given to k24_synth.x (qtfrac/iebfrac are passed to it as arguments) -- an
# option exposed on one side but not the other silently compares two
# different configurations.
#
# The water source is built from terms from the file's volume rate `mdot`
# [m/s]: G = mdot*rho_w*L_w, q_T = qtfrac*G, i_eb = iebfrac*mdot*rho_w. The
# fields some options need are derived from vb exactly as k24_synth.f90 does.

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
xc = Array{Float64}(ds["xc"][:]); yc = Array{Float64}(ds["yc"][:])
dx = xc[2] - xc[1];               dy = yc[2] - yc[1]
h    = Array{Float64}(ds["h"][:, :]);    b    = Array{Float64}(ds["b"][:, :])
mask = Array{Float64}(ds["mask"][:, :]); vb   = Array{Float64}(ds["vb"][:, :])
A    = Array{Float64}(ds["A"][:, :]);    mdot = Array{Float64}(ds["mdot"][:, :])
close(ds)
const Nx, Ny = size(h)

# Terms of the water source and the derived fields, same operations and order
# as k24_synth.f90.
qtfrac  = parse(Float64, get_cfg("qtfrac", "0.0"))
iebfrac = parse(Float64, get_cfg("iebfrac", "0.0"))
ctill   = parse(Float64, get_cfg("ctill", "0.02"))
G      = mdot .* 1000.0 .* 3.34e5
q_T    = qtfrac .* G
i_eb   = iebfrac .* mdot .* 1000.0
ux_b   = vb
uy_b   = 0.5 .* vb
taub   = 2.0e4 .+ vb .* 1.0e10
c_till = ctill .* (1.0 .+ vb .* 1.0e5)

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

psi = get_cfg("psi", "taped")
psi_out_algorithm =
    psi == "recursive"   ? RecursivePsiOut() :
    psi == "iterative"   ? IterativePsiOut() :
    psi == "topological" ? TopologicalPsiOut(allow_cycles = true) :
                           TapedPsiOut()

rs = get_cfg("routing", "warner")
routing_scheme =
    rs == "gdswarner"   ? GDSWarner() :
    rs == "quinn"       ? Quinn() :
    rs == "quinnorig"   ? Quinn(original = true) :
    rs == "tarboton"    ? Tarboton() :
    rs == "modtarboton" ? ModifiedTarboton() :
    rs == "gdstarboton" ? GDSTarboton() :
                          Warner()

fa = get_cfg("fill", "auto")
fill_algorithm =
    fa == "jacobi" ? JacobiFill() :
    fa == "lowest" ? LowestNeighbourFill() :
    fa == "flood"  ? PriorityFloodFill() :
                     nothing

qc = get_cfg("qconv", "auto")
q_conversion = qc == "outflow" ? QFromOutflow() : qc == "face" ? QFromFaceAverage() : nothing

dd = get_cfg("dissdisc", "auto")
dissipation_discretization = dd == "cell" ? CellCentredDissipation() : dd == "face" ? FaceDissipation() : nothing

fr = get_cfg("friction", "cell")
friction_discretization =
    fr == "staggered"     ? StaggeredFriction(ux_b, uy_b) :
    fr == "staggeredquad" ? StaggeredFriction(ux_b, uy_b; quadrature = true) :
                            CellCentredFriction()

dm = get_cfg("drainage", "both")
drainage_mode =
    dm == "efficient"   ? EfficientOnly() :
    dm == "inefficient" ? InefficientOnly() :
                          BothDrainage()

# Half-cell-padded limits so grid.dx/dy come out exactly dx/dy.
xlims = (xc[1] - dx/2, xc[Nx] + dx/2)
ylims = (yc[1] - dy/2, yc[Ny] + dy/2)
grid  = ArrayHydroGrid(Nx, Ny, xlims, ylims; T = Float64)
@assert isapprox(grid.dx, dx; rtol = 1e-12) && isapprox(grid.dy, dy; rtol = 1e-12)

sl = get_cfg("sliding", "none")
sliding_law =
    sl == "weertman"        ? WeertmanSlidingLaw(C = 1.0e5) :
    sl == "powerplastic"    ? PowerPlasticSlidingLaw(c_till = ctill) :
    sl == "regcoulomb"      ? RegularizedCoulombSlidingLaw(c_till = ctill) :
    sl == "field"           ? PrescribedFieldSlidingLaw(grid, taub) :
    sl == "regcoulombfield" ? RegularizedCoulombFieldSlidingLaw(grid, c_till) :
    sl == "shakti"          ? ShaktiRegularizedCoulombSlidingLaw(grid, A, ctill) :
                              NoFrictionSlidingLaw()

model = KazmierczakHydroModel(grid, kappa, vb, A, G, q_T;
    i_eb                       = i_eb,
    coupling_length_kamb86     = parse(Float64, get_cfg("kamb86", "10.0")),
    water_thickness_algorithm  = water_thickness_algorithm,
    psi_out_algorithm          = psi_out_algorithm,
    routing_scheme             = routing_scheme,
    fill_algorithm             = fill_algorithm,
    q_conversion               = q_conversion,
    dissipation_discretization = dissipation_discretization,
    friction_discretization    = friction_discretization,
    drainage_mode              = drainage_mode,
    sliding_law                = sliding_law,
    max_psi_out_calls          = parse(Int, get_cfg("maxpsi", "100000")),
    fill_iters                 = parse(Int, get_cfg("filliters", "10")),
    dissipation_melt           = parse(Bool, get_cfg("dissip", "true")))

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
    defVar(out, "W",      Array(state.W),      ("xc", "yc"))
    defVar(out, "N",      Array(state.N),      ("xc", "yc"))
    defVar(out, "q",      Array(model.q),      ("xc", "yc"))
    defVar(out, "p_w",    p_w,                 ("xc", "yc"))
    defVar(out, "Q_b",    Array(model.Q_b),    ("xc", "yc"))
    defVar(out, "Q_diss", Array(model.Q_diss), ("xc", "yc"))
    defVar(out, "mask",   mask,                ("xc", "yc"))
    defVar(out, "gsx",    Array(model.minus_grad_phi0_sx), ("xc", "yc"))
    defVar(out, "gsy",    Array(model.minus_grad_phi0_sy), ("xc", "yc"))
    defVar(out, "absgs",  Array(model.abs_grad_phi0_s),    ("xc", "yc"))
    defVar(out, "absg",   Array(model.abs_grad_phi0),      ("xc", "yc"))
    defVar(out, "phi0",   Array(model.phi0),               ("xc", "yc"))
end
println("wrote $out_path")
