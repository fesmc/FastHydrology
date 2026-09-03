# Export the Kazmierczak et al 2024 Thwaites-2km dataset to NetCDF so the
# Fortran K24 driver and the Julia reference read bit-identical inputs.
using FastHydrology, NCDatasets, Printf, Statistics

src = ARGS[1]; out = ARGS[2]
Nx, Ny, xlims, ylims, mask, h, b, abs_v_b, A_visc, G, q_T, mdot_mass, kappa =
    FastHydrology.load_Kazmierczak(src; bed_rheology = :hard)

# Cell-centred axes reconstructed from the limits compute_lims produced.
dx = (xlims[2]-xlims[1])/Nx
dy = (ylims[2]-ylims[1])/Ny
xc = [xlims[1] + dx/2 + (i-1)*dx for i in 1:Nx]
yc = [ylims[1] + dy/2 + (j-1)*dy for j in 1:Ny]

mdot_vol = mdot_mass ./ 1000.0    # water-equivalent volume rate [m/s], the Fortran convention

@printf("Nx=%d Ny=%d  dx=%.6f dy=%.6f  grounded=%d\n", Nx, Ny, dx, dy, count(mask .== 1))
@printf("h: %.1f..%.1f   b: %.1f..%.1f\n", minimum(h), maximum(h), minimum(b), maximum(b))
@printf("vb: %.3e..%.3e   A: %.3e..%.3e   mdot_vol: %.3e..%.3e\n",
        minimum(abs_v_b), maximum(abs_v_b), minimum(A_visc), maximum(A_visc),
        minimum(mdot_vol), maximum(mdot_vol))

NCDataset(out, "c") do o
    defDim(o, "xc", Nx); defDim(o, "yc", Ny)
    defVar(o, "xc", collect(xc), ("xc",)); defVar(o, "yc", collect(yc), ("yc",))
    for (n, v) in [("h", h), ("b", b), ("mask", Float64.(mask)),
                   ("vb", abs_v_b), ("A", A_visc), ("mdot", mdot_vol)]
        defVar(o, n, Array{Float64}(v), ("xc", "yc"))
    end
end
println("wrote $out")
