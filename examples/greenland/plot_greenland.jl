# Plot Greenland-example output: side-by-side comparison of BUCKET vs K24
# H_w, N, p_w, and water flux magnitude. Reads the two NetCDF files
# written by examples/greenland/greenland.x.
#
# Usage:
#     julia --project=. plot_greenland.jl                 # uses defaults
#     julia --project=. plot_greenland.jl B_FILE K_FILE OUT_PNG
#
# Dependencies (add to a local Project.toml or your global env):
#     NCDatasets, CairoMakie

using NCDatasets
using CairoMakie

const DEFAULT_BUCKET = joinpath(@__DIR__, "..", "..", "output", "greenland_bucket.nc")
const DEFAULT_K24    = joinpath(@__DIR__, "..", "..", "output", "greenland_k24.nc")
const DEFAULT_OUT    = joinpath(@__DIR__, "..", "..", "output", "greenland_compare.png")

bucket_file = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_BUCKET
k24_file    = length(ARGS) >= 2 ? ARGS[2] : DEFAULT_K24
out_png     = length(ARGS) >= 3 ? ARGS[3] : DEFAULT_OUT

# Read the final time slice of each variable from a NetCDF file
function load_last(ncfile, varname)
    NCDataset(ncfile, "r") do ds
        var = ds[varname]
        nt  = size(var, ndims(var))
        # restart format: (xc, yc, time) → take t=nt slice
        return Array(var[:, :, nt])
    end
end

function load_coords(ncfile)
    NCDataset(ncfile, "r") do ds
        return Array(ds["xc"][:]), Array(ds["yc"][:])
    end
end

println("loading bucket: $bucket_file")
println("loading K24   : $k24_file")

xc, yc = load_coords(bucket_file)
fields = ["H_w", "N", "p_w"]
field_labels = Dict("H_w" => "H_w [m]", "N" => "N [Pa]", "p_w" => "p_w [Pa]")

bucket_data = Dict(f => load_last(bucket_file, f) for f in fields)
k24_data    = Dict(f => load_last(k24_file,    f) for f in fields)

# Flux magnitude as a fourth panel
b_qx = load_last(bucket_file, "q_x")
b_qy = load_last(bucket_file, "q_y")
k_qx = load_last(k24_file,    "q_x")
k_qy = load_last(k24_file,    "q_y")
bucket_data["|q|"] = sqrt.(b_qx .^ 2 .+ b_qy .^ 2)
k24_data["|q|"]    = sqrt.(k_qx .^ 2 .+ k_qy .^ 2)
field_labels["|q|"] = "|q| [m^2 s^-1]"

plot_fields = ["H_w", "N", "p_w", "|q|"]

println("plotting → $out_png")

# 4 rows x 2 cols (bucket | k24), one row per field
fig = Figure(size = (1100, 1800))

for (i, field) in enumerate(plot_fields)
    bd = bucket_data[field]
    kd = k24_data[field]
    vmax = max(maximum(abs, bd), maximum(abs, kd))
    vmin = field == "H_w" ? 0.0 : -vmax
    if field in ("H_w", "|q|")
        crange = (0.0, vmax)
        cmap   = :viridis
    else
        crange = (vmin, vmax)
        cmap   = :balance
    end

    ax_b = Axis(fig[i, 1]; title = "BUCKET $(field_labels[field])",
                xlabel = "x [km]", ylabel = "y [km]", aspect = DataAspect())
    ax_k = Axis(fig[i, 2]; title = "K24    $(field_labels[field])",
                xlabel = "x [km]", ylabel = "y [km]", aspect = DataAspect())

    hm_b = heatmap!(ax_b, xc, yc, bd; colorrange = crange, colormap = cmap)
    hm_k = heatmap!(ax_k, xc, yc, kd; colorrange = crange, colormap = cmap)

    Colorbar(fig[i, 3], hm_k)
end

save(out_png, fig)
println("done.")
