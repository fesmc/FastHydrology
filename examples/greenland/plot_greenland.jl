# Plot Greenland-example output: side-by-side comparison of two NetCDF
# files produced by examples/greenland/greenland.x (W_til, N, p_w, and
# water-flux magnitude). The two files come from two runme invocations
# of the example with different &fhyd { method_transport } values --
# e.g. method_transport=0 saved to output/greenland_bucket/greenland.nc
# and method_transport=1 saved to output/greenland_k24/greenland.nc.
#
# Usage:
#     julia --project=examples/greenland plot_greenland.jl                # mask ocean (default)
#     julia --project=examples/greenland plot_greenland.jl --no-mask      # show ocean
#     julia --project=examples/greenland plot_greenland.jl [opts] B_FILE K_FILE OUT_PNG RESTART
#
# By default cells with f_grnd <= 0 (ocean and partially-floating) are
# masked out so the K24 ice-sheet pattern is visible at its native sheet
# thickness instead of being washed out by the saturated ocean cells.
#
# Dependencies (see Project.toml): NCDatasets, CairoMakie.

using NCDatasets
using CairoMakie

const DEFAULT_BUCKET  = joinpath(@__DIR__, "..", "..", "output", "greenland_bucket", "greenland.nc")
const DEFAULT_K24     = joinpath(@__DIR__, "..", "..", "output", "greenland_k24",    "greenland.nc")
const DEFAULT_OUT     = joinpath(@__DIR__, "..", "..", "output", "greenland_compare.png")
const DEFAULT_RESTART = joinpath(@__DIR__, "..", "..", "input",  "GRL-16KM_yelmo_restart.nc")

# --- Parse args ---
mask_ocean = true
positional = String[]
for a in ARGS
    if a == "--no-mask" || a == "--no-ocean-mask"
        mask_ocean = false
    elseif a == "--mask" || a == "--ocean-mask"
        mask_ocean = true
    else
        push!(positional, a)
    end
end
bucket_file  = length(positional) >= 1 ? positional[1] : DEFAULT_BUCKET
k24_file     = length(positional) >= 2 ? positional[2] : DEFAULT_K24
out_png      = length(positional) >= 3 ? positional[3] : DEFAULT_OUT
restart_file = length(positional) >= 4 ? positional[4] : DEFAULT_RESTART

println("loading bucket: $bucket_file")
println("loading K24   : $k24_file")
mask_ocean && println("loading mask  : $restart_file")

# Read the final time slice of a (xc, yc, time) variable
function load_last(ncfile, varname)
    NCDataset(ncfile, "r") do ds
        var = ds[varname]
        nt  = size(var, ndims(var))
        return Array(var[:, :, nt])
    end
end

function load_coords(ncfile)
    NCDataset(ncfile, "r") do ds
        return Array(ds["xc"][:]), Array(ds["yc"][:])
    end
end

# Apply ocean mask: replace cells with f_grnd <= 0 by NaN so they render
# as the heatmap's nan_color (transparent / grey).
function apply_ocean_mask(data, f_grnd)
    out = float.(data)
    out[f_grnd .<= 0.0] .= NaN
    return out
end

xc, yc = load_coords(bucket_file)
fields = ["W_til", "N", "p_w"]
field_labels = Dict("W_til" => "W_til [m]", "N" => "N [Pa]", "p_w" => "p_w [Pa]")

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

if mask_ocean
    f_grnd = load_last(restart_file, "f_grnd")
    for f in keys(bucket_data)
        bucket_data[f] = apply_ocean_mask(bucket_data[f], f_grnd)
        k24_data[f]    = apply_ocean_mask(k24_data[f],    f_grnd)
    end
end

plot_fields = ["W_til", "N", "p_w", "|q|"]

# BUCKET only writes W_til and (via the closure post-step) N. p_w and |q|
# carry no useful information in BUCKET mode -- skip those panels.
const BUCKET_FIELDS = Set(["W_til", "N"])

println("plotting → $out_png")

# Colorbar range over the *visible* (non-NaN) data so ocean cells don't
# determine the scale. (0, vmax) for nonnegative fields, otherwise (-v, v).
function panel_range(data, field)
    vmask = .!isnan.(data)
    vmax  = any(vmask) ? maximum(abs, data[vmask]) : 0.0
    vmax  = vmax > 0 ? vmax : 1.0
    if field in ("W_til", "|q|")
        return (0.0, vmax)
    else
        return (-vmax, vmax)
    end
end

panel_cmap(field) = field in ("W_til", "|q|") ? :viridis : :balance

# 4 rows x 4 cols (bucket | cbar | k24 | cbar), one row per field
fig = Figure(size = (1400, 1800))

for (i, field) in enumerate(plot_fields)
    kd     = k24_data[field]
    k_rng  = panel_range(kd, field)
    cmap   = panel_cmap(field)

    if field in BUCKET_FIELDS
        bd    = bucket_data[field]
        b_rng = panel_range(bd, field)
        ax_b  = Axis(fig[i, 1]; title = "BUCKET $(field_labels[field])",
                     xlabel = "x [km]", ylabel = "y [km]", aspect = DataAspect())
        hm_b  = heatmap!(ax_b, xc, yc, bd; colorrange = b_rng, colormap = cmap, nan_color = :transparent)
        Colorbar(fig[i, 2], hm_b)
    else
        Label(fig[i, 1], "BUCKET $(field_labels[field])\n(not produced)";
              tellwidth = false, tellheight = false, fontsize = 14, color = :gray)
    end

    ax_k = Axis(fig[i, 3]; title = "K24    $(field_labels[field])",
                xlabel = "x [km]", ylabel = "y [km]", aspect = DataAspect())
    hm_k = heatmap!(ax_k, xc, yc, kd; colorrange = k_rng, colormap = cmap, nan_color = :transparent)
    Colorbar(fig[i, 4], hm_k)
end

save(out_png, fig)
println("done.")
