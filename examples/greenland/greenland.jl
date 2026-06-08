using NCDatasets
using Statistics
using FastHydrology
using CairoMakie
using Oceananigans

function load_data(path::String; bed_rheology = :hard)

    ds = NCDataset(path)

    # grid size
    x = ds["xc"][:]
    y = ds["yc"][:]

    Nx = length(x)
    Ny = length(y)

    xlims = (minimum(x), maximum(x)) .* 1000 # convert to meters
    ylims = (minimum(y), maximum(y)) .* 1000 

    # --- extract fields (NOTE: order is (xc, yc, time)) ---
    t = 1

    mask = Int.(ds["f_ice"][:, :, t] .* ds["f_grnd"][:, :, t] .> 0.0)

    h = ds["H_ice"][:, :, t]
    b = ds["z_bed"][:, :, t]

    abs_v_b = ds["uxy_b"][:, :, t]

    A_visc = ds["ATT_bar"][:, :, t]

    ṁ = -ds["bmb_grnd"][:, :, t] * 1000.0 # 1000 is rho_w

    # --- bed rheology ---
    function initialize_κ!(Nx, Ny, b; bed_rheology = bed_rheology)

        T = eltype(b)

        κ = zeros(T, Nx, Ny)

        if bed_rheology == :hard

            # already zero

        elseif bed_rheology == :soft

            κ .= one(T)

        elseif bed_rheology == :mixed

            κ[b .< -1000] .= one(T)

        elseif bed_rheology == :mixed_smooth

            for i in 1:Nx, j in 1:Ny
                if b[i, j] <= -1500
                    κ[i, j] = one(T)
                elseif b[i, j] <= -500
                    κ[i, j] = (b[i, j] + 500) / (-1000)
                end
            end
        else
            error("Unknown bed_rheology = $bed_rheology")
        end

        return κ
    end

    κ = initialize_κ!(Nx, Ny, b; bed_rheology)

    close(ds)

    return Nx, Ny, xlims, ylims, mask, h, b, abs_v_b, A_visc, ṁ, κ, x, y
end

function plot_results(xc, yc, H_w, N, p_w, q; out_png = "k24_results.png", f_grnd = nothing, mask_ocean = true)

    # ---------------------------
    # fields
    # ---------------------------
    fields = ["H_w", "N", "p_w", "|q|"]

    data = Dict("H_w" => H_w, "N"   => N, "p_w" => p_w, "|q|" => q)

    field_labels = Dict("H_w" => "H_w [m]", "N"   => "N [Pa]", "p_w" => "p_w [Pa]", "|q|" => "|q| [m² s⁻¹]")

    # ---------------------------
    # ocean mask (same logic as script)
    # ---------------------------
    if mask_ocean
        @assert f_grnd !== nothing "f_grnd required for ocean masking"

        for f in keys(data)
            A = float.(data[f])
            A[f_grnd .<= 0.0] .= NaN
            data[f] = A
        end
    end

    # ---------------------------
    # same plotting helpers
    # ---------------------------
    function panel_range(A, field)
        vmask = .!isnan.(A)
        vmax  = any(vmask) ? maximum(abs, A[vmask]) : 1.0
        vmax  = vmax > 0 ? vmax : 1.0

        if field in ("H_w", "|q|")
            return (0.0, vmax)
        else
            return (-vmax, vmax)
        end
    end

    panel_cmap(field) = field in ("H_w", "|q|") ? :viridis : :balance

    # ---------------------------
    # figure
    # ---------------------------
    fig = Figure(size = (1800, 700))

    for (i, f) in enumerate(fields)

        A = data[f]
        rng = panel_range(A, f)
        cmap = panel_cmap(f)

        ax = Axis(fig[1, 2*i-1]; title = "K24 $(field_labels[f])", xlabel = "x [m]", ylabel = "y [m]", aspect = DataAspect())

        hm = heatmap!(ax, xc, yc, A; colorrange = rng, colormap = cmap, nan_color = :transparent)

        Colorbar(fig[1, 2*i], hm)
    end

    save(out_png, fig)
    println("saved → $out_png")

    return fig
end

function main()

    data_dir = "$(@__DIR__)/../../input/GRL-16KM_yelmo_restart.nc"
    longcoupwater = 5.0
    Nx, Ny, xlims, ylims, mask, h, b, abs_v_b, A_visc, ṁ, κ, x, y = load_data(data_dir);
    grid = OGRectHydroGrid(Nx, Ny, xlims, ylims; T = Float64)
    model = KazmierczakHydroModel(grid, κ, abs_v_b, A_visc, ṁ; longcoupwater = longcoupwater);
    state = HydroState(grid, mask, h, b);
    sim = SteadyStateSimulation(model, grid, state);
    FastHydrology.run!(sim)

    H_w = Oceananigans.interior(state.W, :, :, 1)
    N = Oceananigans.interior(state.N, :, :, 1)
    Po = model.rho_i .* model.g .* Oceananigans.interior(state.h, :, :, 1)
    p_w = Po .- N
    q = Oceananigans.interior(model.q, :, :, 1)

    f_grnd = mask

    plot_results(x, y, H_w, N, p_w, q; f_grnd = f_grnd, mask_ocean = true, out_png = "$(@__DIR__)/../../output/k24.jl_results.png")
    
end

main()
