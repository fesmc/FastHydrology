using NCDatasets, Statistics, Printf

f = NCDataset(ARGS[1])   # fortran
j = NCDataset(ARGS[2])   # julia

m = Bool.(j["mask"][:, :] .== 1.0)
n = count(m)
@printf("grounded cells: %d of %d\n\n", n, length(m))
@printf("%-8s %14s %14s %14s %14s %12s\n", "field", "mean(J)", "mean(F)", "max|J|", "max|F-J|", "rel_max")

for v in ["W", "N", "q", "p_w"]
    F = Float64.(f[v][:, :])[m]
    J = Float64.(j[v][:, :])[m]
    bad = .!(isfinite.(F) .& isfinite.(J))
    if any(bad)
        @printf("%-8s  NON-FINITE: F=%d J=%d\n", v, count(.!isfinite.(F)), count(.!isfinite.(J)))
    end
    ok = .!bad
    Fo, Jo = F[ok], J[ok]
    scale = maximum(abs, Jo)
    d = maximum(abs, Fo .- Jo)
    @printf("%-8s %14.6e %14.6e %14.6e %14.6e %12.3e\n",
            v, mean(Jo), mean(Fo), scale, d, scale > 0 ? d/scale : 0.0)
end

println()
println("correlation (grounded cells):")
for v in ["W", "N", "q", "p_w"]
    F = Float64.(f[v][:, :])[m]; J = Float64.(j[v][:, :])[m]
    ok = isfinite.(F) .& isfinite.(J)
    @printf("  %-6s %.10f\n", v, cor(F[ok], J[ok]))
end
close(f); close(j)
