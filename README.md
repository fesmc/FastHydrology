# FastHydrology

> **Compatibility (dev):** must be rebuilt against **fesm-utils dev at `3f415cc`
> (2026-06-26) or later** (the release that folded `coordinates` into `utils/src/coords/`
> and split `constants` out of `precision`). FastHydrology needs **no source changes** — it
> only uses the stable `nml`/`ncio` interfaces — but its `libfasthydro.a` must be relinked
> against the new fesm-utils.

A Fortran library of basal-hydrology models for ice-sheet simulations.
Two orthogonal switches select how till water storage and water transport
are handled, run sequentially each step.

| `method_til`        | what runs on `W_til` (till storage)                     |
|--------------------:|---------------------------------------------------------|
| `0` BUCKET (default)| local mass-balance bucket (van Pelt & Bueler 2015 style)|
| `1` EXTERNAL        | host owns `W_til`; library does not touch it            |

| `method_transport` | what runs on `W` (distributed sheet) and `N`            |
|-------------------:|---------------------------------------------------------|
| `0` NONE           | `W = 0`, `q_x = q_y = 0`; `N` from `bucket%N_closure`   |
| `1` K24            | Kazmierczak 2024 distributed model: `W`, `q_x`, `q_y`, `N`, `p_w`, `Q_b`, `Q_diss` |

The till step runs first and is driven by `mdot`. K24 does not take `mdot`:
its water source is built from the terms of the basal melt rate -- the
geothermal heat `G` and the heat conducted into the ice `q_T`, plus the
frictional heat `Q_b` and the dissipation heat `Q_diss` it computes itself
-- and the water reaching the bed from above, `i_eb` (drained englacial
water, surface input), which is routed but is not melt. The bucket overflow
no longer feeds K24. EXTERNAL is the coupling-friendly mode that lets a host
own `W_til` and use FastHydrology only for `N` and/or transport. Notation
follows van Pelt & Bueler 2015:

- `W_til` : till water storage thickness   [m]
- `W`     : distributed sheet thickness    [m]
- `mdot`  : bucket source rate             [m/s, water-equivalent]
- `G`, `q_T` : geothermal heat into the bed, conductive heat into the ice [W/m2]
- `i_eb`  : water reaching the bed from above [m/s, water-equivalent]

All internal units are SI. The public API takes `time` in years (matching
typical ice-sheet model conventions); namelist `bkt_till_rate` is in m/a
(converted at load time). Output diagnostics (`dW_til_dt`, `overflow`,
`q_x`, `q_y`) are in SI. The grid spacing pair `(dx, dy)` is carried
throughout — no isotropic-grid assumption.

## Build

The build is configured by [`configme`](https://github.com/fesmc/configme),
which fills `config/Makefile` with a machine + compiler fragment and
writes a resolved root `Makefile`:

```sh
configme -m macbook -c gfortran
make fasthydro-static
```

Then `include/libfasthydro.a` is ready to link against (target name
matches yelmo's `yelmo-static` / FastIsostasy's `isostasy-static`). See `config/Makefile`
for the build template and `config/common.mk` for the dependency wiring
(fesm-utils + FFTW + netCDF).

## Greenland example

End-to-end: build a Greenland-16km hydrology field from a yelmo restart,
run it for 1000 a, and inspect `W_til`, `W`, `overflow`, and `N`. The
example is driven by [`runme`](https://github.com/fesmc/runme):

```sh
make greenland
runme -r -e greenland -n examples/greenland/greenland.nml -o output/greenland
```

`runme` stages a clean rundir at `output/greenland/`, symlinks `input/`
for the restart file, and runs the executable from there. Output ends
up at `output/greenland/hydro.nc` with all eight fields (`W_til,
dW_til_dt, overflow, W, N, p_w, q_x, q_y`) on the yelmo `(xc, yc, time)`
grid. Permute configurations by editing the in-tree namelist's
`method_til` / `method_transport` switches between `runme` invocations,
pointing each at a different `-o` rundir; a small wrapper bash script is
the typical way to run a sweep. See
[`examples/greenland/README.md`](examples/greenland/README.md) for the
direct (non-runme) invocation, the plot script, and the namelist
switches.

Plot the two side-by-side (Julia; first-time `Pkg.instantiate()`):

```sh
julia --project=examples/greenland -e 'using Pkg; Pkg.instantiate()'
julia --project=examples/greenland examples/greenland/plot_greenland.jl
# → output/greenland_compare.png
```

See [`examples/greenland/README.md`](examples/greenland/README.md) for
the full setup, including the restart-file conventions and how to switch
the boundary condition.

## The K24 model and FastHydrology.jl

`src/k24.f90` is a port of the `kazmierczak2024` model in
[FastHydrology.jl](https://github.com/TakisAngelides/FastHydrology.jl)
(Kazmierczak et al. 2024, https://doi.org/10.5194/tc-18-5887-2024). Where the
two could differ, the Julia side is the source of truth: every `k24_*` namelist
default reproduces a `KazmierczakHydroModel` constructor keyword, and each
parameter's own comment in `input/yelmo_defaults.nml` names the Julia field it
maps to.

The water source is built exactly as in FastHydrology.jl, as a mass rate
[kg/m2/s]: `mdot_total = (G - q_T + Q_b + Q_diss)/L_w + i_eb`, seeded into the
routing as `mdot_total*dx*dy/rho_w`. Only the public `hydro_update` takes
`i_eb` as a water-equivalent volume rate [m/s], like the rest of the Fortran
API, and converts it.

The routing options of FastHydrology.jl are all available: `k24_routing_scheme`
(Warner by default, GDS-Warner the original K24/KORI scheme, Quinn, Tarboton,
modified and GDS Tarboton; Le Brocq et al. 2006), `k24_fill_algorithm`,
`k24_q_conversion`, `k24_dissipation_discretization`, the taped flow router
(`k24_flux_solver = 3`, the default; the recursive, iterative and topological
routers implement GDS-Warner only), the staggered (C-grid) frictional heat
(`k24_friction_discretization`) and the field-valued sliding laws (prescribed
`tau_b`, a per-cell Coulomb coefficient, Shakti's regularized Coulomb law).

Deliberate deviations:

- **Degenerate cases.** The reference lets IEEE arithmetic produce `Inf`/`NaN`
  at cells where `Q == 0`, `S_inf == 0` or `N_inf == 0` and then overwrites
  them; this library takes the same limits by an explicit branch instead, since
  it is built with `-Ofast` where `Inf`/`NaN` propagation is not dependable.
  Each such site is marked `DEVIATION:` in `src/k24.f90`.

Four clamps that KORI-ULB applies are **off by default**, matching the
reference constructor: `k24_W_min`, `k24_W_max`, `k24_min_pressure_fraction`
and `k24_q_max`. KORI-ULB's own values are noted next to each key in the schema
if you want them back — `k24_q_max` and `k24_min_pressure_fraction` in
particular guard real numerical edge cases, so turning them off is a genuine
tradeoff rather than a free simplification.

The port is verified against the Julia implementation on identical inputs by
`tests/k24_synth.f90` (every configuration in `tests/k24_crossvalidate.sh`,
agreeing to ~1e-15 in double precision) and `tests/k24_greenland.f90` (the real 16 km restart, ~1e-7 through
the `real(sp)` public API). See `tests/README.md`.

## SHMIP driver

`tests/shmip.f90` runs the SHMIP A–D steady-state benchmarks for quick
verification:

```sh
make shmip
./bin/shmip.x par/shmip.nml      # case selected by &shmip { case }
```
