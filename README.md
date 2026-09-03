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
| `1` K24            | Kazmierczak 2024 distributed model: `W`, `q_x`, `q_y`, `N`, `p_w` |

The till step runs first. When BUCKET is on, any source `mdot` that does
not fit under `W_til_max` spills over to feed the transport step as its
source. With `W_til_max = 0` the bucket holds nothing and all source
flows through to transport. EXTERNAL is the coupling-friendly mode that
lets a host own `W_til` and use FastHydrology only for `N` and/or
transport. Notation follows van Pelt & Bueler 2015:

- `W_til` : till water storage thickness   [m]
- `W`     : distributed sheet thickness    [m]
- `mdot`  : source rate from ice base      [m/s, water-equivalent]

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

Two conventions differ deliberately:

- **Units.** FastHydrology.jl carries the melt rate as a mass rate [kg/m2/s]
  and divides by `rho_w` when seeding the flow routing. This library carries it
  as the water-equivalent volume rate `mdot` [m/s] used everywhere else in the
  Fortran API, so every melt-like source term picks up an extra `1/rho_w`
  (`tau_b*v_b/(L_w*rho_w)` rather than `tau_b*v_b/L_w`, and likewise for the
  dissipation term). Everything else is SI and identical.
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
`tests/k24_synth.f90` (18 configurations, agreeing to ~1e-15 in double
precision) and `tests/k24_greenland.f90` (the real 16 km restart, ~1e-7 through
the `real(sp)` public API). See `tests/README.md`.

## SHMIP driver

`tests/shmip.f90` runs the SHMIP A–D steady-state benchmarks for quick
verification:

```sh
make shmip
./bin/shmip.x par/shmip.nml      # case selected by &shmip { case }
```
