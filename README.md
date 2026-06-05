# FastHydrology

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

## SHMIP driver

`tests/shmip.f90` runs the SHMIP A–D steady-state benchmarks for quick
verification:

```sh
make shmip
./bin/shmip.x par/shmip.nml      # case selected by &shmip { case }
```
