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
- `mdot`  : source rate from ice base      [m/a today; m/s after Commit 2]

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
run it for 1000 a, and compare BUCKET-only vs BUCKET+K24 side-by-side.

```sh
make greenland
mkdir -p output

# BUCKET only (method_til=BUCKET, method_transport=NONE)
./bin/greenland.x examples/greenland/greenland_bucket.nml

# BUCKET + K24 (method_til=BUCKET, method_transport=K24)
./bin/greenland.x examples/greenland/greenland_k24.nml
```

Each run writes a NetCDF file (`output/greenland_bucket.nc`,
`output/greenland_k24.nc`) with `W_til, dW_til_dt, overflow, W, N, p_w,
q_x, q_y` on the yelmo `(xc, yc, time)` grid, and prints a one-line
summary per output step.

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
