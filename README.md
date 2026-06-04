# FastHydrology

A Fortran library of basal-hydrology models for ice-sheet simulations.
Selectable at runtime via a single `method` switch:

| `method` | name     | what FastHydrology writes                         |
|---------:|----------|---------------------------------------------------|
|       -1 | EXTERNAL | nothing on `H_w`; closure writes `N`, `p_w`       |
|        0 | NONE     | nothing                                           |
|        1 | BUCKET   | `H_w` (local mass-balance); optional `N`, `p_w`   |
|        2 | K24      | `H_w` (steady-state), `q_x`, `q_y`, `N`, `p_w`    |

EXTERNAL is the coupling-friendly mode: an ice-sheet model owns `H_w`,
and FastHydrology only derives effective pressure `N` (and `p_w`) from
the configured closure on the host's water field. K24 produces an
equilibrium `H_w`; hosts can evolve their own field toward it with the
public elemental helper `relax_H_w(H_w, H_w_eq, tau, dt)`.

## Build

The build is configured by [`configme`](https://github.com/fesmc/configme),
which fills `config/Makefile` with a machine + compiler fragment and
writes a resolved root `Makefile`:

```sh
configme -m macbook -c gfortran
make lib
```

Then `lib/libfasthydro.a` is ready to link against. See `config/Makefile`
for the build template and `config/common.mk` for the dependency wiring
(fesm-utils + FFTW + netCDF).

## Greenland example

End-to-end: build a Greenland-16km hydrology field from a yelmo restart,
run it for 1000 a, and compare BUCKET vs K24 side-by-side.

```sh
make greenland
mkdir -p output

# Bucket model (local mass balance + N closure)
./bin/greenland.x examples/greenland/greenland_bucket.nml

# K24 model (steady-state distributed)
./bin/greenland.x examples/greenland/greenland_k24.nml
```

Each run writes a NetCDF file (`output/greenland_bucket.nc`,
`output/greenland_k24.nc`) with `H_w, dHwdt, N, p_w, q_x, q_y` on the
yelmo `(xc, yc, time)` grid, and prints a one-line summary per output
step:

```
       time      H_w_max     H_w_mean       N_max      N_mean
    100.0000   2.0000E+00   3.3770E-01   2.9729E+07   1.5435E+07
   ...
   1000.0000   2.0000E+00   5.0349E-01   2.9729E+07   1.5435E+07
```

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
