# Greenland example

Loads a Yelmo 16 km Greenland restart from `input/GRL-16KM_yelmo_restart.nc`,
runs `FastHydrology` to a steady-state H_w field using either the BUCKET
or K24 model, and writes a NetCDF output suitable for plotting.

## Build

```sh
make greenland FESMUTILS_DIR=../yelmo/fesm-utils
```

## Run

```sh
mkdir -p output
./bin/greenland.x examples/greenland/greenland_bucket.nml
./bin/greenland.x examples/greenland/greenland_k24.nml
```

Produces `output/greenland_bucket.nc` and `output/greenland_k24.nc`.
Both files have the same variables (H_w, dHwdt, N, p_w, q_x, q_y) on
the (xc, yc, time) Yelmo grid.

## Plot (Julia)

A `Project.toml` + `Manifest.toml` are included for a reproducible
environment. First-time setup:

```sh
julia --project=examples/greenland -e 'using Pkg; Pkg.instantiate()'
```

Then on every run:

```sh
julia --project=examples/greenland examples/greenland/plot_greenland.jl
```

Writes `output/greenland_compare.png` showing H_w, N, p_w, and water
flux magnitude side-by-side for BUCKET and K24.

## Notes

- `bmb_grnd` from the restart is in ice-equivalent m/a with Yelmo's
  sign convention (negative when melting). The driver converts to
  water-equivalent and flips the sign:
  `bmb_w = -bmb_grnd * (rho_ice/rho_w)`.
- `ATT_bar` is used as `A_glen` (depth-averaged rate factor).
- The K24 internal mask is `f_grnd > 0 .and. f_ice > 0`.
- The domain border BC is `MASK_BC_ZERO` by default. Edit
  `&fhyd { mask_bc }` to switch to IMPOSED or MIRROR.
