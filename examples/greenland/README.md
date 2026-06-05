# Greenland example

Loads a Yelmo 16 km Greenland restart from `input/GRL-16KM_yelmo_restart.nc`,
runs `FastHydrology` to a steady-state W_til / W field using the configured
till and transport methods, and writes a NetCDF output suitable for plotting.

## Build

```sh
make greenland
```

## Run

```sh
mkdir -p output
./bin/greenland.x examples/greenland/greenland_bucket.nml
./bin/greenland.x examples/greenland/greenland_k24.nml
```

Produces `output/greenland_bucket.nc` and `output/greenland_k24.nc`.
Both files have the same variables (W_til, dW_til_dt, overflow, W, N,
p_w, q_x, q_y) on the (xc, yc, time) Yelmo grid. Switch between modes by
editing the namelist:

- `method_til`: `0` = BUCKET (default), `1` = EXTERNAL (host owns W_til).
- `method_transport`: `0` = NONE (W=0, no transport; N from
  `bucket%N_closure`), `1` = K24 (writes W, q_x, q_y, N, p_w).

In sequential coupling, the bucket runs first; once `W_til >= W_til_max(i,j)`
the saturation overflow becomes the source term for the K24 transport
model. With `W_til_max = 0` the bucket holds nothing and the full source
`mdot` passes through to the transport model.

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

Writes `output/greenland_compare.png` showing W_til, N, p_w, and water
flux magnitude side-by-side for the BUCKET-only and BUCKET+K24 runs.

## Notes

- `bmb_grnd` from the restart is in ice-equivalent m/a with Yelmo's
  sign convention (negative when melting). The driver converts to
  water-equivalent and flips the sign:
  `mdot = -bmb_grnd * (rho_ice/rho_w)`.
- `ATT_bar` is used as `A_glen` (depth-averaged rate factor).
- The K24 internal mask is `f_grnd > 0 .and. f_ice > 0`.
- The domain border BC is `MASK_BC_ZERO` by default. Edit
  `&fhyd { mask_bc }` to switch to IMPOSED or MIRROR.
