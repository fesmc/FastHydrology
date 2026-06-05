# Greenland example

Loads a Yelmo 16 km Greenland restart from `input/GRL-16KM_yelmo_restart.nc`,
runs `FastHydrology` to a steady-state W_til / W field using the configured
till and transport methods, and writes a NetCDF output suitable for
plotting.

## Build

```sh
make greenland
```

## Run

```sh
mkdir -p output
./bin/greenland.x examples/greenland/greenland.nml
```

Produces `output/greenland.nc` with `W_til, dW_til_dt, overflow, W, N,
p_w, q_x, q_y` on the yelmo `(xc, yc, time)` grid.

The example takes a single namelist. Toggle the configuration by editing
`examples/greenland/greenland.nml`:

- `&fhyd { method_til }`: `0` = BUCKET (default), `1` = EXTERNAL (host
  owns `W_til`).
- `&fhyd { method_transport }`: `0` = NONE (W = 0, no transport; N from
  `bucket%N_closure`), `1` = K24 (writes W, q_x, q_y, N, p_w).
- `&greenland { out_file }`: rename when permuting so successive runs
  don't overwrite each other.

In sequential coupling, the bucket runs first; once `W_til >= W_til_max(i,j)`
the saturation overflow becomes the source term for K24. With
`W_til_max = 0` the bucket holds nothing and the full source `mdot`
passes through to transport.

## Plot (Julia)

A `Project.toml` + `Manifest.toml` are included for a reproducible
environment. First-time setup:

```sh
julia --project=examples/greenland -e 'using Pkg; Pkg.instantiate()'
```

For a side-by-side comparison, run `greenland.x` twice with different
`method_transport` and `out_file` values, then plot:

```sh
julia --project=examples/greenland examples/greenland/plot_greenland.jl
```

The plotter defaults to reading `output/greenland_bucket.nc` and
`output/greenland_k24.nc`; override the filenames as positional args.
Writes `output/greenland_compare.png`.

## Notes

- `bmb_grnd` from the restart is in ice-equivalent m/a with Yelmo's
  sign convention (negative when melting). The driver converts to
  a water-equivalent rate in m/s (SI):
  `mdot = -bmb_grnd * (rho_ice/rho_w) / SEC_PER_YEAR`.
- `ATT_bar` is used as `A_glen` (depth-averaged rate factor).
- The K24 internal mask is `f_grnd > 0 .and. f_ice > 0`.
- The domain border BC is `MASK_BC_ZERO` by default. Edit
  `&fhyd { mask_bc }` to switch to IMPOSED or MIRROR.
