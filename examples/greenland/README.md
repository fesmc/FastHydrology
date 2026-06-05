# Greenland example

Loads a Yelmo 16 km Greenland restart from `input/GRL-16KM_yelmo_restart.nc`,
runs `FastHydrology` to a steady-state W_til / W field using the configured
till and transport methods, and writes a NetCDF output suitable for
plotting.

## Build

```sh
make greenland
```

## Run via `runme` (recommended)

The example is set up to be driven by [`runme`](https://github.com/fesmc/runme):

```sh
runme -r -e greenland -n examples/greenland/greenland.nml -o output/greenland
```

`runme` stages a clean rundir at `output/greenland/`, copies in the
namelist + executable, symlinks `input/` (so the restart file is
reachable), and runs the executable from that directory. Output ends up
at `output/greenland/greenland.nc`.

Permute by hand without storing extra namelist files: edit the in-tree
`greenland.nml` between runs, point each at a different `-o` rundir,
and you get one rundir per permutation. A small bash script that wraps
several `runme` invocations is a good way to automate this.

For example:

```sh
# Pure bucket (no transport)
sed -i.bak 's/method_transport.*/method_transport = 0/' examples/greenland/greenland.nml
runme -r -e greenland -n examples/greenland/greenland.nml -o output/greenland_bucket

# Bucket + K24
sed -i.bak 's/method_transport.*/method_transport = 1/' examples/greenland/greenland.nml
runme -r -e greenland -n examples/greenland/greenland.nml -o output/greenland_k24
```

## Run directly (no runme)

```sh
mkdir -p output/greenland_direct && cd output/greenland_direct
ln -sf ../../input input
../../bin/greenland.x ../../examples/greenland/greenland.nml
```

The example writes `greenland.nc` (and reads `input/GRL-16KM_yelmo_restart.nc`)
relative to the current working directory.

## Output

`greenland.nc` contains the final timestep stack of `W_til, dW_til_dt,
overflow, W, N, p_w, q_x, q_y` on the yelmo `(xc, yc, time)` grid. The
diagnostic line printed each output step shows:

```
       time  W_til_mean       W_mean    overflow       N_mean
```

`W_til` is in m, `W` is in m, `overflow` is in m/s, `N` is in Pa.

## Namelist switches

- `&fhyd { method_til }`: `0` = BUCKET (default), `1` = EXTERNAL (host
  owns `W_til`).
- `&fhyd { method_transport }`: `0` = NONE (W = 0, no transport; N from
  `bucket%N_closure`), `1` = K24.
- `&fhyd { W_til_max }`: scalar default for the per-cell till cap.
  Set to `0` to bypass the bucket entirely (all `mdot` flows through to
  transport).

In sequential coupling, the bucket runs first; once `W_til >= W_til_max(i,j)`
the saturation overflow becomes the source term for K24.

## Plot (Julia)

A `Project.toml` + `Manifest.toml` are included for a reproducible
environment. First-time setup:

```sh
julia --project=examples/greenland -e 'using Pkg; Pkg.instantiate()'
```

For a side-by-side comparison, run `greenland.x` twice (typically with
different `method_transport`), then plot:

```sh
julia --project=examples/greenland examples/greenland/plot_greenland.jl
```

The plotter defaults to `output/greenland_bucket/greenland.nc` and
`output/greenland_k24/greenland.nc` (the rundirs produced by the
example commands above). Override the filenames as positional args.
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
