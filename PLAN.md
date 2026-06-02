# FastHydrology — design & implementation plan

Stand-alone basal-hydrology library wrapping multiple models behind a uniform
API. Two models in v1: a local bucket (yelmo-equivalent) and the K24
steady-state distributed model. Built as a static library against fesm-utils
(FFTW, LIS). Drivable stand-alone via SHMIP test programs; designed to wire
into yelmo without changing the API.

## Goals

- One library, multiple basal-hydrology methods, uniform call signature.
- Each method declares its own H_w ownership semantics; host can opt in or
  out per method.
- Diagnostic outputs (`q_x, q_y, N, p_w, dHwdt`) populated honestly per
  method; nothing pretends to be something it isn't.
- Independent of yelmo at link time, but plug-compatible with yelmo's
  current bucket + effective-pressure-closure usage.

## Method dispatch

`par%method` ∈ {`HYDRO_METHOD_NONE`, `HYDRO_METHOD_BUCKET`, `HYDRO_METHOD_K24`}.
Drop the unused `HYDRO_METHOD_RESERVED = 1`.

| Method | Writes H_w | Writes q_x,q_y | Writes N        | Writes p_w        |
|--------|------------|----------------|-----------------|-------------------|
| NONE   | no         | no             | no              | no                |
| BUCKET | yes        | no             | only if closure | only if closure   |
| K24    | optional¹  | yes            | yes             | yes               |

¹ Gated on `par%k24%update_H_w`; instant (`tau_H = 0`) or relaxation
(`tau_H > 0`) toward `H_w_eq`.

## API surface

```fortran
subroutine hydro_init(hyd, filename, nx, ny)
subroutine hydro_init_state(hyd, z_bed, f_ice, f_grnd, time)
subroutine hydro_update(hyd, H_ice, z_bed, f_ice, f_grnd, mask, &
                       bmb_w, uxy_b, A_glen, time)
```

Public standalone closure subs (callable regardless of `par%method`):

```fortran
hydro_calc_N_overburden(N, ...)
hydro_calc_N_marine    (N, ...)   ! Leguy 2014, geometry only
hydro_calc_N_till      (N, ...)   ! van Pelt & Bueler 2015
hydro_calc_N_two_value (N, ...)
```

## Parameter struct

```
par%method
par%init_method                 ! HYDRO_INIT_ZERO | HYDRO_INIT_EXTERNAL
par%dx, par%dy
par%H_w_max                     ! single value; used by bucket cap + closure saturation

par%bucket%till_rate            ! default 1e-3 m/a
par%bucket%N_closure            ! NONE | OVERBURDEN | MARINE | TILL | TWO_VALUE
par%bucket%n_marine%(...)
par%bucket%n_till%(...)
par%bucket%n_two_value%(...)

par%k24%update_H_w              ! .true. or .false.
par%k24%tau_H                   ! 0 = instant, > 0 = relax
par%k24%(...existing K24 params...)
```

## State struct (flat)

```
hyd%now%time, dt, initialized   ! initialized: bool replaces 1.0e10 sentinel
hyd%now%H_w
hyd%now%dHwdt                   ! (H_w_new - H_w_old)/dt, natural sign
hyd%now%p_w
hyd%now%q_x, q_y
hyd%now%N
hyd%now%kappa                   ! filled only when method = K24
```

All fields allocated unconditionally; kappa is only *filled* in
`hydro_init_state` when `method = K24`.

## Precision

- `wp = sp` (status quo; library-wide constant).
- K24 internals stay `dp`; `hydro_update` does the sp→dp→sp round-trip at the K24 boundary.
- No preprocessor flag for precision in v1.

## Init flow

`hydro_init`:
- loads namelist
- allocates state
- `initialized = .false.`

`hydro_init_state(z_bed, f_ice, f_grnd, time)`:
- if `method = K24`: call `initialize_kappa`
- per `init_method`:
  - `ZERO`: H_w = 0 on grounded cells
  - `EXTERNAL`: leave grounded H_w untouched (host pre-filled it)
- always apply override on floating + adjacent-to-floating: H_w = H_w_max
- zero p_w, q_x, q_y, N, dHwdt
- `time = time`, `dt = 0`, `initialized = .true.`

## Update flow

`hydro_update(...)`:

1. Assert `hyd%now%initialized`.
2. Stash `H_w_old = hyd%now%H_w`.
3. `dt = max(time - hyd%now%time, 0); hyd%now%time = time; hyd%now%dt = dt`.
4. If `dt == 0`: skip to (8) (still recompute dHwdt = 0).
5. Dispatch on `par%method`:
   - `NONE`: no-op.
   - `BUCKET`: yelmo's 4-case logic on `(f_grnd, f_ice, neighbors)`; grounded
     ice-covered interior cells evolve via `H_w += dt·(bmb_w − till_rate)`
     clamped to `[0, H_w_max]`. Floating + adjacent → H_w_max. Grounded
     ice-free → 0.
   - `K24`: `calc_k24(...)` → `q_x, q_y, N, p_w, H_w_eq`. If
     `update_H_w`: apply instant (`H_w := H_w_eq`) or relaxation (`H_w :=
     H_w + (H_w_eq − H_w)·(1 − exp(−dt/tau_H))`). No clamp to H_w_max.
6. If `method = BUCKET` and `par%bucket%N_closure ≠ NONE`: run closure
   post-step → write N; set `p_w = P_o − N`.
7. (For `method = K24`, ignore `par%bucket%N_closure`; K24's N is canonical.)
8. `dHwdt = (H_w − H_w_old) / dt` (zero if dt == 0).

## Mask semantics

- `mask` is host-supplied "active domain" (e.g. basin restriction).
- For BUCKET: mask is intersected with the 4-case logic — `mask = 0` skips the
  cell entirely (no override, no evolution).
- For K24: passes the raw host `mask` through to `calc_k24` (status quo).
  Host may pre-restrict to grounded ice if desired.

## K24 changes

- Add `H_w` (or `H_w_eq`) to `calc_k24`'s OUT signature. Source: existing
  `H_conduit` (computed at k24.f90:286), no new physics.
- Wrap in `hydro_update` per the `update_H_w` / `tau_H` logic.

## Dependencies & build

- Static library, Makefile mirroring yelmo's pattern.
- Hard link dependency on fesm-utils (FFTW, LIS); same `-I.../-L.../-lfesmutils`
  pattern as yelmo.
- No `#ifdef` flags.

## Test programs

- Single `shmip.x` driver in `tests/`.
- Case selector via CLI flag or namelist entry.
- Coverage: SHMIP A–F in scope; A–D required for initial commit (steady-state
  subset). E/F land later for transient validation.

## Closures ported from yelmo

All four ported verbatim from `yelmo/src/physics/basal_dragging.f90`:
`overburden`, `marine` (Leguy 2014, geometry-only — H_w extension deferred),
`till` (van Pelt & Bueler 2015), `two_value`. Parameters re-routed from
yelmo's namelist groups into `par%bucket%n_<closure>%(...)`.

## Implementation milestones

1. **API rename / signature growth.** Add `f_ice, f_grnd` to `hydro_update`
   and `hydro_init_state`. Drop `HYDRO_METHOD_RESERVED`. Replace `1.0e10`
   sentinel with `initialized` bool. Gate kappa init on `method = K24`.
2. **K24 H_w output.** Add `H_w` to `calc_k24` OUT args (= `H_conduit`).
   Wire `par%k24%update_H_w` + `par%k24%tau_H` in `hydro_update`. Add `dHwdt`
   to state and compute via stash-and-diff.
3. **Bucket model.** Implement `hydro_calc_bucket` mirroring yelmo's
   `calc_basal_water_local`. Wire under `case (HYDRO_METHOD_BUCKET)`. Add
   `par%bucket` sub-struct with `till_rate, H_w_max` (top-level), and
   `N_closure` enum.
4. **N closures.** Port the four yelmo subs as public standalone subs. Wire
   post-step dispatch in `hydro_update` when `method = BUCKET` and
   `N_closure ≠ NONE`. Closure-specific parameters in
   `par%bucket%n_<closure>`.
5. **Floating override on init.** Apply H_w = H_w_max on floating +
   adjacent-to-floating in `hydro_init_state`, independent of `init_method`.
6. **Build system.** Makefile linking against fesm-utils (FFTW, LIS).
   Produce `libfasthydro.a`.
7. **SHMIP A–D driver.** Single `shmip.x` test program with case selector.
   Run A–D against K24 and BUCKET, verify against reference outputs.
8. **SHMIP E–F.** Transient cases, after the steady set is locked.

Out of scope for v1:

- Yelmo integration (separate plan).
- H_w extension of the Leguy/marine closure.
- Precision parameterization (deferred to library-wide cross-cut).
- FFTW-optional builds.
- Per-method state sub-structs.
- Replacing FFT smoother with a direct convolution.

## Post-v1 follow-ups (landed)

- NetCDF output via fesm-utils `ncio` (no variable tables).
- SHMIP cases A1..A6, B1..B5, D1..D5 in the driver. Moulin positions
  for B-cases come from a fixed-seed Park-Miller LCG (SHMIP-style but
  not bit-exact). SHMIP C is intentionally not supported (real C
  resolution requires `dt ~ minutes`, out of scope).
- Domain-border BC `par%mask_bc` (+ `par%H_w_bc`) applied on the
  i=1, i=nx, j=1, j=ny rim, for any method that wrote H_w (BUCKET, K24):
  - `MASK_BC_ZERO`    -- H_w = 0 on the rim
  - `MASK_BC_IMPOSED` -- H_w = `par%H_w_bc` on the rim
  - `MASK_BC_MIRROR`  -- H_w on rim copied from inward neighbor (Neumann)
  Yelmo's floating-cell logic (H_w = H_w_max on floating + adjacent,
  H_w = 0 on grounded-ice-free) stays hardcoded, independent of `mask_bc`.

## Resolved: K24 H_w interpretation

K24's `H_w` was initially filled from `H_conduit` (a conduit length
scale, ~O(100 m)). Per the Julia reference
(TakisAngelides/FastHydrology.jl/water_flux.jl, `update_W!`), the
SHMIP-comparable sheet water-layer thickness is

  `W = (12 * eta_w * q_si / mean(|grad phi0_smoothed|))^(1/3)`

clamped to `[W_min, W_max]`. New k24 parameters `eta_w = 1.8e-3 Pa s`,
`W_min = 1e-8 m`, `W_max = 0.015 m` (defaults from Kazmierczak 2022).
SHMIP-A1 now produces H_w ~ 1.2 mm; A5 ~ 9 mm, scaling as melt^(1/3).
