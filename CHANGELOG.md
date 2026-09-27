# Changelog

All notable changes to FastHydrology are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/); versioning is [SemVer](https://semver.org/).

## [Unreleased]
### Added
- Periodic domains: `hydro_init` takes optional `periodic_x`/`periodic_y` (default `.false.`,
  which leaves every result unchanged). A periodic direction wraps with period `nx`/`ny` and no
  halo cells: `mask_bc` is not applied to its rim, and the neighbours in `apply_margin_fill` and
  in every K24 stencil (potential filling, potential gradients, smoothing padding, flow routing)
  wrap instead of being clamped or skipped. `apply_margin_fill`/`apply_mask_bc` take the same
  optional flags; K24 carries them in `k24_param_class` (set by `hydro_init`, not the namelist).
- K24: dissipation-melt source term `|q*grad(phi0)|/(L_w*rho_w)`, resolved by a Picard loop
  (`k24_dissipation_melt`, on by default; `k24_max_dissipation_iters`, `k24_dissipation_rtol`).
- K24: four sliding laws for the frictional-heating term `tau_b*v_b/(L_w*rho_w)`
  (`k24_sliding_law`: none/Weertman/power-plastic/regularized-Coulomb). The two
  pressure-dependent laws make `tau_b` depend on `N`, which turns `q` and `N` into a joint
  fixed point solved by a coupled Picard loop (`k24_max_coupling_iters`, `k24_coupling_rtol`).
- K24: three water-thickness closures (`k24_water_thickness_algorithm`) — Darcy-Weisbach
  (new default), laminar Le Brocq/Weertman (the previous behaviour), areal-conduit — each
  selectable between a domain-mean and a local potential gradient (`k24_gradient_convention`).
- K24: drainage modes (`k24_drainage_mode`) forcing entirely-efficient or
  entirely-inefficient drainage, gating both `update_H`'s `Q_c` and `N_inf`'s opening terms.
- K24: a third flow router, `k24_flux_solver = 1` (iterative, explicit-stack), which
  reproduces the recursive one exactly without its recursion-depth limit.
- K24: `k24_max_psi_out_calls` cap on cells visited per routing sweep, and
  `k24_toposort_allow_cycles` to downgrade the topological router's cycle error to a warning.
- K24: `k24_q_min`/`k24_q_max`/`k24_fill_iters` exposed (they were hard-coded).
- `hyd%now%q`, the distributed flux magnitude, persisted between steps so the Picard loops
  warm-start from it.
- `input/yelmo_defaults.nml`: the canonical `&yhyd` schema, so the library's own examples
  and tests are runnable without a Yelmo checkout.
- `tests/k24_synth.f90` and `tests/k24_greenland.f90`: cross-validation harnesses that run
  the Fortran K24 and FastHydrology.jl on identical inputs and compare them field by field.
  See `tests/README.md`.

### Changed
- K24 is now a port of the `kazmierczak2024` model in FastHydrology.jl, verified against it
  to ~1e-15 (double precision, 18 configurations) and ~1e-7 through the `real(sp)` public API.
- **Numerics (results change).** `corfac` is now anisotropic,
  `(|gsx|*dy + |gsy|*dx)/|grad_s|`, replacing a `sqrt(dx*dy)` placeholder that was only
  correct for square cells. `psi_out` seeds unclamped and applies `max(0, .)` after
  accumulation instead of at seeding. `Po` is built from the raw ice thickness rather than
  the potential-filled one. Potential-filling and the gradient stencil now cover the domain
  border with edge-clamped indices instead of skipping/copying it. The smoothing kernel is
  sized per axis, and `k24_long_coupling_water = 0` now bypasses smoothing entirely.
  `S_inf`/`H`/`N_inf` resolve their degenerate cases the way the reference does (`Q == 0` and
  `S_inf == 0` give exactly 0 and `Po`) rather than flooring at `1e-12`; `H_soft` gained a
  `max(0, .)`.
- **Defaults (results change).** `k24_latent_heat_water` 3.35e5 -> 3.34e5.
  `k24_eta_w` 1.8e-3 -> 5.7039e-11: KORI-ULB's `par.waterviscosity` is a per-year quantity
  and needed the same per-year -> per-second conversion as every other per-year input, so the
  old value was ~3.16e7 too large. The four KORI-ULB clamps are off by default, matching the
  reference constructor: `k24_W_min` 1e-8 -> 0, `k24_W_max` 0.015 -> inf,
  `k24_min_pressure_fraction` 0.02 -> 0, and the new `k24_q_max` -> inf. Pass the old values
  explicitly to restore them.
- **Breaking:** `k24_flux_solver = 1` now selects the iterative router; the topological one
  moved to `2`, matching the reference's ordering.
- **Breaking:** `calc_k24` takes `q` and `N` as `INOUT` (plus optional diagnostic outputs).
- Example namelists use the `&yhyd` group name. They still said `&fhyd`, which was renamed in
  v0.2, so every setting in them had been silently ignored in favour of the schema defaults.

### Fixed
- K24 smoothing convolution: the border was reflected where the reference replicates, and the
  result was cropped at `(i + 2*frb, j + 2*frb)` instead of `(i + frb, j + frb)`, translating
  the smoothed potential gradient by `frb` cells in both directions.
- K24 `N_inf`: the drainage mode's opening coefficients were missing, so entirely-efficient
  and entirely-inefficient drainage did not actually drop their respective opening terms.

## [v0.3] - 2026-07-15
### Changed
- Integrate with the flattened fesm-utils layout (now built against fesm-utils v1.3): update include/library paths in the build config.
- Document fesm-utils compatibility (rebuild against dev ≥ `3f415cc`).

## [v0.2] - 2026-06-18
### Added
- Bucket hydrology model, K24 sheet/conduit model, and multi-method dispatch.
- Configurable N (effective pressure) closures, including constant-N (`N_CLOSURE_CONST`).
- Margin-fill closure for ice-margin cells.
- SHMIP A/B/C/D driver cases and NetCDF output.
- Greenland 16 km end-to-end example from a Yelmo restart, with plotting scripts.
- `runme`/`configme` scaffolding and namelist API for Yelmo coupling.

### Changed
- Switch internals to SI units; carry `dx`/`dy` separately.
- Adopt BvP15 notation; sequential bucket→K24 coupling.
- Rename namelist group `fhyd`→`yhyd`; lift physical constants (rho, g) to top level.
- Thread `input/yelmo_defaults.nml` through every `par_load`.

### Fixed
- K24 `psi_out` uninitialized-memory bug at masked cells.
- K24 H_w interpretation (Kazmierczak 2024 Eq. 8, sheet thickness W).
