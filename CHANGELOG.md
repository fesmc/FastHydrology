# Changelog

All notable changes to FastHydrology are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/); versioning is [SemVer](https://semver.org/).

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
