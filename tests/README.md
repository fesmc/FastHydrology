# Tests

## `shmip.f90`

The SHMIP-style bucket/transport driver. Build and run:

```sh
make shmip
./bin/shmip.x par/shmip.nml
```

## Running the whole cross-validation

`k24_crossvalidate.sh` drives both K24 harnesses below over every supported
configuration and checks each against a tolerance:

```sh
make k24_synth k24_greenland
tests/k24_crossvalidate.sh /path/to/FastHydrology.jl
```

It exits non-zero if any case exceeds `K24_TOL` (default `1e-12`, comfortably
above the ~1e-15 round-off the synthetic case actually achieves). The
Greenland case is reported but not tolerance-checked, since its `real(sp)`
API caps it at ~1e-7.

`K24_INPUT` points the sweep at a different dataset. Both `k24_synth.x` and
`tests/k24_synth.jl` take the grid size and spacing from the input file, so any
NetCDF with `xc`, `yc` and the fields `h`, `b`, `mask`, `vb`, `A`, `mdot`
(water-equivalent volume rate, m/s) works:

```sh
K24_INPUT=/tmp/thwaites_input.nc tests/k24_crossvalidate.sh /path/to/FastHydrology.jl
```

### Thwaites 2 km

`k24_thwaites_export.jl` converts the Kazmierczak et al. 2024 Thwaites dataset
that ships with FastHydrology.jl into that format:

```sh
julia --project=/path/to/FastHydrology.jl tests/k24_thwaites_export.jl \
  "/path/to/FastHydrology.jl/test/Kazmierczak et al 2024/input/Kazmierczak2024/THWAITES2km_m3_HAB_toto.mat" \
  /tmp/thwaites_input.nc
```

345 x 288 with 51945 grounded cells — 25x Greenland-16km and a much harder
test of the flow routing. All 18 configurations agree to ~1e-13 or better.
The `.nc` is ~4.8 MB and derived from FastHydrology.jl's own test data, so it
is generated rather than committed.

Two things this dataset exercises that the others do not:

- **The `max_psi_out_calls` cutoff actually binds.** 51945 grounded cells
  exceeds the default 100000, so the routing is truncated part-way. Raising
  `k24_max_psi_out_calls` to 400000 changes the answer materially (mean `q`
  2.1409e-4 -> 2.1763e-4, mean `N` 1.3807e6 -> 1.3119e6) — and both
  implementations move together to 10 significant figures, so the truncation
  itself lands identically. If you raise it on one side you must raise it on
  the other, or you are comparing two different models.
- **The flow-direction graph is genuinely cyclic.** The topological router
  leaves 26938 of 51945 grounded cells (51.9%) unprocessed, matching
  `TopologicalPsiOut`'s own docstring ("roughly half of all grounded cells at
  the old default `longcoupwater = 5.0`, i.e. `coupling_length_kamb86 = 10`"). Its answer is therefore very
  different from the recursive/iterative one (mean `q` 1.28e-5 vs 2.14e-4) —
  both implementations agree on it exactly, but that agreement is not a reason
  to use it on real data.

## `k24_synth.f90` — K24 cross-validation against FastHydrology.jl

The Fortran K24 model in `src/k24.f90` is a port of the `kazmierczak2024`
model in [FastHydrology.jl](https://github.com/TakisAngelides/FastHydrology.jl).
This test runs both on identical inputs and compares them field by field, so
the port can be re-verified whenever either side changes.

It is a *synthetic* case on purpose. The real Greenland-16km restart in
`input/` cannot reach several code paths: on that dataset `N_inf` saturates at
`Po` on 100% of grounded cells, which hides the `kappa` hard/soft blend, the
drainage modes and the unclamped `N_inf` branch entirely. This case picks
`A_glen` so `N_inf` lands strictly inside `(sigmat*Po, Po)` on most cells (about
14% still clamp, so both branches are covered), and uses `dx != dy` so the
anisotropic `corfac` and the rectangular smoothing kernel are exercised too.

Both sides read their fields from `k24_synth_input.nc` rather than each
building them from the same formulas. Building them independently makes them
differ by ~1 ulp — Intel's `libm` and Julia's do not round `sin`/`cos`
identically — and `potential_filling`'s discrete local-minimum test amplifies
that into a large local difference. That is an artifact of the harness, not of
either model, and reading a shared file removes it.

`k24_synth.f90` calls `calc_k24` directly rather than going through
`hydro_update`, so everything stays in double precision: the public
`hydro_class` state is `real(sp)`, which would cap agreement at ~1e-7.

### Running it

```sh
make k24_synth
./bin/k24_synth.x <namelist> out_fortran.nc [input.nc]
```

The Julia side needs a FastHydrology.jl checkout:

```sh
julia --project=/path/to/FastHydrology.jl tests/k24_synth.jl out_julia.nc
julia --project=/path/to/FastHydrology.jl tests/k24_synth_compare.jl out_fortran.nc out_julia.nc
```

`tests/k24_synth.jl` takes `key=value` overrides matching the namelist
switches: `substrate`, `drainage`, `wthick`, `grad`, `kamb86`, `routing`,
`fill`, `qconv`, `dissdisc`, `friction`, `psi`, `dissip`, `sliding`, `ctill`,
`input`, plus `qtfrac`/`iebfrac` (the source terms), which `k24_synth.x` also
reads from its own arguments. Set the corresponding `k24_*` keys in the
namelist to the same values and the two must agree.

The water source is built from terms on both sides from the fixture's volume
rate `mdot`: `G = mdot*rho_w*L_w`, `q_T = qtfrac*G`, `i_eb = iebfrac*mdot*rho_w`.
The fields some options need (C-grid velocities, a `tau_b` field, a per-cell
Coulomb coefficient) are derived from `vb` the same way on both sides.

To regenerate the input fixture (needs Julia only):

```sh
julia --project=/path/to/FastHydrology.jl tests/k24_synth_gen.jl tests/k24_synth_input.nc 1e-24
```

### Expected result

Every configuration agrees to ~1e-15 relative (double-precision round-off) on
`W`, `N`, `q` and `p_w`, and on the optional intermediate diagnostics
(`phi0`, `abs_g`, `gsx`, `gsy`, `abs_gs`) that `calc_k24` can return. Anything
materially above that is a regression. The configurations covered:

| case | what it exercises |
| --- | --- |
| `hard_both` | baseline: hard bed, both drainage terms |
| `soft_both`, `mixed_both` | the `kappa` hard/soft conduit blend |
| `soft_efficient`, `mixed_eff` | `EfficientOnly`: `Q_c -> 0`, sliding opening term dropped |
| `soft_ineff`, `hard_ineff` | `InefficientOnly`: `Q_c -> Inf`, melt opening term dropped |
| `laminar_mean`, `laminar_local` | the laminar `W` closure, both gradient conventions |
| `areal` | the areal-conduit `W` closure |
| `darcy_local` | Darcy-Weisbach `W` with the local gradient |
| `kamb86_0` | smoothing disabled (bypasses the FFT convolution entirely) |
| `kamb86_4` | a different smoothing kernel size |
| `nodissip`, `warner_nodiss` | dissipation melt off (single routing pass, no Picard loop) |
| `weertman`, `field` | N-independent sliding laws (the second a prescribed `tau_b` field) |
| `powerplastic`, `regcoulomb`, `regcoulfield`, `shakti` | N-dependent sliding laws, i.e. the joint (q, N) Picard loop |
| `terms_qT_ieb` | a nonzero `q_T` and water from above `i_eb` in the source |
| `stagger`, `stagger_quad` | C-grid frictional heat, on faces and at Gauss points |
| `gdswarner`, `gds_recursive`, `gds_iterative`, `gds_face` | the original K24 routing through each flow router, and with face fluxes |
| `quinn`, `quinn_orig`, `tarboton`, `modtarboton`, `gdstarboton` | the other Le Brocq et al. (2006) routing schemes |
| `warner_jacobi`, `warner_lowest`, `warner_outflow` | the default Warner routing with the other fill algorithms and `q` conversion |

Every case except the `gds*` ones runs the default Warner routing with
priority-flood filling, face-average `q` and face dissipation.

## `k24_greenland.f90` — the same cross-validation on real data

Complements the synthetic case above by running both implementations on the
Greenland-16km restart already in `input/`. It covers two branches the
synthetic case does not reach — `S_inf == 0` where `Q == 0` (about 4800 cells
here) and the `N_inf == Po` clamp (100% of grounded cells on this dataset) —
and is the reason the synthetic case exists: those same properties are what
make this one blind to the `kappa` blend and the drainage modes.

It goes through the public `hydro_update` API rather than calling `calc_k24`
directly, so agreement is capped at ~1e-7 by the `real(sp)` state arrays. That
is expected; the synthetic case is the one that resolves to 1e-15.

```sh
make k24_greenland
./bin/k24_greenland.x par/k24_greenland.nml out_fortran.nc 2.90585971523335349e-08
julia --project=/path/to/FastHydrology.jl tests/k24_greenland.jl out_julia.nc \
      mdotscale=2.90585971523335349e-08
julia --project=/path/to/FastHydrology.jl tests/k24_greenland_compare.jl out_fortran.nc out_julia.nc
```

The third argument (and `mdotscale=`) is `(rho_ice/rho_w)/sec_year`, which
converts the restart's `bmb_grnd` from ice-equivalent m/a to a water-equivalent
m/s source rate. Both sides must be given the same value. Passing `1.0`
instead reproduces `examples/greenland/greenland.jl`'s own (unconverted, and
therefore ~3e7 too large) forcing.

Expect ~1e-7 relative agreement on `W`, `N`, `q` and `p_w`, with correlations
of 1.0000000000.
