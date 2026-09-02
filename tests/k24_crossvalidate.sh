#!/bin/bash
# Cross-validate the Fortran K24 (src/k24.f90) against the kazmierczak2024
# model in FastHydrology.jl across every configuration the port supports.
# See tests/README.md for what each case covers and why the two harnesses
# differ in achievable precision.
#
#   usage: tests/k24_crossvalidate.sh /path/to/FastHydrology.jl [workdir]
#
# Run from the repository root, after:  make k24_synth k24_greenland
set -u

JLPROJ="${1:?usage: $0 /path/to/FastHydrology.jl [workdir]}"
WORK="${2:-/tmp/k24_crossvalidate}"
JL="julia --project=$JLPROJ"
NML_SRC="par/k24_greenland.nml"
INPUT="${K24_INPUT:-tests/k24_synth_input.nc}"   # override to run the same sweep on another dataset
TOL="${K24_TOL:-1e-12}"   # relative; override to tighten/loosen

mkdir -p "$WORK"

fail=0

# $1 = case name, $2 = julia overrides, $3.. = sed edits applied to the namelist
run_case () {
    name="$1"; jargs="$2"; shift 2
    nml="$WORK/$name.nml"
    cp "$NML_SRC" "$nml"
    # keep the per-iteration Picard chatter out of the log
    sed -i 's/^    k24_dissipation_verbose          = True/    k24_dissipation_verbose          = False/;
            s/^    k24_coupling_verbose             = True/    k24_coupling_verbose             = False/' "$nml"
    # NML_SRC (par/k24_greenland.nml) sets k24_long_coupling_water = 0.0, correct for its own
    # 16 km grid but not for SYNTH's much finer dx=2000/dy=3000 (k24_synth_gen.jl). Re-normalize
    # to the sweep's own documented baseline (tests/README.md: "the model's default
    # longcoupwater = 5.0") so this file's Greenland-specific value doesn't silently change what
    # every non-longcoup* case below tests. longcoup0/longcoup1 still override this explicitly.
    sed -i 's/^    k24_long_coupling_water          = 0.0/    k24_long_coupling_water          = 5.0/' "$nml"
    for e in "$@"; do sed -i "$e" "$nml"; done

    printf '  %-26s ' "$name"
    if ! ./bin/k24_synth.x "$nml" "$WORK/${name}_f.nc" "$INPUT" > "$WORK/${name}_f.log" 2>&1; then
        echo "FORTRAN FAILED (see $WORK/${name}_f.log)"; fail=1; return
    fi
    if ! $JL tests/k24_synth.jl "$WORK/${name}_j.nc" input="$INPUT" $jargs > "$WORK/${name}_j.log" 2>&1; then
        echo "JULIA FAILED (see $WORK/${name}_j.log)"; fail=1; return
    fi
    # k24_synth_compare.jl columns: field  max|J|  max|F-J|  rel_max  rel_max_all_cells
    $JL tests/k24_synth_compare.jl "$WORK/${name}_f.nc" "$WORK/${name}_j.nc" 2>&1 \
        | awk -v tol="$TOL" '
            /^(W|N|q|p_w) /{
                printf "%s=%s ", $1, $4
                if ($4+0 > tol || $5+0 > tol) bad=1
            }
            END{ if (nomatch) bad=1; printf "%s\n", bad ? "<-- ABOVE TOLERANCE" : "ok"; exit bad }'
    if [ "${PIPESTATUS[1]}" -ne 0 ]; then fail=1; fi
    return 0
}

WT='s/^    k24_water_thickness_algorithm    = 0/    k24_water_thickness_algorithm    = '
GC='s/^    k24_gradient_convention          = 0/    k24_gradient_convention          = '
FS='s/^    k24_flux_solver                  = 0/    k24_flux_solver                  = '
LC='s/^    k24_long_coupling_water          = 5.0/    k24_long_coupling_water          = '
DM='s/^    k24_dissipation_melt             = True/    k24_dissipation_melt             = '
SL='s/^    k24_sliding_law                  = 0/    k24_sliding_law                  = '
ST='s/^    k24_substrate_type               = 0/    k24_substrate_type               = '
DR='s/^    k24_drainage_mode                = 0/    k24_drainage_mode                = '

echo "Case sweep on $INPUT (double precision throughout; expect ~1e-15)"
echo
run_case hard_both      ""
run_case soft_both      "substrate=soft"                        "${ST}1/"
run_case soft_efficient "substrate=soft drainage=efficient"     "${ST}1/" "${DR}1/"
run_case soft_ineff     "substrate=soft drainage=inefficient"   "${ST}1/" "${DR}2/"
run_case mixed_both     "substrate=mixed"                       "${ST}2/"
run_case mixed_eff      "substrate=mixed drainage=efficient"    "${ST}2/" "${DR}1/"
run_case hard_ineff     "drainage=inefficient"                  "${DR}2/"
run_case laminar_mean   "wthick=laminar"                        "${WT}1/"
run_case laminar_local  "wthick=laminar grad=local"             "${WT}1/" "${GC}1/"
run_case areal          "wthick=areal"                          "${WT}2/"
run_case darcy_local    "grad=local"                            "${GC}1/"
run_case longcoup0      "longcoup=0.0"                          "${LC}0.0/"
run_case longcoup1      "longcoup=1.0"                          "${LC}1.0/"
run_case iterative      "psi=iterative"                         "${FS}1/"
run_case nodissip       "dissip=false"                          "${DM}False/"
run_case weertman       "sliding=weertman"                      "${SL}1/" \
    's/^    k24_weertman_C                   = 0.0/    k24_weertman_C                   = 1.0e5/'
run_case powerplastic   "sliding=powerplastic ctill=0.02"       "${SL}2/" \
    's/^    k24_power_plastic_c_till         = 0.0/    k24_power_plastic_c_till         = 0.02/'
run_case regcoulomb     "sliding=regcoulomb ctill=0.02"         "${SL}3/" \
    's/^    k24_reg_coulomb_c_till           = 0.0/    k24_reg_coulomb_c_till           = 0.02/'

echo
echo "Greenland 16 km case (through the real(sp) public API; expect ~1e-7)"
echo
SCALE=2.90585971523335349e-08     # (rho_ice/rho_w)/SEC_PER_YEAR
# NML_SRC's own k24_long_coupling_water = 0.0 (correct at 16 km) is used as-is here, unlike
# in run_case above. tests/k24_greenland.jl's own longcoup default is independently hardcoded
# to "5.0" (get_cfg("longcoup", "5.0")), so it needs the matching override explicitly or the
# two sides would silently compare different physics.
if ./bin/k24_greenland.x "$NML_SRC" "$WORK/grl_f.nc" $SCALE > "$WORK/grl_f.log" 2>&1 &&
   $JL tests/k24_greenland.jl "$WORK/grl_j.nc" mdotscale=$SCALE longcoup=0.0 > "$WORK/grl_j.log" 2>&1; then
    $JL tests/k24_greenland_compare.jl "$WORK/grl_f.nc" "$WORK/grl_j.nc"
else
    echo "  Greenland case FAILED (see $WORK/grl_f.log, $WORK/grl_j.log)"; fail=1
fi

echo
if [ $fail -eq 0 ]; then echo "cross-validation PASSED"; else echo "cross-validation FAILED"; fi
exit $fail
