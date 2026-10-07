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

# Portable in-place sed: BSD sed (macOS) requires an explicit backup suffix
# after -i (even empty), or it silently treats the next argument as that
# suffix instead of as the script -- which then leaves the script string
# itself unconsumed and mistaken for the file operand. Always write a .bak
# and remove it, since that's the one syntax both BSD and GNU sed accept.
sed_i() {
    sed -i.bak "$1" "$2" && rm -f "$2.bak"
}

mkdir -p "$WORK"

fail=0

# $1 = case name, $2 = julia overrides, $3.. = sed edits applied to the namelist
run_case () {
    name="$1"; jargs="$2"; shift 2
    nml="$WORK/$name.nml"
    cp "$NML_SRC" "$nml"
    # keep the per-iteration Picard chatter out of the log
    sed_i 's/^    k24_dissipation_verbose          = True/    k24_dissipation_verbose          = False/;
            s/^    k24_qN_verbose                   = True/    k24_qN_verbose                   = False/' "$nml"
    # NML_SRC (par/k24_greenland.nml) sets k24_coupling_length_kamb86 = 0.0, correct for its own
    # 16 km grid but not for SYNTH's much finer dx=2000/dy=3000 (k24_synth_gen.jl). Re-normalize
    # to the model's default 10.0 so this file's Greenland-specific value doesn't silently change
    # what every non-kamb86* case below tests. kamb86_0/kamb86_4 still override this explicitly.
    sed_i 's/^    k24_coupling_length_kamb86       = 0.0 /    k24_coupling_length_kamb86       = 10.0/' "$nml"
    for e in "$@"; do sed_i "$e" "$nml"; done

    printf '  %-26s ' "$name"
    # The Fortran driver reads qtfrac/iebfrac from the same key=value list and ignores the rest.
    if ! ./bin/k24_synth.x "$nml" "$WORK/${name}_f.nc" "$INPUT" $jargs > "$WORK/${name}_f.log" 2>&1; then
        echo "FORTRAN FAILED (see $WORK/${name}_f.log)"; fail=1; return
    fi
    if ! $JL tests/k24_synth.jl "$WORK/${name}_j.nc" input="$INPUT" $jargs > "$WORK/${name}_j.log" 2>&1; then
        echo "JULIA FAILED (see $WORK/${name}_j.log)"; fail=1; return
    fi
    # k24_synth_compare.jl columns: field  max|J|  max|F-J|  rel_max  rel_max_all_cells
    fields=phi0,absg,gsx,gsy,absgs,q,W,N,p_w,Q_b,Q_diss,C_frz
    case "$jargs" in *ubfac=*) fields=$fields,N_ub ;; esac
    $JL tests/k24_synth_compare.jl "$WORK/${name}_f.nc" "$WORK/${name}_j.nc" $fields 2>&1 \
        | awk -v tol="$TOL" '
            /^(W|N|N_ub|q|p_w|Q_b|Q_diss|C_frz) /{
                printf "%s=%s ", $1, $4
                if ($4+0 > tol || $5+0 > tol) bad=1
            }
            END{ if (nomatch) bad=1; printf "%s\n", bad ? "<-- ABOVE TOLERANCE" : "ok"; exit bad }'
    if [ "${PIPESTATUS[1]}" -ne 0 ]; then fail=1; fi
    return 0
}

WT='s/^    k24_water_thickness_algorithm    = 0/    k24_water_thickness_algorithm    = '
GC='s/^    k24_gradient_convention          = 0/    k24_gradient_convention          = '
FS='s/^    k24_flux_solver                  = 3/    k24_flux_solver                  = '
LC='s/^    k24_coupling_length_kamb86       = 10.0/    k24_coupling_length_kamb86       = '
RT='s/^    k24_routing_scheme               = 0/    k24_routing_scheme               = '
FA='s/^    k24_fill_algorithm               = -1/    k24_fill_algorithm               = '
QC='s/^    k24_q_conversion                 = -1/    k24_q_conversion                 = '
DD='s/^    k24_dissipation_discretization   = -1/    k24_dissipation_discretization   = '
FD='s/^    k24_friction_discretization      = 0/    k24_friction_discretization      = '
QO='s/^    k24_quinn_original               = False/    k24_quinn_original               = True/'
FQ='s/^    k24_friction_quadrature          = False/    k24_friction_quadrature          = True/'
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
run_case kamb86_0       "kamb86=0.0"                            "${LC}0.0/"
run_case kamb86_4       "kamb86=4.0"                            "${LC}4.0/"
run_case nodissip       "dissip=false"                          "${DM}False/"
run_case weertman       "sliding=weertman"                      "${SL}1/" \
    's/^    k24_weertman_C                   = 0.0/    k24_weertman_C                   = 1.0e5/'
run_case powerplastic   "sliding=powerplastic ctill=0.02"       "${SL}2/" \
    's/^    k24_power_plastic_c_till         = 0.0/    k24_power_plastic_c_till         = 0.02/'
run_case regcoulomb     "sliding=regcoulomb ctill=0.02"         "${SL}3/" \
    's/^    k24_reg_coulomb_c_till           = 0.0/    k24_reg_coulomb_c_till           = 0.02/'
run_case field          "sliding=field"                         "${SL}4/"
run_case regcoulfield   "sliding=regcoulombfield ctill=0.02"    "${SL}5/" \
    's/^    k24_reg_coulomb_c_till           = 0.0/    k24_reg_coulomb_c_till           = 0.02/'
run_case shakti         "sliding=shakti ctill=0.02"             "${SL}6/" \
    's/^    k24_shakti_C                     = 0.0/    k24_shakti_C                     = 0.02/'
run_case terms_qT_ieb   "qtfrac=0.3 iebfrac=0.5"

# N from a new sliding speed with the routing held (k24_N_from_ub / N_from_ub!): N_ub is N for
# 0.2 x and 5 x the sliding speed of the solve, compared like every other field.
run_case ub_slower      "sliding=field ubfac=0.2"               "${SL}4/"
run_case ub_faster      "sliding=field ubfac=5.0"               "${SL}4/"
run_case ub_nofric      "ubfac=5.0"
run_case ub_soft        "sliding=field substrate=soft ubfac=5.0" "${SL}4/" "${ST}1/"
run_case ub_mixed_ineff "sliding=field substrate=mixed drainage=inefficient ubfac=5.0" "${SL}4/" "${ST}2/" "${DR}2/"
run_case ub_qT_ieb      "sliding=field qtfrac=0.3 iebfrac=0.5 ubfac=5.0" "${SL}4/"
run_case stagger        "sliding=field friction=staggered"      "${SL}4/" "${FD}1/"
run_case stagger_quad   "sliding=field friction=staggeredquad"  "${SL}4/" "${FD}1/" "$FQ"
run_case gdswarner      "routing=gdswarner"                     "${RT}1/"
run_case gds_recursive  "routing=gdswarner psi=recursive"       "${RT}1/" "${FS}0/"
run_case gds_iterative  "routing=gdswarner psi=iterative"       "${RT}1/" "${FS}1/"
run_case gds_face       "routing=gdswarner qconv=face dissdisc=face" "${RT}1/" "${QC}1/" "${DD}1/"
run_case quinn          "routing=quinn"                         "${RT}2/"
run_case quinn_orig     "routing=quinnorig"                     "${RT}2/" "$QO"
run_case tarboton       "routing=tarboton"                      "${RT}3/"
run_case modtarboton    "routing=modtarboton"                   "${RT}4/"
run_case gdstarboton    "routing=gdstarboton"                   "${RT}5/"
run_case warner_jacobi  "fill=jacobi"                           "${FA}0/"
run_case warner_lowest  "fill=lowest"                           "${FA}1/"
run_case warner_outflow "qconv=outflow dissdisc=cell"           "${QC}0/" "${DD}0/"
run_case warner_nodiss  "dissip=false"                          "${DM}False/"

echo
echo "Greenland 16 km case (through the real(sp) public API; expect ~1e-7)"
echo
SCALE=2.90585971523335349e-08     # (rho_ice/rho_w)/SEC_PER_YEAR
# NML_SRC's own k24_coupling_length_kamb86 = 0.0 (correct at 16 km) is used as-is here, unlike
# in run_case above. tests/k24_greenland.jl's own kamb86 default is 10.0, so it needs the
# matching override explicitly or the two sides would silently compare different physics.
if ./bin/k24_greenland.x "$NML_SRC" "$WORK/grl_f.nc" $SCALE > "$WORK/grl_f.log" 2>&1 &&
   $JL tests/k24_greenland.jl "$WORK/grl_j.nc" mdotscale=$SCALE kamb86=0.0 > "$WORK/grl_j.log" 2>&1; then
    $JL tests/k24_greenland_compare.jl "$WORK/grl_f.nc" "$WORK/grl_j.nc"
else
    echo "  Greenland case FAILED (see $WORK/grl_f.log, $WORK/grl_j.log)"; fail=1
fi

echo
if [ $fail -eq 0 ]; then echo "cross-validation PASSED"; else echo "cross-validation FAILED"; fi
exit $fail
