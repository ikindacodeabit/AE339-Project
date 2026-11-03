#!/bin/bash
#------------------------------------------------------------------------------
# Create ramp-angle copies of rampCase/ in angleSweep/, one per ramp angle, all
# at the freestream Mach number of rampCase. Each copy differs from rampCase
# only in rampAngle and cellSize, and writes the final time only. They place
# several CFD points on the theta-beta-M curve. Run them with  ./run.sh sweep.
#
# Keep the angles below theta_max for the Mach number (28.7 deg at Mach 2.4):
# above it the shock detaches and there is no attached exact solution.
#------------------------------------------------------------------------------
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ANGLES="${ANGLES:-5 10 15 20 25}"   # [deg]
DX="${DX:-0.01}"                    # [m]; 10 mm: beta within ~0.1 deg

for th in $ANGLES; do
    dst="$ROOT/angleSweep/theta_${th}deg"
    rm -rf "$dst"
    mkdir -p "$dst"
    cp -R "$ROOT"/rampCase/{0.orig,constant,system,Allrun,Allclean} "$dst"/

    sed -i '' -E "s/^(rampAngle[[:space:]]+)[^;]+;/\1$th;/" "$dst/system/caseParameters"
    sed -i '' -E "s/^(cellSize[[:space:]]+)[^;]+;/\1$DX;/" "$dst/system/caseParameters"
    sed -i '' -E "s/^(writeInterval[[:space:]]+)[^;]+;/\10.01;/" "$dst/system/controlDict"
    sed -i '' -E "s/^numberOfSubdomains[[:space:]]+[0-9]+;/numberOfSubdomains 1;/" \
        "$dst/system/decomposeParDict"

    echo "angleSweep/theta_${th}deg: rampAngle $th deg, cellSize $DX m"
done
