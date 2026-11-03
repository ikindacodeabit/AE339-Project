#!/bin/bash
#------------------------------------------------------------------------------
# Create mesh-convergence copies of rampCase/ in meshStudy/, one per cell size.
# Each copy differs from rampCase only in cellSize (and writes the final time
# only, to save disk). Run them with  ./run.sh study.
#------------------------------------------------------------------------------
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SIZES="${SIZES:-0.02 0.01 0.005}"   # [m]; each level halves the cell size

for dx in $SIZES; do
    mm=$(python3 -c "print(f'{$dx*1e3:g}')")
    dst="$ROOT/meshStudy/dx_${mm}mm"
    rm -rf "$dst"
    mkdir -p "$dst"
    cp -R "$ROOT"/rampCase/{0.orig,constant,system,Allrun,Allclean} "$dst"/

    sed -i '' -E "s/^(cellSize[[:space:]]+)[^;]+;/\1$dx;/" "$dst/system/caseParameters"
    sed -i '' -E "s/^(writeInterval[[:space:]]+)[^;]+;/\10.01;/" "$dst/system/controlDict"

    # Parallel overhead outweighs the gain below ~50k cells
    nCells=$(python3 -c "print(round(1.25/$dx)*round(2.25/$dx))")
    np=$([ "$nCells" -lt 50000 ] && echo 1 || echo 4)
    sed -i '' -E "s/^numberOfSubdomains[[:space:]]+[0-9]+;/numberOfSubdomains $np;/" \
        "$dst/system/decomposeParDict"

    echo "meshStudy/dx_${mm}mm: cellSize $dx m, ~$nCells cells, $np core(s)"
done
