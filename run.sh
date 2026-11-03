#!/bin/bash
#------------------------------------------------------------------------------
# Driver for the AE339 compression-ramp project (macOS, OpenFOAM-v2606.app).
#
#   ./run.sh            run the main case (rampCase/) and post-process it
#   ./run.sh study      run the mesh-convergence study (meshStudy/) and post-process
#   ./run.sh sweep      run the ramp-angle sweep (angleSweep/) and post-process
#   ./run.sh post       post-processing only
#   ./run.sh transient  run the start-up animation case (transientCase/, ~1 h)
#   ./run.sh movie      transient figures + MP4 animation (results/)
#   ./run.sh view       open the start-up animation in ParaView
#   ./run.sh clean      remove all solution data
#
# OpenFOAM rejects paths that contain spaces, so cases are run through a
# space-free symlink: ~/ae339_project -> this folder.
#------------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
LINK="$HOME/ae339_project"

if [ ! -e "$LINK" ]; then
    ln -s "$HERE" "$LINK"
elif [ "$(cd "$LINK" && pwd -P)" != "$(cd "$HERE" && pwd -P)" ]; then
    echo "error: $LINK exists but does not point to $HERE" >&2
    exit 1
fi

# ~/.OpenFOAM/prefs.sh on this machine selects SYSTEMOPENMPI, but the app's
# MPI library is built as OPENMPI. Override it, or parallel runs load the
# dummy (serial-only) Pstream library.
foam() { openfoam2606 -DWM_MPLIB=OPENMPI -c "$1"; }

PARAVIEW_APP=/Applications/ParaView-6.0.1.app/Contents

runCase() {
    echo "== $1"
    foam "cd '$LINK/$1' && ./Allclean > /dev/null && ./Allrun"
    if ! grep -q "^End" "$HERE/$1/log.rhoCentralFoam"; then
        echo "error: rhoCentralFoam did not finish, see $1/log.rhoCentralFoam" >&2
        exit 1
    fi
}

# Post-process the main case, plus the mesh study and angle sweep if present
postAll() {
    local args=()
    [ -d "$HERE/meshStudy" ] && args+=(--study "$HERE/meshStudy")
    [ -d "$HERE/angleSweep" ] && args+=(--sweep "$HERE/angleSweep")
    python3 "$HERE/scripts/postprocess.py" ${args[@]+"${args[@]}"}
}

case "${1:-main}" in
main)
    runCase rampCase
    postAll
    ;;
study)
    "$HERE/scripts/make_mesh_study.sh"
    for d in "$HERE"/meshStudy/*/; do
        runCase "meshStudy/$(basename "$d")"
    done
    postAll
    ;;
sweep)
    "$HERE/scripts/make_angle_sweep.sh"
    for d in "$HERE"/angleSweep/*/; do
        runCase "angleSweep/$(basename "$d")"
    done
    postAll
    ;;
post)
    postAll
    ;;
transient)
    runCase transientCase
    python3 "$HERE/scripts/transient_figures.py"
    ;;
movie)
    python3 "$HERE/scripts/transient_figures.py"
    "$PARAVIEW_APP/bin/pvbatch" "$HERE/scripts/view_transient.py" \
        --movie "$HERE/results/rampStartup.mp4"
    ;;
view)
    # Started from the project folder: the script falls back to the working
    # directory to find transientCase/ when ParaView does not set __file__.
    (cd "$HERE" && "$PARAVIEW_APP/MacOS/paraview" --script="$HERE/scripts/view_transient.py" &)
    ;;
clean)
    foam "cd '$LINK/rampCase' && ./Allclean"
    foam "cd '$LINK/transientCase' && ./Allclean"
    rm -rf "$HERE/meshStudy" "$HERE/angleSweep"
    ;;
*)
    sed -n '3,13p' "$0"
    exit 1
    ;;
esac
