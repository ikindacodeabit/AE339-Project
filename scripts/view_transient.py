"""
ParaView set-up for the start-up animation of the supersonic ramp flow.

Interactive (sets up a 2x2 layout; then press Play):
    /Applications/ParaView-6.0.1.app/Contents/MacOS/paraview --script=scripts/view_transient.py

Batch (no window):
    PVBATCH=/Applications/ParaView-6.0.1.app/Contents/bin/pvbatch
    $PVBATCH scripts/view_transient.py --movie results/rampStartup.mp4   # needs ffmpeg
    $PVBATCH scripts/view_transient.py --frames results/frames           # PNG sequence
    $PVBATCH scripts/view_transient.py --still 0.0015 out.png            # one time

Panels: Mach number, p/p_inf, T/T_inf and numerical schlieren exp(-k|grad rho|/ref).
Colour ranges are fixed over the whole animation so frames are comparable.
"""

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from paraview.simple import *  # noqa: F401,F403

try:
    ROOT = Path(__file__).resolve().parent.parent
except NameError:  # the GUI's --script does not always define __file__
    ROOT = Path.cwd()

ap = argparse.ArgumentParser()
ap.add_argument("--case", type=Path, default=ROOT / "transientCase")
ap.add_argument("--frames", type=Path, help="write PNG frames of the whole animation here")
ap.add_argument("--movie", type=Path, help="write an H.264 MP4 of the whole animation (uses ffmpeg)")
ap.add_argument("--fps", type=int, default=20, help="movie frame rate")
ap.add_argument("--still", nargs=2, metavar=("TIME", "PNG"), help="render one time to a PNG")
ap.add_argument("--schlieren-ref", type=float, default=25.0,
                help="|grad rho| [kg/m^4] that maps to mid-grey in the schlieren panel "
                     "(~25 for the 2.5 mm mesh; scale with 1/cell size)")
args, _ = ap.parse_known_args(sys.argv[1:])

case = args.case.resolve()
foam = case / f"{case.name}.foam"
if not foam.exists():
    foam = next(case.glob("*.foam"))

prm = {k: float(v) for k, v in re.findall(
    r"^\s*(\w+)\s+([-+0-9.eE]+)\s*;", (case / "system" / "caseParameters").read_text(), re.M)}
pInf, TInf = prm["pInf"], prm["TInf"]

# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

LoadPalette(paletteName="WhiteBackground")

reader = OpenFOAMReader(registrationName=foam.name, FileName=str(foam))
reader.MeshRegions = ["internalMesh"]
reader.CellArrays = ["Ma", "p", "T", "U", "rho", "magGradRho"]
reader.UpdatePipelineInformation()

calc = PythonCalculator(registrationName="nondimensional", Input=reader)
calc.ArrayAssociation = "Cell Data"
calc.Expression = f"inputs[0].CellData['p']/{pInf}"
calc.ArrayName = "p_ratio"

calcT = PythonCalculator(registrationName="T_ratio", Input=calc)
calcT.ArrayAssociation = "Cell Data"
calcT.Expression = f"inputs[0].CellData['T']/{TInf}"
calcT.ArrayName = "T_ratio"

schl = PythonCalculator(registrationName="schlieren", Input=calcT)
schl.ArrayAssociation = "Cell Data"
schl.Expression = f"exp(-0.693*inputs[0].CellData['magGradRho']/{args.schlieren_ref})"
schl.ArrayName = "schlieren"

# Interpolate to points for smooth shading
pts = CellDatatoPointData(registrationName="smooth", Input=schl)

# ---------------------------------------------------------------------------
# Views
# ---------------------------------------------------------------------------

MaMax = round(prm["Mach"] + 0.1, 1)  # just above the freestream Mach number
# p and T ranges cover the M = 2.4 start-up: p/p_inf 0.78-7.0, T/T_inf 0.97-2.01

PANELS = [
    # array,      title,                                          preset,      range,       log
    ("Ma",        "Mach number",                                  "Blues",     (0.0, MaMax), False),
    ("p_ratio",   r"Pressure  $p/p_\infty$ (log scale)",          "Oranges",   (0.7, 8.0),  True),
    ("T_ratio",   r"Temperature  $T/T_\infty$",                   "Reds",      (0.95, 2.05), False),
    ("schlieren", r"Numerical schlieren  $\Vert\nabla\rho\Vert$", "Grayscale", (0.0, 1.0),  False),
]

layout = CreateLayout("Ramp start-up")
views = []
for i, (array, title, preset, rng, log) in enumerate(PANELS):
    v = CreateView("RenderView")
    v.OrientationAxesVisibility = 0
    v.CameraParallelProjection = 1
    v.InteractionMode = "2D"
    views.append(v)

layout.SplitVertical(0, 0.5)
layout.SplitHorizontal(1, 0.5)
layout.SplitHorizontal(2, 0.5)
for cell, v in zip([3, 4, 5, 6], views):
    layout.AssignView(cell, v)

xc = 0.5 * (-prm["xUp"] + prm["rampLength"] + prm["xDown"])
yc = 0.5 * prm["height"]

for v, (array, title, preset, rng, log) in zip(views, PANELS):
    d = Show(pts, v)
    d.Representation = "Surface"
    d.Ambient, d.Diffuse = 1.0, 0.0  # unlit: colours exactly as in the colour map
    ColorBy(d, ("POINTS", array))
    lut = GetColorTransferFunction(array)
    lut.ApplyPreset(preset, True)
    if array != "schlieren":
        lut.InvertTransferFunction()  # light = low, dark = high
    lut.AutomaticRescaleRangeMode = "Never"
    lut.RescaleTransferFunction(*rng)
    if log:
        lut.UseLogScale = 1
    GetOpacityTransferFunction(array).RescaleTransferFunction(*rng)

    if array != "schlieren":
        bar = GetScalarBar(lut, v)
        bar.Title = ""
        bar.ComponentTitle = ""
        bar.Orientation = "Vertical"
        bar.WindowLocation = "Any Location"
        bar.Position = [0.89, 0.12]
        bar.ScalarBarLength = 0.62
        bar.ScalarBarThickness = 14
        bar.LabelFontSize = 15
        bar.RangeLabelFormat = "%.3g"
        d.SetScalarBarVisibility(v, True)
    else:
        d.SetScalarBarVisibility(v, False)

    label = Text(registrationName=f"title {array}", Text=title)
    td = Show(label, v)
    td.WindowLocation = "Upper Left Corner"
    td.FontSize = 20

    # Leave a margin on the right for the colour bar and at the top for the title
    v.CameraPosition = [xc + 0.18, yc + 0.06, 5.0]
    v.CameraFocalPoint = [xc + 0.18, yc + 0.06, 0.0]
    v.CameraViewUp = [0.0, 1.0, 0.0]
    v.CameraParallelScale = 0.75

    # Domain outline, so still air (light) stays distinguishable from the background
    edges = Show(reader, v)
    edges.Representation = "Feature Edges"
    ColorBy(edges, None)
    edges.AmbientColor = edges.DiffuseColor = [0.32, 0.32, 0.31]
    edges.LineWidth = 1.5

# Time stamp in the first panel
tm = AnnotateTimeFilter(registrationName="time", Input=reader)
tm.Format = "t = {time:.2f} ms"
tm.Scale = 1000.0
tmd = Show(tm, views[0])
tmd.WindowLocation = "Upper Right Corner"
tmd.FontSize = 20

for v in views[1:]:
    AddCameraLink(views[0], v, f"link{id(v)}")

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
times = list(reader.TimestepValues)
scene.AnimationTime = times[-1] if times else 0.0
SetActiveView(views[0])

# ---------------------------------------------------------------------------
# Batch output
# ---------------------------------------------------------------------------

if args.still:
    t, png = float(args.still[0]), args.still[1]
    scene.AnimationTime = min(times, key=lambda s: abs(s - t))
    SaveScreenshot(png, layout, SaveAllViews=1, ImageResolution=[1920, 1080],
                   FontScaling="Do not scale fonts")
    print(f"wrote {png} at t = {scene.AnimationTime*1e3:.3f} ms")

def save_frames(folder):
    folder.mkdir(parents=True, exist_ok=True)
    SaveAnimation(str(folder / "frame.png"), layout, SaveAllViews=1,
                  ImageResolution=[1920, 1080], FontScaling="Do not scale fonts",
                  FrameWindow=[0, len(times) - 1])
    print(f"wrote {len(times)} frames to {folder}")


if args.frames:
    save_frames(args.frames)

if args.movie:
    ffmpeg = shutil.which("ffmpeg") or "/opt/homebrew/bin/ffmpeg"
    with tempfile.TemporaryDirectory() as tmp:
        save_frames(Path(tmp))
        pattern = str(Path(tmp) / "frame.%04d.png")
        args.movie.parent.mkdir(parents=True, exist_ok=True)
        # yuv420p + even dimensions: plays in PowerPoint, Keynote and browsers
        subprocess.run([ffmpeg, "-y", "-loglevel", "error", "-framerate", str(args.fps),
                        "-i", pattern, "-c:v", "libx264", "-pix_fmt", "yuv420p",
                        "-crf", "18", "-movflags", "+faststart", str(args.movie)], check=True)
    print(f"wrote {args.movie}")
