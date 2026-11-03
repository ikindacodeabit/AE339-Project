#!/usr/bin/env python3
"""
Figures for the start-up (transient) case: how the steady pattern forms.

    python3 scripts/transient_figures.py                  # transientCase -> results/figures
    python3 scripts/transient_figures.py --case <dir> --out <dir>

Writes
    startup_storyboard.png   numerical schlieren at nine key instants
    startup_xt_wall.png      x-t diagram of the wall pressure (shock paths)
    startup_probes.png       pressure histories at the probes + ramp force
"""

import argparse
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.patches import Polygon
from matplotlib.tri import Triangulation

import theory
from postprocess import (DEPTH, GRID, INK, INK2, ORANGE, ORDINAL, read_field,
                         read_params, read_raw)

ROOT = Path(__file__).resolve().parent.parent
ORANGE_RAMP = LinearSegmentedColormap.from_list(
    "orange", ["#fdeee6", "#f9c4a8", "#f39566", "#eb6834", "#c4481a", "#8f3110"])

# Instants for the storyboard and what happens there [s]. Measured from the
# M = 2.4 run (wall x-t data, probes, ramp force); see README "Sequence of events".
STORY = [
    (0.00030, "Starting shock (right) runs into still air;\nsecond shock follows the Mach 2.4 stream"),
    (0.00090, "Starting shock reflects off the ramp:\nMach reflection with a triple point"),
    (0.00150, "Shock diffracts round the shoulder and rolls\nup a vortex; second shock reaches the corner"),
    (0.00250, "Attached oblique shock grows from the corner,\nbehind the second shock"),
    (0.00400, "Second shock climbs the ramp with a lambda foot;\nvortex carried downstream"),
    (0.00600, "Second shock near the shoulder;\nstarting shock has left (3.3 ms)"),
    (0.00850, "Ramp steady (from 7.0 ms); second shock\nswept downstream past the fan"),
    (0.01200, "Second shock runs out along the\nshock-fan interaction"),
    (0.01800, "Steady state (wall steady from 14.6 ms):\noblique shock, fan and interaction"),
]


class Transient:
    def __init__(self, path):
        self.path = Path(path)
        P = self.prm = read_params(self.path)
        self.times = sorted((float(d.name), d.name) for d in self.path.iterdir()
                            if d.is_dir() and _is_number(d.name) and float(d.name) > 0)
        self.x, self.y = read_field(self.path / "0" / "Cx"), read_field(self.path / "0" / "Cy")
        self.pInf, self.TInf, self.M, self.g = P["pInf"], P["TInf"], P["Mach"], P["gamma"]
        self.a = np.sqrt(self.g * P["R"] * self.TInf)
        self.L, self.H, self.xUp, self.xDown = P["rampLength"], P["height"], P["xUp"], P["xDown"]
        self.th = np.radians(P["rampAngle"])
        self.h = self.L * np.tan(self.th)
        self.sol = theory.ramp_solution(self.M, P["rampAngle"], self.g)
        self.start = theory.starting_shocks(self.M, self.g)

        self.tri = Triangulation(self.x, self.y)
        t = self.tri.triangles
        self.tri.set_mask(self.y[t].mean(1) < self.wall_y(self.x[t].mean(1)))

    def wall_y(self, x):
        return np.clip(x, 0.0, self.L) * np.tan(self.th)

    def nearest(self, t):
        return min(self.times, key=lambda s: abs(s[0] - t))

    def field(self, name, time):
        return read_field(self.path / time / name)

    def body(self):
        xo, yb = self.L + self.xDown, -0.1
        pts = [(-self.xUp, 0), (0, 0), (self.L, self.h), (xo, self.h), (xo, yb), (-self.xUp, yb)]
        return Polygon(pts, closed=True, fc="#d9d8d4", ec=INK2, lw=0.6, zorder=5)


def _is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


# ---------------------------------------------------------------------------

def fig_storyboard(c, out):
    # One reference for all panels so that frames are comparable
    G_ref = np.percentile(c.field("magGradRho", c.times[-1][1]), 99.9)
    fig, axs = plt.subplots(3, 3, figsize=(13.2, 9.4))
    for k, (ax, (t, text)) in enumerate(zip(axs.flat, STORY)):
        tv, tn = c.nearest(t)
        S = np.exp(-15.0 * c.field("magGradRho", tn) / G_ref)
        ax.tripcolor(c.tri, S, cmap="gray", shading="gouraud", vmin=0, vmax=1, rasterized=True)
        ax.add_patch(c.body())
        ax.set_xlim(-c.xUp, c.L + c.xDown)
        ax.set_ylim(-0.03, c.H)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(True)
            sp.set_color(GRID)
        ax.set_title(f"({chr(97 + k)})  t = {tv * 1e3:.2f} ms", loc="left", fontsize=11, color=INK)
        ax.text(0.5, -0.04, text, transform=ax.transAxes, ha="center", va="top",
                fontsize=9, color=INK2)
    fig.suptitle(f"Start-up of Mach {c.M:g} flow over a {np.degrees(c.th):g}° ramp "
                 "(numerical schlieren; domain initially still air)",
                 x=0.01, ha="left", fontsize=13, color=INK)
    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.06, top=0.91, hspace=0.45, wspace=0.04)
    fig.savefig(out / "startup_storyboard.png")
    plt.close(fig)


def fig_xt(c, out):
    base = c.path / "postProcessing" / "wallSamples"
    rows, ts = [], []
    for tv, tn in c.times:
        f = base / tn / "p_walls.raw"
        if f.exists():
            d = read_raw(f)
            rows.append(d[:, 3] / c.pInf)
            ts.append(tv)
            xw = d[:, 0]
    P, ts = np.array(rows), np.array(ts) * 1e3

    fig, ax = plt.subplots(figsize=(9, 6.5))
    pc = ax.pcolormesh(xw, ts, P, cmap=ORANGE_RAMP, norm=LogNorm(0.5, 12), shading="nearest",
                       rasterized=True)
    cb = fig.colorbar(pc, ax=ax, pad=0.02, ticks=[0.5, 1, 2, 3, 5, 8, 12])
    cb.ax.set_yticklabels(["0.5", "1", "2", "3", "5", "8", "12"])
    cb.set_label("wall pressure  p/p∞ (log scale)")
    cb.outline.set_visible(False)

    # 1D Riemann solution, valid on the upstream wall until reflections return
    x = np.array([-c.xUp, 0.0])
    for W, lab, ls in [(c.start.W_primary, "starting shock", "--"),
                       (c.start.W_secondary, "second shock", ":")]:
        ax.plot(x, (x + c.xUp) / (W * c.a) * 1e3, color=INK, ls=ls, lw=1.4,
                label=f"{lab}: {W * c.a:.0f} m/s (exact 1D theory)")
    ax.legend(loc="upper left", fontsize=8.5, frameon=True, facecolor="white", edgecolor=GRID)
    for xv, lab in [(0.0, "corner"), (c.L, "shoulder")]:
        ax.axvline(xv, color=INK2, lw=0.8, ls=":")
        ax.text(xv, 1.005, lab, color=INK2, fontsize=9, ha="center", va="bottom",
                transform=ax.get_xaxis_transform())
    ax.set_ylim(0, ts[-1])
    ax.set_xlabel("x along the wall [m]")
    ax.set_ylabel("t [ms]")
    ax.set_title("x–t diagram of the wall pressure: each line is a wave", loc="left", color=INK,
                 pad=16)
    fig.savefig(out / "startup_xt_wall.png")
    plt.close(fig)


def fig_probes(c, out):
    base = c.path / "postProcessing" / "probes"
    files = sorted(base.glob("*/p"), key=lambda f: float(f.parent.name))
    d = np.vstack([np.loadtxt(f, comments="#") for f in files])
    t, P = d[:, 0] * 1e3, d[:, 1:] / c.pInf
    s, st = c.sol, c.start
    probes = [
        ("upstream wall  (−0.10, 0.01)", 1.0),
        ("on the ramp  (0.25, 0.08)", s.s2.p),
        ("behind the oblique shock  (0.40, 0.20)", s.s2.p),
        ("behind the expansion fan  (1.20, 0.15)", s.s3.p),
        ("shock / fan interaction  (1.60, 0.90)", None),
    ]
    fig, axs = plt.subplots(len(probes) + 1, 1, figsize=(9, 11), sharex=True)
    for i, (ax, (name, p_ss)) in enumerate(zip(axs, probes)):
        ax.plot(t, P[:, i], color=ORDINAL[1], lw=1.2)
        if p_ss is not None:
            ax.axhline(p_ss, color=INK, ls="--", lw=1)
            ax.text(t[-1], p_ss, f"steady theory {p_ss:.2f} ", ha="right", va="bottom",
                    fontsize=8.5, color=INK)
        ax.set_ylabel("p/p∞")
        ax.set_title(name, loc="left", fontsize=10, color=INK)
        ax.grid(axis="y")
    # Exact 1D start-up pulse at the upstream probe
    x0 = -0.10 + c.xUp
    t1, t2 = x0 / (st.W_primary * c.a) * 1e3, x0 / (st.W_secondary * c.a) * 1e3
    axs[0].plot([0, t1, t1, t2, t2], [1, 1, st.p, st.p, 1], color=ORANGE, lw=1.2, ls="--")
    axs[0].text(t2 + 0.1, st.p * 0.9, f"1D Riemann solution: {st.p:.2f} p∞", fontsize=8.5,
                color=ORANGE, va="top")

    f = np.vstack([np.loadtxt(g, comments="#", usecols=(0, 1)) for g in
                   sorted((c.path / "postProcessing" / "rampForce").glob("*/force.dat"),
                          key=lambda g: float(g.parent.name))])
    Fx_th = (s.s2.p - 1.0) * c.pInf * c.h * DEPTH
    axs[-1].plot(f[:, 0] * 1e3, f[:, 1] / Fx_th, color=ORDINAL[2], lw=1.2)
    axs[-1].axhline(1.0, color=INK, ls="--", lw=1)
    axs[-1].set_ylabel("Fₓ / Fₓ,steady")
    axs[-1].set_title("pressure force on the ramp", loc="left", fontsize=10, color=INK)
    axs[-1].grid(axis="y")
    axs[-1].set_xlabel("t [ms]")
    fig.suptitle("Start-up histories (dashed: exact theory)", x=0.01, ha="left",
                 fontsize=13, color=INK)
    fig.tight_layout()
    fig.savefig(out / "startup_probes.png")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", default=ROOT / "transientCase", type=Path)
    ap.add_argument("--out", default=ROOT / "results" / "figures", type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    c = Transient(a.case)
    s = c.start
    print(f"{len(c.times)} frames, t = {c.times[0][0]*1e3:.2f} ... {c.times[-1][0]*1e3:.2f} ms")
    print(f"1D start-up theory: starting shock {s.W_primary*c.a:.1f} m/s, second shock "
          f"{s.W_secondary*c.a:.1f} m/s, p = {s.p:.3f} p_inf, T = {s.T:.3f} T_inf")
    fig_xt(c, a.out)
    fig_probes(c, a.out)
    fig_storyboard(c, a.out)
    print(f"Done: {a.out}")


if __name__ == "__main__":
    main()
