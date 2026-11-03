#!/usr/bin/env python3
"""
Post-processing and validation for the inviscid compression-ramp case.

Reads the latest time of the OpenFOAM case (ASCII fields), compares it with
exact shock-expansion theory (theory.py) and writes figures plus a validation
table to results/.

    python3 scripts/postprocess.py                     # main case
    python3 scripts/postprocess.py --study meshStudy   # + mesh convergence study and GCI
    python3 scripts/postprocess.py --sweep angleSweep  # + ramp-angle sweep
"""

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon
from matplotlib.ticker import MaxNLocator
from matplotlib.tri import LinearTriInterpolator, Triangulation

import theory

ROOT = Path(__file__).resolve().parent.parent
DEPTH = 0.01  # span of the 2D mesh [m], zF - zB in system/blockMeshDict

# Palette: one-hue blue ramp for magnitude, fixed categorical order for series
BLUE_RAMP = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
ORDINAL = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]  # coarse -> fine
ORANGE = "#eb6834"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
CMAP = LinearSegmentedColormap.from_list("blue", BLUE_RAMP)

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.grid": False,
    "grid.color": GRID, "grid.linewidth": 0.6, "legend.frameon": False,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 200, "savefig.bbox": "tight",
})


# ---------------------------------------------------------------------------
# OpenFOAM I/O
# ---------------------------------------------------------------------------

def read_params(case):
    txt = re.sub(r"//.*", "", (case / "system" / "caseParameters").read_text())
    return {k: float(v) for k, v in re.findall(r"^\s*(\w+)\s+([-+0-9.eE]+)\s*;", txt, re.M)}


def latest_time(case):
    times = []
    for d in case.iterdir():
        try:
            times.append((float(d.name), d.name))
        except ValueError:
            pass
    return max(times)[1]


def read_field(path):
    """internalField of an ASCII volScalarField / volVectorField (plain or .gz)."""
    path = Path(path)
    if not path.exists() and path.with_name(path.name + ".gz").exists():
        path = path.with_name(path.name + ".gz")
    txt = gzip.open(path, "rt").read() if path.suffix == ".gz" else path.read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s+(\d+)\s*\(", txt)
    n, start = int(m.group(2)), m.end()
    body = txt[start:txt.index("\n)", start)]
    if m.group(1) == "vector":
        return np.array(body.replace("(", " ").replace(")", " ").split(), float).reshape(n, 3)
    return np.array(body.split(), float)


def read_raw(path):
    """Surface sample in raw format: x y z value."""
    d = np.loadtxt(path, comments="#")
    return d[np.argsort(d[:, 0])]


class Case:
    def __init__(self, path):
        self.path = Path(path)
        P = self.prm = read_params(self.path)
        self.time = latest_time(self.path)
        d = self.path / self.time

        self.x, self.y = read_field(d / "Cx"), read_field(d / "Cy")
        self.p, self.T, self.rho = (read_field(d / f) for f in ("p", "T", "rho"))
        self.U = read_field(d / "U")
        self.Ma = read_field(d / "Ma")
        grad = d / "grad(rho)"
        self.gradRho = read_field(grad) if grad.exists() else None

        self.g, self.R = P["gamma"], P["R"]
        self.pInf, self.TInf, self.M = P["pInf"], P["TInf"], P["Mach"]
        self.rhoInf = self.pInf / (self.R * self.TInf)
        self.UInf = self.M * np.sqrt(self.g * self.R * self.TInf)
        self.L, self.H = P["rampLength"], P["height"]
        self.xUp, self.xDown, self.dx = P["xUp"], P["xDown"], P["cellSize"]
        self.th = np.radians(P["rampAngle"])
        self.h = self.L * np.tan(self.th)
        self.sol = theory.ramp_solution(self.M, P["rampAngle"], self.g)

        self.tri = Triangulation(self.x, self.y)
        t = self.tri.triangles
        self.tri.set_mask(self.y[t].mean(1) < self.wall_y(self.x[t].mean(1)))

    def wall_y(self, x):
        return np.clip(x, 0.0, self.L) * np.tan(self.th)

    def body(self, **kw):
        xo, yb = self.L + self.xDown, -0.1
        pts = [(-self.xUp, 0), (0, 0), (self.L, self.h), (xo, self.h), (xo, yb), (-self.xUp, yb)]
        return Polygon(pts, closed=True, fc="#d9d8d4", ec=INK2, lw=0.8, zorder=5, **kw)

    def wall_samples(self, field="p"):
        return read_raw(self.path / "postProcessing" / "wallSamples" / self.time / f"{field}_walls.raw")

    def forces(self):
        rows = []
        for f in sorted((self.path / "postProcessing" / "rampForce").glob("*/force.dat"),
                        key=lambda f: float(f.parent.name)):
            rows.append(np.loadtxt(f, comments="#", usecols=(0, 1, 2)))
        return np.vstack(rows)

    # -- measurements -------------------------------------------------------

    def shock_points(self, n=15):
        """Shock location along horizontal lines, at the mid-pressure crossing."""
        s = self.sol
        pmid = 0.5 * (1.0 + s.s2.p) * self.pInf
        xi, yi = s.shock_fan_intersection(self.L)
        interp = LinearTriInterpolator(self.tri, self.p)
        lead = np.tan(s.fan_leading_angle())
        pts = []
        for y in np.linspace(0.12, 0.85, n) * min(yi, self.H):
            # stop before the ramp surface (y < h) or the expansion fan (y > h)
            xs = y / np.tan(self.th) if y < self.h else self.L + (y - self.h) / lead
            x = np.linspace(-self.xUp + self.dx, xs - 2 * self.dx, 3000)
            p = np.asarray(interp(x, np.full_like(x, y)).filled(np.nan))
            i = np.flatnonzero(p > pmid)
            if i.size and i[0] > 0:
                k = i[0]
                pts.append((np.interp(pmid, [p[k - 1], p[k]], [x[k - 1], x[k]]), y))
        return np.array(pts)

    def shock_angle(self):
        pts = self.shock_points()
        slope, icpt = np.polyfit(pts[:, 0], pts[:, 1], 1)
        return np.degrees(np.arctan(slope)), -icpt / slope, pts

    def region_masks(self):
        s, L, x, y = self.sol, self.L, self.x, self.y
        # Region 2: core of the wedge between ramp and shock, upstream of the fan
        thick = x * (np.tan(s.beta) - np.tan(self.th))
        f2 = (y - self.wall_y(x)) / np.where(thick > 0, thick, np.nan)
        r2 = (x > 0.3 * L) & (x < 0.9 * L) & (f2 > 0.25) & (f2 < 0.75)
        # Region 3: below the fan's trailing characteristic, downstream of the shoulder
        xr = x - L
        f3 = (y - self.h) / np.where(xr > 0, xr * np.tan(s.mu3), np.nan)
        r3 = (xr > 0.2 * self.xDown) & (xr < 0.9 * self.xDown) & (f3 > 0.2) & (f3 < 0.7)
        return r2, r3

    def region_stats(self):
        out = {}
        angle = np.degrees(np.arctan2(self.U[:, 1], self.U[:, 0]))
        for name, m in zip(("2", "3"), self.region_masks()):
            out[name] = dict(
                M=self.Ma[m].mean(), p=self.p[m].mean() / self.pInf,
                T=self.T[m].mean() / self.TInf, rho=self.rho[m].mean() / self.rhoInf,
                angle=angle[m].mean(), n=int(m.sum()),
            )
        return out

    def ramp_force_pressure(self):
        """Mean ramp pressure p2/pInf implied by the converged x-force."""
        Fx = self.forces()[-1, 1]
        return 1.0 + Fx / (self.pInf * self.h * DEPTH)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def overlay_theory(ax, c, color=ORANGE, label=True):
    s = c.sol
    xi, yi = s.shock_fan_intersection(c.L)
    ax.plot([0, xi], [0, yi], ls="--", lw=1.2, color=color, zorder=6,
            label="exact theory" if label else None)
    ax.plot([c.L, xi], [c.h, yi], ls="--", lw=1.0, color=color, zorder=6)
    r = 0.9 * c.xDown
    ax.plot([c.L, c.L + r], [c.h, c.h + r * np.tan(s.mu3)], ls="--", lw=1.0, color=color, zorder=6)


def domain_axes(ax, c):
    ax.set_xlim(-c.xUp, c.L + c.xDown)
    ax.set_ylim(-0.03, c.H)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.add_patch(c.body())


def fig_field(c, values, label, fname, out, levels=41, lines=True):
    fig, ax = plt.subplots(figsize=(10, 5.6))
    lv = np.linspace(values.min(), values.max(), levels)
    cf = ax.tricontourf(c.tri, values, levels=lv, cmap=CMAP)
    if lines:
        ax.tricontour(c.tri, values, levels=lv[::2], colors="white", linewidths=0.35, alpha=0.6)
    overlay_theory(ax, c)
    domain_axes(ax, c)
    s = c.sol
    for txt, (x, y) in {"1": (-0.15, 0.6), "2": (0.4, 0.17), "3": (1.6, 0.22)}.items():
        ax.text(x, y, txt, fontsize=13, weight="bold", color=INK, ha="center", va="center", zorder=7)
    cb = fig.colorbar(cf, ax=ax, shrink=0.8, pad=0.02, ticks=MaxNLocator(8))
    cb.set_label(label)
    cb.outline.set_visible(False)
    ax.legend(loc="upper left")
    ax.set_title(f"{label}:  M$_\\infty$ = {c.M:g},  ramp {np.degrees(c.th):g}°,  "
                 f"β$_{{theory}}$ = {np.degrees(s.beta):.2f}°", loc="left", color=INK)
    fig.savefig(out / fname)
    plt.close(fig)


def fig_schlieren(c, out):
    if c.gradRho is None:
        print("  (skipping schlieren: run  postProcess -func 'grad(rho)' -latestTime)")
        return
    g = np.linalg.norm(c.gradRho[:, :2], axis=1)
    S = np.exp(-15.0 * g / np.percentile(g, 99.9))
    fig, ax = plt.subplots(figsize=(10, 5.6))
    ax.tripcolor(c.tri, S, cmap="gray", shading="gouraud", vmin=0, vmax=1)
    domain_axes(ax, c)
    ax.set_title("Numerical schlieren  exp(−k|∇ρ|/|∇ρ|max)", loc="left", color=INK)
    fig.savefig(out / "schlieren.png")
    plt.close(fig)


def fig_mesh(c, out):
    """Mesh near the ramp, rebuilt from the block structure (straight-edged blocks)."""
    fig, ax = plt.subplots(figsize=(10, 4.2))
    nY = round(c.H / c.dx)
    xs = np.concatenate([np.linspace(-c.xUp, 0, round(c.xUp / c.dx) + 1),
                         np.linspace(0, c.L, round(c.L / c.dx) + 1),
                         np.linspace(c.L, c.L + c.xDown, round(c.xDown / c.dx) + 1)])
    xs = np.unique(xs)
    yw = c.wall_y(xs)
    eta = np.linspace(0, 1, nY + 1)
    Y = yw[None, :] + eta[:, None] * (c.H - yw[None, :])
    kw = dict(color=INK2, lw=0.25)
    ax.plot(np.tile(xs, (2, 1)), np.vstack([yw, np.full_like(yw, c.H)]), **kw)
    ax.plot(xs, Y.T, **kw)
    ax.add_patch(c.body())
    ax.set_xlim(-0.1, c.L + 0.3)
    ax.set_ylim(-0.02, 0.35)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    nCells = nY * (round(c.xUp / c.dx) + round(c.L / c.dx) + round(c.xDown / c.dx))
    ax.set_title(f"Structured mesh near the ramp: Δ = {c.dx * 1e3:g} mm, {nCells:,} cells",
                 loc="left", color=INK)
    fig.savefig(out / "mesh.png")
    plt.close(fig)


def fig_wall_pressure(cases, out):
    c0 = cases[-1]
    s = c0.sol
    fig, ax = plt.subplots(figsize=(9, 4.5))
    xo = c0.L + c0.xDown
    ax.plot([-c0.xUp, 0, 0, c0.L, c0.L, xo], [1, 1, s.s2.p, s.s2.p, s.s3.p, s.s3.p],
            color=INK, lw=1.6, ls="--", label="exact shock-expansion theory", zorder=5)
    for c, col in zip(cases, ORDINAL[-len(cases):]):
        d = c.wall_samples("p")
        ax.plot(d[:, 0], d[:, 3] / c.pInf, color=col, lw=1.4,
                label=f"CFD, Δ = {c.dx * 1e3:g} mm")
    ax.axvspan(0, c0.L, color=GRID, alpha=0.5, lw=0)
    ax.text(c0.L / 2, 0.3, "ramp", ha="center", color=INK2)
    ax.text(c0.L + 0.04, s.s2.p, f"p₂/p∞ = {s.s2.p:.3f}  (oblique shock)", va="center", color=INK2)
    ax.text(c0.L + 0.6, s.s3.p - 0.3, f"p₃/p∞ = {s.s3.p:.3f}  (Prandtl–Meyer)", va="center", color=INK2)
    ax.set_xlim(-c0.xUp, xo)
    ax.set_ylim(0, None)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("wall pressure  p/p∞")
    ax.grid(axis="y")
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.5))
    ax.set_title("Wall pressure distribution", loc="left", color=INK)
    fig.savefig(out / "wall_pressure.png")
    plt.close(fig)


def fig_shock_angle(c, beta, x0, pts, out):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.tricontourf(c.tri, c.p / c.pInf, levels=41, cmap=CMAP)
    xi, yi = c.sol.shock_fan_intersection(c.L)
    xx = np.linspace(0, xi, 10)
    ax.plot(xx, xx * np.tan(c.sol.beta), color=INK, ls="--", lw=1.2,
            label=f"theory  β = {np.degrees(c.sol.beta):.2f}°")
    ax.plot(xx, (xx - x0) * np.tan(np.radians(beta)), color=ORANGE, lw=1.2,
            label=f"CFD fit  β = {beta:.2f}°")
    ax.plot(pts[:, 0], pts[:, 1], "o", ms=5, mfc=ORANGE, mec="white", mew=1, zorder=7,
            label="CFD shock location (p = (p₁+p₂)/2)")
    ax.add_patch(c.body())
    ax.set_xlim(-0.1, xi + 0.1)
    ax.set_ylim(-0.02, yi + 0.05)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.legend(loc="upper left")
    ax.set_title("Oblique shock angle", loc="left", color=INK)
    fig.savefig(out / "shock_angle.png")
    plt.close(fig)


def theta_beta_curve(M, g, n=600):
    """theta [deg] and beta [deg] along the theta-beta-M curve, with a weak-branch mask."""
    b = np.linspace(np.arcsin(1 / M) + 1e-6, np.pi / 2, n)
    return np.degrees(theory.theta_from_beta(M, b, g)), np.degrees(b), b <= theory.beta_max(M, g)


def fig_theta_beta_m(c, beta_cfd, out, sweep=None):
    """theta-beta-M chart; the case Mach number is drawn in black with the CFD points on it.

    sweep: optional list of (theta_deg, beta_cfd_deg) from the ramp-angle sweep.
    """
    fig, ax = plt.subplots(figsize=(7.5, 5))
    others = [M for M in (1.5, 2.0, 3.0, 5.0) if abs(M - c.M) > 0.05]
    blues = ["#9ec5f4", "#6da7ec", "#3987e5", "#256abf"]
    for M, col, lw, lab in [(M, b, 1.3, None) for M, b in zip(others, blues)] + \
                           [(c.M, INK, 1.8, "this case")]:
        th, b, weak = theta_beta_curve(M, c.g)
        ax.plot(th[weak], b[weak], color=col, lw=lw, zorder=4 if lab else 3)
        ax.plot(th[~weak], b[~weak], color=col, lw=0.9, ls=":", zorder=4 if lab else 3)
        k = np.argmax(th)
        ax.text(th[k], b[k] + 2.0, f"M = {M:g}", color=INK if lab else INK2, ha="center",
                va="bottom", fontsize=9, weight="bold" if lab else "normal", zorder=8,
                bbox=dict(fc="white", ec="none", pad=0.6))
    if sweep:
        ths, bs = np.array(sweep).T
        ax.plot(ths, bs, "x", ms=8, mew=2, color=ORANGE, zorder=7,
                label=f"CFD, ramp-angle sweep (M = {c.M:g})")
    else:
        ax.plot(np.degrees(c.th), beta_cfd, "x", ms=9, mew=2, color=ORANGE, label="CFD", zorder=7)
    ax.plot(np.degrees(c.th), np.degrees(c.sol.beta), "o", ms=10, mfc="none", mec=INK, mew=1.2,
            label=f"main case: θ = {np.degrees(c.th):g}°, β = {np.degrees(c.sol.beta):.2f}°",
            zorder=6)
    ax.set_xlim(0, 50)
    ax.set_ylim(0, 90)
    ax.set_xlabel("deflection angle θ [deg]")
    ax.set_ylabel("shock angle β [deg]")
    ax.grid()
    ax.legend(loc="lower right", title="solid: weak  ·  dotted: strong", title_fontsize=9)
    ax.set_title("θ–β–M relation", loc="left", color=INK)
    fig.savefig(out / "theta_beta_M.png")
    plt.close(fig)


def fig_convergence(c, out):
    f = c.forces()
    tft = (c.xUp + c.L + c.xDown) / c.UInf
    Fx_th = (c.sol.s2.p - 1.0) * c.pInf * c.h * DEPTH
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(8, 5.5), sharex=True)
    a1.plot(f[:, 0] / tft, f[:, 1] / Fx_th, color=ORDINAL[1], lw=1.4)
    a1.axhline(1.0, color=INK, ls="--", lw=1)
    a1.set_ylabel("Fₓ / Fₓ,theory")
    a1.grid(axis="y")
    a1.set_title("Convergence: pressure force on the ramp", loc="left", color=INK)
    t = f[:, 0] / tft
    Fm = f[t > t[-1] - 1.0, 1].mean()  # mean over the last flow-through
    a2.semilogy(t, np.abs(f[:, 1] / Fm - 1.0), color=ORDINAL[1], lw=1.0)
    a2.set_ylabel("|Fₓ / F̄ₓ − 1|")
    a2.text(0.99, 0.95, "F̄ₓ: mean over the last flow-through", transform=a2.transAxes,
            ha="right", va="top", color=INK2, fontsize=9)
    a2.set_xlabel("t / t_flow-through   (t_flow-through = domain length / U∞)")
    a2.grid(axis="y")
    fig.savefig(out / "convergence.png")
    plt.close(fig)


def fig_angle_sweep(cases, metrics, out):
    """CFD vs exact theory over ramp angle: shock angle and ramp pressure."""
    c0 = cases[0]
    th, b, weak = theta_beta_curve(c0.M, c0.g, 4000)
    th_max = th.max()
    tt = np.array([np.degrees(c.th) for c in cases])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4))

    a1.plot(th[weak], b[weak], color=INK, lw=1.6, label="exact, weak shock")
    a1.plot(th[~weak], b[~weak], color=INK, lw=0.9, ls=":", label="exact, strong shock")
    a1.plot(tt, [m["beta"] for m in metrics], "x", ms=8, mew=2, color=ORANGE, zorder=6,
            label=f"CFD, Δ = {c0.dx * 1e3:g} mm")
    a1.axvline(th_max, color=INK2, lw=0.8, ls="--")
    a1.text(th_max - 0.4, 22, f"θ_max = {th_max:.1f}°\n(shock detaches)", ha="right",
            color=INK2, fontsize=9)
    a1.set_xlim(0, 32)
    a1.set_ylim(20, 92)
    a1.set_xlabel("ramp angle θ [deg]")
    a1.set_ylabel("shock angle β [deg]")
    a1.grid()
    a1.legend(loc="center left", bbox_to_anchor=(0.0, 0.62), fontsize=9)
    a1.set_title(f"Shock angle, M∞ = {c0.M:g}", loc="left", color=INK)

    # p2 from the oblique shock alone; p3 only while M2 > 1 (close to theta_max the
    # weak-shock M2 is subsonic and there is no Prandtl-Meyer fan at the shoulder)
    tw = np.linspace(0.01, th_max - 1e-3, 300)
    p2, p3 = [], []
    for t in tw:
        beta = theory.oblique_shock_beta(c0.M, np.radians(t), c0.g)
        p2.append(theory.normal_shock(c0.M * np.sin(beta), c0.g)["p"])
        try:
            with np.errstate(invalid="ignore"):
                p3.append(theory.ramp_solution(c0.M, t, c0.g).s3.p)
        except ValueError:
            p3.append(np.nan)
    a2.plot(tw, p2, color=INK, lw=1.6, label="exact p₂/p∞ (oblique shock)")
    a2.plot(tw, p3, color=INK, lw=1.0, ls="--", label="exact p₃/p∞ (after the expansion fan)")
    a2.plot(tt, [c.ramp_force_pressure() for c in cases], "x", ms=8, mew=2, color=ORANGE,
            zorder=6, label="CFD p₂, from the ramp force")
    a2.plot(tt, [m["stats"]["3"]["p"] for m in metrics], "o", ms=6, mfc="white", mec=ORANGE,
            mew=1.5, zorder=6, label="CFD p₃, region-3 average")
    a2.axvline(th_max, color=INK2, lw=0.8, ls="--")
    a2.set_xlim(0, 32)
    a2.set_ylim(0.8, None)
    a2.set_xlabel("ramp angle θ [deg]")
    a2.set_ylabel("p/p∞")
    a2.grid()
    a2.legend(loc="upper left", fontsize=9)
    a2.set_title("Ramp and shoulder pressure", loc="left", color=INK)
    fig.tight_layout()
    fig.savefig(out / "angle_sweep.png")
    plt.close(fig)


def sweep_table(cases, metrics):
    pct = lambda a, b: 100.0 * (a - b) / b
    lines = ["", f"## Ramp-angle sweep (M∞ = {cases[0].M:g}, Δ = {cases[0].dx*1e3:g} mm)", "",
             "| θ [deg] | β exact | β CFD | β error | p₂/p∞ exact | p₂/p∞ CFD (force) | p₂ error "
             "| M₂ exact | M₂ CFD | M₂ error | p₃/p∞ exact | p₃/p∞ CFD | p₃ error |",
             "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for c, m in zip(cases, metrics):
        s, r2, r3 = c.sol, m["stats"]["2"], m["stats"]["3"]
        pF = c.ramp_force_pressure()
        lines.append(
            f"| {np.degrees(c.th):g} | {np.degrees(s.beta):.3f} | {m['beta']:.3f} | "
            f"{m['beta'] - np.degrees(s.beta):+.3f}° | {s.s2.p:.4f} | {pF:.4f} | "
            f"{pct(pF, s.s2.p):+.2f} % | {s.s2.M:.4f} | {r2['M']:.4f} | {pct(r2['M'], s.s2.M):+.2f} % | "
            f"{s.s3.p:.4f} | {r3['p']:.4f} | {pct(r3['p'], s.s3.p):+.2f} % |")
    return "\n".join(lines)


def mean_cell_size(c):
    """Representative cell size h = sqrt(A/N) of the 2D mesh [m]."""
    area = (c.xUp + c.L + c.xDown) * c.H - 0.5 * c.L * c.h - c.xDown * c.h
    return np.sqrt(area / c.x.size)


def gci(h, phi, Fs=1.25):
    """Three-grid Richardson extrapolation and Grid Convergence Index.

    Celik et al., J. Fluids Eng. 130 (2008) 078001. h and phi are ordered
    fine -> coarse. Returns the observed order p, the extrapolated value and
    the fine-grid GCI (a relative error band on the fine-grid value).
    """
    (h1, h2, h3), (f1, f2, f3) = h, phi
    r21, r32 = h2 / h1, h3 / h2
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0.0 or e32 == 0.0:
        return None
    s = np.sign(e32 / e21)
    p = abs(np.log(abs(e32 / e21))) / np.log(r21)
    for _ in range(200):  # fixed-point iteration for unequal refinement ratios
        q = np.log((r21**p - s) / (r32**p - s))
        p_new = abs(np.log(abs(e32 / e21)) + q) / np.log(r21)
        if abs(p_new - p) < 1e-12:
            break
        p = p_new
    ext = (r21**p * f1 - f2) / (r21**p - 1.0)
    return dict(p=p, ext=ext, gci=Fs * abs(e21 / f1) / (r21**p - 1.0), oscillatory=s < 0)


def gci_table(cases, metrics):
    """GCI for the finest three meshes; the exact value is shown alongside to check it."""
    cs, ms = cases[-3:][::-1], metrics[-3:][::-1]  # fine -> coarse
    h = [mean_cell_size(c) for c in cs]
    s = cs[0].sol
    quantities = [
        ("Shock angle β [deg]", [m["beta"] for m in ms], np.degrees(s.beta)),
        ("p₂/p∞ from the ramp force", [c.ramp_force_pressure() for c in cs], s.s2.p),
        ("p₂/p∞, region-2 core", [m["stats"]["2"]["p"] for m in ms], s.s2.p),
        ("M₂, region-2 core", [m["stats"]["2"]["M"] for m in ms], s.s2.M),
        ("M₃, region-3 average", [m["stats"]["3"]["M"] for m in ms], s.s3.M),
        ("p₃/p∞, region-3 average", [m["stats"]["3"]["p"] for m in ms], s.s3.p),
    ]
    lines = ["", "## Grid Convergence Index (Celik et al. 2008)", "",
             f"Meshes: h = {h[2]*1e3:.2f}, {h[1]*1e3:.2f}, {h[0]*1e3:.2f} mm "
             f"(h = √(A/N)); refinement ratios r₃₂ = {h[2]/h[1]:.3f}, r₂₁ = {h[1]/h[0]:.3f}.", "",
             "| Quantity | coarse φ₃ | medium φ₂ | fine φ₁ | observed order p | extrapolated φ_ext "
             "| GCI_fine | exact | fine-mesh error vs exact | exact inside GCI band? |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|"]
    rows = []
    for name, phi, exact in quantities:
        g = gci(h, phi)
        err = 100.0 * abs(phi[0] - exact) / exact
        if g is None:
            lines.append(f"| {name} | {phi[2]:.4f} | {phi[1]:.4f} | {phi[0]:.4f} | – | – | – | "
                         f"{exact:.4f} | {err:.3f} % | – |")
            continue
        osc = " (oscillatory)" if g["oscillatory"] else ""
        lines.append(f"| {name} | {phi[2]:.4f} | {phi[1]:.4f} | {phi[0]:.4f} | {g['p']:.2f}{osc} | "
                     f"{g['ext']:.4f} | {100*g['gci']:.3f} % | {exact:.4f} | {err:.3f} % | "
                     f"{'yes' if err <= 100*g['gci'] else 'no'} |")
        rows.append((name, g, err))
    lines += ["", "GCI_fine = 1.25 |φ₁ − φ₂| / (|φ₁| (r₂₁ᵖ − 1)). Where the exact value lies "
              "inside the GCI band, the error estimate is conservative."]
    return "\n".join(lines), rows


def fig_mesh_study(cases, metrics, out):
    """Error vs. exact theory against cell size, log-log, with order-1 and order-2 slopes."""
    dx = np.array([c.dx for c in cases]) * 1e3
    s = cases[0].sol
    series = [
        ("wall pressure on the ramp (L1 norm)", [m["L1_wall"] for m in metrics], "#1baf7a"),
        ("p₂ from the integrated ramp force",
         [abs(c.ramp_force_pressure() - s.s2.p) / s.s2.p for c in cases], ORANGE),
        ("M₃ behind the expansion fan",
         [abs(m["stats"]["3"]["M"] - s.s3.M) / s.s3.M for m in metrics], ORDINAL[3]),
        ("shock angle β", [m["beta_rel"] for m in metrics], ORDINAL[1]),
    ]
    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    for lab, err, col in series:
        ax.loglog(dx, err, "o-", color=col, lw=1.4, ms=6, label=lab)
    # Reference slopes through the coarse-mesh point of the first and last series
    for order, ref, ls in [(1, series[0][1][0], ":"), (2, series[-1][1][0], "--")]:
        ax.loglog(dx, ref * (dx / dx[0]) ** order, color=INK2, lw=0.8, ls=ls,
                  label=f"slope {order} (order {order})")
    ax.set_xticks(dx)
    ax.set_xticklabels([f"{d:g}" for d in dx])
    ax.minorticks_off()
    ax.set_xlabel("cell size Δ [mm]")
    ax.set_ylabel("relative error vs. exact theory")
    ax.grid(alpha=0.6)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=9)
    ax.set_title(f"Mesh convergence, M∞ = {cases[0].M:g}", loc="left", color=INK)
    fig.savefig(out / "mesh_convergence.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Validation table
# ---------------------------------------------------------------------------

def metrics_for(c):
    beta, x0, pts = c.shock_angle()
    st = c.region_stats()
    d = c.wall_samples("p")
    on_ramp = (d[:, 0] > 0) & (d[:, 0] < c.L)
    L1 = np.mean(np.abs(d[on_ramp, 3] / c.pInf - c.sol.s2.p)) / c.sol.s2.p
    return dict(dx=c.dx, beta=beta, beta_err=abs(beta - np.degrees(c.sol.beta)),
                beta_rel=abs(np.radians(beta) - c.sol.beta) / c.sol.beta,
                p2_err=abs(st["2"]["p"] - c.sol.s2.p) / c.sol.s2.p, L1_wall=L1,
                peak=d[on_ramp, 3].max() / c.pInf, stats=st, x0=x0, pts=pts)


def validation_table(c, m):
    s = c.sol
    pct = lambda a, b: 100.0 * (a - b) / b
    r2, r3 = m["stats"]["2"], m["stats"]["3"]
    pF = c.ramp_force_pressure()
    rows = [
        ("Shock angle β [deg]", np.degrees(s.beta), m["beta"]),
        ("Region 2  M₂", s.s2.M, r2["M"]),
        ("Region 2  p₂/p∞", s.s2.p, r2["p"]),
        ("Region 2  T₂/T∞", s.s2.T, r2["T"]),
        ("Region 2  ρ₂/ρ∞", s.s2.rho, r2["rho"]),
        ("Region 2  flow angle [deg]", np.degrees(c.th), r2["angle"]),
        ("Ramp pressure from force, p₂/p∞", s.s2.p, pF),
        ("Region 3  M₃", s.s3.M, r3["M"]),
        ("Region 3  p₃/p∞", s.s3.p, r3["p"]),
        ("Region 3  T₃/T∞", s.s3.T, r3["T"]),
        ("Region 3  flow angle [deg]", 0.0, r3["angle"]),
    ]
    lines = [f"# Validation: M∞ = {c.M:g}, ramp angle {np.degrees(c.th):g}°, Δ = {c.dx*1e3:g} mm, "
             f"t = {c.time} s", "",
             "| Quantity | Exact theory | CFD | Error |", "|---|---:|---:|---:|"]
    for name, th, cf in rows:
        err = f"{cf - th:+.3f}°" if "deg" in name else f"{pct(cf, th):+.2f} %"
        lines.append(f"| {name} | {th:.4f} | {cf:.4f} | {err} |")
    lines += ["", f"Region-2 average over {r2['n']} cells, region-3 over {r3['n']} cells.",
              f"Shock line fit virtual origin x₀ = {m['x0']*1e3:+.2f} mm (corner at 0).",
              f"Peak wall pressure on the ramp (corner overshoot): {m['peak']:.3f} p∞."]
    return "\n".join(lines)


# ---------------------------------------------------------------------------

def study_cases(d, key):
    """Cases in directory d, sorted by key(parameters)."""
    dirs = [s for s in d.iterdir() if (s / "system" / "caseParameters").exists()]
    return [Case(s) for s in sorted(dirs, key=lambda s: key(read_params(s)))]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", default=ROOT / "rampCase", type=Path)
    ap.add_argument("--study", default=None, type=Path, help="directory of mesh-study cases")
    ap.add_argument("--sweep", default=None, type=Path, help="directory of ramp-angle sweep cases")
    ap.add_argument("--out", default=ROOT / "results", type=Path)
    a = ap.parse_args()

    figs = a.out / "figures"
    figs.mkdir(parents=True, exist_ok=True)

    print(f"Reading {a.case} ...")
    c = Case(a.case)
    m = metrics_for(c)
    report = [validation_table(c, m), "", "## Exact theory", "", "```", theory.report(c.sol), "```"]
    print(report[0])

    print("Writing figures ...")
    fig_field(c, c.Ma, "Mach number", "mach.png", figs)
    fig_field(c, c.p / c.pInf, "p / p∞", "pressure.png", figs)
    fig_field(c, c.T / c.TInf, "T / T∞", "temperature.png", figs, lines=False)
    fig_schlieren(c, figs)
    fig_mesh(c, figs)
    fig_shock_angle(c, m["beta"], m["x0"], m["pts"], figs)
    fig_convergence(c, figs)

    cases = [c]
    if a.study:
        cases = study_cases(a.study, lambda P: -P["cellSize"])  # coarse -> fine
        ms = [metrics_for(s) for s in cases]
        fig_mesh_study(cases, ms, figs)
        lines = ["", "## Mesh convergence (error vs. exact theory)", "",
                 "| Δ [mm] | cells | β [deg] | β error [deg] | p₂ error (core) | wall L1 error | peak p/p∞ |",
                 "|---:|---:|---:|---:|---:|---:|---:|"]
        for s, mm in zip(cases, ms):
            lines.append(f"| {s.dx*1e3:g} | {s.x.size:,} | {mm['beta']:.3f} | {mm['beta_err']:.3f} | "
                         f"{100*mm['p2_err']:.2f} % | {100*mm['L1_wall']:.2f} % | {mm['peak']:.3f} |")
        report.append("\n".join(lines))
        if len(cases) >= 3:
            report.append(gci_table(cases, ms)[0])
        print("\n".join(report[-2:]))
    fig_wall_pressure(cases, figs)

    sweep = None
    if a.sweep:
        sc = study_cases(a.sweep, lambda P: P["rampAngle"])
        sm = [metrics_for(s) for s in sc]
        fig_angle_sweep(sc, sm, figs)
        report.append(sweep_table(sc, sm))
        print(report[-1])
        sweep = [(np.degrees(s.th), mm["beta"]) for s, mm in zip(sc, sm)]
    fig_theta_beta_m(c, m["beta"], figs, sweep)

    (a.out / "validation.md").write_text("\n".join(report) + "\n")
    print(f"Done: {figs}")


if __name__ == "__main__":
    main()
