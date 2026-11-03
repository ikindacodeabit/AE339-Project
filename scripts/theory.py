"""
Exact inviscid theory for supersonic flow over a compression ramp with a
shoulder (shock-expansion theory), for a calorically perfect gas.

    Region 1 : freestream
    Region 2 : behind the attached oblique shock from the ramp corner
    Region 3 : behind the Prandtl-Meyer expansion fan from the ramp shoulder

Run directly to print the solution for the default case:
    python3 theory.py [Mach] [rampAngleDeg] [gamma]
"""

import sys
from dataclasses import dataclass

import numpy as np


# ---------------------------------------------------------------------------
# Oblique shock
# ---------------------------------------------------------------------------

def theta_from_beta(M1, beta, g=1.4):
    """Flow deflection angle theta [rad] for shock angle beta [rad] (theta-beta-M relation)."""
    num = M1**2 * np.sin(beta)**2 - 1.0
    den = M1**2 * (g + np.cos(2.0 * beta)) + 2.0
    return np.arctan(2.0 / np.tan(beta) * num / den)


def beta_max(M1, g=1.4):
    """Shock angle [rad] at maximum deflection; separates weak and strong branches."""
    b = np.linspace(np.arcsin(1.0 / M1) + 1e-9, np.pi / 2 - 1e-9, 20001)
    return b[np.argmax(theta_from_beta(M1, b, g))]


def _bisect(f, a, b, tol=1e-13, maxit=200):
    fa = f(a)
    for _ in range(maxit):
        m = 0.5 * (a + b)
        fm = f(m)
        if fa * fm <= 0.0:
            b = m
        else:
            a, fa = m, fm
        if b - a < tol:
            break
    return 0.5 * (a + b)


def oblique_shock_beta(M1, theta, g=1.4, strong=False):
    """Shock angle beta [rad] for deflection theta [rad]. Weak solution by default."""
    mu = np.arcsin(1.0 / M1)
    bm = beta_max(M1, g)
    if theta > theta_from_beta(M1, bm, g):
        raise ValueError("theta exceeds theta_max: the shock detaches (no attached solution)")
    f = lambda b: theta_from_beta(M1, b, g) - theta
    return _bisect(f, bm, np.pi / 2) if strong else _bisect(f, mu + 1e-12, bm)


def normal_shock(Mn1, g=1.4):
    """Ratios across a normal shock for upstream normal Mach number Mn1."""
    p = 1.0 + 2.0 * g / (g + 1.0) * (Mn1**2 - 1.0)
    rho = (g + 1.0) * Mn1**2 / ((g - 1.0) * Mn1**2 + 2.0)
    T = p / rho
    Mn2 = np.sqrt((1.0 + 0.5 * (g - 1.0) * Mn1**2) / (g * Mn1**2 - 0.5 * (g - 1.0)))
    p0 = rho**(g / (g - 1.0)) * (1.0 / p)**(1.0 / (g - 1.0))
    return dict(p=p, rho=rho, T=T, Mn2=Mn2, p0=p0)


# ---------------------------------------------------------------------------
# Prandtl-Meyer expansion
# ---------------------------------------------------------------------------

def prandtl_meyer(M, g=1.4):
    """Prandtl-Meyer function nu(M) [rad]."""
    a = np.sqrt((g + 1.0) / (g - 1.0))
    return a * np.arctan(np.sqrt((M**2 - 1.0) / a**2)) - np.arctan(np.sqrt(M**2 - 1.0))


def inverse_prandtl_meyer(nu, g=1.4):
    """Mach number for a given Prandtl-Meyer angle nu [rad]."""
    nu_max = 0.5 * np.pi * (np.sqrt((g + 1.0) / (g - 1.0)) - 1.0)
    if not 0.0 <= nu < nu_max:
        raise ValueError("nu outside [0, nu_max)")
    return _bisect(lambda M: prandtl_meyer(M, g) - nu, 1.0, 1e3)


def isentropic_p(M, g=1.4):
    """Static-to-total pressure ratio p/p0."""
    return (1.0 + 0.5 * (g - 1.0) * M**2) ** (-g / (g - 1.0))


def isentropic_T(M, g=1.4):
    """Static-to-total temperature ratio T/T0."""
    return 1.0 / (1.0 + 0.5 * (g - 1.0) * M**2)


# ---------------------------------------------------------------------------
# Ramp solution
# ---------------------------------------------------------------------------

@dataclass
class State:
    M: float
    p: float     # p/p_inf
    T: float     # T/T_inf
    rho: float   # rho/rho_inf


@dataclass
class RampSolution:
    M1: float
    theta: float           # ramp angle [rad]
    gamma: float
    beta: float            # oblique shock angle [rad]
    s1: State
    s2: State
    s3: State
    mu2: float             # Mach angle in region 2 [rad]
    mu3: float             # Mach angle in region 3 [rad]
    nu2: float
    nu3: float
    p02_p01: float         # total-pressure ratio across the shock

    def fan_leading_angle(self):
        """Angle [rad] of the fan's leading characteristic to the horizontal."""
        return self.theta + self.mu2

    def fan_trailing_angle(self):
        """Angle [rad] of the fan's trailing characteristic to the horizontal."""
        return self.mu3

    def shock_fan_intersection(self, L):
        """Point where the corner shock meets the fan's leading characteristic.

        Corner at (0, 0), shoulder at (L, L tan(theta)). Above this point the
        expansion fan starts to weaken and curve the shock.
        """
        tb, tl, h = np.tan(self.beta), np.tan(self.fan_leading_angle()), L * np.tan(self.theta)
        x = (h - L * tl) / (tb - tl)
        return x, x * tb


def ramp_solution(M1, theta_deg, g=1.4):
    th = np.radians(theta_deg)
    beta = oblique_shock_beta(M1, th, g)

    ns = normal_shock(M1 * np.sin(beta), g)
    M2 = ns["Mn2"] / np.sin(beta - th)
    s2 = State(M=M2, p=ns["p"], T=ns["T"], rho=ns["rho"])

    # Isentropic turn back through theta at the shoulder
    nu2 = prandtl_meyer(M2, g)
    nu3 = nu2 + th
    M3 = inverse_prandtl_meyer(nu3, g)
    p3 = s2.p * isentropic_p(M3, g) / isentropic_p(M2, g)
    T3 = s2.T * isentropic_T(M3, g) / isentropic_T(M2, g)
    s3 = State(M=M3, p=p3, T=T3, rho=p3 / T3)

    return RampSolution(
        M1=M1, theta=th, gamma=g, beta=beta,
        s1=State(M=M1, p=1.0, T=1.0, rho=1.0), s2=s2, s3=s3,
        mu2=np.arcsin(1.0 / M2), mu3=np.arcsin(1.0 / M3),
        nu2=nu2, nu3=nu3, p02_p01=ns["p0"],
    )


# ---------------------------------------------------------------------------
# Start-up: a uniform stream entering still gas (transientCase)
# ---------------------------------------------------------------------------

@dataclass
class StartingShocks:
    Ms: float          # shock Mach number relative to the gas ahead of it
    W_primary: float   # lab speed of the shock running into the still gas [/a_inf]
    W_secondary: float # lab speed of the upstream-facing shock (swept downstream) [/a_inf]
    u_star: float      # gas velocity between the shocks [/a_inf]
    p: float           # p/p_inf between the shocks
    T: float           # T/T_inf between the shocks


def starting_shocks(M1, g=1.4):
    """1D Riemann problem: gas at Mach M1 colliding with the same gas at rest.

    Both initial states have the same p and T, so the solution is symmetric
    about the frame moving at U/2: two shocks of equal strength, a stagnant
    central state in that frame, and no visible contact surface. Valid along
    the upstream wall until the waves reflected from the ramp arrive.
    """
    u_star = 0.5 * M1
    k = 0.5 * (g + 1.0) * u_star          # Ms - 1/Ms, from u_p/a = 2/(g+1) (Ms - 1/Ms)
    Ms = 0.5 * (k + np.sqrt(k**2 + 4.0))
    ns = normal_shock(Ms, g)
    return StartingShocks(Ms=Ms, W_primary=Ms, W_secondary=M1 - Ms,
                          u_star=u_star, p=ns["p"], T=ns["T"])


def report(sol):
    d = np.degrees
    lines = [
        f"M1 = {sol.M1:g}, theta = {d(sol.theta):g} deg, gamma = {sol.gamma:g}",
        f"Oblique shock angle beta      = {d(sol.beta):.3f} deg",
        f"Total pressure ratio p02/p01  = {sol.p02_p01:.4f}",
        "",
        f"{'Region':<26}{'M':>8}{'p/p1':>9}{'T/T1':>9}{'rho/rho1':>10}",
    ]
    for name, s in [("1 freestream", sol.s1), ("2 behind oblique shock", sol.s2),
                    ("3 behind expansion fan", sol.s3)]:
        lines.append(f"{name:<26}{s.M:8.4f}{s.p:9.4f}{s.T:9.4f}{s.rho:10.4f}")
    lines += [
        "",
        f"nu2 = {d(sol.nu2):.3f} deg, nu3 = {d(sol.nu3):.3f} deg",
        f"Fan leading characteristic  = {d(sol.fan_leading_angle()):.3f} deg to horizontal",
        f"Fan trailing characteristic = {d(sol.fan_trailing_angle()):.3f} deg to horizontal",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    M = float(sys.argv[1]) if len(sys.argv) > 1 else 3.0
    th = float(sys.argv[2]) if len(sys.argv) > 2 else 15.0
    g = float(sys.argv[3]) if len(sys.argv) > 3 else 1.4
    print(report(ramp_solution(M, th, g)))
