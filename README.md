# AE339 Project: inviscid supersonic flow over a compression ramp

Inviscid (Euler) **Mach 2.4** flow of air (the group's Mach number) over a 15°
compression ramp followed by a flat shoulder, solved with OpenFOAM
`rhoCentralFoam`. The flow has three features. Exact theory predicts the first
two:

1. an **attached oblique shock** from the ramp corner (θ–β–M relation),
2. a centred **Prandtl–Meyer expansion fan** from the shoulder,
3. the **shock–expansion interaction** downstream, where the fan weakens and
   curves the shock. This part has no simple closed-form solution.

What the project shows:

- **Validation:** the CFD is compared with exact shock-expansion theory, for the
  main 15° case and for a sweep of ramp angles from 5° to 25° at Mach 2.4.
- **Grid convergence:** three meshes (20, 10 and 5 mm), with the error against
  the exact solution and a Grid Convergence Index (GCI) estimate.

```
                         top (supersonic far field)
     +-----------------------------------------------------------+
     |                                                           |
  M  |                                                           |  outlet
 --> |                          shoulder                         |
     |                       ___________________________________ |
     |                  ___/  (Prandtl-Meyer expansion fan)      |
     |_____________ ___/  ramp, 15 deg                           |
     upstream     corner  (oblique shock)
     |<- 0.25 m ->|<- 0.5 m ->|<------------ 1.5 m ------------>|
```

## Layout

| Path | Contents |
|---|---|
| `rampCase/system/caseParameters` | **All inputs in one file**: Mach, p∞, T∞, γ, R, geometry, cell size |
| `rampCase/0.orig/` | Initial and boundary conditions (p, T, U) |
| `rampCase/constant/thermophysicalProperties` | Calorically perfect air, μ = 0 (inviscid) |
| `rampCase/system/blockMeshDict` | Parametric 3-block structured mesh |
| `rampCase/system/{controlDict,fvSchemes,fvSolution}` | Solver, numerics, monitoring |
| `rampCase/Allrun`, `Allclean` | Mesh → solve → reconstruct → post-process |
| `scripts/theory.py` | Exact oblique-shock + Prandtl–Meyer solution |
| `scripts/postprocess.py` | Figures and validation table |
| `scripts/make_mesh_study.sh` | Builds the mesh-convergence cases (`meshStudy/`) |
| `scripts/make_angle_sweep.sh` | Builds the ramp-angle validation cases (`angleSweep/`) |
| `transientCase/` | Start-up animation case: same flow, starts from still air, 240 frames |
| `scripts/transient_figures.py` | Storyboard, x–t diagram and probe histories of the start-up |
| `scripts/view_transient.py` | ParaView set-up (GUI) and MP4 rendering (batch) for the start-up |
| `run.sh` | One-command driver |
| `results/` | `validation.md` and `figures/` |

## Running

Requires the OpenFOAM-v2606 macOS app (`openfoam2606` on the PATH) and Python 3
with NumPy and Matplotlib.

```bash
./run.sh            # main case (5 mm mesh, 112,500 cells, 4 cores), then post-process
./run.sh study      # 20 / 10 / 5 mm mesh-convergence study + GCI, then post-process
./run.sh sweep      # ramp angles 5-25 deg at 10 mm (validation on the theta-beta-M curve)
./run.sh post       # post-processing only (includes the study and sweep if present)
./run.sh clean      # delete all solution data
python3 scripts/theory.py 2.4 15   # exact solution for any M, theta
```

To change the flow or geometry, edit `rampCase/system/caseParameters` only.
The mesh, boundary conditions, gas properties and analytical comparison all
read from it. If θ exceeds θ_max for the chosen Mach number, the shock
detaches and `theory.py` reports that there is no attached solution.

To view the fields in ParaView, open `rampCase/rampCase.foam`.

### Notes for this machine

- **Paths with spaces.** OpenFOAM refuses to run in a directory whose path
  contains spaces (`College stoof/Autumn 2026`). `run.sh` runs everything
  through the symlink `~/ae339_project` → this folder. For parallel runs,
  `Allrun` passes `mpirun -x PWD` so that the MPI ranks keep that
  space-free path.
- **MPI library.** `~/.OpenFOAM/prefs.sh` sets `WM_MPLIB=SYSTEMOPENMPI`.
  The v2606 app is built with `OPENMPI`, so `run.sh` overrides it with
  `-DWM_MPLIB=OPENMPI`. Without the override, the parallel solver loads a
  serial-only stub library and exits.

## Start-up animation (`transientCase/`)

The same flow and geometry, run as a **time-accurate start-up** for viewing in
ParaView. At t = 0 the domain is filled with still air at (p∞, T∞), and the
Mach 2.4 stream enters through the inlet, as when a supersonic wind tunnel
starts. The steady shock pattern of `rampCase/` is fully established on the wall
by 14.6 ms.

| | |
|---|---|
| Mesh | 2.5 mm (900 × 500 = 450,000 cells), twice as fine as `rampCase/` |
| Output | 240 frames, one every 75 µs, from 0 to 18 ms (compressed ASCII, 1.6 GB) |
| Fields | `Ma`, `p`, `T`, `U`, `rho`, `magGradRho` (\|∇ρ\| for numerical schlieren) |
| Extras | Pressure probes at 5 points, ramp force, wall p/T/Ma at every frame |
| Run time | About 2 h 20 min on 4 cores (`./run.sh transient`) |

```bash
./run.sh transient   # run it (about 2.3 h), then make the start-up figures
./run.sh view        # open it in ParaView with a 2x2 layout (Mach, p, T, schlieren); press Play
./run.sh movie       # render results/rampStartup.mp4 (1920x1080, H.264)
```

To view it by hand, open `transientCase/transientCase.foam` in ParaView, tick
the fields, click **Apply**, and press **Play**. The `.foam` file is an empty
marker that tells ParaView which case folder to read.

### Sequence of events

See `results/figures/startup_storyboard.png`. The times below were measured
from the wall-pressure x–t data, the probes and the ramp force.

1. **Starting shock (0–0.37 ms).** The incoming Mach 2.4 stream acts like a
   piston on the still air. A shock runs ahead into the still air at 678 m/s.
   A second, upstream-facing shock is swept downstream at only 155 m/s. Between
   them the gas moves at 417 m/s with p = 4.28 p∞ and T = 1.65 T∞. Both initial
   states have the same p and T, so this is an exact, symmetric 1D Riemann
   problem with no visible contact surface. The probe on the upstream wall
   (x = −0.10 m) reproduces the exact pulse: 4.280 p∞ from 0.221 to 0.956 ms,
   against 4.280 p∞ from 0.221 to 0.965 ms exact (`startup_probes.png`).
2. **Reflection off the ramp (0.37–1.1 ms).** The starting shock reaches the
   corner at about 0.37 ms and the shoulder at 1.1 ms. The reflection is a Mach
   reflection: an incident shock, a reflected shock and a Mach stem meet at a
   triple point. The ramp force peaks at 3.2× its steady value at 1.4 ms.
3. **Diffraction at the shoulder (1.1–2.5 ms).** The shock turns round the
   convex corner, and the flow rolls up a vortex at the sharp edge. The Euler
   equations have no viscosity, so the vortex comes from the sharp corner and
   the scheme's numerical dissipation. It is carried downstream (x ≈ 1.0 m at
   5 ms, 1.85 m at 7 ms).
4. **The second shock climbs the ramp (1.6–7 ms).** The second shock reaches
   the corner at about 1.6 ms and moves slowly up the ramp with a lambda foot.
   Behind it is the real Mach 2.4 stream, so the attached oblique shock grows
   from the corner. The starting shock leaves through the outlet at 3.3 ms. The
   ramp pressure is steady (within 3 %) from 7.0 ms.
5. **Sweep-out (7–14.6 ms).** The second shock passes the shoulder, crosses the
   expansion fan and runs out along the shock–fan interaction. The wall is steady
   up to x = 1.0 m from 8.9 ms, to 1.5 m from 11.3 ms, and everywhere from
   14.6 ms. The region-3 probe is steady from 9.6 ms. The 1D estimate for the
   second shock to cross the 2.25 m domain at 155 m/s is 14.5 ms.
6. **Steady state (from 14.6 ms).** What remains is the steady pattern of
   `rampCase/`: the oblique shock, the expansion fan, and the shock curving where
   they interact. The final probe pressures match theory: 2.403 and 2.404 p∞
   behind the shock (exact 2.405) and 1.011 p∞ behind the fan (exact 1.011). The
   final frame agrees with the steady 5 mm solution with a median pressure
   difference of 0.01 %. The probe in the interaction zone, (1.60, 0.90), still
   varies by ±0.6 % at 18 ms.

`startup_xt_wall.png` shows the same story as an x–t diagram of the wall
pressure: each wave is a line, and its slope gives the wave speed.

## Numerical method

| Item | Choice |
|---|---|
| Solver | `rhoCentralFoam`: density-based, explicit, finite-volume |
| Convective flux | Kurganov–Tadmor central-upwind (Kurganov, Noelle & Petrova 2001) |
| Reconstruction | Second-order, van Leer TVD limiter on ρ, U, T |
| Time integration | Forward Euler, adaptive Δt at Courant number 0.2. Marched to steady state over 3.7 flow-through times (10 ms) |
| Gas | Calorically perfect air: γ = 1.4, R = 287.1 J/(kg K), μ = 0 |
| Mesh | Structured, uniform Δ = 5 mm, 450 × 250 cells |

| Boundary | p | T | U | Reason |
|---|---|---|---|---|
| inlet | fixed | fixed | fixed (M = 2.4) | Supersonic inflow: all characteristics enter |
| outlet | zeroGradient | zeroGradient | zeroGradient | Supersonic outflow: all characteristics leave |
| top | zeroGradient | zeroGradient | zeroGradient | Flow is tangential or leaving |
| walls | zeroGradient | zeroGradient | slip | Inviscid wall: no normal velocity |
| front/back | empty | empty | empty | 2D |

Free stream: p∞ = 101 325 Pa, T∞ = 300 K, U∞ = 833.40 m/s (M = 2.4). One
flow-through time (2.25 m / U∞) is 2.70 ms.

## Results (5 mm mesh, t = 0.01 s = 3.7 flow-throughs)

Exact theory for M∞ = 2.4, θ = 15°: β = 38.214°, θ_max = 28.7°, so the shock is
attached and well below detachment.

| Quantity | Exact theory | CFD | Error |
|---|---:|---:|---:|
| Shock angle β | 38.214° | 38.224° | +0.010° |
| p₂/p∞ behind shock | 2.4048 | 2.4052 | +0.02 % |
| M₂ behind shock | 1.7927 | 1.7913 | −0.08 % |
| T₂/T∞ behind shock | 1.3100 | 1.3106 | +0.04 % |
| ρ₂/ρ∞ behind shock | 1.8357 | 1.8352 | −0.03 % |
| Flow angle behind shock | 15.000° | 15.002° | +0.002° |
| p₂/p∞ from integrated ramp force | 2.4048 | 2.4045 | −0.01 % |
| M₃ after expansion fan | 2.3495 | 2.3483 | −0.05 % |
| p₃/p∞ after expansion fan | 1.0113 | 1.0117 | +0.04 % |
| T₃/T∞ after expansion fan | 1.0228 | 1.0233 | +0.04 % |

### Mesh convergence

Three meshes, each with half the cell size of the previous one. The error is
measured against the exact solution.

| Δ | Cells | β error | p₂ error (region-2 core) | Wall-pressure L1 error | Peak corner p/p∞ | Solver time |
|---:|---:|---:|---:|---:|---:|---:|
| 20 mm | 7,119 | 0.121° | 0.55 % | 6.99 % | 3.00 | 10 s (1 core) |
| 10 mm | 28,125 | 0.030° | 0.00 % | 3.59 % | 2.99 | 87 s (1 core) |
| 5 mm | 112,500 | 0.012° | 0.02 % | 1.86 % | 2.98 | 8 min (4 cores) |

- The **wall-pressure L1 error halves each time Δ halves** (ratios 1.95 and
  1.93). That is first-order convergence, which is expected for any
  shock-capturing scheme: at a discontinuity, TVD limiters drop to first order.
- The **shock angle** converges at between first and second order (error ratios
  4.0 and 2.5).
- In the smooth **region-2 core**, p₂ is already within 0.02 % on the 10 mm
  mesh. Below that the error is set by the residual unsteadiness (≈ 2 × 10⁻⁴,
  see below), not by the mesh.

### Grid Convergence Index (GCI)

Richardson extrapolation and GCI by the three-grid procedure of Celik et al.,
*J. Fluids Eng.* 130 (2008) 078001, with h = √(A/N) (refinement ratio 2.0). Here
the exact answer is known, so it is possible to check whether the GCI error band
contains it:

| Quantity | φ (20 / 10 / 5 mm) | Observed order p | Extrapolated | GCI (5 mm) | Exact | Actual error (5 mm) | Exact inside band? |
|---|---|---:|---:|---:|---:|---:|:---:|
| p₂/p∞ from ramp force | 2.4022 / 2.4035 / 2.4045 | 0.40 | 2.4077 | 0.16 % | 2.4048 | 0.012 % | yes |
| M₂, region-2 core | 1.7875 / 1.7887 / 1.7913 | 1.12 | 1.7935 | 0.16 % | 1.7927 | 0.076 % | yes |
| M₃, region-3 average | 2.3407 / 2.3456 / 2.3479 | 1.10 | 2.3499 | 0.11 % | 2.3495 | 0.067 % | yes |
| Shock angle β | 38.334 / 38.244 / 38.226 | 2.38 | 38.222 | 0.014 % | 38.214 | 0.032 % | no |

For the integrated force and the regional averages, the observed order is
about 1 and the GCI band contains the exact value, so it is a safe error
estimate. For β it is not. The fitted shock angle keeps a small offset of
≈ 0.01° that does not shrink with the mesh (from the half-pressure criterion on
a smeared, slightly unsteady shock). Richardson extrapolation cannot see an
error that does not change with the mesh. This is why validation against an exact solution is needed
as well as a grid-convergence study. The full table, including two quantities
with oscillatory convergence at the noise floor, is in `results/validation.md`.

### Validation over ramp angle

The same case at θ = 5°, 10°, 15°, 20° and 25° (10 mm mesh, 75–90 s each in
serial), up to just below θ_max = 28.7°:

| θ | β exact | β CFD | p₂/p∞ exact | p₂/p∞ CFD (force) | M₂ exact | M₂ CFD | p₃/p∞ exact | p₃/p∞ CFD |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5° | 28.528° | 28.558° | 1.3661 | 1.3659 | 2.1986 | 2.1985 | 1.0002 | 1.0003 |
| 10° | 33.023° | 33.049° | 1.8292 | 1.8286 | 1.9994 | 1.9979 | 1.0025 | 1.0028 |
| 15° | 38.214° | 38.244° | 2.4048 | 2.4035 | 1.7927 | 1.7887 | 1.0113 | 1.0122 |
| 20° | 44.336° | 44.392° | 3.1155 | 3.1120 | 1.5689 | 1.5639 | 1.0333 | 1.0350 |
| 25° | 52.172° | 52.193° | 4.0257 | 4.0249 | 1.3061 | 1.3035 | 1.0772 | 1.0777 |

β is within 0.06° and p₂ within 0.11 % at every angle, including θ = 25°, where
M₂ = 1.31 and the flow behind the shock is close to sonic. Above θ_max the
shock would detach and stand ahead of the corner as a curved bow shock, so
there is no exact attached solution to compare with.

Full tables: `results/validation.md`. Figures in
`results/figures/`:

| Figure | Shows |
|---|---|
| `mach.png`, `pressure.png`, `temperature.png` | Fields with the exact shock and fan overlaid |
| `schlieren.png` | Numerical schlieren: shock, fan, curved shock after the interaction |
| `shock_angle.png` | Shock located along 15 lines and fitted, compared with θ–β–M |
| `theta_beta_M.png` | θ–β–M chart with the M = 2.4 curve and the CFD points of the angle sweep |
| `angle_sweep.png` | β and p₂, p₃ against ramp angle: CFD vs. exact, up to θ_max |
| `wall_pressure.png` | Wall pressure vs. shock-expansion theory, all meshes |
| `convergence.png` | Ramp force history: steady after ~0.6 flow-throughs |
| `mesh.png`, `mesh_convergence.png` | Mesh, and error vs. cell size with order-1 and order-2 slopes |

### Points to discuss

- **Shock–expansion interaction.** Above (1.22, 0.96) m, the fan's leading
  characteristic (48.9° to the horizontal) reaches the shock. The shock then
  weakens and curves toward the free-stream Mach angle (24.6°). Theory gives
  the intersection point; the curved part needs the method of characteristics
  or CFD. The domain (1.25 m high) is tall enough to show the start of the curved part.
- **Corner overshoot and undershoot.** The wall pressure spikes to about
  3.0 p∞ just after the ramp corner and dips to about 0.78 p∞ just after the
  shoulder. These are numerical artifacts: the exact solution is a step. The
  wall cells at a corner have to turn the flow 15° within one cell. The
  spike's width shrinks with the cell size (see the mesh study), but its
  height does not (3.00, 2.99, 2.98 p∞). That is typical behaviour of a
  captured discontinuity, as with the Gibbs phenomenon. The integrated force
  is still within 0.02 %.
- **Residual unsteadiness.** The ramp force reaches 10⁻³ of its final value
  after 0.57 flow-throughs and then levels off at about 2 × 10⁻⁴ (r.m.s.) instead of
  converging further. This is a small limit cycle from the captured shock
  moving across a grid that is not aligned with it. It also causes the faint
  striations behind the shock in the schlieren.
- **Time step.** The Courant number had to be 0.2. At 0.5, the run diverged
  at the outlet when the start-up shock left the domain.
