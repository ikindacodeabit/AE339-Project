# Validation: M∞ = 2.4, ramp angle 15°, Δ = 5 mm, t = 0.01 s

| Quantity | Exact theory | CFD | Error |
|---|---:|---:|---:|
| Shock angle β [deg] | 38.2137 | 38.2240 | +0.010° |
| Region 2  M₂ | 1.7927 | 1.7913 | -0.08 % |
| Region 2  p₂/p∞ | 2.4048 | 2.4052 | +0.02 % |
| Region 2  T₂/T∞ | 1.3100 | 1.3106 | +0.04 % |
| Region 2  ρ₂/ρ∞ | 1.8357 | 1.8352 | -0.03 % |
| Region 2  flow angle [deg] | 15.0000 | 15.0021 | +0.002° |
| Ramp pressure from force, p₂/p∞ | 2.4048 | 2.4045 | -0.01 % |
| Region 3  M₃ | 2.3495 | 2.3483 | -0.05 % |
| Region 3  p₃/p∞ | 1.0113 | 1.0117 | +0.04 % |
| Region 3  T₃/T∞ | 1.0228 | 1.0233 | +0.04 % |
| Region 3  flow angle [deg] | 0.0000 | 0.0006 | +0.001° |

Region-2 average over 1005 cells, region-3 over 9131 cells.
Shock line fit virtual origin x₀ = -0.24 mm (corner at 0).
Peak wall pressure on the ramp (corner overshoot): 2.982 p∞.

## Exact theory

```
M1 = 2.4, theta = 15 deg, gamma = 1.4
Oblique shock angle beta      = 38.214 deg
Total pressure ratio p02/p01  = 0.9346

Region                           M     p/p1     T/T1  rho/rho1
1 freestream                2.4000   1.0000   1.0000    1.0000
2 behind oblique shock      1.7927   2.4048   1.3100    1.8357
3 behind expansion fan      2.3495   1.0113   1.0228    0.9888

nu2 = 20.513 deg, nu3 = 35.513 deg
Fan leading characteristic  = 48.906 deg to horizontal
Fan trailing characteristic = 25.190 deg to horizontal
```

## Mesh convergence (error vs. exact theory)

| Δ [mm] | cells | β [deg] | β error [deg] | p₂ error (core) | wall L1 error | peak p/p∞ |
|---:|---:|---:|---:|---:|---:|---:|
| 20 | 7,119 | 38.334 | 0.121 | 0.55 % | 6.99 % | 3.000 |
| 10 | 28,125 | 38.244 | 0.030 | 0.00 % | 3.59 % | 2.987 |
| 5 | 112,500 | 38.226 | 0.012 | 0.02 % | 1.86 % | 2.982 |

## Grid Convergence Index (Celik et al. 2008)

Meshes: h = 19.03, 9.57, 4.79 mm (h = √(A/N)); refinement ratios r₃₂ = 1.988, r₂₁ = 2.000.

| Quantity | coarse φ₃ | medium φ₂ | fine φ₁ | observed order p | extrapolated φ_ext | GCI_fine | exact | fine-mesh error vs exact | exact inside GCI band? |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| Shock angle β [deg] | 38.3342 | 38.2438 | 38.2261 | 2.38 | 38.2218 | 0.014 % | 38.2137 | 0.032 % | no |
| p₂/p∞ from the ramp force | 2.4022 | 2.4035 | 2.4045 | 0.40 | 2.4077 | 0.164 % | 2.4048 | 0.012 % | yes |
| p₂/p∞, region-2 core | 2.4181 | 2.4047 | 2.4052 | 4.84 (oscillatory) | 2.4052 | 0.001 % | 2.4048 | 0.017 % | no |
| M₂, region-2 core | 1.7875 | 1.7887 | 1.7913 | 1.12 | 1.7935 | 0.155 % | 1.7927 | 0.076 % | yes |
| M₃, region-3 average | 2.3407 | 2.3456 | 2.3479 | 1.10 | 2.3499 | 0.107 % | 2.3495 | 0.067 % | yes |
| p₃/p∞, region-3 average | 1.0114 | 1.0122 | 1.0119 | 1.29 (oscillatory) | 1.0116 | 0.029 % | 1.0113 | 0.055 % | no |

GCI_fine = 1.25 |φ₁ − φ₂| / (|φ₁| (r₂₁ᵖ − 1)). Where the exact value lies inside the GCI band, the error estimate is conservative.

## Ramp-angle sweep (M∞ = 2.4, Δ = 10 mm)

| θ [deg] | β exact | β CFD | β error | p₂/p∞ exact | p₂/p∞ CFD (force) | p₂ error | M₂ exact | M₂ CFD | M₂ error | p₃/p∞ exact | p₃/p∞ CFD | p₃ error |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 28.528 | 28.558 | +0.031° | 1.3661 | 1.3659 | -0.01 % | 2.1986 | 2.1985 | -0.00 % | 1.0002 | 1.0003 | +0.01 % |
| 10 | 33.023 | 33.049 | +0.026° | 1.8292 | 1.8286 | -0.03 % | 1.9994 | 1.9979 | -0.07 % | 1.0025 | 1.0028 | +0.04 % |
| 15 | 38.214 | 38.244 | +0.030° | 2.4048 | 2.4035 | -0.05 % | 1.7927 | 1.7887 | -0.22 % | 1.0113 | 1.0122 | +0.09 % |
| 20 | 44.336 | 44.392 | +0.056° | 3.1155 | 3.1120 | -0.11 % | 1.5689 | 1.5639 | -0.31 % | 1.0333 | 1.0350 | +0.17 % |
| 25 | 52.172 | 52.193 | +0.022° | 4.0257 | 4.0249 | -0.02 % | 1.3061 | 1.3035 | -0.21 % | 1.0772 | 1.0777 | +0.05 % |
