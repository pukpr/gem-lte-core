# The stellarator part: wavenumber > 0 (2026-10-22)

Companion to `qbo_switching_plasma.md` (the QBO as the reversed-field-pinch modification of the tokamak). Code:
`enso_stellarator.py`, `qbo_low_stellarator.py`; figures `enso_stellarator.png`, `qbo_low_stellarator.png`.

## Mathematical form

On a common manifold F(t) with base winding k:

  ζ(λ, t) = Σ_h Σ_n B_hn · Y_n(λ) · sin(2π·h·k·F(t) + φ_hn)

- **Tokamak / reversed-field pinch (axisymmetric):** n = 0 only; the zonal-mean flow (the QBO) carries all windings.
- **Stellarator:** n > 0 modes Y_n(λ) appear, and each winding h may select its own modes (a matrix B_hn that is not
  rank 1). Y_n is cos(nπx/L) in a bounded basin (standing) and e^{inλ} around the full equatorial ring (helical).

## 1. ENSO: winding-selective modes do not help

Eight equatorial Kaplan quads, 130°E–90°W, fitted 1950–2005, scored 2006–2022, common manifold (base 0.2017), windings
1–9× and the slow winding.

| Structure | Parameters | Held-out mean r | Held-out pattern r | Training mean r |
|---|---|---|---|---|
| one shared pattern for all windings (rank 1) | 28 | **+0.31** | **+0.15** | +0.42 |
| free profile per winding | 160 | +0.27 | +0.11 | +0.50 |
| basin modes n = 0, 1, 2 mixed per winding | 60 | +0.24 | +0.02 | +0.47 |
| one basin mode per winding | 20 | +0.17 | +0.10 | +0.38 |

Within the Pacific basin, every winding shares one standing pattern. The winding-specific shapes (the eastern vs
central-Pacific split in `enso_spatial_mode.py`) fit the training years but do not predict. In the basin, ENSO acts
as one mode.

## 2. Equatorial winds: wavenumber > 0 emerges below 70 hPa, and it follows the ENSO manifold

ERA5 monthly zonal wind, 5°S–5°N at each longitude, 1979–2025, monthly climatology removed per longitude.

**Share of variance by zonal wavenumber:**

| hPa | 1–10 | 20–50 | 70 | 100 | 125–150 | 200 |
|---|---|---|---|---|---|---|
| n = 0 | 0.90–0.99 | 0.99–1.00 | 0.92 | 0.23 | 0.17 | 0.20 |
| n = 1 | 0.01–0.07 | ≤ 0.01 | 0.04 | 0.34 | 0.46–0.48 | 0.40 |
| n = 2 | ≤ 0.02 | 0.00 | 0.02 | 0.24 | 0.20–0.21 | 0.22 |

The transition is abrupt: from 70 to 100 hPa the zonal-mean share falls from 0.92 to 0.23.

**Which manifold drives each wavenumber** (regression on nine windings of each manifold, fitted 1979–2010; held-out
correlation 2011–2022; in brackets the training fit's percentile among 300 phase-randomized surrogates):

| Level | n = 0, QBO windings | n = 0, ENSO windings | n = 1–2, QBO windings | n = 1–2, ENSO windings |
|---|---|---|---|---|
| 5 hPa | **+0.46 (100)** | −0.01 (36) | +0.08 | −0.05 |
| 30 hPa | **+0.67 (100)** | −0.01 (70) | −0.06 | −0.09 |
| 50 hPa | **+0.34 (83)** | +0.26 (2) | −0.02 | +0.03 |
| 100 hPa | −0.21 (41) | +0.10 (76) | −0.02 | **+0.19 to +0.22 (68–98)** |
| 125–150 hPa | −0.02 to −0.09 | +0.03 to +0.05 | 0.00 to +0.02 | **+0.24 to +0.26 (79–94)** |
| 200 hPa | 0.00 | +0.09 | −0.01 | **+0.25 to +0.28 (80–87)** |

- **Above 70 hPa:** the flow is axisymmetric (n = 0), and its time history follows the QBO manifold, out of sample
  (+0.67 at 30 hPa). The ENSO windings do not explain it.
- **Below 70 hPa:** wavenumbers 1–2 carry most of the variance, and they follow the ENSO manifold out of sample
  (+0.19 to +0.28 at every level from 100 to 200 hPa); the QBO windings give about zero there.
- So the vertical stack separates cleanly: a tokamak / reversed-field-pinch layer (n = 0, QBO manifold) over a
  stellarator layer (n > 0, ENSO manifold), with the boundary between 70 and 100 hPa, at the tropical tropopause.

## Caveats

- The n > 0 result says that the longitudinally varying winds near the tropopause are ENSO-driven, which is already
  well known (the Walker circulation). What is new is only the framing and that the tidally forced ENSO windings
  predict it on held-out years. Held-out correlations of 0.2–0.3 are modest, and the levels are not independent.
- Surrogate percentiles for n > 0 (68–98) are suggestive, not strong.
- Monthly data cannot see the fast waves (Kelvin and mixed Rossby-gravity waves, days to weeks) that also carry n > 0
  in the lower stratosphere; daily ERA5 at several longitudes would.
- In the Pacific SST itself, the stellarator-style per-winding selectivity was not supported (section 1); the n > 0
  structure there is one fixed standing pattern, not winding-selective.
