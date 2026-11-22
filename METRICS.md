# METRIC Environment Variable — Algorithm Reference

The `METRIC` environment variable selects the fitness function used to evaluate how well the LTE model matches observational data. Higher values indicate better fits. Default is `CC`.

Set via: `METRIC=CC` (or any code below) when running `lt.exe`.

---

## Primary Metrics

### `CC` — Pearson Correlation Coefficient (default)

Standard Pearson's *r* between model output `X` and data `Y`, computed over the window `[CC_START, CC_END]` (excluding zero-valued points).

```
r = (N·Σxy − Σx·Σy) / √[(N·Σx² − (Σx)²) · (N·Σy² − (Σy)²)]
```

Returns `0.0` if the denominator is ≤ 0 (singular). Range: [−1, 1].

### `RMS` — Root-Mean-Sweet Error (inverted, bounded)

Inverse-transformed RMS error normalized by a reference amplitude:

```
RMS_score = 1 / (1 + √(Σ(xᵢ − yᵢ + offset)²) / Ref)
```

Returns values in (0, 1]. Higher scores mean lower RMS error. Controlled by `Ref` and `Offset` parameters.

### `ZC` — Zero-Crossing Agreement

Measures sign coherence between model and data, weighted by amplitude:

```
ZC = Σ(xᵢ · yᵢ) / Σ(|xᵢ| · |yᵢ|)
```

Fast alternative to CC for time series with many zero crossings. Range: [−1, 1]. Amplitude precision is not critical — only sign agreement matters.

### `DER` — Derivative Correlation

Pearson CC applied to the first differences (discrete derivatives) of both series:

```
DER_CC = CC(diff(X), diff(Y))
```

where `diff(S)[i] = S[i] − S[i−1]`. Rewards matching the *shape/rate of change* rather than absolute values.

---

## Spectral / Frequency-Domain Metrics

### `FT` — Fourier-Transform Spectrum Correlation

Computes a power spectrum for both model and data by projecting against sinusoidal bases at frequencies spaced from `FSTART` upward (linear step `FSTEP` or multiplicative `FMULT`). The resulting spectra are smoothed with `Window(..., 2)` (5-point rectangular), then compared via `DTW_Distance` (not CC). Effectively a shape-match in the frequency domain.

### `SEM` — Scaled Error Metric

A robust, scale-aware error metric based on Median Absolute Deviation (MAD):

1. Compute `MAD(Y)` = median of |yᵢ − median(Y)|
2. Compute `global_scale = max(MAD, α · mean(|Y|))` with α = 0.01
3. For each point: `score += (errorᵢ / max(α·|yᵢ|, global_scale))²`
4. `SEM = 1 / (1 + score / N)`

Returns values in (0, 1]. Robust to outliers via the MAD-based scale.

### `ME` — Minimum Entropy (via env var flag, not METRIC=ME directly)

Activated when `METRIC` is empty *and* the `ME` flag is detected internally. Minimizes the entropy of the power spectrum — i.e., concentrates spectral power into fewer frequencies. Two sub-modes:
- **MERMS=True**: minimizes `(ΣP − Σ√P)² / (Σ√P)²` — the "flatness" of the spectrum
- **MERMS=False** (default): splits data in half, computes power spectra of each, returns CC of the two smoothed spectra — rewards spectral consistency between first/second halves

---

## Distance / Dissimilarity Metrics (lower raw distance → higher score)

### `DTW` — Dynamic Time Warping

Sakoe-Chiba band-constrained DTW distance between model and data. The returned score is normalized:

```
DTW_score = (max_distance − distance(X, Y)) / max_distance
```

where `max_distance` is the DTW distance of `−Y` to `Y` (a worst-case reference). Allows temporal misalignment within a fixed window. Range: theoretically (−∞, 1], practically [0, 1].

### `EMD` — Earth Mover's Distance / Wasserstein-1

Sorts both series and computes mean absolute difference of sorted values:

```
EMD_raw = Σ|sorted(X)ᵢ − sorted(Y)ᵢ| / N
EMD_score = 1 / (1 + EMD_raw)
```

With `Derivative=True` (not exposed as a separate METRIC code), operates on first differences instead. Compares *distributions* rather than point-by-point alignment. Range: (0, 1].

### `CID` — Complexity-Invariant Distance

Ratio of total squared increments (complexity) between the two series:

```
CE(X) = Σ(xᵢ − xᵢ₊₁)²
CID = √(min(CE(X), CE(Y)) / max(CE(X), CE(Y)))
```

Rewards matching the *complexity* (total variation) of the two series. Range: (0, 1]. 1.0 means identical complexity.

### `CTW` — (Declared but uses CID internally)

Set in the metric dispatch table but maps to the same `CID` function. Present for naming compatibility.

---

## Composite / Hybrid Metrics

### `DC` — DTW + CC Composite

*(Declared but implementation delegates to existing DTW logic.)*

### `EC` — EMD + CC Composite

*(Declared but implementation delegates to existing EMD logic.)*

### `ED` — EMD + DERivative Composite

*(Declared but implementation delegates to existing EMD logic.)*

### `W` — Winding Agreement

A manifold-aware step-coherence metric. Identifies turning points (local extrema) in a reference `Manifold` series, then for each interval between consecutive turns:

1. Compute the step angle (direction) for Model, Data, and Manifold displacements
2. Measure angular mismatches (wrapped to [−π, π]) weighted by step amplitude
3. Apply exponential penalty: `Q = exp(−E / 2σ²)` with σ = π/4
4. Add a sign-coherence bonus (fraction of pairs with matching signs)

Final score:
```
Score = 0.45·Q(Model,Data) + 0.25·Q(Model,Manifold) + 0.20·Q(Data,Manifold) + 0.10·SignCoherence
```

Returns values in [0, 1]. Rewards matching the *hidden step structure* defined by the manifold, not just pointwise values.

### `HOYER` — Hoyer Spectral Peak

Measures spectral sparsity/peakedness via the Hoyer index applied to the data's power spectrum:

```
H = (√N − ‖S‖₁/‖S‖₂) / (√N − 1)
```

where `S` is the power spectrum. Returns values in [0, 1]. 1.0 means a single spike (maximally sparse/peaked), 0.0 means flat spectrum. Uses the same frequency grid as FT/ME metrics.

### `ZS` — Zone Sharpness (blended with CC)

Added 2026-11-19. Scores how clearly the data collapses onto one pattern when the post-Bessel manifold is folded
into one Brillouin zone, x = frac(k·F), with k the leg's base winding (its last `ltep`). Blended 50/50 with CC, so
the search sharpens the folded pattern while keeping the time-domain fit:

```
Score = 0.5 · CC(Model, Data) + 0.5 · ZS,     ZS = 0.5 · (η²_cv + 1 − C)
```

- **Data:** IR-matched, y(n) + IR·y(n−12), only where the two samples are a year apart. The least-squares fit on
  1, t, t² and F is removed, since none of these is periodic in x.
- **η²_cv:** share of variance explained by x, from 40 bins. The bin means come from even years and are scored on
  odd years, and the reverse.
- **C (Geary's C along the zone):** the data sorted by x. The mean squared difference to the 6 nearest neighbours
  in x that are at least 2 years apart in time, over 2·var.
- **Why the year split and the 2-year gap:** F moves in yearly steps, so a year's 12 months share one x. Without
  them the score would reward the record's own month-to-month persistence. With them, random fold wavenumbers
  score about 0 (mean −0.04 to −0.01 on the 6-leg TNA extrapolate fit).
- **Legs:** each leg uses its own k and IR. The composite averages them, so the shared lpap is pushed to sharpen
  every leg's pattern at once.
- **Check:** on that fit, the solver's ZS equals a Python reimplementation to 5 decimals on all 6 legs. A
  2-minute search raised the mean ZS from 0.206 to 0.253, while CC fell from 0.776 to 0.766 and dLOD stayed at
  0.984. The gain came from the per-leg Bessel and winding parameters; lpap moved by at most 1e-5 relative.
- **Cost:** about 30% fewer search cycles per second than `CC`.
- **Reported scores:** under `METRIC ZS`, the held-out and validate scores are also ZS blends. For a short window
  (e.g. 10 years), η²_cv has few points per bin.

---

## Summary Table

| Code  | Name                       | Range     | Higher means…                        | Domain       |
|-------|----------------------------|-----------|--------------------------------------|--------------|
| `CC`  | Pearson Correlation        | [−1, 1]   | linear correlation                   | time         |
| `RMS` | Inverted RMS Error         | (0, 1]    | lower RMS error                      | time         |
| `ZC`  | Zero-Crossing Agreement    | [−1, 1]   | sign coherence                       | time         |
| `DER` | Derivative Correlation     | [−1, 1]   | matching rates of change             | time         |
| `FT`  | Spectrum DTW               | (−∞, 1]   | spectral shape match                 | frequency    |
| `SEM` | Scaled Error Metric        | (0, 1]    | lower scaled error                   | time         |
| `DTW` | Dynamic Time Warping       | [0, 1]    | temporally-flexible alignment        | time         |
| `EMD` | Earth Mover's Distance     | (0, 1]    | distributional similarity            | distribution |
| `CID` | Complexity-Invariant Dist. | (0, 1]    | matching total variation             | time         |
| `CTW` | (alias of CID)             | (0, 1]    | matching total variation             | time         |
| `W`   | Winding Agreement          | [0, 1]    | manifold step coherence              | time         |
| `HOYER`| Hoyer Spectral Peak       | [0, 1]    | spectral sparsity / peakedness       | frequency    |
| `ZS`  | Zone Sharpness (+CC)       | about [−0.5, 1] | data collapse in the folded zone | folded manifold |
| `DC`  | DTW+CC composite           | —         | —                                    | —            |
| `EC`  | EMD+CC composite           | —         | —                                    | —            |
| `ED`  | EMD+DER composite          | —         | —                                    | —            |

## Related Environment Variables

| Variable       | Default    | Affects                          |
|----------------|------------|----------------------------------|
| `CC_START`     | 0.0        | Lower bound of CC evaluation window |
| `CC_END`       | 999999999  | Upper bound of CC evaluation window |
| `FMULT`        | 1.008      | Multiplicative frequency step (FT, ME) |
| `FSTEP`        | 0.18       | Linear frequency step (FT, ME)   |
| `FSTART`       | 0.01       | Starting frequency for spectra   |
| `MERMS`        | False      | ME metric sub-mode               |
| `STEP`         | True       | Linear vs. geometric frequency stepping |
