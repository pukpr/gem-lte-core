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
