# Winding Harmonics in LTE

How the harmonics of the base winding enter a fit, and how Brillouin-zone folding shows them,
using the PDO fit on the dLOD-calibrated tidal table as the worked example
(`../Oct2026_generator/indep/pdo`, record 1880–2022, tidal table 0.9995 correlated with
measured dLOD). Statements under "Established" are computed from that fit or from the tidal
periods. The "Interpretation" section is labelled as such.

Figures, both in `../Oct2026_generator/indep/pdo/`:

- `pdosite1940-1970.png`, from `lte_gui.py` / `plot.py` with the Brillouin-zone folding option on (the `ZONE` setting; called BZONE here to keep it apart from the ZONAL tide mode) (sections 5 to 7);
- `winding_scalogram.png`, from `winding_scalogram.py` (section 8).

Numbers for sections 5 to 7: `../Oct2026_generator/windharm/zonefold_pdo.py`.

## Summary

- The annual impulse samples the fast lunar tides once a year, so each tide is seen at its
  alias. Mf and Mm alias to about 3.8 and 3.9 years; Mt, their sum frequency, aliases to
  about 127 years.
- Integrating those yearly samples gives the manifold F(t), a staircase with one step per year.
- Folding the manifold into one Brillouin zone, θ = k·F mod 1 with k the base winding, turns
  the response into a single curve g(θ). One zone is 1/k = 4.96 units of F.
- On average the manifold advances by about one zone per year (k × mean annual step = 0.95),
  with wide scatter, so the zone is covered.
- The harmonics h·k are the Fourier series of g(θ). For PDO the sixth harmonic dominates, in
  the model and in the data folded the same way.
- The winding scalogram scans all winding numbers without assuming harmonics of k. Its
  strongest ridge in the PDO data is at M = 1.21 = 6.00·k, present through the whole record.

## Established

### 1. Tidal aliases under annual sampling

With the model year of 365.2463 days, a line of period P has f = 365.2463 / P cycles per year
and is seen at f − round(f).

| Constituent | Period (days) | Cycles per year | Alias (cycles/yr) | Alias period (yr) |
|-------------|---------------|-----------------|-------------------|-------------------|
| Mf | 13.6608 | 26.737 | −0.2632 | 3.80 |
| Mm | 27.5546 | 13.255 | +0.2554 | 3.92 |
| Mt | 9.1330 | 39.992 | −0.0079 | 127 |

Mt is the sum frequency of Mf and Mm (1/13.6608 + 1/27.5546 = 1/9.1330), and the aliases obey
the same sum: −0.2632 + 0.2554 = −0.0079. The 127-year term is therefore the beat between the
3.80- and 3.92-year aliases.

### 2. Origin of the three lines

The forcing is the degree-2 zonal potential,

```
V(t) ∝ (a/r)³ · P₂(sin δ),      P₂(x) = (3x² − 1) / 2
```

with δ the Moon's declination and r its distance.

- **Mf** comes from the square of the declination. For sin δ ∝ sin(ωt),
  `P₂ = 1/4 − (3/4)·cos(2ωt)`, a line at twice the Moon's mean-longitude rate.
- **Mm** comes from the distance factor (a/r)³, at the anomalistic rate.
- **Mt** comes from the product of the two, and from the elliptic inequality in the Moon's
  longitude, which puts the argument 2L + M′ inside sin²δ.

### 3. The manifold

The tidal forcing is multiplied by the annual impulse and integrated. The manifold changes in
about one month per year and is flat otherwise.

| Quantity | Value |
|----------|-------|
| Base winding k | 0.201774 |
| Zone width 1/k | 4.96 |
| Manifold range, 1880–2022 | 39.3 (about 8 zones) |
| Mean absolute annual step | 4.73 |
| Spread of the annual step (sd of the absolute value) | 3.21 |
| k × mean absolute annual step | 0.95 zones |
| Annual step in zones, 10th / 50th / 90th percentile | −1.57 / −0.05 / +1.51 |

The same product, k × mean annual step, is 0.93–0.96 in the AMO, NINO4, IODE and NAO fits on
the same table, which share this manifold and nearly this k.

### 4. The response

```
model(t) = level + k0·F(t) + trend·t + Σ_j A_j · sin(2π·k_j·F(t) + φ_j)
```

followed by the lag-12 term, `model(t) − IR·model(t − 12 months)`, with IR = −0.434. The
windings k_j are one low winding (−0.064) and the base winding k with its harmonics 2 to 9.
The amplitudes and phases come from a linear regression at each evaluation; k, the low winding
and the list of harmonics are searched.

The fit scores 0.772 over the training record and +0.545 over 2000–2005, which the search did
not see.

### 5. Brillouin-zone folding

Every term h·k with integer h has the same value at F and at F + 1/k. The sum over harmonics
therefore depends on F only through the zone phase,

```
θ(t) = k·F(t) mod 1,        g(θ) = Σ_h A_h · sin(2π·h·θ + φ_h)
```

so the eight zones the manifold crosses in 143 years collapse onto one curve. Two panels of the
figure show this:

- **"Latent/hidden forcing layer"** plots θ(t). The yearly steps of about one zone, with their
  scatter, appear as jumps that land throughout the interval 0 to 1. Over the record, six equal
  bins of θ hold 26, 19, 20, 21, 23 and 34 of the 143 yearly values.
- **"LTE Modulation"** plots the response against θ, for the model and for the data. Both have
  the lag-12 term put on the same footing: the model with it undone, the data as
  D(t) + IR·D(t − 12), which is the quantity the regression is fitted to. The model traces
  g(θ) directly; the data scatter about it.

Because the zone is covered, the harmonic amplitudes are the Fourier coefficients of g, and the
same decomposition can be applied to the data without reference to which harmonics the model
was given.

### 6. Harmonic spectrum of the folded response

Each series is regressed on a background (constant, F, t, t²) plus sin and cos of 2π·h·k·F for
h = 1 to 12, as in the "Harmonic spectrum of the zone-folded modulation" panel.

| h | Data, folded | Model regression |
|---|--------------|------------------|
| 1 | 0.179 | 0.153 |
| 2 | 0.110 | 0.120 |
| 3 | 0.233 | 0.218 |
| 4 | 0.200 | 0.226 |
| 5 | 0.269 | 0.283 |
| **6** | **0.748** | **0.776** |
| 7 | 0.159 | 0.179 |
| 8 | 0.355 | 0.397 |
| 9 | 0.112 | 0.146 |
| 10 | 0.034 | not in the model |
| 11 | 0.071 | not in the model |
| 12 | 0.053 | not in the model |

- The sixth harmonic is the largest in both, by about a factor of two over the next (h = 8).
- It is an interior peak. In the data, harmonics 10 to 12, which the model does not contain,
  are the smallest of the twelve.
- The earlier PDO fit in this directory (`pdo/`, a different tidal table) also has h = 6
  largest.

### 7. Stability of the sixth harmonic across the record

The data-side decomposition repeated on separate parts of the record, with the manifold and k
unchanged:

| Era | Largest harmonic | h = 6 amplitude | h = 6 phase (rad) | h = 6 relative to the mean of the other eleven |
|-----|------------------|-----------------|-------------------|-----------------------------------------------|
| 1880–1940 | 6 | 0.865 | −1.41 | 4.5 |
| 1940–1970 | 6 | 0.711 | −1.65 | 3.0 |
| 1970–2022 | 6 | 0.526 | −1.11 | 2.8 |

The sixth harmonic is the largest in each era and its phase stays within about 0.5 rad. Its
amplitude declines through the record.

For scale, phase-randomized surrogates of the PDO series (same power spectrum), folded on the
same manifold, give an h = 6 amplitude of 0.130 ± 0.069, and their largest harmonic of the
twelve averages 0.228; none of 2,000 reaches 0.748.

**Limit of both checks.** The manifold and k were searched against the whole PDO record, so the
eras are not independent of the fit, and the surrogates were folded on a manifold tuned to the
real series, not searched themselves. These show that the six-lobed curve is consistent through
the record and far above what an untuned series gives. They do not replace a fit made with part
of the record withheld.

### 8. The winding scalogram

Zone folding assumes the response is built from integer harmonics of one base winding. The
winding scalogram does not. For every winding number M on a grid from 0 to 5 in steps of 0.01,
and every 15-year Gaussian window centred at t0, it computes

```
G_X(t0, M) = < w(t − t0) · X(t) · exp(−i·2π·M·F(t)) > / < w(t − t0) >
```

for X the standardized data and the standardized model (the model with the lag-12 term undone).
This is a windowed Fourier transform taken along the manifold F instead of along time, so M is
conjugate to the forcing level, and |G|² shows when in the record a winding at rate M is
present. The reference curve in the left panel is the same transform applied to AR(1) red noise
(lag-1 coefficient 0.97), which measures how much power any persistent series shows at each M
on this manifold.

The zone-folded harmonic spectrum of section 6 is the time-average of this transform read off
at M = h·k only. The scalogram adds two things: whether the ridges fall on multiples of k when
M is free, and how each ridge behaves through time.

**Where the ridges fall.** Time-mean power, strongest first:

| Ridge in the data (M) | M / k | log2 power | Bits above the AR(1) reference | Same ridge in the model (bits above) |
|-----------------------|-------|------------|--------------------------------|--------------------------------------|
| 1.21 | 6.00 | −2.9 | +2.8 | +3.6 |
| 0.61 | 3.02 | −4.4 | +1.2 | +1.1 |
| 1.62 | 8.03 | −4.6 | +2.2 | +3.0 |
| 0.82 | 4.06 | −5.2 | +1.0 | +1.9 |

- The four strongest ridges in the data sit at 6, 3, 8 and 4 times the base winding, to within
  the grid spacing. The scan was not told about k.
- The sixth harmonic is the strongest ridge anywhere on the grid, in data and model.
- Only 1.6% of the grid exceeds the red-noise reference by more than 2 bits, and the median
  excess is zero, so the spectrum is not broadly elevated.
- Harmonics 1, 2 and 9 to 12 are at or below the reference in the data, in line with their
  small amplitudes in section 6.

**Through time.** Along the M = 6·k ridge, for the 22 window centres clear of the record edges
(1897 to 2002), with the red-noise reference at −5.7:

| Window centre | Data, log2 power | Model, log2 power | Data phase (rad) |
|---------------|------------------|-------------------|------------------|
| 1897 | −2.7 | −1.9 | −3.02 |
| 1922 | −2.8 | −2.3 | −2.77 |
| 1947 | −2.7 | −2.3 | −2.67 |
| 1972 | −2.8 | −2.0 | −2.56 |
| 1987 | −3.2 | −1.8 | −2.58 |
| 2002 | −3.6 | −1.8 | −2.82 |

- The ridge is continuous across the record and stays 2 to 3 bits above the reference.
- Its phase moves by less than 0.5 rad in 105 years.
- The data ridge weakens after about 1975 while the model ridge does not, which is the same
  decline seen in section 7.

**Limits.** Adjacent windows overlap heavily (15-year width, 5-year spacing), so the 22 values
are not independent. With a manifold range of 39.3, winding numbers closer than about 0.03
cannot be separated. As in section 7, the manifold was searched against the whole record.

## Interpretation

These points are not established by the above.

- **The harmonics as a mode shape.** If g(θ) is the response of the basin to the forcing level,
  a dominant h = 6 means that response has six lobes across one zone. That would make the
  dominant harmonic a property of the basin, which is consistent with it differing between
  indices (6 for PDO, 4 for NAO).
- **One zone per year.** k × mean annual step ≈ 1 may be a condition the search converges to,
  or it may be inherited from the shared starting values. A search started from a different k
  would distinguish the two.
- **The highest allowed harmonic.** In the fits where the largest harmonic is the top of the
  allowed range (9 for AMO, NINO4, IODE and many SST quads), the amplitudes rise steadily toward
  the limit. That pattern is what truncation of a longer series would produce, so those cases
  do not identify a dominant mode until higher harmonics are allowed. Folding the data, as in
  section 6, is the direct way to see whether higher harmonics are present.
- **The declining amplitude** of the sixth harmonic from 0.87 to 0.53 across the record, also
  visible as the weakening data ridge after about 1975 in the scalogram, could be a real change
  in the response, a change in data quality, or the manifold drifting out of register with the
  true forcing toward one end. The model ridge does not weaken, so the fit does not reproduce
  it, and the present fit does not separate these causes.

## Open checks

1. Refit PDO with 1940–1970 withheld from the search and repeat sections 6 and 7 on the
   withheld years.
2. Start a PDO search from a base winding well away from 0.2018 and see whether k × mean annual
   step returns to about 1.
3. Apply the data-side folding of section 6 and the scalogram of section 8 to the other
   indices and the SST quads, to find which have an interior dominant harmonic, which are
   limited by the allowed range, and whether their ridges also fall on multiples of k.
4. Run the scalogram with the manifold of a fit that withheld part of the record, and read the
   ridge inside the withheld years.
