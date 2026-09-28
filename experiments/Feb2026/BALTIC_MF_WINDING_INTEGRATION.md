# Does Baltic's M=1.245 Winding Correspond to a Unity Integral Winding of the Mf-Driven Manifold Swing?

A follow-on investigation to `WINDING_NARRATIVE.md` section (8) and
`AMO_FORMULATION.md`'s Coulomb-friction resonance work, asking a more
specific question: can Baltic's dominant winding number, `M=1.24538`
(`kN_baltic`'s own independent fit: `M=1.24674`) — the 6th harmonic of the
shared `M(NM)≈0.2076` backbone found across seven climate indices — be
derived, even approximately, from the raw Mf (13.66083d fortnightly lunar
tide) forcing mechanics, rather than only cross-validated empirically as
in section (8)?

**Verdict up front: no.** Along the way, a real winding-convention bug was
found and corrected (section 2) — the true Ada argument is
`sin(2π·M·Forcing)`, not bare `sin(M·Forcing)` — and even after fixing
that, a striking-looking 0.064%-from-`2π` match (section 4) turned out to
be a false positive once properly sensitivity-checked, not a genuine
resonance. A follow-on lattice/Brillouin-zone hypothesis (section 5) —
does the manifold's discrete staircase form a regular lattice whose
reciprocal-space zone edge mode-locks to `M=1.245`, Arnold-tongue style —
is rejected even more decisively: the staircase isn't a regular lattice
in the first place, and the direct aliasing test (comparing the actual
basis vectors, not just magnitudes) shows no relationship at all. This
document records the full investigation, since the convention correction
and both false-positive/rejection catches are informative about how this
manifold's scale actually works and how easy it is to mistake a
coincidence for a finding here.

## 1. The mechanism confirmed first (solid ground)

Two Doodson constituents identified exactly by period, out of the 42-term
set (`DOODSON_ARGS` in `lte_forward.py`):

| name | row | (s,h,p,N) | period |
|---|---|---|---|
| Mf (fortnightly) | 4 | (2,0,0,0) | 13.66083d (exact machine match to `JERK_REFERENCE_PERIOD`) |
| Mt | 19 | (3,0,-1,0) | 9.13295d (matches `AMO_FORMULATION.md`'s Mt=9.133d) |

Baltic's own calibrated impulse comb (`impulse_delta`, sampling=12/yr)
fires **exactly once per year**, at the *"asym"* branch position (month
index 8 of 12), not the nominal "delA" branch — the calibrated `delB`
places that branch's position (`dpos=14`) outside the valid `0..12`
range entirely, so only the semi-annual companion pulse ever fires in
practice. This is a real, literal **annual stroboscope**: once a year, at
a fixed calendar phase, the impulse comb samples the instantaneous value
of the fast Mf-dominated tide sum and feeds it into a near-perfect
integrator (`lag_a = 1 - ma ≈ 0.99999`).

Because Mf does not complete an integer number of cycles per year
(`year_length/Mf = 26.7368`, fractional part `0.7368`), that annual sample
sweeps slowly through a full Mf cycle — a classic stroboscopic alias, with
predicted beat period

```
T = 1 / |1 - frac(year_length/Mf)| = 3.7988 yr
```

This is not a guess: the actual periodogram of the real, stroboscopically
-sampled Mf-only signal (145 discrete annual values, isolated by zeroing
every `lpap` row except row 4) peaks at **T=3.794yr** (best-fit, fine
grid) against the **T=3.7988yr** predicted from Mf's period and the
calibrated year length alone — a 0.13% match. This part of the mechanism
is solid and already load-bearing evidence elsewhere in this project
(`WINDING_NARRATIVE.md` section 8's SST-persistence cross-validation).

## 2. A convention correction, made mid-investigation

**The real Ada winding argument includes an explicit `2π`.** The actual
formula, `gem-lte-primitives.adb:533` (`function LTE`):

```ada
SW := Sin (2.0 * Pi * Wave_Numbers (J) * Res (I).Value + M.Phase) * ...
```

— confirmed identical in the verified Python port, `lte_forward.py`'s
`lte_response()` and `regression_factors()` (both use
`2 * math.pi * m[k] * fv`). `winding_scalogram.py`'s own transform basis
(`exp(-i·2π·M·Forcing)`) matches this. Everything in section 1 above (the
Mf/Mt identification, the beat-period confirmation) was built on
`cv_rolling_blocked.prepare()` + `regression_factors`/`lte_response`, so
it already used the correct convention throughout and is unaffected by
what follows.

**A separate, unrelated tool in this project — `grid_lasso_optimizer.py`
and `pysr_optimizer.py`, written fresh this session — used bare
`sin(M·Forcing)`, missing the `2π` factor entirely.** Since `6 ≈ 2π`
(4.7% apart), a bare-convention fit that actually recovers the
*fundamental* `M(NM)≈0.2076` produces a bare-M value of
`M(NM)×2π≈1.304` — close enough to the *real* printed
`K6 = M(NM)×6 = 1.245` to look like a match, but for the wrong reason.
This means that investigation's earlier "recovery" of M≈1.245 (including
recovering it exactly on a synthetic calibration target built with the
same bare formula) never actually tested against the true Ada convention
at all — it was self-consistent with its own bug, not a confirmation.
This is a real, useful catch but a separate one from the analysis below,
which uses the *correct* convention throughout: since the true argument
is `2π·M·Forcing`, the number of full rotations swept as `Forcing`
changes by `ΔF` is `2π·M·ΔF / (2π) = M·ΔF` — **no additional division by
`2π`** (an error present in an earlier draft of this document, corrected
below).

## 3. The hypothesis under test

If the manifold's local swing over one beat cycle (`T ≈ 3.8yr`) is some
characteristic peak-to-peak magnitude `ΔF`, does `M · ΔF ≈ 2π` — i.e.,
does the fitted winding number correspond to almost exactly **one full
2π rotation** of the true `sin(2π·M·Forcing)` argument per beat cycle? A
clean match here would be a genuinely strong, falsifiable claim, in the
same spirit as the exact-harmonic (`K6 = 6×M(NM)`) and exact-period
(Mf/Mt row identification) matches already confirmed elsewhere.

## 4. Method 1 — direct measurement, and why the first version was a false positive

Reconstruct the Mf-only forcing (via `tide_sum` → `impulse_delta` → `iir`,
Baltic's own calibrated `delA`/`asym`/`ma`/`mp`, all other `lpap` rows
zeroed). At the specific choice of a 3.7988yr window, non-overlapping,
starting at the record's first sample:

```
mean local p2p (Mf-only, 3.7988yr windows) = 5.042
M(K6=1.24538) × 5.042 = 6.2798   vs   2π = 6.2832   (0.064% off)
```

That is a striking number — but a sensitivity sweep (the same discipline
`AMO_FORMULATION.md` applies to its own year-length and `mp` findings)
shows it does not survive:

| window width | mean p2p | M×p2p |
|---|---|---|
| 3.0yr | 4.471 | 5.568 |
| 3.4yr | 5.084 | 6.331 |
| **3.7988yr** | 5.042 | **6.280** |
| 4.0yr | 5.365 | 6.682 |
| 4.5yr | 5.540 | 6.900 |
| 5.0yr | 6.061 | 7.548 |
| 6.0yr | 6.504 | 8.100 |

`M×p2p` climbs **smoothly and monotonically with window width** — exactly
as expected for a near-integrator accumulating a noisy input (wider
windows capture more accumulated wander, full stop) — with no plateau,
kink, or distinguishing feature at 3.7988yr specifically. The apparent
match to `2π` is one point on a smooth, unremarkable curve, not a
resonance.

Worse, at **fixed** width (3.7988yr) the result is not even stable to
which sample the windows start on:

| window start offset | mean p2p | M×p2p |
|---|---|---|
| 0.0yr | 5.042 | 6.280 |
| 0.5yr | 5.209 | 6.487 |
| 1.0yr | 5.092 | 6.341 |
| 1.5yr | 4.744 | 5.908 |
| 1.9yr | 4.966 | 6.184 |
| 2.5yr | 5.001 | 6.228 |
| 3.0yr | 5.256 | 6.546 |

A ±5% spread just from window placement, and the finer-grained,
most-averaged estimate (1695 overlapping windows, 0.1yr step) lands at
`M×p2p = 6.317` — 0.5% from `2π`, not the original 0.064%. **This is the
same false-positive pattern already caught once earlier in this
investigation (a single-lag check that looked like a π/2 resonance and
wasn't): a smoothly-varying quantity happening to pass near a round
number at one specific, otherwise-unremarkable parameter choice.**

## 5. A lattice/Brillouin-zone hypothesis, tested and rejected

A further, more structural hypothesis: `Forcing(t)` is literally a
staircase (constant between the once-a-year impulse updates). Does the
*discrete set* of 145 plateau values it visits form a regular lattice —
and if so, does the associated reciprocal-space zone boundary
(`M_BZ = π/Δ`, the Brillouin-zone analogue of a lattice with spacing `Δ`)
line up with the fitted `M=1.245`, in the sense of genuine mode-locking
(an Arnold-tongue-style phenomenon)?

**First check — is there a regular lattice at all?** Sorting the 145
plateau values and measuring nearest-neighbor spacing in *value* space:
mean spacing `0.373`, but `std/mean = 1.24` — wildly irregular (a real
lattice has `std/mean ≈ 0`). The histogram is visibly multi-modal/clumpy
(heavy concentration around [-47,-42] and [-9,0], sparser elsewhere),
consistent with a random-walk-like process that revisits some regions
far more than others, not one that marches through evenly-spaced levels.
**The literal crystallographic Brillouin-zone concept does not apply —
there is no regular lattice for it to act on.** Naively computing
`M_BZ = π / (value-sorted spacing)` confirms this is the wrong quantity:
it gives `8.4`–`13.6`, off from `M=1.245` by 500-1000%.

A more physically defensible version uses the *temporal* (sequential,
year-to-year) step size as the effective lattice constant instead of the
value-sorted gap — `Δ = mean|step| = 2.6446` (§ previous discussion) —
giving `M_BZ = π/Δ = 1.190`, only 4.4% from `M=1.245`. Suggestive, but
this project has already caught two separate ~4-6% "matches" that turned
out to be coincidental (§4 here, and the earlier 6-vs-2π confusion), so
a magnitude comparison alone is not enough to trust it.

**The decisive test**: aliasing has a precise meaning — two winding
numbers `M` and `M'` are aliased under a discrete sampling iff
`sin(2π·M·x)` and `sin(2π·M'·x)` are the *same vector* at the sampled
points, not merely similar in magnitude. Testing this directly at the
candidate alias `M' = M + 2π/Δ = 3.6212`, evaluated at the 145 real
discrete plateau values:

```
corr(sin(2π·M·plateau), sin(2π·M'·plateau))  = +0.112
corr(cos(2π·M·plateau), cos(2π·M'·plateau))  = -0.159
```

Nowhere near `±1`. Worse, an arbitrary control with no lattice
justification at all (`M+1`) gives an equally weak correlation (`-0.054`)
— **the candidate alias is not behaving any differently from a random
other frequency.** Single-mode fit quality against real Baltic MSL data,
compared at `M`, at the candidate alias, and at that same arbitrary
control, on both the 145 discrete points and the full monthly record,
shows the same story: no distinguishing relationship anywhere.

This decisively rejects the hypothesis, more conclusively than the
magnitude-based checks elsewhere in this document — it fails the actual
defining test of aliasing, not just a numerical-closeness heuristic. It
converges with the first check's finding: there is no regular lattice in
Baltic's staircase for a Brillouin-zone/mode-locking argument to act on
in the first place.

## 6. Verdict

| Claim | Status |
|---|---|
| Mf (row 4) and Mt (row 19) identified exactly by period | Confirmed |
| Annual-stroboscope aliasing predicts a 3.7988yr beat period | Confirmed — matches the real periodogram peak to 0.13% |
| That beat signature persists in real SST data before Baltic MSL data exists | Confirmed (`WINDING_NARRATIVE.md` §8) |
| The real Ada/`lte_forward.py` winding argument is `2π·M·Forcing`, not bare `M·Forcing` | Confirmed directly from source (`gem-lte-primitives.adb:533`) |
| `grid_lasso_optimizer.py`/`pysr_optimizer.py`'s bare-`M` convention doesn't match the real model | Confirmed — a real bug in those two tools, unrelated to the mechanism above |
| M=1.245 corresponds to a unity (`2π`) winding per local 3.8yr manifold swing | **Not confirmed** — the number moves smoothly and substantially with both window width and window placement; the apparent 0.064% match was one point on that unremarkable curve, not a distinguishing feature |
| The manifold's discrete plateau values form a regular lattice | **Rejected** — nearest-neighbor spacing is wildly irregular (std/mean=1.24), histogram is multi-modal |
| M=1.245 mode-locks to a Brillouin-zone edge of that (candidate) lattice | **Rejected** — the decisive aliasing test (correlation of the actual basis vectors at the discrete points) shows no relationship, indistinguishable from an arbitrary control frequency |

The honest reading is unchanged from the previous draft's conclusion,
now reached the right way: the annual-stroboscope/beat-period *mechanism*
is real and precisely verified, but that does not make every quantity one
can construct from it — including "windings per local swing," under
either the correct or the previously-erroneous formula — similarly
precise. `M`'s specific value is fit jointly with the other free
parameters shaping `Forcing`'s scale (it has no fixed reference dimension
of its own); nothing in this integration-theory route pins that scale
down independently of the regression that already produced it. The
`mp`-sensitivity resonance (a genuinely sharp, ~1%-wide, reproducible
peak that does *not* move around under a matching sensitivity check)
remains the strongest confirmed link between a fixed-scale parameter and
the M=1.245 winding's fit quality — this document's own attempted
mechanism does not add to that.
