# The QBO read as a control problem (2026-10-26)

Code: `qbo_control.py` (figure `qbo_control.png`, numbers `qbo_control.json`). Plant: the linear-edge QBO device of
`qbo_tokamak.md` (`COMB|linear`, parameters from `qbo_tokamak.json`). Data: KIT radiosonde winds at 30, 50, 70 hPa.

This is a supplementary view, borrowed from the Chandler-wobble counter-torque analysis
([ChandlerWobble docs](https://github.com/pukpr/ChandlerWobble/blob/main/docs/README.md), `cw_counter_torque.py`).
That analysis runs the plant backwards to find the excitation the observed wobble required, budgets it against the
lunar torque, and lets a pulsed controller at the comb instants remove the wobble. In steady state the controller fires
the aliased lunar pulse train, reversed. The question here is what the same steps say about the QBO, and whether a
QBO could be controlled the way tokamak plasmas are.

## 1. Running the plant backwards

The device is linear in its pulses. The comb fires twice a year, on day 239 and half a year later. Each pulse aₙ
passes through a leaky integrator (τᵢ = 0.76 yr), reaches the edge, and penetrates downward by advection–diffusion.
So each level's wind is uₗ(t) = cₗ Σ aₙ kₗ(t − tₙ).

Solving for the 173 pulse amplitudes (1940–2026) that reproduce all three levels together, with ridge regularisation
chosen by generalised cross-validation (127 effective degrees of freedom), reconstructs the winds at r = 0.92 / 0.90 /
0.82. That is expected of an inversion and is not evidence by itself. The content is in the **required pulse train
a_req** (panels a, b).

## 2. Required vs astronomical pulses

The astronomical pulses at the same instants are a_ast,n = sin β_moon · sin δ_sun · (solar distance)³, with Meeus
elements as in `qbo_tokamak.py`.

| | 1955–2005 | 2006–2026 |
|---|---|---|
| pulse-by-pulse r(a_req, a_ast) | **+0.60** (98th percentile vs shifted trains) | −0.09 (20th) |
| 2.37-yr alias component, a_req / a_ast | ratio 1.00, phase +5° | ratio 0.59, phase **+101°** |

- **To 2005 the input nature "applied" resembles the lunar comb** at the pulse level and is in phase at the alias.
  This is **not independent**: the plant's parameters, including the trigger day, were fitted on 1953–2005 to this same
  comb.
- **After 2006 it does not.** The alias component has turned by about 100°, which is the half-cycle phase slip of
  `qbo_tokamak.md` §4 expressed as a change in the required input.

## 3. The extra input (panel c)

The residual rₙ = a_req,n − κ a_ast,n is what had to be added beyond the lunar comb, with κ fitted on 1955–2005.

- Its rms, relative to 1955–2005, is 1.39 for 2006–14, 1.28 for 2015–21 and **1.99 for 2021–26**.
- The pulses beyond 2σ after 1955 cluster late: 2014.65 (+2.6), **2016.15 (+2.7, the 2016 disruption)**, 2022.65
  (−2.4), 2023.65 (+3.5), 2024.65 (−2.8), 2025.65 (+2.4). Before 2006 only four pulses exceed 2σ, none beyond 2.2.
- The 2019–20 disruption does not stand out at the pulse level.

Read as control, the 2016 disruption is a **single large out-of-schedule input**, consistent with its usual attribution
to Rossby-wave momentum from the Northern Hemisphere. Since 2021 the input has stayed noisier than at any time before.

## 4. Budget: a gate, not a momentum source

- **What the QBO needs:** at 30 hPa, |du/dt| has a median of 2.0 and a 90th percentile of 10.2 m s⁻¹ per month (at the
  transitions).
- **What the Moon supplies:** the horizontal lunar tidal acceleration peaks at 8.2 × 10⁻⁷ m s⁻². That would be
  2.2 m s⁻¹ per month only if it never reversed, and it reverses twice a day. The **zonal-mean zonal component is
  exactly zero**: the degree-2 force is the gradient of a potential whose zonal-mean part has no longitude dependence.

So the lunar input cannot supply the QBO's momentum, which comes from the vertically propagating waves. If the
astronomical comb matters, it acts as a **gate or trigger on that wave forcing**, setting when the reversals happen.
That is the tokamak control situation:
- the plasma's energy comes from its own heating and current;
- the control system (coils, pellets) supplies small inputs whose timing paces or steers the instabilities. ELM pacing
  by pellet injection is the closest example: the pellets set when edge-localised modes fire, not their energy.

This is the same conclusion as for the Chandler wobble, where the counter-torque is 10⁻⁴–10⁻³ of the lunar torque.

## 5. Can feedback control it? Not like the Chandler wobble (panels d, e)

The Chandler controller fires u_n = −g × (state just before the pulse) at the comb instants. Applied to the QBO device
from 1990:

| controller | residual QBO after 2000 (rms) |
|---|---|
| feedback, gain 0.2 / 0.5 / 0.8 / 1.0 / 1.5 / 2.0 | 1.02 / 1.04 / 1.03 / 1.01 / 0.93 / 0.82 |
| feedback, gain 3.0 | 3.39 (unstable) |
| feedforward, 25 / 50 / 75 / 100 % of the coming lunar pulse | 0.75 / 0.50 / 0.25 / 0.00 |

Feedback barely acts. **The QBO device is kick-dominated:** only 52% of its state survives each half year
(e^{−0.5/τᵢ}), so the signal at any time is mostly the newest pulse, which a controller acting at the same instant
cannot see. The Chandler wobble is the opposite: a resonance (Q ≈ 30–85) that accumulates many small kicks, which a
sampled feedback controller can damp.

Only feedforward, a controller that knows the coming pulse from the ephemeris, cancels the device QBO. Exact
cancellation at 100% is an identity for a linear plant; the point is the contrast with feedback. Tokamak disruption
and ELM control are mostly predictive or feedforward for the same reason: by the time the state shows the event, it is
too late to act.

## 6. What this view adds

- **The QBO's input is a schedule, not a state.** In this device the QBO is reproduced by an input known in advance
  (the ephemeris and the solar comb), not by a feedback loop on its own state. Its departures from that schedule
  (2016, and the noisier input since 2021) are directly measurable as extra input.
- **The lunar role is necessarily a gate.** Any version of the astronomical hypothesis must be a trigger or gate on
  wave-driven momentum (§4). The model's gated penetration (`qbo_tokamak.md` §7) is the right kind of mechanism.
- **The Chandler analogy holds for the budget, not for the dynamics.** The Chandler wobble is a resonator that
  integrates its forcing; the QBO device is a kicked, short-memory relaxer. A QBO "controller" would have to be
  predictive.

## Caveats

- The plant is the linear-edge device fitted to 1953–2005. The threshold (reversed-field-pinch) edge is not linear, so
  it cannot be inverted this way.
- The 1955–2005 agreement in §2 is in-sample with respect to the plant fit. Only the post-2006 disagreement is
  out-of-sample.
- The inversion is ridge-regularised. The 2σ residual pulses depend somewhat on the regularisation; the era rms
  values are more robust than individual pulses.
- The budget compares orders of magnitude and uses no wave-forcing data. An ERA5 EP-flux budget across the 2016
  disruption would put a number on the "extra input".
