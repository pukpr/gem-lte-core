#!/usr/bin/env python3
"""signal_operators.py — small, genuinely cross-index operators, factored
out here (rather than duplicated per-script) as the first concrete step
toward unifying the resolved shallow-water line of work across indices.

Every resolved model this session (baltic/amo/nino4/qbo/pdo) has built
its own bespoke fixed/moving-gauge extraction inline; those stay
per-script because each domain's geometry differs. What DOESN'T need to
differ is what you do to a gauge's OUTPUT once you have it -- starting
with the delayed-difference operator this module exists for.
"""
from __future__ import annotations

import numpy as np


def delayed_difference(x: np.ndarray, lag: int) -> np.ndarray:
    """D(t) = x(t) - x(t-lag). A first-difference-at-lag operator -- its
    transfer function has magnitude |2*sin(pi*f*lag*dt)|, which is EXACTLY
    zero at f=0 and at every integer multiple of 1/(lag*dt), and reaches
    its own maximum (2x) at half-integer multiples. For lag=12 on monthly
    data (lag*dt=1yr), this kills DC and every exact-integer-year period,
    while passing periods well off those multiples largely unchanged --
    see `transfer_response` for the exact attenuation at any given period.
    """
    return x[lag:] - x[:-lag]


def transfer_response(period_yr: float, lag_dt_yr: float) -> float:
    """Closed-form fractional amplitude passed by a delayed-difference
    operator (lag expressed in the same time units as period_yr) at a
    given period: |2*sin(pi*lag/period)| / 2, normalized to 1.0 at the
    operator's peak response (half-integer multiples of the lag) so the
    return value reads directly as "fraction of amplitude passed" in
    [0, 1]. E.g. transfer_response(60, 1.0) ~= 0.05 -- at a 12-month
    (1-year) lag, a 60-year period signal keeps only ~5% of its
    amplitude, which is the quantitative form of "a 12-month delayed
    difference removes a 60-year cycle."
    """
    return float(abs(np.sin(np.pi * lag_dt_yr / period_yr)))


def dipole(field_a: np.ndarray, field_b: np.ndarray) -> np.ndarray:
    """A spatial dipole difference between two sampled locations (e.g.
    Azores-minus-Iceland for NAO) -- the spatial analogue of
    delayed_difference: instead of differencing a single location across
    TIME, this differences two locations at the SAME time. Both operators
    remove whatever component is common/slowly-varying across the thing
    being differenced (past self, or the other pole), which is exactly
    why either one can plausibly explain why a derived index (NAO) lacks
    a slow cycle (AMO's ~60yr) that the underlying field still carries."""
    return field_a - field_b
