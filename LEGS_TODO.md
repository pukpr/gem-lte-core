# Leg groups: to do

Companion to `LEGS.md`. Written 2026-10-12.

## Open issue: the legs share windings but not the response

In a two-leg run the manifold parameters, the windings (`ltep`) and the harmonic list are common, but
each leg regresses its own amplitudes and phases. For two samplings of the same index this leaves the
interpretation uncertain: a good composite does not show that one response fits both records.

Measured on the 5-day + monthly SOI result (`experiments/Oct2026_generator/legs/soi_5_m`,
`leg_agreement.py`):

| Windings | Agreement of the two legs' amplitudes and phases |
|---|---|
| All 15 | +0.33 |
| Harmonics up to 30 | +0.52 |
| Harmonics above 30 | −0.07 |

For two different indices (AMO and PDO, say) different amplitudes and phases are expected and this is
not an issue.

## To do

1. **Common response option for same-index legs.** A switch (suggested `LEG_COMMON TRUE`) with one of:
   - *Joint regression (preferred).* Each leg builds its normal equations; the `Exchange` sums them
     with `LEG_WEIGHT` and solves once; both legs use the same amplitudes and phases. Needs
     `Regression_Factors` split into a build step and a solve step.
   - *Leader's amplitudes (quick check).* The follower skips its regression and applies the leader's.
     Leg 2 then only tests leg 1's response.
2. **Per-leg scale factor** for the common response, so the two series need not be in identical units.
3. **Report the agreement at save.** Print the leg-to-leg amplitude agreement next to `COMPOSITE`, so
   every run states how common the response is.
4. **Composite options.** `min` of the two legs as an alternative to the weighted mean, so one leg
   cannot be traded for the other.
5. **Cross-interval scores in the status line** (each leg's score on the interval only the other leg
   trains), for display only.
6. ~~**Save race** between pairs finishing together (see Limits in `LEGS.md`): have the best pair's
   follower save under the leader's control.~~ Done 2026-10-12: the saving pair is chosen centrally.
7. ~~**Test-mode exit when standard input is an idle pipe.** The keyboard read blocks; exits normally
   from a terminal or with `< /dev/null`.~~ Done 2026-10-12: a test run no longer reads the keyboard.
8. **`VALIDATE`, `ALTERNATE`, `SPLIT_TRAINING`** with legs; **`TIDES ZONAL`** forcing cache with two
   time grids.
9. **Plot scripts and `lte_gui.py`**: accept `lte_results.1.csv` / `.2.csv`.
10. Done 2026-10-12 (`LEGS N`, 2 to 9). ~~**More than two legs** (three or more indices on one manifold): the `Exchange` generalizes to one
    leader and several followers.~~

## Next experiment

Two legs on two different indices, monthly AMO and PDO, replacing the earlier round robin: seed with
the slow winding, the base winding (~0.2) and its harmonics 2 to 9; common manifold from 1880.
Run directory: `experiments/Oct2026_generator/legs/amo_pdo`.
