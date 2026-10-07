# AMO/PDO common-manifold round-robin, v2: declaration (2026-10-06, before running)

Replicates ../_rr/round_robin.py (2026-10-02) with explicit acceptance, tolerance and exit criteria.

Shared manifold (as before): lpap, delA, delB, asym, ma, mp, shfT, init, year, impA, impB, offs, bg,
ltep[0] (layer k), ltep[1] (Bessel k). Per index: ltep[2..] + harm, IR, impC, annual terms and the
linear outer coefficients. Starting fits: Feb2026/amo and Feb2026/pdo as of this declaration
(copied to start/). Each index uses its own resp and lte_run.sh settings unchanged (METRIC=DTW,
EXCLUDE 2000-05, VALIDATE). Frozen binary ../../Oct2026_layer_ir/lt.exe.

Phase 1, noise floor: 3 replicate 300 s searches of each index from its own start. eps_i = the
largest absolute difference in VALIDATE score between replicates and the start re-score.
eps = max(eps_amo, eps_pdo, 0.005).
Reference: s_i(own) = VALIDATE score of index i's own start (TEST_ONLY).

Phase 2, round-robin (budget 12 rounds), alternating AMO, PDO, ...:
  search index X 300 s from (current common manifold + X's own outer terms), with MANIFOLD
  anchored to the current common column 4; carry X's shared parameters to every index; re-score
  each with TEST_ONLY (outer linear terms re-solved). Sharing penalty
  Delta_i = s_i(own) - s_i(common), on VALIDATE.
  Accept the round iff dLOD >= 0.99 for both AND max_i Delta_i decreases by more than eps.
Exits:
  SUCCESS: max_i Delta_i <= eps (with dLOD >= 0.99 for both).
  PLATEAU: 2 consecutive rejected rounds -> "not shareable at tolerance eps"; report Delta_i.
  BUDGET: 12 rounds.
Train and test (2000-05) are reported but play no part in the decisions.

Certification (separate; never used for decisions): on SUCCESS, the 1880-1949 back-cast and the
frozen-manifold nulls (shifted clock, fake Moon) are run on the common manifold.

AMENDMENT (2026-10-06 10:38, after the own-fit re-scores and before any replicate or round result):
the own fits score dLOD 0.9880 (AMO) and 0.9889 (PDO), below the declared 0.99, so no round could
ever pass. The gate is changed to dLOD >= 0.985, i.e. within 0.003 of the own fits. Nothing else
changes.
