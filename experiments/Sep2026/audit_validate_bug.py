#!/usr/bin/env python3
"""Audit which VALIDATE=TRUE cells were actually hit by the Forcing/harm
cross-thread resync bug fixed in gem-lte-primitives-solution.adb on
2026-09-29 (see HANDOFF.md).

For each cell, replays the CURRENTLY SAVED lt.exe.p through
Feb2026/lte_forward.py's forward() directly (no lt.exe launch) and
compares the resulting Forcing (pre-regression, depends only on lt.exe.p's
scalars + tidal periods -- NOT on TRAIN_START/END/EXCLUDE) against the
cell's own saved lte_results.csv column 4.

For an UNCORRUPTED save, this should match to float64 noise (~1e-9 relative,
as independently verified on amo before/after the Ada fix). A cell hit by
the bug will show a much larger, structural mismatch, because its saved
Forcing came from a different thread's trial than its own Model/lt.exe.p.

Usage: python3 audit_validate_bug.py [--threshold 1e-3]
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FEB2026 = HERE.parent / "Feb2026"
sys.path.insert(0, str(FEB2026))
import lte_forward as L  # noqa: E402

import numpy as np


def audit_cell(cell_dir: Path):
    resp_path = cell_dir / "lt.exe.resp"
    p_path = cell_dir / "lt.exe.p"
    csv_path = cell_dir / "lte_results.csv"
    if not (resp_path.is_file() and p_path.is_file() and csv_path.is_file()):
        return None

    resp = L.read_resp(resp_path)
    v = resp.get("VALIDATE", "FALSE").strip().upper()
    if v not in ("TRUE", "T", "1"):
        return None  # only auditing VALIDATE=TRUE cells

    cell = cell_dir.name
    resp.setdefault("CLIMATE_INDEX", str(cell_dir / f"{cell}.dat"))
    params = json.loads(p_path.read_text())

    cfg = L.Config(resp, {})
    try:
        r = L.forward(cfg, params)
    except Exception as e:  # noqa: BLE001
        return dict(cell=cell, error=str(e))

    exp = np.loadtxt(csv_path, delimiter=",")
    got_forcing = r["forcing"]
    if exp.shape[0] != len(got_forcing):
        return dict(cell=cell, error=f"row count mismatch "
                                      f"{exp.shape[0]} vs {len(got_forcing)}")

    exp_f = exp[:, 3]
    d = exp_f - got_forcing
    max_abs = float(np.max(np.abs(d)))
    scale = max(float(np.max(np.abs(exp_f))), 1e-30)
    rel = max_abs / scale
    corr = float(np.corrcoef(exp_f, got_forcing)[0, 1]) if np.std(got_forcing) > 0 else 1.0
    return dict(cell=cell, max_abs=max_abs, rel=rel, corr=corr)


def main():
    threshold = 1e-3
    if "--threshold" in sys.argv:
        threshold = float(sys.argv[sys.argv.index("--threshold") + 1])

    results = []
    seen = set()
    for d in sorted(HERE.glob("k*")):
        if not d.is_dir():
            continue
        # Sep2026 dirs with VALIDATE=TRUE in resp but no own lt.exe.p/
        # lte_results.csv yet keep their real, finalized result in the
        # Feb2026 companion dir instead (same pattern as
        # plot_hidden_latent_forcing_all_quads.py's freshest_results_csv).
        candidate = d
        if not ((d / "lt.exe.p").is_file() and (d / "lte_results.csv").is_file()):
            alt = FEB2026 / d.name
            if (alt / "lt.exe.p").is_file() and (alt / "lte_results.csv").is_file():
                candidate = alt
        res = audit_cell(candidate)
        if res is not None:
            results.append(res)
            seen.add(d.name)

    ok, suspect, errored = [], [], []
    for res in results:
        if "error" in res:
            errored.append(res)
        elif res["rel"] > threshold or res["corr"] < 0.999:
            suspect.append(res)
        else:
            ok.append(res)

    print(f"Audited {len(results)} VALIDATE=TRUE cells "
          f"(threshold rel={threshold})\n")

    print(f"=== OK ({len(ok)}) -- Forcing matches saved lt.exe.p ===")
    for res in sorted(ok, key=lambda r: -r["rel"]):
        print(f"  {res['cell']:<14} rel={res['rel']:.2e}  corr={res['corr']:.9f}")

    print(f"\n=== SUSPECT ({len(suspect)}) -- likely hit by the resync bug ===")
    for res in sorted(suspect, key=lambda r: -r["rel"]):
        print(f"  {res['cell']:<14} rel={res['rel']:.2e}  corr={res['corr']:.9f}")

    if errored:
        print(f"\n=== ERRORED ({len(errored)}) -- could not replay ===")
        for res in errored:
            print(f"  {res['cell']:<14} {res['error']}")

    out = HERE / "audit_validate_bug_results.json"
    out.write_text(json.dumps(dict(ok=ok, suspect=suspect, errored=errored),
                               indent=2))
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
