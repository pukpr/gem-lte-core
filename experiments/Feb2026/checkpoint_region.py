#!/usr/bin/env python3
"""checkpoint_region.py <region> <seed> -- snapshot a region's current
lt.exe.p/lt.exe.resp/lte_results.csv into checkpoints/<region>/, compute
its train/holdout R (2000-2005 block), and append to the running ledger
(checkpoints/ledger.json) -- so a known-good, CV'd fit is never lost to
a subsequent rerun, and there's a single source of truth for which
regions are validated and what seeded what.
"""
import json, shutil, sys, datetime
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent
CKPT = ROOT / "checkpoints"

def main():
    region, seed = sys.argv[1], sys.argv[2]
    rdir = ROOT / region
    raw = np.loadtxt(rdir / "lte_results.csv", delimiter=",")
    dates, model, data = raw[:,0], raw[:,1], raw[:,2]
    train_m = (dates < 2000) | (dates >= 2005)
    hold_m = (dates >= 2000) & (dates < 2005)
    train_r = float(np.corrcoef(model[train_m], data[train_m])[0,1])
    hold_r = float(np.corrcoef(model[hold_m], data[hold_m])[0,1])

    dest = CKPT / region
    dest.mkdir(parents=True, exist_ok=True)
    for fn in ["lt.exe.p", "lt.exe.resp", "lte_results.csv"]:
        shutil.copy(rdir / fn, dest / fn)

    ledger_path = CKPT / "ledger.json"
    ledger = json.loads(ledger_path.read_text()) if ledger_path.exists() else []
    ledger = [e for e in ledger if e["region"] != region]  # replace prior entry for this region
    ledger.append(dict(region=region, seed=seed, train_r=train_r, holdout_r=hold_r,
                        timestamp=datetime.datetime.now().isoformat(timespec="seconds")))
    ledger_path.write_text(json.dumps(ledger, indent=1))
    print(f"[checkpoint] {region} (seed={seed}): train_r={train_r:+.4f} holdout_r={hold_r:+.4f} -> saved to {dest}")

if __name__ == "__main__":
    main()
