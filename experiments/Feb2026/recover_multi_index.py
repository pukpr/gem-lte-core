#!/usr/bin/env python3
"""Recover per-index params from multi_index/ and write to each index's subdir.

Combines shared params + per-index params into a full lt.exe.p for each index.

Usage:
    python3 recover_multi_index.py                    # from multi_index/
    python3 recover_multi_index.py --suffix multi     # saves as lt.exe.multi.p
    python3 recover_multi_index.py --overwrite        # overwrites lt.exe.p
"""
from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser(description="Recover multi-index params")
    ap.add_argument("--suffix", default="p", help="output suffix (default: p, overwrites original)")
    ap.add_argument("--overwrite", action="store_true", help="overwrite lt.exe.p directly")
    args = ap.parse_args()

    multi_dir = ROOT / "multi_index"
    if not multi_dir.exists():
        print(f"ERROR: {multi_dir} not found. Run multi-index optimization first.")
        return 2

    shared = json.loads((multi_dir / "shared_params.p").read_text())
    per_index = json.loads((multi_dir / "per_index_params.json").read_text())

    print(f"Loaded shared params from: {multi_dir / 'shared_params.p'}")
    print(f"Loaded per-index params for: {list(per_index.keys())}")

    for idx, pip in per_index.items():
        idx_dir = ROOT / idx
        if not idx_dir.is_dir():
            print(f"  SKIP {idx}: directory not found")
            continue

        # Build full params: shared + per-index overrides
        params = dict(shared)
        params["ltep"] = pip["ltep"]
        params["IR"] = pip["IR"]
        if "harm" in pip:
            params["harm"] = pip["harm"]

        # Determine output file
        if args.overwrite:
            out_path = idx_dir / "lt.exe.p"
            # Backup original
            backup = idx_dir / "lt.exe.p.bak"
            if out_path.exists():
                shutil.copy2(out_path, backup)
                print(f"  Backed up {idx}/lt.exe.p -> lt.exe.p.bak")
        else:
            out_path = idx_dir / f"lt.exe.{args.suffix}.p"

        with open(out_path, "w") as f:
            json.dump(params, f, indent=2)

        # Verify
        orig_path = idx_dir / "lt.exe.p"
        if orig_path.exists():
            orig = json.loads(orig_path.read_text())
            print(f"  {idx}: ltep={params['ltep']} (was {orig.get('ltep', 'N/A')})")
            print(f"       IR={params['IR']:.6f} (was {orig.get('IR', 'N/A'):.6f})")
        else:
            print(f"  {idx}: ltep={params['ltep']}, IR={params['IR']:.6f}")
        print(f"       -> {out_path}")

    print(f"\nDone. Recovered {len(per_index)} index params.")


if __name__ == "__main__":
    main()
