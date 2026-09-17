#!/usr/bin/env python3
"""Trial tracker for LTE parameter optimization runs.

Records each run with its configuration, starting params, and results.
Enables comparison across indices and seeding from previous best runs.

Usage:
    ./trial_tracker.py record <index> <subset> <cc_before> <cc_after> [--notes "..."]
    ./trial_tracker.py list [--index amo] [--sort cc_gain]
    ./trial_tracker.py best [--index nino4]  # best seed for an index
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TRACKER_FILE = ROOT / "trials.jsonl"


def record(index: str, subset: str, cc_before: float, cc_after: float,
           params_file: str | None = None, notes: str = "",
           seed_from: str = ""):
    """Record a trial run."""
    entry = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "index": index,
        "subset": subset,
        "cc_before": round(cc_before, 6),
        "cc_after": round(cc_after, 6),
        "cc_gain": round(cc_after - cc_before, 6),
        "seed_from": seed_from,
        "params_file": params_file or "",
        "notes": notes,
    }
    with open(TRACKER_FILE, "a") as f:
        f.write(json.dumps(entry) + "\n")
    print(f"Recorded: {index}/{subset} CC {cc_before:.6f} -> {cc_after:.6f} "
          f"(gain {entry['cc_gain']:+.6f})")
    return entry


def list_trials(index: str | None = None, sort_by: str = "timestamp"):
    """List recorded trials."""
    if not TRACKER_FILE.exists():
        print("No trials recorded yet.")
        return []

    trials = []
    for line in TRACKER_FILE.read_text().strip().split("\n"):
        if line.strip():
            trials.append(json.loads(line))

    if index:
        trials = [t for t in trials if t["index"] == index]

    if sort_by == "cc_gain":
        trials.sort(key=lambda t: t["cc_gain"], reverse=True)
    elif sort_by == "cc_after":
        trials.sort(key=lambda t: t["cc_after"], reverse=True)
    else:
        trials.sort(key=lambda t: t["timestamp"])

    print(f"{'Timestamp':<22} {'Index':<8} {'Subset':<12} "
          f"{'CC_before':>10} {'CC_after':>10} {'Gain':>10} {'Seed':<8}")
    print("-" * 90)
    for t in trials:
        print(f"{t['timestamp']:<22} {t['index']:<8} {t['subset']:<12} "
              f"{t['cc_before']:>10.6f} {t['cc_after']:>10.6f} "
              f"{t['cc_gain']:>+10.6f} {t.get('seed_from',''):<8}")

    if trials:
        gains = [t["cc_gain"] for t in trials]
        print(f"\nTotal trials: {len(trials)}")
        print(f"Total gain: {sum(gains):+.6f}")
        print(f"Best gain:  {max(gains):+.6f}")

    return trials


def best_seed(target_index: str):
    """Find the best cumulative CC for a target index.

    Shows the chain of optimizations and the total gain from the original baseline.
    """
    if not TRACKER_FILE.exists():
        print("No trials recorded yet.")
        return None

    trials = []
    for line in TRACKER_FILE.read_text().strip().split("\n"):
        if line.strip():
            trials.append(json.loads(line))

    # Filter to target index
    idx_trials = [t for t in trials if t["index"] == target_index]
    if not idx_trials:
        print(f"No trials for {target_index}.")
        return None

    # Sort by timestamp to show the chain
    idx_trials.sort(key=lambda t: t["timestamp"])

    print(f"\n=== Optimization chain for {target_index} ===")
    cumulative_before = idx_trials[0]["cc_before"]
    print(f"Original baseline: {cumulative_before:.6f}")
    for i, t in enumerate(idx_trials):
        print(f"  Step {i+1}: {t['subset']:12s} {t['cc_before']:.6f} -> "
              f"{t['cc_after']:.6f} (gain {t['cc_gain']:+.6f}) "
              f"[seed: {t.get('seed_from','')}]")
        cumulative_before = t["cc_after"]

    print(f"\nFinal CC: {idx_trials[-1]['cc_after']:.6f}")
    original_baseline = idx_trials[0]["cc_before"]
    # For nino4, the original baseline is the NINO4 own params (0.388)
    # The first trial's cc_before is the PDO-seeded value (0.347)
    # So total gain vs original = final - original_own
    print(f"Gain from first seed: {idx_trials[-1]['cc_after'] - original_baseline:+.6f}")

    return idx_trials[-1]


def main():
    ap = argparse.ArgumentParser(description="LTE trial tracker")
    sub = ap.add_subparsers(dest="command")

    rec = sub.add_parser("record", help="Record a trial run")
    rec.add_argument("index")
    rec.add_argument("subset")
    rec.add_argument("cc_before", type=float)
    rec.add_argument("cc_after", type=float)
    rec.add_argument("--params-file", default="")
    rec.add_argument("--seed-from", default="")
    rec.add_argument("--notes", default="")

    lst = sub.add_parser("list", help="List trials")
    lst.add_argument("--index", default=None)
    lst.add_argument("--sort", default="timestamp",
                     choices=["timestamp", "cc_gain", "cc_after"])

    bst = sub.add_parser("best", help="Find best seed for an index")
    bst.add_argument("index")

    args = ap.parse_args()

    if args.command == "record":
        record(args.index, args.subset, args.cc_before, args.cc_after,
               params_file=args.params_file, seed_from=args.seed_from,
               notes=args.notes)
    elif args.command == "list":
        list_trials(args.index, args.sort)
    elif args.command == "best":
        best_seed(args.index)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
