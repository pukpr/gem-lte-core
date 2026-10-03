#!/usr/bin/env python3
"""Track the COLLECTIVE optimality measure over time (2026-09-28,
per the user's framing: the objective is to move the whole 89-quad
collective toward a better shared parameterization, not to chase every
individual cell's own CC/DTW/dLOD -- a few weak individual cells are
fine if the collective measure is improving; a uniform shared bias in
any parameter (dLOD, year-length, any .p scalar) is fine too, as long as
it's uniform across the collective.

This is a thin wrapper around mlr_shared_windings.py: run it, capture
its overall/mean R^2 numbers, append one line to
collective_r2_log.jsonl with a timestamp and a free-text label for
whatever just happened (a sweep batch, a LOCKW change, a manual fix
round). This log is the actual go/no-go signal for future work --
compare consecutive entries, not per-cell scores, to judge whether a
change helped or hurt.

Usage: python3 collective_r2_log.py "label describing what just happened" [--root DIR]
(--root is passed through to mlr_shared_windings.py and recorded in the log)
"""
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, "/home/paul/eval/gem-lte-core/experiments/Sep2026")
import sweep

LOG_PATH = sweep.HERE / "collective_r2_log.jsonl"


def main() -> None:
    label = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("--") else "unlabeled checkpoint"
    extra = [a for a in sys.argv[1:] if a != label]
    root = extra[extra.index("--root") + 1] if "--root" in extra else None

    proc = subprocess.run(
        [sys.executable, str(sweep.HERE / "mlr_shared_windings.py")] + extra,
        capture_output=True, text=True, cwd=str(sweep.HERE),
    )
    out = proc.stdout

    def grab(prefix: str) -> str | None:
        for line in out.splitlines():
            if line.strip().startswith(prefix):
                return line
        return None

    overall_line = grab("OVERALL R^2")
    indep_line = grab("Per-quad-independent baseline")
    summary_lines = [l for l in out.splitlines() if l.strip().startswith(
        ("Geo-smooth per-quad", "Independent per-quad", "Mean gap", "Parameter count"))]

    entry = dict(
        timestamp=time.strftime("%Y-%m-%dT%H:%M:%S"),
        label=label,
        root=root,
        overall_geo_smooth_r2=overall_line,
        independent_baseline=indep_line,
        summary=summary_lines,
    )

    with open(LOG_PATH, "a") as f:
        f.write(json.dumps(entry) + "\n")

    print(f"Logged to {LOG_PATH}:")
    print(f"  label: {label}")
    print(f"  {overall_line}")
    print(f"  {indep_line}")

    # Show the trend if there's history
    if LOG_PATH.exists():
        lines = LOG_PATH.read_text().splitlines()
        if len(lines) > 1:
            print(f"\nHistory ({len(lines)} entries):")
            for line in lines:
                e = json.loads(line)
                print(f"  {e['timestamp']}  {e['label']:<50} {e.get('overall_geo_smooth_r2','?')}")


if __name__ == "__main__":
    main()
