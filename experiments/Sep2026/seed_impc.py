#!/usr/bin/env python3
"""Seed ImpC=0.01 into a cell's lt.exe.p (2026-09-29, per user request).

gem-random_descent.adb's Markov procedure never perturbs a parameter
sitting at exactly 0.0 (Set(I)=0.0 branch just recurses to pick a
different parameter, unless FLIP<0.0, which no existing resp sets) --
D.B.ImpC has stayed 0.0 in every fit across this whole project not
because it doesn't help, but because the search structurally can never
discover otherwise. Seeding a small nonzero value lets ordinary
relative-perturbation random descent actually explore it (Markov's
elsif branch: Set(I) := Set(I) * (1 + Spread*...), with a small FLIP
probability of a full sign flip each iteration).

Only ever seeds cells whose CURRENT impC is exactly 0.0 -- never
overwrites an already-fitted nonzero value.

IMPORTANT: edits the file as TEXT (a targeted regex substitution of just
the "impC" line), not via json.loads/json.dumps -- Ada's JSON reader
choked on a full round-trip through Python's compact single-line
json.dumps output ("2:9: comma expected"), corrupting two cells'
checkpoints the first time this was tried (2026-09-29, recovered via
`git checkout`). Ada's own writer (GNATCOLL.JSON, Write(Compact=>False))
always produces multi-line, 2-space-indented, trailing-newline-terminated
JSON, and its reader appears to expect that shape back -- preserving the
file byte-for-byte except the one changed line is the only edit
confirmed safe.

Usage: python3 seed_impc.py CELL [CELL ...]
"""
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEED_VALUE = 0.01
IMPC_RE = re.compile(r'("impC"\s*:\s*)([^,\n]+)(,?)')


def seed_one(cell: str) -> str:
    p_path = HERE / cell / "lt.exe.p"
    if not p_path.is_file():
        return f"{cell}: no lt.exe.p, skipped"
    text = p_path.read_text()
    params = json.loads(text)  # read-only: validate + inspect current value
    current = params.get("impC", 0.0)
    if current != 0.0:
        return f"{cell}: impC already {current}, left unchanged"
    new_text, n = IMPC_RE.subn(rf"\g<1>{SEED_VALUE}\g<3>", text, count=1)
    if n != 1:
        return f"{cell}: could not locate a single \"impC\" line, skipped"
    p_path.write_text(new_text)
    return f"{cell}: impC 0.0 -> {SEED_VALUE}"


def main() -> None:
    cells = sys.argv[1:]
    if not cells:
        print("usage: seed_impc.py CELL [CELL ...]", file=sys.stderr)
        sys.exit(2)
    for cell in cells:
        print(seed_one(cell))


if __name__ == "__main__":
    main()
