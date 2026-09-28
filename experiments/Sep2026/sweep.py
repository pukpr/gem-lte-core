#!/usr/bin/env python3
"""Pilot: iterative flood-fill sweep of climate-index "quads" (Kaplan SST
grid cells) using the real lt.exe stochastic optimizer, seeded
geographically from whichever already-good result (a flagship, or another
already-solved cell) is nearest -- see /home/paul/.claude/plans/
twinkly-yawning-candy.md for the full design.

Fully isolated in experiments/Sep2026/ -- never writes into experiments/
Feb2026/. Flagship directories there are read-only donor sources.

Usage:
    python3 sweep.py --pilot [--n-cells 10] [--timeout 300] [--force]
"""
from __future__ import annotations

import argparse
import json
import math
import re
import shlex
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
FEB = HERE.parent / "Feb2026"
LEDGER_PATH = HERE / "sweep_ledger.json"
WATCH_EXP = HERE / "watch_run.exp"

FLAGSHIPS = ["nino4", "amo", "pdo", "baltic", "sam",
             "nino34", "nao", "tna", "iode", "iodw"]
# nino4/amo/pdo from k_sst_seeded_fit.py's own DONOR_COORDS; baltic/sam
# are best-guess placements (no coordinate exists on record for them) --
# confirmed with the user as acceptable for this pilot. nino34/nao/tna/
# iode/iodw added later (already mature, independently-fitted indices --
# 12-92 backup snapshots each -- used here purely as ADDITIONAL donor
# nodes, never re-fit themselves): nino34 is the standard Nino3.4 box
# center; nao/tna are basin-representative points (NAO has no true single
# point -- a station-difference index -- so this is a best-guess anchor
# for spatial donor-proximity only); iode/iodw are the Indian Ocean
# Dipole's east/west poles -- the first flagships in that basin at all,
# which the original 6 left completely uncovered.
#
# "pna" REMOVED 2026-09-27 (was (50.0, -165.0), also a best-guess
# placement with no real coordinate on record) -- the user flagged it as
# a bad donor for Pacific cells after it produced a failed seed for
# kN060_W170, and real, manually-verified Pacific grid cells now exist
# to seed from instead: kN000_W170, kN020_W170, kN040_W170 (see
# HANDOFF.md's gold-standard section) are real fitted quads at the same
# longitude (W170), not a guess, and should be preferred for any Pacific
# donor lookup. They aren't flagships (they're regular grid cells with
# `manual-gold-standard` ledger entries), so `current_solved_nodes()`
# already picks them up automatically without needing a FLAGSHIPS entry.
FLAGSHIP_COORDS = {
    "nino4": (0.0, -175.0),
    "amo": (40.0, -40.0),
    "pdo": (40.0, -180.0),
    "baltic": (58.0, 20.0),
    "sam": (-60.0, 0.0),
    "nino34": (0.0, -145.0),
    "nao": (45.0, -30.0),
    "tna": (12.0, -27.0),
    "iode": (-5.0, 100.0),
    "iodw": (-8.0, 60.0),
}

GRID_RE = re.compile(r"^k([NS])(\d+)_([EW])(\d+)$")

# Base env, matching the convention already established across every
# generated lte_run.sh in this project (see e.g. Feb2026/kN000_W170/
# lte_run.sh, Feb2026/kS020_E170/lte_run.sh) -- METRIC and TIMEOUT are
# overridden per attempt.
BASE_OVERRIDES = {
    "DLOD_REF": "TRUE",
    "EXCLUDE": "true",
    "F9": "1",
    "IDATE": "1920.9",
    "LOCKA": "FALSE",
    "LOCKT": "FALSE",
    "TEST_ONLY": "false",
    "TRAIN_END": "2005",
    "TRAIN_START": "2000",
    "TREND": "true",
    "ZONE": "FALSE",
    "NUMBER_OF_PROCESSORS": "2",  # REVERTED to 2 (was dropped to 1 out of
                                  # caution about a memory-kill issue that
                                  # turned out to be caused by
                                  # `ulimit -s unlimited`, already fixed
                                  # separately -- see write_lte_run_sh).
                                  # N=1 was very likely why the backbone
                                  # refix pass had a near-0% success rate:
                                  # a single random-descent thread has far
                                  # less chance to escape a bad-alias local
                                  # optimum within TIMEOUT than 2+ parallel
                                  # threads exploring independently.
}

DLOD_FLOOR = 0.994  # restored: the earlier ~0.968 plateau (kN000_W150,
                     # kN060_E010) was traced to wrong-donor-seeding
                     # (verbatim init copy from a donor whose own record
                     # starts on a different date), not a genuine
                     # convergence/timeout limit -- 0.994 is comfortably
                     # achievable once seeding derives the correct value.


# ---------------------------------------------------------------------
# Geo helpers
# ---------------------------------------------------------------------

def parse_grid_name(name: str) -> tuple[float, float] | None:
    m = GRID_RE.match(name)
    if not m:
        return None
    ns, lat, ew, lon = m.groups()
    lat_v = float(lat) * (1.0 if ns == "N" else -1.0)
    lon_v = float(lon) * (1.0 if ew == "E" else -1.0)
    return lat_v, lon_v


def haversine_km(lat1, lon1, lat2, lon2) -> float:
    R = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlmb = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlmb / 2) ** 2
    return 2 * R * math.asin(math.sqrt(a))


# ---------------------------------------------------------------------
# Ledger (this sweep's own -- separate from Feb2026/checkpoints/ledger.json)
# ---------------------------------------------------------------------

def load_sweep_ledger() -> list[dict]:
    if LEDGER_PATH.exists():
        return json.loads(LEDGER_PATH.read_text())
    return []


def append_ledger(entry: dict) -> None:
    entries = load_sweep_ledger()
    entries.append(entry)
    LEDGER_PATH.write_text(json.dumps(entries, indent=2))


def latest_entry_per_cell(entries: list[dict]) -> dict[str, dict]:
    """Last entry per cell wins -- needed so an explicit invalidation
    entry (accepted=False, reason='invalidated-...') correctly overrides
    an earlier accepted=True entry for the same cell, rather than that
    stale acceptance leaking back in because SOME entry for the cell has
    accepted=True somewhere in history."""
    by_cell: dict[str, dict] = {}
    for e in entries:
        by_cell[e["cell"]] = e
    return by_cell


def solved_cells(entries: list[dict]) -> dict[str, dict]:
    return {c: e for c, e in latest_entry_per_cell(entries).items() if e.get("accepted")}


# ---------------------------------------------------------------------
# Candidate selection (from Feb2026's EXISTING one-shot-regression ledger)
# ---------------------------------------------------------------------

def load_poor_cells(threshold: float = 0.0) -> list[dict]:
    ledger = json.loads((FEB / "checkpoints" / "ledger.json").read_text())
    poor = [e for e in ledger if e["holdout_r"] < threshold]
    poor.sort(key=lambda e: e["holdout_r"])
    return poor


def current_solved_nodes() -> dict[str, tuple[float, float]]:
    """Flagships plus every cell this sweep has already solved -- the
    live donor pool for flood-fill proximity, not just the 6 flagships.
    Grows every time a new cell is accepted."""
    nodes = dict(FLAGSHIP_COORDS)
    for cell, e in latest_entry_per_cell(load_sweep_ledger()).items():
        if e.get("accepted"):
            coord = parse_grid_name(cell)
            if coord:
                nodes[cell] = coord
    return nodes


def pick_pilot_cells(n: int, exclude_attempted: bool = True) -> list[str]:
    """Poor-holdout cells nearest to the CURRENT solved set (flagships +
    anything already solved this sweep) -- proper outward flood-fill
    growth, not just proximity to the original 6 flagships."""
    poor = load_poor_cells()
    nodes = current_solved_nodes()
    attempted = {e["cell"] for e in load_sweep_ledger()} if exclude_attempted else set()
    scored = []
    for e in poor:
        c = e["region"]
        if c in attempted:
            continue
        coord = parse_grid_name(c)
        if coord is None:
            continue
        lat, lon = coord
        dist = min(haversine_km(lat, lon, *co) for co in nodes.values())
        scored.append((dist, c))
    scored.sort()
    return [name for _, name in scored[:n]]


# ---------------------------------------------------------------------
# resp editing -- VALIDATE must be set IN the resp (GEM.Getenv precedence
# is command-line > resp > env var > default, so an env override of
# VALIDATE is silently ignored whenever a resp already sets it actively,
# as several grid cells' resp files do). METRIC is safe via env override
# -- confirmed no grid-cell resp mentions METRIC at all.
# ---------------------------------------------------------------------

def set_resp_key(resp_path: Path, key: str, value: str) -> None:
    lines = resp_path.read_text().splitlines()
    pattern = re.compile(rf"^-?{re.escape(key)}\b")
    out = [ln for ln in lines if not pattern.match(ln.strip())]
    out.append(f"{key:<30s} {value}")
    resp_path.write_text("\n".join(out) + "\n")


def write_lte_run_sh(run_dir: Path, overrides: dict[str, str]) -> None:
    """Self-contained -- deliberately does NOT import lte_gui.py (per
    user direction: lte_run.sh already carries everything needed for
    non-interactive use). Matches the exact format lte_gui.py's own
    write_run_scripts produces, for consistency with every other
    generated lte_run.sh in this project."""
    lines = [
        "#!/usr/bin/env bash",
        "# Generated by sweep.py for one optimization attempt.",
        'cd "$(dirname "$0")"',
    ]
    for key in sorted(overrides):
        lines.append(f"export {key}={shlex.quote(overrides[key])}")
    # 64MB, not "unlimited" -- every project lte_run.sh uses "unlimited",
    # but every sweep.py attempt using it this session got killed by the
    # harness's own memory-safety monitor within seconds, before any
    # output at all; a direct side-by-side test with a bounded stack ran
    # fine (modest ~58MB RSS) in the exact same system state. 64MB is far
    # more than any observed real usage, just not "request everything".
    lines.append("ulimit -s 65536")
    lines.append("../lt.exe -j")
    sh_path = run_dir / "lte_run.sh"
    sh_path.write_text("\n".join(lines) + "\n")
    sh_path.chmod(0o755)


# ---------------------------------------------------------------------
# Per-cell working directory setup
# ---------------------------------------------------------------------

def cell_source_dir(cell: str) -> Path:
    return FEB / cell


def data_start_date(dat_path: Path) -> float:
    with dat_path.open() as f:
        first_line = f.readline()
    return float(first_line.split()[0])


def read_forcing_at(csv_path: Path, date: float, tol: float = 0.001) -> float | None:
    """Column 4 (post-Bessel Forcing) of an lte_results.csv at a given date."""
    import csv as csv_mod
    with csv_path.open(newline="") as f:
        for row in csv_mod.reader(f):
            if not row:
                continue
            try:
                d = float(row[0])
            except ValueError:
                continue
            if abs(d - date) < tol:
                return float(row[3])
    return None


def invert_bessel(target_post_bessel: float, impA: float, impB: float,
                   offs: float, bg: float, k: float, guess: float = 0.0) -> float:
    """Newton's-method inverse of Ada's Bessel(v, impA, impB, k, offs, bg)
    -- see gem-lte-primitives.adb. Confirmed exact (residual < 1e-9) every
    time it's been checked against a real Ada-computed CSV value this
    session."""
    def bessel(v):
        return (v + impA * math.sin(2 * math.pi * k * v) + impB * math.cos(2 * math.pi * k * v)
                + offs * impA * math.sin(4 * math.pi * k * v) + bg * impB * math.cos(4 * math.pi * k * v))

    def dbessel(v, h=1e-6):
        return (bessel(v + h) - bessel(v - h)) / (2 * h)

    v = guess if guess else target_post_bessel
    for _ in range(60):
        d = dbessel(v)
        if d == 0:
            break
        v_new = v - (bessel(v) - target_post_bessel) / d
        if abs(v_new - v) < 1e-10:
            v = v_new
            break
        v = v_new
    return v


def derive_seed_init(donor_dir: Path, target_start: float) -> float | None:
    """The donor's own value of `init` AT target_start, derived from the
    donor's own real (Ada-computed) lte_results.csv by inverting the
    Bessel nonlinearity -- NOT a verbatim copy of the donor's own `init`
    field, which means "value at the DONOR's own start" and is wrong
    (often catastrophically -- confirmed this session: seeding
    kN000_W150 with nino4's raw init sent dLOD to -0.48; the derived
    value at 1950 gave dLOD~0.995) whenever the donor's own record starts
    on a different date than the target's. Returns None if the donor's
    own CSV doesn't reach back to target_start (nothing to derive from)."""
    donor_p = json.loads((donor_dir / "lt.exe.p").read_text())
    donor_resp = read_resp_dict(donor_dir / "lt.exe.resp")
    nm = int(float(donor_resp.get("NM", len(donor_p["ltep"]))))
    k = donor_p["ltep"][nm - 1]
    csv_path = donor_dir / "lte_results.csv"
    if not csv_path.exists():
        return None
    post_bessel = read_forcing_at(csv_path, target_start)
    if post_bessel is None:
        return None
    return invert_bessel(post_bessel, donor_p["impA"], donor_p["impB"],
                          donor_p["offs"], donor_p["bg"], k,
                          guess=donor_p["init"])


def read_resp_dict(resp_path: Path) -> dict[str, str]:
    out = {}
    for line in resp_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("-"):
            continue
        parts = line.split(None, 1)
        if len(parts) == 2:
            out[parts[0]] = parts[1].strip().strip('"')
    return out


def clear_stale_checkpoint(target: Path, cell: str) -> None:
    """lt.exe writes its own per-CLIMATE_INDEX checkpoint file,
    `lt.exe.<cell>.dat.p`, alongside the generic `lt.exe.p` at the end of
    every run -- and, per the user directly (2026-09-27), PREFERS loading
    that per-index file over `lt.exe.p` on its next run if both exist.
    That means every setup_cell_dir* call that writes a fresh donor seed
    to `lt.exe.p` was being silently ignored on any cell that had already
    been attempted once before (which left its own per-index checkpoint
    behind) -- the run would just reload its own prior state instead of
    the intended fresh seed. This is almost certainly why the 24-cell
    backbone-refix backlog showed ~0% success across BOTH N=1 and N=2:
    only each cell's very first-ever attempt could have actually used the
    intended donor seed; every retry after that silently no-opped. Must
    be called every time `lt.exe.p` is freshly written for a cell."""
    stale = target / f"lt.exe.{cell}.dat.p"
    if stale.exists():
        stale.unlink()


def setup_cell_dir(cell: str, donor_dir: Path) -> Path:
    """Copy the cell's own data/resp template from Feb2026 (first time
    only) and seed lt.exe.p from the donor. If the donor's own record
    starts on (approximately) the same date as this cell's real data,
    its `init` is already valid there -- wholesale copy, the kS040_W050
    precedent. Otherwise derive the donor's own value AT this cell's own
    start date (see derive_seed_init) -- confirmed this session to be
    the difference between a catastrophic fit (dLOD ~ -0.48) and a good
    one (dLOD ~ 0.995), far outweighing METRIC/VALIDATE choice. IDATE is
    then set explicitly to the cell's own real start too, rather than
    relying on IIR's own implicit clamping when a stale/inherited IDATE
    from a shared template precedes the real data (every grid cell's
    resp defaults to IDATE=1880 even though the data itself starts 1950)."""
    target = HERE / cell
    target.mkdir(exist_ok=True)
    src = cell_source_dir(cell)
    dat_name = f"{cell}.dat"
    for fname in (dat_name, "dlod_ref.dat", "lt.exe.resp"):
        src_f = src / fname
        if src_f.exists() and not (target / fname).exists():
            (target / fname).write_bytes(src_f.read_bytes())

    target_start = data_start_date(target / dat_name)
    donor_dat = donor_dir / f"{donor_dir.name}.dat"
    donor_start = data_start_date(donor_dat) if donor_dat.exists() else None

    p = json.loads((donor_dir / "lt.exe.p").read_text())
    if donor_start is not None and abs(donor_start - target_start) > 0.5:
        seed = derive_seed_init(donor_dir, target_start)
        if seed is not None:
            p = dict(p)
            p["init"] = seed
    clear_stale_checkpoint(target, cell)
    (target / "lt.exe.p").write_text(json.dumps(p, indent=2))
    set_resp_key(target / "lt.exe.resp", "IDATE", f"{target_start:.4f}")
    return target


# ---------------------------------------------------------------------
# Triangulated (multi-donor, inverse-distance-weighted) seeding -- for
# stragglers a single nearest donor couldn't resolve, or a result whose
# scores look like a fluke of one particular donor's own quirks.
# ---------------------------------------------------------------------

SCALAR_FIELDS = ["offs", "bg", "impA", "impB", "impC", "delA", "delB", "asym",
                 "ann1", "ann2", "IR", "ma", "mp", "shfT"]


def find_k_nearest_solved(cell: str, k: int = 3, exclude: set[str] = frozenset()
                           ) -> list[tuple[str, float]]:
    lat, lon = parse_grid_name(cell)
    nodes = current_solved_nodes()
    dists = sorted(((n, haversine_km(lat, lon, *co)) for n, co in nodes.items()
                     if n not in exclude),
                    key=lambda x: x[1])
    return dists[:k]


def donor_dir_for(name: str) -> Path:
    return FEB / name if name in FLAGSHIPS else HERE / name


def triangulate_seed_p(cell: str, target_start: float, k: int = 3,
                        exclude: set[str] = frozenset()
                        ) -> tuple[dict, list[str]]:
    """Inverse-distance-weighted (IDW, weight=1/d^2) blend of the K
    nearest solved donors' shape parameters, rather than a single
    nearest-neighbor copy -- a real kriging-style triangulation. Each
    donor's OWN `init` is first translated to "value at target_start" via
    derive_seed_init (same Bessel-inversion as the single-donor path)
    before being blended, so the average is physically meaningful rather
    than averaging two numbers anchored at different dates. `lpap`/`ltep`/
    `harm` are blended element-wise only when every donor shares the same
    array length; otherwise falls back to the nearest donor's own arrays
    (blending across genuinely different harmonic *schemes* would produce
    a meaningless result, unlike blending nearby values of the same
    quantity)."""
    neighbors = find_k_nearest_solved(cell, k, exclude=exclude)
    weights = [1.0 / max(d, 1.0) ** 2 for _, d in neighbors]
    wsum = sum(weights)
    weights = [w / wsum for w in weights]

    donor_ps, donor_inits = [], []
    for name, _dist in neighbors:
        ddir = donor_dir_for(name)
        p = json.loads((ddir / "lt.exe.p").read_text())
        donor_ps.append(p)
        donor_dat = ddir / f"{ddir.name}.dat"
        d_start = data_start_date(donor_dat) if donor_dat.exists() else None
        if d_start is not None and abs(d_start - target_start) > 0.5:
            init_val = derive_seed_init(ddir, target_start)
        else:
            init_val = p["init"]
        donor_inits.append(init_val if init_val is not None else p["init"])

    blended = {f: sum(p.get(f, 0.0) * w for p, w in zip(donor_ps, weights))
               for f in SCALAR_FIELDS}
    blended["init"] = sum(v * w for v, w in zip(donor_inits, weights))

    base = donor_ps[0]
    if len({len(p["ltep"]) for p in donor_ps}) == 1:
        n = len(base["ltep"])
        blended["ltep"] = [sum(p["ltep"][i] * w for p, w in zip(donor_ps, weights))
                            for i in range(n)]
    else:
        blended["ltep"] = base["ltep"]

    harm_lens = {len(p.get("harm", [])) for p in donor_ps}
    if len(harm_lens) == 1 and harm_lens != {0}:
        n = harm_lens.pop()
        blended["harm"] = [sum(p["harm"][i] * w for p, w in zip(donor_ps, weights))
                            for i in range(n)]
    else:
        blended["harm"] = base.get("harm", [])

    if len({len(p["lpap"]) for p in donor_ps}) == 1:
        n = len(base["lpap"])
        blended["lpap"] = [[sum(p["lpap"][i][j] * w for p, w in zip(donor_ps, weights))
                             for j in range(3)] for i in range(n)]
    else:
        blended["lpap"] = base["lpap"]

    return blended, [n for n, _ in neighbors]


def setup_cell_dir_triangulated(cell: str, k: int = 3,
                                 exclude: set[str] = frozenset()
                                 ) -> tuple[Path, list[str]]:
    target = HERE / cell
    target.mkdir(exist_ok=True)
    src = cell_source_dir(cell)
    dat_name = f"{cell}.dat"
    for fname in (dat_name, "dlod_ref.dat", "lt.exe.resp"):
        src_f = src / fname
        if src_f.exists() and not (target / fname).exists():
            (target / fname).write_bytes(src_f.read_bytes())
    target_start = data_start_date(target / dat_name)
    blended, donor_names = triangulate_seed_p(cell, target_start, k=k, exclude=exclude)
    clear_stale_checkpoint(target, cell)
    (target / "lt.exe.p").write_text(json.dumps(blended, indent=2))
    set_resp_key(target / "lt.exe.resp", "IDATE", f"{target_start:.4f}")
    return target, donor_names


# ---------------------------------------------------------------------
# Watcher invocation + output parsing
# ---------------------------------------------------------------------

TRIPLET_RE = re.compile(
    r"(CC|DTW)\s+([-\d.]+)\s+V:\s+([-\d.]+)\s+([-\d.]+)\s+(\d+)\s+(\d+)")
PAIR_RE = re.compile(r"(CC|DTW)\s+([-\d.]+)\s+([-\d.]+)\s+(\d+)\s+(\d+)")
DLOD_RE = re.compile(r"([-\d.]+):dLOD:")


def run_attempt(run_dir: Path, internal_timeout_s: int) -> dict:
    """Invoke watch_run.exp, collect every tagged line, and boil it down
    to: the best (highest-train) triplet/pair seen, the best (highest)
    dLOD seen, and whether any thread reported a deadlock.

    internal_timeout_s is lt.exe's OWN TIMEOUT env var (how long its
    search loop runs before Report_Final fires). The Expect watcher's
    OWN external patience must be meaningfully LONGER than that -- Ada's
    timeout check is periodic, not instantaneous, and Report_Final itself
    (writing CSV, printing several lines per thread) takes real wall-clock
    time on top of it. Confirmed by direct test: an external timeout equal
    to the internal one raced Report_Final and lost, seeing nothing."""
    watcher_timeout_s = internal_timeout_s + 90
    proc = subprocess.run(
        ["expect", "-f", str(WATCH_EXP), str(run_dir), str(watcher_timeout_s)],
        capture_output=True, text=True, timeout=watcher_timeout_s + 30,
    )
    out = proc.stdout
    triplets, pairs, dlods = [], [], []
    deadlocked = False
    for line in out.splitlines():
        if line.startswith("MATCHED_TRIPLET:"):
            m = TRIPLET_RE.search(line)
            if m:
                triplets.append({
                    "metric": m.group(1), "train": float(m.group(2)),
                    "validate": float(m.group(3)), "test": float(m.group(4)),
                })
        elif line.startswith("MATCHED_PAIR:"):
            m = PAIR_RE.search(line)
            if m:
                pairs.append({
                    "metric": m.group(1), "val1": float(m.group(2)),
                    "val2": float(m.group(3)),
                })
        elif line.startswith("MATCHED_DEADLOCK:"):
            deadlocked = True
        elif line.startswith("MATCHED_DLOD:"):
            m = DLOD_RE.search(line)
            if m:
                dlods.append(float(m.group(1)))
    best_triplet = max(triplets, key=lambda t: t["train"]) if triplets else None
    best_pair = max(pairs, key=lambda p: p["val1"]) if pairs else None
    best_dlod = max(dlods) if dlods else None
    return dict(deadlocked=deadlocked, triplet=best_triplet, pair=best_pair,
                dlod=best_dlod, raw_tail=out[-2000:])


# ---------------------------------------------------------------------
# Per-cell attempt state machine
# ---------------------------------------------------------------------

def stiff(result: dict) -> bool:
    """'Too stiff' = deadlocked, or a VALIDATE=TRUE triplet whose train
    score is fine but validate collapses (large train-validate gap)."""
    if result["deadlocked"]:
        return True
    t = result["triplet"]
    if t is not None and t["train"] > 0.3 and t["validate"] < 0.05:
        return True
    return False


def gate_ok(result: dict) -> bool:
    if result["dlod"] is None or result["dlod"] < DLOD_FLOOR:
        return False
    return result["triplet"] is not None or result["pair"] is not None


# Empirically-grounded range: every confirmed-good fit across this whole
# project (nino4, pdo, every sane grid cell) has its backbone winding
# (ltep[nm-1]) somewhere in ~0.2-1.0 (the shared tidal backbone and its
# first couple of dyadic doublings); every confirmed-bad one found by
# comparing kS020_E050's manual "gold standard" against its automated
# neighbors landed at >=10 (18.68, 100, 150, 200, 600, 678 recurring
# identically across unrelated cells) -- a clean bimodal split, no
# borderline cases observed, so a wide margin is safe.
BACKBONE_MIN, BACKBONE_MAX = 0.1, 2.5


def read_backbone(cell_dir: Path) -> float | None:
    try:
        p = json.loads((cell_dir / "lt.exe.p").read_text())
        resp = read_resp_dict(cell_dir / "lt.exe.resp")
        nm = int(float(resp.get("NM", len(p["ltep"]))))
        return p["ltep"][nm - 1]
    except (FileNotFoundError, KeyError, IndexError, ValueError):
        return None


def backbone_ok(cell_dir: Path) -> bool:
    """Catches the failure mode dLOD/CC alone miss entirely: a real
    optimizer run, given a wide SPREAD_MAX, occasionally drifts a
    perfectly good seed to a spurious high-frequency alias that still
    scores well on CC/dLOD by coincidence -- confirmed directly this
    session (kS040_E090, backbone 0.4455, correctly seeded kS040_E070,
    which drifted to 18.68 on its own and then propagated that bad value
    to kS040_E050 via flood-fill donor inheritance)."""
    bb = read_backbone(cell_dir)
    return bb is not None and BACKBONE_MIN <= abs(bb) <= BACKBONE_MAX


def attempt(cell_dir: Path, validate: bool, metric: str, timeout_s: int) -> dict:
    set_resp_key(cell_dir / "lt.exe.resp", "VALIDATE", "TRUE" if validate else "FALSE")
    overrides = dict(BASE_OVERRIDES)
    overrides["CLIMATE_INDEX"] = f"{cell_dir.name}.dat"
    overrides["METRIC"] = metric
    overrides["TIMEOUT"] = str(timeout_s)
    write_lte_run_sh(cell_dir, overrides)
    # Each of the 4 cascade attempts (DTW/CC x VALIDATE TRUE/FALSE) must
    # start from the SAME donor-derived lt.exe.p, not from whatever
    # per-CLIMATE_INDEX checkpoint the PREVIOUS attempt in this same
    # cascade just wrote -- see clear_stale_checkpoint's docstring. Without
    # this, attempts 2-4 silently continue from attempt 1's ending state
    # instead of being independent fresh trials.
    clear_stale_checkpoint(cell_dir, cell_dir.name)
    return run_attempt(cell_dir, timeout_s)


def solve_cell(cell: str, donor_dir: Path, timeout_s: int) -> dict:
    """Attempt order revised from the original CC-first design based on
    direct evidence this session (kN000_W150, both un-seeded historical
    logs and fresh nino4-seeded runs): DTW consistently generalizes
    better than CC for these short-record cells -- CC's held-out score
    collapsed (0.26, 0.155) even with good train/validate, while DTW
    stayed consistent (0.55-0.68) across all three splits every time it
    was tried. DTW now goes first; CC is the fallback, not the reverse."""
    cell_dir = setup_cell_dir(cell, donor_dir)
    return _run_cascade(cell, cell_dir, timeout_s)


def solve_cell_triangulated(cell: str, timeout_s: int, k: int = 3,
                             exclude: set[str] = frozenset()) -> dict:
    """Same attempt cascade, but seeded from an inverse-distance-weighted
    blend of the K nearest solved donors instead of a single one -- for
    stragglers a lone donor couldn't resolve, or a result whose scores
    look like a fluke of one particular donor's own quirks."""
    cell_dir, donor_names = setup_cell_dir_triangulated(cell, k=k, exclude=exclude)
    result = _run_cascade(cell, cell_dir, timeout_s)
    result["triangulated_from"] = donor_names
    return result


def _run_cascade(cell: str, cell_dir: Path, timeout_s: int) -> dict:
    log = []

    def record(tag, res):
        log.append((tag, res))

    # Attempt 1: VALIDATE=TRUE, METRIC=DTW
    r1 = attempt(cell_dir, validate=True, metric="DTW", timeout_s=timeout_s)
    record("validate_true_dtw", r1)
    if not stiff(r1) and gate_ok(r1) and backbone_ok(cell_dir):
        return finalize(cell, "validate_true_dtw", r1, log)

    # Attempt 2: VALIDATE=FALSE, METRIC=DTW
    r2 = attempt(cell_dir, validate=False, metric="DTW", timeout_s=timeout_s)
    record("validate_false_dtw", r2)
    if not stiff(r2) and gate_ok(r2) and backbone_ok(cell_dir):
        return finalize(cell, "validate_false_dtw", r2, log)

    # Attempt 3: VALIDATE=TRUE, METRIC=CC -- fallback if DTW itself is stiff
    r3 = attempt(cell_dir, validate=True, metric="CC", timeout_s=timeout_s)
    record("validate_true_cc", r3)
    if not stiff(r3) and gate_ok(r3) and backbone_ok(cell_dir):
        return finalize(cell, "validate_true_cc", r3, log)

    # Attempt 4 (final fallback): VALIDATE=FALSE, METRIC=CC
    r4 = attempt(cell_dir, validate=False, metric="CC", timeout_s=timeout_s)
    record("validate_false_cc", r4)
    if gate_ok(r4) and backbone_ok(cell_dir):
        return finalize(cell, "validate_false_cc", r4, log)

    return dict(cell=cell, accepted=False, reason="unsolved-needs-review",
                attempts=[{"tag": t, "deadlocked": r["deadlocked"],
                           "dlod": r["dlod"]} for t, r in log])


def solve_cell_with_retry(cell: str, donor_dir: Path, timeout_s: int,
                           retries: int = 1) -> dict:
    """A full re-attempt (fresh random search from the same seed) is cheap
    relative to giving up -- confirmed directly this session: kN060_E010
    failed all 4 attempts once (stuck at dLOD~0.92, a bad-luck local
    optimum) then solved cleanly on a completely fresh call with the
    identical method and donor. Per the user's own guidance, the goal is
    a good comprehensive set, not a perfect run on every single cell --
    retry a bounded number of times, then move on and leave it for later
    triangulation/kriging from whichever of its neighbors DID solve."""
    result = None
    for i in range(retries + 1):
        result = solve_cell(cell, donor_dir, timeout_s)
        result["retry_index"] = i
        if result.get("accepted"):
            return result
    return result


def solve_cell_triangulated_with_retry(cell: str, timeout_s: int, k: int = 3,
                                        retries: int = 1,
                                        exclude: set[str] = frozenset()) -> dict:
    # Always exclude the cell itself -- without this, a cell that
    # already has an accepted ledger entry (e.g. a refit attempt on an
    # already-solved cell) can select ITSELF as a triangulation donor,
    # since find_k_nearest_solved/current_solved_nodes has no other way
    # to know it shouldn't. Confirmed directly this session: refitting
    # kN020_E150 without this exclusion silently included kN020_E150 as
    # one of its own 3 donors and made the fit worse, not better.
    exclude = set(exclude) | {cell}
    result = None
    for i in range(retries + 1):
        result = solve_cell_triangulated(cell, timeout_s, k=k, exclude=exclude)
        result["retry_index"] = i
        if result.get("accepted"):
            return result
    return result


def finalize(cell: str, tag: str, result: dict, log) -> dict:
    scores = result["triplet"] or result["pair"]
    return dict(
        cell=cell, accepted=True, attempt_tag=tag,
        validate_mode=(result["triplet"] is not None),
        metric=scores["metric"], scores=scores, dlod=result["dlod"],
        n_attempts=len(log),
    )


# ---------------------------------------------------------------------
# Main pilot loop
# ---------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", action="store_true", required=True)
    ap.add_argument("--n-cells", type=int, default=10)
    ap.add_argument("--timeout", type=int, default=300)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    entries = load_sweep_ledger()
    solved = solved_cells(entries)

    pilot_cells = pick_pilot_cells(args.n_cells)
    print(f"[sweep] pilot cells (poor-holdout, nearest a flagship): {pilot_cells}")

    # solved-set for flood fill: name -> (coords, donor directory)
    solved_nodes = {f: (FLAGSHIP_COORDS[f], FEB / f) for f in FLAGSHIPS}
    for cell, e in solved.items():
        coord = parse_grid_name(cell)
        if coord:
            solved_nodes[cell] = (coord, HERE / cell)

    remaining = [c for c in pilot_cells if c not in solved or args.force]
    while remaining:
        # pick the remaining cell nearest to any currently-solved node
        best_cell, best_dist, best_donor = None, None, None
        for c in remaining:
            coord = parse_grid_name(c)
            if coord is None:
                continue
            lat, lon = coord
            for name, (dc, dp) in solved_nodes.items():
                d = haversine_km(lat, lon, *dc)
                if best_dist is None or d < best_dist:
                    best_cell, best_dist, best_donor = c, d, dp
        if best_cell is None:
            break
        remaining.remove(best_cell)

        print(f"\n[sweep] === {best_cell}  (seed: nearest solved node, "
              f"{best_dist:.0f} km, from {best_donor}) ===")
        t0 = time.time()
        result = solve_cell(best_cell, best_donor, args.timeout)
        result["seed_from_path"] = str(best_donor)
        result["timestamp"] = datetime.now(timezone.utc).isoformat()
        result["elapsed_s"] = round(time.time() - t0, 1)
        append_ledger(result)

        if result.get("accepted"):
            coord = parse_grid_name(best_cell)
            solved_nodes[best_cell] = (coord, HERE / best_cell)
            print(f"[sweep] SOLVED {best_cell}: {result['attempt_tag']} "
                  f"dLOD={result['dlod']:.6f} scores={result['scores']}")
        else:
            print(f"[sweep] UNSOLVED {best_cell}: {result.get('reason')}")

    print("\n[sweep] pilot pass complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
