#!/usr/bin/env python3
"""Refinement sweep: one long CC search per quad from its own current
parameters, kept only if it is better (2026-10-01).

Recipe (the user's kS040_W010 run after the harmonic-save fix): METRIC=CC,
VALIDATE=TRUE, F9=1, the quad's own lt.exe.resp untouched. The user runs
this with 16 threads for ~200 s; the default here is 8 threads x 250 s with 2
quads at a time (same 16 cores, similar thread-seconds per quad).
On Oct2026_altquad/kS040_W010 that took the fit CC from 0.55 to 0.72 with
the manifold still locked. Earlier sweeps (300 s, pre-fix binary) lost any
harmonic change at save, so there is room to gain across the set.

Per quad, in a scratch copy under <root>/_refine/:
  1. baseline: TEST_ONLY score of the current parameters (CC, VALIDATE=TRUE)
  2. search:   lt.exe, METRIC=CC, TIMEOUT s, NUMBER_OF_PROCESSORS=2
  3. re-score: TEST_ONLY of what the search saved
  4. accept if train >= baseline's, validate >= baseline's - 0.01, test
     (the 5-year excluded gap, a noisy CC) >= baseline's - --test-tol,
     dLOD >= floor, backbone in range and (with --ref) manifold CC vs the
     reference >= --min-mcc; then the quad's files are replaced (old ones kept in
     _refine/backup/<quad>/). Otherwise the quad is left untouched.

lt.exe never exits after its final report, and `script` puts it in its own
session, so each run is ended by killing every lt.exe whose cwd is the run
folder (the Oct2026_altquad sweep leaked 118 processes without this).

Resumable (quads already in _refine/results.jsonl are skipped).

Usage:
  python3 sweep_cc_refine.py --root ../Oct2026_altquad \\
      --ref ../Oct2026_altquad/altMedianManifold.dat [--ALL | -set Q ...]
"""
import argparse
import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sweep  # noqa: E402

FILES = ("lt.exe.p", "lt.exe.resp", "lte_run.sh", "lte_results.csv",
         "lt.exe.windings.json", "dlod_compare.csv", "dlod_ref.dat")
LOCK = threading.Lock()


def log(msg):
    with LOCK:
        print(f"{time.strftime('%H:%M:%S')} {msg}", flush=True)


def kill_lt_in(d: Path) -> None:
    for pid in os.listdir("/proc"):
        if not pid.isdigit():
            continue
        try:
            if (os.readlink(f"/proc/{pid}/cwd") == str(d)
                    and open(f"/proc/{pid}/cmdline", "rb").read().split(b"\0")[0].endswith(b"lt.exe")):
                os.kill(int(pid), signal.SIGKILL)
        except OSError:
            continue


def run_lt(d: Path, ov: dict, test_only: bool, timeout: int) -> str:
    sweep.write_lte_run_sh(d, ov)
    sweep.clear_stale_checkpoint(d, d.name)
    res = d / "lte_results.csv"
    res.unlink(missing_ok=True)
    (d / "run.out").unlink(missing_ok=True)
    pr = subprocess.Popen(["bash", "-c", "script -qfc ./lte_run.sh /dev/null > run.out 2>&1 < /dev/null"],
                          cwd=str(d), start_new_session=True)
    t0 = time.time()
    limit = (60 if test_only else timeout) + 300
    while pr.poll() is None and time.time() - t0 < limit:
        time.sleep(3)
        out = (d / "run.out").read_text(errors="ignore") if (d / "run.out").exists() else ""
        settled = res.exists() and time.time() - res.stat().st_mtime > 4
        quiet = time.time() - (d / "run.out").stat().st_mtime > 6 if (d / "run.out").exists() else False
        # one dLOD line is enough (only the saving thread prints it); then wait for quiet output
        if ":dLOD:" in out and settled and quiet and (test_only or time.time() - t0 > timeout * 0.5):
            break
    kill_lt_in(d)
    if pr.poll() is None:
        try:
            os.killpg(pr.pid, signal.SIGKILL)
        except OSError:
            pass
        pr.wait()
    return (d / "run.out").read_text(errors="ignore") if (d / "run.out").exists() else ""


def parse(out: str) -> dict:
    trip, dl = [], []
    for line in out.splitlines():
        m = sweep.TRIPLET_RE.search(line)
        if m:
            trip.append(dict(metric=m.group(1), train=float(m.group(2)),
                             validate=float(m.group(3)), test=float(m.group(4))))
        m = sweep.DLOD_RE.search(line)
        if m:
            dl.append(float(m.group(1)))
    t = trip[-1] if trip else None
    return dict(triplet=t, dlod=(max(dl) if dl else None),
                min=(min(t["train"], t["validate"], t["test"]) if t else None))


def load_results(csv: Path) -> np.ndarray | None:
    """lte_results.csv, retried until complete: in TEST_ONLY both threads
    write it, so a read can land mid-rewrite."""
    n_expected = sum(1 for _ in open(csv.parent / f"{csv.parent.name}.dat"))
    for _ in range(20):
        try:
            a = np.loadtxt(csv, delimiter=",", ndmin=2)
            if a.shape[0] == n_expected and a.shape[1] >= 4:
                return a
        except (OSError, ValueError):
            pass
        time.sleep(0.5)
    return None


def fit_cc(csv: Path) -> float | None:
    a = load_results(csv)
    return None if a is None else float(np.corrcoef(a[:, 1], a[:, 2])[0, 1])


def mcc(csv: Path, ref: np.ndarray | None) -> float | None:
    a = load_results(csv)
    if ref is None or a is None or len(a) != len(ref):
        return None
    return float(np.corrcoef(a[:, 3], ref)[0, 1])


def overrides(d: Path, args, test_only: bool) -> dict:
    ov = dict(sweep.BASE_OVERRIDES)
    ov.update(CLIMATE_INDEX=f"{d.name}.dat", METRIC="CC", F9="1", ACCEL="FALSE",
              NUMBER_OF_PROCESSORS=str(args.threads),
              TEST_ONLY="true" if test_only else "false",
              TIMEOUT="30" if test_only else str(args.timeout))
    return ov


def process(q: str, args, ref, prior: dict | None = None) -> dict:
    """prior: an earlier record whose final re-score failed to report (lt.exe
    aborted mid-save); its search result is still in _refine/<q>/lt.exe.p, so
    only the re-score and the accept decision are redone."""
    root = args.root
    src = root / q
    d = root / "_refine" / q
    t0 = time.time()
    forced = q in set(args.force_quads or [])
    try:
        if prior is None:
            if d.exists():
                shutil.rmtree(d)
            d.mkdir(parents=True)
            for f in ("lt.exe.p", "lt.exe.resp", "dlod_ref.dat"):
                shutil.copy(src / f, d / f)
            shutil.copy(src / f"{q}.dat", d / f"{q}.dat")
            sweep.set_resp_key(d / "lt.exe.resp", "VALIDATE", "TRUE")
            if forced:
                for kv in args.force_set:
                    k, _, v = kv.partition("=")
                    sweep.set_resp_key(d / "lt.exe.resp", k, v)

            base = parse(run_lt(d, overrides(d, args, True), True, 30))
            base_fit = fit_cc(d / "lte_results.csv")
            if forced:   # keep the baseline-with-settings state in case the search is rejected
                (d / "_base").mkdir(exist_ok=True)
                for f in FILES:
                    if (d / f).exists():
                        shutil.copy2(d / f, d / "_base" / f)
            seed_p = (d / "lt.exe.p").read_bytes()
            run_lt(d, overrides(d, args, False), False, args.timeout)
        else:
            base, base_fit, seed_p = prior["base"], prior["base_fit"], b""
        new = parse(run_lt(d, overrides(d, args, True), True, 30))
        new_fit = fit_cc(d / "lte_results.csv")
        m = mcc(d / "lte_results.csv", ref)
        bt, nt = base["triplet"], new["triplet"]
        fails = []
        if bt is None or nt is None:
            fails.append("no scores")
        else:
            if nt["train"] < bt["train"]: fails.append("train")
            if nt["validate"] < bt["validate"] - 0.01: fails.append("validate")
            if nt["test"] < bt["test"] - args.test_tol: fails.append("test")
        if new["dlod"] is None or new["dlod"] < sweep.DLOD_FLOOR: fails.append(f"dLOD {new['dlod']}")
        if not sweep.backbone_ok(d): fails.append("backbone")
        if ref is not None and (m is None or m < args.min_mcc): fails.append("manifold")
        ok = not fails
        rec = dict(fails=fails, quad=q, base=base, new=new, base_fit=base_fit, new_fit=new_fit, mcc=m,
                   backbone=sweep.read_backbone(d), changed=(d / "lt.exe.p").read_bytes() != seed_p,
                   action="kept_new" if ok else "kept_old")
        if not ok and forced and prior is None and base["triplet"] is not None:
            #  --force-set quads must end up carrying the settings: install the
            #  baseline evaluated with them (same parameters as before).
            for f in FILES:
                if (d / "_base" / f).exists():
                    shutil.copy2(d / "_base" / f, d / f)
            rec["action"] = "kept_old_with_settings"
            ok = True
        if ok:
            bk = root / "_refine" / "backup" / q
            bk.mkdir(parents=True, exist_ok=True)
            for f in FILES:
                if (src / f).exists():
                    shutil.copy2(src / f, bk / f)
            sweep.write_lte_run_sh(d, overrides(d, args, False))   # leave the search invocation
            for f in FILES:
                if (d / f).exists():
                    shutil.copy2(d / f, src / f)
    except Exception as exc:
        rec = dict(quad=q, action="exception", error=repr(exc)[:300])
    finally:
        kill_lt_in(d)
    rec["elapsed_s"] = round(time.time() - t0)
    with LOCK, (args.root / "_refine" / "results.jsonl").open("a") as fh:
        fh.write(json.dumps(rec) + "\n")
    if rec["action"] == "exception":
        log(f"{q}: EXCEPTION {rec['error']}")
    else:
        b, n = rec["base"]["triplet"] or {}, rec["new"]["triplet"] or {}
        log(f"{q}: {rec['action']}  train/val/test {b.get('train', 0):.3f}/{b.get('validate', 0):.3f}/{b.get('test', 0):.3f}"
            f" -> {n.get('train', 0):.3f}/{n.get('validate', 0):.3f}/{n.get('test', 0):.3f}  "
            f"fit CC {rec['base_fit']} -> {rec['new_fit']}  dLOD {rec['new']['dlod']}  "
            f"{'rejected: ' + ', '.join(rec['fails']) if rec['fails'] else ''}  ({rec['elapsed_s']}s)")
    return rec


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    which = ap.add_mutually_exclusive_group()
    which.add_argument("--ALL", action="store_true", help="every grid quad (default)")
    which.add_argument("-set", nargs="+", metavar="QUAD")
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--ref", type=Path, default=None, help="reference manifold .dat for the quadrature gate")
    ap.add_argument("--min-mcc", type=float, default=0.5)
    ap.add_argument("--timeout", type=int, default=250)
    ap.add_argument("--threads", type=int, default=8,
                    help="search threads per lt.exe (NUMBER_OF_PROCESSORS)")
    ap.add_argument("--test-tol", type=float, default=0.10,
                    help="allowed drop in the 2000-05 test CC (a 5-year CC is noisy)")
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--force-set", action="append", default=[], metavar="KEY=VALUE",
                    help="resp setting applied to the --force-quads (e.g. ALPHA=0.25); those quads always "
                         "end up with it: the searched fit if it beats the baseline scored WITH the setting, "
                         "else that baseline")
    ap.add_argument("--force-quads", nargs="+", metavar="QUAD", default=[])
    ap.add_argument("--rescore-failed", action="store_true",
                    help="redo only the final re-score + decision for quads recorded with 'dLOD None'")
    args = ap.parse_args()
    args.root = args.root.resolve()
    ref = np.loadtxt(args.ref)[:, 1] if args.ref else None

    rdir = args.root / "_refine"
    rdir.mkdir(exist_ok=True)
    for link, tgt in (("lt.exe", args.root / "lt.exe"), ("dlod3.dat", args.root / "dlod3.dat")):
        if not (rdir / link).exists():
            (rdir / link).symlink_to(tgt)
    cells = sorted(p.name for p in args.root.iterdir() if p.is_dir() and sweep.parse_grid_name(p.name))
    targets = [c for c in (args.set or cells) if c in cells]
    res_path = rdir / "results.jsonl"
    done = {json.loads(l)["quad"] for l in res_path.read_text().splitlines() if l} if res_path.exists() else set()
    todo = [c for c in targets if c not in done]
    if args.rescore_failed:
        recs = [json.loads(l) for l in res_path.read_text().splitlines() if l]
        bad = {r["quad"]: r for r in recs if any(str(f).startswith("dLOD None") for f in r.get("fails", []))}
        res_path.write_text("".join(json.dumps(r) + "\n" for r in recs if r["quad"] not in bad))
        log(f"re-scoring {len(bad)} quads whose final re-score did not report")
        with ThreadPoolExecutor(args.workers) as ex:
            list(ex.map(lambda q: process(q, args, ref, prior=bad[q]), sorted(bad)))
        todo = []
    log(f"{len(targets)} quads, {len(done)} already done, {len(todo)} to run; "
        f"TIMEOUT {args.timeout}s, {args.threads} threads x {args.workers} workers, root {args.root.name}")
    with ThreadPoolExecutor(args.workers) as ex:
        list(ex.map(lambda q: process(q, args, ref), todo))
    recs = [json.loads(l) for l in res_path.read_text().splitlines() if l]
    log(f"DONE: {sum(r['action'] == 'kept_new' for r in recs)} improved, "
        f"{sum(r['action'].startswith('kept_old') for r in recs)} unchanged, "
        f"{sum(r['action'] == 'exception' for r in recs)} exceptions")
    with open(rdir / "global_sst_from_quads.txt", "w") as fh:
        subprocess.run([sys.executable, str(HERE / "global_sst_from_quads.py"), "--root", str(args.root),
                        "--out", str(rdir / "global_sst_from_quads.png")], stdout=fh, stderr=subprocess.STDOUT)
    (rdir / "REFINE_DONE").write_text(time.strftime("%Y-%m-%dT%H:%M:%S") + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
