#!/usr/bin/env python3
"""Command-line equivalent of the lte_gui.py "Run lt" button.

Runs the LTE calibration (lt.exe -> enso_opt) in a selected index
subdirectory. Each setting's default is the SELECTED SUBDIR's lt.exe.resp
value when it defines one (e.g. most index dirs set IDATE=1880.0, not the
1920.9 below); otherwise it falls back to the GUI defaults:

    JSON mode on (-j), METRIC=CC, TIMEOUT=1000, IDATE=1920.9,
    EXCLUDE=true, TREND=true, F9=1, TEST_ONLY=false,
    DLOD_REF=TRUE, LOCKT=FALSE, LOCKA=FALSE, ZONE=FALSE, RIDGE=0.0, sim off.

An explicit CLI flag always wins over both the resp file and the GUI
default. The only thing you normally change is the CV interval
(TRAIN_START / TRAIN_END), which defaults to 1880 .. 1885 (or the resp
file's TRAIN_START/TRAIN_END, if it sets them).

Usage:
    ./lte_run.py amo
    ./lte_run.py amo --cv 1880 1885
    ./lte_run.py nino4 --timeout 3600 --metric DTW
    ./lte_run.py baltic --no-exclude --cv 1940 2020 --ridge 0.1
    ./lte_run.py baltic --no-exclude --cv 1940 2020 --coverage 0.5
    ./lte_run.py baltic --cv 1940 2020 --enclosed
    ./lte_run.py baltic --cv 1940 2020 --enclosed u

Produces (in the index dir): lte_results.csv, updated lt.exe.p, lt.exe.resp.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

from lte_forward import read_resp

ROOT = Path(__file__).resolve().parent          # .../experiments/Feb2026

# GUI defaults (App.__init__ + run_lt) — used only when the selected
# subdir's lt.exe.resp doesn't define the corresponding setting.
DEFAULTS = {
    "cv_begin":      "1880",
    "cv_end":        "1885",
    "metric":        "CC",        # metric_var, upper() applied
    "timeout":       1000.0,      # time_var (seconds)
    "idate":         "1920.9",    # hardcoded in run_lt
    "exclude":       "true",      # exclude_var
    "trend":         "true",      # trend_var
    "f9":            "1",         # filter_var
    "test_only":     "false",     # test_only_var
    "dlod_ref":      "TRUE",      # lod_var ("no LOD cal" checked)
    "lockt":         "FALSE",     # tidal_phase_var
    "locka":         "FALSE",     # tidal_amplitude_var
    "zone":          "FALSE",     # zone_var
    "use_json":      True,        # use_json_var -> "-j" flag
    "ridge":         0.0,         # not in the GUI; matches lt.exe's own
                                   # GEM.Getenv("RIDGE", 0.0) default (unregularized OLS)
    "coverage":      1.0,         # not in the GUI; matches lt.exe's own
                                   # GEM.Getenv("COVERAGE", 1.0) default (whole interval)
    "enclosing":     "FALSE",     # not in the GUI; matches lt.exe's own
                                   # GEM.Getenv("ENCLOSING", False) default
    "enclosed":      "UL",        # not in the GUI; matches lt.exe's own
                                   # GEM.Getenv("ENCLOSED", "UL") default
}


def _peek_index(argv: list[str]) -> str | None:
    """Find the positional 'index' arg without fully parsing argv — used to
    locate and load the selected subdir's lt.exe.resp BEFORE the real
    argparser is built, so its values can become that parser's defaults."""
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("index", nargs="?")
    pre_args, _ = pre.parse_known_args(argv)
    return pre_args.index


def _resp_float(resp: dict, key: str, fallback: float) -> float:
    try:
        return float(resp[key])
    except (KeyError, ValueError):
        return fallback


def _bool_str(raw, true_str: str, false_str: str) -> str:
    """Interpret a resp/env boolean string the same way Ada's GEM.Getenv
    does (case-insensitive TRUE/T/1), rendered in the caller's casing
    convention (e.g. "true"/"false" for EXCLUDE, "TRUE"/"FALSE" for ZONE)."""
    return true_str if str(raw).strip().upper() in ("TRUE", "T", "1") else false_str


def resolve_lt_cmd(run_dir: Path, use_json: bool) -> list[str]:
    """Same as lte_gui.resolve_lt_cmd: look for lt.exe / lt in the parent."""
    json_flag = ["-j"] if use_json else []
    for name in ("lt.exe", "lt"):
        lt_path = ROOT / name
        if lt_path.exists():
            return [str(lt_path)] + json_flag
    raise FileNotFoundError(f"No lt or lt.exe found in {ROOT}")


def build_env(index: str, args: argparse.Namespace, resp: dict) -> dict:
    """Build the subprocess env. --cv/--metric/--timeout/--idate already
    carry resp-derived defaults from argparse (see main()); the store_true
    flags below default to None (not passed) so an explicit flag can be
    told apart from "use the resp file's value, or the GUI default"."""
    env = os.environ.copy()
    env["METRIC"] = args.metric.strip().upper()
    env["TIMEOUT"] = f"{args.timeout:.6f}"
    env["TRAIN_START"] = str(args.cv[0]).strip()
    env["TRAIN_END"] = str(args.cv[1]).strip()
    env["CLIMATE_INDEX"] = f"{index}.dat"
    env["IDATE"] = args.idate
    env["RIDGE"] = f"{args.ridge:.10g}"
    env["COVERAGE"] = f"{args.coverage:.10g}"
    env["EXCLUDE"] = ("false" if args.no_exclude else
                       _bool_str(resp.get("EXCLUDE", DEFAULTS["exclude"]), "true", "false"))
    env["TREND"] = ("false" if args.no_trend else
                     _bool_str(resp.get("TREND", DEFAULTS["trend"]), "true", "false"))
    env["F9"] = "0" if args.no_filter else str(resp.get("F9", DEFAULTS["f9"]))
    env["TEST_ONLY"] = ("true" if args.test_only else
                         _bool_str(resp.get("TEST_ONLY", DEFAULTS["test_only"]), "true", "false"))
    env["DLOD_REF"] = ("FALSE" if args.no_lod else
                        _bool_str(resp.get("DLOD_REF", DEFAULTS["dlod_ref"]), "TRUE", "FALSE"))
    env["LOCKT"] = ("TRUE" if args.lock_tidal_phase else
                     _bool_str(resp.get("LOCKT", DEFAULTS["lockt"]), "TRUE", "FALSE"))
    env["LOCKA"] = ("TRUE" if args.lock_tidal_amplitude else
                     _bool_str(resp.get("LOCKA", DEFAULTS["locka"]), "TRUE", "FALSE"))
    env["ZONE"] = ("TRUE" if args.zone else
                    _bool_str(resp.get("ZONE", DEFAULTS["zone"]), "TRUE", "FALSE"))
    env["ENCLOSING"] = ("TRUE" if args.enclosed is not None else
                         _bool_str(resp.get("ENCLOSING", DEFAULTS["enclosing"]), "TRUE", "FALSE"))
    env["ENCLOSED"] = (args.enclosed.upper() if args.enclosed is not None
                        else resp.get("ENCLOSED", DEFAULTS["enclosed"]))
    return env


def main() -> int:
    # Peek the positional 'index' before the real parser exists, so the
    # selected subdir's lt.exe.resp can supply this parser's defaults
    # instead of the hardcoded GUI DEFAULTS.
    index_guess = _peek_index(sys.argv[1:])
    resp = {}
    if index_guess:
        guess_dir = (ROOT / index_guess).resolve()
        if guess_dir.is_dir():
            resp = read_resp(guess_dir / "lt.exe.resp")

    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("index", help="subdirectory to run, e.g. amo")
    ap.add_argument("--cv", nargs=2, metavar=("BEGIN", "END"),
                    default=[resp.get("TRAIN_START", DEFAULTS["cv_begin"]),
                             resp.get("TRAIN_END", DEFAULTS["cv_end"])],
                    help="CV interval TRAIN_START TRAIN_END "
                         "(default: from lt.exe.resp, else 1880 1885)")
    ap.add_argument("--metric", default=resp.get("METRIC", DEFAULTS["metric"]),
                    choices=["CC", "cc", "DTW", "dtw", "H", "HOYER", "hoyer"],
                    help="objective metric (default: from lt.exe.resp, else CC)")
    ap.add_argument("--timeout", type=float,
                    default=_resp_float(resp, "TIMEOUT", DEFAULTS["timeout"]),
                    help="TIMEOUT seconds passed to the model "
                         "(default: from lt.exe.resp, else 1000)")
    # store_true with default=None (not False) so build_env can tell "flag
    # explicitly passed" apart from "not passed, use the resp/GUI default".
    ap.add_argument("--test-only", action="store_true", default=None,
                    help="force TEST_ONLY=true (default: from lt.exe.resp, else false)")
    ap.add_argument("--zone", action="store_true", default=None,
                    help="force ZONE=TRUE (default: from lt.exe.resp, else FALSE)")
    ap.add_argument("--enclosed", choices=["u", "ul"], nargs="?", const="ul",
                    default=None,
                    help="force ENCLOSING=TRUE — fit only on the data BEFORE "
                         "TRAIN_START, and gate accept/reject on ENCLOSED's "
                         "choice of region(s): 'ul' (default when --enclosed "
                         "is passed with no value) scores LOWER+UPPER "
                         "combined, so the search can't drift away from "
                         "LOWER's own pattern without penalty; 'u' scores "
                         "UPPER (after TRAIN_END) only, a stricter pure-"
                         "generalization signal. Either way, the enclosed "
                         "[TRAIN_START,TRAIN_END] interval is reported as "
                         "the true held-out validation. Standalone, like "
                         "--coverage — takes precedence over EXCLUDE/"
                         "SPLIT_TRAINING, and --coverage wins if both are "
                         "set (default: from lt.exe.resp, else off)")
    ap.add_argument("--no-lod", action="store_true", default=None,
                    help="force DLOD_REF=FALSE (default: from lt.exe.resp, else TRUE)")
    ap.add_argument("--lock-tidal-phase", action="store_true", default=None,
                    help="force LOCKT=TRUE (default: from lt.exe.resp, else FALSE)")
    ap.add_argument("--lock-tidal-amplitude", action="store_true", default=None,
                    help="force LOCKA=TRUE (default: from lt.exe.resp, else FALSE)")
    ap.add_argument("--no-exclude", action="store_true", default=None,
                    help="force EXCLUDE=false (default: from lt.exe.resp, else true)")
    ap.add_argument("--no-trend", action="store_true", default=None,
                    help="force TREND=false (default: from lt.exe.resp, else true)")
    ap.add_argument("--no-filter", action="store_true", default=None,
                    help="force F9=0 (default: from lt.exe.resp, else 1)")
    ap.add_argument("--no-json", action="store_true",
                    help="omit the -j JSON flag (GUI checks JSON by default)")
    ap.add_argument("--idate", default=resp.get("IDATE", DEFAULTS["idate"]),
                    help="default: from lt.exe.resp, else 1920.9")
    ap.add_argument("--ridge", type=float,
                    default=_resp_float(resp, "RIDGE", DEFAULTS["ridge"]),
                    help="L2/ridge penalty on Regression_Coefficients' normal "
                         "equations — shrinks the OLS-fit level/k0/amp-phase/"
                         "trend/accel/annual coefficients (0.0 = unregularized "
                         "OLS, lt.exe's own default; default: from lt.exe.resp, "
                         "else 0.0)")
    ap.add_argument("--coverage", type=float,
                    default=_resp_float(resp, "COVERAGE", DEFAULTS["coverage"]),
                    help="fraction (0.0-1.0) of the training interval actually "
                         "used to fit the regression; the resulting fit is "
                         "scored on the WHOLE training interval, and that "
                         "whole-interval score directly gates each search "
                         "thread's own accept/reject — rejects parameter sets "
                         "that only fit the coverage slice, not the withheld "
                         "remainder (1.0 = whole interval, lt.exe's own "
                         "default, i.e. no effect; default: from lt.exe.resp, "
                         "else 1.0)")
    args = ap.parse_args()

    run_dir = (ROOT / args.index).resolve()
    if not run_dir.is_dir():
        print(f"error: index directory not found: {run_dir}", file=sys.stderr)
        return 2
    if not (run_dir / f"{args.index}.dat").exists():
        print(f"error: {run_dir / f'{args.index}.dat'} not found "
              "(CLIMATE_INDEX target)", file=sys.stderr)
        return 2

    # Reload from the now-authoritative run_dir (args.index), rather than
    # trusting the pre-parse guess above, before it feeds build_env.
    resp = read_resp(run_dir / "lt.exe.resp")
    index = run_dir.name.strip().strip('"').strip("'")

    # GUI run_lt: delete stale param file so the model regenerates it.
    par = run_dir / f"lt.exe.{index}.dat.par"
    if par.exists():
        par.unlink()

    try:
        lt_cmd = resolve_lt_cmd(run_dir, use_json=not args.no_json)
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    env = build_env(index, args, resp)

    print(f"run dir : {run_dir}")
    print(f"cmd     : {' '.join(lt_cmd)}")
    print(f"cv      : {env['TRAIN_START']} .. {env['TRAIN_END']}   "
          f"metric={env['METRIC']}  timeout={env['TIMEOUT']}  ridge={env['RIDGE']}  "
          f"coverage={env['COVERAGE']}  enclosing={env['ENCLOSING']}"
          f"{'/' + env['ENCLOSED'] if env['ENCLOSING'] == 'TRUE' else ''}")

    # GUI _open_terminal on headless Linux: bash -c 'ulimit -s unlimited && <cmd>'
    bare_cmd = " ".join(_shq(a) for a in lt_cmd)
    proc = subprocess.run(
        ["bash", "-c", f"ulimit -s unlimited && {bare_cmd}"],
        cwd=str(run_dir), env=env,
    )
    if proc.returncode != 0:
        print(f"error: {lt_cmd[0]} exited with code {proc.returncode}",
              file=sys.stderr)
        return proc.returncode

    results = run_dir / "lte_results.csv"
    print(f"OK: wrote {results}" if results.exists()
          else "finished, but no lte_results.csv was produced")
    return 0


def _shq(s: str) -> str:
    import shlex
    return shlex.quote(s)


if __name__ == "__main__":
    sys.exit(main())
