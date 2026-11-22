#!/usr/bin/env python3
"""Wisdom-of-the-crowds tidal set from every saved lpap.

Scans for lt.exe.p (in use) and lt.exe.p.* (backups) files, reads each one's
"lpap" list of [period, amplitude, phase] triples, groups matching periods
across files, and writes one averaged triple per period.

  * Periods: sorted and split wherever the relative gap to the previous
    period exceeds --tol (default 2e-4). The real tidal lines in these files
    are at least ~1e-3 apart (e.g. 9.5435 vs 9.5568 d), while one line's
    refits differ by ~1e-5. The output period is the group's median.
  * Fast lines: rows with |period| < --min-period (default 2 d: the diurnal
    and semidiurnal lines, 0.50-1.12 d) are dropped before anything else.
  * Variance: each set's amplitudes are then scaled so the set's own series
    has unit variance, sum(A^2)/2 = 1, so a fit with runaway amplitudes
    weighs no more than any other (off with --no-normalize). Amplitudes in
    the output are then relative.
  * Sign: A*f(wt+phi) = |A|*f(wt+phi+pi), so a negative amplitude is made
    positive and its phase shifted by pi (turn off with --keep-sign).
  * Phases are folded into [0, 2*pi) and averaged as a circular mean. R in
    [0,1] measures agreement: 1 = identical phases, ~0 = uniformly spread.
  * Amplitude: median of |A| by default (--amp median). The saved fits have
    heavy tails -- means run 5-50x the medians, and one runaway fit puts the
    182.6 d mean at ~1e4 -- so --amp mean is the plain average and --amp trim
    a 10%-each-end trimmed mean. All three go to the CSV.
  * Exclusions: a set is dropped when its folder path, or any lt.exe.resp*
    in its folder, matches --exclude (default 'qbo', case-insensitive), and so
    is every copy of that lpap saved elsewhere. --exclude '' keeps them all.
  * Row order matters: the model assigns each lpap row, in file order, to the
    nearest unmatched period of its built-in table within 1%. Sorted by
    period, a rare period missing from the table (e.g. 9.0829 d) takes a
    real line's slot and that line (e.g. 9.1330 d) is dropped. So rows are
    written most-common first, and with --template REF.lt.exe.p they follow
    REF's lpap order exactly: REF periods without a crowd value get amplitude
    0 (e.g. the dropped diurnal lines), and crowd periods outside REF's table
    are left out of the lpap (they stay in the CSV, marked in_lpap=0).
  * Selection: --include REGEX keeps only files whose full path (file name
    or any folder name) matches it, case-insensitive; e.g. --include wow
    keeps the fits saved as exceptional (lt.exe.p.*wow*, *WOW*, or any
    lt.exe.p* in a *.wow* folder). --skip REGEX drops files by path alone,
    without the copy rule of --exclude (use it for scratch copies of fits
    that should still count where they were first saved).
  * Duplicates: identical lpap lists (copies of the same fit) count once
    unless --no-dedupe is given.

usage: lpap_wisdom.py [ROOT ...] [-o wisdom.lt.exe.p] [--csv stats.csv]
       default ROOT is / (skipping /proc /sys /dev /run /mnt)
"""
import argparse, csv, hashlib, json, math, os, re, statistics, sys

SKIP = {'/proc', '/sys', '/dev', '/run', '/mnt'}
TWO_PI = 2 * math.pi


def find_files(roots):
    for root in roots:
        for d, dirs, files in os.walk(root, onerror=lambda e: None):
            dirs[:] = [x for x in dirs if os.path.join(d, x) not in SKIP]
            for f in files:
                if f == 'lt.exe.p' or f.startswith('lt.exe.p.'):
                    p = os.path.join(d, f)
                    if os.path.isfile(p):
                        yield p, ('used' if f == 'lt.exe.p' else 'backup')


def load_lpap(path):
    try:
        with open(path) as fh:
            p = json.load(fh)
    except (OSError, ValueError, UnicodeDecodeError):
        return None
    l = p.get('lpap') if isinstance(p, dict) else None
    if not isinstance(l, list):
        return None
    rows = [r for r in l if isinstance(r, list) and len(r) >= 3
            and all(isinstance(v, (int, float)) and math.isfinite(v) for v in r[:3])]
    return rows or None


def excluded_dir(d, pat, cache={}):
    if d not in cache:
        hit = bool(pat.search(d))
        if not hit:
            try:
                for f in os.listdir(d):
                    if f.startswith('lt.exe.resp'):
                        with open(os.path.join(d, f), errors='ignore') as fh:
                            if pat.search(fh.read()):
                                hit = True
                                break
            except OSError:
                pass
        cache[d] = hit
    return cache[d]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('roots', nargs='*', default=['/'])
    ap.add_argument('-o', '--out', default='wisdom.lt.exe.p', help='output JSON with the averaged lpap')
    ap.add_argument('--csv', default='wisdom_lpap_stats.csv', help='per-period statistics')
    ap.add_argument('--tol', type=float, default=2e-4, help='relative period gap that starts a new group')
    ap.add_argument('--min-frac', type=float, default=0.0,
                    help='drop periods present in less than this fraction of the sets')
    ap.add_argument('--exclude', default='qbo',
                    help="regex (case-insensitive) on folder path / lt.exe.resp text; '' = none")
    ap.add_argument('--include', default='',
                    help="regex (case-insensitive) a file's full path must match; '' = all")
    ap.add_argument('--skip', default='',
                    help="regex (case-insensitive) on the file path: drop that file only")
    ap.add_argument('--no-dedupe', action='store_true', help='count copies of the same lpap separately')
    ap.add_argument('--min-period', type=float, default=2.0,
                    help='drop rows with |period| below this many days (default 2)')
    ap.add_argument('--no-normalize', action='store_true',
                    help='do not scale each set to unit variance first')
    ap.add_argument('--amp', choices=['median', 'mean', 'trim'], default='median',
                    help='amplitude estimator written to the lpap (default median)')
    ap.add_argument('--template', help='lt.exe.p whose lpap order/period table the output follows')
    ap.add_argument('--keep-sign', action='store_true', help='do not fold negative amplitudes into the phase')
    a = ap.parse_args()

    pat = re.compile(a.exclude, re.I) if a.exclude else None
    inc = re.compile(a.include, re.I) if a.include else None
    skp = re.compile(a.skip, re.I) if a.skip else None
    nsel = 0
    found, bad = [], set()          # found: (path, kind, rows, hash)
    nfiles = {'used': 0, 'backup': 0}
    nolpap = 0
    for path, kind in find_files(a.roots):
        if (inc and not inc.search(path)) or (skp and skp.search(path)):
            continue
        nsel += 1
        rows = load_lpap(path)
        if rows is None:
            nolpap += 1
            continue
        h = hashlib.md5(json.dumps([r[:3] for r in rows]).encode()).hexdigest()
        if pat and excluded_dir(os.path.dirname(path), pat):
            bad.add(h)
        found.append((path, kind, rows, h))
    seen, sets, nexcl = set(), [], 0
    for path, kind, rows, h in found:
        if h in bad:
            nexcl += 1
            continue
        nfiles[kind] += 1
        if not a.no_dedupe:
            if h in seen:
                continue
            seen.add(h)
        sets.append((path, kind, rows))
    if not sets:
        sys.exit('no lpap found')

    # every triple, tagged with its set index
    trip = []
    nzero = 0
    for i, (_, _, rows) in enumerate(sets):
        rows = [r[:3] for r in rows if abs(r[0]) >= a.min_period]
        scale = 1.0
        if not a.no_normalize:
            var = sum(r[1] ** 2 for r in rows) / 2
            if not var > 0 or not math.isfinite(var):
                nzero += 1
                continue
            scale = 1 / math.sqrt(var)
        for per, amp, ph in rows:
            amp *= scale
            if amp < 0 and not a.keep_sign:
                amp, ph = -amp, ph + math.pi
            trip.append((per, amp, ph % TWO_PI, i))
    trip.sort(key=lambda t: t[0])

    groups, cur = [], [trip[0]]
    for t in trip[1:]:
        prev = cur[-1][0]
        if abs(t[0] - prev) > a.tol * max(abs(prev), 1e-12) or (t[0] < 0) != (prev < 0):
            groups.append(cur)
            cur = [t]
        else:
            cur.append(t)
    groups.append(cur)

    out_rows, stats, warn = [], [], 0
    for g in groups:
        frac = len({t[3] for t in g}) / len(sets)
        if frac < a.min_frac:
            continue
        if len({t[3] for t in g}) < len(g):
            warn += 1               # one set has two rows inside the tolerance
        per = statistics.median(t[0] for t in g)
        amps = [t[1] for t in g]
        c = sum(math.cos(t[2]) for t in g) / len(g)
        s = sum(math.sin(t[2]) for t in g) / len(g)
        R = math.hypot(c, s)
        phase = math.atan2(s, c) % TWO_PI
        srt = sorted(amps)
        k = len(srt) // 10
        est = {'median': statistics.median(srt), 'mean': statistics.fmean(srt),
               'trim': statistics.fmean(srt[k:len(srt) - k] or srt)}
        amp = est[a.amp]
        out_rows.append([per, amp, phase, frac])
        stats.append({'period': per, 'period_min': g[0][0], 'period_max': g[-1][0],
                      'n': len(g), 'frac_of_sets': round(frac, 4),
                      'amp_median': est['median'], 'amp_trim10': est['trim'], 'amp_mean': est['mean'],
                      'amp_sd': statistics.pstdev(amps), 'phase_mean': phase, 'phase_R': R})

    for st in stats:
        st['in_lpap'] = 1
    if a.template:
        tmpl = [r[0] for r in load_lpap(a.template)]
        lp, used = [], set()
        for tp in tmpl:
            best = min((i for i in range(len(out_rows)) if i not in used),
                       key=lambda i: abs(abs(out_rows[i][0]) - abs(tp)), default=None)
            if best is not None and abs(abs(out_rows[best][0]) - abs(tp)) <= a.tol * abs(tp):
                used.add(best)
                lp.append([tp] + out_rows[best][1:3])
            else:
                lp.append([tp, 0.0, 0.0])
        for i, st in enumerate(stats):
            st['in_lpap'] = int(i in used)
        nofit = sum(1 for r in lp if r[1] == 0.0)
        print(f"template {a.template}: {len(lp)} rows, {nofit} set to 0, "
              f"{len(out_rows) - len(used)} crowd periods not in its table")
    else:
        lp = [r[:3] for r in sorted(out_rows, key=lambda r: -r[3])]
    with open(a.out, 'w') as fh:
        json.dump({'lpap': lp,
                   '_wisdom': {'sets': len(sets), 'files_used': nfiles['used'],
                               'files_backup': nfiles['backup'], 'deduped': not a.no_dedupe,
                               'sign_folded': not a.keep_sign, 'min_period': a.min_period,
                               'normalized': not a.no_normalize, 'tol': a.tol, 'amp': a.amp,
                               'roots': a.roots, 'exclude': a.exclude,
                               'include': a.include, 'skip': a.skip,
                               'files_excluded': nexcl, 'template': a.template}}, fh, indent=2)
    with open(a.csv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(stats[0]))
        w.writeheader()
        w.writerows(stats)

    if inc or skp:
        print(f"selected {nsel} files (include {a.include!r}, skip {a.skip!r})")
    if nexcl:
        print(f"excluded {nexcl} files matching {a.exclude!r} (or copies of their lpap)")
    print(f"files with lpap: {nfiles['used']} used + {nfiles['backup']} backup "
          f"({nolpap} without); {len(sets)} {'distinct ' if not a.no_dedupe else ''}sets")
    print(f"{len(out_rows)} periods -> {a.out}, stats -> {a.csv}")
    if nzero:
        print(f"skipped {nzero} sets with zero or non-finite variance")
    if warn:
        print(f"warning: {warn} groups hold two rows from the same set; try a smaller --tol")


if __name__ == '__main__':
    main()
