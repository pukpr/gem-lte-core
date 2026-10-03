#!/usr/bin/env python3
"""NM=3 layered AMO: ltep = [layer k1, Bessel kB, regressed kW] (2026-10-02 20:07 build).
TEST_ONLY scores for a grid of (start parameters, kB, kW); usage: three.py"""
import json, re, shutil, sys, itertools
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
sys.argv = [sys.argv[0]]; from single_k import B, src, ov0, score, sweep
STARTS = {"amo": src / "lt.exe.p", "s25": B / "search_k2.5" / "lt.exe.p"}
def make3(name, start, kB, kW):
    d = B / name
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.resp", "dlod_ref.dat"): shutil.copy(src / f, d / f)
    shutil.copy(src / "amo.dat", d / f"{name}.dat")
    p = json.loads(STARTS[start].read_text()); p["ltep"] = [p["ltep"][0], kB, kW]; (d / "lt.exe.p").write_text(json.dumps(p, indent=2))
    rp = d / "lt.exe.resp"; sweep.set_resp_key(rp, "NM", "3"); rp.write_text(re.sub(r'^LAYER_K\s.*\n', '', rp.read_text(), flags=re.M))
    return d
if __name__ == "__main__":
    # regression check: installed amo (NM 2 + LAYER_K) must still score 0.7450
    d = B / "check_installed"
    if d.exists(): shutil.rmtree(d)
    d.mkdir()
    for f in ("lt.exe.resp", "lt.exe.p", "dlod_ref.dat"): shutil.copy(src / f, d / f)
    shutil.copy(src / "amo.dat", d / "check_installed.dat")
    jobs = [("check", None, None, None)] + [(f"t3_{s}_B{kB:g}_W{kW:g}", s, kB, kW) for s, kB, kW in itertools.product(("amo", "s25"), (0.231889340633804, 2.5, 2.6), (4.0, 4.2, 4.3, 4.5))]
    def run(j):
        name, s, kB, kW = j
        return score(d) if s is None else score(make3(name, s, kB, kW))
    with ThreadPoolExecutor(8) as ex: recs = list(ex.map(run, jobs))
    for x in recs:
        t = x["triplet"]; print(f"{x['name']:24s} ltep {[round(v,4) for v in x['ltep']]}: {t['train']:.4f}/{t.get('validate',0):.4f}/{t['test']:.4f} dLOD {x['dlod']:.4f} amps {[round(a[1],3) for a in x['kap']]}")
    json.dump(recs, open(B / "three.json", "w"), indent=1)
