#!/usr/bin/env python3
"""NAO-family legs for the directed gradient run (2026-11-21), all on nao.dat's footing: monthly, 13-month boxcar.
nao.dat is nao_station_monthly.dat through a 13-month boxcar (CC 1.00), so the station series is not used again.
  nao.dat        copied from experiments/Feb2026/nao (already smoothed)
  nao_gib.dat    Gibraltar-Iceland index (GEM-LTE/run/nao_gib.dat), 13-month boxcar
  azores.dat, gibraltar.dat, reykjavik.dat   Hadley station MSLP (~/eval/hadley/stations, 1850-2004): monthly anomalies
                 (calendar-month means removed), 13-month boxcar
  azores_reykjavik.dat   Azores minus Reykjavik anomaly, 13-month boxcar
Boxcar: centred, full windows only (6 months dropped at each end and around gaps). Dates as year + (month-1)/12.
    python3 make_data.py"""
import os, numpy as np
H = os.path.expanduser("~/eval"); OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
def fts(f):
    out = {}
    for l in open(f).readlines()[1:]:
        x = l.split(); y = int(x[0])
        for m, v in enumerate(x[1:13]):
            if float(v) > -900: out[y * 12 + m] = float(v)
    return out
def anomalies(d):
    ks = np.array(sorted(d)); v = np.array([d[k] for k in ks])
    for m in range(12): v[ks % 12 == m] -= v[ks % 12 == m].mean()
    return dict(zip(ks, v))
def dat(f):
    d = np.loadtxt(f); return {int(round(t * 12)): v for t, v in d[:, :2]}
def box13(d):
    out = {}
    for k in d:
        w = [d.get(k + j) for j in range(-6, 7)]
        if all(x is not None for x in w): out[k] = sum(w) / 13
    return out
def write(name, d, note):
    ks = sorted(d)
    with open(os.path.join(OUT, name), "w") as fh:
        for k in ks: fh.write(f"{k // 12 + (k % 12) / 12:.6f}\t{d[k]:.6f}\n")
    print(f"{name:22s} {ks[0] // 12}-{ks[-1] // 12}  {len(ks)} months  {note}")
S = os.path.join(H, "hadley/stations")
az, gb, rk = (anomalies(fts(os.path.join(S, f))) for f in ("08506_Azores_MSLP_18502004_v2.fts", "08495_Gibraltar_MSLP_18502004_v2.fts", "04030_Reykjavik_MSLP_18502004_v2.fts"))
import shutil; shutil.copy(os.path.join(H, "GEM-LTE/experiments/Feb2026/nao/nao.dat"), os.path.join(OUT, "nao.dat")); print("nao.dat                copied")
write("nao_gib.dat", box13(dat(os.path.join(H, "GEM-LTE/run/nao_gib.dat"))), "Gibraltar-Iceland index")
write("azores.dat", box13(az), "Azores MSLP anomaly (hPa)")
write("gibraltar.dat", box13(gb), "Gibraltar MSLP anomaly (hPa)")
write("reykjavik.dat", box13(rk), "Reykjavik MSLP anomaly (hPa)")
write("azores_reykjavik.dat", box13({k: az[k] - rk[k] for k in az if k in rk}), "Azores - Reykjavik (hPa)")
