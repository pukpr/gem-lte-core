#!/usr/bin/env python3
"""Monthly axial atmospheric angular momentum (AAM) as a climate-index .dat.

Source: GFZ ESMGFZ AAM v1.0, 3-hourly effective angular momentum functions from
ECMWF operational analysis, 1976 onward (atmospheric pressure tides removed,
inverted barometer over the ocean), downloaded 2026-11-12 from
https://rz-vm480.gfz.de/files/ESMGFZ/EAM/operational_AAM/ into gfz_raw/.

The axial (z) excitation chi3 = mass + motion term; the length-of-day change it
excites is dLOD = LOD0 * chi3, LOD0 = 86400 s, written here in ms.

Writes (date = year + (month - 1) / 12, as in the Feb2026 .dat files; complete
months only):
    aam.dat        total chi3 (mass + motion) in ms of LOD
    aam_anom.dat   aam.dat minus its 1976-2025 monthly climatology
    aam_motion.dat, aam_mass.dat   the two terms separately (ms)

Also the longer NOAA PSL globally integrated (relative, wind) AAM, NCEP/NCAR
reanalysis, 1958 - Feb 2014, 3-month running means (Weickmann et al. 2000),
https://psl.noaa.gov/data/correlation/glaam.data -> psl_glaam.data:
    glaam.dat      in ms of LOD (dLOD = LOD0 * dH / (C * Omega)), its small
                   remaining monthly climatology removed
"""
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LOD0_MS = 86400.0e3


def read_3h():
    rows = []
    for f in sorted((HERE / "gfz_raw").glob("ESMGFZ_AAM_v1.0_03h_[12][0-9][0-9][0-9].asc")):
        for line in f.read_text().splitlines():
            p = line.split()
            if len(p) == 11 and p[0].isdigit() and len(p[0]) == 4:
                rows.append([float(v) for v in p])
    a = np.array(rows)
    # columns: Y M D H MJD  mass x y z  motion x y z
    return a[:, 0].astype(int), a[:, 1].astype(int), a[:, 7], a[:, 10]


def monthly(y, m, v):
    keys = y * 12 + (m - 1)
    out = {}
    for k in np.unique(keys):
        sel = keys == k
        yy, mm = divmod(int(k), 12)
        days = (np.datetime64(f"{yy + (mm + 1) // 12:04d}-{(mm + 1) % 12 + 1:02d}-01")
                - np.datetime64(f"{yy:04d}-{mm + 1:02d}-01")).astype(int)
        if sel.sum() == 8 * days:  # complete month of 3-hourly records
            out[yy + mm / 12.0] = v[sel].mean()
    d = np.array(sorted(out))
    return d, np.array([out[k] for k in d])


def write(name, d, v):
    with open(HERE / name, "w") as f:
        for t, x in zip(d, v):
            f.write(f"{t:.6f}\t{x:.6f}\n")


C_OMEGA = 8.0365e37 * 7.292115e-5   # Earth's polar moment of inertia (kg m^2) x rotation rate


def glaam():
    d, v = [], []
    for line in (HERE / "psl_glaam.data").read_text().splitlines()[1:]:
        p = line.split()
        if not p or not re.fullmatch(r"\d{4}", p[0]):
            continue
        for k, x in enumerate(re.findall(r"-?\d+\.\d+E[+-]\d+", line[5:])):
            x = float(x)
            if abs(x + 99.99) > 1e-6 and x != 0.0:   # -99.99 missing; Mar 2014 is a 0.0 placeholder
                d.append(int(p[0]) + k / 12.0)
                v.append(x * LOD0_MS / C_OMEGA)
    d, v = np.array(d), np.array(v)
    month = np.rint((d % 1) * 12).astype(int)
    clim = np.array([v[month == k].mean() for k in range(12)])
    write("glaam.dat", d, v - clim[month])
    print(f"glaam: {len(d)} months {d[0]:.3f}-{d[-1]:.3f}; sd {v.std():.3f} ms, "
          f"remaining seasonal range {np.ptp(clim):.3f} ms removed")


def main():
    glaam()
    y, m, mass, motion = read_3h()
    d, tot = monthly(y, m, (mass + motion) * LOD0_MS)
    _, ma = monthly(y, m, mass * LOD0_MS)
    _, mo = monthly(y, m, motion * LOD0_MS)
    month = np.rint((d % 1) * 12).astype(int)
    base = (d >= 1976) & (d < 2026)
    clim = np.array([tot[base & (month == k)].mean() for k in range(12)])
    write("aam.dat", d, tot)
    write("aam_anom.dat", d, tot - clim[month])
    write("aam_mass.dat", d, ma)
    write("aam_motion.dat", d, mo)
    print(f"{len(d)} months {d[0]:.3f}-{d[-1]:.3f}; total sd {tot.std():.3f} ms, "
          f"anomaly sd {(tot - clim[month]).std():.3f} ms, seasonal range {np.ptp(clim):.3f} ms")


if __name__ == "__main__":
    main()
