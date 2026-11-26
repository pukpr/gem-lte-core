#!/usr/bin/env python3
"""Summarize the 30-quad basin-seeded fast-refine Kaplan sweep.

Writes quad5_besk_summary.csv and quad5_besk_summary.png.
"""
import csv
import json
import re
import statistics as st
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
sample = json.loads((HERE / "sample.json").read_text())
log_txt = (HERE / "sweep.log").read_text()

pat_run = re.compile(r"^(\S+): \S+ \| train ([+-][\d.]+).*?held out \S+ ([+-][\d.]+)", re.M)
pat_ref = re.compile(
    r"REFINE summary: mean ([\d.]+) -> ([\d.]+); dLOD ([\d.]+) -> ([\d.]+)\s+\[total ([\d.]+) s\]"
)
pat_name = re.compile(r"^k([NS])(\d{3})_([EW])(\d{3})$")


def latlon(q):
    m = pat_name.match(q)
    return (
        int(m.group(2)) * (1 if m.group(1) == "N" else -1),
        int(m.group(4)) * (1 if m.group(3) == "E" else -1),
    )


def cvec(path):
    return np.array([x[1] * np.exp(1j * x[2]) for x in json.load(open(path))["k_amp_phase"]])


def agree(a, b):
    return float(np.real(np.vdot(a, b)) / np.sqrt(np.sum(abs(a) ** 2) * np.sum(abs(b) ** 2)))


seed_scores = {}
cur_q = None
for line in log_txt.splitlines():
    m = re.search(r" (k[NS]\d{3}_[EW]\d{3}) \((atlantic|pacific|indian)\): set up", line)
    if m:
        cur_q = m.group(1)
        seed_scores[cur_q] = []
    elif cur_q and line.strip().startswith("seed,"):
        ms = re.search(r"train ([+-][\d.]+)\s+held out \S+ ([+-][\d.]+)", line)
        if ms:
            seed_scores[cur_q].append((float(ms.group(1)), float(ms.group(2))))

rows = []
for q, basin in sorted(sample.items()):
    d = HERE / q
    if not (d / "legs.png").exists():
        continue
    s0 = json.loads((d / "seed.p").read_text())
    p = json.loads((d / "result.p").read_text())
    new = [(float(a), float(b)) for _, a, b in pat_run.findall((d / "run.out").read_text())]
    old = seed_scores.get(q, [])
    m_ref = pat_ref.search((d / "refine.out").read_text())
    dl0, dl1, rsec = (
        (float(m_ref.group(3)), float(m_ref.group(4)), float(m_ref.group(5)))
        if m_ref
        else (float("nan"), float("nan"), float("nan"))
    )
    C = {k: cvec(d / f"lt.exe.windings.{k}.json") for k in range(1, 6)}
    mean_ag = float(np.mean([agree(C[1], C[k]) for k in range(2, 6)]))
    lat, lon = latlon(q)
    rows.append(
        dict(
            quad=q,
            basin=basin,
            lat=lat,
            lon=lon,
            besk_seed=s0["besk"],
            besk=p["besk"],
            base_seed=s0["ltep"][-1],
            base_k=p["ltep"][-1],
            train_seed=st.mean(x for x, _ in old),
            train_new=st.mean(x for x, _ in new),
            held_seed=st.mean(y for _, y in old),
            held_new=st.mean(y for _, y in new),
            centre_held_seed=old[0][1],
            centre_held_new=new[0][1],
            dlod0=dl0,
            dlod1=dl1,
            refine_s=rsec,
            mean_agree=mean_ag,
        )
    )

with open(HERE / "quad5_besk_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

print(
    f"{len(rows)} of {len(sample)} quads finished | "
    f"mean REFINE time: {st.mean(r['refine_s'] for r in rows):.2f} s"
)
print(
    f"{'quad':11s} {'basin':8s} {'besk s->f':>15s} | {'train s->f':>13s} | "
    f"{'held5 s->f':>15s} | {'ctr_ho s->f':>15s} | {'dLOD':>6s} | {'agree':>6s}"
)
for r in rows:
    print(
        f"{r['quad']:11s} {r['basin']:8s} {r['besk_seed']:.4f}->{r['besk']:.4f} | "
        f"{r['train_seed']:.3f}->{r['train_new']:.3f} | "
        f"{r['held_seed']:+.3f}->{r['held_new']:+.3f} | "
        f"{r['centre_held_seed']:+.3f}->{r['centre_held_new']:+.3f} | "
        f"{r['dlod1']:.4f} | {r['mean_agree']:+.2f}"
    )

print("\n=== BASIN SUMMARY ===")
for b in ("atlantic", "pacific", "indian", "ALL"):
    v = rows if b == "ALL" else [r for r in rows if r["basin"] == b]
    print(
        f"{b:8s} (n={len(v):2d}): besk med {st.median(r['besk'] for r in v):.4f} "
        f"[{min(r['besk'] for r in v):.4f}..{max(r['besk'] for r in v):.4f}] | "
        f"train {st.mean(r['train_seed'] for r in v):.3f}->{st.mean(r['train_new'] for r in v):.3f} | "
        f"held5 mean {st.mean(r['held_seed'] for r in v):+.3f}->{st.mean(r['held_new'] for r in v):+.3f} "
        f"(med {st.median(r['held_seed'] for r in v):+.3f}->{st.median(r['held_new'] for r in v):+.3f}) | "
        f"ctr_ho mean {st.mean(r['centre_held_seed'] for r in v):+.3f}->{st.mean(r['centre_held_new'] for r in v):+.3f} | "
        f"improved {sum(r['held_new'] > r['held_seed'] for r in v)}/{len(v)} | "
        f"dLOD {st.mean(r['dlod1'] for r in v):.4f} | agree {st.median(r['mean_agree'] for r in v):+.2f}"
    )

fig, axs = plt.subplots(1, 2, figsize=(15, 5.2), facecolor="#fcfcfb")
lons = [r["lon"] for r in rows]
lats = [r["lat"] for r in rows]
for ax, key, title, cmap, vmin, vmax in [
    (axs[0], "held_new", "5-leg mean held-out CC (2000-2005)", "RdBu_r", -0.8, 0.8),
    (axs[1], "besk", "Final shared besk by basin", "viridis", 0.195, 0.236),
]:
    ax.set_facecolor("#fcfcfb")
    vals = [r[key] for r in rows]
    sc = ax.scatter(
        lons, lats, c=vals, cmap=cmap, vmin=vmin, vmax=vmax,
        s=460, marker="s", edgecolor="#222", lw=0.7,
    )
    for r in rows:
        lbl = f"{r[key]:+.2f}" if key == "held_new" else f"{r[key]:.3f}"
        ax.text(
            r["lon"], r["lat"], lbl, ha="center", va="center", fontsize=6.5, color="#111",
            bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.55, lw=0),
        )
    ax.set_xlim(-185, 185)
    ax.set_ylim(-75, 80)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-60, 81, 20))
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.set_title(title, loc="left", fontsize=11)
    ax.grid(True, color="#e4e3df", lw=0.5)
    fig.colorbar(sc, ax=ax, shrink=0.82)

fig.suptitle(
    "Kaplan 30-Quad Fast-Refine Sweep (Basin-Seeded, Shared Params, 2000-2005 Held Out)",
    fontsize=12, x=0.02, ha="left",
)
fig.tight_layout()
fig.savefig(HERE / "quad5_besk_summary.png", dpi=130, facecolor="#fcfcfb")
print("Saved:", HERE / "quad5_besk_summary.csv", "and", HERE / "quad5_besk_summary.png")
