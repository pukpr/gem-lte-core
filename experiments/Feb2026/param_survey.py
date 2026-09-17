#!/usr/bin/env python3
"""param_survey.py — exploratory look at the regression-fit parameters
(ann1, ann2, sem1, sem2) plus IR, across a set of index subdirectories.

Scope note: of the "applied after the forcing manifold" group (level, k0,
per-mode amp/phase, trend, accel, ann1, ann2, sem1, sem2), only
ann1/ann2/sem1/sem2 are actually persisted in lt.exe.p — level/k0/trend/
accel/per-mode amp-phase are recomputed fresh from whatever CV window is
active on each forward() call and never written to disk. This script only
looks at what's actually on disk: ann1/ann2/sem1/sem2 + IR (IR isn't a
regression output — it's a free/searched Param_B scalar — but is included
per request since it's likewise applied to Model after Forcing exists).

Usage:
    ./param_survey.py
    ./param_survey.py nino4 pdo amo          # subset of indices
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lte_forward import Config, read_resp, forward  # noqa: E402

ROOT = Path(__file__).resolve().parent

DEFAULT_INDICES = ["nino4", "pdo", "amo", "tna", "nao", "baltic", "iode"]
PARAMS = ["ann1", "ann2", "sem1", "sem2", "IR"]

# Same project-wide defaults lte_run.py falls back to when an index's own
# lt.exe.resp doesn't set TRAIN_START/TRAIN_END/EXCLUDE.
CV_DEFAULT = ("1880", "1885")
EXCLUDE_DEFAULT = "true"


def load_params(indices: list[str]) -> dict[str, dict[str, float]]:
    data = {}
    for idx in indices:
        p_path = ROOT / idx / "lt.exe.p"
        if not p_path.exists():
            print(f"  (skipping {idx}: {p_path} not found)", file=sys.stderr)
            continue
        raw = json.loads(p_path.read_text())
        data[idx] = {k: float(raw[k]) for k in PARAMS if k in raw}
        missing = [k for k in PARAMS if k not in raw]
        if missing:
            print(f"  (warning: {idx} missing {missing})", file=sys.stderr)
    return data


def print_table(data: dict[str, dict[str, float]]) -> None:
    indices = list(data.keys())
    col_w = 13
    header = "index".ljust(10) + "".join(p.rjust(col_w) for p in PARAMS)
    print(header)
    print("-" * len(header))
    for idx in indices:
        row = idx.ljust(10)
        for p in PARAMS:
            v = data[idx].get(p)
            row += (f"{v:+.6f}".rjust(col_w) if v is not None
                    else "n/a".rjust(col_w))
        print(row)


def print_similarity_summary(data: dict[str, dict[str, float]]) -> None:
    print("\n--- cross-index spread per parameter (similar -> different) ---")
    stats = []
    for p in PARAMS:
        vals = np.array([data[idx][p] for idx in data if p in data[idx]])
        if len(vals) < 2:
            continue
        mean, std = vals.mean(), vals.std()
        scale = np.abs(vals).mean() + 1e-12  # avoid div-by-zero near 0
        rel_spread = std / scale
        stats.append((rel_spread, p, mean, std, vals.min(), vals.max()))
    stats.sort()  # most similar (lowest relative spread) first
    print(f"{'param':<8}{'rel.spread':>12}{'mean':>14}{'std':>14}"
          f"{'min':>14}{'max':>14}")
    for rel_spread, p, mean, std, vmin, vmax in stats:
        print(f"{p:<8}{rel_spread:>12.3f}{mean:>+14.6f}{std:>14.6f}"
              f"{vmin:>+14.6f}{vmax:>+14.6f}")
    print("\n(rel.spread = std / mean(|value|) across indices — lower means")
    print(" the parameter lands in a similar place for every index; higher")
    print(" means it's doing something more index-specific.)")


def print_outliers(data: dict[str, dict[str, float]], z_thresh: float = 1.5) -> None:
    print(f"\n--- indices that stand out per parameter (|z| > {z_thresh}) ---")
    any_flagged = False
    for p in PARAMS:
        items = [(idx, data[idx][p]) for idx in data if p in data[idx]]
        if len(items) < 3:
            continue
        vals = np.array([v for _, v in items])
        mean, std = vals.mean(), vals.std()
        if std == 0:
            continue
        for idx, v in items:
            z = (v - mean) / std
            if abs(z) > z_thresh:
                any_flagged = True
                print(f"  {p:<6} {idx:<10} value={v:+.6f}  z={z:+.2f}")
    if not any_flagged:
        print("  (none)")


def compute_winding(idx: str) -> dict | None:
    """Run forward() once for `idx` and pull out the Forcing->Model mapping's
    parameters: the winding frequencies (wave-numbers M — ltep's base modes
    plus their integer harmonics) with each one's fitted amplitude/phase,
    plus the "independent" (forcing-blind, calendar/date-only) terms:
    annual/semi-annual (ann1/ann2/sem1/sem2) and trend/accel. Returns None
    (with a printed warning) if the index can't be evaluated at all —
    exploratory tool, one bad index shouldn't kill the rest.
    """
    idx_dir = ROOT / idx
    p_path = idx_dir / "lt.exe.p"
    if not p_path.exists():
        return None
    try:
        params = json.loads(p_path.read_text())
        resp = read_resp(idx_dir / "lt.exe.resp")
        resp["CLIMATE_INDEX"] = str(idx_dir / f"{idx}.dat")
        overrides = {
            "TRAIN_START": resp.get("TRAIN_START", CV_DEFAULT[0]),
            "TRAIN_END": resp.get("TRAIN_END", CV_DEFAULT[1]),
            "EXCLUDE": resp.get("EXCLUDE", EXCLUDE_DEFAULT),
        }
        if "harm" in params:
            overrides["NH"] = " ".join(
                str(int(float(h))) for h in params["harm"])
        cfg = Config(resp, overrides)
        r = forward(cfg, params)
    except Exception as exc:  # noqa: BLE001 — exploratory, keep going
        print(f"  (skipping {idx} for winding analysis: {exc})",
              file=sys.stderr)
        return None

    semi1, semi2, ann1, ann2 = r["reg"]["annual"]
    return {
        "m": np.asarray(r["m"], dtype=float),
        "amp": np.asarray(r["reg"]["amp"], dtype=float),
        "phase": np.asarray(r["reg"]["phase"], dtype=float),
        "nm": r["nm"],
        "nh": r["nh"],
        "k0": r["reg"]["k0"],  # the "wavenumber=0" term — plotted at M=0
        "ann1": ann1, "ann2": ann2, "semi1": semi1, "semi2": semi2,
        "trend": r["reg"]["trend"], "accel": r["reg"]["accel"],
        "cc": r["cc"],
        "cv": (overrides["TRAIN_START"], overrides["TRAIN_END"]),
        "exclude": overrides["EXCLUDE"],
    }


def winding_xy(w: dict) -> tuple[np.ndarray, np.ndarray]:
    """(x, y) for plotting: winding frequencies as |M| (a winding number is
    a rate, not a signed direction — the negative-M and positive-M fits are
    the same physical winding), plus k0 — the "wavenumber=0" term, i.e. the
    plain-linear-in-Forcing response rather than an oscillatory one — added
    at x=0 alongside them, NOT absolute-valued away like the others (zero
    is already its own position, not a sign to discard)."""
    x = np.concatenate([np.abs(w["m"]), [0.0]])
    y = np.concatenate([np.abs(w["amp"]), [abs(w["k0"])]])
    return x, y


def print_winding_tables(winding: dict[str, dict]) -> None:
    print("\n=== Forcing -> Model mapping: winding frequencies (ltep + "
          "harmonics) ===")
    print("(fit window per index — resp's TRAIN_START/END if set, else "
          f"{CV_DEFAULT[0]}-{CV_DEFAULT[1]}; EXCLUDE from resp if set, "
          f"else {EXCLUDE_DEFAULT})\n")
    for idx, w in winding.items():
        print(f"-- {idx}  (cv={w['cv'][0]}-{w['cv'][1]}  "
              f"exclude={w['exclude']}  cc={w['cc']:.4f}, "
              f"{w['nm']} base mode(s) + {w['nh']} harmonic(s)) --")
        order = np.argsort(-np.abs(w["amp"]))
        for i in order:
            tag = "base" if i < w["nm"] else "harm"
            print(f"    [{tag}] M={w['m'][i]:+.6f}   "
                  f"amp={w['amp'][i]:+.6f}   phase={w['phase'][i]:+.6f}")
        print(f"    [k0, M=0]  k0={w['k0']:+.6f}")
        print(f"    [date-only] ann1={w['ann1']:+.6f}  ann2={w['ann2']:+.6f}"
              f"  sem1={w['semi1']:+.6f}  sem2={w['semi2']:+.6f}  "
              f"trend={w['trend']:+.6f}  accel={w['accel']:+.6f}")


def plot_periodogram(winding: dict[str, dict], out_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, (ax_raw, ax_norm) = plt.subplots(
        2, 1, figsize=(12, 10), sharex=True)
    colors = plt.cm.tab10(np.linspace(0, 1, max(len(winding), 1)))
    n = len(winding)

    for i, (idx, w) in enumerate(winding.items()):
        # Small per-index x-offset so overlapping winding frequencies
        # across indices stay visually distinguishable as separate stems
        # rather than perfectly overlapping.
        jitter = (i - (n - 1) / 2) * 0.0015
        x_raw, y = winding_xy(w)
        x = x_raw + jitter
        y_norm = y / (y.max() + 1e-12)  # each index's own peak -> 1.0

        for ax, yy in ((ax_raw, y), (ax_norm, y_norm)):
            ax.vlines(x, 0, yy, color=colors[i], linewidth=2, alpha=0.85,
                      label=idx if ax is ax_raw else None)
            ax.scatter(x, yy, color=colors[i], s=18, zorder=3)

    ax_raw.axhline(0, color="k", lw=0.6)
    ax_raw.set_ylabel("|amplitude| (raw)")
    ax_raw.set_title("Forcing -> Model winding-frequency periodogram, by "
                      "index — raw magnitude")
    ax_raw.legend(loc="upper right", ncol=2, fontsize=9)
    ax_raw.grid(True, ls=":", lw=0.4, alpha=0.6)

    ax_norm.axhline(0, color="k", lw=0.6)
    ax_norm.set_xlabel("Winding frequency |M| (ltep base modes + integer "
                        "harmonics; k0, the wavenumber=0 term, at x=0)")
    ax_norm.set_xlim(left=0)
    ax_norm.set_ylabel("|amplitude| / index's own max")
    ax_norm.set_title("Same data, normalized per index — for comparing "
                       "SHAPE regardless of scale (top panel can be "
                       "dominated by one index's current absolute scale)")
    ax_norm.grid(True, ls=":", lw=0.4, alpha=0.6)

    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def load_locations(indices: list[str]) -> dict[str, tuple[float, float]]:
    """Read (lat, lon) in signed decimal degrees for each index from
    ID.yml (the same file lte_gui.py uses for its Name/Country fields)."""
    import yaml
    id_path = ROOT / "ID.yml"
    raw = yaml.safe_load(id_path.read_text())
    locs = {}
    for idx in indices:
        rec = raw.get(idx)
        if rec is None:
            print(f"  (no ID.yml entry for {idx}, skipping on the map)",
                  file=sys.stderr)
            continue
        lat = float(rec["Latitude"]) * (-1.0 if rec["N/S"] == "S" else 1.0)
        lon = float(rec["Longitude"]) * (-1.0 if rec["E/W"] == "W" else 1.0)
        locs[idx] = (lat, lon)
    return locs


def plot_geo_fingerprints(winding: dict[str, dict], out_path: Path) -> None:
    """Each index's winding-frequency fingerprint (|amplitude| vs M,
    normalized to that index's own peak — comparing SHAPE, not the
    currently-uneven absolute scale) as a small bar chart placed at its
    real-world location on a Mercator map.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    locs = load_locations(list(winding.keys()))
    if not locs:
        print("(no locations available, skipping geo-fingerprint map)",
              file=sys.stderr)
        return

    proj = ccrs.Mercator()
    fig = plt.figure(figsize=(16, 9))
    ax = plt.axes(projection=proj)
    ax.set_extent([-180, 180, -40, 75], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="0.92")
    ax.add_feature(cfeature.OCEAN, facecolor="0.98")
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, alpha=0.4)
    ax.gridlines(draw_labels=False, linewidth=0.3, alpha=0.3)
    ax.set_title("Winding-frequency fingerprint by index location "
                  "(x = |winding frequency|, k0 at x=0; "
                  "bar height = |amplitude| / that index's own peak)")

    colors = plt.cm.tab10(np.linspace(0, 1, max(len(winding), 1)))
    inset_w, inset_h = 0.13, 0.16  # figure-fraction size of each mini-chart

    # Shared winding-frequency (x) axis across every inset, so position
    # along x is directly comparable index to index — only the amplitude
    # (y) stays independently normalized per index, since that depends on
    # each climate index's own physical scale, not something to compare
    # directly. Winding frequency is a rate, not a signed direction, so
    # it's plotted as |M| — never negative; k0 (the wavenumber=0, plain-
    # linear-in-Forcing term) sits at x=0 rather than being dropped.
    all_m = np.concatenate([winding_xy(w)[0] for w in winding.values()])
    m_pad = all_m.max() * 0.05 + 1e-6
    shared_xlim = (0.0, all_m.max() + m_pad)  # never negative
    bar_width = (all_m.max() + 1e-6) * 0.02

    # True (unmoved) marker positions, in figure-fraction coordinates —
    # kept separate from the (possibly nudged) inset centers below so a
    # leader line can always point back to the real location.
    idx_order = list(locs.keys())
    true_xy = np.zeros((len(idx_order), 2))
    for i, idx in enumerate(idx_order):
        lat, lon = locs[idx]
        x_data, y_data = proj.transform_point(lon, lat, ccrs.PlateCarree())
        ax.plot(x_data, y_data, marker="o", markersize=5,
                color=colors[i], transform=None, zorder=5)
        x_disp, y_disp = ax.transData.transform((x_data, y_data))
        true_xy[i] = fig.transFigure.inverted().transform((x_disp, y_disp))

    # Simple iterative pairwise repulsion so insets whose true locations are
    # close together (e.g. amo/tna/nao/baltic all sit in/around the N.
    # Atlantic) don't render on top of each other — push apart only enough
    # to clear inset_w x inset_h, symmetrically, a few passes.
    inset_xy = true_xy.copy()
    min_dx, min_dy = inset_w * 1.05, inset_h * 1.05
    for _ in range(200):
        moved = False
        for i in range(len(idx_order)):
            for j in range(i + 1, len(idx_order)):
                dx = inset_xy[j, 0] - inset_xy[i, 0]
                dy = inset_xy[j, 1] - inset_xy[i, 1]
                overlap_x = min_dx - abs(dx)
                overlap_y = min_dy - abs(dy)
                if overlap_x > 0 and overlap_y > 0:
                    moved = True
                    # Push apart along whichever axis has less overlap to
                    # resolve, so pairs mostly separate vertically or
                    # horizontally rather than diagonally drifting far.
                    if overlap_x < overlap_y:
                        shift = (overlap_x / 2 + 1e-4) * (1 if dx >= 0 else -1)
                        inset_xy[j, 0] += shift
                        inset_xy[i, 0] -= shift
                    else:
                        shift = (overlap_y / 2 + 1e-4) * (1 if dy >= 0 else -1)
                        inset_xy[j, 1] += shift
                        inset_xy[i, 1] -= shift
        if not moved:
            break
    # Keep every inset's center within the figure so none run off the edge.
    margin_x, margin_y = inset_w / 2 + 0.01, inset_h / 2 + 0.01
    inset_xy[:, 0] = np.clip(inset_xy[:, 0], margin_x, 1 - margin_x)
    inset_xy[:, 1] = np.clip(inset_xy[:, 1], margin_y, 1 - margin_y)

    for i, idx in enumerate(idx_order):
        w = winding[idx]
        x_fig, y_fig = inset_xy[i]
        tx_fig, ty_fig = true_xy[i]
        if np.hypot(x_fig - tx_fig, y_fig - ty_fig) > 1e-3:
            fig.add_artist(Line2D(
                [tx_fig, x_fig], [ty_fig, y_fig], transform=fig.transFigure,
                color=colors[i], linewidth=0.7, alpha=0.6, zorder=4))

        iax = fig.add_axes(
            [x_fig - inset_w / 2, y_fig - inset_h / 2, inset_w, inset_h])
        m, y = winding_xy(w)
        y = y / (y.max() + 1e-12)
        iax.bar(m, y, width=bar_width,
                color=colors[i], edgecolor="k", linewidth=0.3)
        iax.set_xlim(*shared_xlim)
        iax.set_ylim(0, 1.15)
        iax.text(0.5, 0.95, idx, transform=iax.transAxes,
                 ha="center", va="top", fontsize=9, fontweight="bold",
                 bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                           edgecolor="none", alpha=0.75))
        iax.tick_params(labelsize=6)
        iax.set_yticks([])
        for spine in ("top", "right"):
            iax.spines[spine].set_visible(False)
        for spine in iax.spines.values():
            spine.set_edgecolor(colors[i])

    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def main() -> int:
    indices = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_INDICES
    data = load_params(indices)
    if not data:
        print("No lt.exe.p files found.", file=sys.stderr)
        return 2

    print(f"Loaded {len(data)}/{len(indices)} indices: {', '.join(data)}\n")
    print_table(data)
    print_similarity_summary(data)
    print_outliers(data)

    print("\nNote: level, k0, trend, accel, and per-mode amp/phase are NOT")
    print("in lt.exe.p — they're recomputed fresh per forward() call from")
    print("whichever CV window is active, never persisted to disk. The")
    print("winding-frequency analysis below recomputes them via one")
    print("forward() pass per index (see the window/exclude printed for")
    print("each) — those numbers will move if you re-run with a different")
    print("--cv, unlike the ann1/ann2/sem1/sem2/IR table above.")

    winding = {}
    for idx in data:  # only indices that loaded successfully above
        w = compute_winding(idx)
        if w is not None:
            winding[idx] = w
    if winding:
        print_winding_tables(winding)
        out_path = ROOT / "winding_periodogram.png"
        plot_periodogram(winding, out_path)
        print(f"\nSaved periodogram chart to {out_path}")

        geo_path = ROOT / "geo_fingerprints.png"
        plot_geo_fingerprints(winding, geo_path)
        print(f"Saved geographic fingerprint map to {geo_path}")
    else:
        print("\n(no indices could be evaluated for the winding analysis)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
