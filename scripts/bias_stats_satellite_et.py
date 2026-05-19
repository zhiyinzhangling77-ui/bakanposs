"""
scripts/bias_stats_satellite_et.py
==================================

Section 10 of SATELLITE_ET_NOTES.md を再現する集計スクリプト。

入力:
  --master  data/master_full_v2.csv
            列: date, site, ET_mm (EC), ET_metv3_mm (sat), metv3_qflag, Ta_min, ...

処理 (fair comparison, per site):
  1. qflag フィルタ:    metv3_qflag <= QFLAG_MAX (default 1)
  2. 雪/凍結日除外:     Ta_min > FREEZE_TA_MIN_THRESHOLD (default 0 °C)
  3. 両方が値を持つ日のみ (paired comparison)
  4. 統計量を算出:      n, mean_ec, mean_sat, bias, bias_pct,
                         rmse, r (Pearson), KGE
  5. CSV と scatter 図 (pair plot) を出力

出力: ./output_bias_stats/
  bias_stats_summary.csv         サイト別 1 行
  fig_scatter_ec_vs_metv3.png    poster Methods 用 scatter

使い方:
  python3 scripts/bias_stats_satellite_et.py
  python3 scripts/bias_stats_satellite_et.py --qflag-max 1 --freeze-ta 0
"""

from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["axes.unicode_minus"] = False

SITES = ["Oran", "Tarazona"]
SITE_COL = {"Oran": "#E85D04", "Tarazona": "#1D9E75"}


def parse_args():
    p = argparse.ArgumentParser(description="EC vs METv3 bias stats")
    p.add_argument("--master", default="data/master_full_v2.csv")
    p.add_argument("--qflag-max", type=float, default=1.0,
                   help="metv3_qflag の上限 (default 1 = good+nominal)")
    p.add_argument("--freeze-ta", type=float, default=0.0,
                   help="Ta_min がこの値以下の日を除外 [°C]")
    p.add_argument("--out", default="output_bias_stats")
    return p.parse_args()


def kge(obs: np.ndarray, sim: np.ndarray) -> float:
    """Kling-Gupta Efficiency.  r, α, β の3成分."""
    o, s = obs, sim
    if len(o) < 3 or np.std(o) == 0 or np.std(s) == 0:
        return np.nan
    r = float(np.corrcoef(o, s)[0, 1])
    alpha = float(np.std(s) / np.std(o))
    beta  = float(np.mean(s) / np.mean(o)) if np.mean(o) != 0 else np.nan
    if np.isnan(beta):
        return np.nan
    return float(1 - np.sqrt((r - 1) ** 2 + (alpha - 1) ** 2 + (beta - 1) ** 2))


def stats_per_site(df: pd.DataFrame, site: str,
                    qflag_max: float, freeze_ta: float) -> dict:
    sub = df[df["site"] == site].copy()
    n_raw = len(sub)

    sub = sub[sub["metv3_qflag"] <= qflag_max]
    n_after_qflag = len(sub)

    if "Ta_min" in sub.columns:
        sub = sub[sub["Ta_min"].fillna(99) > freeze_ta]
    n_after_freeze = len(sub)

    sub = sub.dropna(subset=["ET_mm", "ET_metv3_mm"])
    n_paired = len(sub)

    if n_paired < 5:
        return dict(site=site, n_raw=n_raw, n_after_qflag=n_after_qflag,
                      n_after_freeze=n_after_freeze, n_paired=n_paired,
                      valid=False, reason="too few paired days")

    ec  = sub["ET_mm"].values
    sat = sub["ET_metv3_mm"].values

    mean_ec  = float(np.mean(ec))
    mean_sat = float(np.mean(sat))
    bias     = mean_sat - mean_ec
    bias_pct = 100.0 * bias / mean_ec if mean_ec != 0 else np.nan
    rmse     = float(np.sqrt(np.mean((sat - ec) ** 2)))
    r        = float(np.corrcoef(ec, sat)[0, 1])
    kge_val  = kge(ec, sat)

    date_min = pd.to_datetime(sub["date"]).min().date()
    date_max = pd.to_datetime(sub["date"]).max().date()

    return dict(site=site, valid=True,
                  n_raw=n_raw, n_after_qflag=n_after_qflag,
                  n_after_freeze=n_after_freeze, n_paired=n_paired,
                  date_min=str(date_min), date_max=str(date_max),
                  mean_ec_mmd=round(mean_ec, 3),
                  mean_metv3_mmd=round(mean_sat, 3),
                  bias_mmd=round(bias, 3),
                  bias_pct=round(bias_pct, 1),
                  rmse_mmd=round(rmse, 3),
                  pearson_r=round(r, 3),
                  kge=round(kge_val, 3))


def draw_scatter(df: pd.DataFrame, results: list[dict], out_png: Path,
                  qflag_max: float, freeze_ta: float):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), facecolor="white")
    for ax, site in zip(axes, SITES):
        sub = df[df["site"] == site].copy()
        sub = sub[sub["metv3_qflag"] <= qflag_max]
        if "Ta_min" in sub.columns:
            sub = sub[sub["Ta_min"].fillna(99) > freeze_ta]
        sub = sub.dropna(subset=["ET_mm", "ET_metv3_mm"])
        ec  = sub["ET_mm"].values
        sat = sub["ET_metv3_mm"].values

        ax.scatter(ec, sat, s=14, alpha=0.45, color=SITE_COL[site],
                     edgecolor="none")
        lim = max(ec.max() if len(ec) else 1,
                    sat.max() if len(sat) else 1) * 1.05
        ax.plot([0, lim], [0, lim], color="#444", lw=1, ls="--", alpha=0.8)

        res = next((r for r in results if r["site"] == site), None)
        if res and res["valid"]:
            txt = (f"n = {res['n_paired']}\n"
                     f"mean EC  = {res['mean_ec_mmd']:.2f} mm/d\n"
                     f"mean Sat = {res['mean_metv3_mmd']:.2f} mm/d\n"
                     f"bias = {res['bias_mmd']:+.2f}  ({res['bias_pct']:+.0f} %)\n"
                     f"RMSE = {res['rmse_mmd']:.2f} mm/d\n"
                     f"r = {res['pearson_r']:.2f},  KGE = {res['kge']:.2f}")
            ax.text(0.04, 0.96, txt, transform=ax.transAxes,
                      fontsize=9, va="top", ha="left", family="monospace",
                      bbox=dict(boxstyle="round,pad=0.4", fc="white",
                                  ec="#888", alpha=0.92))

        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.set_xlabel("EC tower ET [mm day$^{-1}$]", fontsize=11)
        ax.set_ylabel("Meteosat ETv3 [mm day$^{-1}$]", fontsize=11)
        ax.set_title(site, fontsize=12, weight="bold", loc="left")
        ax.grid(alpha=0.25)

    fig.suptitle(
        f"EC vs Meteosat ETv3 (qflag ≤ {qflag_max:.0f}, Ta_min > {freeze_ta:.0f} °C)",
        fontsize=12, weight="bold", y=0.99)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  [save] {out_png}")


def main():
    args = parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("EC vs Meteosat ETv3 — bias statistics")
    print(f"  qflag_max = {args.qflag_max}, freeze_ta = {args.freeze_ta} °C")
    print(f"  master    = {args.master}")
    print(f"  out       = {out}")
    print("=" * 60)

    df = pd.read_csv(args.master)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    print(f"  raw rows: {len(df):,}")

    needed = ["ET_mm", "ET_metv3_mm", "metv3_qflag"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise SystemExit(f"[FATAL] master file missing columns: {missing}")

    results = []
    for s in SITES:
        r = stats_per_site(df, s, args.qflag_max, args.freeze_ta)
        results.append(r)
        if r["valid"]:
            print(f"\n  {s}:")
            print(f"    n_raw → qflag → freeze → paired :  "
                    f"{r['n_raw']} → {r['n_after_qflag']} → "
                    f"{r['n_after_freeze']} → {r['n_paired']}")
            print(f"    period   : {r['date_min']} → {r['date_max']}")
            print(f"    mean EC  = {r['mean_ec_mmd']:.2f} mm/d")
            print(f"    mean Sat = {r['mean_metv3_mmd']:.2f} mm/d")
            print(f"    bias     = {r['bias_mmd']:+.2f} mm/d "
                    f"({r['bias_pct']:+.0f} %)")
            print(f"    RMSE     = {r['rmse_mmd']:.2f} mm/d")
            print(f"    Pearson r= {r['pearson_r']:.2f}")
            print(f"    KGE      = {r['kge']:.2f}")
        else:
            print(f"\n  {s}: INVALID ({r['reason']})")

    summary = pd.DataFrame(results)
    csv_fp = out / "bias_stats_summary.csv"
    summary.to_csv(csv_fp, index=False)
    print(f"\n  [save] {csv_fp}")

    fig_fp = out / "fig_scatter_ec_vs_metv3.png"
    draw_scatter(df, results, fig_fp, args.qflag_max, args.freeze_ta)

    print(f"\n[done] {out}/")


if __name__ == "__main__":
    main()
