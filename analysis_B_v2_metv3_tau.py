"""
analysis_B_v2_metv3_tau.py
==========================

解析B v2: LSA SAF Meteosat ETv3 (30-min) で 解析A の τ ≈ 3-4 d を独立検証する。

v1 (MOD16A2) の致命的欠点:
  MOD16 は 8-day composite。τ ~ 3 d の指数減衰は 8-day 窓内で平滑化されて
  原理的に観測不能 (v1 で実証: Tarazona R² = -0.82, τ pegged at floor)。

v2 が解決するもの:
  METv3 は 30-min 瞬時 ET → daily 化しても 48 個の半時間値の合計。
  forward-fill 不要、3-4 d スケールの decay を直接観測できる。
  かつ 2018-2024 フル期間カバー → Oran (2018-2020) と Tarazona (2020-2024) の
  両方を比較可能になり、v1 で抜けていた Oran 検証が復活。

入力:
  --metv3-csv : metv3_daily_all.csv  (scripts/load_metv3.py が事前生成)
                列: date, site, ET_mm, n_obs, qflag_mean
                site は "Oran" / "TzM" (TzM → Tarazona に内部マッピング)
  --tara-csv  : Tarazona EC daily CSV (要 Irrig_mm)
  --oran-csv  : Oran     EC daily CSV (要 Rain_mm)

  EC reference (解析A v27):
    Oran     : τ = 2.82 d (SE 0.90)
    Tarazona : τ = 3.36 d (SE 0.62)

処理:
  1) metv3_daily_all.csv 読込 → site rename → 単位はそのまま mm/day
  2) Water-input 読込 + event 検出 (v1 と同じ閾値)
  3) 指数減衰 fit (inline; v1 と同じロジック) + bootstrap
  4) τ_metv3 vs τ_EC 比較 (両サイト)

出力先: ./output_analysis_B_v2/
  v2_metv3_tau_summary.csv
  v2_metv3_daily_extract.csv          (date,site,ET,LE,Irrig/Rain,events flag)
  fig01_metv3_timeseries.png
  fig02_metv3_recovery.png
  fig03_metv3_tau_comparison.png

Usage:
  python3 analysis_B_v2_metv3_tau.py
  python3 analysis_B_v2_metv3_tau.py --quick
  # METv3 daily がまだ無い場合は事前に:
  #   python3 scripts/load_metv3.py
"""

from __future__ import annotations
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Reuse v1 fit + helpers (event detection, exponential fit, validation)
sys.path.insert(0, str(Path(__file__).parent))
from analysis_B_v1_mod16_tau import (
    SITES, EC_TAU_REFERENCE, LE_PER_MM,
    IRRIG_THRESHOLD, RAIN_THRESHOLD,
    load_tara_irrig, load_oran_rain,
    detect_events, estimate_le_inf, fit_with_validation,
    _find_col,
)

warnings.filterwarnings("ignore")


# ================================================================
# CLI
# ================================================================

# METv3 daily uses "TzM" for Tarazona
METV3_SITE_RENAME = {"TzM": "Tarazona", "Oran": "Oran"}


def parse_args():
    p = argparse.ArgumentParser(
        description="解析B v2 — Meteosat ETv3 (30-min) で τ 独立検証",
    )
    p.add_argument("--metv3-csv",
                   default="/home/shion-nagamine/bakanposs/metv3_daily_all.csv",
                   help="METv3 daily aggregated CSV "
                        "(scripts/load_metv3.py の出力)")
    p.add_argument("--tara-csv",
                   default="/home/shion-nagamine/Dataset/Eddy data in Spain/"
                           "Daily_Summary_Filtered_forPred_ActEne26.csv")
    p.add_argument("--oran-csv",
                   default="/home/shion-nagamine/Dataset/Eddy data in Spain/"
                           "Oran_EddyDaily_MASTER_2018_2020_correct.csv")
    p.add_argument("--out", default="./output_analysis_B_v2")
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--quick", action="store_true",
                   help="n_boot=300, faster sanity-check run")
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args()


# ================================================================
# METv3 daily loader
# ================================================================

def load_metv3_daily(csv_path: Path) -> pd.DataFrame:
    """Read scripts/load_metv3.py の出力 (date,site,ET_mm,n_obs,qflag_mean)
    → unified long DataFrame (date, site, ET_mm_per_day, LE_Wm2)."""
    if not csv_path.exists():
        sys.exit(f"[FATAL] METv3 daily CSV not found: {csv_path}\n"
                  f"  → 先に scripts/load_metv3.py を実行してください "
                  f"(~3.5h, 120k NetCDF を処理)")
    df = pd.read_csv(csv_path)
    print(f"opened {csv_path.name}")
    print(f"  raw rows: {len(df):,}")
    print(f"  raw cols: {list(df.columns)}")

    date_col = _find_col(df, ["date", "Date"])
    site_col = _find_col(df, ["site", "name"])
    et_col   = _find_col(df, ["ET_mm", "ET_mm_per_day", "ET"])
    if not all([date_col, site_col, et_col]):
        sys.exit(f"[FATAL] required cols missing: "
                 f"date={date_col}, site={site_col}, et={et_col}")

    df = df.rename(columns={date_col: "date",
                              site_col: "site_raw",
                              et_col: "ET_mm_per_day"})
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df["site"] = df["site_raw"].map(lambda x: METV3_SITE_RENAME.get(x, x))
    df = df[df["site"].isin(SITES)].dropna(subset=["date", "ET_mm_per_day"])

    out = pd.DataFrame({
        "date": df["date"].values,
        "site": df["site"].values,
        "ET_mm_per_day": pd.to_numeric(df["ET_mm_per_day"],
                                          errors="coerce").values,
    })
    out["LE_Wm2"] = out["ET_mm_per_day"] * LE_PER_MM
    out = (out.sort_values(["site", "date"])
              .drop_duplicates(subset=["site", "date"])
              .reset_index(drop=True))

    print(f"  retained: {len(out):,} (site,date) rows")
    for s in SITES:
        sub = out[out["site"] == s]
        if len(sub):
            print(f"    {s:9s}: {len(sub):,} days  "
                  f"({sub['date'].min().date()} → {sub['date'].max().date()})  "
                  f"mean ET = {sub['ET_mm_per_day'].mean():.2f} mm/d")
    return out


# ================================================================
# Per-site driver
# ================================================================

def fit_metv3_tau(daily, water_df, site, water_col, threshold,
                    n_boot=2000):
    if daily.empty:
        return dict(site=site, valid=False, reason="empty METv3 series")
    df = (daily[daily["site"] == site][["date", "LE_Wm2"]]
            .rename(columns={"LE_Wm2": "LE"}))
    if water_df is not None and not water_df.empty:
        df = df.merge(water_df, on="date", how="left")
    if water_col not in df.columns:
        df[water_col] = 0.0
    df[water_col] = df[water_col].fillna(0)

    events = detect_events(df, water_col=water_col, threshold=threshold,
                              min_window=4, max_window=14)
    if not events:
        return dict(site=site, valid=False, reason="no events found", n_events=0)
    le_inf = estimate_le_inf(events)
    res = fit_with_validation(events, le_inf, n_boot=n_boot)
    res.update(site=site, n_events=len(events), le_inf=le_inf,
                  water_col=water_col, threshold=threshold)
    return res


# ================================================================
# Figures
# ================================================================

def fig_timeseries(daily, water_dfs, out_path):
    fig, axes = plt.subplots(2, 1, figsize=(13, 6), sharex=True)
    for ax, site in zip(axes, SITES):
        sub = daily[daily["site"] == site].sort_values("date")
        ax.plot(sub["date"], sub["ET_mm_per_day"], "-", color="tab:orange",
                  lw=0.7, label=f"METv3 ET ({site})")
        ax.set_ylabel("ET [mm/day]")
        ax.set_title(f"{site}: Meteosat ETv3 (30-min → daily aggregated)")
        ax.grid(alpha=0.3)
        wcol = "Irrig_mm" if site == "Tarazona" else "Rain_mm"
        thr  = IRRIG_THRESHOLD if site == "Tarazona" else RAIN_THRESHOLD
        wdf = water_dfs.get(site)
        if wdf is not None and not wdf.empty and wcol in wdf.columns:
            ev = wdf[wdf[wcol] > thr]
            ax2 = ax.twinx()
            ax2.bar(ev["date"], ev[wcol], width=1, alpha=0.4,
                       color="tab:blue", label=wcol)
            ax2.set_ylabel(f"{wcol} [mm/day]", color="tab:blue")
        ax.legend(loc="upper left")
    axes[-1].set_xlabel("date")
    fig.suptitle("METv3 ET time series with water events", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def fig_recovery(daily, water_dfs, fit_results, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, site in zip(axes, ["Oran", "Tarazona"]):
        res = fit_results.get(site, {})
        wcol = "Irrig_mm" if site == "Tarazona" else "Rain_mm"
        thr  = IRRIG_THRESHOLD if site == "Tarazona" else RAIN_THRESHOLD
        df = (daily[daily["site"] == site][["date", "LE_Wm2"]]
                .rename(columns={"LE_Wm2": "LE"}))
        wdf = water_dfs.get(site)
        if wdf is not None and not wdf.empty:
            df = df.merge(wdf, on="date", how="left")
        if wcol not in df.columns:
            df[wcol] = 0.0
        df[wcol] = df[wcol].fillna(0)
        events = detect_events(df, water_col=wcol, threshold=thr,
                                  min_window=4, max_window=14)
        for ev in events:
            evd = ev["data"]
            ax.plot(evd["days_since_event"], evd["LE"], "o", ms=3,
                      alpha=0.3, color="gray")
        if res.get("valid") and np.isfinite(res.get("tau", np.nan)):
            tau, le0, le_inf = res["tau"], res["le0"], res["le_inf"]
            xx = np.linspace(0, 14, 200)
            yy = le_inf + (le0 - le_inf) * np.exp(-xx / tau)
            ax.plot(xx, yy, "-", color="tab:red", lw=2,
                      label=(f"fit: τ={tau:.2f} d, LE0={le0:.1f}, "
                             f"LE∞={le_inf:.1f}\n"
                             f"amp={res['amplitude']:.1f} W/m², "
                             f"n_ev={res['n_events']}, R²={res['r2']:.2f}"))
        else:
            ax.text(0.05, 0.95,
                      f"[no valid fit]\n{res.get('reason', '?')}",
                      transform=ax.transAxes, va="top")
        ax.set_xlabel("days since water event")
        ax.set_ylabel("LE [W/m²] (METv3)")
        ax.set_title(site)
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right")
    fig.suptitle("METv3 LE recovery curves", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def fig_tau_compare(sat_fits, ec_taus, out_path):
    fig, ax = plt.subplots(figsize=(8, 5))
    sites = list(SITES)
    x = np.arange(len(sites))
    w = 0.35
    ec_v = [ec_taus.get(s, {}).get("tau", np.nan) for s in sites]
    ec_e = [ec_taus.get(s, {}).get("se",  np.nan) for s in sites]
    sat_v = [sat_fits.get(s, {}).get("tau",   np.nan) for s in sites]
    sat_e = [sat_fits.get(s, {}).get("tau_se", np.nan) for s in sites]
    ax.bar(x - w/2, ec_v, w, yerr=ec_e, capsize=4,
            color="tab:blue", alpha=0.8, label="EC tower (analysis A v27)")
    ax.bar(x + w/2, sat_v, w, yerr=sat_e, capsize=4,
            color="tab:orange", alpha=0.8, label="METv3 satellite (this work)")
    for i, (e, s) in enumerate(zip(ec_v, sat_v)):
        if np.isfinite(e):
            ax.text(i - w/2, e + 0.3, f"{e:.2f}", ha="center", fontsize=10)
        if np.isfinite(s):
            ax.text(i + w/2, s + 0.3, f"{s:.2f}", ha="center", fontsize=10)
    ax.set_xticks(x); ax.set_xticklabels(sites)
    ax.set_ylabel("τ [days]")
    ax.set_title("Recovery time τ : EC tower vs Meteosat ETv3")
    ax.axhspan(3.0, 4.0, color="green", alpha=0.1,
                 label="analysis A target: τ ≈ 3-4 d")
    ax.legend()
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ================================================================
# MAIN
# ================================================================

def main():
    args = parse_args()
    n_boot = 300 if args.quick else args.n_boot
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) METv3 daily → site x date long table
    print("\n=== 1) reading METv3 daily ===")
    daily = load_metv3_daily(Path(args.metv3_csv))

    # 2) water-input series
    print("\n=== 2) loading water-input series ===")
    tara_water = load_tara_irrig(Path(args.tara_csv))
    oran_water = load_oran_rain(Path(args.oran_csv))
    print(f"  Tarazona Irrig events (>{IRRIG_THRESHOLD}): "
          f"{(tara_water.get('Irrig_mm', pd.Series([], dtype=float)) > IRRIG_THRESHOLD).sum()}")
    print(f"  Oran     Rain  events (>{RAIN_THRESHOLD}):  "
          f"{(oran_water.get('Rain_mm', pd.Series([], dtype=float)) > RAIN_THRESHOLD).sum()}")
    water_dfs = {"Tarazona": tara_water, "Oran": oran_water}

    # 3) τ-fit per site
    print("\n=== 3) fitting τ on METv3 LE ===")
    sat_fits = {}
    for site in SITES:
        wcol = "Irrig_mm" if site == "Tarazona" else "Rain_mm"
        thr  = IRRIG_THRESHOLD if site == "Tarazona" else RAIN_THRESHOLD
        print(f"  -- {site} ({wcol} > {thr}) --")
        res = fit_metv3_tau(daily, water_dfs[site], site, wcol, thr,
                              n_boot=n_boot)
        sat_fits[site] = res
        if res.get("valid"):
            print(f"    τ = {res['tau']:.2f} d  "
                   f"CI [{res['ci_lo']:.2f}, {res['ci_hi']:.2f}]  "
                   f"amp = {res['amplitude']:.1f} W/m²  "
                   f"R² = {res['r2']:.3f}  "
                   f"n_events = {res['n_events']}")
        else:
            print(f"    [invalid] {res.get('reason', 'unknown')}")

    # 4) EC reference + summary
    print("\n=== 4) EC reference (analysis A v27) ===")
    for s, v in EC_TAU_REFERENCE.items():
        print(f"  {s}: τ = {v['tau']:.2f} d (SE {v['se']:.2f}) — {v['source']}")

    rows = []
    for site in SITES:
        sf = sat_fits[site]
        ef = EC_TAU_REFERENCE[site]
        row = dict(
            site=site,
            tau_ec=ef["tau"], se_ec=ef["se"],
            tau_metv3=sf.get("tau", np.nan),
            tau_se_metv3=sf.get("tau_se", np.nan),
            ci_lo_metv3=sf.get("ci_lo", np.nan),
            ci_hi_metv3=sf.get("ci_hi", np.nan),
            n_events_metv3=sf.get("n_events", 0),
            amp_metv3=sf.get("amplitude", np.nan),
            le_inf_metv3=sf.get("le_inf", np.nan),
            r2_metv3=sf.get("r2", np.nan),
            valid_metv3=sf.get("valid", False),
            reason=sf.get("reason", ""),
        )
        if sf.get("valid") and np.isfinite(sf.get("tau_se", np.nan)):
            mde = 1.96 * np.sqrt(ef["se"]**2 + sf["tau_se"]**2)
            row["MDE"] = mde
            row["diff"] = abs(ef["tau"] - sf["tau"])
            row["significant_diff"] = row["diff"] > mde
        else:
            row["MDE"] = np.nan
            row["diff"] = np.nan
            row["significant_diff"] = None
        rows.append(row)
    summary = pd.DataFrame(rows)

    summary.to_csv(out_dir / "v2_metv3_tau_summary.csv", index=False)
    daily.to_csv(out_dir / "v2_metv3_daily_extract.csv", index=False)

    pd.set_option("display.float_format", "{:+.3f}".format)
    print("\n=== SUMMARY ===")
    cols_show = ["site", "tau_ec", "se_ec", "tau_metv3", "tau_se_metv3",
                   "n_events_metv3", "amp_metv3", "r2_metv3",
                   "MDE", "diff", "significant_diff", "reason"]
    print(summary[cols_show].to_string(index=False))
    print(f"\nwrote {out_dir / 'v2_metv3_tau_summary.csv'}")

    if args.no_plots:
        print("[skipped plots --no-plots]")
        return

    print("\n=== 5) generating figures ===")
    fig_timeseries(daily, water_dfs, out_dir / "fig01_metv3_timeseries.png")
    print("  wrote fig01_metv3_timeseries.png")
    fig_recovery(daily, water_dfs, sat_fits,
                   out_dir / "fig02_metv3_recovery.png")
    print("  wrote fig02_metv3_recovery.png")
    fig_tau_compare(sat_fits, EC_TAU_REFERENCE,
                       out_dir / "fig03_metv3_tau_comparison.png")
    print("  wrote fig03_metv3_tau_comparison.png")

    print(f"\n[DONE] outputs in {out_dir.absolute()}")


if __name__ == "__main__":
    main()
