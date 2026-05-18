"""
analysis_B_v5_oran_lst_response.py
====================================

解析B v5: Oran rainfed で METv3 LE が雨直後にスパイクするか直接観測する。

v3 で発見した不一致:
  τ_bias (EC − METv3) = 0.83 d  <<  τ_EC = 2.82 d
仮説:
  雨 → 表層湿潤 → LST 冷却 (蒸発冷却) → METv3 アルゴリズムが ET 上昇
  → satellite ET が短時間スパイク → bias が EC より速く decay
TzM drip では:
  表層 5 cm は乾いたまま (wet-bulb at 10-30 cm) → LST 不変
  → satellite flat → bias decay = EC decay shape

このスクリプトでは:
  1) Oran で各 rain event 後 0-10 d の METv3 LE 単体を per-event 重ね plot
  2) 同じ window の EC LE を並べる
  3) METv3 LE と EC LE それぞれに exp decay fit を試みる
     - METv3 LE が雨後 spike + fast decay → 仮説支持
     - METv3 LE が flat → 別の bias mechanism
  4) TzM Tarazona で同じ図を作って物理コントラストを示す

入力:
  --metv3-csv : metv3_daily_all.csv
  --tara-csv  : Tarazona EC daily
  --oran-csv  : Oran     EC daily

出力先: ./output_analysis_B_v5/
  v5_oran_lst_response.csv
  fig01_oran_metv3_vs_ec_recovery.png  (per-event panel)
  fig02_tara_metv3_vs_ec_recovery.png
  fig03_contrast_decay_speeds.png      (oran fast vs tara slow)

Usage:
  python3 analysis_B_v5_oran_lst_response.py
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
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).parent))
from analysis_B_v1_mod16_tau import (
    SITES, LE_PER_MM,
    IRRIG_THRESHOLD, RAIN_THRESHOLD,
)
from analysis_B_v2_metv3_tau import load_metv3_daily
from analysis_B_v3_bias_tau import load_tara_ec, load_oran_ec

warnings.filterwarnings("ignore")


TAU_FLOOR = 0.5
TAU_CEIL  = 60.0


def parse_args():
    p = argparse.ArgumentParser(
        description="解析B v5 — Oran rain → METv3 LST cooling response",
    )
    p.add_argument("--metv3-csv",
                   default="/home/shion-nagamine/bakanposs/metv3_daily_all.csv")
    p.add_argument("--tara-csv",
                   default="/home/shion-nagamine/Dataset/Eddy data in Spain/"
                           "Daily_Summary_Filtered_forPred_ActEne26.csv")
    p.add_argument("--oran-csv",
                   default="/home/shion-nagamine/Dataset/Eddy data in Spain/"
                           "Oran_EddyDaily_MASTER_2018_2020_correct.csv")
    p.add_argument("--out", default="./output_analysis_B_v5")
    p.add_argument("--max-window", type=int, default=10)
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args()


# ================================================================
# Build per-event windows for both EC LE and METv3 LE
# ================================================================

def build_events_double(ec_df: pd.DataFrame, metv3_df: pd.DataFrame,
                         site: str, water_col: str, threshold: float,
                         max_window: int = 10):
    """Return long DataFrame: event_start, days_since_event, LE_EC, LE_METv3."""
    metv3_site = (metv3_df[metv3_df["site"] == site][["date", "LE_Wm2"]]
                       .rename(columns={"LE_Wm2": "LE_METv3"}))
    merged = ec_df.merge(metv3_site, on="date", how="inner").rename(
        columns={"LE_EC_Wm2": "LE_EC"})
    event_dates = merged[merged[water_col].fillna(0) > threshold]["date"].tolist()
    # remove events within max_window of next (avoid overlap)
    keep = []
    last = None
    for d in event_dates:
        d = pd.Timestamp(d)
        if last is None or (d - last).days >= 4:
            keep.append(d)
            last = d
    pieces = []
    for ed in keep:
        end = ed + pd.Timedelta(days=max_window)
        ev = merged[(merged["date"] >= ed) & (merged["date"] <= end)].copy()
        ev["days_since_event"] = (ev["date"] - ed).dt.days
        ev["event_start"] = ed
        ev["pulse_mm"] = float(merged.loc[
            merged["date"] == ed, water_col].iloc[0])
        pieces.append(ev[["event_start", "date", "days_since_event",
                              "LE_EC", "LE_METv3", "pulse_mm"]])
    if not pieces:
        return pd.DataFrame()
    return pd.concat(pieces, ignore_index=True)


# ================================================================
# Recovery fit (2-param exp with LE_inf as 7-d tail median)
# ================================================================

def exp_2p(d, le0, tau, le_inf):
    return le_inf + (le0 - le_inf) * np.exp(-d / tau)


def fit_recovery(arr_d, arr_y):
    """Median-per-day NLS fit. Returns (le0, tau, le_inf, r2) or None."""
    df = pd.DataFrame({"d": arr_d, "y": arr_y}).dropna()
    if df.empty:
        return None
    g = (df.groupby("d")["y"].agg(["median", "count"]).reset_index())
    g = g[g["count"] >= 3]
    if len(g) < 3:
        return None
    x = g["d"].values.astype(float)
    y = g["median"].values
    tail = df[df["d"] >= 7]["y"]
    le_inf = float(np.median(tail)) if len(tail) >= 5 else float(
        np.percentile(df["y"], 25))
    le0_floor = max(le_inf + 1, 10.0)
    le0_ceil  = max(y.max() * 1.5, 350.0)

    def model(d, le0, tau):
        return exp_2p(d, le0, tau, le_inf)

    for p0 in [[y.max(), 3.0], [y.max() * 1.1, 5.0],
                 [(y.max() + le_inf) / 2, 2.0]]:
        try:
            popt, _ = curve_fit(model, x, y, p0=p0,
                                  bounds=([le0_floor, TAU_FLOOR],
                                          [le0_ceil, TAU_CEIL]),
                                  maxfev=8000)
            yp = model(x, *popt)
            r2 = (1 - np.sum((y - yp)**2) / np.sum((y - y.mean())**2)
                  if np.var(y) > 0 else np.nan)
            return dict(le0=popt[0], tau=popt[1], le_inf=le_inf,
                          r2=float(r2), amplitude=popt[0] - le_inf,
                          n_bins=len(g))
        except Exception:
            continue
    return None


# ================================================================
# Figures
# ================================================================

def fig_recovery_panel(events_long: pd.DataFrame, site: str,
                        fit_ec: dict | None, fit_metv3: dict | None,
                        out_path: Path):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    if events_long.empty:
        for ax in axes:
            ax.text(0.5, 0.5, "no events", ha="center", transform=ax.transAxes)
        fig.savefig(out_path, dpi=150)
        plt.close(fig)
        return

    for ax, col, label, c, fit in zip(
            axes,
            ["LE_EC", "LE_METv3"],
            ["EC LE", "METv3 LE"],
            ["tab:blue", "tab:orange"],
            [fit_ec, fit_metv3]):
        for ed, sub in events_long.groupby("event_start"):
            ax.plot(sub["days_since_event"], sub[col], "-o", ms=3,
                      lw=0.7, alpha=0.35, color=c)
        med = (events_long.groupby("days_since_event")[col]
                       .median().reset_index())
        ax.plot(med["days_since_event"], med[col], "o", ms=8,
                  color=c, markeredgecolor="k", label="median per t")
        if fit and np.isfinite(fit.get("tau", np.nan)):
            xx = np.linspace(0, events_long["days_since_event"].max(), 200)
            yy = exp_2p(xx, fit["le0"], fit["tau"], fit["le_inf"])
            ax.plot(xx, yy, "-", color="tab:red", lw=2,
                      label=(f"fit: τ={fit['tau']:.2f} d, "
                             f"LE0={fit['le0']:.1f}, LE∞={fit['le_inf']:.1f}\n"
                             f"amp={fit['amplitude']:.1f}, R²={fit['r2']:.2f}"))
        ax.set_xlabel("days since event")
        ax.set_ylabel(f"{label}  [W/m²]")
        ax.set_title(f"{site}: {label}")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right", fontsize=8)
    fig.suptitle(f"{site}: per-event EC vs METv3 LE recovery",
                  fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def fig_contrast(summary: pd.DataFrame, out_path: Path):
    """Side-by-side bar: tau_EC vs tau_METv3 for each site."""
    fig, ax = plt.subplots(figsize=(8, 5))
    sites = summary["site"].tolist()
    x = np.arange(len(sites))
    w = 0.35
    tau_ec    = summary["tau_EC"].values
    tau_metv3 = summary["tau_METv3"].values
    ax.bar(x - w/2, tau_ec, w, color="tab:blue", alpha=0.8,
            label="EC LE direct recovery")
    ax.bar(x + w/2, tau_metv3, w, color="tab:orange", alpha=0.8,
            label="METv3 LE direct recovery")
    for i, (e, m) in enumerate(zip(tau_ec, tau_metv3)):
        if np.isfinite(e):
            ax.text(i - w/2, e + 0.3, f"{e:.2f}", ha="center", fontsize=10)
        if np.isfinite(m):
            ax.text(i + w/2, m + 0.3, f"{m:.2f}", ha="center", fontsize=10)
    ax.set_xticks(x); ax.set_xticklabels(sites)
    ax.set_ylabel("τ [days]  (NaN = no valid fit)")
    ax.set_title("Per-event LE recovery τ : EC vs METv3 satellite")
    ax.grid(axis="y", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ================================================================
# MAIN
# ================================================================

def main():
    args = parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("\n=== 1) loading data ===")
    metv3 = load_metv3_daily(Path(args.metv3_csv))
    tara_ec = load_tara_ec(Path(args.tara_csv))
    oran_ec = load_oran_ec(Path(args.oran_csv))

    rows = []
    for site, ec_df, wcol, thr in [
            ("Oran",     oran_ec, "Rain_mm",  RAIN_THRESHOLD),
            ("Tarazona", tara_ec, "Irrig_mm", IRRIG_THRESHOLD)]:
        print(f"\n=== 2) {site}: build per-event windows ===")
        events_long = build_events_double(ec_df, metv3, site, wcol, thr,
                                              max_window=args.max_window)
        n_events = events_long["event_start"].nunique() if len(events_long) else 0
        n_obs    = len(events_long)
        print(f"  {site}: n_events={n_events}  n_obs={n_obs}")
        if events_long.empty:
            rows.append(dict(site=site, n_events=0, n_obs=0,
                                 tau_EC=np.nan, le0_EC=np.nan,
                                 le_inf_EC=np.nan, r2_EC=np.nan,
                                 tau_METv3=np.nan, le0_METv3=np.nan,
                                 le_inf_METv3=np.nan, r2_METv3=np.nan))
            continue
        fit_ec    = fit_recovery(events_long["days_since_event"].values,
                                   events_long["LE_EC"].values)
        fit_metv3 = fit_recovery(events_long["days_since_event"].values,
                                   events_long["LE_METv3"].values)
        print(f"  EC fit   : "
              f"{(fit_ec or {}).get('tau', 'None')}  "
              f"amp={(fit_ec or {}).get('amplitude', float('nan'))}  "
              f"R²={(fit_ec or {}).get('r2', float('nan')):.2f}"
              if fit_ec else "  EC fit   : None")
        print(f"  METv3 fit: "
              f"{(fit_metv3 or {}).get('tau', 'None')}  "
              f"amp={(fit_metv3 or {}).get('amplitude', float('nan'))}  "
              f"R²={(fit_metv3 or {}).get('r2', float('nan')):.2f}"
              if fit_metv3 else "  METv3 fit: None")
        rows.append(dict(
            site=site, n_events=n_events, n_obs=n_obs,
            tau_EC    =(fit_ec or {}).get("tau",    np.nan),
            le0_EC    =(fit_ec or {}).get("le0",    np.nan),
            le_inf_EC =(fit_ec or {}).get("le_inf", np.nan),
            amp_EC    =(fit_ec or {}).get("amplitude", np.nan),
            r2_EC     =(fit_ec or {}).get("r2",     np.nan),
            tau_METv3   =(fit_metv3 or {}).get("tau",    np.nan),
            le0_METv3   =(fit_metv3 or {}).get("le0",    np.nan),
            le_inf_METv3=(fit_metv3 or {}).get("le_inf", np.nan),
            amp_METv3   =(fit_metv3 or {}).get("amplitude", np.nan),
            r2_METv3    =(fit_metv3 or {}).get("r2",     np.nan),
        ))

        if not args.no_plots:
            out_fig = out_dir / f"fig0{1 if site == 'Oran' else 2}_" \
                          f"{site.lower()}_metv3_vs_ec_recovery.png"
            fig_recovery_panel(events_long, site, fit_ec, fit_metv3, out_fig)
            print(f"  wrote {out_fig.name}")

        # save the long pool
        events_long.to_csv(
            out_dir / f"v5_{site.lower()}_event_pool.csv", index=False)

    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "v5_oran_lst_response.csv", index=False)
    pd.set_option("display.float_format", "{:+.3f}".format)
    print("\n=== SUMMARY ===")
    print(summary[["site", "n_events", "n_obs",
                       "tau_EC", "amp_EC", "r2_EC",
                       "tau_METv3", "amp_METv3", "r2_METv3"]].to_string(
                       index=False))
    print(f"\nwrote {out_dir / 'v5_oran_lst_response.csv'}")

    # ── interpretation ─────────────────────────────────────────────────
    print("\n=== INTERPRETATION ===")
    oran = summary[summary.site == "Oran"].iloc[0]
    tara = summary[summary.site == "Tarazona"].iloc[0]
    if np.isfinite(oran["tau_METv3"]) and oran["amp_METv3"] > 5:
        print("  Oran:  METv3 LE shows a measurable post-rain "
               f"recovery (amp={oran['amp_METv3']:.1f} W/m², "
               f"τ={oran['tau_METv3']:.2f} d). Hypothesis SUPPORTED: rain → "
               "LST cooling → METv3 ET briefly elevated → bias decays "
               "faster than EC LE alone.")
    else:
        print("  Oran:  METv3 LE does not exhibit a clear rain response "
               f"(amp={oran.get('amp_METv3', float('nan')):.1f} W/m²). "
               "Fast bias decay must come from another mechanism — most "
               "likely EC LE itself decaying faster after rain than after "
               "irrigation (different surface state).")
    if not np.isfinite(tara["tau_METv3"]) or tara["amp_METv3"] < 10:
        print("  Tarazona: METv3 LE is approximately flat after irrigation "
               f"(amp={tara.get('amp_METv3', float('nan')):.1f} W/m²) — "
               "confirms 'satellite does not see drip wet-bulb' premise. "
               "Bias decay shape is therefore inherited from EC dynamics.")

    if not args.no_plots:
        fig_contrast(summary, out_dir / "fig03_contrast_decay_speeds.png")
        print(f"  wrote fig03_contrast_decay_speeds.png")

    print(f"\n[DONE] outputs in {out_dir.absolute()}")


if __name__ == "__main__":
    main()
