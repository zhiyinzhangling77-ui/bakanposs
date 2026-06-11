"""
analysis_B_v3_bias_tau.py
=========================

解析B v3: Bias-based τ — Δ(t) = LE_EC(t) - LE_METv3(t) を days_since_event
の関数として fit する。

v2 で判明したこと:
  satellite LE 自体は灌漑シグナルがほぼゼロ (Tarazona amp=1 W/m²) で、
  「satellite LE の τ」では検証不能。これは衛星アルゴリズムが灌漑直後の
  LE パルスを補足できないことの帰結 (TzM MBE = -2.34 mm/d と整合)。

v3 のアプローチ:
  「satellite-EC bias がどう減衰するか」を直接 fit する。
    bias(t=0)  : EC up, satellite flat  → 大きな bias
    bias(t→∞)  : EC が baseline へ      → 小さな bias (or systematic offset)
  これは解析C の tau_fit.py と同じ手法だが、event-based stratification は
  ここで初めて適用される (解析C は days_since_irrig を直接特徴量に)。

入力:
  --metv3-csv : metv3_daily_all.csv
  --tara-csv  : Tarazona EC daily (LE_avg を使う)
  --oran-csv  : Oran     EC daily (Eddy_ET_mm_d × 28.36 を LE とする)

モデル (解析A v27 と同形式):
  bias(t) = a · exp(-t/τ) + c
    a  : 灌漑直後の transient bias    [W/m²]
    τ  : decay time constant          [d]   ← これが目的
    c  : permanent (systematic) bias  [W/m²]

出力先: ./output/analysis_B/v3/
  v3_bias_tau_summary.csv
  v3_bias_perevent_pooled.csv
  fig01_bias_recovery.png
  fig02_bias_tau_comparison.png

Usage:
  python3 analysis_B_v3_bias_tau.py
  python3 analysis_B_v3_bias_tau.py --quick
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

# Reuse v1 helpers
sys.path.insert(0, str(Path(__file__).parent))
from analysis_B_v1_mod16_tau import (
    SITES, EC_TAU_REFERENCE, LE_PER_MM,
    IRRIG_THRESHOLD, RAIN_THRESHOLD,
    _find_col,
)
from analysis_B_v2_metv3_tau import load_metv3_daily

warnings.filterwarnings("ignore")


# ================================================================
# CLI
# ================================================================

TAU_FLOOR = 0.5
TAU_CEIL  = 60.0
TAU_PEG_RATIO = 0.9
R2_MIN_VALID = 0.30
N_EVENTS_MIN_VALID = 3
MIN_PER_BIN = 3
CI_PCT = (2.5, 97.5)


def parse_args():
    p = argparse.ArgumentParser(
        description="解析B v3 — bias-based τ (Δ = EC - METv3) decay",
    )
    p.add_argument("--metv3-csv",
                   default="/home/shion-nagamine/bakanposs/metv3_daily_all.csv")
    p.add_argument("--tara-csv",
                   default="/home/shion-nagamine/Dataset/Eddy data in Spain/"
                           "Daily_Summary_Filtered_forPred_ActEne26.csv")
    p.add_argument("--oran-csv",
                   default="/home/shion-nagamine/Dataset/Eddy data in Spain/"
                           "Oran_EddyDaily_MASTER_2018_2020_correct.csv")
    p.add_argument("--out", default="./output/analysis_B/v3")
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--quick", action="store_true")
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args()


# ================================================================
# EC daily loaders (returning date,LE_EC,water)
# ================================================================

def load_tara_ec(path: Path) -> pd.DataFrame:
    """Tarazona daily EC: returns date, LE_EC_Wm2, Irrig_mm."""
    if not path.exists():
        sys.exit(f"[FATAL] Tarazona EC CSV not found: {path}")
    df = pd.read_csv(path)
    date_col = _find_col(df, ["date", "Date"])
    le_col   = _find_col(df, ["LE_avg", "LE_Wm2", "LE_corr", "LE"])
    irr_col  = _find_col(df, ["Irrig_mm", "Irrigation_mm"])
    if not all([date_col, le_col, irr_col]):
        sys.exit(f"[FATAL] Tarazona EC cols missing: "
                 f"date={date_col}, LE={le_col}, irrig={irr_col}")
    out = pd.DataFrame({
        "date": pd.to_datetime(df[date_col], errors="coerce"),
        "LE_EC_Wm2": pd.to_numeric(df[le_col], errors="coerce"),
        "Irrig_mm":  pd.to_numeric(df[irr_col], errors="coerce").fillna(0),
    }).dropna(subset=["date"]).drop_duplicates("date").sort_values("date")
    print(f"  Tarazona EC : {len(out):,} days  "
          f"LE non-null: {out['LE_EC_Wm2'].notna().sum():,}")
    return out.reset_index(drop=True)


def load_oran_ec(path: Path) -> pd.DataFrame:
    """Oran daily EC: returns date, LE_EC_Wm2 (= ET_mm × 28.36), Rain_mm."""
    if not path.exists():
        sys.exit(f"[FATAL] Oran EC CSV not found: {path}")
    df = pd.read_csv(path)
    date_col = _find_col(df, ["Date", "date"])
    et_col   = _find_col(df, ["Eddy_ET_mm_d", "ET_mm_d", "ET"])
    rain_col = _find_col(df, ["Rain_mm", "P_mm", "Precip_mm"])
    nobs_col = _find_col(df, ["n_obs_le"])
    if not all([date_col, et_col, rain_col]):
        sys.exit(f"[FATAL] Oran EC cols missing: "
                 f"date={date_col}, ET={et_col}, rain={rain_col}")
    out = pd.DataFrame({
        "date":      pd.to_datetime(df[date_col], errors="coerce",
                                      dayfirst=True),
        "LE_EC_Wm2": pd.to_numeric(df[et_col], errors="coerce") * LE_PER_MM,
        "Rain_mm":   pd.to_numeric(df[rain_col], errors="coerce").fillna(0),
    })
    if nobs_col:
        n = pd.to_numeric(df[nobs_col], errors="coerce")
        out.loc[n.fillna(0) < 24, "LE_EC_Wm2"] = np.nan
    out = (out.dropna(subset=["date"])
              .drop_duplicates("date").sort_values("date").reset_index(drop=True))
    print(f"  Oran     EC : {len(out):,} days  "
          f"LE non-null: {out['LE_EC_Wm2'].notna().sum():,}")
    return out


# ================================================================
# Event detection on bias series
# ================================================================

def build_bias_events(merged: pd.DataFrame, water_col: str, threshold: float,
                       min_window: int = 4, max_window: int = 14):
    """merged must have: date, LE_EC_Wm2, LE_METv3_Wm2, <water_col>.

    Returns list of dicts {event_start, data: DataFrame with days_since_event
    and bias = LE_EC - LE_METv3, pulse_mm}.
    """
    merged = merged.sort_values("date").reset_index(drop=True).copy()
    merged["bias"] = merged["LE_EC_Wm2"] - merged["LE_METv3_Wm2"]
    event_dates = merged[merged[water_col].fillna(0) > threshold]["date"].values
    events = []
    for i, ed in enumerate(event_dates):
        ed = pd.Timestamp(ed)
        if i + 1 < len(event_dates):
            next_ed = pd.Timestamp(event_dates[i + 1])
            window_end = min(next_ed - pd.Timedelta(days=1),
                               ed + pd.Timedelta(days=max_window))
        else:
            window_end = ed + pd.Timedelta(days=max_window)
        if (window_end - ed).days < min_window:
            continue
        ev = merged[(merged["date"] >= ed) & (merged["date"] <= window_end)].copy()
        ev["days_since_event"] = (ev["date"] - ed).dt.days
        ev = ev[ev["bias"].notna()]
        if len(ev) < min_window:
            continue
        events.append(dict(
            event_start=ed, data=ev,
            pulse_mm=float(merged.loc[merged["date"] == ed, water_col].iloc[0]),
        ))
    return events


# ================================================================
# Bias decay fit:  Δ(t) = a · exp(-t/τ) + c    (full 3-parameter model)
# ================================================================

def bias_model(t, a, tau, c):
    return a * np.exp(-t / tau) + c


def fit_bias_decay(arr_t, arr_y, min_per_bin: int = MIN_PER_BIN):
    """Median-per-day NLS fit of bias(t) = a*exp(-t/tau) + c."""
    df = pd.DataFrame({"t": arr_t, "y": arr_y}).dropna()
    if df.empty:
        return None
    g = df.groupby("t")["y"].agg(["median", "count"]).reset_index()
    g = g[g["count"] >= min_per_bin]
    if len(g) < 4:    # need >= 4 bins for 3-param model
        return None
    x = g["t"].values.astype(float)
    y = g["median"].values
    y0_guess = y[x == x.min()].mean() if (x == x.min()).any() else y.max()
    y_inf_guess = y[x == x.max()].mean() if (x == x.max()).any() else y.min()
    a_guess = y0_guess - y_inf_guess
    a_lo, a_hi = -500.0, 500.0
    c_lo, c_hi = -500.0, 500.0

    for p0 in [[a_guess, 3.0, y_inf_guess],
                 [a_guess * 1.5, 5.0, y_inf_guess],
                 [a_guess * 0.8, 7.0, y_inf_guess]]:
        try:
            popt, _ = curve_fit(bias_model, x, y, p0=p0,
                                  bounds=([a_lo, TAU_FLOOR, c_lo],
                                          [a_hi, TAU_CEIL, c_hi]),
                                  maxfev=8000)
            yp = bias_model(x, *popt)
            sse = float(np.sum((y - yp) ** 2))
            r2 = (float(1 - sse / np.sum((y - y.mean()) ** 2))
                  if np.var(y) > 0 else np.nan)
            return dict(a=popt[0], tau=popt[1], c=popt[2],
                          r2=r2, sse=sse, n_bins=len(g),
                          amplitude=popt[0])
        except Exception:
            continue
    return None


def bootstrap_bias_taus(arr_t, arr_y, n_boot=2000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(arr_t)
    taus, a_arr, c_arr = [], [], []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        res = fit_bias_decay(arr_t[idx], arr_y[idx])
        if res is None:
            continue
        if res["tau"] >= TAU_CEIL * TAU_PEG_RATIO:
            continue
        if res["tau"] <= TAU_FLOOR * 1.5:
            continue
        if res["r2"] < R2_MIN_VALID:
            continue
        taus.append(res["tau"])
        a_arr.append(res["a"])
        c_arr.append(res["c"])
    return np.array(taus), np.array(a_arr), np.array(c_arr)


def is_bias_fit_valid(res, n_events):
    if n_events < N_EVENTS_MIN_VALID:
        return False, f"n_events={n_events} < {N_EVENTS_MIN_VALID}"
    if res is None:
        return False, "fit returned None"
    if np.isnan(res.get("tau", np.nan)):
        return False, "tau is NaN"
    tau = res["tau"]
    if tau >= TAU_CEIL * TAU_PEG_RATIO:
        return False, f"tau={tau:.1f} pegged at ceiling"
    if tau <= TAU_FLOOR * 1.5:
        return False, f"tau={tau:.2f} pegged at floor"
    if res.get("r2", -np.inf) < R2_MIN_VALID:
        return False, f"R2={res['r2']:.3f} < {R2_MIN_VALID}"
    return True, "OK"


def fit_bias_with_validation(events, n_boot=2000):
    if len(events) < N_EVENTS_MIN_VALID:
        return dict(valid=False, reason=f"n_events={len(events)}<3",
                     tau=np.nan, a=np.nan, c=np.nan, r2=np.nan,
                     ci_lo=np.nan, ci_hi=np.nan, tau_se=np.nan,
                     n_events=len(events), n_data=0, n_boot_valid=0,
                     a_ci_lo=np.nan, a_ci_hi=np.nan,
                     c_ci_lo=np.nan, c_ci_hi=np.nan)
    pooled = pd.concat([ev["data"] for ev in events])
    arr_t = pooled["days_since_event"].values.astype(float)
    arr_y = pooled["bias"].values
    res = fit_bias_decay(arr_t, arr_y)
    is_valid, reason = is_bias_fit_valid(res, len(events))
    out = dict(
        n_events=len(events), n_data=len(pooled),
        tau=res["tau"] if res else np.nan,
        a=res["a"] if res else np.nan,
        c=res["c"] if res else np.nan,
        amplitude=res["a"] if res else np.nan,
        r2=res["r2"] if res else np.nan,
        valid=is_valid, reason=reason,
    )
    if is_valid:
        boot_tau, boot_a, boot_c = bootstrap_bias_taus(arr_t, arr_y,
                                                          n_boot=n_boot)
        if len(boot_tau) >= 100:
            out["ci_lo"], out["ci_hi"] = np.percentile(boot_tau, CI_PCT)
            out["tau_se"] = float(np.std(boot_tau))
            out["a_ci_lo"], out["a_ci_hi"] = np.percentile(boot_a, CI_PCT)
            out["c_ci_lo"], out["c_ci_hi"] = np.percentile(boot_c, CI_PCT)
            out["n_boot_valid"] = len(boot_tau)
        else:
            for k in ("ci_lo", "ci_hi", "tau_se",
                       "a_ci_lo", "a_ci_hi", "c_ci_lo", "c_ci_hi"):
                out[k] = np.nan
            out["n_boot_valid"] = len(boot_tau)
    else:
        for k in ("ci_lo", "ci_hi", "tau_se",
                   "a_ci_lo", "a_ci_hi", "c_ci_lo", "c_ci_hi"):
            out[k] = np.nan
        out["n_boot_valid"] = 0
    return out


# ================================================================
# Per-site driver
# ================================================================

def run_site(ec_df, metv3_df, site, water_col, threshold, n_boot):
    metv3_site = metv3_df[metv3_df["site"] == site][["date", "LE_Wm2"]].rename(
        columns={"LE_Wm2": "LE_METv3_Wm2"})
    merged = ec_df.merge(metv3_site, on="date", how="inner")
    print(f"  {site}: EC × METv3 inner-join = {len(merged):,} days  "
          f"(EC: {len(ec_df)}, METv3: {len(metv3_site)})")
    events = build_bias_events(merged, water_col=water_col, threshold=threshold,
                                  min_window=4, max_window=14)
    if not events:
        return dict(site=site, valid=False, reason="no events",
                      n_events=0), pd.DataFrame()
    pooled = pd.concat([ev["data"] for ev in events])
    pooled["site"] = site
    res = fit_bias_with_validation(events, n_boot=n_boot)
    res.update(site=site, water_col=water_col, threshold=threshold,
                  n_events=len(events))
    return res, pooled


# ================================================================
# Figures
# ================================================================

def fig_bias_recovery(per_event_pooled: dict, fit_results: dict,
                        out_path: Path):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, site in zip(axes, ["Oran", "Tarazona"]):
        pooled = per_event_pooled.get(site, pd.DataFrame())
        res = fit_results.get(site, {})
        if not pooled.empty:
            ax.plot(pooled["days_since_event"], pooled["bias"], "o",
                      ms=3, alpha=0.25, color="gray", label="all event-days")
            g = (pooled.groupby("days_since_event")["bias"]
                          .agg(["median", "count"]).reset_index())
            ax.plot(g["days_since_event"], g["median"], "o", ms=8,
                      color="tab:blue", label="median per t")
        if res.get("valid") and np.isfinite(res.get("tau", np.nan)):
            tau, a, c = res["tau"], res["a"], res["c"]
            xx = np.linspace(0, 14, 200)
            yy = a * np.exp(-xx / tau) + c
            ax.plot(xx, yy, "-", color="tab:red", lw=2,
                      label=(f"fit: τ={tau:.2f} d  "
                             f"a={a:+.1f} W/m²  c={c:+.1f}\n"
                             f"n_ev={res['n_events']}  R²={res['r2']:.2f}"))
        else:
            ax.text(0.05, 0.95,
                      f"[no valid fit]\n{res.get('reason', '?')}",
                      transform=ax.transAxes, va="top")
        ax.axhline(0, color="k", lw=0.8, ls=":")
        ax.set_xlabel("days since water event")
        ax.set_ylabel("bias = LE_EC − LE_METv3  [W/m²]")
        ax.set_title(site)
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right", fontsize=8)
    fig.suptitle("EC − METv3 bias decay after water events", fontsize=12)
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
            color="tab:blue", alpha=0.8,
            label="EC-only (analysis A v27 — LE recovery)")
    ax.bar(x + w/2, sat_v, w, yerr=sat_e, capsize=4,
            color="tab:green", alpha=0.8,
            label="bias decay τ (this work — EC − METv3)")
    for i, (e, s) in enumerate(zip(ec_v, sat_v)):
        if np.isfinite(e):
            ax.text(i - w/2, e + 0.3, f"{e:.2f}", ha="center", fontsize=10)
        if np.isfinite(s):
            ax.text(i + w/2, s + 0.3, f"{s:.2f}", ha="center", fontsize=10)
    ax.set_xticks(x); ax.set_xticklabels(sites)
    ax.set_ylabel("τ [days]")
    ax.set_title("Recovery τ : EC-only vs (EC − METv3) bias decay")
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

    print("\n=== 1) loading METv3 daily ===")
    metv3 = load_metv3_daily(Path(args.metv3_csv))

    print("\n=== 2) loading EC daily ===")
    tara_ec = load_tara_ec(Path(args.tara_csv))
    oran_ec = load_oran_ec(Path(args.oran_csv))

    print("\n=== 3) fitting bias decay τ ===")
    fits = {}
    pooled_by_site = {}
    # Tarazona (drip)
    print("  -- Tarazona (Irrig > 0.5) --")
    res, pooled = run_site(tara_ec, metv3, "Tarazona",
                              "Irrig_mm", IRRIG_THRESHOLD, n_boot)
    fits["Tarazona"] = res
    pooled_by_site["Tarazona"] = pooled
    if res.get("valid"):
        print(f"    τ = {res['tau']:.2f} d  "
               f"CI [{res['ci_lo']:.2f}, {res['ci_hi']:.2f}]  "
               f"a = {res['a']:+.1f} W/m²  c = {res['c']:+.1f}  "
               f"R² = {res['r2']:.3f}  n_ev = {res['n_events']}")
    else:
        print(f"    [invalid] {res.get('reason', '?')}")

    # Oran (rainfed)
    print("  -- Oran (Rain > 3.0) --")
    res, pooled = run_site(oran_ec, metv3, "Oran",
                              "Rain_mm", RAIN_THRESHOLD, n_boot)
    fits["Oran"] = res
    pooled_by_site["Oran"] = pooled
    if res.get("valid"):
        print(f"    τ = {res['tau']:.2f} d  "
               f"CI [{res['ci_lo']:.2f}, {res['ci_hi']:.2f}]  "
               f"a = {res['a']:+.1f} W/m²  c = {res['c']:+.1f}  "
               f"R² = {res['r2']:.3f}  n_ev = {res['n_events']}")
    else:
        print(f"    [invalid] {res.get('reason', '?')}")

    # 4) summary
    print("\n=== 4) EC reference τ (analysis A v27) ===")
    for s, v in EC_TAU_REFERENCE.items():
        print(f"  {s}: τ = {v['tau']:.2f} d (SE {v['se']:.2f})")

    rows = []
    for site in SITES:
        sf = fits[site]
        ef = EC_TAU_REFERENCE[site]
        row = dict(
            site=site,
            tau_ec=ef["tau"], se_ec=ef["se"],
            tau_bias=sf.get("tau", np.nan),
            tau_se_bias=sf.get("tau_se", np.nan),
            tau_ci_lo=sf.get("ci_lo", np.nan),
            tau_ci_hi=sf.get("ci_hi", np.nan),
            a_Wm2=sf.get("a", np.nan),
            a_ci_lo=sf.get("a_ci_lo", np.nan),
            a_ci_hi=sf.get("a_ci_hi", np.nan),
            c_Wm2=sf.get("c", np.nan),
            c_ci_lo=sf.get("c_ci_lo", np.nan),
            c_ci_hi=sf.get("c_ci_hi", np.nan),
            n_events=sf.get("n_events", 0),
            r2=sf.get("r2", np.nan),
            valid=sf.get("valid", False),
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
    summary.to_csv(out_dir / "v3_bias_tau_summary.csv", index=False)

    # also save per-event pooled bias (for downstream plotting)
    if any(len(p) for p in pooled_by_site.values()):
        all_pooled = pd.concat(
            [p for p in pooled_by_site.values() if len(p)], ignore_index=True)
        all_pooled.to_csv(out_dir / "v3_bias_perevent_pooled.csv", index=False)

    pd.set_option("display.float_format", "{:+.3f}".format)
    print("\n=== SUMMARY ===")
    cols_show = ["site", "tau_ec", "se_ec",
                   "tau_bias", "tau_se_bias", "tau_ci_lo", "tau_ci_hi",
                   "a_Wm2", "c_Wm2", "n_events", "r2",
                   "MDE", "diff", "significant_diff", "reason"]
    print(summary[cols_show].to_string(index=False))
    print(f"\nwrote {out_dir / 'v3_bias_tau_summary.csv'}")

    if args.no_plots:
        print("[skipped plots --no-plots]")
        return

    print("\n=== 5) generating figures ===")
    fig_bias_recovery(pooled_by_site, fits,
                          out_dir / "fig01_bias_recovery.png")
    print("  wrote fig01_bias_recovery.png")
    fig_tau_compare(fits, EC_TAU_REFERENCE,
                       out_dir / "fig02_bias_tau_comparison.png")
    print("  wrote fig02_bias_tau_comparison.png")

    print(f"\n[DONE] outputs in {out_dir.absolute()}")


if __name__ == "__main__":
    main()
