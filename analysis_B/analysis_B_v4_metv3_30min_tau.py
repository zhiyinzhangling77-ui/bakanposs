"""
analysis_B_v4_metv3_30min_tau.py
================================

解析B v4: METv3 を **30-min cadence** (= 0.0208 d 解像度) でそのまま使い、
event 後の LE 回復曲線を sub-daily で fit する。

v2 で見えたこと:
  METv3 daily aggregation では amplitude が壊滅的に小さい (Tarazona 1 W/m²)。
  daily 平均化で 30-min のピークが消える可能性 → 高解像度で再挑戦。

v4 のロジック:
  1) 各 event 日について、event_date から +max_window 日分の NetCDF を読む
     (1 day x 48 半時間値 = ~672 ファイル / event window / site)
  2) 30-min 瞬時 ET [mm/h] → LE [W/m²] = ET_mm_per_h * LE_PER_MM
  3) days_since_event を連続 (0.0, 0.021, 0.042, ...) で持つ
  4) median-per-bin (例えば 0.25-day bin) で fit
  5) 解析A v27 と同じ exp_model_2p で τ 推定

入力:
  --metv3-base : /home/shion-nagamine/Dataset/METv3/
                 配下に YYYY/MM/DDMM/YYYYMMDD_HHMM.nc が並ぶ
  --tara-csv   : Tarazona Irrig 用
  --oran-csv   : Oran     Rain  用
  --bin-hours  : days_since_event のビン幅 (デフォルト 6h = 0.25 d)
  --max-events : サブセット高速化用 (デフォルト = 全 event)

出力先: ./output_analysis_B_v4/
  v4_metv3_30min_tau_summary.csv
  v4_metv3_30min_event_pool.csv     (raw 30-min event-pool data)
  fig01_30min_recovery.png

Usage:
  python3 analysis_B_v4_metv3_30min_tau.py --quick
  python3 analysis_B_v4_metv3_30min_tau.py --max-events 30
  python3 analysis_B_v4_metv3_30min_tau.py             # 全 event, n_boot=2000

NOTE: フル実行で event x 14 日 x 48 ファイル x 2 site = 数万 NetCDF read。
      初回 30-60 分かかる可能性あり (キャッシュなし)。--max-events で短縮可。
"""

from __future__ import annotations
import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Reuse helpers from v1 / v2
sys.path.insert(0, str(Path(__file__).parent))
from analysis_B_v1_mod16_tau import (
    SITES, EC_TAU_REFERENCE, LE_PER_MM,
    IRRIG_THRESHOLD, RAIN_THRESHOLD,
    load_tara_irrig, load_oran_rain,
)

# Sub-daily unit conversion:
#   ET [mm/h]  ×  λρ_w/3600  =  LE [W/m²]
#   λ ≈ 2.45e6 J/kg, ρ_w = 1000 kg/m³  →  680.6 W/m² per (mm/h)
# (= 24 × daily LE_PER_MM, because 1 mm/h sustained for 24h = 24 mm/day)
LE_PER_MM_PER_HOUR = 24.0 * LE_PER_MM   # ≈ 680.6 W/m² per (mm/h)

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
DT_HOURS = 0.5


def parse_args():
    p = argparse.ArgumentParser(
        description="解析B v4 — METv3 30-min τ (sub-daily resolution)",
    )
    p.add_argument("--metv3-base",
                   default="/home/shion-nagamine/Dataset/METv3",
                   help="METv3 NetCDF base directory")
    p.add_argument("--tara-csv",
                   default="/home/shion-nagamine/Dataset/Eddy data in Spain/"
                           "Daily_Summary_Filtered_forPred_ActEne26.csv")
    p.add_argument("--oran-csv",
                   default="/home/shion-nagamine/Dataset/Eddy data in Spain/"
                           "Oran_EddyDaily_MASTER_2018_2020_correct.csv")
    p.add_argument("--out", default="./output_analysis_B_v4")
    p.add_argument("--max-window", type=int, default=10,
                   help="post-event window in days (default 10; 8-14 typical)")
    p.add_argument("--bin-hours", type=float, default=6.0,
                   help="days_since_event bin width in hours (only used "
                        "when --aggregate raw)")
    p.add_argument("--aggregate", choices=["noon", "daymax", "raw"],
                   default="noon",
                   help=("noon: daily mean over 10:00-14:00 UTC -> one obs "
                         "per integer day per event (cancels diurnal cycle, "
                         "recommended). daymax: daily max LE. raw: pool all "
                         "30-min obs (dominated by diurnal cycle, low R^2)."))
    p.add_argument("--max-events", type=int, default=0,
                   help="cap events per site (0 = no cap)")
    p.add_argument("--daytime-only", action="store_true",
                   help="restrict to 09:00-15:00 UTC (peak LE; reduces nights "
                        "that bias the fit toward 0)")
    p.add_argument("--n-boot", type=int, default=500,
                   help="bootstrap iters (default 500 — slower than v2 because "
                        "n_obs per event is ~480 vs ~10)")
    p.add_argument("--quick", action="store_true",
                   help="n_boot=100 + max_events=10 (sanity test)")
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args()


# ================================================================
# NetCDF read at single pixel for time window
# ================================================================

def _file_path(base: Path, ts: pd.Timestamp) -> Path:
    """YYYY/MM/MMDD/YYYYMMDD_HHMM.nc"""
    mmdd = f"{ts.month:02d}{ts.day:02d}"
    name = f"{ts.year:04d}{ts.month:02d}{ts.day:02d}_{ts.hour:02d}{ts.minute:02d}.nc"
    return base / f"{ts.year}" / f"{ts.month:02d}" / mmdd / name


def _resolve_indices(sample_file: Path) -> dict:
    """Find nearest (i_lat, i_lon) for each site."""
    import xarray as xr
    with xr.open_dataset(sample_file, engine="h5netcdf") as ds:
        lats = ds["lat"].values
        lons = ds["lon"].values
    out = {}
    for name, coords in SITES.items():
        i_lat = int(np.argmin(np.abs(lats - coords["lat"])))
        i_lon = int(np.argmin(np.abs(lons - coords["lon"])))
        out[name] = (i_lat, i_lon)
    return out


def read_window(base: Path, t_start: pd.Timestamp, t_end: pd.Timestamp,
                 site: str, indices: dict) -> pd.DataFrame:
    """Read all 30-min METv3 files between t_start and t_end (inclusive),
    return DataFrame with ts, ET_mm_per_h."""
    import xarray as xr
    i_lat, i_lon = indices[site]
    rows = []
    ts = t_start.floor("30min")
    while ts <= t_end:
        p = _file_path(base, ts)
        if p.exists():
            try:
                with xr.open_dataset(p, engine="h5netcdf") as ds:
                    v = float(ds["ET"].values[0, i_lat, i_lon])
                rows.append({"ts": ts, "ET_mm_per_h": v})
            except Exception:
                pass
        ts += pd.Timedelta(minutes=30)
    return pd.DataFrame(rows)


# ================================================================
# Event detection (daily water table → event dates)
# ================================================================

def detect_event_dates(water_df: pd.DataFrame, water_col: str,
                          threshold: float, min_gap_days: int = 4):
    """Return event dates spaced at least min_gap_days apart."""
    df = water_df[water_df[water_col].fillna(0) > threshold].sort_values("date")
    keep = []
    last = None
    for d in df["date"]:
        d = pd.Timestamp(d)
        if last is None or (d - last).days >= min_gap_days:
            keep.append(d)
            last = d
    return keep


# ================================================================
# Sub-daily τ fit (mirrors v1 but on continuous days_since_event)
# ================================================================

def exp_model_2p(d, le0, tau, le_inf_fixed):
    return le_inf_fixed + (le0 - le_inf_fixed) * np.exp(-d / tau)


def fit_subdaily_tau(arr_d, arr_y, le_inf_fixed, bin_d):
    """Bin arr_d at bin_d (days) resolution, take median, fit exponential."""
    df = pd.DataFrame({"d": arr_d, "y": arr_y}).dropna()
    if df.empty:
        return None
    df["bin"] = (df["d"] / bin_d).round() * bin_d
    g = df.groupby("bin")["y"].agg(["median", "count"]).reset_index()
    g = g[g["count"] >= MIN_PER_BIN]
    if len(g) < 4:
        return None
    x = g["bin"].values.astype(float)
    y = g["median"].values
    le0_floor = max(le_inf_fixed + 1, 10.0)
    le0_ceil = max(y.max() * 1.5, 350.0)

    def model(d, le0, tau):
        return exp_model_2p(d, le0, tau, le_inf_fixed)

    for p0 in [[y.max(), 5.0],
                 [y.max() * 1.1, 3.0],
                 [(y.max() + le_inf_fixed) / 2, 7.0]]:
        try:
            popt, _ = curve_fit(model, x, y, p0=p0,
                                  bounds=([le0_floor, TAU_FLOOR],
                                          [le0_ceil, TAU_CEIL]),
                                  maxfev=8000)
            yp = model(x, *popt)
            sse = float(np.sum((y - yp) ** 2))
            r2 = (float(1 - sse / np.sum((y - y.mean()) ** 2))
                  if np.var(y) > 0 else np.nan)
            return dict(le0=popt[0], tau=popt[1], r2=r2, sse=sse,
                          n_bins=len(g), amplitude=popt[0] - le_inf_fixed)
        except Exception:
            continue
    return None


def bootstrap_subdaily(arr_d, arr_y, le_inf, bin_d, n_boot=500, seed=42):
    rng = np.random.default_rng(seed)
    n = len(arr_d)
    taus = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        res = fit_subdaily_tau(arr_d[idx], arr_y[idx], le_inf, bin_d)
        if res is None:
            continue
        if res["tau"] >= TAU_CEIL * TAU_PEG_RATIO:
            continue
        if res["tau"] <= TAU_FLOOR * 1.5:
            continue
        if res["r2"] < R2_MIN_VALID:
            continue
        taus.append(res["tau"])
    return np.array(taus)


def estimate_le_inf(pooled: pd.DataFrame, tail_d: float = 7.0) -> float:
    tail = pooled[pooled["days_since_event"] >= tail_d]["LE"]
    if tail.dropna().shape[0] < 5:
        return float(np.percentile(pooled["LE"].dropna(), 25))
    return float(np.median(tail.dropna()))


def aggregate_per_day(pool: pd.DataFrame, mode: str) -> pd.DataFrame:
    """Collapse 30-min observations to one obs per (event, integer day).

    mode='noon'   : mean LE over 10:00-14:00 UTC  (peak transpiration window)
    mode='daymax' : daily max LE within the integer day
    mode='raw'    : return pool unchanged (caller does the binning)
    """
    if mode == "raw":
        return pool
    df = pool.copy()
    df["day"] = np.floor(df["days_since_event"]).astype(int)
    if mode == "noon":
        df = df[df["ts"].dt.hour.between(10, 13)]
        agg = (df.groupby(["event_start", "day"])
                  ["LE"].mean().reset_index())
    elif mode == "daymax":
        agg = (df.groupby(["event_start", "day"])
                  ["LE"].max().reset_index())
    else:
        return pool
    agg = agg.rename(columns={"day": "days_since_event"})
    return agg


# ================================================================
# MAIN per-site driver
# ================================================================

def collect_events(base: Path, water_df: pd.DataFrame, water_col: str,
                     threshold: float, site: str, indices: dict,
                     max_window: int, max_events: int,
                     daytime_only: bool) -> pd.DataFrame:
    event_dates = detect_event_dates(water_df, water_col, threshold)
    if max_events > 0:
        event_dates = event_dates[:max_events]
    print(f"    {len(event_dates)} events detected (window={max_window} d)")
    pieces = []
    t0 = time.time()
    for k, ed in enumerate(event_dates):
        t_end = ed + pd.Timedelta(days=max_window) - pd.Timedelta(minutes=30)
        win = read_window(base, ed, t_end, site, indices)
        if win.empty:
            continue
        win["days_since_event"] = (win["ts"] - ed).dt.total_seconds() / 86400.0
        win["LE"] = win["ET_mm_per_h"].clip(lower=0) * LE_PER_MM_PER_HOUR
        if daytime_only:
            win = win[win["ts"].dt.hour.between(9, 14)]
        win["event_start"] = ed
        win["pulse_mm"] = float(water_df.loc[
            water_df["date"] == ed, water_col].iloc[0])
        pieces.append(win)
        if (k + 1) % 10 == 0:
            dt = time.time() - t0
            eta = dt / (k + 1) * (len(event_dates) - k - 1)
            print(f"      [{k + 1}/{len(event_dates)}] elapsed {dt:.0f}s  "
                   f"ETA {eta:.0f}s")
    if not pieces:
        return pd.DataFrame()
    return pd.concat(pieces, ignore_index=True)


# ================================================================
# Figures
# ================================================================

def fig_recovery(pool_by_site: dict, fits: dict, bin_d: float,
                   out_path: Path):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, site in zip(axes, ["Oran", "Tarazona"]):
        pool = pool_by_site.get(site, pd.DataFrame())
        res = fits.get(site, {})
        if not pool.empty:
            ax.plot(pool["days_since_event"], pool["LE"], ".", ms=2,
                      alpha=0.15, color="gray")
            tmp = pool.copy()
            tmp["bin"] = (tmp["days_since_event"] / bin_d).round() * bin_d
            g = (tmp.groupby("bin")["LE"]
                          .agg(["median", "count"]).reset_index())
            g = g[g["count"] >= MIN_PER_BIN]
            ax.plot(g["bin"], g["median"], "o", ms=6,
                      color="tab:purple",
                      label=f"median per {bin_d*24:.0f}h bin")
        if res.get("valid") and np.isfinite(res.get("tau", np.nan)):
            tau, le0, le_inf = res["tau"], res["le0"], res["le_inf"]
            xx = np.linspace(0, 14, 400)
            yy = le_inf + (le0 - le_inf) * np.exp(-xx / tau)
            ax.plot(xx, yy, "-", color="tab:red", lw=2,
                      label=(f"fit: τ={tau:.2f} d, LE0={le0:.1f}, "
                             f"LE∞={le_inf:.1f}\n"
                             f"amp={res['amplitude']:.1f} W/m², "
                             f"R²={res['r2']:.2f}"))
        else:
            ax.text(0.05, 0.95,
                      f"[no valid fit]\n{res.get('reason', '?')}",
                      transform=ax.transAxes, va="top")
        ax.set_xlabel("days since event")
        ax.set_ylabel("LE [W/m²] (METv3 30-min)")
        ax.set_title(site)
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right", fontsize=8)
    fig.suptitle("METv3 30-min LE recovery after water events", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


# ================================================================
# MAIN
# ================================================================

def main():
    args = parse_args()
    if args.quick:
        args.n_boot = 100
        if args.max_events == 0:
            args.max_events = 10
    bin_d = args.bin_hours / 24.0
    base = Path(args.metv3_base)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not base.exists():
        sys.exit(f"[FATAL] METv3 base not found: {base}")

    # discover nearest-pixel indices from any existing file
    print("\n=== 1) resolving nearest METv3 pixel for each site ===")
    sample = next((p for p in base.rglob("*.nc")), None)
    if sample is None:
        sys.exit(f"[FATAL] no NetCDFs under {base}")
    indices = _resolve_indices(sample)
    for s, (ila, ilo) in indices.items():
        print(f"  {s}: i_lat={ila}, i_lon={ilo}")

    print("\n=== 2) loading water-input series ===")
    tara_water = load_tara_irrig(Path(args.tara_csv))
    oran_water = load_oran_rain(Path(args.oran_csv))
    print(f"  Tarazona Irrig events (>{IRRIG_THRESHOLD}): "
          f"{(tara_water.get('Irrig_mm', pd.Series([], dtype=float)) > IRRIG_THRESHOLD).sum()}")
    print(f"  Oran     Rain  events (>{RAIN_THRESHOLD}):  "
          f"{(oran_water.get('Rain_mm', pd.Series([], dtype=float)) > RAIN_THRESHOLD).sum()}")

    print("\n=== 3) collecting 30-min event windows ===")
    pool_by_site = {}
    print("  -- Tarazona --")
    pool_by_site["Tarazona"] = collect_events(
        base, tara_water, "Irrig_mm", IRRIG_THRESHOLD, "Tarazona",
        indices, args.max_window, args.max_events, args.daytime_only)
    print("  -- Oran --")
    pool_by_site["Oran"] = collect_events(
        base, oran_water, "Rain_mm", RAIN_THRESHOLD, "Oran",
        indices, args.max_window, args.max_events, args.daytime_only)

    print(f"\n=== 4) aggregating per (event,day) with mode={args.aggregate!r} ===")
    fits = {}
    for site in SITES:
        pool = pool_by_site[site]
        if pool.empty:
            fits[site] = dict(valid=False, reason="empty pool", n_events=0)
            print(f"  -- {site} -- [invalid] empty pool")
            continue
        agg = aggregate_per_day(pool, args.aggregate)
        n_events = agg["event_start"].nunique() if "event_start" in agg else \
                       pool["event_start"].nunique()
        n_obs = len(agg)
        # bin width: integer-day for noon/daymax, sub-daily for raw
        eff_bin_d = 1.0 if args.aggregate in ("noon", "daymax") else bin_d
        le_inf = estimate_le_inf(agg)
        print(f"  -- {site} -- n_events={n_events}, n_obs={n_obs}, "
               f"le_inf={le_inf:.1f} W/m²  "
               f"(LE range: {agg['LE'].min():.1f} → {agg['LE'].max():.1f})")
        res = fit_subdaily_tau(
            agg["days_since_event"].values.astype(float),
            agg["LE"].values, le_inf, eff_bin_d)
        if res is None:
            fits[site] = dict(valid=False,
                                reason="fit returned None",
                                n_events=n_events, n_obs=n_obs,
                                le_inf=le_inf,
                                tau=np.nan, le0=np.nan,
                                amplitude=np.nan, r2=np.nan,
                                ci_lo=np.nan, ci_hi=np.nan, tau_se=np.nan)
            print(f"  -- {site} -- [invalid] fit None  "
                   f"(n_obs={n_obs}, n_events={n_events})")
            continue
        # validation
        valid = (res["r2"] >= R2_MIN_VALID
                  and TAU_FLOOR * 1.5 < res["tau"] < TAU_CEIL * TAU_PEG_RATIO
                  and n_events >= N_EVENTS_MIN_VALID)
        out = dict(site=site, valid=valid,
                       tau=res["tau"], le0=res["le0"], le_inf=le_inf,
                       amplitude=res["amplitude"], r2=res["r2"],
                       n_events=n_events, n_obs=n_obs, n_bins=res["n_bins"])
        if valid:
            boot = bootstrap_subdaily(
                agg["days_since_event"].values.astype(float),
                agg["LE"].values, le_inf, eff_bin_d, n_boot=args.n_boot)
            if len(boot) >= 50:
                out["ci_lo"], out["ci_hi"] = np.percentile(boot, CI_PCT)
                out["tau_se"] = float(np.std(boot))
                out["n_boot_valid"] = len(boot)
            else:
                out["ci_lo"] = out["ci_hi"] = out["tau_se"] = np.nan
                out["n_boot_valid"] = len(boot)
            out["reason"] = "OK"
        else:
            out["ci_lo"] = out["ci_hi"] = out["tau_se"] = np.nan
            out["n_boot_valid"] = 0
            out["reason"] = (f"R2={res['r2']:.3f} <0.3 "
                                f"or tau={res['tau']:.2f} pegged")
        fits[site] = out
        if valid:
            print(f"  -- {site} --  τ = {out['tau']:.2f} d  "
                   f"CI [{out['ci_lo']:.2f}, {out['ci_hi']:.2f}]  "
                   f"amp = {out['amplitude']:.1f} W/m²  "
                   f"R² = {out['r2']:.3f}  n_ev={n_events}, n_obs={n_obs}")
        else:
            print(f"  -- {site} -- [invalid] {out['reason']}  "
                   f"(τ={out['tau']:.2f}, R²={out['r2']:.3f}, "
                   f"n_obs={n_obs}, n_events={n_events})")

    # 5) summary CSV
    rows = []
    for site in SITES:
        sf = fits[site]
        ef = EC_TAU_REFERENCE[site]
        row = dict(
            site=site,
            tau_ec=ef["tau"], se_ec=ef["se"],
            tau_sat=sf.get("tau", np.nan),
            tau_se_sat=sf.get("tau_se", np.nan),
            ci_lo_sat=sf.get("ci_lo", np.nan),
            ci_hi_sat=sf.get("ci_hi", np.nan),
            n_events=sf.get("n_events", 0),
            n_obs=sf.get("n_obs", 0),
            amp_sat=sf.get("amplitude", np.nan),
            le_inf=sf.get("le_inf", np.nan),
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
            row["MDE"] = row["diff"] = np.nan
            row["significant_diff"] = None
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "v4_metv3_30min_tau_summary.csv", index=False)

    if any(len(p) for p in pool_by_site.values()):
        all_pool = pd.concat([p.assign(site=s) for s, p in pool_by_site.items()
                                if len(p)], ignore_index=True)
        all_pool.to_csv(out_dir / "v4_metv3_30min_event_pool.csv",
                            index=False)

    pd.set_option("display.float_format", "{:+.3f}".format)
    print("\n=== SUMMARY ===")
    cols_show = ["site", "tau_ec", "se_ec", "tau_sat", "tau_se_sat",
                   "n_events", "n_obs", "amp_sat", "r2",
                   "MDE", "diff", "significant_diff", "reason"]
    print(summary[cols_show].to_string(index=False))
    print(f"\nwrote {out_dir / 'v4_metv3_30min_tau_summary.csv'}")

    if args.no_plots:
        print("[skipped plots --no-plots]")
        return

    print("\n=== 5) generating figures ===")
    fig_recovery(pool_by_site, fits, bin_d,
                   out_dir / "fig01_30min_recovery.png")
    print("  wrote fig01_30min_recovery.png")
    print(f"\n[DONE] outputs in {out_dir.absolute()}")


if __name__ == "__main__":
    main()
