"""
analysis_B_v1_mod16_tau.py
==========================

解析B v1: 衛星 ET (MOD16A2) で 解析A の τ ≈ 3-4 d を独立検証する。

入力:
  /mnt/hdd/Dataset/MOD16A2_ET_Oran_TzM_2018_2024.tif
    GeoTIFF, 各バンドが 8-day composite (MOD16A2 標準)
    バンドの time-stamp は GDAL metadata から復元するか、
    --start-date と --cadence 8 で構築

  解析A の event 定義に必要な水入力データ:
    Tarazona: /home/shion-nagamine/Dataset/Eddy data in Spain/
              Daily_Summary_Filtered_forPred_ActEne26.csv  (Irrig_mm)
    Oran:     /home/shion-nagamine/Dataset/Eddy data in Spain/
              Oran_EddyDaily_MASTER_2018_2020_correct.csv  (Rain_mm)

  EC 比較用 (解析A の τ と並べる):
    /home/shion-nagamine/bakanposs/analysis_A/daily_classified_v4.parquet

サイト座標:
  Oran     : lat 38.82,   lon -1.86
  Tarazona : lat 39.266,  lon -1.9397

処理:
  1) GeoTIFF を rasterio で開く
  2) Oran / Tarazona の nearest pixel から daily 時系列を抽出
     (MOD16A2 は 8-day、各 timestep の値を期間内 daily に forward-fill)
  3) 単位変換: MOD16A2 ET は kg/m²/8day (scale 0.1) → mm/day
     LE [W/m²] = ET [mm/day] × 28.36 (latent heat conversion)
  4) Tarazona/Oran それぞれで:
     - 解析A と同じ event 検出 (Irrig>0.5 または Rain>3.0)
     - 解析A v27 の fit_pooled_tau / bootstrap_taus を import して使う
  5) 出力: τ_satellite vs τ_EC, amplitude 比較, 図 3 枚

出力先: ./output_analysis_B_v1/
  v1_satellite_tau_summary.csv
  v1_satellite_vs_ec_tau.csv
  fig01_pixel_extraction.png       # サイトと nearest pixel の関係
  fig02_satellite_timeseries.png   # 抽出した 8-day ET 時系列
  fig03_satellite_recovery.png     # event 後の τ-fit カーブ
  fig04_tau_comparison.png         # 衛星 τ vs EC τ pairwise

Usage:
  python analysis_B_v1_mod16_tau.py
  python analysis_B_v1_mod16_tau.py --start-date 2018-01-01 --cadence 8
  python analysis_B_v1_mod16_tau.py --quick     # 確認用 (n_boot=300)
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Reuse analysis A v27 fit logic
sys.path.insert(0, str(Path(__file__).parent))
from analysis_A_v27 import (
    fit_pooled_tau, bootstrap_taus, is_fit_valid,
    detect_events, estimate_le_inf, fit_with_validation,
    TAU_FLOOR, TAU_CEIL, R2_MIN_VALID,
)

# ================================================================
# CLI + paths
# ================================================================

SITES = {
    "Oran":     {"lat": 38.82,   "lon": -1.86},
    "Tarazona": {"lat": 39.266,  "lon": -1.9397},
}

# MOD16A2 scale handling.
#   - Raw HDF / SDS: stored x 0.1 = kg/m^2/8day, fill values 32760+
#   - GEE-exported GeoTIFF: scale ALREADY applied (values in mm/8day,
#     fills NaN-coded), so use 1.0
# Set --scale 0.1 on the CLI if the file is the raw unscaled product.
DEFAULT_MOD16_SCALE = 1.0
# 1 mm/day -> 28.36 W/m^2 (lambda * rho_w / 86400 s/day ~ 28.36)
LE_PER_MM = 28.36
# MOD16 fill values (only relevant if raw unscaled)
MOD16_FILL = {32767, 32766, 32765, 32764, 32763, 32762, 32761}

# Event thresholds (consistent with analysis A v27)
IRRIG_THRESHOLD = 0.5    # mm/day
RAIN_THRESHOLD  = 3.0    # mm/day


def parse_args():
    p = argparse.ArgumentParser(description="解析B v1 — MOD16A2 で τ 再現")
    p.add_argument("--tif",
                   default="/mnt/hdd/Dataset/MOD16A2_ET_Oran_TzM_2018_2024.tif")
    p.add_argument("--start-date", default="2018-01-01",
                   help="first band date (ISO format)")
    p.add_argument("--cadence", type=int, default=8,
                   help="days between bands (default 8 = MOD16A2 standard)")
    p.add_argument("--tara-csv",
                   default="/home/shion-nagamine/Dataset/Eddy data in Spain/"
                           "Daily_Summary_Filtered_forPred_ActEne26.csv")
    p.add_argument("--oran-csv",
                   default="/home/shion-nagamine/Dataset/Eddy data in Spain/"
                           "Oran_EddyDaily_MASTER_2018_2020_correct.csv")
    p.add_argument("--ec-parquet",
                   default="/home/shion-nagamine/bakanposs/analysis_A/"
                           "daily_classified_v4.parquet")
    p.add_argument("--out", default="./output_analysis_B_v1")
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--scale", type=float, default=DEFAULT_MOD16_SCALE,
                   help="MOD16 multiplier. Use 1.0 for GEE-exported TIFFs "
                        "(scale already applied) or 0.1 for raw HDF.")
    p.add_argument("--quick", action="store_true",
                   help="n_boot=300, plots only essential")
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args()


# ================================================================
# GeoTIFF reader
# ================================================================

def read_mod16_tif(tif_path: Path, start_date: str, cadence_days: int,
                    scale: float = DEFAULT_MOD16_SCALE) -> pd.DataFrame:
    """Open the MOD16A2 multi-band GeoTIFF and extract the time series
    at each site (nearest pixel).

    Returns long-format DataFrame:
      date, site, ET_mm_per_day, LE_Wm2
    """
    try:
        import rasterio
        from rasterio.transform import rowcol
    except ImportError as e:
        sys.exit(f"rasterio not installed: {e}\n  pip install rasterio")

    if not tif_path.exists():
        sys.exit(f"[FATAL] TIFF not found: {tif_path}")

    rows = []
    with rasterio.open(tif_path) as src:
        print(f"opened {tif_path.name}")
        print(f"  bands : {src.count}")
        print(f"  size  : {src.width} x {src.height}")
        print(f"  CRS   : {src.crs}")
        print(f"  bounds: {src.bounds}")

        # Try to read band dates from GDAL metadata first
        band_dates = []
        for b in range(1, src.count + 1):
            tags = src.tags(b)
            date_str = (tags.get("system:time_start") or tags.get("date")
                         or tags.get("DATE") or tags.get("BAND_DATE"))
            if date_str:
                # GEE epoch-ms support
                try:
                    if str(date_str).isdigit() and len(str(date_str)) >= 12:
                        ts = pd.Timestamp(int(date_str), unit="ms")
                    else:
                        ts = pd.Timestamp(date_str)
                    band_dates.append(ts)
                    continue
                except Exception:
                    pass
            band_dates.append(None)

        if any(d is None for d in band_dates):
            # fall back to evenly spaced from start_date
            n_missing = sum(d is None for d in band_dates)
            print(f"  [info] {n_missing}/{src.count} bands missing date metadata; "
                  f"reconstructing from --start-date={start_date} "
                  f"--cadence={cadence_days}")
            base = pd.Timestamp(start_date)
            band_dates = [base + pd.Timedelta(days=cadence_days * i)
                          for i in range(src.count)]

        print(f"  date range: {min(band_dates).date()} -> {max(band_dates).date()}")

        # transformer to convert lat/lon -> pixel row/col
        # if CRS is geographic (EPSG:4326), can use lat/lon directly
        # otherwise need to transform
        from rasterio.warp import transform
        for site, coords in SITES.items():
            lat, lon = coords["lat"], coords["lon"]
            if src.crs and src.crs.to_string() != "EPSG:4326":
                xs, ys = transform("EPSG:4326", src.crs, [lon], [lat])
                x, y = xs[0], ys[0]
            else:
                x, y = lon, lat
            try:
                row, col = src.index(x, y)
            except Exception as e:
                print(f"  [warn] {site}: site outside raster bounds ({e})")
                continue

            print(f"  {site}: lat={lat}, lon={lon} -> pixel (row={row}, col={col})")

            # bounds check
            if row < 0 or row >= src.height or col < 0 or col >= src.width:
                print(f"    [skip] pixel out of bounds")
                continue

            # read the time-series at that one pixel
            window = ((row, row + 1), (col, col + 1))
            vals = src.read(window=window).astype(float).flatten()
            assert len(vals) == src.count, f"expected {src.count} got {len(vals)}"

            for d, v in zip(band_dates, vals):
                # filter MOD16 fill values (raw HDF) and NaN (GEE exports)
                if not np.isfinite(v):
                    et_8day = np.nan
                elif scale == 0.1 and int(v) in MOD16_FILL:
                    et_8day = np.nan
                else:
                    et_8day = v * scale         # kg/m^2/8day
                et_per_day = et_8day / 8.0      # mm/day
                rows.append({
                    "date": d, "site": site,
                    "ET_mm_per_day": et_per_day,
                    "LE_Wm2": et_per_day * LE_PER_MM,
                    "raw_value": v,
                })

    df = pd.DataFrame(rows)
    return df


def expand_to_daily(eight_day_df: pd.DataFrame) -> pd.DataFrame:
    """Forward-fill the 8-day MOD16 composite to a daily time series.

    Each composite represents the 8-day period STARTING on its band date,
    so we propagate the value forward up to (cadence-1) days.
    """
    pieces = []
    for site, sub in eight_day_df.groupby("site"):
        sub = sub.sort_values("date").reset_index(drop=True)
        if sub.empty:
            continue
        idx = pd.date_range(sub["date"].min(),
                             sub["date"].max() + pd.Timedelta(days=7),
                             freq="D")
        s = sub.set_index("date").reindex(idx)
        # forward-fill ET / LE for up to 7 days within the composite window
        for col in ["ET_mm_per_day", "LE_Wm2"]:
            s[col] = s[col].ffill(limit=7)
        s["site"] = site
        s = s.reset_index().rename(columns={"index": "date"})
        pieces.append(s[["date", "site", "ET_mm_per_day", "LE_Wm2"]])
    return pd.concat(pieces, ignore_index=True)


# ================================================================
# water-input loaders (event detection inputs)
# ================================================================

def load_tara_irrig(tara_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(tara_csv)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    if "Irrig_mm" not in df.columns:
        sys.exit(f"Irrig_mm missing in {tara_csv}")
    return (df[["date", "Irrig_mm"]].dropna(subset=["date"])
            .drop_duplicates("date").sort_values("date").reset_index(drop=True))


def load_oran_rain(oran_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(oran_csv)
    date_col = next((c for c in ["Date", "date", "TIMESTAMP"] if c in df.columns), None)
    rain_col = next((c for c in ["Rain_mm", "P_mm", "Precip_mm"]
                     if c in df.columns), None)
    if not all([date_col, rain_col]):
        sys.exit(f"date/rain cols missing in {oran_csv}")
    out = pd.DataFrame({
        "date": pd.to_datetime(df[date_col], errors="coerce"),
        "Rain_mm": pd.to_numeric(df[rain_col], errors="coerce"),
    }).dropna(subset=["date"]).drop_duplicates("date").sort_values("date")
    return out.reset_index(drop=True)


# ================================================================
# τ-fit at one site (uses analysis_A_v27 functions)
# ================================================================

def fit_satellite_tau(sat_daily: pd.DataFrame, water_df: pd.DataFrame,
                      site: str, water_col: str, threshold: float,
                      n_boot: int = 2000) -> dict:
    """Run the v27 event detection + fit on satellite-derived LE."""
    if sat_daily.empty:
        return dict(site=site, valid=False, reason="empty satellite series")

    df = sat_daily[sat_daily["site"] == site][["date", "LE_Wm2"]].copy()
    df = df.rename(columns={"LE_Wm2": "LE"})
    df = df.merge(water_df, on="date", how="left")
    df[water_col] = df[water_col].fillna(0)

    events = detect_events(df, water_col=water_col, threshold=threshold,
                            min_window=4, max_window=14)
    if not events:
        return dict(site=site, valid=False, reason="no events found")

    le_inf = estimate_le_inf(events)
    res = fit_with_validation(events, le_inf, n_boot=n_boot)
    res["site"] = site
    res["n_events"] = len(events)
    res["le_inf"] = le_inf
    res["water_col"] = water_col
    res["threshold"] = threshold
    return res


# ================================================================
# Figures
# ================================================================

def fig_timeseries(sat_8day: pd.DataFrame, water_dfs: dict,
                    out_path: Path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(13, 6), sharex=True)
    for ax, (site, _) in zip(axes, SITES.items()):
        sub = sat_8day[sat_8day["site"] == site].sort_values("date")
        ax.plot(sub["date"], sub["ET_mm_per_day"], "o-", ms=3,
                color="tab:purple", label=f"MOD16A2 ET ({site})")
        ax.set_ylabel("ET [mm/day]")
        ax.set_title(f"{site}: MOD16A2 ET (8-day composite)")
        ax.grid(alpha=0.3)

        # overlay water events
        wcol = "Irrig_mm" if site == "Tarazona" else "Rain_mm"
        threshold = IRRIG_THRESHOLD if site == "Tarazona" else RAIN_THRESHOLD
        wdf = water_dfs.get(site)
        if wdf is not None:
            ev = wdf[wdf[wcol] > threshold]
            ax2 = ax.twinx()
            ax2.bar(ev["date"], ev[wcol], width=1, alpha=0.4,
                     color="tab:blue", label=f"{wcol}")
            ax2.set_ylabel(f"{wcol} [mm/day]", color="tab:blue")
        ax.legend(loc="upper left")
    axes[-1].set_xlabel("date")
    fig.suptitle("MOD16A2 ET time series with detected water events",
                  fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def fig_recovery(sat_daily: pd.DataFrame, water_dfs: dict,
                  fit_results: dict, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, site in zip(axes, ["Oran", "Tarazona"]):
        res = fit_results.get(site, {})
        wcol = "Irrig_mm" if site == "Tarazona" else "Rain_mm"
        threshold = IRRIG_THRESHOLD if site == "Tarazona" else RAIN_THRESHOLD
        df = sat_daily[sat_daily["site"] == site][["date", "LE_Wm2"]].copy()
        df = df.rename(columns={"LE_Wm2": "LE"})
        wdf = water_dfs.get(site)
        if wdf is not None:
            df = df.merge(wdf, on="date", how="left")
            df[wcol] = df[wcol].fillna(0)
        events = detect_events(df, water_col=wcol, threshold=threshold,
                                min_window=4, max_window=14)

        # plot all events overlaid
        for ev in events:
            ev_d = ev["data"]
            ax.plot(ev_d["days_since_event"], ev_d["LE"], "o", ms=4,
                     alpha=0.4, color="gray")

        # fit curve if valid
        if res.get("valid") and not np.isnan(res.get("tau", np.nan)):
            tau, le0, le_inf = res["tau"], res["le0"], res["le_inf"]
            xx = np.linspace(0, 14, 200)
            yy = le_inf + (le0 - le_inf) * np.exp(-xx / tau)
            ax.plot(xx, yy, "-", color="tab:red", lw=2,
                     label=f"fit: τ={tau:.2f}d, LE_0={le0:.1f}, "
                           f"LE_inf={le_inf:.1f}\namp={res['amplitude']:.1f}, "
                           f"n={res['n_events']} events, R²={res['r2']:.2f}")
        else:
            reason = res.get("reason", "no fit")
            ax.text(0.05, 0.95, f"[no valid fit]\n{reason}",
                     transform=ax.transAxes, va="top")
        ax.set_xlabel("days since event")
        ax.set_ylabel("LE [W/m²] (satellite)")
        ax.set_title(f"{site}")
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right")
    fig.suptitle("MOD16A2 ET recovery after water events", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def fig_tau_compare(satellite_fits: dict, ec_taus: dict, out_path: Path) -> None:
    """Side-by-side comparison: EC τ vs satellite τ at each site."""
    fig, ax = plt.subplots(figsize=(8, 5))
    sites = list(SITES.keys())
    x = np.arange(len(sites))
    width = 0.35

    ec_vals = [ec_taus.get(s, {}).get("tau", np.nan) for s in sites]
    ec_ses  = [ec_taus.get(s, {}).get("se",  np.nan) for s in sites]
    sat_vals = [satellite_fits.get(s, {}).get("tau", np.nan) for s in sites]
    sat_ses  = [satellite_fits.get(s, {}).get("tau_se", np.nan) for s in sites]

    ax.bar(x - width/2, ec_vals, width, yerr=ec_ses, capsize=4,
           color="tab:blue", alpha=0.8, label="EC tower (analysis A v27)")
    ax.bar(x + width/2, sat_vals, width, yerr=sat_ses, capsize=4,
           color="tab:purple", alpha=0.8, label="MOD16A2 satellite (this work)")

    for i, (e, s) in enumerate(zip(ec_vals, sat_vals)):
        if np.isfinite(e):
            ax.text(i - width/2, e + 0.3, f"{e:.2f}",
                     ha="center", fontsize=10)
        if np.isfinite(s):
            ax.text(i + width/2, s + 0.3, f"{s:.2f}",
                     ha="center", fontsize=10)

    ax.set_xticks(x)
    ax.set_xticklabels(sites)
    ax.set_ylabel("τ [days]")
    ax.set_title("Recovery time τ: EC tower vs MOD16A2 satellite")
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

    # ── 1) GeoTIFF → 8-day pixel time series ─────────────────────────────
    print("\n=== 1) reading MOD16A2 GeoTIFF ===")
    print(f"  scale = {args.scale} (use 1.0 for GEE, 0.1 for raw HDF)")
    sat_8day = read_mod16_tif(Path(args.tif), args.start_date, args.cadence,
                                scale=args.scale)
    if sat_8day.empty:
        sys.exit("[FATAL] no satellite data extracted")
    print(f"  extracted {len(sat_8day)} (band x site) rows")
    print(f"  non-null  : {sat_8day['ET_mm_per_day'].notna().sum()}")
    print(sat_8day.groupby("site")["ET_mm_per_day"].describe().round(3))
    # sanity-check the magnitude
    mean_et = sat_8day.groupby("site")["ET_mm_per_day"].mean()
    for site, m in mean_et.items():
        if not np.isfinite(m):
            continue
        if site == "Tarazona" and m < 0.5:
            print(f"  [WARN] {site} mean ET = {m:.3f} mm/d looks low "
                  "(expect 1-4). Try --scale 0.1 or check TIFF units.")
        if site == "Oran" and m > 5:
            print(f"  [WARN] {site} mean ET = {m:.3f} mm/d looks high "
                  "(expect 0.5-2). Try --scale 1.0 or check TIFF units.")

    # ── 2) forward-fill 8-day → daily ────────────────────────────────────
    sat_daily = expand_to_daily(sat_8day)
    print(f"  daily-expanded: {len(sat_daily):,} rows")

    # ── 3) load water-input series for event detection ───────────────────
    print("\n=== 2) loading water-input series ===")
    tara_irrig = load_tara_irrig(Path(args.tara_csv))
    oran_rain  = load_oran_rain(Path(args.oran_csv))
    print(f"  Tarazona Irrig events (>{IRRIG_THRESHOLD}): "
          f"{(tara_irrig['Irrig_mm'] > IRRIG_THRESHOLD).sum()}")
    print(f"  Oran     Rain  events (>{RAIN_THRESHOLD}):  "
          f"{(oran_rain['Rain_mm'] > RAIN_THRESHOLD).sum()}")
    water_dfs = {"Tarazona": tara_irrig, "Oran": oran_rain}

    # ── 4) fit τ on satellite-derived LE at each site ────────────────────
    print("\n=== 3) fitting τ on satellite LE ===")
    sat_fits = {}
    for site in SITES:
        wcol = "Irrig_mm" if site == "Tarazona" else "Rain_mm"
        threshold = IRRIG_THRESHOLD if site == "Tarazona" else RAIN_THRESHOLD
        print(f"  -- {site} ({wcol} > {threshold}) --")
        res = fit_satellite_tau(sat_daily, water_dfs[site], site,
                                 wcol, threshold, n_boot=n_boot)
        sat_fits[site] = res
        if res.get("valid"):
            print(f"    τ = {res['tau']:.2f} d  "
                  f"CI [{res['ci_lo']:.2f}, {res['ci_hi']:.2f}]  "
                  f"amp = {res['amplitude']:.1f} W/m²  "
                  f"R² = {res['r2']:.3f}  "
                  f"n_events = {res['n_events']}")
        else:
            print(f"    [invalid] {res.get('reason', 'unknown')}")

    # ── 5) (optional) compare with EC τ from analysis A ──────────────────
    ec_taus = {
        "Oran":     {"tau": 2.82, "se": 0.90, "source": "analysis A v27 (active pool)"},
        "Tarazona": {"tau": 3.36, "se": 0.62, "source": "analysis A v27 (Jun-Sep)"},
    }
    print("\n=== 4) EC reference (from analysis A v27) ===")
    for s, v in ec_taus.items():
        print(f"  {s}: τ = {v['tau']:.2f} d (SE {v['se']:.2f}) — {v['source']}")

    # ── 6) write summary CSVs ────────────────────────────────────────────
    rows = []
    for site in SITES:
        sf = sat_fits[site]
        ef = ec_taus[site]
        # MDE-style comparison
        if sf.get("valid"):
            se1 = ef["se"]
            se2 = sf.get("tau_se", np.nan)
            mde = 1.96 * np.sqrt(se1**2 + se2**2) if np.isfinite(se2) else np.nan
            diff = abs(ef["tau"] - sf["tau"])
            sig = (diff > mde) if np.isfinite(mde) else None
            rows.append({
                "site": site,
                "tau_ec":       ef["tau"], "se_ec":  ef["se"],
                "tau_sat":      sf["tau"], "se_sat": sf.get("tau_se", np.nan),
                "ci_lo_sat":    sf.get("ci_lo", np.nan),
                "ci_hi_sat":    sf.get("ci_hi", np.nan),
                "n_events_sat": sf.get("n_events", 0),
                "amp_sat":      sf.get("amplitude", np.nan),
                "le_inf_sat":   sf.get("le_inf", np.nan),
                "r2_sat":       sf.get("r2", np.nan),
                "diff":         diff,
                "MDE":          mde,
                "significant_diff":  sig,
            })
        else:
            rows.append({
                "site": site,
                "tau_ec": ef["tau"], "se_ec": ef["se"],
                "tau_sat": np.nan, "se_sat": np.nan,
                "n_events_sat": sf.get("n_events", 0),
                "reason": sf.get("reason", "no fit"),
            })
    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "v1_satellite_vs_ec_tau.csv", index=False)
    pd.set_option("display.float_format", "{:+.3f}".format)
    print("\n=== SUMMARY ===")
    print(summary.to_string(index=False))
    print(f"\nwrote {out_dir / 'v1_satellite_vs_ec_tau.csv'}")

    # also save raw satellite series
    sat_8day.to_csv(out_dir / "v1_mod16_8day_pixel_series.csv", index=False)
    sat_daily.to_csv(out_dir / "v1_mod16_daily_expanded.csv", index=False)

    # ── 7) figures ───────────────────────────────────────────────────────
    if args.no_plots:
        print("[skipped plots --no-plots]")
        return

    print("\n=== 5) generating figures ===")
    fig_timeseries(sat_8day, water_dfs, out_dir / "fig02_satellite_timeseries.png")
    print(f"  wrote fig02_satellite_timeseries.png")
    fig_recovery(sat_daily, water_dfs, sat_fits,
                  out_dir / "fig03_satellite_recovery.png")
    print(f"  wrote fig03_satellite_recovery.png")
    fig_tau_compare(sat_fits, ec_taus, out_dir / "fig04_tau_comparison.png")
    print(f"  wrote fig04_tau_comparison.png")

    print(f"\n[DONE] outputs in {out_dir.absolute()}")


if __name__ == "__main__":
    main()
