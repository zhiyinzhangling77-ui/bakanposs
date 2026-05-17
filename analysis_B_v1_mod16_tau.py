"""
analysis_B_v1_mod16_tau.py
==========================

解析B v1: 衛星 ET (MOD16A2 v061) で 解析A の τ ≈ 3-4 d を独立検証する。
旧 v1 の 2 つの問題を修正:
  (1) スケール過剰適用 → AppEEARS CSV は raw int (要 *0.1)、
                          GEE TIFF は scale 適用済み (要 *1.0)
                          → --scale フラグで明示制御 + magnitude sanity-check
  (2) 旧 TIFF は 2018-2021 のみ → 新 AppEEARS CSV は 2018-2024 フルカバー

入力:
  --csv : NASA AppEEARS Sample 抽出 CSV
          /mnt/hdd/Dataset/MOD16A2_2018_2024_Oran_TzM/
            MOD16A2-2018-2024-Oran-TzM-MOD16A2-061-results.csv

          想定カラム (AppEEARS 標準):
            ID                            : サイト名 (Oran / Tarazona, etc.)
            Latitude, Longitude
            Date                          : YYYY-MM-DD
            MOD16A2_061_ET_500m           : raw ET (要 *0.1 で kg/m^2/8day)
            MOD16A2_061_ET_QC_500m        : QC bit-field
            MOD16A2_061_PET_500m          : PET (使わない)

  --tara-csv :  Tarazona EC daily CSV (要 Irrig_mm)
  --oran-csv :  Oran     EC daily CSV (要 Rain_mm)

  EC reference (解析A v27 の τ をベタ書き、比較用):
    Oran     : τ = 2.82 d (SE 0.90)
    Tarazona : τ = 3.36 d (SE 0.62)

処理:
  1) AppEEARS CSV 読込 → サイト別に長表に
  2) 単位変換: raw * scale / 8.0 = mm/day, LE = ET * 28.36 [W/m^2]
     - --scale 0.1 で AppEEARS raw を補正 (デフォルト)
     - QC で品質悪 (cloud / fill) をマスク
  3) Magnitude sanity-check: Tarazona summer mean が < 0.5 なら WARN
  4) 8-day → daily forward-fill (上限 7 日)
  5) Water-input (Irrig / Rain) を読み event 検出
  6) 自前の inline 指数減衰 fit + bootstrap
  7) τ_satellite vs τ_EC pairwise 比較
  8) 出力: CSV 1 + 図 3 枚

出力先: ./output_analysis_B_v1/
  v1_satellite_vs_ec_tau.csv
  v1_mod16_8day_pixel_series.csv
  v1_mod16_daily_expanded.csv
  fig01_satellite_timeseries.png
  fig02_satellite_recovery.png
  fig03_tau_comparison.png

Usage:
  python3 analysis_B_v1_mod16_tau.py
  python3 analysis_B_v1_mod16_tau.py --quick
  python3 analysis_B_v1_mod16_tau.py --scale 0.1   # AppEEARS raw (default)
  python3 analysis_B_v1_mod16_tau.py --scale 1.0   # already-scaled file
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

warnings.filterwarnings("ignore")


# ================================================================
# CONFIG
# ================================================================

# サイト座標 (Oran rainfed / Tarazona drip-irrigated almond)
SITES = {
    "Oran":     {"lat": 38.82,   "lon": -1.86},
    "Tarazona": {"lat": 39.266,  "lon": -1.9397},
}

# Site name aliases (AppEEARS ID column may differ in case/format)
SITE_ALIASES = {
    "Oran":     ["Oran", "oran", "ORAN"],
    "Tarazona": ["Tarazona", "TzM", "tzm", "TARAZONA", "tarazona",
                 "Tarazona_de_la_Mancha", "TarazonaDeLaMancha",
                 "Tarazona-de-la-Mancha"],
}

# AppEEARS で raw → mm/8day に変換するスケール (HDF native = 0.1)
# GEE 経由ですでに 0.1 適用済みなら 1.0 に
DEFAULT_MOD16_SCALE = 0.1

# kg/m^2/day = mm/day  →  W/m^2 (λρw / 86400 ≈ 28.36)
LE_PER_MM = 28.36

# AppEEARS は QC を生 bit-field の int で渡す。MOD16 ET QA の下位 2 bit が
# Mandatory QC: 0=good, 1=marginal, 2=cloud, 3=missing
MOD16_QA_GOOD = {0, 1}    # 1 は marginal だが運用上採用 (NASA 推奨)

# fill values (raw int)
MOD16_FILL = {32767, 32766, 32765, 32764, 32763, 32762, 32761}

# Event thresholds (解析A v27 と整合)
IRRIG_THRESHOLD = 0.5    # mm/day  → Tarazona drip event
RAIN_THRESHOLD  = 3.0    # mm/day  → Oran rain event

# τ fit constraints
TAU_FLOOR = 0.5
TAU_CEIL  = 60.0
TAU_PEG_RATIO = 0.9
R2_MIN_VALID = 0.30
N_EVENTS_MIN_VALID = 3
MIN_PER_BIN = 3
CI_PCT = (2.5, 97.5)

# 解析A v27 の最終 τ (このスクリプトと並べる)
EC_TAU_REFERENCE = {
    "Oran":     {"tau": 2.82, "se": 0.90,
                 "source": "analysis A v27 (active pool, n=9 rain events)"},
    "Tarazona": {"tau": 3.36, "se": 0.62,
                 "source": "analysis A v27 (Jun-Sep, n=41 irrig events)"},
}


def parse_args():
    p = argparse.ArgumentParser(
        description="解析B v1 — MOD16A2 (AppEEARS CSV) で τ 独立検証",
    )
    p.add_argument("--csv",
                   default="/mnt/hdd/Dataset/MOD16A2_2018_2024_Oran_TzM/"
                           "MOD16A2-2018-2024-Oran-TzM-MOD16A2-061-results.csv",
                   help="AppEEARS results CSV (point sample export)")
    p.add_argument("--tara-csv",
                   default="/home/shion-nagamine/Dataset/Eddy data in Spain/"
                           "Daily_Summary_Filtered_forPred_ActEne26.csv")
    p.add_argument("--oran-csv",
                   default="/home/shion-nagamine/Dataset/Eddy data in Spain/"
                           "Oran_EddyDaily_MASTER_2018_2020_correct.csv")
    p.add_argument("--out", default="./output_analysis_B_v1")
    p.add_argument("--scale", type=float, default=DEFAULT_MOD16_SCALE,
                   help=("MOD16 ET scale. 0.1 for AppEEARS raw int (default), "
                         "1.0 if values are already in mm/8day."))
    p.add_argument("--no-qc-filter", action="store_true",
                   help="skip QC masking even if QC column exists")
    p.add_argument("--n-boot", type=int, default=2000)
    p.add_argument("--quick", action="store_true",
                   help="n_boot=300, faster sanity-check run")
    p.add_argument("--no-plots", action="store_true")
    return p.parse_args()


# ================================================================
# AppEEARS CSV reader
# ================================================================

def _find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    """Case-insensitive lookup of the first column whose name matches any
    candidate (substring match)."""
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    # substring fallback
    for cand in candidates:
        for lc, orig in lower.items():
            if cand.lower() in lc:
                return orig
    return None


def _resolve_site(label: str) -> str | None:
    """Map AppEEARS ID label to one of SITES keys."""
    if not isinstance(label, str):
        return None
    for canonical, aliases in SITE_ALIASES.items():
        for a in aliases:
            if a.lower() == label.strip().lower():
                return canonical
    return None


def read_mod16_appeears(csv_path: Path, scale: float,
                          apply_qc: bool) -> pd.DataFrame:
    """Read NASA AppEEARS point-sample CSV → long DataFrame (date, site,
    ET_mm_per_day, LE_Wm2, qc, raw_value)."""
    if not csv_path.exists():
        sys.exit(f"[FATAL] CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)
    print(f"opened {csv_path.name}")
    print(f"  raw rows  : {len(df):,}")
    print(f"  raw cols  : {list(df.columns)[:6]}{'...' if len(df.columns) > 6 else ''}")

    id_col   = _find_col(df, ["ID", "Sample", "Name"])
    date_col = _find_col(df, ["Date"])
    et_col   = _find_col(df, ["MOD16A2_061_ET_500m", "ET_500m",
                                "MOD16A2_006_ET_500m"])
    qc_col   = _find_col(df, ["MOD16A2_061_ET_QC_500m", "ET_QC_500m"])

    if not all([id_col, date_col, et_col]):
        sys.exit(f"[FATAL] required cols missing: "
                 f"id={id_col}, date={date_col}, et={et_col}\n"
                 f"  columns present: {list(df.columns)}")
    print(f"  using cols: id={id_col!r}, date={date_col!r}, "
          f"et={et_col!r}, qc={qc_col!r}")

    df = df.rename(columns={id_col: "raw_id", date_col: "raw_date",
                              et_col: "raw_et"})
    df["date"] = pd.to_datetime(df["raw_date"], errors="coerce")
    df["site"] = df["raw_id"].map(_resolve_site)
    df = df.dropna(subset=["date", "site"])

    # date range diagnostic
    print(f"  date range : {df['date'].min().date()} → "
          f"{df['date'].max().date()}")
    print(f"  site counts: "
          f"{dict(df['site'].value_counts())}")

    # value processing
    raw = pd.to_numeric(df["raw_et"], errors="coerce")
    et_mask = raw.isin(MOD16_FILL) | raw.isna()
    et_8day = raw * scale                # raw → mm/8day (or already, if scale=1.0)
    et_8day = et_8day.where(~et_mask, np.nan)

    # QC masking
    if qc_col and apply_qc:
        qc = pd.to_numeric(df[qc_col], errors="coerce")
        qa_main = (qc.astype("Int64") & 0b11).astype("Int64")
        bad = ~qa_main.isin(MOD16_QA_GOOD)
        nbad = int(bad.fillna(True).sum())
        et_8day = et_8day.where(~bad.fillna(True), np.nan)
        print(f"  QC masked  : {nbad}/{len(df)} obs flagged")
    else:
        if qc_col:
            print(f"  [note] QC column present but --no-qc-filter set")

    out = pd.DataFrame({
        "date": df["date"].values,
        "site": df["site"].values,
        "ET_mm_per_day": et_8day.values / 8.0,
        "LE_Wm2":        et_8day.values / 8.0 * LE_PER_MM,
        "raw_value":     raw.values,
    })

    # dedup on (date,site): AppEEARS may have multiple lat/lons per site
    out = (out.sort_values(["site", "date"])
              .drop_duplicates(subset=["site", "date"])
              .reset_index(drop=True))

    print(f"  retained   : {len(out):,} (site,date) rows  "
          f"valid ET: {out['ET_mm_per_day'].notna().sum():,}")
    print("  per-site summary (mm/day):")
    print(out.groupby("site")["ET_mm_per_day"].describe().round(3))
    return out


def expand_to_daily(eight_day_df: pd.DataFrame) -> pd.DataFrame:
    """8-day composite を 8 日間 forward-fill (上限 7 日)。
    各 composite はその開始日からの 8 日間値とみなす。"""
    pieces = []
    for site, sub in eight_day_df.groupby("site"):
        sub = sub.sort_values("date").reset_index(drop=True)
        if sub.empty:
            continue
        idx = pd.date_range(sub["date"].min(),
                              sub["date"].max() + pd.Timedelta(days=7),
                              freq="D")
        s = sub.set_index("date").reindex(idx)
        for col in ["ET_mm_per_day", "LE_Wm2"]:
            s[col] = s[col].ffill(limit=7)
        s["site"] = site
        s = s.reset_index().rename(columns={"index": "date"})
        pieces.append(s[["date", "site", "ET_mm_per_day", "LE_Wm2"]])
    return pd.concat(pieces, ignore_index=True)


# ================================================================
# Water-input loaders (event detection inputs)
# ================================================================

def load_tara_irrig(path: Path) -> pd.DataFrame:
    if not path.exists():
        print(f"  [warn] Tarazona water CSV missing: {path}")
        return pd.DataFrame(columns=["date", "Irrig_mm"])
    df = pd.read_csv(path)
    date_col = _find_col(df, ["date", "Date", "TIMESTAMP"])
    irr_col  = _find_col(df, ["Irrig_mm", "Irrigation_mm", "Irrig"])
    if not all([date_col, irr_col]):
        print(f"  [warn] Tarazona Irrig cols missing in {path}")
        return pd.DataFrame(columns=["date", "Irrig_mm"])
    out = pd.DataFrame({
        "date":     pd.to_datetime(df[date_col], errors="coerce"),
        "Irrig_mm": pd.to_numeric(df[irr_col], errors="coerce"),
    }).dropna(subset=["date"]).drop_duplicates("date").sort_values("date")
    return out.reset_index(drop=True)


def load_oran_rain(path: Path) -> pd.DataFrame:
    if not path.exists():
        print(f"  [warn] Oran water CSV missing: {path}")
        return pd.DataFrame(columns=["date", "Rain_mm"])
    df = pd.read_csv(path)
    date_col = _find_col(df, ["Date", "date", "TIMESTAMP"])
    rain_col = _find_col(df, ["Rain_mm", "P_mm", "Precip_mm", "Precipitation"])
    if not all([date_col, rain_col]):
        print(f"  [warn] Oran rain cols missing in {path}")
        return pd.DataFrame(columns=["date", "Rain_mm"])
    out = pd.DataFrame({
        "date":    pd.to_datetime(df[date_col], errors="coerce"),
        "Rain_mm": pd.to_numeric(df[rain_col], errors="coerce"),
    }).dropna(subset=["date"]).drop_duplicates("date").sort_values("date")
    return out.reset_index(drop=True)


# ================================================================
# Inline exponential decay fit (mirrors analysis_A_v27 logic)
# ================================================================

def exp_model_2p(d, le0, tau, le_inf_fixed):
    return le_inf_fixed + (le0 - le_inf_fixed) * np.exp(-d / tau)


def fit_pooled_tau(arr_d, arr_y, le_inf_fixed,
                    min_per_bin: int = MIN_PER_BIN):
    """Median-per-day NLS fit with τ ∈ [TAU_FLOOR, TAU_CEIL]."""
    df = pd.DataFrame({"d": arr_d, "y": arr_y}).dropna()
    if df.empty:
        return None
    g = df.groupby("d")["y"].agg(["median", "count"]).reset_index()
    g = g[g["count"] >= min_per_bin]
    if len(g) < 3:
        return None
    x = g["d"].values.astype(float)
    y = g["median"].values
    le0_floor = max(le_inf_fixed + 1, 10.0)
    le0_ceil  = max(y.max() * 1.5, 350)

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
                          n_bins=len(g),
                          amplitude=popt[0] - le_inf_fixed)
        except Exception:
            continue
    return None


def bootstrap_taus(arr_d, arr_y, le_inf, n_boot=2000, seed=42):
    rng = np.random.default_rng(seed)
    n = len(arr_d)
    taus = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        res = fit_pooled_tau(arr_d[idx], arr_y[idx], le_inf)
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


def is_fit_valid(res, n_events):
    if n_events < N_EVENTS_MIN_VALID:
        return False, f"insufficient n_events ({n_events} < {N_EVENTS_MIN_VALID})"
    if res is None:
        return False, "fit returned None"
    if np.isnan(res.get("tau", np.nan)):
        return False, "tau is NaN"
    tau = res["tau"]
    if tau >= TAU_CEIL * TAU_PEG_RATIO:
        return False, f"tau={tau:.1f} pegged at ceiling ({TAU_CEIL})"
    if tau <= TAU_FLOOR * 1.5:
        return False, f"tau={tau:.2f} pegged at floor ({TAU_FLOOR})"
    if res.get("r2", -np.inf) < R2_MIN_VALID:
        return False, f"R2={res['r2']:.3f} < {R2_MIN_VALID}"
    return True, "OK"


def detect_events(df, water_col, threshold, min_window=4, max_window=14):
    df = df.sort_values("date").reset_index(drop=True).copy()
    event_dates = df[df[water_col].fillna(0) > threshold]["date"].values
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
        evd = df[(df["date"] >= ed) & (df["date"] <= window_end)].copy()
        evd["days_since_event"] = (evd["date"] - ed).dt.days
        evd = evd[evd["LE"].notna()]
        if len(evd) < min_window:
            continue
        events.append(dict(
            event_start=ed, data=evd,
            pulse_mm=float(df.loc[df["date"] == ed, water_col].iloc[0]),
        ))
    return events


def estimate_le_inf(events, tail_threshold=7):
    if not events:
        return np.nan
    pooled = pd.concat([ev["data"] for ev in events])
    mask = pooled["days_since_event"] >= tail_threshold
    if mask.sum() < 5:
        return float(np.percentile(pooled["LE"], 25))
    return float(np.median(pooled.loc[mask, "LE"]))


def fit_with_validation(events, le_inf, n_boot=2000):
    if len(events) < N_EVENTS_MIN_VALID:
        return dict(valid=False, success=False,
                     reason=f"n_events={len(events)} < {N_EVENTS_MIN_VALID}",
                     n_events=len(events), n_data=0,
                     tau=np.nan, le0=np.nan, amplitude=np.nan, r2=np.nan,
                     le_inf=le_inf,
                     ci_lo=np.nan, ci_hi=np.nan, tau_se=np.nan,
                     n_boot_valid=0)
    pooled = pd.concat([ev["data"] for ev in events])
    arr_d = pooled["days_since_event"].values.astype(float)
    arr_y = pooled["LE"].values

    res = fit_pooled_tau(arr_d, arr_y, le_inf)
    is_valid, reason = is_fit_valid(res, len(events))
    out = dict(n_events=len(events), n_data=len(pooled), le_inf=le_inf,
                tau=res["tau"] if res else np.nan,
                le0=res["le0"] if res else np.nan,
                amplitude=res["amplitude"] if res else np.nan,
                r2=res["r2"] if res else np.nan,
                valid=is_valid, reason=reason)
    if is_valid:
        boot_taus = bootstrap_taus(arr_d, arr_y, le_inf, n_boot=n_boot)
        if len(boot_taus) >= 100:
            out["ci_lo"], out["ci_hi"] = np.percentile(boot_taus, CI_PCT)
            out["tau_se"] = float(np.std(boot_taus))
            out["n_boot_valid"] = len(boot_taus)
        else:
            out["ci_lo"] = out["ci_hi"] = out["tau_se"] = np.nan
            out["n_boot_valid"] = len(boot_taus)
    else:
        out["ci_lo"] = out["ci_hi"] = out["tau_se"] = np.nan
        out["n_boot_valid"] = 0
    return out


# ================================================================
# Per-site driver
# ================================================================

def fit_satellite_tau(sat_daily, water_df, site, water_col, threshold,
                       n_boot=2000):
    if sat_daily.empty:
        return dict(site=site, valid=False, reason="empty satellite series")
    df = (sat_daily[sat_daily["site"] == site][["date", "LE_Wm2"]]
            .rename(columns={"LE_Wm2": "LE"}))
    if water_df is not None and not water_df.empty:
        df = df.merge(water_df, on="date", how="left")
    if water_col not in df.columns:
        df[water_col] = 0.0
    df[water_col] = df[water_col].fillna(0)

    events = detect_events(df, water_col=water_col, threshold=threshold,
                             min_window=4, max_window=14)
    if not events:
        return dict(site=site, valid=False, reason="no events found",
                     n_events=0)
    le_inf = estimate_le_inf(events)
    res = fit_with_validation(events, le_inf, n_boot=n_boot)
    res.update(site=site, n_events=len(events), le_inf=le_inf,
                 water_col=water_col, threshold=threshold)
    return res


# ================================================================
# Figures
# ================================================================

def fig_timeseries(sat_8day, water_dfs, out_path):
    fig, axes = plt.subplots(2, 1, figsize=(13, 6), sharex=True)
    for ax, site in zip(axes, SITES):
        sub = sat_8day[sat_8day["site"] == site].sort_values("date")
        ax.plot(sub["date"], sub["ET_mm_per_day"], "o-", ms=3,
                 color="tab:purple", label=f"MOD16A2 ET ({site})")
        ax.set_ylabel("ET [mm/day]")
        ax.set_title(f"{site}: MOD16A2 ET (8-day composite)")
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
    fig.suptitle("MOD16A2 ET time series with water events", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def fig_recovery(sat_daily, water_dfs, fit_results, out_path):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, site in zip(axes, ["Oran", "Tarazona"]):
        res = fit_results.get(site, {})
        wcol = "Irrig_mm" if site == "Tarazona" else "Rain_mm"
        thr  = IRRIG_THRESHOLD if site == "Tarazona" else RAIN_THRESHOLD
        df = (sat_daily[sat_daily["site"] == site][["date", "LE_Wm2"]]
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
            ax.plot(evd["days_since_event"], evd["LE"], "o", ms=4,
                     alpha=0.4, color="gray")
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
        ax.set_ylabel("LE [W/m²] (satellite-derived)")
        ax.set_title(site)
        ax.grid(alpha=0.3)
        ax.legend(loc="upper right")
    fig.suptitle("MOD16A2 LE recovery curves", fontsize=12)
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
            color="tab:purple", alpha=0.8, label="MOD16A2 satellite (this work)")
    for i, (e, s) in enumerate(zip(ec_v, sat_v)):
        if np.isfinite(e):
            ax.text(i - w/2, e + 0.3, f"{e:.2f}", ha="center", fontsize=10)
        if np.isfinite(s):
            ax.text(i + w/2, s + 0.3, f"{s:.2f}", ha="center", fontsize=10)
    ax.set_xticks(x); ax.set_xticklabels(sites)
    ax.set_ylabel("τ [days]")
    ax.set_title("Recovery time τ : EC tower vs MOD16A2 satellite")
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

    # 1) AppEEARS CSV → site x date long table
    print("\n=== 1) reading MOD16A2 (AppEEARS CSV) ===")
    print(f"  scale = {args.scale}  "
          f"(0.1 = AppEEARS raw int, 1.0 = pre-scaled)")
    sat_8day = read_mod16_appeears(Path(args.csv),
                                       scale=args.scale,
                                       apply_qc=not args.no_qc_filter)
    if sat_8day.empty:
        sys.exit("[FATAL] no satellite rows extracted")

    # magnitude sanity-check
    mean_et = sat_8day.groupby("site")["ET_mm_per_day"].mean()
    for site, m in mean_et.items():
        if not np.isfinite(m):
            continue
        if site == "Tarazona" and m < 0.5:
            print(f"  [WARN] Tarazona mean ET = {m:.3f} mm/d looks low. "
                   f"Expected 1-4 mm/d. Try --scale 0.1 (currently {args.scale}).")
        if site == "Tarazona" and m > 8:
            print(f"  [WARN] Tarazona mean ET = {m:.3f} mm/d looks high. "
                   f"Try --scale 1.0 (currently {args.scale}).")
        if site == "Oran" and m > 5:
            print(f"  [WARN] Oran mean ET = {m:.3f} mm/d looks high. "
                   f"Expected 0.3-2 mm/d.")

    # 2) 8-day → daily ffill
    sat_daily = expand_to_daily(sat_8day)
    print(f"  daily-expanded: {len(sat_daily):,} rows")

    # 3) water-input series
    print("\n=== 2) loading water-input series ===")
    tara_water = load_tara_irrig(Path(args.tara_csv))
    oran_water = load_oran_rain(Path(args.oran_csv))
    print(f"  Tarazona Irrig events (>{IRRIG_THRESHOLD}): "
          f"{(tara_water.get('Irrig_mm', pd.Series([], dtype=float)) > IRRIG_THRESHOLD).sum()}")
    print(f"  Oran     Rain  events (>{RAIN_THRESHOLD}):  "
          f"{(oran_water.get('Rain_mm', pd.Series([], dtype=float)) > RAIN_THRESHOLD).sum()}")
    water_dfs = {"Tarazona": tara_water, "Oran": oran_water}

    # 4) τ-fit per site
    print("\n=== 3) fitting τ on satellite LE ===")
    sat_fits = {}
    for site in SITES:
        wcol = "Irrig_mm" if site == "Tarazona" else "Rain_mm"
        thr  = IRRIG_THRESHOLD if site == "Tarazona" else RAIN_THRESHOLD
        print(f"  -- {site} ({wcol} > {thr}) --")
        res = fit_satellite_tau(sat_daily, water_dfs[site], site,
                                  wcol, thr, n_boot=n_boot)
        sat_fits[site] = res
        if res.get("valid"):
            print(f"    τ = {res['tau']:.2f} d  "
                   f"CI [{res['ci_lo']:.2f}, {res['ci_hi']:.2f}]  "
                   f"amp = {res['amplitude']:.1f} W/m²  "
                   f"R² = {res['r2']:.3f}  "
                   f"n_events = {res['n_events']}")
        else:
            print(f"    [invalid] {res.get('reason', 'unknown')}")

    # 5) EC reference + summary
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
            tau_sat=sf.get("tau", np.nan),
            tau_se_sat=sf.get("tau_se", np.nan),
            ci_lo_sat=sf.get("ci_lo", np.nan),
            ci_hi_sat=sf.get("ci_hi", np.nan),
            n_events_sat=sf.get("n_events", 0),
            amp_sat=sf.get("amplitude", np.nan),
            le_inf_sat=sf.get("le_inf", np.nan),
            r2_sat=sf.get("r2", np.nan),
            valid_sat=sf.get("valid", False),
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

    summary.to_csv(out_dir / "v1_satellite_vs_ec_tau.csv", index=False)
    sat_8day.to_csv(out_dir / "v1_mod16_8day_pixel_series.csv", index=False)
    sat_daily.to_csv(out_dir / "v1_mod16_daily_expanded.csv", index=False)

    pd.set_option("display.float_format", "{:+.3f}".format)
    print("\n=== SUMMARY ===")
    cols_show = ["site", "tau_ec", "se_ec", "tau_sat", "tau_se_sat",
                  "n_events_sat", "amp_sat", "r2_sat",
                  "MDE", "diff", "significant_diff", "reason"]
    print(summary[cols_show].to_string(index=False))
    print(f"\nwrote {out_dir / 'v1_satellite_vs_ec_tau.csv'}")

    if args.no_plots:
        print("[skipped plots --no-plots]")
        return

    print("\n=== 5) generating figures ===")
    fig_timeseries(sat_8day, water_dfs,
                     out_dir / "fig01_satellite_timeseries.png")
    print("  wrote fig01_satellite_timeseries.png")
    fig_recovery(sat_daily, water_dfs, sat_fits,
                   out_dir / "fig02_satellite_recovery.png")
    print("  wrote fig02_satellite_recovery.png")
    fig_tau_compare(sat_fits, EC_TAU_REFERENCE,
                      out_dir / "fig03_tau_comparison.png")
    print("  wrote fig03_tau_comparison.png")

    print(f"\n[DONE] outputs in {out_dir.absolute()}")


if __name__ == "__main__":
    main()
