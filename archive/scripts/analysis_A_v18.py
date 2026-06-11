"""
================================================================
解析A v18 : III + 3種降水データ比較 + τ Diagnostic 拡張
================================================================

【v17 からの追加点】
  1. Inter-Irrigation Interval (III) の年別計算
       - 灌漑イベント間の日数を年別に集計
       - III_mean / III_median / III_p75 を算出
       - τ × III の散布図で「τ は灌漑間隔に規定される」仮説を検証

  2. 3種降水データの Tarazona 座標抽出 + rain gauge との比較
       - CHIRPS v2   (daily, 0.05°) : 衛星+rain gauge 融合
       - TerraClimate (monthly, ~4km): 長期気候ベースライン
       - ERA5         (hourly, ~31km): 大気再解析

  3. fig06_tau_iii_scatter.png  : τ × III 散布図 (年ラベル付き)
  4. fig07_precip_compare.png   : 3種降水 + rain gauge の年別比較
  5. fig08_tau_drivers.png      : τ × III × 降水量 の統合診断

【Tarazona 座標】
  lat = 39.085, lon = -1.358  (論文記載値を使用)
  → NC ファイルから最近傍グリッドを sel(method="nearest") で抽出

【ERA5 hourly → daily 変換】
  ERA5 の tp (total precipitation) は hourly 積算値 [m]。
  → 日次合計 [mm] = sum(tp * 1000) per day

【科学的位置付け】
  rain gauge が ground truth。
  CHIRPS / TerraClimate / ERA5 は:
    (a) rain gauge の空間代表性を評価する補助データ
    (b) 欠損年の補完候補
  として使う。τ × III の相関が強ければ、
  「τ は灌漑管理(スケジュール)に規定される」という主張が成立する。

【入力】
  daily_classified_v4.parquet
  Tarazona daily CSV (Irrig_mm, Rain_mm)
  CHIRPS  : chirps-v2.0.{year}.days_p05.nc  (2020–2024)
  TerraClimate: TerraClimate_ppt_{year}.nc   (2020–2024)
  ERA5    : ERA5_*_hourly_ppt.nc             (複数ファイル)

【出力】./output_analysis_A_v18/
  fig01–fig05 : v17 と同一
  fig06_tau_iii_scatter.png
  fig07_precip_compare.png
  fig08_tau_drivers.png
  v18_iii_by_year.csv
  v18_precip_comparison.csv
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from scipy import stats
from scipy.optimize import curve_fit

try:
    import xarray as xr
    HAS_XARRAY = True
except ImportError:
    HAS_XARRAY = False
    print("  [warn] xarray not found — NC extraction skipped")

warnings.filterwarnings("ignore")

# ================================================================
# パス設定
# ================================================================
INPUT_DIR      = Path("/home/shion-nagamine/bakanposs/analysis_A")
PARQUET        = INPUT_DIR / "daily_classified_v4.parquet"
TARA_DAILY_CSV = Path(
    "/home/shion-nagamine/Dataset/Eddy data in Spain/"
    "Daily_Summary_Filtered_forPred_ActEne26.csv")

CHIRPS_DIR   = Path("/home/shion-nagamine/Dataset/Chirps v2 for persipitation")
TERRA_DIR    = Path("/home/shion-nagamine/Dataset/TerraClimate_ppt")
ERA5_DIR     = Path("/home/shion-nagamine/Dataset/ERA5_ppt")

OUT_DIR = Path("./output_analysis_A_v18")
OUT_DIR.mkdir(exist_ok=True)

# ================================================================
# Tarazona 座標 (論文記載値)
# ================================================================
TARA_LAT =  39.085
TARA_LON =  -1.358
STUDY_YEARS = [2020, 2021, 2022, 2023, 2024]

# ================================================================
# パラメータ (v17 と同じ)
# ================================================================
IRRIG_THRESHOLD      = 0.5
MIN_IRRIG_PER_MONTH  = 2

SEASONS = {
    "spring (1-4月)"      : [1, 2, 3, 4],
    "shoulder (5,6,10月)" : [5, 6, 10],
    "summer (7-9月)"      : [7, 8, 9],
}

N_BOOT            = 5000
CI_PCT            = (2.5, 97.5)
MIN_N_PER_CLASS   = 5
DENOM_FLOORS      = {"LE_corr": 5.0, "EF_corr": 0.05, "ET": 0.3}
RATIO_CLIP        = 10.0
N_MIN_VERDICT     = 30
CI_WIDTH_MAX      = 0.5
DEEP_ROOT_BAND    = 0.15
R2_MIN_TRUST      = 0.70

SITE_COL = {"Oran": "#E85D04", "Tarazona": "#1D9E75"}
FIG_BG   = "#F8F9FA"
VERDICT_COL = {
    "deep_root"         : "#1D9E75",
    "shallow"           : "#E85D04",
    "negative_anomaly"  : "#7C3AED",
    "uncertain"         : "#9E9E9E",
    "uncertain_wide_CI" : "#BDBDBD",
    "insufficient_data" : "#E8E8E8",
}
VERDICT_JP = {
    "deep_root"         : "深根候補",
    "shallow"           : "浅根(表層感受性)",
    "negative_anomaly"  : "異常値",
    "uncertain"         : "判定不能",
    "uncertain_wide_CI" : "CI幅超過→保留",
    "insufficient_data" : "n不足→保留",
}

# ================================================================
# 1. データ読込 (v17 と同じ)
# ================================================================

def load_and_merge() -> pd.DataFrame:
    daily = pd.read_parquet(PARQUET)
    daily["date"] = pd.to_datetime(daily["date"])
    csv = pd.read_csv(TARA_DAILY_CSV)
    csv["date"] = pd.to_datetime(csv["date"], errors="coerce")
    irrig_cols = [c for c in ["Irrig_mm", "Rain_mm", "IrrigRain_mm"] if c in csv.columns]
    extra = (csv[["date"] + irrig_cols]
             .dropna(subset=["date"]).drop_duplicates("date"))
    tara = daily[daily["site"] == "Tarazona"].merge(extra, on="date", how="left")
    oran = daily[daily["site"] == "Oran"].copy()
    for c in irrig_cols:
        oran[c] = 0.0
    df = pd.concat([oran, tara]).sort_values(["site", "date"]).reset_index(drop=True)
    df["year"]  = df["date"].dt.year
    df["month"] = df["date"].dt.month
    return df


def add_days_since_irrig(df):
    df = df.copy()
    df["days_since_irrig"] = np.nan
    for site, idx in df.groupby("site").groups.items():
        sub = df.loc[idx].sort_values("date")
        last, out = None, []
        for _, row in sub.iterrows():
            irrig = row.get("Irrig_mm", 0) or 0
            if pd.notna(irrig) and irrig > IRRIG_THRESHOLD:
                last = row["date"]; out.append(0)
            elif last is None:
                out.append(np.nan)
            else:
                out.append((row["date"] - last).days)
        df.loc[sub.index, "days_since_irrig"] = out
    return df


def add_irrig_active_month(df):
    df = df.copy()
    df["_ym"] = df["date"].dt.to_period("M")
    df["_is_irrig"] = df["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD
    monthly = (df.groupby(["site", "_ym"])["_is_irrig"]
                 .sum().reset_index(name="n_irrig"))
    active_set = set(
        zip(monthly.loc[monthly["n_irrig"] >= MIN_IRRIG_PER_MONTH, "site"],
            monthly.loc[monthly["n_irrig"] >= MIN_IRRIG_PER_MONTH, "_ym"]))
    df["irrig_active_month"] = df.apply(
        lambda r: (r["site"], r["_ym"]) in active_set, axis=1)
    df.drop(columns=["_ym", "_is_irrig"], inplace=True)
    return df


# ================================================================
# 2. Inter-Irrigation Interval (III) の計算
# ================================================================

def calc_iii_by_year(df: pd.DataFrame) -> pd.DataFrame:
    """
    灌漑イベント間の日数 (Inter-Irrigation Interval) を年別に計算する。

    τ と III の相関を検証することで
    「τ は灌漑管理スケジュール(間隔)に規定される」仮説を評価できる。

    仮説:
      III が短い(頻繁に灌漑) → 植物が常に水分補給される
        → 灌漑停止後のLE低下が急峻 → τ 短い
      III が長い(間隔が大きい) → 植物が乾燥に適応/深層水に依存
        → 灌漑停止後も LE を維持 → τ 長い
    """
    tara = df[(df["site"] == "Tarazona") &
              (df["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD)].copy()
    rows = []
    for yr in STUDY_YEARS:
        sub = tara[tara["year"] == yr].sort_values("date")
        dates = sub["date"].values
        if len(dates) < 2:
            rows.append(dict(year=yr, iii_mean=np.nan, iii_median=np.nan,
                             iii_p25=np.nan, iii_p75=np.nan, n_events=len(dates)))
            continue
        intervals = np.array([(dates[i+1] - dates[i]) / np.timedelta64(1, "D")
                               for i in range(len(dates) - 1)])
        # 異常に長い間隔(冬季の灌漑停止期間)を除外: 60日超はシーズン境界とみなす
        intervals = intervals[intervals <= 60]
        if len(intervals) == 0:
            rows.append(dict(year=yr, iii_mean=np.nan, iii_median=np.nan,
                             iii_p25=np.nan, iii_p75=np.nan, n_events=len(dates)))
            continue
        rows.append(dict(
            year=yr,
            iii_mean=float(np.mean(intervals)),
            iii_median=float(np.median(intervals)),
            iii_p25=float(np.percentile(intervals, 25)),
            iii_p75=float(np.percentile(intervals, 75)),
            n_events=len(dates),
        ))
    return pd.DataFrame(rows)


# ================================================================
# 3. 3種降水データの抽出
# ================================================================

def _nearest_val(ds, lat_name, lon_name, var_name, lat, lon):
    """xarray Dataset から最近傍グリッドの時系列を返す"""
    return ds[var_name].sel(
        {lat_name: lat, lon_name: lon}, method="nearest")


def extract_chirps(years=STUDY_YEARS) -> pd.DataFrame:
    """
    CHIRPS v2 daily (0.05°) から Tarazona の年次降水量を抽出。
    変数名: precip [mm/day], 座標: latitude / longitude
    """
    if not HAS_XARRAY:
        return pd.DataFrame({"year": years, "chirps_mm": [np.nan]*len(years)})
    rows = []
    for yr in years:
        fp = CHIRPS_DIR / f"chirps-v2.0.{yr}.days_p05.nc"
        if not fp.exists():
            rows.append(dict(year=yr, chirps_mm=np.nan)); continue
        try:
            ds  = xr.open_dataset(fp)
            # 変数名の自動検出 (precip が標準)
            var = "precip" if "precip" in ds else list(ds.data_vars)[0]
            lat_name = "latitude"  if "latitude"  in ds.coords else "lat"
            lon_name = "longitude" if "longitude" in ds.coords else "lon"
            ts  = _nearest_val(ds, lat_name, lon_name, var, TARA_LAT, TARA_LON)
            # 年間合計 [mm]; CHIRPS は mm/day なのでそのまま sum
            annual = float(ts.values[ts.values >= 0].sum())
            rows.append(dict(year=yr, chirps_mm=round(annual, 1)))
            ds.close()
        except Exception as e:
            print(f"  [warn] CHIRPS {yr}: {e}")
            rows.append(dict(year=yr, chirps_mm=np.nan))
    return pd.DataFrame(rows)


def extract_terraclimate(years=STUDY_YEARS) -> pd.DataFrame:
    """
    TerraClimate monthly ppt [mm/month] から年次合計を抽出。
    変数名: ppt, 座標: lat / lon
    """
    if not HAS_XARRAY:
        return pd.DataFrame({"year": years, "terra_mm": [np.nan]*len(years)})
    rows = []
    for yr in years:
        fp = TERRA_DIR / f"TerraClimate_ppt_{yr}.nc"
        if not fp.exists():
            rows.append(dict(year=yr, terra_mm=np.nan)); continue
        try:
            ds  = xr.open_dataset(fp)
            var = "ppt" if "ppt" in ds else list(ds.data_vars)[0]
            lat_name = "lat" if "lat" in ds.coords else "latitude"
            lon_name = "lon" if "lon" in ds.coords else "longitude"
            ts  = _nearest_val(ds, lat_name, lon_name, var, TARA_LAT, TARA_LON)
            annual = float(ts.values[ts.values >= 0].sum())
            rows.append(dict(year=yr, terra_mm=round(annual, 1)))
            ds.close()
        except Exception as e:
            print(f"  [warn] TerraClimate {yr}: {e}")
            rows.append(dict(year=yr, terra_mm=np.nan))
    return pd.DataFrame(rows)


def extract_era5(years=STUDY_YEARS) -> pd.DataFrame:
    """
    ERA5 hourly total precipitation [m] → daily sum → 年次合計 [mm] を抽出。
    複数ファイルをまとめて処理し、対象年のみ抽出する。

    ERA5 の tp は hourly 積算値 [m]。
    日次合計 [mm] = sum(tp * 1000) per UTC day。
    座標名: latitude / longitude
    """
    if not HAS_XARRAY:
        return pd.DataFrame({"year": years, "era5_mm": [np.nan]*len(years)})

    era5_files = sorted(ERA5_DIR.glob("ERA5_*_hourly_ppt.nc"))
    if not era5_files:
        print("  [warn] ERA5 files not found")
        return pd.DataFrame({"year": years, "era5_mm": [np.nan]*len(years)})

    try:
        ds = xr.open_mfdataset(era5_files, combine="by_coords")
        var = "tp" if "tp" in ds else list(ds.data_vars)[0]
        lat_name = "latitude"  if "latitude"  in ds.coords else "lat"
        lon_name = "longitude" if "longitude" in ds.coords else "lon"
        ts = _nearest_val(ds, lat_name, lon_name, var, TARA_LAT, TARA_LON)

        # hourly [m] → daily [mm]
        daily_mm = (ts * 1000).resample(time="1D").sum()
        daily_df = daily_mm.to_dataframe(name="era5_mm").reset_index()
        daily_df["year"] = pd.to_datetime(daily_df["time"]).dt.year
        annual = (daily_df[daily_df["era5_mm"] >= 0]
                  .groupby("year")["era5_mm"].sum()
                  .reset_index(name="era5_mm"))
        annual["era5_mm"] = annual["era5_mm"].round(1)
        ds.close()
    except Exception as e:
        print(f"  [warn] ERA5 extraction failed: {e}")
        return pd.DataFrame({"year": years, "era5_mm": [np.nan]*len(years)})

    # 対象年のみ返す
    result = pd.DataFrame({"year": years})
    return result.merge(annual[annual["year"].isin(years)], on="year", how="left")


# ================================================================
# 4. 三点パッケージ + 判定 (v17 と同じ)
# ================================================================

def safe_ratio(numer, denom, var):
    floor = DENOM_FLOORS.get(var, 0.05)
    if abs(denom) < floor: return np.nan
    r = numer / denom
    return np.nan if abs(r) > RATIO_CLIP else r


def three_point_package(arr_n, arr_s, var, seed=42):
    n_n, n_s = len(arr_n), len(arr_s)
    empty = dict(sds=np.nan, sds_lo=np.nan, sds_hi=np.nan,
                 abs_diff=np.nan, abs_lo=np.nan, abs_hi=np.nan,
                 rb=np.nan, p=np.nan, med_n=np.nan, med_s=np.nan,
                 n_n=n_n, n_s=n_s)
    if n_n < MIN_N_PER_CLASS or n_s < MIN_N_PER_CLASS: return empty
    med_n, med_s = float(np.median(arr_n)), float(np.median(arr_s))
    sds_pt, abs_pt = safe_ratio(med_n - med_s, med_n, var), med_n - med_s
    rng = np.random.default_rng(seed)
    boot_sds, boot_abs = [], []
    for _ in range(N_BOOT):
        nb = rng.choice(arr_n, size=n_n, replace=True)
        sb = rng.choice(arr_s, size=n_s, replace=True)
        r = safe_ratio(np.median(nb) - np.median(sb), np.median(nb), var)
        if not np.isnan(r): boot_sds.append(r)
        boot_abs.append(np.median(nb) - np.median(sb))
    sds_lo = sds_hi = np.nan
    if len(boot_sds) >= N_BOOT * 0.1:
        sds_lo, sds_hi = np.percentile(boot_sds, CI_PCT)
    abs_lo, abs_hi = np.percentile(boot_abs, CI_PCT)
    try:
        u, p = stats.mannwhitneyu(arr_n, arr_s, alternative="two-sided")
        rb = 1 - 2 * u / (n_n * n_s)
    except ValueError:
        u = p = rb = np.nan
    return dict(sds=sds_pt, sds_lo=sds_lo, sds_hi=sds_hi,
                abs_diff=abs_pt, abs_lo=abs_lo, abs_hi=abs_hi,
                rb=rb, p=p, med_n=med_n, med_s=med_s, n_n=n_n, n_s=n_s)


def verdict(pkg):
    n_n, n_s = pkg.get("n_n", 0), pkg.get("n_s", 0)
    sds, lo, hi = pkg.get("sds", np.nan), pkg.get("sds_lo", np.nan), pkg.get("sds_hi", np.nan)
    if np.isnan(sds) or n_n < N_MIN_VERDICT or n_s < N_MIN_VERDICT:
        return "insufficient_data", "weak"
    if np.isnan(lo): return "uncertain", "weak"
    if hi - lo > CI_WIDTH_MAX: return "uncertain_wide_CI", "weak"
    if lo <= 0 <= hi and abs(sds) <= DEEP_ROOT_BAND: return "deep_root", "strong"
    if lo > 0: return "shallow", "strong"
    if hi < 0: return "negative_anomaly", "weak"
    return "uncertain", "weak"


# ================================================================
# 5. Recovery τ (v17 と同じ)
# ================================================================

def exp_model(d, le_inf, le0, tau):
    return le_inf + (le0 - le_inf) * np.exp(-d / tau)

def lin_model(d, a, b): return a + b * d
def log_model(d, a, b): return a + b * np.log(np.where(d > 0, d, 0.5))


def _bin_median(sub, var, min_per_bin=5):
    g = (sub.groupby("days_since_irrig")[var]
            .agg(["median", "count"]).reset_index())
    return g[g["count"] >= min_per_bin].reset_index(drop=True)


def fit_one(x, y, model_fn, p0, bounds):
    try:
        popt, _ = curve_fit(model_fn, x, y, p0=p0, bounds=bounds, maxfev=8000)
        y_pred = model_fn(x, *popt)
        ss_res = np.sum((y - y_pred) ** 2)
        ss_tot = np.sum((y - y.mean()) ** 2)
        r2  = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
        k, n = len(popt), len(y)
        aic = n * np.log(ss_res / n + 1e-12) + 2 * k if n > 0 else np.nan
        return popt, r2, aic
    except Exception:
        return [np.nan] * 3, np.nan, np.nan


def fit_recovery(df_active, var="LE_corr"):
    sub     = df_active[df_active[var].notna() & df_active["days_since_irrig"].notna()]
    grouped = _bin_median(sub, var)
    base    = dict(tau=np.nan, le0=np.nan, le_inf=np.nan,
                   r2_exp=np.nan, r2_lin=np.nan, r2_log=np.nan,
                   aic_exp=np.nan, aic_lin=np.nan, aic_log=np.nan,
                   n_points=len(grouped), grouped=grouped)
    if len(grouped) < 4: return base
    x, y = grouped["days_since_irrig"].values.astype(float), grouped["median"].values
    popt_e, r2_e, aic_e = fit_one(x, y, exp_model,
        p0=[y.min(), y.max(), 5.0], bounds=([0,0,0.5],[np.inf,np.inf,100]))
    popt_l, r2_l, aic_l = fit_one(x, y, lin_model,
        p0=[y.mean(),-1.0], bounds=([-np.inf,-np.inf],[np.inf,np.inf]))
    popt_g, r2_g, aic_g = fit_one(x, y, log_model,
        p0=[y.mean(),-5.0], bounds=([-np.inf,-np.inf],[np.inf,np.inf]))
    base.update(dict(
        tau=float(popt_e[2]) if not np.isnan(popt_e[2]) else np.nan,
        le0=float(popt_e[1]) if not np.isnan(popt_e[1]) else np.nan,
        le_inf=float(popt_e[0]) if not np.isnan(popt_e[0]) else np.nan,
        popt_exp=popt_e, popt_lin=popt_l, popt_log=popt_g,
        r2_exp=r2_e, r2_lin=r2_l, r2_log=r2_g,
        aic_exp=aic_e, aic_lin=aic_l, aic_log=aic_g,
    ))
    return base


def fit_tau_by_year(df_active, var="LE_corr"):
    rows = []
    for yr in STUDY_YEARS:
        sub = df_active[df_active["year"] == yr]
        res = fit_recovery(sub, var)
        trusted = (not np.isnan(res["r2_exp"]) and res["r2_exp"] >= R2_MIN_TRUST)
        reason  = ("OK" if trusted else
                   "fit失敗" if np.isnan(res["r2_exp"]) else
                   f"R²={res['r2_exp']:.2f} < {R2_MIN_TRUST}")
        rows.append(dict(year=yr, tau=res["tau"], r2=res["r2_exp"],
                         n_points=res["n_points"], tau_reliable=trusted, reason=reason))
    return pd.DataFrame(rows)


def calc_irrig_stats_by_year(df):
    tara = df[df["site"] == "Tarazona"].copy()
    tara["is_irrig"] = tara["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD
    rows = []
    for yr in STUDY_YEARS:
        sub = tara[tara["year"] == yr]
        n_ev    = int(sub["is_irrig"].sum())
        tot_mm  = float(sub["Irrig_mm"].fillna(0).sum())
        avg_mm  = tot_mm / n_ev if n_ev > 0 else np.nan
        rain    = float(sub["Rain_mm"].fillna(0).sum()) if "Rain_mm" in sub.columns else np.nan
        rows.append(dict(year=yr, n_irrig_events=n_ev,
                         total_irrig_mm=round(tot_mm, 1),
                         avg_mm_per_event=round(avg_mm, 1) if not np.isnan(avg_mm) else np.nan,
                         gauge_rain_mm=round(rain, 1) if not np.isnan(rain) else np.nan))
    return pd.DataFrame(rows)


# ================================================================
# 6. 可視化ヘルパー
# ================================================================

def _sig_star(p):
    if np.isnan(p): return ""
    return "★★★" if p<0.001 else "★★" if p<0.01 else "★" if p<0.05 else "n.s."

def _ax_style(ax, title="", xlabel="", ylabel=""):
    ax.set_title(title, fontweight="bold", fontsize=10, pad=6)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.grid(axis="both", alpha=0.2, lw=0.7)
    ax.spines[["top", "right"]].set_visible(False)


# ================================================================
# 7. fig03 — 三点指標間整合性
# ================================================================

def plot_three_point_matrix(pkg_df, save_dir):
    sub = pkg_df[pkg_df["var"]=="LE_corr"].dropna(subset=["sds","abs_diff","rb"])
    if len(sub) < 2:
        print("  [skip] three-point matrix"); return
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.patch.set_facecolor(FIG_BG)
    fig.suptitle("fig03 — 三点指標間の整合性 (LE_corr)\n"
                 "3指標が同じ方向を向いているほど主張が Robust",
                 fontsize=12, fontweight="bold")
    pairs = [
        ("sds","abs_diff","SDS","Abs diff [W/m²]","SDS vs 絶対差"),
        ("sds","rb","SDS","rank-biserial r","SDS vs 効果量"),
        ("abs_diff","rb","Abs diff [W/m²]","rank-biserial r","絶対差 vs 効果量"),
    ]
    for ax, (xk,yk,xl,yl,title) in zip(axes, pairs):
        for _, r in sub.iterrows():
            v, strength = verdict({"sds":r["sds"],"sds_lo":r["sds_lo"],
                                    "sds_hi":r["sds_hi"],"n_n":r["n_n"],"n_s":r["n_s"]})
            col   = SITE_COL.get(r["site"], "#888")
            mk    = "o" if r["site"]=="Oran" else "^"
            alpha = 0.9 if strength=="strong" else 0.4
            ax.scatter(r[xk], r[yk], c=col, s=90, alpha=alpha,
                       edgecolors="black", lw=1, marker=mk)
            ax.annotate(f"{r['site'][:3]}\n{r['season'][:6]}", (r[xk],r[yk]),
                        fontsize=7.5, xytext=(5,5), textcoords="offset points")
        ax.axhline(0, color="gray", lw=0.6, alpha=0.5)
        ax.axvline(0, color="gray", lw=0.6, alpha=0.5)
        _ax_style(ax, title, xl, yl)
    leg = [mpatches.Patch(color=SITE_COL["Oran"],     label="Oran (rainfed)"),
           mpatches.Patch(color=SITE_COL["Tarazona"], label="Tarazona (irrigated)"),
           Line2D([0],[0], marker="o", color="gray", ms=8, lw=0,
                  markeredgecolor="black", label="strong verdict"),
           Line2D([0],[0], marker="o", color="gray", ms=8, lw=0,
                  alpha=0.35, markeredgecolor="black", label="weak verdict")]
    fig.legend(handles=leg, loc="lower center", ncol=4, fontsize=9, frameon=True)
    plt.tight_layout(rect=[0,0.07,1,1])
    fp = save_dir / "fig03_three_point_matrix.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 8. fig04 — 年別 τ Standalone
# ================================================================

def plot_tau_sensitivity(tau_by_year, tau_all, save_dir):
    valid = tau_by_year.dropna(subset=["tau"])
    if len(valid) < 2:
        print("  [skip] tau sensitivity"); return
    fig, ax = plt.subplots(figsize=(10, 5.5))
    fig.patch.set_facecolor(FIG_BG)
    tmin, tmax = valid["tau"].min(), valid["tau"].max()
    colors = ["#1D9E75" if r.tau_reliable else "#BDBDBD" for r in valid.itertuples()]
    bars = ax.bar(valid["year"], valid["tau"], color=colors, alpha=0.80,
                  edgecolor="black", lw=1.2, width=0.6)
    ax.axhline(tau_all, color="#E85D04", ls="--", lw=2.2,
               label=f"全期間 τ={tau_all:.1f}d")
    ax.axhline(5, color="gray", ls=":", lw=1.0, alpha=0.7, label="τ=5d 閾値")
    trusted = valid[valid["tau_reliable"]]
    if len(trusted) >= 2:
        ax.axhspan(trusted["tau"].min(), trusted["tau"].max(),
                   color="#1D9E75", alpha=0.08,
                   label=f"信頼年: {trusted['tau'].min():.1f}–{trusted['tau'].max():.1f}d")
    for bar, (_, row) in zip(bars, valid.iterrows()):
        ax.text(bar.get_x()+bar.get_width()/2, row["tau"]+0.15,
                f"{row['tau']:.1f}d\n(n={int(row['n_points'])})\n{row['reason']}",
                ha="center", va="bottom", fontsize=8.5,
                color="black" if row["tau_reliable"] else "#999")
    ax.set_ylim(0, max(valid["tau"].max()*1.45, 8))
    ax.set_xticks(valid["year"])
    leg_patches = [mpatches.Patch(color="#1D9E75", alpha=0.85, label="信頼年 (R²≥0.7)"),
                   mpatches.Patch(color="#BDBDBD", alpha=0.85, label="低信頼年")]
    h, _ = ax.get_legend_handles_labels()
    ax.legend(handles=h+leg_patches, fontsize=9.5)
    r_flag = ("★ 信頼年 τ 変動小 → Robust"
              if len(trusted)>=2 and trusted["tau"].max()-trusted["tau"].min()<3
              else "△ 信頼年でも変動あり → 要注意")
    r_col = "#1D9E75" if "★" in r_flag else "#E85D04"
    body = (f"{r_flag}\n信頼年 τ 平均={trusted['tau'].mean():.1f}±{trusted['tau'].std():.1f}d"
            if len(trusted)>=2 else r_flag)
    ax.text(0.98, 0.98, body, transform=ax.transAxes, va="top", ha="right",
            fontsize=10, fontweight="bold", color=r_col,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=r_col, alpha=0.9))
    _ax_style(ax, f"fig04 — 年別 τ Sensitivity (Tarazona)\n全期間 τ={tau_all:.1f}d",
              "Year", "τ [days]")
    plt.tight_layout()
    fp = save_dir / "fig04_tau_sensitivity.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 9. fig05 — τ Diagnostic (灌漑統計との4段パネル)
# ================================================================

def plot_tau_diagnostic(tau_by_year, irrig_stats, save_dir):
    """
    年別 τ と灌漑統計を4段に並べて τ 変動の要因を診断する。
    (A) τ  (B) 灌漑イベント数  (C) mm/event  (D) rain gauge 降水量
    """
    merged = tau_by_year.merge(irrig_stats, on="year", how="left")
    years  = merged["year"].values
    x      = np.arange(len(years))
    w      = 0.6

    fig = plt.figure(figsize=(14, 12))
    fig.patch.set_facecolor(FIG_BG)
    fig.suptitle(
        "fig05 — τ Diagnostic: 年別 τ × 灌漑統計\n"
        "τ の年間変動の要因を診断する",
        fontsize=12, fontweight="bold", y=0.98)
    gs = GridSpec(4, 1, hspace=0.55, figure=fig)

    # (A) τ
    ax0 = fig.add_subplot(gs[0])
    tau_vals = merged["tau"].values
    colors0  = ["#1D9E75" if r.tau_reliable else "#BDBDBD" for r in merged.itertuples()]
    ax0.bar(x, tau_vals, color=colors0, alpha=0.85, edgecolor="black", lw=1, width=w)
    for i, (v, r) in enumerate(zip(tau_vals, merged.itertuples())):
        if not np.isnan(v):
            ax0.text(i, v+0.3, f"{v:.1f}d\n({r.reason})",
                     ha="center", va="bottom", fontsize=8,
                     color="black" if r.tau_reliable else "#888")
    ax0.axhline(5, color="gray", ls=":", lw=1, alpha=0.7, label="τ=5d 閾値")
    ax0.set_ylim(0, np.nanmax(tau_vals)*1.4+2)
    ax0.set_xticks(x); ax0.set_xticklabels(years)
    leg_a = [mpatches.Patch(color="#1D9E75", alpha=0.85, label="信頼年 (R²≥0.7)"),
             mpatches.Patch(color="#BDBDBD", alpha=0.85, label="低信頼年 (R²<0.7)")]
    ax0.legend(handles=leg_a+ax0.get_legend_handles_labels()[0], fontsize=8, loc="upper right")
    _ax_style(ax0, "(A) 年別 τ [days]", "", "τ [days]")

    # (B) 灌漑イベント数
    ax1 = fig.add_subplot(gs[1])
    n_ev = merged["n_irrig_events"].values
    ax1.bar(x, n_ev, color="#1A73E8", alpha=0.75, edgecolor="black", lw=1, width=w)
    for i, v in enumerate(n_ev):
        if not np.isnan(v):
            ax1.text(i, v+0.5, f"{int(v)}", ha="center", va="bottom", fontsize=9)
    ax1.set_xticks(x); ax1.set_xticklabels(years)
    _ax_style(ax1, "(B) 年別 灌漑イベント数\n多い → τ 短い可能性", "", "イベント数")

    # (C) mm/event
    ax2 = fig.add_subplot(gs[2])
    mm_ev = merged["avg_mm_per_event"].values
    ax2.bar(x, mm_ev, color="#FFA000", alpha=0.80, edgecolor="black", lw=1, width=w)
    for i, v in enumerate(mm_ev):
        if not np.isnan(v):
            ax2.text(i, v+0.3, f"{v:.1f}", ha="center", va="bottom", fontsize=9)
    ax2.set_xticks(x); ax2.set_xticklabels(years)
    _ax_style(ax2, "(C) 年別 灌漑量 [mm/event]\n多い → 1回で多くの水 → τ 長い可能性",
              "", "mm/event")

    # (D) 降水量
    ax3 = fig.add_subplot(gs[3])
    rain_col = "gauge_rain_mm" if "gauge_rain_mm" in merged.columns else "Rain_mm"
    rain = merged[rain_col].values if rain_col in merged.columns else np.full(len(years), np.nan)
    ax3.bar(x, rain, color="#9575CD", alpha=0.75, edgecolor="black", lw=1, width=w)
    for i, v in enumerate(rain):
        if not np.isnan(v):
            ax3.text(i, v+2, f"{v:.0f}", ha="center", va="bottom", fontsize=9)
    ax3.set_xticks(x); ax3.set_xticklabels(years)
    _ax_style(ax3, "(D) 年別 降水量 [mm] (rain gauge)\n多い → 表層 SWC 高", "Year", "Rain [mm]")

    plt.subplots_adjust(top=0.94)
    fp = save_dir / "fig05_tau_diagnostic.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 10. fig06 — τ × III 散布図
# ================================================================

def plot_tau_iii_scatter(tau_df: pd.DataFrame, iii_df: pd.DataFrame, save_dir: Path):
    """
    τ と III_mean の散布図。
    相関が強ければ「τ は灌漑スケジュールに規定される」仮説が支持される。
    """
    merged = tau_df.merge(iii_df, on="year", how="inner").dropna(subset=["tau","iii_mean"])
    if len(merged) < 3:
        print("  [skip] tau-III scatter — not enough data"); return

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.patch.set_facecolor(FIG_BG)
    fig.suptitle(
        "fig06 — τ × 灌漑間隔 (III) 散布図\n"
        "「τ は灌漑スケジュールに規定される」仮説の検証",
        fontsize=12, fontweight="bold")

    for ax, (x_col, x_label) in zip(axes, [
        ("iii_mean",   "III_mean [days] (灌漑イベント平均間隔, 60日超除外)"),
        ("iii_median", "III_median [days]"),
    ]):
        x = merged[x_col].values
        y = merged["tau"].values
        colors = ["#1D9E75" if r.tau_reliable else "#BDBDBD"
                  for r in merged.itertuples()]

        for i, row in merged.iterrows():
            col = "#1D9E75" if row["tau_reliable"] else "#BDBDBD"
            ax.scatter(row[x_col], row["tau"], c=col, s=100,
                       edgecolors="black", lw=1.2, zorder=5)
            ax.annotate(str(int(row["year"])),
                        (row[x_col], row["tau"]),
                        xytext=(5, 5), textcoords="offset points",
                        fontsize=9, fontweight="bold")

        # 信頼年のみで回帰直線
        trusted = merged[merged["tau_reliable"]]
        if len(trusted) >= 3:
            xv = trusted[x_col].values
            yv = trusted["tau"].values
            slope, intercept, r, p_val, _ = stats.linregress(xv, yv)
            x_line = np.linspace(xv.min()-0.5, xv.max()+0.5, 100)
            ax.plot(x_line, intercept + slope * x_line,
                    "r--", lw=1.8, alpha=0.8,
                    label=f"信頼年回帰  r={r:.2f}, p={p_val:.2f}")
            ax.legend(fontsize=9)
            corr_txt = (f"r = {r:.2f}\np = {p_val:.3f}\n"
                        f"{'→ 正の相関: τ∝III → 灌漑間隔が支配的' if r > 0.5 else '→ 相関弱い'}")
            ax.text(0.97, 0.05, corr_txt,
                    transform=ax.transAxes, va="bottom", ha="right",
                    fontsize=9, color="red" if r > 0.5 else "gray",
                    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="red" if r>0.5 else "gray", alpha=0.85))

        ax.set_ylim(bottom=0)
        _ax_style(ax,
                  f"τ vs {x_col.replace('_',' ').upper()}",
                  x_label, "τ [days]")
        leg = [mpatches.Patch(color="#1D9E75", alpha=0.85, label="信頼年 (R²≥0.7)"),
               mpatches.Patch(color="#BDBDBD", alpha=0.85, label="低信頼年")]
        ax.legend(handles=leg + ax.get_legend_handles_labels()[0], fontsize=8.5)

    plt.tight_layout()
    fp = save_dir / "fig06_tau_iii_scatter.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 11. fig07 — 3種降水データ + rain gauge 比較
# ================================================================

def plot_precip_compare(precip_df: pd.DataFrame, save_dir: Path):
    """
    rain gauge / CHIRPS / TerraClimate / ERA5 の年別降水量を比較する。

    目的:
      (1) rain gauge の空間代表性を評価
      (2) データソース間のバイアスを把握
      (3) 2021年の高降水量が異常値か否かを複数データで確認
    """
    years = precip_df["year"].values
    x     = np.arange(len(years))
    w     = 0.2

    cols = {
        "gauge_rain_mm": ("#2196F3", "Rain gauge (現地実測)"),
        "chirps_mm":     ("#4CAF50", "CHIRPS v2 (衛星+gauge融合)"),
        "terra_mm":      ("#FF9800", "TerraClimate (月次, ~4km)"),
        "era5_mm":       ("#9C27B0", "ERA5 (hourly, ~31km)"),
    }

    fig, ax = plt.subplots(figsize=(12, 6))
    fig.patch.set_facecolor(FIG_BG)

    for j, (col, (color, label)) in enumerate(cols.items()):
        if col not in precip_df.columns: continue
        vals = precip_df[col].values
        offset = (j - 1.5) * w
        bars = ax.bar(x + offset, vals, width=w, color=color, alpha=0.80,
                      edgecolor="black", lw=0.8, label=label)
        for bar, v in zip(bars, vals):
            if not np.isnan(v):
                ax.text(bar.get_x() + bar.get_width()/2, v + 3,
                        f"{v:.0f}", ha="center", va="bottom",
                        fontsize=7.5, color=color)

    ax.set_xticks(x); ax.set_xticklabels(years)
    ax.legend(fontsize=9, loc="upper right")

    # 2021 年を強調
    idx_2021 = list(years).index(2021) if 2021 in years else None
    if idx_2021 is not None:
        ax.axvspan(idx_2021 - 0.5, idx_2021 + 0.5,
                   color="#FFF9C4", alpha=0.5, zorder=0)
        ax.text(idx_2021, ax.get_ylim()[1] * 0.97,
                "2021\n(τ短い年)", ha="center", va="top",
                fontsize=8.5, color="#F57F17", fontweight="bold")

    _ax_style(ax,
              "fig07 — 年別降水量: 4データソース比較 (Tarazona)\n"
              "rain gauge が ground truth。他3つは空間代表性と補完性の評価",
              "Year", "年間降水量 [mm]")
    plt.tight_layout()
    fp = save_dir / "fig07_precip_compare.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 12. fig08 — τ × III × 降水量 統合診断
# ================================================================

def plot_tau_drivers(tau_df, iii_df, precip_df, save_dir):
    """
    τ の年別変動を III・降水量・灌漑量の3軸で同時に診断する。

    レイアウト:
      上段: τ (bar) + III_mean (折れ線, 右軸)
      下段: 降水量 (gauge + CHIRPS) + 灌漑総量 (右軸)

    読み方:
      τ と III が同じ年に高/低なら「τ ∝ III」支持
      τ が降水量とも連動するなら「大気水分も影響」
    """
    merged = tau_df.merge(iii_df, on="year", how="left")
    merged = merged.merge(precip_df, on="year", how="left")
    years  = merged["year"].values
    x      = np.arange(len(years))
    w      = 0.5

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    fig.patch.set_facecolor(FIG_BG)
    fig.suptitle(
        "fig08 — τ の年別変動ドライバー統合診断\n"
        "τ × III × 降水量 × 灌漑量を1枚で俯瞰",
        fontsize=12, fontweight="bold")

    # ── 上段: τ (bar) + III_mean (折れ線) ──
    colors_tau = ["#1D9E75" if r.tau_reliable else "#BDBDBD"
                  for r in merged.itertuples()]
    ax1.bar(x, merged["tau"], color=colors_tau, alpha=0.80,
            edgecolor="black", lw=1, width=w, label="τ [days]")
    for xi, row in zip(x, merged.itertuples()):
        if not np.isnan(row.tau):
            ax1.text(xi, row.tau + 0.4, f"{row.tau:.1f}d",
                     ha="center", va="bottom", fontsize=8.5,
                     color="black" if row.tau_reliable else "#999")

    ax1r = ax1.twinx()
    iii_vals = merged["iii_mean"].values
    ax1r.plot(x, iii_vals, "o-", color="#E85D04", lw=2, ms=8,
              label="III_mean [days]")
    for xi, v in zip(x, iii_vals):
        if not np.isnan(v):
            ax1r.text(xi, v + 0.2, f"{v:.1f}d",
                      ha="center", va="bottom", fontsize=8, color="#E85D04")

    ax1.set_ylabel("τ [days]", fontsize=9)
    ax1r.set_ylabel("III_mean [days]", fontsize=9, color="#E85D04")
    ax1r.tick_params(axis="y", colors="#E85D04")
    ax1.set_ylim(bottom=0)
    ax1r.set_ylim(bottom=0)

    # 凡例を統合
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax1r.get_legend_handles_labels()
    leg_p = [mpatches.Patch(color="#1D9E75", alpha=0.85, label="τ 信頼年"),
             mpatches.Patch(color="#BDBDBD", alpha=0.85, label="τ 低信頼年")]
    ax1.legend(handles=h1+h2+leg_p, fontsize=8.5, loc="upper left")
    _ax_style(ax1,
              "(A) τ (bar) と灌漑間隔 III_mean (折れ線)\n"
              "τ と III が連動 → 灌漑スケジュールが τ を規定",
              "", "τ [days]")

    # ── 下段: 降水量 + 灌漑総量 ──
    ax2.bar(x - 0.15, merged.get("gauge_rain_mm", np.nan), width=0.3,
            color="#2196F3", alpha=0.75, edgecolor="black", lw=0.8,
            label="Rain gauge [mm]")
    if "chirps_mm" in merged.columns:
        ax2.bar(x + 0.15, merged["chirps_mm"], width=0.3,
                color="#4CAF50", alpha=0.75, edgecolor="black", lw=0.8,
                label="CHIRPS [mm]")

    ax2r = ax2.twinx()
    if "total_irrig_mm" in merged.columns:
        ax2r.plot(x, merged["total_irrig_mm"], "s--", color="#FF9800",
                  lw=2, ms=8, label="灌漑総量 [mm]")
        ax2r.set_ylabel("灌漑総量 [mm]", fontsize=9, color="#FF9800")
        ax2r.tick_params(axis="y", colors="#FF9800")

    ax2.set_xticks(x); ax2.set_xticklabels(years)
    ax2.set_ylabel("降水量 [mm]", fontsize=9)
    ax2.set_ylim(bottom=0)
    h3, l3 = ax2.get_legend_handles_labels()
    h4, l4 = ax2r.get_legend_handles_labels()
    ax2.legend(handles=h3+h4, fontsize=8.5, loc="upper left")
    _ax_style(ax2,
              "(B) 降水量 (gauge + CHIRPS) と灌漑総量\n"
              "2021年の高降水量が複数データで確認されるか確認",
              "Year", "降水量 [mm]")

    plt.tight_layout(h_pad=2)
    fp = save_dir / "fig08_tau_drivers.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 10. fig01–fig05 (v17 と同じ — 省略なし)
# ================================================================

def plot_verdict_panel(pkg_df, save_dir):
    le_df = pkg_df[pkg_df["var"]=="LE_corr"].sort_values(["site","season"]).reset_index(drop=True)
    n_rows = len(le_df)
    if n_rows == 0: return
    fig, axes = plt.subplots(1, 3, figsize=(22, max(5, n_rows*1.1)))
    fig.patch.set_facecolor(FIG_BG)
    fig.suptitle("fig01 — 三点パッケージ統合判定 (LE_corr)\n"
                 "判定ガード: n_n,n_s≥30 かつ 95%CI幅≤0.50",
                 fontsize=13, fontweight="bold", y=1.01)
    y_pos  = np.arange(n_rows)
    labels = [f"{r.site}\n{r.season}" for r in le_df.itertuples()]
    ax = axes[0]
    for i, r in enumerate(le_df.itertuples()):
        pkg = {"sds":r.sds,"sds_lo":r.sds_lo,"sds_hi":r.sds_hi,"n_n":r.n_n,"n_s":r.n_s}
        v, strength = verdict(pkg)
        col = VERDICT_COL[v]; alpha = 0.90 if strength=="strong" else 0.45
        if not np.isnan(r.sds):
            ax.barh(i, r.sds, color=col, alpha=alpha, edgecolor="black", lw=0.8, height=0.65)
        lo_ = getattr(r,"sds_lo",np.nan); hi_ = getattr(r,"sds_hi",np.nan)
        if not (np.isnan(lo_) or np.isnan(hi_)):
            ax.errorbar(r.sds, i, xerr=[[r.sds-lo_],[hi_-r.sds]],
                        color="black", capsize=5, lw=1.2, fmt="none")
        ax.text(-0.55, i, f"n=({int(r.n_n)},{int(r.n_s)})",
                va="center", ha="left", fontsize=7.5, color="#555")
        if not np.isnan(r.sds):
            ax.text(r.sds+0.03, i, f"{r.sds:+.2f}  [{VERDICT_JP[v]}]",
                    va="center", fontsize=8.5,
                    fontweight="bold" if strength=="strong" else "normal")
    ax.axvline(0, color="gray", lw=1.0)
    ax.axvspan(-DEEP_ROOT_BAND, DEEP_ROOT_BAND, color="#1D9E75", alpha=0.07)
    ax.set_xlim(-0.6, 1.0); ax.set_yticks(y_pos); ax.set_yticklabels(labels, fontsize=9)
    patches = [mpatches.Patch(color=c, label=VERDICT_JP[k], alpha=0.8)
               for k,c in VERDICT_COL.items()]
    ax.legend(handles=patches, loc="lower right", fontsize=7, title="判定", title_fontsize=8)
    _ax_style(ax, "(1) SDS\nSDS≈0→深根的  SDS>0→浅根的", "SDS")
    ax = axes[1]
    for i, r in enumerate(le_df.itertuples()):
        if not np.isnan(r.abs_diff):
            ax.barh(i, r.abs_diff, color="#5B4FCF", alpha=0.75,
                    edgecolor="black", lw=0.8, height=0.65)
            if not (np.isnan(r.abs_lo) or np.isnan(r.abs_hi)):
                ax.errorbar(r.abs_diff, i,
                            xerr=[[r.abs_diff-r.abs_lo],[r.abs_hi-r.abs_diff]],
                            color="black", capsize=5, lw=1.2, fmt="none")
            ax.text(r.abs_diff+(1.5 if r.abs_diff>=0 else -1.5), i,
                    f"{r.abs_diff:+.1f} W/m²", va="center", fontsize=8.5)
    ax.axvline(0, color="gray", lw=1.0); ax.set_yticks(y_pos); ax.set_yticklabels([])
    _ax_style(ax, "(2) 絶対差 = LE_n − LE_s", "LE_n − LE_s [W/m²]")
    ax = axes[2]
    for i, r in enumerate(le_df.itertuples()):
        if not np.isnan(r.rb):
            col = "#E85D04" if abs(r.rb)>0.3 else "#FFA000" if abs(r.rb)>0.1 else "#9E9E9E"
            ax.barh(i, r.rb, color=col, alpha=0.80, edgecolor="black", lw=0.8, height=0.65)
            ax.text(r.rb+0.03, i, f"{r.rb:+.2f}  p={r.p:.1e} {_sig_star(r.p)}",
                    va="center", fontsize=8.5)
    ax.axvline(0, color="gray", lw=1.0); ax.set_xlim(-1, 1.5)
    ax.set_yticks(y_pos); ax.set_yticklabels([])
    _ax_style(ax, "(3) rank-biserial r + 有意性", "rank-biserial r")
    plt.tight_layout()
    plt.savefig(save_dir/"fig01_verdict_panel.png", dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close(); print(f"  [save] {save_dir}/fig01_verdict_panel.png")


def plot_recovery_curve(fit_res, tau_by_year, save_dir):
    grouped = fit_res.get("grouped")
    if grouped is None or len(grouped) < 2: return
    tau, le0, le_inf = fit_res["tau"], fit_res["le0"], fit_res["le_inf"]
    r2_exp, aic_e, aic_l, aic_g = fit_res["r2_exp"], fit_res["aic_exp"], fit_res["aic_lin"], fit_res["aic_log"]
    fig, axes = plt.subplots(2, 2, figsize=(16, 12)); fig.patch.set_facecolor(FIG_BG)
    x_obs = grouped["days_since_irrig"].values.astype(float)
    y_obs = grouped["median"].values; n_obs = grouped["count"].values
    x_fit = np.linspace(0, x_obs.max()+2, 300); sizes = 40 + n_obs*4
    ax = axes[0,0]
    ax.scatter(x_obs, y_obs, s=sizes, c="#1D9E75", edgecolors="black", lw=1, zorder=5, label="観測中央値")
    if not np.isnan(tau):
        ax.plot(x_fit, exp_model(x_fit, le_inf, le0, tau), "k-", lw=2.5,
                label=f"Exp fit: τ={tau:.1f}d, R²={r2_exp:.3f}")
        ax.axhline(le_inf, color="#E85D04", ls="--", lw=1.5,
                   label=f"LE_∞={le_inf:.0f} W/m²")
        ax.axhline(le0, color="#1A73E8", ls="--", lw=1.5, label=f"LE_0={le0:.0f} W/m²")
        ax.axvline(tau, color="gray", ls=":", lw=1.5, label=f"τ={tau:.1f}d")
        ic = "#E85D04" if tau<5 else "#FFA000" if tau<=14 else "#1D9E75"
        it = (f"τ={tau:.1f}d<5d → 灌漑依存" if tau<5 else
              f"τ={tau:.1f}d>14d → 深層水?" if tau>14 else f"τ={tau:.1f}d → 中間")
        ax.text(0.97, 0.97, it, transform=ax.transAxes, va="top", ha="right",
                fontsize=10, fontweight="bold", color=ic,
                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=ic, alpha=0.9))
    ax.set_ylim(bottom=0); ax.legend(fontsize=8.5)
    _ax_style(ax, f"(A) Recovery curve — τ={tau:.1f}d",
              "Days since last irrigation", "LE_corr 中央値 [W/m²]")
    ax = axes[0,1]
    ax.scatter(x_obs, y_obs, s=sizes, c="#1D9E75", edgecolors="black", lw=1, zorder=5, label="観測中央値")
    best_aic = min(v for v in [aic_e,aic_l,aic_g] if not np.isnan(v)) if any(
        not np.isnan(v) for v in [aic_e,aic_l,aic_g]) else np.nan
    for name, fn, popt, r2, aic, ls, lw_ in [
        ("Exponential", exp_model, fit_res.get("popt_exp"), r2_exp,             aic_e, "k-",  2.5),
        ("Linear",      lin_model, fit_res.get("popt_lin"), fit_res["r2_lin"],   aic_l, "b--", 1.8),
        ("Logarithmic", log_model, fit_res.get("popt_log"), fit_res["r2_log"],   aic_g, "r-.", 1.8),
    ]:
        if popt is not None and not np.isnan(popt[0]):
            da = aic-best_aic if not np.isnan(aic) and not np.isnan(best_aic) else np.nan
            bm = "  ← BEST" if not np.isnan(da) and abs(da)<0.1 else ""
            lbl = f"{name}  R²={r2:.3f}, ΔAIC={da:.1f}{bm}" if not np.isnan(da) else f"{name}  R²={r2:.3f}"
            ax.plot(x_fit, fn(x_fit, *popt), ls, lw=lw_, label=lbl)
    ax.set_ylim(bottom=0); ax.legend(fontsize=8.5)
    _ax_style(ax, "(B) モデル比較", "Days since last irrigation", "LE_corr 中央値 [W/m²]")
    ax = axes[1,0]
    valid = tau_by_year.dropna(subset=["tau"])
    if len(valid) >= 2:
        colors_bar = ["#1D9E75" if r.tau_reliable else "#BDBDBD" for r in valid.itertuples()]
        ax.bar(valid["year"], valid["tau"], color=colors_bar, alpha=0.85,
               edgecolor="black", lw=1, width=0.6)
        ax.axhline(tau, color="#E85D04", ls="--", lw=2, label=f"全期間 τ={tau:.1f}d")
        ax.axhline(5, color="gray", ls=":", lw=1, alpha=0.7, label="τ=5d 閾値")
        trusted = valid[valid["tau_reliable"]]
        if len(trusted)>=2:
            ax.axhspan(trusted["tau"].min(), trusted["tau"].max(), color="#1D9E75", alpha=0.08,
                       label=f"信頼年: {trusted['tau'].min():.1f}–{trusted['tau'].max():.1f}d")
        for _, row in valid.iterrows():
            ax.text(row["year"], row["tau"]+0.3, f"{row['tau']:.1f}d",
                    ha="center", va="bottom", fontsize=8,
                    color="black" if row["tau_reliable"] else "#999")
        leg_p = [mpatches.Patch(color="#1D9E75", alpha=0.85, label="信頼年 (R²≥0.7)"),
                 mpatches.Patch(color="#BDBDBD", alpha=0.85, label="低信頼年")]
        h, _ = ax.get_legend_handles_labels()
        ax.legend(handles=h+leg_p, fontsize=9)
    _ax_style(ax, "(C) 年別 τ — 信頼性フラグ付き", "Year", "τ [days]")
    ax = axes[1,1]; ax.set_xlim(0,10); ax.set_ylim(0,10); ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(mpatches.Rectangle((0,0),10,10,fc=FIG_BG,ec="none"))
    ax.text(5,9.5,"(D) 二層モデル解釈",ha="center",va="top",fontsize=11,fontweight="bold")
    def box(x,y,w,h,fc,text="",fs=9):
        ax.add_patch(mpatches.FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.15",
                                             fc=fc,ec="black",lw=1.5,alpha=0.85))
        if text: ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=fs,fontweight="bold")
    def arrow(x1,y1,x2,y2,col="black"):
        ax.annotate("",xy=(x2,y2),xytext=(x1,y1),
                    arrowprops=dict(arrowstyle="-|>",color=col,lw=1.8))
    box(0.5,7.5,2.5,1.2,fc="#AED6F1",text="灌漑パルス\n(drip)")
    arrow(3.0,8.1,4.5,8.1)
    box(4.5,7.2,4.5,1.8,fc="#D5EAF5",text=f"Fast pool\n(灌漑水→τ≈{tau:.1f}d)")
    arrow(6.75,7.2,6.75,5.8,col="#1A73E8")
    box(4.5,4.2,4.5,1.5,fc="#D5F5E3",text=f"Slow pool\n(深層水?→LE_∞≈{le_inf:.0f}W/m²)")
    arrow(6.75,4.2,6.75,2.8,col="#1D9E75")
    box(4.0,1.5,5.5,1.2,fc="#FDE8D8",text=f"蒸散 LE: {le0:.0f}→{le_inf:.0f} W/m²")
    ax.text(5,0.5,"SWC-ET decoupling は灌漑水が主要因\n(深根を棄却せず、但し支配的でない)",
            ha="center",va="bottom",fontsize=8.5,style="italic",color="#333",
            bbox=dict(boxstyle="round,pad=0.3",fc="white",ec="#999",alpha=0.8))
    plt.tight_layout(h_pad=3,w_pad=3)
    plt.savefig(save_dir/"fig02_recovery_curve.png", dpi=150, bbox_inches="tight", facecolor=FIG_BG)
    plt.close(); print(f"  [save] {save_dir}/fig02_recovery_curve.png")


# ================================================================
# 11. MAIN
# ================================================================

def main():
    print("="*60)
    print("解析A v18 — III + 3種降水データ + τ Diagnostic 拡張")
    print("="*60)

    df = load_and_merge()
    df = add_days_since_irrig(df)
    df = add_irrig_active_month(df)

    # ── 三点パッケージ ──
    print(f"\n{'='*60}\n三点パッケージ計算\n{'='*60}")
    rows = []
    for season_label, months in SEASONS.items():
        for site in ["Oran","Tarazona"]:
            sub = df[(df["site"]==site) & df["is_growing"] & df["month"].isin(months)]
            if site=="Tarazona": sub = sub[sub["irrig_active_month"]]
            for var in ["LE_corr","EF_corr","ET"]:
                if var not in sub.columns: continue
                n_arr = sub.loc[sub["drought_type"]=="normal",   var].dropna().values
                s_arr = sub.loc[sub["drought_type"]=="soil dry", var].dropna().values
                pkg   = three_point_package(n_arr, s_arr, var)
                v, strength = verdict(pkg)
                pkg.update(dict(site=site, season=season_label, var=var,
                                verdict=v, verdict_strength=strength))
                rows.append(pkg)
                if var=="LE_corr":
                    print(f"  [{site:9s}/{season_label:22s}] "
                          f"SDS={pkg['sds']:+.3f} abs={pkg['abs_diff']:+.1f}W/m² "
                          f"rb={pkg['rb']:+.3f} p={pkg['p']:.1e} "
                          f"n=({pkg['n_n']},{pkg['n_s']}) → [{VERDICT_JP[v]}] ({strength})")
    pkg_df = pd.DataFrame(rows)
    pkg_df.to_csv(OUT_DIR/"v18_three_point_package.csv", index=False)

    # ── Recovery τ ──
    print(f"\n{'='*60}\nRecovery τ\n{'='*60}")
    tara_active = df[(df["site"]=="Tarazona") & df["irrig_active_month"] & df["is_growing"]]
    fit_res     = fit_recovery(tara_active)
    tau_by_year = fit_tau_by_year(tara_active)
    irrig_stats = calc_irrig_stats_by_year(df)
    print(f"  τ(全期間)={fit_res['tau']:.2f}d  LE_0={fit_res['le0']:.1f}  LE_∞={fit_res['le_inf']:.1f} W/m²")
    print(f"  R² exp={fit_res['r2_exp']:.3f} lin={fit_res['r2_lin']:.3f} log={fit_res['r2_log']:.3f}")
    print(f"  AIC exp={fit_res['aic_exp']:.1f} lin={fit_res['aic_lin']:.1f} log={fit_res['aic_log']:.1f}")
    print(tau_by_year.to_string(index=False))

    # ── III ──
    print(f"\n{'='*60}\nInter-Irrigation Interval (III)\n{'='*60}")
    iii_df = calc_iii_by_year(df)
    print(iii_df.to_string(index=False))
    merged_check = tau_by_year.merge(iii_df, on="year").dropna(subset=["tau","iii_mean"])
    if len(merged_check) >= 3:
        trusted = merged_check[merged_check["tau_reliable"]]
        if len(trusted) >= 3:
            r, p = stats.pearsonr(trusted["tau"], trusted["iii_mean"])
            print(f"\n  τ × III_mean 相関 (信頼年): r={r:.3f}, p={p:.3f}")
            if r > 0.5:
                print(f"  → 正の相関 → 「τ は灌漑間隔に規定される」仮説を支持")
            else:
                print(f"  → 相関が弱い → τ に他の要因も関与")

    # ── 3種降水データ抽出 ──
    print(f"\n{'='*60}\n降水データ抽出\n{'='*60}")
    chirps = extract_chirps()
    terra  = extract_terraclimate()
    era5   = extract_era5()
    precip_df = (irrig_stats[["year","gauge_rain_mm"]]
                 .merge(chirps, on="year", how="left")
                 .merge(terra,  on="year", how="left")
                 .merge(era5,   on="year", how="left"))
    print(precip_df.to_string(index=False))
    precip_df.to_csv(OUT_DIR/"v18_precip_comparison.csv", index=False)
    iii_df.to_csv(OUT_DIR/"v18_iii_by_year.csv", index=False)

    # ── 可視化 ──
    print(f"\n--- 可視化 ---")
    plot_verdict_panel(pkg_df, OUT_DIR)
    plot_recovery_curve(fit_res, tau_by_year, OUT_DIR)
    plot_three_point_matrix(pkg_df, OUT_DIR)
    plot_tau_sensitivity(tau_by_year, fit_res["tau"], OUT_DIR)
    plot_tau_diagnostic(tau_by_year, irrig_stats, OUT_DIR)
    plot_tau_iii_scatter(tau_by_year, iii_df, OUT_DIR)
    plot_precip_compare(precip_df, OUT_DIR)
    drivers_df = precip_df.merge(
        irrig_stats[["year", "total_irrig_mm", "n_irrig_events"]], on="year", how="left")
    plot_tau_drivers(tau_by_year, iii_df, drivers_df, OUT_DIR)

    # ── 最終サマリー ──
    print(f"\n{'='*60}\n★ 最終サマリー\n{'='*60}")
    for _, row in pkg_df[pkg_df["var"]=="LE_corr"].iterrows():
        v, strength = row["verdict"], row["verdict_strength"]
        lo, hi = row.get("sds_lo",np.nan), row.get("sds_hi",np.nan)
        ci = f"[{lo:+.2f},{hi:+.2f}]" if not np.isnan(lo) else "[CI算出不可]"
        print(f"  [{row['site']:9s}/{row['season']:22s}] "
              f"SDS={row['sds']:+.3f} {ci} n=({row['n_n']:.0f},{row['n_s']:.0f}) "
              f"→ {VERDICT_JP[v]} ({strength})")
    tau = fit_res["tau"]
    print(f"\n  Recovery τ={tau:.1f}d  主張: 'Irrigation-driven decoupling (τ≈{tau:.0f}d)'")
    print(f"\n[done] → {OUT_DIR}/")


if __name__ == "__main__":
    main()
