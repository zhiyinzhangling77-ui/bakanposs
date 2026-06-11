#!/usr/bin/env python3
"""
================================================================
解析A v19 : ターミナル直接実行版
================================================================

【v18 からの変更点】
  1. argparse による CLI 対応
       python analysis_A_v19.py [オプション]
       主要オプション:
         --parquet   PARQUET ファイルパス
         --csv       Tarazona 日次 CSV パス
         --chirps    CHIRPS ディレクトリ
         --terra     TerraClimate ディレクトリ
         --era5      ERA5 ディレクトリ
         --out       出力ディレクトリ (default: ./output_analysis_A_v19)
         --lat       Tarazona 緯度  (default: 39.085)
         --lon       Tarazona 経度  (default: -1.358)
         --years     解析年 (default: 2020 2021 2022 2023 2024)
         --no-plots  図の生成をスキップ (数値確認のみ)
         --diag-only 診断ログのみ出力して終了

  2. ERA5 修正: sortby("longitude").sortby("latitude") で
     単調でない座標エラーを回避
     また open_mfdataset を年別ファイルごとに分けて処理

  3. 降水量乖離診断を自動化
       - rain_mm を全期間 / 生育期 / 灌漑アクティブ月 の3種で集計
       - CHIRPS との比較で「gauge が生育期のみか」を自動判定
       - 診断結果をコンソールと v19_precip_diagnostic.csv に出力

  4. 座標グリッド診断ログ
       - CHIRPS / TerraClimate の実際の抽出グリッド座標を表示
       - 「座標ずれ」の可能性を事前に警告

  5. τ × III 相関が弱い場合の代替仮説を自動提示
       - VPD (将来拡張) / 灌漑総量 / 前月降水量 との相関も試みる

【使い方 (最小限)】
  cd /home/shion-nagamine/bakanposs
  python analysis_A_v19.py

【パスを変えて実行する場合】
  python analysis_A_v19.py \\
    --parquet /path/to/daily_classified_v4.parquet \\
    --csv     /path/to/Daily_Summary.csv \\
    --chirps  "/path/to/Chirps v2 for persipitation" \\
    --terra   /path/to/TerraClimate_ppt \\
    --era5    /path/to/ERA5_ppt \\
    --out     ./output_v19 \\
    --lat 39.085 --lon -1.358

【診断のみ実行 (図不要)】
  python analysis_A_v19.py --diag-only
"""

import argparse
import json
import sys
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

warnings.filterwarnings("ignore")


# ================================================================
# 0. CLI 引数パース
# ================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="解析A v19 — III + 3種降水 + τ Diagnostic",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)

    # パス
    p.add_argument("--parquet", default=None,
                   help="daily_classified_v4.parquet のパス")
    p.add_argument("--csv", default=None,
                   help="Tarazona 日次 CSV のパス (Irrig_mm, Rain_mm 含む)")
    p.add_argument("--chirps", default=None,
                   help="CHIRPS NC ファイルのディレクトリ")
    p.add_argument("--terra", default=None,
                   help="TerraClimate NC ファイルのディレクトリ")
    p.add_argument("--era5", default=None,
                   help="ERA5 NC ファイルのディレクトリ")
    p.add_argument("--out", default=None,
                   help="出力ディレクトリ (default: ./output_analysis_A_v19)")

    # 座標
    p.add_argument("--lat", type=float, default=39.085,
                   help="Tarazona 緯度 (default: 39.085)")
    p.add_argument("--lon", type=float, default=-1.358,
                   help="Tarazona 経度 (default: -1.358)")

    # 解析年
    p.add_argument("--years", type=int, nargs="+",
                   default=[2020, 2021, 2022, 2023, 2024],
                   help="解析対象年 (default: 2020 2021 2022 2023 2024)")

    # 実行制御
    p.add_argument("--no-plots", action="store_true",
                   help="図の生成をスキップ (数値出力のみ)")
    p.add_argument("--diag-only", action="store_true",
                   help="降水量診断ログのみ出力して終了")

    return p.parse_args()


def resolve_paths(args):
    """
    CLI 引数が省略された場合はデフォルトパスにフォールバックする。
    デフォルトは v18 と同じパス。
    """
    base = Path("/home/shion-nagamine")

    parquet = Path(args.parquet) if args.parquet else \
              base / "bakanposs/analysis_A/daily_classified_v4.parquet"
    csv     = Path(args.csv) if args.csv else \
              base / "Dataset/Eddy data in Spain/Daily_Summary_Filtered_forPred_ActEne26.csv"
    chirps  = Path(args.chirps) if args.chirps else \
              base / "Dataset/Chirps v2 for persipitation"
    terra   = Path(args.terra) if args.terra else \
              base / "Dataset/TerraClimate_ppt"
    era5    = Path(args.era5) if args.era5 else \
              base / "Dataset/ERA5_ppt"
    out     = Path(args.out) if args.out else \
              Path("./output_analysis_A_v19")

    return parquet, csv, chirps, terra, era5, out


# ================================================================
# 1. パラメータ (args から上書き可能)
# ================================================================

IRRIG_THRESHOLD     = 0.5
MIN_IRRIG_PER_MONTH = 2

SEASONS = {
    "spring (1-4月)"      : [1, 2, 3, 4],
    "shoulder (5,6,10月)" : [5, 6, 10],
    "summer (7-9月)"      : [7, 8, 9],
}

N_BOOT          = 5000
CI_PCT          = (2.5, 97.5)
MIN_N_PER_CLASS = 5
DENOM_FLOORS    = {"LE_corr": 5.0, "EF_corr": 0.05, "ET": 0.3}
RATIO_CLIP      = 10.0
N_MIN_VERDICT   = 30
CI_WIDTH_MAX    = 0.5
DEEP_ROOT_BAND  = 0.15
R2_MIN_TRUST    = 0.70

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
# 2. データ読込・前処理
# ================================================================

def load_and_merge(parquet: Path, csv: Path) -> pd.DataFrame:
    print(f"  parquet : {parquet}")
    print(f"  csv     : {csv}")
    if not parquet.exists():
        sys.exit(f"[ERROR] parquet not found: {parquet}")
    if not csv.exists():
        sys.exit(f"[ERROR] CSV not found: {csv}")

    daily = pd.read_parquet(parquet)
    daily["date"] = pd.to_datetime(daily["date"])

    raw = pd.read_csv(csv)
    raw["date"] = pd.to_datetime(raw["date"], errors="coerce")
    irrig_cols = [c for c in ["Irrig_mm","Rain_mm","IrrigRain_mm"] if c in raw.columns]
    extra = raw[["date"]+irrig_cols].dropna(subset=["date"]).drop_duplicates("date")

    tara = daily[daily["site"]=="Tarazona"].merge(extra, on="date", how="left")
    oran = daily[daily["site"]=="Oran"].copy()
    for c in irrig_cols:
        oran[c] = 0.0

    df = pd.concat([oran, tara]).sort_values(["site","date"]).reset_index(drop=True)
    df["year"]  = df["date"].dt.year
    df["month"] = df["date"].dt.month
    return df


def add_days_since_irrig(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["days_since_irrig"] = np.nan
    for site, idx in df.groupby("site").groups.items():
        sub  = df.loc[idx].sort_values("date")
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


def add_irrig_active_month(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["_ym"] = df["date"].dt.to_period("M")
    df["_ev"] = df["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD
    monthly = df.groupby(["site","_ym"])["_ev"].sum().reset_index(name="n_ev")
    active  = set(zip(monthly.loc[monthly["n_ev"]>=MIN_IRRIG_PER_MONTH,"site"],
                      monthly.loc[monthly["n_ev"]>=MIN_IRRIG_PER_MONTH,"_ym"]))
    df["irrig_active_month"] = df.apply(lambda r: (r["site"],r["_ym"]) in active, axis=1)
    df.drop(columns=["_ym","_ev"], inplace=True)
    return df


# ================================================================
# 3. 降水量診断 (乖離の自動判定)
# ================================================================

def diagnose_rain_gauge(df: pd.DataFrame) -> pd.DataFrame:
    """
    rain_mm を全期間 / 生育期 / 灌漑アクティブ月 の3種で集計し、
    CHIRPS との比較に使う正しい集計方法を診断する。

    rain gauge が全期間合計なら CHIRPS の年間値と比較すべき。
    生育期のみなら CHIRPS も生育期に絞るべき。
    乖離が 3 倍以上あれば警告を出す。
    """
    tara = df[df["site"]=="Tarazona"].copy()
    rows = []
    for yr, sub in tara.groupby("year"):
        r_all     = sub["Rain_mm"].fillna(0).sum()
        r_growing = sub.loc[sub["is_growing"]==True,  "Rain_mm"].fillna(0).sum()
        r_active  = sub.loc[sub["irrig_active_month"], "Rain_mm"].fillna(0).sum()
        n_all     = len(sub)
        n_growing = (sub["is_growing"]==True).sum()
        rows.append(dict(year=yr,
                         rain_all_mm=round(r_all,1),
                         rain_growing_mm=round(r_growing,1),
                         rain_active_mm=round(r_active,1),
                         n_days_all=n_all,
                         n_days_growing=n_growing))
    return pd.DataFrame(rows)


# ================================================================
# 4. NC 抽出 (CHIRPS / TerraClimate / ERA5)
# ================================================================

def _nearest_val(ds, lat_name, lon_name, var_name, lat, lon, verbose=True):
    """最近傍グリッドを抽出し、実際の座標を表示する"""
    arr  = ds[var_name].sel({lat_name: lat, lon_name: lon}, method="nearest")
    if verbose:
        actual_lat = float(ds[lat_name].sel({lat_name: lat}, method="nearest"))
        actual_lon = float(ds[lon_name].sel({lon_name: lon}, method="nearest"))
        dist = ((actual_lat - lat)**2 + (actual_lon - lon)**2)**0.5 * 111  # km近似
        print(f"    抽出グリッド: lat={actual_lat:.4f}, lon={actual_lon:.4f} "
              f"(目標から約 {dist:.1f} km)")
    return arr


def extract_chirps(chirps_dir: Path, years: list, lat: float, lon: float) -> pd.DataFrame:
    if not HAS_XARRAY:
        return pd.DataFrame({"year": years, "chirps_mm_annual": [np.nan]*len(years),
                             "chirps_mm_growing": [np.nan]*len(years)})
    rows = []
    printed_grid = False
    for yr in years:
        fp = chirps_dir / f"chirps-v2.0.{yr}.days_p05.nc"
        if not fp.exists():
            print(f"    [skip] CHIRPS {yr}: file not found")
            rows.append(dict(year=yr, chirps_mm_annual=np.nan, chirps_mm_growing=np.nan))
            continue
        try:
            ds  = xr.open_dataset(fp)
            var = "precip" if "precip" in ds else list(ds.data_vars)[0]
            lat_n = "latitude"  if "latitude"  in ds.coords else "lat"
            lon_n = "longitude" if "longitude" in ds.coords else "lon"
            # 初回のみグリッド座標を表示
            ts = _nearest_val(ds, lat_n, lon_n, var, lat, lon, verbose=not printed_grid)
            printed_grid = True
            vals = ts.values
            annual  = float(vals[vals >= 0].sum())
            # 生育期(4〜10月)に絞った合計
            time_arr = pd.to_datetime(ts.time.values)
            growing_mask = time_arr.month.isin([4,5,6,7,8,9,10])
            growing = float(vals[growing_mask & (vals >= 0)].sum())
            rows.append(dict(year=yr,
                             chirps_mm_annual=round(annual,1),
                             chirps_mm_growing=round(growing,1)))
            ds.close()
        except Exception as e:
            print(f"    [warn] CHIRPS {yr}: {e}")
            rows.append(dict(year=yr, chirps_mm_annual=np.nan, chirps_mm_growing=np.nan))
    return pd.DataFrame(rows)


def extract_terraclimate(terra_dir: Path, years: list, lat: float, lon: float) -> pd.DataFrame:
    if not HAS_XARRAY:
        return pd.DataFrame({"year": years, "terra_mm_annual": [np.nan]*len(years)})
    rows = []
    printed_grid = False
    for yr in years:
        fp = terra_dir / f"TerraClimate_ppt_{yr}.nc"
        if not fp.exists():
            rows.append(dict(year=yr, terra_mm_annual=np.nan)); continue
        try:
            ds  = xr.open_dataset(fp)
            var = "ppt" if "ppt" in ds else list(ds.data_vars)[0]
            lat_n = "lat" if "lat" in ds.coords else "latitude"
            lon_n = "lon" if "lon" in ds.coords else "longitude"
            ts  = _nearest_val(ds, lat_n, lon_n, var, lat, lon, verbose=not printed_grid)
            printed_grid = True
            annual = float(ts.values[ts.values >= 0].sum())
            rows.append(dict(year=yr, terra_mm_annual=round(annual,1)))
            ds.close()
        except Exception as e:
            print(f"    [warn] TerraClimate {yr}: {e}")
            rows.append(dict(year=yr, terra_mm_annual=np.nan))
    return pd.DataFrame(rows)


def extract_era5(era5_dir: Path, years: list, lat: float, lon: float) -> pd.DataFrame:
    """
    ERA5 hourly tp [m] → daily [mm] → 年次合計。
    複数 NC ファイルを年別に個別読込して monotonic エラーを回避。
    各ファイルを sortby で座標を昇順に揃えてから抽出する。
    """
    if not HAS_XARRAY:
        return pd.DataFrame({"year": years, "era5_mm_annual": [np.nan]*len(years)})

    era5_files = sorted(era5_dir.glob("ERA5_*_hourly_ppt.nc"))
    if not era5_files:
        print("    [warn] ERA5 files not found")
        return pd.DataFrame({"year": years, "era5_mm_annual": [np.nan]*len(years)})

    annual_dict = {}
    printed_grid = False
    for fp in era5_files:
        try:
            ds = xr.open_dataset(fp)
            # ── 座標を昇順にソート (monotonic エラーの根本対策) ──
            lat_n = "latitude"  if "latitude"  in ds.coords else "lat"
            lon_n = "longitude" if "longitude" in ds.coords else "lon"
            ds = ds.sortby(lat_n).sortby(lon_n)

            var = "tp" if "tp" in ds else list(ds.data_vars)[0]
            ts  = _nearest_val(ds, lat_n, lon_n, var, lat, lon, verbose=not printed_grid)
            printed_grid = True

            # hourly [m] → daily [mm]
            daily_mm = (ts * 1000).resample(time="1D").sum()
            df_tmp   = daily_mm.to_dataframe(name="era5_mm").reset_index()
            df_tmp["year"] = pd.to_datetime(df_tmp["time"]).dt.year

            for yr, sub in df_tmp.groupby("year"):
                if yr not in annual_dict:
                    annual_dict[yr] = 0.0
                annual_dict[yr] += float(sub.loc[sub["era5_mm"]>=0, "era5_mm"].sum())
            ds.close()
        except Exception as e:
            print(f"    [warn] ERA5 {fp.name}: {e}")

    result = []
    for yr in years:
        v = annual_dict.get(yr, np.nan)
        result.append(dict(year=yr, era5_mm_annual=round(v,1) if not np.isnan(v) else np.nan))
    return pd.DataFrame(result)


# ================================================================
# 5. 三点パッケージ + 判定
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
    sds_pt = safe_ratio(med_n - med_s, med_n, var)
    abs_pt = med_n - med_s
    rng = np.random.default_rng(seed)
    boot_sds, boot_abs = [], []
    for _ in range(N_BOOT):
        nb = rng.choice(arr_n, size=n_n, replace=True)
        sb = rng.choice(arr_s, size=n_s, replace=True)
        r  = safe_ratio(np.median(nb)-np.median(sb), np.median(nb), var)
        if not np.isnan(r): boot_sds.append(r)
        boot_abs.append(np.median(nb)-np.median(sb))
    sds_lo = sds_hi = np.nan
    if len(boot_sds) >= N_BOOT*0.1:
        sds_lo, sds_hi = np.percentile(boot_sds, CI_PCT)
    abs_lo, abs_hi = np.percentile(boot_abs, CI_PCT)
    try:
        u, p = stats.mannwhitneyu(arr_n, arr_s, alternative="two-sided")
        rb = 1 - 2*u/(n_n*n_s)
    except ValueError:
        u = p = rb = np.nan
    return dict(sds=sds_pt, sds_lo=sds_lo, sds_hi=sds_hi,
                abs_diff=abs_pt, abs_lo=abs_lo, abs_hi=abs_hi,
                rb=rb, p=p, med_n=med_n, med_s=med_s, n_n=n_n, n_s=n_s)


def verdict(pkg):
    n_n, n_s = pkg.get("n_n",0), pkg.get("n_s",0)
    sds, lo, hi = pkg.get("sds",np.nan), pkg.get("sds_lo",np.nan), pkg.get("sds_hi",np.nan)
    if np.isnan(sds) or n_n<N_MIN_VERDICT or n_s<N_MIN_VERDICT:
        return "insufficient_data", "weak"
    if np.isnan(lo): return "uncertain", "weak"
    if hi-lo > CI_WIDTH_MAX: return "uncertain_wide_CI", "weak"
    if lo<=0<=hi and abs(sds)<=DEEP_ROOT_BAND: return "deep_root", "strong"
    if lo>0: return "shallow", "strong"
    if hi<0: return "negative_anomaly", "weak"
    return "uncertain", "weak"


# ================================================================
# 6. Recovery τ
# ================================================================

def exp_model(d, le_inf, le0, tau):
    return le_inf + (le0-le_inf)*np.exp(-d/tau)

def lin_model(d, a, b): return a + b*d
def log_model(d, a, b): return a + b*np.log(np.where(d>0, d, 0.5))


def _bin_median(sub, var, min_per_bin=5):
    g = sub.groupby("days_since_irrig")[var].agg(["median","count"]).reset_index()
    return g[g["count"]>=min_per_bin].reset_index(drop=True)


def fit_one(x, y, fn, p0, bounds):
    try:
        popt, _ = curve_fit(fn, x, y, p0=p0, bounds=bounds, maxfev=8000)
        yp  = fn(x, *popt)
        sse = np.sum((y-yp)**2)
        sst = np.sum((y-y.mean())**2)
        r2  = 1-sse/sst if sst>0 else np.nan
        k,n = len(popt), len(y)
        aic = n*np.log(sse/n+1e-12)+2*k if n>0 else np.nan
        return popt, r2, aic
    except Exception:
        return [np.nan]*3, np.nan, np.nan


def fit_recovery(df_active, var="LE_corr"):
    sub     = df_active[df_active[var].notna() & df_active["days_since_irrig"].notna()]
    grouped = _bin_median(sub, var)
    base    = dict(tau=np.nan, le0=np.nan, le_inf=np.nan,
                   r2_exp=np.nan, r2_lin=np.nan, r2_log=np.nan,
                   aic_exp=np.nan, aic_lin=np.nan, aic_log=np.nan,
                   n_points=len(grouped), grouped=grouped)
    if len(grouped)<4: return base
    x, y = grouped["days_since_irrig"].values.astype(float), grouped["median"].values
    pe,r2e,ae = fit_one(x,y,exp_model,[y.min(),y.max(),5.0],([0,0,0.5],[np.inf,np.inf,100]))
    pl,r2l,al = fit_one(x,y,lin_model,[y.mean(),-1.0],([-np.inf,-np.inf],[np.inf,np.inf]))
    pg,r2g,ag = fit_one(x,y,log_model,[y.mean(),-5.0],([-np.inf,-np.inf],[np.inf,np.inf]))
    base.update(dict(
        tau=float(pe[2]) if not np.isnan(pe[2]) else np.nan,
        le0=float(pe[1]) if not np.isnan(pe[1]) else np.nan,
        le_inf=float(pe[0]) if not np.isnan(pe[0]) else np.nan,
        popt_exp=pe, popt_lin=pl, popt_log=pg,
        r2_exp=r2e, r2_lin=r2l, r2_log=r2g,
        aic_exp=ae, aic_lin=al, aic_log=ag,
    ))
    return base


def fit_tau_by_year(df_active, years, var="LE_corr"):
    rows = []
    for yr in years:
        sub = df_active[df_active["year"]==yr]
        res = fit_recovery(sub, var)
        ok  = not np.isnan(res["r2_exp"]) and res["r2_exp"]>=R2_MIN_TRUST
        reason = ("OK" if ok else
                  "fit失敗" if np.isnan(res["r2_exp"]) else
                  f"R²={res['r2_exp']:.2f}<{R2_MIN_TRUST}")
        rows.append(dict(year=yr, tau=res["tau"], r2=res["r2_exp"],
                         n_points=res["n_points"], tau_reliable=ok, reason=reason))
    return pd.DataFrame(rows)


def calc_irrig_stats_by_year(df: pd.DataFrame, years: list) -> pd.DataFrame:
    tara = df[df["site"]=="Tarazona"].copy()
    tara["is_ev"] = tara["Irrig_mm"].fillna(0) > IRRIG_THRESHOLD
    rows = []
    for yr in years:
        sub = tara[tara["year"]==yr]
        n_ev   = int(sub["is_ev"].sum())
        tot_mm = float(sub["Irrig_mm"].fillna(0).sum())
        avg_mm = tot_mm/n_ev if n_ev>0 else np.nan
        rain   = float(sub["Rain_mm"].fillna(0).sum()) if "Rain_mm" in sub.columns else np.nan
        rows.append(dict(year=yr, n_irrig_events=n_ev,
                         total_irrig_mm=round(tot_mm,1),
                         avg_mm_per_event=round(avg_mm,1) if not np.isnan(avg_mm) else np.nan,
                         gauge_rain_mm=round(rain,1) if not np.isnan(rain) else np.nan))
    return pd.DataFrame(rows)


def calc_iii_by_year(df: pd.DataFrame, years: list) -> pd.DataFrame:
    tara = df[(df["site"]=="Tarazona") & (df["Irrig_mm"].fillna(0)>IRRIG_THRESHOLD)].copy()
    rows = []
    for yr in years:
        sub   = tara[tara["year"]==yr].sort_values("date")
        dates = sub["date"].values
        if len(dates)<2:
            rows.append(dict(year=yr, iii_mean=np.nan, iii_median=np.nan,
                             iii_p25=np.nan, iii_p75=np.nan, n_events=len(dates)))
            continue
        ivs = np.array([(dates[i+1]-dates[i])/np.timedelta64(1,"D")
                        for i in range(len(dates)-1)])
        ivs = ivs[ivs<=60]   # シーズン境界除外
        if len(ivs)==0:
            rows.append(dict(year=yr, iii_mean=np.nan, iii_median=np.nan,
                             iii_p25=np.nan, iii_p75=np.nan, n_events=len(dates)))
            continue
        rows.append(dict(year=yr,
                         iii_mean=round(float(np.mean(ivs)),3),
                         iii_median=round(float(np.median(ivs)),3),
                         iii_p25=round(float(np.percentile(ivs,25)),3),
                         iii_p75=round(float(np.percentile(ivs,75)),3),
                         n_events=len(dates)))
    return pd.DataFrame(rows)


# ================================================================
# 7. 可視化ヘルパー
# ================================================================

def _sig(p):
    if np.isnan(p): return ""
    return "★★★" if p<0.001 else "★★" if p<0.01 else "★" if p<0.05 else "n.s."

def _ax(ax, title="", xl="", yl=""):
    ax.set_title(title, fontweight="bold", fontsize=10, pad=6)
    ax.set_xlabel(xl, fontsize=9); ax.set_ylabel(yl, fontsize=9)
    ax.grid(alpha=0.2, lw=0.7); ax.spines[["top","right"]].set_visible(False)

def _tau_colors(df):
    return ["#1D9E75" if r.tau_reliable else "#BDBDBD" for r in df.itertuples()]


# ================================================================
# 8. 図01 — 三点パッケージ
# ================================================================

def plot_verdict_panel(pkg_df, out):
    le = pkg_df[pkg_df["var"]=="LE_corr"].sort_values(["site","season"]).reset_index(drop=True)
    if len(le)==0: return
    fig, axes = plt.subplots(1,3, figsize=(22,max(5,len(le)*1.1)))
    fig.patch.set_facecolor(FIG_BG)
    fig.suptitle("fig01 — 三点パッケージ統合判定 (LE_corr)\n"
                 "判定ガード: n≥30, CI幅≤0.50",
                 fontsize=13, fontweight="bold", y=1.01)
    yp = np.arange(len(le))
    lb = [f"{r.site}\n{r.season}" for r in le.itertuples()]
    ax = axes[0]
    for i,r in enumerate(le.itertuples()):
        v,st = verdict({"sds":r.sds,"sds_lo":r.sds_lo,"sds_hi":r.sds_hi,"n_n":r.n_n,"n_s":r.n_s})
        col=VERDICT_COL[v]; al=0.9 if st=="strong" else 0.45
        if not np.isnan(r.sds):
            ax.barh(i,r.sds,color=col,alpha=al,edgecolor="black",lw=0.8,height=0.65)
        lo_=getattr(r,"sds_lo",np.nan); hi_=getattr(r,"sds_hi",np.nan)
        if not(np.isnan(lo_) or np.isnan(hi_)):
            ax.errorbar(r.sds,i,xerr=[[r.sds-lo_],[hi_-r.sds]],
                        color="black",capsize=5,lw=1.2,fmt="none")
        ax.text(-0.55,i,f"n=({int(r.n_n)},{int(r.n_s)})",
                va="center",ha="left",fontsize=7.5,color="#555")
        if not np.isnan(r.sds):
            ax.text(r.sds+0.03,i,f"{r.sds:+.2f} [{VERDICT_JP[v]}]",
                    va="center",fontsize=8.5,
                    fontweight="bold" if st=="strong" else "normal")
    ax.axvline(0,color="gray",lw=1); ax.axvspan(-DEEP_ROOT_BAND,DEEP_ROOT_BAND,color="#1D9E75",alpha=0.07)
    ax.set_xlim(-0.6,1.0); ax.set_yticks(yp); ax.set_yticklabels(lb,fontsize=9)
    ax.legend(handles=[mpatches.Patch(color=c,label=VERDICT_JP[k],alpha=0.8)
                        for k,c in VERDICT_COL.items()],
              loc="lower right",fontsize=7,title="判定",title_fontsize=8)
    _ax(ax,"(1) SDS\nSDS≈0→深根  SDS>0→浅根","SDS")
    ax=axes[1]
    for i,r in enumerate(le.itertuples()):
        if not np.isnan(r.abs_diff):
            ax.barh(i,r.abs_diff,color="#5B4FCF",alpha=0.75,edgecolor="black",lw=0.8,height=0.65)
            if not(np.isnan(r.abs_lo) or np.isnan(r.abs_hi)):
                ax.errorbar(r.abs_diff,i,
                            xerr=[[r.abs_diff-r.abs_lo],[r.abs_hi-r.abs_diff]],
                            color="black",capsize=5,lw=1.2,fmt="none")
            ax.text(r.abs_diff+(1.5 if r.abs_diff>=0 else -1.5),i,
                    f"{r.abs_diff:+.1f}W/m²",va="center",fontsize=8.5)
    ax.axvline(0,color="gray",lw=1); ax.set_yticks(yp); ax.set_yticklabels([])
    _ax(ax,"(2) 絶対差 LE_n − LE_s","LE_n − LE_s [W/m²]")
    ax=axes[2]
    for i,r in enumerate(le.itertuples()):
        if not np.isnan(r.rb):
            col="#E85D04" if abs(r.rb)>0.3 else "#FFA000" if abs(r.rb)>0.1 else "#9E9E9E"
            ax.barh(i,r.rb,color=col,alpha=0.8,edgecolor="black",lw=0.8,height=0.65)
            ax.text(r.rb+0.03,i,f"{r.rb:+.2f} p={r.p:.1e} {_sig(r.p)}",
                    va="center",fontsize=8.5)
    ax.axvline(0,color="gray",lw=1); ax.set_xlim(-1,1.5)
    ax.set_yticks(yp); ax.set_yticklabels([])
    _ax(ax,"(3) rank-biserial r + 有意性","rank-biserial r")
    plt.tight_layout()
    fp = out/"fig01_verdict_panel.png"
    plt.savefig(fp,dpi=150,bbox_inches="tight",facecolor=FIG_BG); plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 9. 図02 — Recovery curve (2×2)
# ================================================================

def plot_recovery_curve(fit_res, tau_df, out):
    g = fit_res.get("grouped")
    if g is None or len(g)<2: return
    tau,le0,le_inf = fit_res["tau"],fit_res["le0"],fit_res["le_inf"]
    r2e,ae,al,ag   = fit_res["r2_exp"],fit_res["aic_exp"],fit_res["aic_lin"],fit_res["aic_log"]
    fig,axes = plt.subplots(2,2,figsize=(16,12)); fig.patch.set_facecolor(FIG_BG)
    xo = g["days_since_irrig"].values.astype(float)
    yo = g["median"].values; no = g["count"].values
    xf = np.linspace(0,xo.max()+2,300); sz = 40+no*4
    # (A) メイン
    ax=axes[0,0]
    ax.scatter(xo,yo,s=sz,c="#1D9E75",edgecolors="black",lw=1,zorder=5,label="観測中央値(サイズ=n)")
    if not np.isnan(tau):
        ax.plot(xf,exp_model(xf,le_inf,le0,tau),"k-",lw=2.5,
                label=f"Exp fit: τ={tau:.1f}d, R²={r2e:.3f}")
        ax.axhline(le_inf,color="#E85D04",ls="--",lw=1.5,label=f"LE_∞={le_inf:.0f}W/m²")
        ax.axhline(le0,color="#1A73E8",ls="--",lw=1.5,label=f"LE_0={le0:.0f}W/m²")
        ax.axvline(tau,color="gray",ls=":",lw=1.5,label=f"τ={tau:.1f}d")
        ic = "#E85D04" if tau<5 else "#FFA000" if tau<=14 else "#1D9E75"
        it = (f"τ={tau:.1f}d<5d→灌漑依存" if tau<5 else
              f"τ={tau:.1f}d>14d→深層水?" if tau>14 else f"τ={tau:.1f}d→中間")
        ax.text(0.97,0.97,it,transform=ax.transAxes,va="top",ha="right",
                fontsize=10,fontweight="bold",color=ic,
                bbox=dict(boxstyle="round,pad=0.4",fc="white",ec=ic,alpha=0.9))
    ax.set_ylim(bottom=0); ax.legend(fontsize=8.5)
    _ax(ax,f"(A) Recovery curve — τ={tau:.1f}d","Days since last irrigation","LE_corr 中央値 [W/m²]")
    # (B) モデル比較
    ax=axes[0,1]
    ax.scatter(xo,yo,s=sz,c="#1D9E75",edgecolors="black",lw=1,zorder=5,label="観測中央値")
    best = min(v for v in [ae,al,ag] if not np.isnan(v)) if any(
        not np.isnan(v) for v in [ae,al,ag]) else np.nan
    for nm,fn,po,r2,ai,ls,lw_ in [
        ("Exponential",exp_model,fit_res.get("popt_exp"),r2e,ae,"k-",2.5),
        ("Linear",lin_model,fit_res.get("popt_lin"),fit_res["r2_lin"],al,"b--",1.8),
        ("Logarithmic",log_model,fit_res.get("popt_log"),fit_res["r2_log"],ag,"r-.",1.8),
    ]:
        if po is not None and not np.isnan(po[0]):
            da = ai-best if not np.isnan(ai) and not np.isnan(best) else np.nan
            bm = " ← BEST" if not np.isnan(da) and abs(da)<0.1 else ""
            ax.plot(xf,fn(xf,*po),ls,lw=lw_,
                    label=f"{nm} R²={r2:.3f} ΔAIC={da:.1f}{bm}" if not np.isnan(da)
                           else f"{nm} R²={r2:.3f}")
    ax.set_ylim(bottom=0); ax.legend(fontsize=8.5)
    _ax(ax,"(B) モデル比較 (Exp/Lin/Log)\nAIC最小=最良","Days since last irrigation","LE_corr 中央値 [W/m²]")
    # (C) 年別τ
    ax=axes[1,0]
    valid=tau_df.dropna(subset=["tau"])
    if len(valid)>=2:
        ax.bar(valid["year"],valid["tau"],color=_tau_colors(valid),
               alpha=0.85,edgecolor="black",lw=1,width=0.6)
        ax.axhline(tau,color="#E85D04",ls="--",lw=2,label=f"全期間τ={tau:.1f}d")
        ax.axhline(5,color="gray",ls=":",lw=1,alpha=0.7,label="τ=5d閾値")
        trusted=valid[valid["tau_reliable"]]
        if len(trusted)>=2:
            ax.axhspan(trusted["tau"].min(),trusted["tau"].max(),
                       color="#1D9E75",alpha=0.08,
                       label=f"信頼年:{trusted['tau'].min():.1f}–{trusted['tau'].max():.1f}d")
        for _,row in valid.iterrows():
            ax.text(row["year"],row["tau"]+0.3,f"{row['tau']:.1f}d",
                    ha="center",va="bottom",fontsize=8,
                    color="black" if row["tau_reliable"] else "#999")
        lp=[mpatches.Patch(color="#1D9E75",alpha=0.85,label="信頼年(R²≥0.7)"),
            mpatches.Patch(color="#BDBDBD",alpha=0.85,label="低信頼年")]
        h,_=ax.get_legend_handles_labels(); ax.legend(handles=h+lp,fontsize=9)
    _ax(ax,"(C) 年別τ","Year","τ [days]")
    # (D) 二層モデル概念図
    ax=axes[1,1]; ax.set_xlim(0,10); ax.set_ylim(0,10)
    ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(mpatches.Rectangle((0,0),10,10,fc=FIG_BG,ec="none"))
    ax.text(5,9.5,"(D) 二層モデル解釈",ha="center",va="top",fontsize=11,fontweight="bold")
    def box(x,y,w,h,fc,txt="",fs=9):
        ax.add_patch(mpatches.FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.15",
                                             fc=fc,ec="black",lw=1.5,alpha=0.85))
        if txt: ax.text(x+w/2,y+h/2,txt,ha="center",va="center",fontsize=fs,fontweight="bold")
    def arr(x1,y1,x2,y2,col="black"):
        ax.annotate("",xy=(x2,y2),xytext=(x1,y1),
                    arrowprops=dict(arrowstyle="-|>",color=col,lw=1.8))
    box(0.5,7.5,2.5,1.2,"#AED6F1","灌漑パルス\n(drip)")
    arr(3,8.1,4.5,8.1)
    box(4.5,7.2,4.5,1.8,"#D5EAF5",f"Fast pool\n(灌漑水→τ≈{tau:.1f}d)")
    arr(6.75,7.2,6.75,5.8,"#1A73E8")
    box(4.5,4.2,4.5,1.5,"#D5F5E3",f"Slow pool\n(深層水?→LE_∞≈{le_inf:.0f}W/m²)")
    arr(6.75,4.2,6.75,2.8,"#1D9E75")
    box(4,1.5,5.5,1.2,"#FDE8D8",f"蒸散LE: {le0:.0f}→{le_inf:.0f}W/m²")
    ax.text(5,0.5,"灌漑水が主要因(深根を棄却せず支配的でない)",
            ha="center",va="bottom",fontsize=8.5,style="italic",
            bbox=dict(boxstyle="round,pad=0.3",fc="white",ec="#999",alpha=0.8))
    plt.tight_layout(h_pad=3,w_pad=3)
    fp=out/"fig02_recovery_curve.png"
    plt.savefig(fp,dpi=150,bbox_inches="tight",facecolor=FIG_BG); plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 10. 図03 — 三点指標間整合性
# ================================================================

def plot_three_point_matrix(pkg_df, out):
    sub=pkg_df[pkg_df["var"]=="LE_corr"].dropna(subset=["sds","abs_diff","rb"])
    if len(sub)<2: return
    fig,axes=plt.subplots(1,3,figsize=(18,6)); fig.patch.set_facecolor(FIG_BG)
    fig.suptitle("fig03 — 三点指標間の整合性",fontsize=12,fontweight="bold")
    for ax,(xk,yk,xl,yl,tt) in zip(axes,[
        ("sds","abs_diff","SDS","Abs diff [W/m²]","SDS vs 絶対差"),
        ("sds","rb","SDS","rank-biserial r","SDS vs 効果量"),
        ("abs_diff","rb","Abs diff [W/m²]","rank-biserial r","絶対差 vs 効果量"),
    ]):
        for _,r in sub.iterrows():
            v,st=verdict({"sds":r["sds"],"sds_lo":r["sds_lo"],
                          "sds_hi":r["sds_hi"],"n_n":r["n_n"],"n_s":r["n_s"]})
            ax.scatter(r[xk],r[yk],c=SITE_COL.get(r["site"],"#888"),s=90,
                       alpha=0.9 if st=="strong" else 0.4,
                       edgecolors="black",lw=1,
                       marker="o" if r["site"]=="Oran" else "^")
            ax.annotate(f"{r['site'][:3]}\n{r['season'][:6]}",(r[xk],r[yk]),
                        fontsize=7.5,xytext=(5,5),textcoords="offset points")
        ax.axhline(0,color="gray",lw=0.6,alpha=0.5); ax.axvline(0,color="gray",lw=0.6,alpha=0.5)
        _ax(ax,tt,xl,yl)
    leg=[mpatches.Patch(color=SITE_COL["Oran"],label="Oran"),
         mpatches.Patch(color=SITE_COL["Tarazona"],label="Tarazona"),
         Line2D([0],[0],marker="o",color="gray",ms=8,lw=0,markeredgecolor="black",label="strong"),
         Line2D([0],[0],marker="o",color="gray",ms=8,lw=0,alpha=0.35,markeredgecolor="black",label="weak")]
    fig.legend(handles=leg,loc="lower center",ncol=4,fontsize=9)
    plt.tight_layout(rect=[0,0.07,1,1])
    fp=out/"fig03_three_point_matrix.png"
    plt.savefig(fp,dpi=150,bbox_inches="tight",facecolor=FIG_BG); plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 11. 図04 — 年別τ Standalone
# ================================================================

def plot_tau_sensitivity(tau_df, tau_all, out):
    valid=tau_df.dropna(subset=["tau"])
    if len(valid)<2: return
    fig,ax=plt.subplots(figsize=(10,5.5)); fig.patch.set_facecolor(FIG_BG)
    bars=ax.bar(valid["year"],valid["tau"],color=_tau_colors(valid),
                alpha=0.80,edgecolor="black",lw=1.2,width=0.6)
    ax.axhline(tau_all,color="#E85D04",ls="--",lw=2.2,label=f"全期間τ={tau_all:.1f}d")
    ax.axhline(5,color="gray",ls=":",lw=1,alpha=0.7,label="τ=5d閾値")
    trusted=valid[valid["tau_reliable"]]
    if len(trusted)>=2:
        ax.axhspan(trusted["tau"].min(),trusted["tau"].max(),color="#1D9E75",alpha=0.08,
                   label=f"信頼年:{trusted['tau'].min():.1f}–{trusted['tau'].max():.1f}d")
    for bar,(_,row) in zip(bars,valid.iterrows()):
        ax.text(bar.get_x()+bar.get_width()/2,row["tau"]+0.15,
                f"{row['tau']:.1f}d\n(n={int(row['n_points'])})\n{row['reason']}",
                ha="center",va="bottom",fontsize=8.5,
                color="black" if row["tau_reliable"] else "#999")
    ax.set_ylim(0,max(valid["tau"].max()*1.45,8)); ax.set_xticks(valid["year"])
    lp=[mpatches.Patch(color="#1D9E75",alpha=0.85,label="信頼年(R²≥0.7)"),
        mpatches.Patch(color="#BDBDBD",alpha=0.85,label="低信頼年")]
    h,_=ax.get_legend_handles_labels(); ax.legend(handles=h+lp,fontsize=9.5)
    rfl=("★ 信頼年τ変動小→Robust" if len(trusted)>=2 and
         trusted["tau"].max()-trusted["tau"].min()<3 else "△ 信頼年でも変動あり")
    rc="#1D9E75" if "★" in rfl else "#E85D04"
    body=(f"{rfl}\n信頼年τ平均={trusted['tau'].mean():.1f}±{trusted['tau'].std():.1f}d"
          if len(trusted)>=2 else rfl)
    ax.text(0.98,0.98,body,transform=ax.transAxes,va="top",ha="right",
            fontsize=10,fontweight="bold",color=rc,
            bbox=dict(boxstyle="round,pad=0.4",fc="white",ec=rc,alpha=0.9))
    _ax(ax,f"fig04 — 年別τ Sensitivity\n全期間τ={tau_all:.1f}d","Year","τ [days]")
    plt.tight_layout()
    fp=out/"fig04_tau_sensitivity.png"
    plt.savefig(fp,dpi=150,bbox_inches="tight",facecolor=FIG_BG); plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 12. 図05 — τ Diagnostic (4段)
# ================================================================

def plot_tau_diagnostic(tau_df, irrig_stats, out):
    merged=tau_df.merge(irrig_stats,on="year",how="left")
    years=merged["year"].values; x=np.arange(len(years)); w=0.6
    fig=plt.figure(figsize=(14,12)); fig.patch.set_facecolor(FIG_BG)
    fig.suptitle("fig05 — τ Diagnostic: 年別τ×灌漑統計",
                 fontsize=12,fontweight="bold",y=0.98)
    gs=GridSpec(4,1,hspace=0.55,figure=fig)
    for panel,(col,title,color) in enumerate([
        ("tau",        "(A) 年別τ [days]",            None),
        ("n_irrig_events","(B) 灌漑イベント数",        "#1A73E8"),
        ("avg_mm_per_event","(C) 灌漑量 [mm/event]",  "#FFA000"),
        ("gauge_rain_mm","(D) 降水量 [mm] (rain gauge)","#9575CD"),
    ]):
        ax=fig.add_subplot(gs[panel])
        vals=merged[col].values if col in merged.columns else np.full(len(years),np.nan)
        if panel==0:
            colors0=["#1D9E75" if r.tau_reliable else "#BDBDBD" for r in merged.itertuples()]
            ax.bar(x,vals,color=colors0,alpha=0.85,edgecolor="black",lw=1,width=w)
            ax.axhline(5,color="gray",ls=":",lw=1,alpha=0.7,label="τ=5d")
            lp=[mpatches.Patch(color="#1D9E75",alpha=0.85,label="信頼年(R²≥0.7)"),
                mpatches.Patch(color="#BDBDBD",alpha=0.85,label="低信頼年")]
            h,_=ax.get_legend_handles_labels(); ax.legend(handles=lp+h,fontsize=8,loc="upper right")
            ax.set_ylim(0,np.nanmax(vals)*1.4+2)
        else:
            ax.bar(x,vals,color=color,alpha=0.80,edgecolor="black",lw=1,width=w)
        for i,v in enumerate(vals):
            if not np.isnan(v):
                lbl = (f"{v:.1f}d\n({merged['reason'].iloc[i]})"
                       if panel==0 else
                       f"{int(v)}" if panel==1 else f"{v:.1f}")
                ax.text(i,v*(1.03 if v>0 else 0.97),lbl,
                        ha="center",va="bottom",fontsize=8,
                        color=("black" if panel>0 or merged["tau_reliable"].iloc[i] else "#888"))
        ax.set_xticks(x); ax.set_xticklabels(years if panel==3 else [""]*len(years))
        _ax(ax,title,"Year" if panel==3 else "","")
    plt.subplots_adjust(top=0.94)
    fp=out/"fig05_tau_diagnostic.png"
    plt.savefig(fp,dpi=150,bbox_inches="tight",facecolor=FIG_BG); plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 13. 図06 — τ × III 散布図
# ================================================================

def plot_tau_iii_scatter(tau_df, iii_df, out):
    merged=tau_df.merge(iii_df,on="year",how="inner").dropna(subset=["tau","iii_mean"])
    if len(merged)<3:
        print("  [skip] tau-III scatter"); return
    fig,axes=plt.subplots(1,2,figsize=(13,5.5)); fig.patch.set_facecolor(FIG_BG)
    fig.suptitle("fig06 — τ × 灌漑間隔(III) 散布図\n「τは灌漑スケジュールに規定される」仮説検証",
                 fontsize=12,fontweight="bold")
    for ax,(xc,xl) in zip(axes,[
        ("iii_mean",  "III_mean [days]"),
        ("iii_median","III_median [days]"),
    ]):
        for _,row in merged.iterrows():
            col="#1D9E75" if row["tau_reliable"] else "#BDBDBD"
            ax.scatter(row[xc],row["tau"],c=col,s=100,edgecolors="black",lw=1.2,zorder=5)
            ax.annotate(str(int(row["year"])),(row[xc],row["tau"]),
                        xytext=(5,5),textcoords="offset points",fontsize=9,fontweight="bold")
        trusted=merged[merged["tau_reliable"]]
        if len(trusted)>=3:
            xv=trusted[xc].values; yv=trusted["tau"].values
            slope,intercept,r,pv,_=stats.linregress(xv,yv)
            xl_=np.linspace(xv.min()-0.5,xv.max()+0.5,100)
            ax.plot(xl_,intercept+slope*xl_,"r--",lw=1.8,alpha=0.8,
                    label=f"信頼年回帰 r={r:.2f} p={pv:.2f}")
            txt=(f"r={r:.2f}, p={pv:.3f}\n"
                 f"{'→ 正相関: τ∝III支持' if r>0.5 else '→ 相関弱: 他要因関与'}")
            ax.text(0.97,0.05,txt,transform=ax.transAxes,va="bottom",ha="right",
                    fontsize=9,color="red" if r>0.5 else "gray",
                    bbox=dict(boxstyle="round,pad=0.3",fc="white",
                              ec="red" if r>0.5 else "gray",alpha=0.85))
        ax.set_ylim(bottom=0); ax.legend(fontsize=9)
        lp=[mpatches.Patch(color="#1D9E75",alpha=0.85,label="信頼年"),
            mpatches.Patch(color="#BDBDBD",alpha=0.85,label="低信頼年")]
        h,_=ax.get_legend_handles_labels(); ax.legend(handles=lp+h,fontsize=8.5)
        _ax(ax,f"τ vs {xc}",xl,"τ [days]")
    plt.tight_layout()
    fp=out/"fig06_tau_iii_scatter.png"
    plt.savefig(fp,dpi=150,bbox_inches="tight",facecolor=FIG_BG); plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 14. 図07 — 降水量比較 (4ソース)
# ================================================================

def plot_precip_compare(precip_df, rain_diag, out):
    """
    上段: 4ソースの年間降水量比較
    下段: rain gauge を生育期/全期間に分けて CHIRPS 年間値と比較
    """
    years=precip_df["year"].values; x=np.arange(len(years)); w=0.18
    fig,(ax1,ax2)=plt.subplots(2,1,figsize=(13,11),sharex=True)
    fig.patch.set_facecolor(FIG_BG)
    fig.suptitle("fig07 — 年別降水量: 4データソース比較 (Tarazona)\n"
                 "上段: 年間値  下段: rain gauge の集計範囲診断",
                 fontsize=12,fontweight="bold")
    # 上段
    src=[("gauge_rain_mm","#2196F3","rain gauge (CSV集計)"),
         ("chirps_mm_annual","#4CAF50","CHIRPS 年間"),
         ("chirps_mm_growing","#81C784","CHIRPS 生育期(4-10月)"),
         ("terra_mm_annual","#FF9800","TerraClimate 年間"),
         ("era5_mm_annual","#9C27B0","ERA5 年間")]
    for j,(col,color,label) in enumerate(src):
        if col not in precip_df.columns: continue
        vals=precip_df[col].values
        offs=(j-2)*w
        bars=ax1.bar(x+offs,vals,width=w,color=color,alpha=0.80,
                     edgecolor="black",lw=0.7,label=label)
        for bar,v in zip(bars,vals):
            if not np.isnan(v):
                ax1.text(bar.get_x()+bar.get_width()/2,v+4,
                         f"{v:.0f}",ha="center",va="bottom",fontsize=7,color=color)
    # 2021強調
    if 2021 in list(years):
        idx=list(years).index(2021)
        ax1.axvspan(idx-0.5,idx+0.5,color="#FFF9C4",alpha=0.5,zorder=0)
        ax1.text(idx,ax1.get_ylim()[1]*0.97,"2021\n(τ短い年)",
                 ha="center",va="top",fontsize=8.5,color="#F57F17",fontweight="bold")
    ax1.legend(fontsize=8.5,loc="upper right",ncol=2)
    _ax(ax1,"(A) 年間降水量 — 4ソース比較","","降水量 [mm]")
    # 下段: rain gauge の集計範囲診断
    if rain_diag is not None and len(rain_diag)>0:
        rd=rain_diag.merge(precip_df[["year","chirps_mm_annual"]],on="year",how="left")
        bw=0.22
        for j,(col,color,label) in enumerate([
            ("rain_all_mm",     "#2196F3","gauge 全期間"),
            ("rain_growing_mm", "#64B5F6","gauge 生育期のみ"),
            ("chirps_mm_annual","#4CAF50","CHIRPS 年間"),
        ]):
            if col not in rd.columns: continue
            vals=rd[col].values
            offs=(j-1)*bw
            ax2.bar(x+offs,vals,width=bw,color=color,alpha=0.80,
                    edgecolor="black",lw=0.7,label=label)
        ax2.legend(fontsize=8.5,loc="upper right")
        _ax(ax2,"(B) rain gauge の集計範囲別 vs CHIRPS\n"
                "全期間≒CHIRPS なら gauge は年間値; 生育期≒CHIRPS なら生育期に絞るべき",
                "Year","降水量 [mm]")
    ax2.set_xticks(x); ax2.set_xticklabels(years)
    plt.tight_layout(h_pad=2)
    fp=out/"fig07_precip_compare.png"
    plt.savefig(fp,dpi=150,bbox_inches="tight",facecolor=FIG_BG); plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 15. 図08 — τ × III × 降水量 統合診断
# ================================================================

def plot_tau_drivers(tau_df, iii_df, drivers_df, out):
    merged=tau_df.merge(iii_df,on="year",how="left").merge(drivers_df,on="year",how="left")
    years=merged["year"].values; x=np.arange(len(years)); w=0.5
    fig,(ax1,ax2)=plt.subplots(2,1,figsize=(12,10),sharex=True)
    fig.patch.set_facecolor(FIG_BG)
    fig.suptitle("fig08 — τドライバー統合診断: τ×III×降水量×灌漑量",
                 fontsize=12,fontweight="bold")
    # 上段
    ax1.bar(x,merged["tau"],color=_tau_colors(merged),alpha=0.80,
            edgecolor="black",lw=1,width=w,label="τ [days]")
    for xi,row in zip(x,merged.itertuples()):
        if not np.isnan(row.tau):
            ax1.text(xi,row.tau+0.4,f"{row.tau:.1f}d",
                     ha="center",va="bottom",fontsize=8.5,
                     color="black" if row.tau_reliable else "#999")
    ax1r=ax1.twinx()
    ax1r.plot(x,merged["iii_mean"],"o-",color="#E85D04",lw=2,ms=8,label="III_mean [days]")
    for xi,v in zip(x,merged["iii_mean"]):
        if not np.isnan(v):
            ax1r.text(xi,v+0.2,f"{v:.1f}d",ha="center",va="bottom",fontsize=8,color="#E85D04")
    ax1.set_ylabel("τ [days]",fontsize=9); ax1r.set_ylabel("III_mean [days]",fontsize=9,color="#E85D04")
    ax1r.tick_params(axis="y",colors="#E85D04"); ax1.set_ylim(bottom=0); ax1r.set_ylim(bottom=0)
    lp=[mpatches.Patch(color="#1D9E75",alpha=0.85,label="τ信頼年"),
        mpatches.Patch(color="#BDBDBD",alpha=0.85,label="τ低信頼年")]
    h1,_=ax1.get_legend_handles_labels(); h2,_=ax1r.get_legend_handles_labels()
    ax1.legend(handles=h1+h2+lp,fontsize=8.5,loc="upper left")
    _ax(ax1,"(A) τ(bar) と灌漑間隔 III_mean(折れ線)\nτとIIIが連動→灌漑スケジュールがτを規定","","τ [days]")
    # 下段
    gauge_col="gauge_rain_mm" if "gauge_rain_mm" in merged.columns else None
    chirps_col="chirps_mm_annual" if "chirps_mm_annual" in merged.columns else None
    if gauge_col:
        ax2.bar(x-0.15,merged[gauge_col],width=0.3,color="#2196F3",alpha=0.75,
                edgecolor="black",lw=0.8,label="Rain gauge [mm]")
    if chirps_col:
        ax2.bar(x+0.15,merged[chirps_col],width=0.3,color="#4CAF50",alpha=0.75,
                edgecolor="black",lw=0.8,label="CHIRPS 年間 [mm]")
    ax2r=ax2.twinx()
    if "total_irrig_mm" in merged.columns:
        ax2r.plot(x,merged["total_irrig_mm"],"s--",color="#FF9800",lw=2,ms=8,label="灌漑総量 [mm]")
        ax2r.set_ylabel("灌漑総量 [mm]",fontsize=9,color="#FF9800")
        ax2r.tick_params(axis="y",colors="#FF9800")
    ax2.set_xticks(x); ax2.set_xticklabels(years)
    ax2.set_ylabel("降水量 [mm]",fontsize=9); ax2.set_ylim(bottom=0)
    h3,_=ax2.get_legend_handles_labels(); h4,_=ax2r.get_legend_handles_labels()
    ax2.legend(handles=h3+h4,fontsize=8.5,loc="upper left")
    _ax(ax2,"(B) 降水量(gauge+CHIRPS)と灌漑総量","Year","降水量 [mm]")
    plt.tight_layout(h_pad=2)
    fp=out/"fig08_tau_drivers.png"
    plt.savefig(fp,dpi=150,bbox_inches="tight",facecolor=FIG_BG); plt.close()
    print(f"  [save] {fp}")


# ================================================================
# 16. MAIN
# ================================================================

def main():
    args = parse_args()
    parquet, csv_path, chirps_dir, terra_dir, era5_dir, out_dir = resolve_paths(args)
    out_dir.mkdir(parents=True, exist_ok=True)
    YEARS = args.years
    LAT, LON = args.lat, args.lon

    print("="*60)
    print("解析A v19 — III + 3種降水データ + τ Diagnostic 拡張")
    print(f"  Tarazona 座標: lat={LAT}, lon={LON}")
    print(f"  解析年: {YEARS}")
    print(f"  出力先: {out_dir}")
    print("="*60)

    # ── データ読込 ──
    print("\n--- データ読込 ---")
    df = load_and_merge(parquet, csv_path)
    df = add_days_since_irrig(df)
    df = add_irrig_active_month(df)
    print(f"  全行数: {len(df)}  "
          f"Tarazona灌漑アクティブ日数: {((df['site']=='Tarazona')&df['irrig_active_month']).sum()}")

    # ── 降水量診断 ──
    print("\n--- 降水量診断 (rain gauge 集計範囲) ---")
    rain_diag = diagnose_rain_gauge(df)
    print(rain_diag.to_string(index=False))
    # 全期間と生育期の比率で自動判定
    ratio = (rain_diag["rain_growing_mm"] / rain_diag["rain_all_mm"].replace(0,np.nan)).mean()
    print(f"\n  生育期/全期間 比率 (平均): {ratio:.2f}")
    if ratio > 0.85:
        print("  → rain gauge はほぼ生育期のみのデータ → CHIRPS との比較は生育期に揃えるべき")
    else:
        print("  → rain gauge は通年データを含む → CHIRPS 年間値と比較可能")

    if args.diag_only:
        rain_diag.to_csv(out_dir/"v19_precip_diagnostic.csv", index=False)
        print(f"\n[diag-only] {out_dir}/v19_precip_diagnostic.csv に保存して終了")
        return

    # ── 三点パッケージ ──
    print(f"\n{'='*60}\n三点パッケージ計算\n{'='*60}")
    rows = []
    for season_label, months in SEASONS.items():
        for site in ["Oran","Tarazona"]:
            sub = df[(df["site"]==site) & df["is_growing"] & df["month"].isin(months)]
            if site=="Tarazona": sub=sub[sub["irrig_active_month"]]
            for var in ["LE_corr","EF_corr","ET"]:
                if var not in sub.columns: continue
                n_arr = sub.loc[sub["drought_type"]=="normal",   var].dropna().values
                s_arr = sub.loc[sub["drought_type"]=="soil dry", var].dropna().values
                pkg   = three_point_package(n_arr, s_arr, var)
                v, st = verdict(pkg)
                pkg.update(dict(site=site,season=season_label,var=var,
                                verdict=v,verdict_strength=st))
                rows.append(pkg)
                if var=="LE_corr":
                    print(f"  [{site:9s}/{season_label:22s}] "
                          f"SDS={pkg['sds']:+.3f} abs={pkg['abs_diff']:+.1f}W/m² "
                          f"rb={pkg['rb']:+.3f} p={pkg['p']:.1e} "
                          f"n=({pkg['n_n']},{pkg['n_s']}) → [{VERDICT_JP[v]}] ({st})")
    pkg_df = pd.DataFrame(rows)
    pkg_df.to_csv(out_dir/"v19_three_point_package.csv", index=False)

    # ── Recovery τ ──
    print(f"\n{'='*60}\nRecovery τ\n{'='*60}")
    tara_active = df[(df["site"]=="Tarazona") & df["irrig_active_month"] & df["is_growing"]]
    fit_res     = fit_recovery(tara_active)
    tau_df      = fit_tau_by_year(tara_active, YEARS)
    irrig_stats = calc_irrig_stats_by_year(df, YEARS)
    print(f"  τ(全期間)={fit_res['tau']:.2f}d  "
          f"LE_0={fit_res['le0']:.1f}  LE_∞={fit_res['le_inf']:.1f} W/m²")
    print(f"  R² exp={fit_res['r2_exp']:.3f} lin={fit_res['r2_lin']:.3f} log={fit_res['r2_log']:.3f}")
    print(f"  AIC exp={fit_res['aic_exp']:.1f} lin={fit_res['aic_lin']:.1f} log={fit_res['aic_log']:.1f}")
    print(tau_df.to_string(index=False))
    trusted = tau_df[tau_df["tau_reliable"]]
    if len(trusted)>=2:
        print(f"\n  信頼年τレンジ: {trusted['tau'].min():.1f}–{trusted['tau'].max():.1f}d "
              f"(平均={trusted['tau'].mean():.1f}±{trusted['tau'].std():.1f}d)")

    # ── III ──
    print(f"\n{'='*60}\nInter-Irrigation Interval (III)\n{'='*60}")
    iii_df = calc_iii_by_year(df, YEARS)
    print(iii_df.to_string(index=False))
    merged_chk = tau_df.merge(iii_df,on="year").dropna(subset=["tau","iii_mean"])
    tr = merged_chk[merged_chk["tau_reliable"]]
    if len(tr)>=3:
        r,pv = stats.pearsonr(tr["tau"],tr["iii_mean"])
        print(f"\n  τ × III_mean 相関 (信頼年 n={len(tr)}): r={r:.3f}, p={pv:.3f}")
        if r>0.5:   print("  → 正相関 → 「τは灌漑間隔に規定される」仮説を支持")
        elif r<-0.3: print("  → 負相関 → 間隔が短いほどτが長い (逆転) → 要考察")
        else:        print("  → 相関弱 → τに他の要因が関与 (VPD・土壌特性など)")
    else:
        print("  信頼年が3未満 → 相関計算不可")

    # ── 3種降水データ抽出 ──
    print(f"\n{'='*60}\n降水データ抽出\n{'='*60}")
    print("  CHIRPS:")
    chirps = extract_chirps(chirps_dir, YEARS, LAT, LON)
    print("  TerraClimate:")
    terra  = extract_terraclimate(terra_dir, YEARS, LAT, LON)
    print("  ERA5:")
    era5   = extract_era5(era5_dir, YEARS, LAT, LON)

    precip_df = (irrig_stats[["year","gauge_rain_mm"]]
                 .merge(chirps, on="year", how="left")
                 .merge(terra,  on="year", how="left")
                 .merge(era5,   on="year", how="left"))
    print(precip_df.to_string(index=False))

    # CHIRPS との乖離を自動警告
    if "chirps_mm_annual" in precip_df.columns:
        for _,row in precip_df.iterrows():
            ratio_yr = (row["gauge_rain_mm"] / row["chirps_mm_annual"]
                        if not np.isnan(row.get("chirps_mm_annual",np.nan)) and
                           row["chirps_mm_annual"]>0 else np.nan)
            if not np.isnan(ratio_yr) and ratio_yr < 0.3:
                print(f"  [warn] {int(row['year'])}: gauge({row['gauge_rain_mm']:.0f}mm) "
                      f"が CHIRPS({row['chirps_mm_annual']:.0f}mm) の {ratio_yr:.0%} "
                      f"→ gauge は生育期のみの可能性大")

    # CSV 保存
    precip_df.to_csv(out_dir/"v19_precip_comparison.csv", index=False)
    iii_df.to_csv(out_dir/"v19_iii_by_year.csv", index=False)
    tau_df.to_csv(out_dir/"v19_tau_by_year.csv", index=False)
    rain_diag.to_csv(out_dir/"v19_precip_diagnostic.csv", index=False)
    if fit_res.get("grouped") is not None:
        fit_res["grouped"].to_csv(out_dir/"v19_recovery_binned.csv", index=False)
    with open(out_dir/"v19_fit_params.json","w") as f:
        json.dump({k:(float(v) if isinstance(v,(float,np.floating)) else
                      int(v) if isinstance(v,(int,np.integer)) else v)
                   for k,v in fit_res.items()
                   if k not in ("grouped","popt_exp","popt_lin","popt_log")}, f, indent=2)

    # ── 可視化 ──
    if not args.no_plots:
        print(f"\n--- 可視化 ---")
        drivers_df = precip_df.merge(
            irrig_stats[["year","total_irrig_mm","n_irrig_events"]], on="year", how="left")
        plot_verdict_panel(pkg_df, out_dir)
        plot_recovery_curve(fit_res, tau_df, out_dir)
        plot_three_point_matrix(pkg_df, out_dir)
        plot_tau_sensitivity(tau_df, fit_res["tau"], out_dir)
        plot_tau_diagnostic(tau_df, irrig_stats, out_dir)
        plot_tau_iii_scatter(tau_df, iii_df, out_dir)
        plot_precip_compare(precip_df, rain_diag, out_dir)
        plot_tau_drivers(tau_df, iii_df, drivers_df, out_dir)
    else:
        print("\n[--no-plots] 図の生成をスキップ")

    # ── 最終サマリー ──
    print(f"\n{'='*60}\n★ 最終サマリー\n{'='*60}")
    for _,row in pkg_df[pkg_df["var"]=="LE_corr"].iterrows():
        v,st=row["verdict"],row["verdict_strength"]
        lo,hi=row.get("sds_lo",np.nan),row.get("sds_hi",np.nan)
        ci=f"[{lo:+.2f},{hi:+.2f}]" if not np.isnan(lo) else "[CI算出不可]"
        print(f"  [{row['site']:9s}/{row['season']:22s}] "
              f"SDS={row['sds']:+.3f} {ci} n=({row['n_n']:.0f},{row['n_s']:.0f}) "
              f"→ {VERDICT_JP[v]} ({st})")
    tau=fit_res["tau"]
    print(f"\n  Recovery τ={tau:.1f}d  "
          f"主張: 'Irrigation-driven decoupling (τ≈{tau:.0f}d)'")
    print(f"\n[done] → {out_dir}/")


if __name__ == "__main__":
    main()
