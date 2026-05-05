"""
解析A v9: 深層水アクセス検出 — 再定義版

v8 → v9 の変更点
================
[改1] deep_access を LE ベースで再定義
        旧: SWC<thr ∧ NDWI>thr ∧ VPD>thr  → 同時成立がほぼゼロ
        新: SWC<p25 ∧ LE>median(LE)        → 物理的に直接的

[改2] NDWI をアノマリ(月別中央値からの残差)で評価
        サイト依存・季節依存を除去

[改3] サイト間比較を主分析に昇格
        "dry SWC 日"の LE/EF を Oran vs Tarazona で Mann-Whitney 検定
        → ここで仮説を直接検証する

[改4] GRACE 経度を 0–360° に変換して再試行
        JPL GRACE-FO は 0–360 表記なので lon<0 は 360+lon に変換

[改5] SWC→NDWI のラグ相関を ±30 日で計算
        深根なら NDWI は SWC 低下後も "遅れて" しか落ちない

[改6] 補助的に v8 と同じ 3 軸分類も残す(参考用)
        ただし NDWI は anomaly>0 を採用してサイト依存性を除去
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats, signal
from scipy.stats import gaussian_kde
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")


# ================================================================
# 0. 設定
# ================================================================
BASE_EC  = Path("/home/shion-nagamine/Dataset/Eddy data in Spain")
BASE_NC  = Path("/mnt/hdd/Dataset")
ERA5_DIR = BASE_NC / "ERA5_2m_Temperature"

PATHS = {
    "oran_ec"  : BASE_EC / "Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.csv",
    "tara_ec"  : BASE_EC / "Daily_Summary_Filtered_forPred_ActEne26.csv",
    "oran_ndwi": BASE_NC / "Sentinel2_NDWI/Oran_NDWI_Export.csv",
    "tara_ndwi": BASE_NC / "Sentinel2_NDWI/TzM_NDWI_Export.csv",
    "grace"    : BASE_NC / "GRACE-FO_TWL/GRCTellus.JPL.200204_202602.GLO.RL06.3M.MSCNv04CRI.nc",
}

SITES = {
    "Oran"    : {"lat": 38.82,   "lon": -1.86},
    "Tarazona": {"lat": 39.266,  "lon": -1.9397},
}

GROWING_MONTHS = {
    "Oran"    : [11, 12, 1, 2, 3, 4, 5, 6],
    "Tarazona": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
}

EF_DENOM_MIN = 10
ERA5_YEARS   = list(range(2018, 2025))

SAVE_DIR = Path("./output_analysis_A_v9")
SAVE_DIR.mkdir(exist_ok=True)


# ================================================================
# 1. ERA5 VPD (v8 と同じ)
# ================================================================

def magnus_e(temp_c):
    return 6.112 * np.exp(17.67 * temp_c / (temp_c + 243.5))


def load_era5_vpd(lat, lon, years=ERA5_YEARS, era5_dir=ERA5_DIR):
    records = []
    t_vars  = ["t2m", "T2M", "temperature_2m", "t"]
    d_vars  = ["d2m", "D2M", "dewpoint_2m", "td"]
    rh_vars = ["r", "RH", "relative_humidity", "rh"]

    for year in years:
        nc_path = era5_dir / f"{year}.nc"
        if not nc_path.exists():
            print(f"  [ERA5] {year}: ファイルなし")
            continue
        try:
            ds = xr.open_dataset(nc_path, engine="netcdf4")
            drop_coords = [c for c in ds.coords
                           if c not in ("valid_time","time","latitude","longitude","lat","lon")]
            ds = ds.drop_vars(drop_coords, errors="ignore")

            def get_var(cands):
                for n in cands:
                    if n in ds: return ds[n]
                raise KeyError(cands)

            def sel_pt(da):
                lat_n = next(c for c in da.dims if "lat" in c.lower())
                lon_n = next(c for c in da.dims if "lon" in c.lower())
                return da.sel({lat_n: lat, lon_n: lon}, method="nearest")

            def to_daily(da, name):
                ser = sel_pt(da).to_series()
                ser.index = pd.to_datetime(ser.index)
                df_ = ser.reset_index()
                df_.columns = ["time", name]
                return df_

            t2m_df = to_daily(get_var(t_vars), "t2m")
            if t2m_df["t2m"].median() > 100: t2m_df["t2m"] -= 273.15
            try:
                d2m_df = to_daily(get_var(d_vars), "d2m")
                if d2m_df["d2m"].median() > 100: d2m_df["d2m"] -= 273.15
                mh = pd.merge(t2m_df, d2m_df, on="time")
                mh["VPD"] = np.maximum(magnus_e(mh["t2m"]) - magnus_e(mh["d2m"]), 0)
            except KeyError:
                rh_df = to_daily(get_var(rh_vars), "rh")
                if rh_df["rh"].median() < 2: rh_df["rh"] *= 100
                mh = pd.merge(t2m_df, rh_df, on="time")
                es = magnus_e(mh["t2m"])
                mh["VPD"] = np.maximum(es - es * mh["rh"] / 100, 0)

            mh["date"] = mh["time"].dt.normalize()
            daily = mh.groupby("date")["VPD"].mean().reset_index()
            daily.columns = ["date", "VPD_era5"]
            records.append(daily)
            print(f"  [ERA5] {year}: {len(daily)}日  median={daily['VPD_era5'].median():.2f}hPa")
            ds.close()
        except Exception as e:
            print(f"  [ERA5] {year}: error {e}")

    if not records:
        raise RuntimeError(f"ERA5 が読めません: {era5_dir}")
    out = pd.concat(records, ignore_index=True).sort_values("date")
    return out


def merge_era5_vpd(ec_df, era5_vpd):
    merged = ec_df.copy()
    if "VPD" in merged.columns:
        merged = merged.rename(columns={"VPD": "VPD_ec"})
    merged = pd.merge(merged, era5_vpd, on="date", how="left")
    merged = merged.rename(columns={"VPD_era5": "VPD"})
    print(f"  ERA5 VPD マッチ: {merged['VPD'].notna().sum()}/{len(merged)} 日")
    return merged


# ================================================================
# 2. EC 読み込み (v8 と同じ)
# ================================================================

def load_oran_ec(filepath):
    df = pd.read_csv(filepath)
    dt_col = next((c for c in df.columns if c.lower() in ("datetime","timestamp")), None)
    df["datetime"] = pd.to_datetime(df[dt_col], errors="coerce")
    df = df.dropna(subset=["datetime"]).sort_values("datetime").reset_index(drop=True)
    for col in ["SWC_1_1_1","LE","H","G","NETRAD","VPD","ET"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["date"] = df["datetime"].dt.normalize()
    daily = df.groupby("date", as_index=False).agg(
        SWC=("SWC_1_1_1","mean"), LE=("LE","mean"), H=("H","mean"),
        G=("G","mean"), Rn=("NETRAD","mean"), VPD=("VPD","mean"), ET=("ET","sum"))
    denom = daily["Rn"] - daily["G"]
    valid = denom > EF_DENOM_MIN
    daily["EF"] = np.nan
    daily.loc[valid, "EF"] = (daily.loc[valid,"LE"]/denom[valid]).clip(0,1.5)
    daily = daily[(daily["SWC"]>0) & (daily["SWC"]<100)].reset_index(drop=True)
    daily["site"] = "Oran"
    print(f"[Oran EC] {len(daily)} 日分")
    return daily


def load_tarazona_ec(filepath):
    df = pd.read_csv(filepath)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df = df.dropna(subset=["date"]).sort_values("date").reset_index(drop=True)
    df = df.rename(columns={"SWC_avg":"SWC","LE_avg":"LE","H_avg":"H","G_avg":"G",
                             "NetRad_avg":"Rn","VPD_mean":"VPD","ET_avg":"ET"})
    for col in ["SWC","LE","H","G","Rn","VPD","ET"]:
        if col in df.columns: df[col] = pd.to_numeric(df[col], errors="coerce")
    denom = df["Rn"] - df["G"]
    valid = denom > EF_DENOM_MIN
    df["EF"] = np.nan
    df.loc[valid,"EF"] = (df.loc[valid,"LE"]/denom[valid]).clip(0,1.5)
    df = df[(df["SWC"]>0) & (df["SWC"]<100)].reset_index(drop=True)
    df["site"] = "Tarazona"
    print(f"[Tarazona EC] {len(df)} 日分")
    cols = [c for c in ["date","site","SWC","LE","H","G","Rn","EF","VPD","ET"] if c in df.columns]
    return df[cols]


def normalize_swc(df, site):
    if df["SWC"].dropna().max() <= 1.0:
        df = df.copy(); df["SWC"] *= 100
        print(f"  SWC [{site}]: m³/m³ → %")
    return df


# ================================================================
# 3. NDWI 読み込み・マッチング・[改2] アノマリ計算
# ================================================================

def load_ndwi_csv(filepath, site_name):
    df = pd.read_csv(filepath)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df = df.rename(columns={"NDWI":"NDWI_s2"})
    df = df[["date","NDWI_s2"]].dropna().sort_values("date").reset_index(drop=True)
    print(f"[NDWI {site_name}] {len(df)} シーン")
    return df


def match_ndwi(ec_df, ndwi_df, max_days=5):
    ec   = ec_df.sort_values("date").reset_index(drop=True).copy()
    ndwi = ndwi_df.sort_values("date").reset_index(drop=True).copy()
    tol  = pd.Timedelta(days=max_days)
    fwd = pd.merge_asof(ec, ndwi, on="date", direction="backward", tolerance=tol)
    bwd = pd.merge_asof(ec, ndwi, on="date", direction="forward",  tolerance=tol,
                        suffixes=("","_bwd"))
    fwd_v = fwd["NDWI_s2"].values
    bwd_c = "NDWI_s2_bwd" if "NDWI_s2_bwd" in bwd.columns else "NDWI_s2"
    bwd_v = bwd[bwd_c].values
    ec["NDWI_s2"] = np.where(np.isnan(fwd_v) & np.isnan(bwd_v), np.nan,
                     np.where(np.isnan(fwd_v), bwd_v,
                     np.where(np.isnan(bwd_v), fwd_v, (fwd_v+bwd_v)/2)))
    print(f"  NDWI マッチ: {ec['NDWI_s2'].notna().sum()}/{len(ec)} 日")
    return ec


def add_ndwi_anomaly(df):
    """月別中央値からの残差として NDWI アノマリを計算 [改2]"""
    df = df.copy()
    df["month"] = df["date"].dt.month
    monthly_med = df.groupby("month")["NDWI_s2"].transform("median")
    df["NDWI_anom"] = df["NDWI_s2"] - monthly_med
    df = df.drop(columns="month")
    print(f"  NDWI_anom: median={df['NDWI_anom'].median():.4f}  "
          f"std={df['NDWI_anom'].std():.4f}")
    return df


# ================================================================
# 4. GRACE [改4: 経度を 0–360 に変換]
# ================================================================

def load_grace(filepath, lat, lon):
    ds = xr.open_dataset(filepath)
    lon_grid = ds["lon"].values
    if lon_grid.min() >= 0 and lon < 0:
        lon_use = 360 + lon
        print(f"  [GRACE] 経度変換: {lon} → {lon_use} (0–360°系)")
    else:
        lon_use = lon
    lwe = ds["lwe_thickness"].sel(lat=lat, lon=lon_use, method="nearest")
    sf  = float(ds["scale_factor"].sel(lat=lat, lon=lon_use, method="nearest").values)
    df  = (lwe * sf).to_dataframe(name="GWS_cm").reset_index()
    df  = df.rename(columns={"time":"date"})
    df["date"] = pd.to_datetime(df["date"]).astype("datetime64[ns]")
    df  = df[["date","GWS_cm"]].dropna()
    print(f"[GRACE] {len(df)} ヶ月分  (lat={lat}, lon={lon_use})")
    return df


def merge_grace(ec_df, grace_df):
    ec_s = ec_df.copy(); ec_s["date"] = ec_s["date"].astype("datetime64[ns]")
    grace_daily = (grace_df.set_index("date").resample("D")
                   .interpolate("linear").reset_index())
    merged = pd.merge_asof(ec_s.sort_values("date"), grace_daily.sort_values("date"),
                           on="date", direction="nearest", tolerance=pd.Timedelta("32D"))
    print(f"  GRACE マッチ: {merged['GWS_cm'].notna().sum()}/{len(merged)} 日")
    return merged


# ================================================================
# 5. 生育期フィルタ
# ================================================================

def filter_growing_season(df, site):
    months = GROWING_MONTHS[site]
    out = df[df["date"].dt.month.isin(months)].copy().reset_index(drop=True)
    print(f"  [生育期] {site}: {len(df)} → {len(out)} 日")
    return out


# ================================================================
# 6. [改1] LE ベースの deep_access 再定義
# ================================================================

def classify_le_based(df, le_q=0.5, swc_q=0.25):
    """
    deep_access = SWC<p_swc_q AND LE>p_le_q
        → 「土壌は乾いているのに蒸散は活発」= 物理的に直接的な深層水アクセス指標

    Returns: df with column 'group_le' ∈ {deep_access, surface_dependent, wet_active, wet_inactive}
    """
    sub = df.dropna(subset=["SWC","LE"]).copy().reset_index(drop=True)
    swc_thr = sub["SWC"].quantile(swc_q)
    le_thr  = sub["LE"].quantile(le_q)

    cond_dry  = sub["SWC"] < swc_thr
    cond_high = sub["LE"]  > le_thr

    sub["group_le"] = "wet_inactive"
    sub.loc[~cond_dry &  cond_high, "group_le"] = "wet_active"
    sub.loc[ cond_dry &  cond_high, "group_le"] = "deep_access"
    sub.loc[ cond_dry & ~cond_high, "group_le"] = "surface_dependent"
    return sub, swc_thr, le_thr


# ================================================================
# 7. [改6] NDWIアノマリ版の3軸分類
# ================================================================

def classify_3axis_anom(df, swc_thr, vpd_thr, ndwi_anom_thr=0.0, use_grace=False):
    sub = df.dropna(subset=["SWC","NDWI_anom","VPD"]).copy().reset_index(drop=True)
    cond_swc_dry  = sub["SWC"] < swc_thr
    cond_ndwi_hi  = sub["NDWI_anom"] > ndwi_anom_thr
    cond_vpd_hi   = sub["VPD"] > vpd_thr

    if use_grace and "GWS_cm" in sub.columns:
        gws_med = sub["GWS_cm"].median()
        cond_deep = cond_ndwi_hi & (sub["GWS_cm"] > gws_med)
        note = f"+GRACE(GWS>{gws_med:.2f}cm)"
    else:
        cond_deep, note = cond_ndwi_hi, ""

    sub["group"] = "wet_shallow"
    sub.loc[~cond_swc_dry &  cond_ndwi_hi,                   "group"] = "wet_deep"
    sub.loc[ cond_swc_dry &  cond_deep    &  cond_vpd_hi,   "group"] = "deep_access"
    sub.loc[ cond_swc_dry &  cond_deep    & ~cond_vpd_hi,   "group"] = "atm_driven"
    sub.loc[ cond_swc_dry & ~cond_ndwi_hi &  cond_vpd_hi,   "group"] = "true_drought"
    sub.loc[ cond_swc_dry & ~cond_ndwi_hi & ~cond_vpd_hi,   "group"] = "compound_dry"
    print(sub["group"].value_counts().to_string())
    return sub, note


# ================================================================
# 8. [改5] SWC→NDWI ラグ相関
# ================================================================

def lag_correlation(df, max_lag=30):
    s = df[["date","SWC","NDWI_s2"]].dropna().sort_values("date").reset_index(drop=True)
    s = s.set_index("date").asfreq("D").interpolate("linear", limit=10)
    lags, rs = [], []
    for L in range(-max_lag, max_lag+1):
        if L < 0:
            r = s["SWC"].corr(s["NDWI_s2"].shift(L))
        else:
            r = s["SWC"].shift(L).corr(s["NDWI_s2"])
        lags.append(L); rs.append(r)
    rs = np.array(rs); lags = np.array(lags)
    best_lag = int(lags[np.nanargmax(rs)])
    return lags, rs, best_lag


# ================================================================
# 9. 閾値決定 (v8 と同じ、簡略版)
# ================================================================

def _kde(s):
    s = s.dropna()
    k = gaussian_kde(s, bw_method="silverman")
    x = np.linspace(s.min(), s.max(), 500)
    return x, k(x), s


def determine_swc_threshold(series, site):
    s = series.dropna()
    p25 = np.percentile(s, 25)
    print(f"  SWC [{site}]: 採用 p25={p25:.3f}%  (n={len(s)})")
    return float(p25)


def determine_vpd_threshold(series, site):
    s = series.dropna()
    p50 = float(np.median(s))
    print(f"  VPD [{site}]: 採用 p50={p50:.2f}hPa  (n={len(s)})")
    return p50


# ================================================================
# 10. 統計検定: [改3] サイト間比較を主分析に
# ================================================================

def cross_site_test(oran_df, tara_df):
    """dry SWC 日における LE/EF を Oran vs Tarazona で比較"""
    print(f"\n{'='*60}\n[★主分析] サイト間比較 — dry SWC 日の LE/EF\n{'='*60}")
    results = {}
    for var in ["LE","EF","ET"]:
        if var not in oran_df.columns or var not in tara_df.columns: continue
        # dry SWC = サイト内 p25 以下
        oq = oran_df["SWC"].quantile(0.25)
        tq = tara_df["SWC"].quantile(0.25)
        d_oran = oran_df[oran_df["SWC"] < oq][var].dropna()
        d_tara = tara_df[tara_df["SWC"] < tq][var].dropna()
        if len(d_oran) < 3 or len(d_tara) < 3:
            print(f"  {var}: サンプル不足"); continue
        u, p = stats.mannwhitneyu(d_tara, d_oran, alternative="greater")
        r = u / (len(d_oran)*len(d_tara))
        results[var] = dict(p=p, r=r,
                            oran_med=d_oran.median(), tara_med=d_tara.median(),
                            n_o=len(d_oran), n_t=len(d_tara))
        sig = "★★★ 有意" if p<0.001 else ("★ 有意" if p<0.05 else "n.s.")
        print(f"  {var:3s}: Oran={d_oran.median():7.2f}(n={len(d_oran)})  "
              f"Tarazona={d_tara.median():7.2f}(n={len(d_tara)})  "
              f"p={p:.2e}  r={r:.3f}  {sig}")
    return results


def within_site_test_le_based(df, site):
    """[改1] LE ベースの群間比較"""
    print(f"\n{'='*60}\n[副分析] LEベース群比較 [{site}]\n{'='*60}")
    g = {k: df[df["group_le"]==k] for k in ["deep_access","surface_dependent",
                                              "wet_active","wet_inactive"]}
    res = {}
    for var in ["LE","EF","ET","NDWI_s2","VPD"]:
        if var not in df.columns: continue
        d1 = g["deep_access"][var].dropna()
        d2 = g["surface_dependent"][var].dropna()
        if len(d1)<3 or len(d2)<3:
            print(f"  {var}: deep={len(d1)}, surf={len(d2)} 不足"); continue
        u, p = stats.mannwhitneyu(d1, d2, alternative="two-sided")
        delta = d1.median() - d2.median()
        res[var] = dict(p=p, delta=delta, deep_med=d1.median(),
                        surf_med=d2.median(), n_d=len(d1), n_s=len(d2))
        sig = "★" if p<0.05 else "n.s."
        print(f"  {var:8s}: deep={d1.median():7.3f}(n={len(d1)})  "
              f"surf_dep={d2.median():7.3f}(n={len(d2)})  Δ={delta:+.3f}  "
              f"p={p:.4f}  {sig}")
    return res


# ================================================================
# 11. 可視化
# ================================================================

G_COL = {"deep_access":"#7C3AED","surface_dependent":"#DC2626",
         "wet_active":"#1D9E75","wet_inactive":"#3B8BD4"}


def plot_le_based_classification(df, swc_thr, le_thr, site, save_dir):
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    fig.patch.set_facecolor("#F8F8F8")

    # Panel1: SWC vs LE 散布図
    ax = axes[0]
    for grp, col in G_COL.items():
        sub = df[df["group_le"]==grp]
        ax.scatter(sub["SWC"], sub["LE"], c=col, alpha=0.6, s=30, label=grp,
                   edgecolors="none")
    ax.axvline(swc_thr, color="#555", lw=1.5, ls="--", label=f"SWC p25={swc_thr:.2f}%")
    ax.axhline(le_thr,  color="#555", lw=1.5, ls=":",  label=f"LE p50={le_thr:.1f}W/m²")
    ax.set_xlabel("SWC (%)"); ax.set_ylabel("LE (W/m²)")
    ax.set_title(f"[{site}] SWC vs LE — LE-based classification")
    ax.legend(fontsize=8); ax.grid(alpha=0.25)

    # Panel2: LE boxplot
    ax = axes[1]
    order = ["surface_dependent","deep_access","wet_inactive","wet_active"]
    data = [df[df["group_le"]==g]["LE"].dropna().values for g in order]
    bp = ax.boxplot(data, patch_artist=True, widths=0.55,
                    medianprops=dict(color="white", lw=2.2))
    for p, g in zip(bp["boxes"], order):
        p.set_facecolor(G_COL[g]); p.set_alpha(0.78)
    ax.set_xticklabels([g.replace("_","\n") for g in order], fontsize=8)
    ax.set_ylabel("LE (W/m²)")
    ax.set_title(f"[{site}] LE by LE-based group")
    ax.grid(axis="y", alpha=0.25)

    # Panel3: サンプル数
    ax = axes[2]
    cnt = df["group_le"].value_counts()
    ax.barh(cnt.index, cnt.values, color=[G_COL.get(g,"#aaa") for g in cnt.index],
            alpha=0.85)
    ax.set_xlabel("日数"); ax.set_title(f"[{site}] Sample counts")
    ax.grid(axis="x", alpha=0.25)

    plt.tight_layout()
    fp = Path(save_dir) / f"v9_le_based_{site}.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [保存] {fp}")


def plot_cross_site(oran_df, tara_df, results, save_dir):
    """[改3] サイト間比較プロット — 主分析"""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    fig.patch.set_facecolor("#F8F8F8")
    fig.suptitle("Cross-site comparison on dry SWC days (SWC < within-site p25)",
                 fontsize=13, fontweight="bold")

    oq = oran_df["SWC"].quantile(0.25)
    tq = tara_df["SWC"].quantile(0.25)
    o_dry = oran_df[oran_df["SWC"] < oq]
    t_dry = tara_df[tara_df["SWC"] < tq]

    for ax, var, ylab in zip(axes, ["LE","EF","ET"],
                              ["LE (W/m²)","EF","ET"]):
        if var not in oran_df.columns or var not in tara_df.columns:
            ax.set_title(f"{var}: 列なし"); continue
        d1 = o_dry[var].dropna(); d2 = t_dry[var].dropna()
        bp = ax.boxplot([d1, d2], patch_artist=True, widths=0.55,
                        medianprops=dict(color="white", lw=2.5),
                        labels=[f"Oran\n(n={len(d1)})", f"Tarazona\n(n={len(d2)})"])
        for p, c in zip(bp["boxes"], ["#E85D04","#1D9E75"]):
            p.set_facecolor(c); p.set_alpha(0.8)
        ax.set_ylabel(ylab); ax.set_title(f"{var} on dry SWC days")
        ax.grid(axis="y", alpha=0.25)
        if var in results:
            r = results[var]
            sig = "★★★" if r["p"]<0.001 else ("★" if r["p"]<0.05 else "n.s.")
            ax.annotate(f"p={r['p']:.2e}  {sig}\n"
                        f"Δ(med)={r['tara_med']-r['oran_med']:+.2f}",
                        xy=(0.5,0.97), xycoords="axes fraction", ha="center",va="top",
                        fontsize=9, color="#7C3AED", fontweight="bold")

    plt.tight_layout()
    fp = Path(save_dir) / "v9_cross_site_dry_SWC.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [保存] {fp}")


def plot_lag_correlation(lag_o, r_o, best_o, lag_t, r_t, best_t, save_dir):
    fig, ax = plt.subplots(figsize=(10, 5.5))
    fig.patch.set_facecolor("#F8F8F8")
    ax.plot(lag_o, r_o, "-o", color="#E85D04", ms=3, label=f"Oran  (peak@{best_o}d)")
    ax.plot(lag_t, r_t, "-o", color="#1D9E75", ms=3, label=f"Tarazona (peak@{best_t}d)")
    ax.axvline(0, color="#555", lw=1, ls="--", alpha=0.6)
    ax.axhline(0, color="#555", lw=0.8, alpha=0.3)
    ax.set_xlabel("Lag (days)  [SWC leads when lag>0]")
    ax.set_ylabel("Pearson r (SWC vs NDWI)")
    ax.set_title("Lag correlation — SWC → NDWI response\n"
                 "(深根なら NDWI の応答が遅れる/弱い)", fontweight="bold")
    ax.legend(); ax.grid(alpha=0.25)
    plt.tight_layout()
    fp = Path(save_dir) / "v9_lag_correlation.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [保存] {fp}")


# ================================================================
# 12. CSV書き出し
# ================================================================

def export_csv(df, site, save_dir):
    df.to_csv(Path(save_dir)/f"v9_{site}_merged.csv", index=False)
    if "group_le" in df.columns:
        cols = [c for c in ["date","SWC","NDWI_s2","NDWI_anom","VPD","LE","EF","ET",
                            "GWS_cm","group_le","group"] if c in df.columns]
        df[df["group_le"]=="deep_access"][cols].to_csv(
            Path(save_dir)/f"v9_{site}_deep_access_LE_based.csv", index=False)
        n = (df["group_le"]=="deep_access").sum()
        print(f"  [保存] LEベース deep_access = {n} 日")


# ================================================================
# 13. サイト実行ラッパー
# ================================================================

def run_site_v9(site_name, ec_df, ndwi_path, lat, lon, swc_thr, vpd_thr):
    print(f"\n{'#'*60}\n# サイト: {site_name}\n{'#'*60}")
    ec_gs  = filter_growing_season(ec_df, site_name)
    ndwi   = load_ndwi_csv(ndwi_path, site_name)
    merged = match_ndwi(ec_gs, ndwi)
    merged = add_ndwi_anomaly(merged)

    grace_used = False
    try:
        grace_df = load_grace(PATHS["grace"], lat, lon)
        merged   = merge_grace(merged, grace_df)
        grace_used = merged["GWS_cm"].notna().sum() > 0
    except Exception as e:
        print(f"  [GRACE] スキップ: {e}")

    # [改1] LEベース分類
    le_classified, swc_dyn, le_dyn = classify_le_based(merged, le_q=0.5, swc_q=0.25)
    print(f"  [LEベース] SWC<{swc_dyn:.2f}% & LE>{le_dyn:.2f}W/m² → "
          f"deep_access={(le_classified['group_le']=='deep_access').sum()} 日")

    # [改6] アノマリ版3軸分類 (補助)
    classified_anom, _ = classify_3axis_anom(le_classified, swc_thr, vpd_thr,
                                              ndwi_anom_thr=0.0, use_grace=grace_used)

    # 検定
    res_le = within_site_test_le_based(classified_anom, site_name)

    # 可視化
    plot_le_based_classification(classified_anom, swc_dyn, le_dyn, site_name, SAVE_DIR)
    export_csv(classified_anom, site_name, SAVE_DIR)
    return classified_anom, res_le


# ================================================================
# MAIN
# ================================================================

if __name__ == "__main__":
    print("="*60)
    print("解析A v9: 深層水アクセス検出 — 再定義版")
    print("="*60)

    # Step1: EC
    oran_ec = load_oran_ec(PATHS["oran_ec"])
    tara_ec = load_tarazona_ec(PATHS["tara_ec"])

    # Step2: SWC単位
    print(f"\n--- Step2: SWC単位統一 ---")
    oran_ec = normalize_swc(oran_ec, "Oran")
    tara_ec = normalize_swc(tara_ec, "Tarazona")

    # Step3: ERA5 VPD
    print(f"\n--- Step3: ERA5 VPD ---")
    era5_o = load_era5_vpd(SITES["Oran"]["lat"],     SITES["Oran"]["lon"])
    era5_t = load_era5_vpd(SITES["Tarazona"]["lat"], SITES["Tarazona"]["lon"])
    oran_ec = merge_era5_vpd(oran_ec, era5_o)
    tara_ec = merge_era5_vpd(tara_ec, era5_t)

    # Step4: 閾値
    print(f"\n--- Step4: 閾値 ---")
    swc_o = determine_swc_threshold(oran_ec["SWC"], "Oran")
    swc_t = determine_swc_threshold(tara_ec["SWC"], "Tarazona")
    vpd_o = determine_vpd_threshold(oran_ec["VPD"], "Oran")
    vpd_t = determine_vpd_threshold(tara_ec["VPD"], "Tarazona")

    # Step5: サイト別解析
    oran_res, oran_stats = run_site_v9("Oran", oran_ec, PATHS["oran_ndwi"],
                                        SITES["Oran"]["lat"], SITES["Oran"]["lon"],
                                        swc_o, vpd_o)
    tara_res, tara_stats = run_site_v9("Tarazona", tara_ec, PATHS["tara_ndwi"],
                                        SITES["Tarazona"]["lat"], SITES["Tarazona"]["lon"],
                                        swc_t, vpd_t)

    # Step6: [★改3] サイト間比較 — 主分析
    cross_results = cross_site_test(oran_res, tara_res)
    plot_cross_site(oran_res, tara_res, cross_results, SAVE_DIR)

    # Step7: [改5] ラグ相関
    print(f"\n--- Step7: SWC→NDWI ラグ相関 ---")
    lag_o, r_o, best_o = lag_correlation(oran_res, max_lag=30)
    lag_t, r_t, best_t = lag_correlation(tara_res, max_lag=30)
    print(f"  Oran     最大相関 lag = {best_o:+d}日  (r_max={np.nanmax(r_o):.3f})")
    print(f"  Tarazona 最大相関 lag = {best_t:+d}日  (r_max={np.nanmax(r_t):.3f})")
    print(f"  → 深根なら lag が大きい / 相関が低い (NDWI が SWC に追従しない)")
    plot_lag_correlation(lag_o, r_o, best_o, lag_t, r_t, best_t, SAVE_DIR)

    # Step8: 最終サマリー
    print(f"\n{'='*60}\n★ 最終サマリー (v9)\n{'='*60}")
    print("\n[1] LEベース deep_access 検出数")
    for name, df_ in [("Oran",oran_res),("Tarazona",tara_res)]:
        n_da = (df_["group_le"]=="deep_access").sum()
        n_sd = (df_["group_le"]=="surface_dependent").sum()
        print(f"   {name:9s}: deep_access={n_da:3d}日   surface_dependent={n_sd:3d}日")

    print("\n[2] 仮説検証 (Tarazona > Oran on dry SWC days)")
    for var in ["LE","EF","ET"]:
        if var in cross_results:
            r = cross_results[var]
            sig = "★★★" if r["p"]<0.001 else ("★" if r["p"]<0.05 else "n.s.")
            print(f"   {var:3s}: Tarazona={r['tara_med']:7.2f}  Oran={r['oran_med']:7.2f}  "
                  f"p={r['p']:.2e}  {sig}")

    print(f"\n[完了] {SAVE_DIR}/ に保存")
