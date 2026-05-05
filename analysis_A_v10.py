"""
解析A v10: 深層水アクセス検出 — 交絡対策版

v9 の根本的問題への対応
=========================
[F1] LE閾値を絶対値化 + アクティブ植被日のみ採用
        旧: LE > サイト内p50  → Oranでは2.49 W/m² (= 植物がある日かないかの判別)
        新: LE > LE_ABS_MIN(=50 W/m²)  かつ  NDWI > NDWI_ACTIVE_MIN
            → 「実質的に蒸散している日」のみで議論する

[F2] サイト間比較を "drought response ratio" に変更
        旧: LE_dry を直接比較 (Oran=3.4 vs Tarazona=195)
            → 季節・作物・LAI差の交絡で因果性なし
        新: LE_dry / LE_wet の比 を比較
            → 各サイトが「自分のbaselineに対してどれだけ落ちたか」=
              本物の表層SWC感受性を比較できる

[F3] ラグ相関を「30日移動平均からの偏差(アノマリ)」で計算
        旧: 季節周期に支配されてフラット → 解釈不能
        新: 短期偏差で見る → 真のSWC→NDWI 応答ラグを抽出

[F4] 負LE/異常値フィルタ
        新: LE < -20 W/m² の日を除外 (ノイズ・夜間反転の典型値)

[F5] 月別 phenologically-matched 比較
        新: 同月内で Oran と Tarazona の dry SWC 日 LE を比較
            → 「春季限定」など phenology を揃えた比較を追加

[F6] Tarazona 内 deep_access 群比較を主結果に格上げ
        v9 で示された「LE↑ EF↑ NDWI↑ VPD↑」が真の主張
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
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

# [F1] 物理的に意味のあるしきい値
LE_ABS_MIN       = 50.0    # W/m²   実質的な蒸散があるとみなす最小LE
LE_NEG_FILTER    = -20.0   # W/m²   この値より低いLEは異常値として除外 [F4]
NDWI_ACTIVE_MIN  = -0.05   # 植被が活動中とみなす最小NDWI(サイト依存性低減のため緩い値)
EF_DENOM_MIN     = 10
ERA5_YEARS       = list(range(2018, 2025))

SAVE_DIR = Path("./output_analysis_A_v10")
SAVE_DIR.mkdir(exist_ok=True)


# ================================================================
# 1. ERA5 VPD (v9 と同じ)
# ================================================================

def magnus_e(t):
    return 6.112 * np.exp(17.67 * t / (t + 243.5))


def load_era5_vpd(lat, lon, years=ERA5_YEARS, era5_dir=ERA5_DIR):
    records = []
    t_vars  = ["t2m", "T2M", "temperature_2m", "t"]
    d_vars  = ["d2m", "D2M", "dewpoint_2m", "td"]
    rh_vars = ["r", "RH", "relative_humidity", "rh"]
    for year in years:
        nc = era5_dir / f"{year}.nc"
        if not nc.exists(): continue
        try:
            ds = xr.open_dataset(nc, engine="netcdf4")
            drop = [c for c in ds.coords if c not in
                    ("valid_time","time","latitude","longitude","lat","lon")]
            ds = ds.drop_vars(drop, errors="ignore")
            def get(cands):
                for n in cands:
                    if n in ds: return ds[n]
                raise KeyError(cands)
            def sel(da):
                ln = next(c for c in da.dims if "lat" in c.lower())
                lo = next(c for c in da.dims if "lon" in c.lower())
                return da.sel({ln:lat, lo:lon}, method="nearest")
            def to_d(da, name):
                s = sel(da).to_series(); s.index = pd.to_datetime(s.index)
                d_ = s.reset_index(); d_.columns = ["time", name]; return d_
            t = to_d(get(t_vars), "t2m")
            if t["t2m"].median() > 100: t["t2m"] -= 273.15
            try:
                d = to_d(get(d_vars), "d2m")
                if d["d2m"].median() > 100: d["d2m"] -= 273.15
                m = pd.merge(t, d, on="time")
                m["VPD"] = np.maximum(magnus_e(m["t2m"]) - magnus_e(m["d2m"]), 0)
            except KeyError:
                rh = to_d(get(rh_vars), "rh")
                if rh["rh"].median() < 2: rh["rh"] *= 100
                m = pd.merge(t, rh, on="time")
                es = magnus_e(m["t2m"])
                m["VPD"] = np.maximum(es - es * m["rh"]/100, 0)
            m["date"] = m["time"].dt.normalize()
            d_ = m.groupby("date")["VPD"].mean().reset_index()
            d_.columns = ["date", "VPD_era5"]
            records.append(d_); ds.close()
        except Exception as e:
            print(f"  [ERA5] {year}: {e}")
    return pd.concat(records, ignore_index=True).sort_values("date")


def merge_era5(ec, era5):
    m = ec.copy()
    if "VPD" in m.columns: m = m.rename(columns={"VPD":"VPD_ec"})
    m = pd.merge(m, era5, on="date", how="left").rename(columns={"VPD_era5":"VPD"})
    return m


# ================================================================
# 2. EC 読み込み
# ================================================================

def load_oran_ec(fp):
    df = pd.read_csv(fp)
    dt = next((c for c in df.columns if c.lower() in ("datetime","timestamp")), None)
    df["datetime"] = pd.to_datetime(df[dt], errors="coerce")
    df = df.dropna(subset=["datetime"]).sort_values("datetime").reset_index(drop=True)
    for c in ["SWC_1_1_1","LE","H","G","NETRAD","VPD","ET"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["date"] = df["datetime"].dt.normalize()
    daily = df.groupby("date", as_index=False).agg(
        SWC=("SWC_1_1_1","mean"), LE=("LE","mean"), H=("H","mean"),
        G=("G","mean"), Rn=("NETRAD","mean"), VPD=("VPD","mean"), ET=("ET","sum"))
    denom = daily["Rn"] - daily["G"]; v = denom > EF_DENOM_MIN
    daily["EF"] = np.nan
    daily.loc[v,"EF"] = (daily.loc[v,"LE"]/denom[v]).clip(0,1.5)
    daily = daily[(daily["SWC"]>0) & (daily["SWC"]<100)].reset_index(drop=True)
    daily["site"] = "Oran"
    print(f"[Oran EC] {len(daily)} 日分")
    return daily


def load_tarazona_ec(fp):
    df = pd.read_csv(fp)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df = df.dropna(subset=["date"]).sort_values("date").reset_index(drop=True)
    df = df.rename(columns={"SWC_avg":"SWC","LE_avg":"LE","H_avg":"H","G_avg":"G",
                             "NetRad_avg":"Rn","VPD_mean":"VPD","ET_avg":"ET"})
    for c in ["SWC","LE","H","G","Rn","VPD","ET"]:
        if c in df.columns: df[c] = pd.to_numeric(df[c], errors="coerce")
    denom = df["Rn"] - df["G"]; v = denom > EF_DENOM_MIN
    df["EF"] = np.nan
    df.loc[v,"EF"] = (df.loc[v,"LE"]/denom[v]).clip(0,1.5)
    df = df[(df["SWC"]>0) & (df["SWC"]<100)].reset_index(drop=True)
    df["site"] = "Tarazona"
    print(f"[Tarazona EC] {len(df)} 日分")
    return df[[c for c in ["date","site","SWC","LE","H","G","Rn","EF","VPD","ET"]
               if c in df.columns]]


def normalize_swc(df, site):
    if df["SWC"].dropna().max() <= 1.0:
        df = df.copy(); df["SWC"] *= 100
        print(f"  SWC [{site}]: m³/m³ → %")
    return df


# ================================================================
# 3. NDWI / GRACE
# ================================================================

def load_ndwi(fp, site):
    df = pd.read_csv(fp)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df = df.rename(columns={"NDWI":"NDWI_s2"})[["date","NDWI_s2"]].dropna()
    df = df.sort_values("date").reset_index(drop=True)
    print(f"[NDWI {site}] {len(df)} シーン")
    return df


def match_ndwi(ec, ndwi, max_days=5):
    ec = ec.sort_values("date").reset_index(drop=True).copy()
    ndwi = ndwi.sort_values("date").reset_index(drop=True).copy()
    tol = pd.Timedelta(days=max_days)
    fwd = pd.merge_asof(ec, ndwi, on="date", direction="backward", tolerance=tol)
    bwd = pd.merge_asof(ec, ndwi, on="date", direction="forward",  tolerance=tol,
                        suffixes=("","_bwd"))
    fv = fwd["NDWI_s2"].values
    bc = "NDWI_s2_bwd" if "NDWI_s2_bwd" in bwd.columns else "NDWI_s2"
    bv = bwd[bc].values
    ec["NDWI_s2"] = np.where(np.isnan(fv) & np.isnan(bv), np.nan,
                     np.where(np.isnan(fv), bv,
                     np.where(np.isnan(bv), fv, (fv+bv)/2)))
    print(f"  NDWI マッチ: {ec['NDWI_s2'].notna().sum()}/{len(ec)}")
    return ec


def add_ndwi_anomaly(df):
    df = df.copy()
    df["month"] = df["date"].dt.month
    med = df.groupby("month")["NDWI_s2"].transform("median")
    df["NDWI_anom"] = df["NDWI_s2"] - med
    return df.drop(columns="month")


def load_grace(fp, lat, lon):
    ds = xr.open_dataset(fp)
    if ds["lon"].values.min() >= 0 and lon < 0:
        lon = 360 + lon
    lwe = ds["lwe_thickness"].sel(lat=lat, lon=lon, method="nearest")
    sf  = float(ds["scale_factor"].sel(lat=lat, lon=lon, method="nearest").values)
    df  = (lwe*sf).to_dataframe(name="GWS_cm").reset_index().rename(columns={"time":"date"})
    df["date"] = pd.to_datetime(df["date"]).astype("datetime64[ns]")
    return df[["date","GWS_cm"]].dropna()


def merge_grace(ec, grace):
    ec = ec.copy(); ec["date"] = ec["date"].astype("datetime64[ns]")
    g = grace.set_index("date").resample("D").interpolate("linear").reset_index()
    return pd.merge_asof(ec.sort_values("date"), g.sort_values("date"),
                         on="date", direction="nearest", tolerance=pd.Timedelta("32D"))


def filter_growing_season(df, site):
    out = df[df["date"].dt.month.isin(GROWING_MONTHS[site])].copy().reset_index(drop=True)
    print(f"  [生育期] {site}: {len(df)} → {len(out)}")
    return out


# ================================================================
# 4. [F4] LE/NDWI フィルタ
# ================================================================

def filter_quality(df, site):
    """LE 異常値除外 + アクティブ植被日のみ採用"""
    n0 = len(df)
    df = df[df["LE"].notna()].copy()
    df = df[df["LE"] > LE_NEG_FILTER]
    n1 = len(df)
    if "NDWI_s2" in df.columns:
        df = df[df["NDWI_s2"].fillna(-99) > NDWI_ACTIVE_MIN]
    n2 = len(df)
    print(f"  [品質F {site}] {n0} → {n1} (LE>{LE_NEG_FILTER}) → {n2} (NDWI>{NDWI_ACTIVE_MIN})")
    return df.reset_index(drop=True)


# ================================================================
# 5. [F1] 改良版 LE-based 分類
# ================================================================

def classify_le_absolute(df, swc_q=0.25, le_min=LE_ABS_MIN):
    """deep_access = SWC<p25 AND LE>LE_ABS_MIN(絶対値)"""
    sub = df.dropna(subset=["SWC","LE"]).copy().reset_index(drop=True)
    swc_thr = sub["SWC"].quantile(swc_q)
    cond_dry  = sub["SWC"] < swc_thr
    cond_high = sub["LE"]  > le_min

    sub["group_le"] = "wet_inactive"
    sub.loc[~cond_dry &  cond_high, "group_le"] = "wet_active"
    sub.loc[ cond_dry &  cond_high, "group_le"] = "deep_access"
    sub.loc[ cond_dry & ~cond_high, "group_le"] = "surface_dependent"
    return sub, swc_thr, le_min


# ================================================================
# 6. [F2] Drought Response Ratio
# ================================================================

def drought_response_ratio(df, site, swc_q=0.25, le_min=LE_ABS_MIN, n_boot=2000):
    """
    DRR = LE_dry_median / LE_wet_median
    DRR≈1 なら "土壌乾燥に非依存"(深根サイン)
    DRR≪1 なら "土壌乾燥で激減"(浅根)
    ブートストラップで95%CIも返す
    """
    sub = df.dropna(subset=["SWC","LE"]).copy()
    sub = sub[sub["LE"] > le_min]  # アクティブ日のみ
    if len(sub) < 20:
        print(f"  [DRR {site}] アクティブ日不足 (n={len(sub)})")
        return None
    swc_thr = sub["SWC"].quantile(swc_q)
    le_dry = sub.loc[sub["SWC"] <  swc_thr, "LE"].dropna()
    le_wet = sub.loc[sub["SWC"] >= swc_thr, "LE"].dropna()
    if len(le_dry) < 5 or len(le_wet) < 5:
        print(f"  [DRR {site}] 群サンプル不足 dry={len(le_dry)} wet={len(le_wet)}")
        return None
    drr_point = le_dry.median() / le_wet.median()
    rng = np.random.default_rng(42)
    boots = []
    for _ in range(n_boot):
        b_d = rng.choice(le_dry.values, size=len(le_dry), replace=True)
        b_w = rng.choice(le_wet.values, size=len(le_wet), replace=True)
        if np.median(b_w) > 0:
            boots.append(np.median(b_d) / np.median(b_w))
    boots = np.array(boots)
    ci_lo, ci_hi = np.percentile(boots, [2.5, 97.5])
    print(f"  [DRR {site}] {drr_point:.3f}  [95%CI {ci_lo:.3f}, {ci_hi:.3f}]  "
          f"(dry n={len(le_dry)}, wet n={len(le_wet)})")
    return dict(drr=drr_point, ci_lo=ci_lo, ci_hi=ci_hi,
                le_dry=le_dry.median(), le_wet=le_wet.median(),
                n_dry=len(le_dry), n_wet=len(le_wet))


# ================================================================
# 7. [F3] アノマリベースのラグ相関
# ================================================================

def lag_correlation_anomaly(df, max_lag=30, smooth_window=30):
    """30日移動平均からの偏差で計算 → 季節周期を除去"""
    s = df[["date","SWC","NDWI_s2"]].dropna().sort_values("date")
    s = s.set_index("date").asfreq("D").interpolate("linear", limit=10)
    sw = s["SWC"]    - s["SWC"].rolling(smooth_window, min_periods=10, center=True).mean()
    nd = s["NDWI_s2"] - s["NDWI_s2"].rolling(smooth_window, min_periods=10, center=True).mean()
    lags, rs = [], []
    for L in range(-max_lag, max_lag+1):
        if L < 0:
            r = sw.corr(nd.shift(L))
        else:
            r = sw.shift(L).corr(nd)
        lags.append(L); rs.append(r)
    rs = np.array(rs); lags = np.array(lags)
    best = int(lags[np.nanargmax(np.abs(rs))])
    return lags, rs, best


# ================================================================
# 8. [F5] Phenology-matched 月別比較
# ================================================================

def monthly_matched_comparison(oran, tara, swc_q=0.25):
    """両サイトに共通する月だけで dry SWC 日 LE を比較"""
    print(f"\n{'='*60}\n[★F5] Phenology-matched 月別比較\n{'='*60}")
    common = sorted(set(GROWING_MONTHS["Oran"]) & set(GROWING_MONTHS["Tarazona"]))
    print(f"  共通月: {common}")
    rows = []
    for m in common:
        o = oran[(oran["date"].dt.month == m) & oran["LE"].notna() &
                 (oran["LE"] > LE_NEG_FILTER)]
        t = tara[(tara["date"].dt.month == m) & tara["LE"].notna() &
                 (tara["LE"] > LE_NEG_FILTER)]
        if len(o) < 10 or len(t) < 10: continue
        o_q = o["SWC"].quantile(swc_q); t_q = t["SWC"].quantile(swc_q)
        o_dry = o[o["SWC"] < o_q]["LE"].dropna()
        t_dry = t[t["SWC"] < t_q]["LE"].dropna()
        if len(o_dry) < 3 or len(t_dry) < 3: continue
        try:
            _, p = stats.mannwhitneyu(t_dry, o_dry, alternative="greater")
        except ValueError:
            p = np.nan
        rows.append(dict(month=m, n_o=len(o_dry), n_t=len(t_dry),
                         le_o=o_dry.median(), le_t=t_dry.median(), p=p,
                         ratio=t_dry.median()/max(o_dry.median(), 1e-3)))
    if not rows:
        print("  比較可能な月なし"); return pd.DataFrame()
    res = pd.DataFrame(rows)
    print(res.to_string(index=False))
    return res


# ================================================================
# 9. 統計検定
# ================================================================

def within_site_test(df, site):
    print(f"\n{'='*60}\n[主分析] LEベース群比較 [{site}]\n{'='*60}")
    g = {k: df[df["group_le"]==k] for k in
         ["deep_access","surface_dependent","wet_active","wet_inactive"]}
    res = {}
    for var in ["LE","EF","ET","NDWI_s2","NDWI_anom","VPD"]:
        if var not in df.columns: continue
        d1 = g["deep_access"][var].dropna()
        d2 = g["surface_dependent"][var].dropna()
        if len(d1)<3 or len(d2)<3:
            print(f"  {var:9s}: n不足 deep={len(d1)} surf={len(d2)}"); continue
        u, p = stats.mannwhitneyu(d1, d2, alternative="two-sided")
        res[var] = dict(p=p, deep_med=d1.median(), surf_med=d2.median(),
                        n_d=len(d1), n_s=len(d2))
        sig = "★★★" if p<0.001 else ("★" if p<0.05 else "n.s.")
        print(f"  {var:9s}: deep={d1.median():8.3f}(n={len(d1)})  "
              f"surf_dep={d2.median():8.3f}(n={len(d2)})  "
              f"Δ={d1.median()-d2.median():+.3f}  p={p:.2e}  {sig}")
    return res


# ================================================================
# 10. 可視化
# ================================================================

G_COL = {"deep_access":"#7C3AED","surface_dependent":"#DC2626",
         "wet_active":"#1D9E75","wet_inactive":"#3B8BD4"}


def plot_drr(drr_o, drr_t, save_dir):
    if drr_o is None or drr_t is None: return
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor("#F8F8F8")
    sites  = ["Oran\n(cereal\nshallow)", "Tarazona\n(almond\ndeep)"]
    vals   = [drr_o["drr"], drr_t["drr"]]
    los    = [drr_o["ci_lo"], drr_t["ci_lo"]]
    his    = [drr_o["ci_hi"], drr_t["ci_hi"]]
    err_lo = [v-l for v,l in zip(vals, los)]
    err_hi = [h-v for h,v in zip(his, vals)]
    bars = ax.bar(sites, vals, yerr=[err_lo, err_hi], capsize=12,
                   color=["#E85D04","#1D9E75"], alpha=0.85, edgecolor="black", lw=1.2)
    ax.axhline(1.0, color="#555", lw=1.2, ls="--",
               label="DRR=1.0 (no drought sensitivity = 深根仮説)")
    ax.set_ylabel("Drought Response Ratio  (LE_dry / LE_wet)", fontsize=11)
    ax.set_title("★ Drought Response Ratio (95% CI bootstrap)\n"
                 "高いほど表層SWC低下に対して LE を維持 = 深根",
                 fontsize=12, fontweight="bold")
    ax.set_ylim(0, max(1.4, max(his)*1.15))
    for b, v in zip(bars, vals):
        ax.text(b.get_x()+b.get_width()/2, v+0.03, f"{v:.2f}",
                ha="center", fontsize=12, fontweight="bold")
    ax.legend(fontsize=10); ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    fp = Path(save_dir) / "v10_DRR.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close(); print(f"  [保存] {fp}")


def plot_within_site(df, site, swc_thr, le_thr, results, save_dir):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.patch.set_facecolor("#F8F8F8")
    fig.suptitle(f"[{site}] LE-based deep_access — within-site comparison",
                 fontsize=13, fontweight="bold")

    # SWC vs LE 散布図
    ax = axes[0,0]
    for grp, col in G_COL.items():
        sub = df[df["group_le"]==grp]
        ax.scatter(sub["SWC"], sub["LE"], c=col, alpha=0.6, s=30,
                   label=f"{grp} (n={len(sub)})", edgecolors="none")
    ax.axvline(swc_thr, color="#555", lw=1.5, ls="--",
               label=f"SWC p25={swc_thr:.2f}%")
    ax.axhline(le_thr,  color="#555", lw=1.5, ls=":",
               label=f"LE_min={le_thr:.0f} W/m² (絶対値)")
    ax.set_xlabel("SWC (%)"); ax.set_ylabel("LE (W/m²)")
    ax.set_title("SWC vs LE")
    ax.legend(fontsize=8); ax.grid(alpha=0.25)

    # NDWI vs VPD 散布(deep_accessで右上に集まれば仮説支持)
    ax = axes[0,1]
    for grp in ["deep_access","surface_dependent","wet_active","wet_inactive"]:
        sub = df[df["group_le"]==grp].dropna(subset=["NDWI_s2","VPD"])
        ax.scatter(sub["VPD"], sub["NDWI_s2"], c=G_COL[grp], alpha=0.6, s=30,
                   label=grp, edgecolors="none")
    ax.set_xlabel("VPD (hPa)"); ax.set_ylabel("NDWI (S2)")
    ax.set_title("NDWI vs VPD\n→ deep_access が右上 = 仮説支持")
    ax.legend(fontsize=8); ax.grid(alpha=0.25)

    # 4変数 boxplot
    for ax_, var, ylab in [(axes[1,0], "EF", "EF"), (axes[1,1], "NDWI_s2", "NDWI")]:
        order = ["surface_dependent","deep_access","wet_inactive","wet_active"]
        data = [df[df["group_le"]==g][var].dropna().values for g in order]
        bp = ax_.boxplot(data, patch_artist=True, widths=0.55,
                         medianprops=dict(color="white", lw=2.2))
        for p, g in zip(bp["boxes"], order):
            p.set_facecolor(G_COL[g]); p.set_alpha(0.78)
        ax_.set_xticklabels([g.replace("_","\n") for g in order], fontsize=8)
        ax_.set_ylabel(ylab); ax_.set_title(f"{var} by group")
        ax_.grid(axis="y", alpha=0.25)
        if var in results:
            r = results[var]
            sig = "★★★" if r["p"]<0.001 else ("★" if r["p"]<0.05 else "n.s.")
            ax_.annotate(f"deep vs surf: p={r['p']:.2e} {sig}",
                         xy=(0.5,0.97), xycoords="axes fraction",
                         ha="center", va="top", fontsize=9,
                         color="#7C3AED", fontweight="bold")

    plt.tight_layout()
    fp = Path(save_dir) / f"v10_within_{site}.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close(); print(f"  [保存] {fp}")


def plot_lag_anomaly(lag_o, r_o, b_o, lag_t, r_t, b_t, save_dir):
    fig, ax = plt.subplots(figsize=(10, 5.5))
    fig.patch.set_facecolor("#F8F8F8")
    ax.plot(lag_o, r_o, "-o", color="#E85D04", ms=3,
            label=f"Oran  (peak@{b_o:+d}d, |r|={np.nanmax(np.abs(r_o)):.3f})")
    ax.plot(lag_t, r_t, "-o", color="#1D9E75", ms=3,
            label=f"Tarazona (peak@{b_t:+d}d, |r|={np.nanmax(np.abs(r_t)):.3f})")
    ax.axvline(0, color="#555", lw=1, ls="--", alpha=0.6)
    ax.axhline(0, color="#555", lw=0.8, alpha=0.3)
    ax.set_xlabel("Lag (days)  [SWC leads NDWI when lag>0]")
    ax.set_ylabel("Pearson r  (anomaly, 30d-detrended)")
    ax.set_title("[F3] Lag correlation on ANOMALIES (deseasonalized)\n"
                 "深根なら lag>0 で正のピーク・かつピーク位置が大きい",
                 fontweight="bold")
    ax.legend(); ax.grid(alpha=0.25)
    plt.tight_layout()
    fp = Path(save_dir) / "v10_lag_anomaly.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close(); print(f"  [保存] {fp}")


def plot_monthly(monthly_df, save_dir):
    if monthly_df is None or monthly_df.empty: return
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.patch.set_facecolor("#F8F8F8")
    ax = axes[0]
    x = np.arange(len(monthly_df))
    w = 0.35
    ax.bar(x-w/2, monthly_df["le_o"], w, color="#E85D04", alpha=0.85,
           label="Oran (cereal)")
    ax.bar(x+w/2, monthly_df["le_t"], w, color="#1D9E75", alpha=0.85,
           label="Tarazona (almond)")
    ax.set_xticks(x); ax.set_xticklabels([f"M{m}" for m in monthly_df["month"]])
    ax.set_ylabel("LE on dry SWC days (W/m²)")
    ax.set_title("[F5] Monthly phenology-matched comparison",
                 fontweight="bold")
    ax.legend(); ax.grid(axis="y", alpha=0.25)
    ax = axes[1]
    ax.bar(x, monthly_df["ratio"], color="#7C3AED", alpha=0.85)
    ax.axhline(1.0, color="#555", lw=1, ls="--", label="ratio=1")
    ax.set_xticks(x); ax.set_xticklabels([f"M{m}" for m in monthly_df["month"]])
    ax.set_ylabel("LE_Tarazona / LE_Oran  (dry SWC days)")
    ax.set_title("Tarazona / Oran ratio")
    ax.legend(); ax.grid(axis="y", alpha=0.25)
    plt.tight_layout()
    fp = Path(save_dir) / "v10_monthly_matched.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close(); print(f"  [保存] {fp}")


# ================================================================
# 11. CSV
# ================================================================

def export_csv(df, site, save_dir):
    df.to_csv(Path(save_dir)/f"v10_{site}_merged.csv", index=False)
    if "group_le" in df.columns:
        for grp in ["deep_access","surface_dependent","wet_active","wet_inactive"]:
            sub = df[df["group_le"]==grp]
            if len(sub):
                sub.to_csv(Path(save_dir)/f"v10_{site}_{grp}.csv", index=False)


# ================================================================
# 12. サイト実行
# ================================================================

def run_site(site, ec, ndwi_path, lat, lon):
    print(f"\n{'#'*60}\n# {site}\n{'#'*60}")
    ec_gs = filter_growing_season(ec, site)
    ndwi  = load_ndwi(ndwi_path, site)
    m = match_ndwi(ec_gs, ndwi)
    m = add_ndwi_anomaly(m)
    try:
        g = load_grace(PATHS["grace"], lat, lon)
        m = merge_grace(m, g)
    except Exception as e:
        print(f"  GRACE skip: {e}")
    # [F4] フィルタ
    m_f = filter_quality(m, site)
    # [F1] 絶対LEで分類
    classified, swc_thr, le_thr = classify_le_absolute(m_f)
    print(f"  [LE-abs] SWC<{swc_thr:.2f}% & LE>{le_thr:.0f}W/m² → "
          f"deep_access={(classified['group_le']=='deep_access').sum()} 日")
    # 検定
    res = within_site_test(classified, site)
    plot_within_site(classified, site, swc_thr, le_thr, res, SAVE_DIR)
    export_csv(classified, site, SAVE_DIR)
    return m_f, classified, res


# ================================================================
# MAIN
# ================================================================

if __name__ == "__main__":
    print("="*60)
    print("解析A v10: 深層水アクセス検出 — 交絡対策版")
    print("="*60)

    # Step1-3
    oran_ec = load_oran_ec(PATHS["oran_ec"])
    tara_ec = load_tarazona_ec(PATHS["tara_ec"])
    print("\n--- SWC 単位統一 ---")
    oran_ec = normalize_swc(oran_ec, "Oran")
    tara_ec = normalize_swc(tara_ec, "Tarazona")
    print("\n--- ERA5 VPD ---")
    e_o = load_era5_vpd(SITES["Oran"]["lat"], SITES["Oran"]["lon"])
    e_t = load_era5_vpd(SITES["Tarazona"]["lat"], SITES["Tarazona"]["lon"])
    oran_ec = merge_era5(oran_ec, e_o)
    tara_ec = merge_era5(tara_ec, e_t)

    # Step4: サイト解析
    o_full, o_cls, o_res = run_site("Oran", oran_ec, PATHS["oran_ndwi"],
                                     SITES["Oran"]["lat"], SITES["Oran"]["lon"])
    t_full, t_cls, t_res = run_site("Tarazona", tara_ec, PATHS["tara_ndwi"],
                                     SITES["Tarazona"]["lat"], SITES["Tarazona"]["lon"])

    # Step5: [F2] DRR
    print(f"\n{'='*60}\n[★F2] Drought Response Ratio (主結果)\n{'='*60}")
    drr_o = drought_response_ratio(o_full, "Oran")
    drr_t = drought_response_ratio(t_full, "Tarazona")
    plot_drr(drr_o, drr_t, SAVE_DIR)

    # Step6: [F3] アノマリラグ相関
    print(f"\n{'='*60}\n[★F3] Anomaly-based Lag correlation\n{'='*60}")
    lag_o, r_o, b_o = lag_correlation_anomaly(o_full)
    lag_t, r_t, b_t = lag_correlation_anomaly(t_full)
    print(f"  Oran     peak lag={b_o:+d}日  |r|max={np.nanmax(np.abs(r_o)):.3f}")
    print(f"  Tarazona peak lag={b_t:+d}日  |r|max={np.nanmax(np.abs(r_t)):.3f}")
    plot_lag_anomaly(lag_o, r_o, b_o, lag_t, r_t, b_t, SAVE_DIR)

    # Step7: [F5] phenology-matched
    monthly = monthly_matched_comparison(o_full, t_full)
    if not monthly.empty:
        monthly.to_csv(SAVE_DIR/"v10_monthly_matched.csv", index=False)
        plot_monthly(monthly, SAVE_DIR)

    # Step8: 最終サマリー
    print(f"\n{'='*60}\n★★★ 最終サマリー (v10) ★★★\n{'='*60}")
    print("\n[1] Tarazona 内 deep_access vs surface_dependent (★主結果)")
    if "LE" in t_res:
        r = t_res["LE"]
        sig = "★★★" if r["p"]<0.001 else ("★" if r["p"]<0.05 else "n.s.")
        print(f"   LE      : deep={r['deep_med']:.1f}  surf_dep={r['surf_med']:.1f}  "
              f"p={r['p']:.2e}  {sig}")
    if "EF" in t_res:
        r = t_res["EF"]
        sig = "★★★" if r["p"]<0.001 else ("★" if r["p"]<0.05 else "n.s.")
        print(f"   EF      : deep={r['deep_med']:.3f} surf_dep={r['surf_med']:.3f} "
              f"p={r['p']:.2e}  {sig}")
    if "NDWI_s2" in t_res:
        r = t_res["NDWI_s2"]
        sig = "★★★" if r["p"]<0.001 else ("★" if r["p"]<0.05 else "n.s.")
        print(f"   NDWI    : deep={r['deep_med']:.3f} surf_dep={r['surf_med']:.3f} "
              f"p={r['p']:.2e}  {sig}")

    print("\n[2] Drought Response Ratio (LE_dry/LE_wet, 95%CI)")
    if drr_o:
        print(f"   Oran     : {drr_o['drr']:.3f}  [{drr_o['ci_lo']:.3f}, {drr_o['ci_hi']:.3f}]"
              f"  → {'高感受' if drr_o['drr']<0.7 else '中間' if drr_o['drr']<0.9 else '非感受(深根)'}")
    if drr_t:
        print(f"   Tarazona : {drr_t['drr']:.3f}  [{drr_t['ci_lo']:.3f}, {drr_t['ci_hi']:.3f}]"
              f"  → {'高感受' if drr_t['drr']<0.7 else '中間' if drr_t['drr']<0.9 else '非感受(深根)'}")
    if drr_o and drr_t:
        if drr_t["drr"] > drr_o["drr"]:
            print(f"   ⇒ Tarazona の DRR が高い = 表層SWC低下に対する LE 維持力が強い"
                  f" = 深根仮説を支持")

    print("\n[3] アノマリラグ相関 (季節除去後)")
    print(f"   Oran     peak lag={b_o:+d}日  |r|max={np.nanmax(np.abs(r_o)):.3f}")
    print(f"   Tarazona peak lag={b_t:+d}日  |r|max={np.nanmax(np.abs(r_t)):.3f}")
    print(f"   ⇒ 浅根なら lag≈0 付近で正の鋭いピーク, 深根なら遅延 or 弱相関")

    print(f"\n[完了] {SAVE_DIR}/")
