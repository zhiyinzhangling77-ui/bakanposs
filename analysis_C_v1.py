"""
解析C v1: フェノロジー × フラックス
====================================
解析Aで「深根アクセス」仮説は完全には立証されなかった。
解析Cでは植生フェノロジー(NDVI)を独立軸として導入し、
A/B の文脈を補強する:

  [目的1] 生育期/非生育期を NDVI から客観的に定義し、
          v9 で使っていた手動の月フィルター(GROWING_MONTHS)と比較
  [目的2] NDVI vs EF / LE / GPP_proxy の関係を Oran vs Tarazona で比較
  [目的3] サイト間でフェノロジー位相がどう違うかを定量化
          (ピーク NDVI 日、立ち上がり/枯れ落ちの傾き)
  [目的4] [A 連結] 低 SWC × 高 NDVI 期で Tarazona が EF を保てるか
                   = 深根アクセスの状況証拠の再評価
  [目的5] [B 連結] EC LE と NDVI×Rn(簡易 GPP プロキシ)の整合性チェック

入力:
  EC daily : 解析Aと同じ Oran / Tarazona の CSV
  NDVI     : MODIS MOD13Q1 (16-day, 250 m) を CSV で受け取る
             scripts/gee_extract.js を拡張して GEE から書き出す
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from scipy.signal import savgol_filter
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# 解析Aの読み込み関数を再利用 (重複定義を避ける)
from analysis_A_v9 import (
    load_oran_ec, load_tarazona_ec, normalize_swc,
    SITES, GROWING_MONTHS, PATHS as A_PATHS,
)


# ================================================================
# 0. 設定
# ================================================================
BASE_NC = Path("/mnt/hdd/Dataset")
# AppEEARS で抽出した単一 CSV (Oran/TzM 両方を含む)
NDVI_APPEEARS_CSV = BASE_NC / "MOD13Q1_NDVI_EVI/MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv"
SITE_ID_IN_CSV = {"Oran": "Oran", "Tarazona": "TzM"}

SAVE_DIR = Path("./output_analysis_C_v1")
SAVE_DIR.mkdir(exist_ok=True)

SITE_COL = {"Oran": "#E85D04", "Tarazona": "#1D9E75"}


# ================================================================
# 1. NDVI 読み込み・スムージング
# ================================================================

def load_ndvi_csv(filepath, site_name, site_id_in_csv=None,
                   reliability_max=1):
    """AppEEARS の MOD13Q1 NDVI/EVI 抽出 CSV (両サイト同梱) を読む.

    AppEEARS 列名例:
      ID, Latitude, Longitude, Date,
      MOD13Q1_061__250m_16_days_NDVI,
      MOD13Q1_061__250m_16_days_EVI,
      MOD13Q1_061__250m_16_days_pixel_reliability,
      MOD13Q1_061__250m_16_days_VI_Quality, ...

    pixel_reliability: 0=good, 1=marginal, 2=snow/ice, 3=cloudy, -1=fill
    AppEEARS は通常 scale factor 適用済みなので NDVI は [-0.2, 1.0] のはず。
    """
    df = pd.read_csv(filepath)

    # 列名を正規化
    cols = {c.lower(): c for c in df.columns}
    def col(*needles):
        for k, orig in cols.items():
            if all(n in k for n in needles):
                return orig
        return None

    c_id  = col("id") or col("name") or col("category")
    c_dt  = col("date")
    c_ndvi = col("ndvi") and next(c for c in df.columns
                                   if "ndvi" in c.lower()
                                   and "quality" not in c.lower()
                                   and "reliability" not in c.lower())
    c_rel  = col("pixel", "reliability")
    c_qa   = col("vi", "quality")
    if not all([c_id, c_dt, c_ndvi]):
        raise KeyError(f"必須列が不足: {df.columns.tolist()}")

    if site_id_in_csv is None:
        site_id_in_csv = SITE_ID_IN_CSV.get(site_name, site_name)
    sub = df[df[c_id].astype(str).str.strip() == site_id_in_csv].copy()
    if sub.empty:
        ids_seen = sorted(df[c_id].dropna().unique().tolist())
        raise ValueError(f"site_id={site_id_in_csv} が CSV に無い. 候補: {ids_seen}")

    sub = sub.rename(columns={c_dt: "date", c_ndvi: "NDVI"})
    sub["date"] = pd.to_datetime(sub["date"], errors="coerce")
    sub["NDVI"] = pd.to_numeric(sub["NDVI"], errors="coerce")

    # スケールファクタ未適用なら適用
    if sub["NDVI"].abs().max() > 2.0:
        sub["NDVI"] = sub["NDVI"] * 1e-4

    n_raw = len(sub)
    if c_rel and c_rel in sub.columns:
        sub[c_rel] = pd.to_numeric(sub[c_rel], errors="coerce")
        sub = sub[(sub[c_rel] >= 0) & (sub[c_rel] <= reliability_max)]
    sub = sub[sub["NDVI"].between(-0.2, 1.0)]
    sub = sub.dropna(subset=["date", "NDVI"]).sort_values("date").reset_index(drop=True)

    print(f"[NDVI {site_name}] AppEEARS site_id={site_id_in_csv}  "
          f"{len(sub)}/{n_raw} シーン (reliability<={reliability_max})  "
          f"range=[{sub['NDVI'].min():.3f}, {sub['NDVI'].max():.3f}]")
    return sub[["date", "NDVI"]]


def smooth_ndvi(ndvi_df, window=5, poly=2):
    """サブピクセル雲混入を抑えるための Savitzky-Golay."""
    s = ndvi_df.set_index("date")["NDVI"].asfreq("16D").interpolate("linear", limit=3)
    if len(s.dropna()) >= window:
        s_sg = pd.Series(savgol_filter(s.fillna(s.median()).values, window, poly),
                         index=s.index, name="NDVI_sg")
    else:
        s_sg = s.rename("NDVI_sg")
    return s_sg.reset_index()


def match_ndvi(ec_df, ndvi_df, max_days=10):
    """前後最近傍の平均 (NDWI と同パターン)."""
    ec   = ec_df.sort_values("date").reset_index(drop=True).copy()
    ndvi = ndvi_df.sort_values("date").reset_index(drop=True).copy()
    tol = pd.Timedelta(days=max_days)
    fwd = pd.merge_asof(ec, ndvi, on="date", direction="backward", tolerance=tol)
    bwd = pd.merge_asof(ec, ndvi, on="date", direction="forward",  tolerance=tol,
                        suffixes=("", "_bwd"))
    fwd_v = fwd["NDVI"].values
    bwd_c = "NDVI_bwd" if "NDVI_bwd" in bwd.columns else "NDVI"
    bwd_v = bwd[bwd_c].values
    ec["NDVI"] = np.where(np.isnan(fwd_v) & np.isnan(bwd_v), np.nan,
                  np.where(np.isnan(fwd_v), bwd_v,
                  np.where(np.isnan(bwd_v), fwd_v, (fwd_v + bwd_v) / 2)))
    print(f"  NDVI マッチ: {ec['NDVI'].notna().sum()}/{len(ec)} 日")
    return ec


# ================================================================
# 2. NDVI ベース生育期判定
# ================================================================

def ndvi_growing_threshold(ndvi_series, method="otsu"):
    """NDVI 分布から生育期/非生育期の閾値を決める."""
    s = pd.Series(ndvi_series).dropna().values
    if len(s) < 30:
        return float(np.nanmedian(s))
    if method == "otsu":
        # 1次元 Otsu
        hist, edges = np.histogram(s, bins=64)
        prob = hist / hist.sum()
        omega = np.cumsum(prob)
        mu = np.cumsum(prob * 0.5 * (edges[:-1] + edges[1:]))
        mu_t = mu[-1]
        denom = omega * (1 - omega)
        denom[denom == 0] = np.nan
        sigma_b2 = (mu_t * omega - mu) ** 2 / denom
        idx = int(np.nanargmax(sigma_b2))
        thr = float(0.5 * (edges[idx] + edges[idx + 1]))
    else:
        thr = float(np.percentile(s, 33))
    return thr


def classify_phenology(df, ndvi_thr):
    df = df.copy()
    df["phen"] = np.where(df["NDVI"] >= ndvi_thr, "growing", "non_growing")
    return df


def compare_growing_definitions(df, site):
    """月ベース vs NDVI ベースの生育期定義の一致度."""
    months = GROWING_MONTHS[site]
    df = df.copy()
    df["month_growing"]  = df["date"].dt.month.isin(months)
    df["ndvi_growing"]   = df["phen"] == "growing"
    n = len(df.dropna(subset=["NDVI"]))
    if n == 0:
        return None
    agree = (df["month_growing"] == df["ndvi_growing"]).sum()
    only_month = ( df["month_growing"] & ~df["ndvi_growing"]).sum()
    only_ndvi  = (~df["month_growing"] &  df["ndvi_growing"]).sum()
    print(f"  [{site}] 生育期定義の一致 {agree}/{n} ({100*agree/n:.1f}%)  "
          f"month-only={only_month}  ndvi-only={only_ndvi}")
    return dict(n=n, agree=agree, only_month=only_month, only_ndvi=only_ndvi)


# ================================================================
# 3. GPP プロキシと EF
# ================================================================

def add_gpp_proxy(df):
    """GPP_proxy = NDVI × Rn (Monteith 型の簡易版).
    ε × FPAR × PAR の代わり。サイト間の絶対値ではなく相対比較に使う。"""
    df = df.copy()
    if "Rn" in df.columns:
        df["GPP_proxy"] = df["NDVI"] * df["Rn"].clip(lower=0)
    else:
        df["GPP_proxy"] = np.nan
    return df


# ================================================================
# 4. 解析: NDVI と各フラックスの関係
# ================================================================

def correlate_ndvi_flux(df, site):
    out = {}
    print(f"\n[{site}] NDVI vs フラックス Spearman 相関")
    for var in ["LE", "EF", "ET", "GPP_proxy", "H"]:
        if var not in df.columns:
            continue
        sub = df[["NDVI", var]].dropna()
        if len(sub) < 10:
            continue
        r, p = stats.spearmanr(sub["NDVI"], sub[var])
        out[var] = dict(rho=r, p=p, n=len(sub))
        sig = "★★★" if p < 0.001 else ("★" if p < 0.05 else "n.s.")
        print(f"  NDVI ~ {var:9s}: ρ={r:+.3f}  p={p:.2e}  n={len(sub):4d}  {sig}")
    return out


def phenology_phase(ndvi_sg_df, site):
    """ピーク NDVI 日、立ち上がり/枯れ落ち傾き(日あたり)."""
    s = ndvi_sg_df.dropna(subset=["NDVI_sg"]).copy()
    s["doy"] = s["date"].dt.dayofyear
    yearly = []
    for y, g in s.groupby(s["date"].dt.year):
        if len(g) < 12:
            continue
        peak_idx = g["NDVI_sg"].idxmax()
        peak_doy = int(g.loc[peak_idx, "doy"])
        peak_val = float(g.loc[peak_idx, "NDVI_sg"])
        before = g[g["doy"] < peak_doy].tail(8)
        after  = g[g["doy"] > peak_doy].head(8)
        slope_up   = np.polyfit(before["doy"], before["NDVI_sg"], 1)[0] if len(before) >= 3 else np.nan
        slope_down = np.polyfit(after ["doy"], after ["NDVI_sg"], 1)[0] if len(after)  >= 3 else np.nan
        yearly.append((y, peak_doy, peak_val, slope_up, slope_down))
    if not yearly:
        return None
    out = pd.DataFrame(yearly, columns=["year", "peak_doy", "peak_ndvi",
                                        "slope_up", "slope_down"])
    print(f"\n[{site}] フェノロジー位相 (年別)")
    print(out.to_string(index=False, float_format=lambda x: f"{x:7.4f}"))
    return out


# ================================================================
# 5. A 連結: 低 SWC × 高 NDVI 期の Tarazona vs Oran
# ================================================================

def cross_site_dry_canopy(oran, tara):
    """
    深根仮説の状況証拠:
        SWC 低 (サイト内 p25 未満) かつ NDVI 高 (サイト内 p67 以上)
        の日の LE / EF を Oran と Tarazona で比較。
    Tarazona が有意に高ければ、緑の維持を支える追加水源 = 深根アクセス
    の状況証拠になる。
    """
    print(f"\n{'='*60}\n[A 連結] 低 SWC × 高 NDVI 日のサイト間比較\n{'='*60}")

    def filt(df):
        sub = df.dropna(subset=["SWC", "NDVI"])
        sw = sub["SWC"].quantile(0.25)
        nv = sub["NDVI"].quantile(0.67)
        return sub[(sub["SWC"] < sw) & (sub["NDVI"] >= nv)]

    do = filt(oran)
    dt = filt(tara)
    print(f"  Oran     n={len(do)}    Tarazona n={len(dt)}")
    res = {}
    for var in ["LE", "EF", "ET"]:
        if var not in oran.columns or var not in tara.columns:
            continue
        a = do[var].dropna()
        b = dt[var].dropna()
        if len(a) < 3 or len(b) < 3:
            print(f"  {var}: サンプル不足 (Oran={len(a)} Tarazona={len(b)})")
            continue
        u, p = stats.mannwhitneyu(b, a, alternative="greater")
        r = u / (len(a) * len(b))
        res[var] = dict(p=p, r=r, oran_med=float(a.median()),
                        tara_med=float(b.median()), n_o=len(a), n_t=len(b))
        sig = "★★★" if p < 0.001 else ("★" if p < 0.05 else "n.s.")
        print(f"  {var:3s}: Oran={a.median():7.2f}(n={len(a)})  "
              f"Tarazona={b.median():7.2f}(n={len(b)})  p={p:.2e}  {sig}")
    return res


# ================================================================
# 6. 可視化
# ================================================================

def plot_ndvi_timeseries(per_site, save_dir):
    fig, ax = plt.subplots(figsize=(12, 5))
    fig.patch.set_facecolor("#F8F8F8")
    for site, (df, sg) in per_site.items():
        ax.plot(sg["date"], sg["NDVI_sg"], "-",
                color=SITE_COL[site], lw=2, label=f"{site} (SG)")
        raw = df.dropna(subset=["NDVI"])
        ax.scatter(raw["date"], raw["NDVI"], s=10, alpha=0.35,
                   color=SITE_COL[site], edgecolors="none")
    ax.set_xlabel("Date"); ax.set_ylabel("NDVI")
    ax.set_title("MOD13Q1 NDVI time series (16-day, smoothed)", fontweight="bold")
    ax.grid(alpha=0.25); ax.legend()
    plt.tight_layout()
    fp = Path(save_dir) / "C_ndvi_timeseries.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [保存] {fp}")


def plot_ndvi_vs_flux(per_site, save_dir):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.patch.set_facecolor("#F8F8F8")
    for ax, var, ylab in zip(axes, ["EF", "LE", "GPP_proxy"],
                              ["EF", "LE (W/m²)", "GPP proxy = NDVI×Rn"]):
        for site, (df, _) in per_site.items():
            sub = df[["NDVI", var]].dropna()
            if len(sub) < 5:
                continue
            ax.scatter(sub["NDVI"], sub[var], s=12, alpha=0.4,
                       color=SITE_COL[site], edgecolors="none", label=site)
            # 1次回帰
            try:
                m, c = np.polyfit(sub["NDVI"], sub[var], 1)
                xs = np.linspace(sub["NDVI"].min(), sub["NDVI"].max(), 100)
                ax.plot(xs, m * xs + c, color=SITE_COL[site], lw=1.5)
            except Exception:
                pass
        ax.set_xlabel("NDVI"); ax.set_ylabel(ylab)
        ax.set_title(f"NDVI vs {var}")
        ax.grid(alpha=0.25); ax.legend(fontsize=9)
    plt.tight_layout()
    fp = Path(save_dir) / "C_ndvi_vs_flux.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [保存] {fp}")


def plot_growing_definition(df, site, save_dir):
    fig, ax = plt.subplots(figsize=(12, 4.5))
    fig.patch.set_facecolor("#F8F8F8")
    sub = df.dropna(subset=["NDVI"]).copy()
    ax.scatter(sub["date"], sub["NDVI"], c=np.where(sub["phen"]=="growing", "#1D9E75", "#888"),
               s=14, alpha=0.7, edgecolors="none")
    months = GROWING_MONTHS[site]
    in_m = sub["date"].dt.month.isin(months)
    for d in sub.loc[in_m, "date"]:
        ax.axvline(d, color="#3B8BD4", alpha=0.05, lw=0.5)
    ax.set_title(f"[{site}] NDVI-based growing-season classification "
                 f"(green=growing) vs month filter (blue band)")
    ax.set_ylabel("NDVI"); ax.grid(alpha=0.25)
    plt.tight_layout()
    fp = Path(save_dir) / f"C_growing_def_{site}.png"
    plt.savefig(fp, dpi=150, bbox_inches="tight", facecolor="#F8F8F8")
    plt.close()
    print(f"  [保存] {fp}")


# ================================================================
# 7. サイト実行
# ================================================================

def run_site_C(site, ec_df, ndvi_path):
    print(f"\n{'#'*60}\n# 解析C: {site}\n{'#'*60}")
    ndvi_df = load_ndvi_csv(ndvi_path, site)
    sg_df   = smooth_ndvi(ndvi_df)

    merged = match_ndvi(ec_df, ndvi_df)
    merged = add_gpp_proxy(merged)

    thr = ndvi_growing_threshold(merged["NDVI"].values, method="otsu")
    print(f"  [{site}] NDVI 生育期閾値 (Otsu) = {thr:.3f}")
    merged = classify_phenology(merged, thr)
    compare_growing_definitions(merged, site)

    corr = correlate_ndvi_flux(merged, site)
    phase = phenology_phase(sg_df, site)

    out_csv = SAVE_DIR / f"C_{site}_merged.csv"
    merged.to_csv(out_csv, index=False)
    print(f"  [保存] {out_csv}")

    plot_growing_definition(merged, site, SAVE_DIR)
    return merged, sg_df, dict(thr=thr, corr=corr, phase=phase)


# ================================================================
# MAIN
# ================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("解析C v1: フェノロジー × フラックス")
    print("=" * 60)

    oran_ec = normalize_swc(load_oran_ec(A_PATHS["oran_ec"]), "Oran")
    tara_ec = normalize_swc(load_tarazona_ec(A_PATHS["tara_ec"]), "Tarazona")

    oran_m, oran_sg, oran_meta = run_site_C("Oran",     oran_ec, NDVI_APPEEARS_CSV)
    tara_m, tara_sg, tara_meta = run_site_C("Tarazona", tara_ec, NDVI_APPEEARS_CSV)

    per_site = {"Oran": (oran_m, oran_sg), "Tarazona": (tara_m, tara_sg)}
    plot_ndvi_timeseries(per_site, SAVE_DIR)
    plot_ndvi_vs_flux(per_site, SAVE_DIR)

    a_results = cross_site_dry_canopy(oran_m, tara_m)

    print(f"\n{'='*60}\n★ 解析C 最終サマリー\n{'='*60}")
    for site, meta in [("Oran", oran_meta), ("Tarazona", tara_meta)]:
        if meta["phase"] is not None:
            mp = meta["phase"]["peak_doy"].mean()
            print(f"  {site:9s}: NDVI 閾値={meta['thr']:.3f}  "
                  f"平均ピーク DOY={mp:.0f}")

    print("\n[A 連結] 低SWC × 高NDVI 日の LE/EF (Tarazona > Oran か?)")
    for var in ["LE", "EF", "ET"]:
        if var in a_results:
            r = a_results[var]
            sig = "★★★" if r["p"] < 0.001 else ("★" if r["p"] < 0.05 else "n.s.")
            print(f"   {var:3s}: Oran={r['oran_med']:7.2f}  "
                  f"Tarazona={r['tara_med']:7.2f}  p={r['p']:.2e}  {sig}")

    print(f"\n[完了] {SAVE_DIR}/ に保存")
