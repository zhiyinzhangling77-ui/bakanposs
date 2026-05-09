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

# 解析Aの設定値・補助関数を再利用 (重複定義を避ける)
from analysis_A_v9 import (
    load_oran_ec, load_tarazona_ec,
    SITES, GROWING_MONTHS, PATHS as A_PATHS,
)
# クリーンローダ (バグ修正版) は共通モジュールから
from data_loaders import (
    load_oran_ec_clean, load_tarazona_ec_clean, normalize_swc,
    EF_DENOM_MIN, SENTINEL_THR,
)


def diagnose_oran_parse_failures(path):
    """Oran CSV のうち TIMESTAMP がパース不能だった行の中身を覗く."""
    df = pd.read_csv(path)
    ts = df["TIMESTAMP"].astype(str).str.strip()
    dt = pd.to_datetime(ts, format="%Y/%m/%d %H:%M:%S", errors="coerce")
    miss = dt.isna()
    if miss.any():
        dt2 = pd.to_datetime(ts[miss], format="%Y/%m/%d", errors="coerce")
        dt = dt.where(~miss, dt2)
    failed = dt.isna()
    n_fail = int(failed.sum())
    print(f"\n[Oran parse 失敗診断] {n_fail}/{len(ts)} 行 ({100 * n_fail / len(ts):.1f}%)")
    if n_fail == 0:
        return
    uniq = ts[failed].drop_duplicates().tolist()
    print(f"  ユニーク失敗値: {len(uniq)} 種")
    for s in uniq[:20]:
        print(f"    {repr(s)}")


def inspect_tarazona_units(path):
    """Tarazona の ET_avg vs ET_sum, VPD_mean vs VPD_kPa を比較し
    単位を判定."""
    df = pd.read_csv(path)
    print("\n[Tarazona 単位検証]")
    if {"ET_avg", "ET_sum"}.issubset(df.columns):
        a = pd.to_numeric(df["ET_avg"], errors="coerce")
        s = pd.to_numeric(df["ET_sum"], errors="coerce")
        ok = a.notna() & s.notna() & (a > 0)
        if ok.sum() > 10:
            ratio = (s[ok] / a[ok]).median()
            print(f"  ET_avg med={a.median():.3f}  ET_sum med={s.median():.3f}  "
                  f"sum/avg ratio≈{ratio:.1f} (24=hourly, 48=halfhourly)")
    if {"VPD_mean", "VPD_kPa"}.issubset(df.columns):
        m = pd.to_numeric(df["VPD_mean"], errors="coerce")
        k = pd.to_numeric(df["VPD_kPa"], errors="coerce")
        ok = m.notna() & k.notna() & (k > 0)
        if ok.sum() > 10:
            ratio = (m[ok] / k[ok]).median()
            print(f"  VPD_mean med={m.median():.1f}  VPD_kPa med={k.median():.3f}  "
                  f"mean/kPa ratio≈{ratio:.0f} (1000=Pa, 10=hPa)")


# ================================================================
# 0. 設定
#    (loader 本体は data_loaders.py に切り出した。
#     A/B でも同じローダを使うため共有モジュール化)
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
    c_evi = next((c for c in df.columns
                  if "evi" in c.lower()
                  and "quality" not in c.lower()
                  and "reliability" not in c.lower()), None)
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

    rename_map = {c_dt: "date", c_ndvi: "NDVI"}
    if c_evi:
        rename_map[c_evi] = "EVI"
    sub = sub.rename(columns=rename_map)
    sub["date"] = pd.to_datetime(sub["date"], errors="coerce")
    sub["NDVI"] = pd.to_numeric(sub["NDVI"], errors="coerce")
    if "EVI" in sub.columns:
        sub["EVI"] = pd.to_numeric(sub["EVI"], errors="coerce")

    if sub["NDVI"].abs().max() > 2.0:
        sub["NDVI"] = sub["NDVI"] * 1e-4
    if "EVI" in sub.columns and sub["EVI"].abs().max() > 2.0:
        sub["EVI"] = sub["EVI"] * 1e-4

    n_raw = len(sub)
    if c_rel and c_rel in sub.columns:
        sub[c_rel] = pd.to_numeric(sub[c_rel], errors="coerce")
        sub = sub[(sub[c_rel] >= 0) & (sub[c_rel] <= reliability_max)]
    sub = sub[sub["NDVI"].between(-0.2, 1.0)]
    sub = sub.dropna(subset=["date", "NDVI"]).sort_values("date").reset_index(drop=True)

    print(f"[NDVI {site_name}] AppEEARS site_id={site_id_in_csv}  "
          f"{len(sub)}/{n_raw} シーン (reliability<={reliability_max})  "
          f"NDVI=[{sub['NDVI'].min():.3f}, {sub['NDVI'].max():.3f}]"
          + (f"  EVI=[{sub['EVI'].min():.3f}, {sub['EVI'].max():.3f}]"
             if "EVI" in sub.columns else ""))
    keep = ["date", "NDVI"] + (["EVI"] if "EVI" in sub.columns else [])
    return sub[keep]


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
    """前後最近傍の平均で NDVI (および EVI) を結合."""
    ec   = ec_df.sort_values("date").reset_index(drop=True).copy()
    ndvi = ndvi_df.sort_values("date").reset_index(drop=True).copy()
    tol = pd.Timedelta(days=max_days)

    vi_cols = [c for c in ["NDVI", "EVI"] if c in ndvi.columns]
    fwd = pd.merge_asof(ec, ndvi, on="date", direction="backward", tolerance=tol)
    bwd = pd.merge_asof(ec, ndvi, on="date", direction="forward",  tolerance=tol,
                        suffixes=("", "_bwd"))
    for col in vi_cols:
        fwd_v = fwd[col].values
        bwd_c = f"{col}_bwd" if f"{col}_bwd" in bwd.columns else col
        bwd_v = bwd[bwd_c].values
        ec[col] = np.where(np.isnan(fwd_v) & np.isnan(bwd_v), np.nan,
                    np.where(np.isnan(fwd_v), bwd_v,
                    np.where(np.isnan(bwd_v), fwd_v, (fwd_v + bwd_v) / 2)))
    print(f"  NDVI マッチ: {ec['NDVI'].notna().sum()}/{len(ec)} 日"
          + (f"  EVI マッチ: {ec['EVI'].notna().sum()}" if "EVI" in ec.columns else ""))
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

def cross_site_dry_canopy(oran, tara, mode="ndvi_p67"):
    """
    深根仮説の状況証拠:
        SWC 低 (サイト内 p25 未満) かつ 「植生が活発」な日の LE / EF を
        Oran と Tarazona で比較。
    Tarazona が有意に高ければ、緑の維持を支える追加水源 = 深根アクセスの状況証拠。

    mode:
      'ndvi_p67'     : NDVI ≥ p67  (元仕様, 厳しい)
      'ndvi_p50'     : NDVI ≥ p50  (緩和)
      'growing_flag' : phen=='growing'  (Otsu 由来)
    """
    print(f"\n{'='*60}\n[A 連結 / mode={mode}] 低 SWC × 緑 のサイト間比較\n{'='*60}")

    def filt(df):
        sub = df.dropna(subset=["SWC", "NDVI"])
        sw = sub["SWC"].quantile(0.25)
        if mode == "growing_flag":
            return sub[(sub["SWC"] < sw) & (sub.get("phen") == "growing")]
        q = 0.50 if mode == "ndvi_p50" else 0.67
        nv = sub["NDVI"].quantile(q)
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
# 5b. 診断: EC データの実態を可視化
# ================================================================

def diagnose_ec(df, site):
    """Rn / G / denom / EF などの欠損・分布を吐く."""
    print(f"\n[EC診断 {site}]  rows={len(df)}")
    print(f"  columns: {list(df.columns)}")
    for c in ["LE", "H", "G", "Rn", "SWC", "VPD", "ET", "EF", "NDVI"]:
        if c in df.columns:
            s = pd.to_numeric(df[c], errors="coerce")
            n = int(s.notna().sum())
            if n == 0:
                print(f"  {c:5s}: all NaN")
                continue
            print(f"  {c:5s}: n={n:4d}/{len(df)}  "
                  f"min={s.min():>8.3f}  med={s.median():>8.3f}  "
                  f"p75={s.quantile(0.75):>8.3f}  max={s.max():>8.3f}")
    if all(c in df.columns for c in ["Rn", "G"]):
        denom = df["Rn"] - df["G"]
        thr = 10
        print(f"  denom=Rn-G : n_valid(>{thr})={(denom > thr).sum()}/{len(df)}  "
              f"med={denom.median():.2f}  p25={denom.quantile(0.25):.2f}  "
              f"p75={denom.quantile(0.75):.2f}")


def inspect_raw_oran_csv(path):
    """v9 が読む前の生 CSV を覗く → EF 列名や G 列名のズレを発見."""
    print(f"\n[Raw Oran CSV inspect] {path}")
    head = pd.read_csv(path, nrows=3)
    print(f"  shape head: {head.shape}")
    print(f"  columns ({len(head.columns)}): {list(head.columns)}")


def dump_oran_first_rows(path, n=10):
    """先頭 n 行で DateTime と TIMESTAMP の値を確認 — どの列が半時間刻みかを目視."""
    df = pd.read_csv(path, nrows=n)
    cols = [c for c in ["DateTime", "TIMESTAMP", "year", "Julian",
                        "Time_hours", "Time_days"] if c in df.columns]
    print(f"\n[Oran 最初 {n} 行 timestamp 関連列]")
    for c in cols:
        vals = df[c].tolist()
        print(f"  {c:11s}: {vals}")


def deep_inspect_oran_raw(path):
    """
    Oran 生 CSV を半時間 / 日次レベルで徹底診断:
      - 行頻度 (半時間? 1時間? 日?)
      - NETRAD の真の分布 (符号反転? 単位ズレ?)
      - 日中(SW_IN>100)時の Rn 符号
      - SW_IN/OUT, LW_IN/OUT から Rn を再計算した値と比較
    """
    print(f"\n[Oran 生 CSV 徹底診断]")
    df = pd.read_csv(path)
    print(f"  total rows: {len(df)}")
    df["dt"] = pd.to_datetime(df["DateTime"], errors="coerce")
    df = df.dropna(subset=["dt"]).sort_values("dt").reset_index(drop=True)
    df["date"] = df["dt"].dt.date
    rpd = df.groupby("date").size()
    print(f"  date range: {df['dt'].min()} ~ {df['dt'].max()}")
    print(f"  unique days: {df['date'].nunique()}")
    print(f"  rows/day: median={rpd.median():.0f}  min={rpd.min()}  max={rpd.max()}")
    delta = df["dt"].diff().dropna().mode()
    if len(delta):
        print(f"  most common Δt: {delta.iloc[0]}")

    nr = pd.to_numeric(df["NETRAD"], errors="coerce")
    print(f"\n  NETRAD raw (n={nr.notna().sum()}/{len(nr)}):")
    print(f"    min={nr.min():.2f}  p1={nr.quantile(0.01):.2f}  "
          f"p25={nr.quantile(0.25):.2f}  med={nr.median():.2f}  "
          f"p75={nr.quantile(0.75):.2f}  p99={nr.quantile(0.99):.2f}  max={nr.max():.2f}")
    print(f"    n(==-9999)={int((nr == -9999).sum())}  "
          f"n(<=-500)={int((nr <= -500).sum())}  "
          f"n(==0)={int((nr == 0).sum())}  "
          f"n(>0)={int((nr > 0).sum())}  "
          f"n(>200)={int((nr > 200).sum())}")

    if "SW_IN" in df.columns:
        sw = pd.to_numeric(df["SW_IN"], errors="coerce")
        print(f"\n  SW_IN raw: med={sw.median():.2f}  "
              f"p99={sw.quantile(0.99):.2f}  max={sw.max():.2f}  "
              f"n(>100)={int((sw > 100).sum())}")
        daytime = (sw > 100) & nr.notna()
        if daytime.sum() > 0:
            nr_day = nr[daytime]
            print(f"    daytime(SW_IN>100, n={daytime.sum()}) NETRAD: "
                  f"med={nr_day.median():.2f}  "
                  f"p25={nr_day.quantile(0.25):.2f}  "
                  f"p75={nr_day.quantile(0.75):.2f}  "
                  f"frac<0={(nr_day < 0).mean() * 100:.1f}%")

    cols_lw = ["SW_IN", "SW_OUT", "LW_IN", "LW_OUT"]
    if all(c in df.columns for c in cols_lw):
        for c in cols_lw:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        rn_calc = (df["SW_IN"] - df["SW_OUT"]) + (df["LW_IN"] - df["LW_OUT"])
        print(f"\n  Rn 再計算 (SW_IN-SW_OUT+LW_IN-LW_OUT):")
        print(f"    min={rn_calc.min():.2f}  med={rn_calc.median():.2f}  "
              f"max={rn_calc.max():.2f}  n_valid={rn_calc.notna().sum()}")
        ok = rn_calc.notna() & nr.notna()
        if ok.sum() > 100:
            r = np.corrcoef(rn_calc[ok], nr[ok])[0, 1]
            ratio = (nr[ok] / rn_calc[ok].replace(0, np.nan)).median()
            diff = (nr[ok] - rn_calc[ok]).median()
            print(f"    corr(NETRAD, Rn_calc) = {r:+.3f}  "
                  f"median(NETRAD/Rn_calc) = {ratio:+.3f}  "
                  f"median(NETRAD-Rn_calc) = {diff:+.2f}")

    for c in ["LE", "H", "G", "ET", "VPD"]:
        if c in df.columns:
            s = pd.to_numeric(df[c], errors="coerce")
            print(f"  {c:5s} raw: med={s.median():>8.3f}  "
                  f"p25={s.quantile(0.25):>8.3f}  p75={s.quantile(0.75):>8.3f}  "
                  f"min={s.min():>8.3f}  max={s.max():>8.3f}")


# ================================================================
# 5c. Tarazona の強い NDVI~GPP_proxy 信号を深掘り
# ================================================================

def partial_correlation_analysis(df, site):
    """
    LE = a + b·NDVI + c·Rn の標準化重回帰 + 部分相関で
    NDVI と Rn の独立寄与を切り分ける.
    """
    g = df[df.get("phen") == "growing"].dropna(
        subset=["NDVI", "Rn", "LE"]).copy()
    if len(g) < 30:
        print(f"\n[Partial corr {site}] n<30, skip")
        return
    le = g["LE"].values
    nv = g["NDVI"].values
    rn = g["Rn"].values
    z_le = (le - le.mean()) / le.std()
    z_nv = (nv - nv.mean()) / nv.std()
    z_rn = (rn - rn.mean()) / rn.std()

    X = np.column_stack([np.ones(len(z_le)), z_nv, z_rn])
    beta, *_ = np.linalg.lstsq(X, z_le, rcond=None)
    pred = X @ beta
    r2 = 1.0 - ((z_le - pred) ** 2).sum() / ((z_le - z_le.mean()) ** 2).sum()

    def partial_r(y, x, z):
        ry = y - z * (np.dot(z, y) / np.dot(z, z))
        rx = x - z * (np.dot(z, x) / np.dot(z, z))
        return float(np.corrcoef(ry, rx)[0, 1])

    pr_nv = partial_r(z_le, z_nv, z_rn)
    pr_rn = partial_r(z_le, z_rn, z_nv)

    print(f"\n[Partial corr {site}] 生育期 n={len(g)}")
    print(f"  標準化回帰  z(LE) = {beta[0]:+.3f} "
          f"+ ({beta[1]:+.3f})·z(NDVI) + ({beta[2]:+.3f})·z(Rn)")
    print(f"  R²            = {r2:.3f}")
    print(f"  partial r(LE, NDVI | Rn)   = {pr_nv:+.3f}  "
          f"(NDVI の独立寄与)")
    print(f"  partial r(LE, Rn   | NDVI) = {pr_rn:+.3f}  "
          f"(Rn の独立寄与)")
    if abs(pr_rn) > abs(pr_nv) * 1.5:
        print(f"  → Rn 主導: {site} の LE は放射駆動")
    elif abs(pr_nv) > abs(pr_rn) * 1.5:
        print(f"  → NDVI 主導: {site} の LE は植生主導")
    else:
        print(f"  → NDVI と Rn が同等寄与")
    return dict(beta=beta.tolist(), r2=r2, pr_nv=pr_nv, pr_rn=pr_rn, n=len(g))


def _sig_stars(p):
    if not np.isfinite(p): return "n.s."
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "n.s."


PAPER_RC = {
    "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
    "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.9,
}
SITE_COLOR = {"Oran": "#c75b3c", "Tarazona": "#2f7a7a"}


def plot_final_summary(oran, tara, save_dir):
    """
    解析C の主結果を 1 枚で:
      左 : 低 SWC × 高 NDVI (p50) 期の EF 分布 (box) — 深根仮説の主証拠
      右 : 生育期 NDVI~LE を VPD 三分位で層別したサイト比較
    """
    save_dir = Path(save_dir)

    def filt(df):
        sub = df.dropna(subset=["SWC", "NDVI"])
        return sub[(sub["SWC"] < sub["SWC"].quantile(0.25)) &
                   (sub["NDVI"] >= sub["NDVI"].quantile(0.50))]

    oran_dry = filt(oran)
    tara_dry = filt(tara)

    with plt.rc_context(PAPER_RC):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5),
                                 gridspec_kw=dict(wspace=0.3))

        ax = axes[0]
        data, labels, ns = [], [], []
        for name, sub in [("Oran", oran_dry), ("Tarazona", tara_dry)]:
            v = sub["EF"].dropna().values
            if len(v):
                data.append(v); labels.append(name); ns.append(len(v))
        bp = ax.boxplot(data, labels=labels, patch_artist=True, widths=0.55,
                        showmeans=True, showfliers=False,
                        medianprops=dict(color="black", lw=1.4),
                        meanprops=dict(marker="D", markerfacecolor="white",
                                       markeredgecolor="black", markersize=6))
        for patch, name in zip(bp["boxes"], labels):
            patch.set_facecolor(SITE_COLOR[name]); patch.set_alpha(0.55)
            patch.set_edgecolor(SITE_COLOR[name])
        # data points overlay
        for i, (v, name) in enumerate(zip(data, labels)):
            jitter = (np.random.RandomState(0).rand(len(v)) - 0.5) * 0.15
            ax.scatter(np.full_like(v, i + 1, dtype=float) + jitter, v,
                       s=14, color=SITE_COLOR[name], alpha=0.55,
                       edgecolors="white", linewidths=0.4, zorder=3)

        for i, n in enumerate(ns):
            ax.text(i + 1, -0.03, f"n={n}", ha="center", va="top", fontsize=9)

        ax.set_ylabel("Evaporative Fraction (EF)")
        ax.set_title("(a)  EF under low SWC × high NDVI",
                     loc="left", fontweight="bold")
        ax.grid(True, alpha=0.25, axis="y")
        ax.set_ylim(-0.08, 1.55)

        if len(data) == 2 and min(len(data[0]), len(data[1])) >= 3:
            _, p = stats.mannwhitneyu(data[1], data[0], alternative="greater")
            star = _sig_stars(p)
            y_top = 1.4
            ax.plot([1, 1, 2, 2], [y_top - 0.04, y_top, y_top, y_top - 0.04],
                    color="black", lw=1.0)
            ax.text(1.5, y_top + 0.02, f"{star}  p = {p:.1e}",
                    ha="center", va="bottom", fontsize=10)

        ax = axes[1]
        qlabels = ["low", "mid", "high"]
        bar = {}
        for site, df in [("Oran", oran), ("Tarazona", tara)]:
            sub = df[df.get("phen") == "growing"].dropna(
                subset=["NDVI", "LE", "VPD"]).copy()
            rho = [np.nan] * 3
            if len(sub) >= 30:
                sub["q"] = pd.qcut(sub["VPD"], q=3, labels=qlabels)
                for j, q in enumerate(qlabels):
                    ss = sub[sub["q"] == q]
                    if len(ss) >= 10:
                        rho[j], _ = stats.spearmanr(ss["NDVI"], ss["LE"])
            bar[site] = rho

        x = np.arange(3); width = 0.36
        ax.bar(x - width / 2, bar["Oran"],     width, label="Oran",
               color=SITE_COLOR["Oran"],     alpha=0.85, edgecolor="white")
        ax.bar(x + width / 2, bar["Tarazona"], width, label="Tarazona",
               color=SITE_COLOR["Tarazona"], alpha=0.85, edgecolor="white")
        for j, q in enumerate(qlabels):
            for k, site in enumerate(["Oran", "Tarazona"]):
                v = bar[site][j]
                if np.isfinite(v):
                    off = -width / 2 if site == "Oran" else width / 2
                    ax.text(j + off, v + (0.03 if v >= 0 else -0.05),
                            f"{v:+.2f}", ha="center", fontsize=9)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels(qlabels)
        ax.set_xlabel("VPD tertile (within site)")
        ax.set_ylabel(r"Spearman $\rho$ (NDVI, LE)")
        ax.set_title("(b)  NDVI–LE coupling vs VPD stress",
                     loc="left", fontweight="bold")
        ax.legend(frameon=False, loc="upper right")
        ax.grid(True, alpha=0.25, axis="y")
        ax.set_ylim(-0.45, 0.9)

        out = save_dir / "C_final_summary.png"
        plt.savefig(out, dpi=180, bbox_inches="tight")
        plt.close()
        print(f"  [保存] {out}")


def plot_monthly_ndvi(per_site, save_dir):
    """月別 NDVI 中央値 ± IQR をサイトで重ねる (フェノロジー位相比較)."""
    save_dir = Path(save_dir)
    with plt.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for site, (df, _) in per_site.items():
            sub = df.dropna(subset=["NDVI"]).copy()
            if "date" not in sub.columns or len(sub) == 0:
                continue
            sub["month"] = pd.to_datetime(sub["date"]).dt.month
            grp = sub.groupby("month")["NDVI"]
            med = grp.median()
            q25 = grp.quantile(0.25)
            q75 = grp.quantile(0.75)
            x = med.index.values
            ax.fill_between(x, q25.values, q75.values,
                            color=SITE_COLOR[site], alpha=0.20)
            ax.plot(x, med.values, "-o", color=SITE_COLOR[site],
                    label=site, lw=2, markersize=5)
        ax.set_xticks(range(1, 13))
        ax.set_xlabel("Month"); ax.set_ylabel("NDVI")
        ax.set_title("Monthly NDVI (median ± IQR)",
                     loc="left", fontweight="bold")
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.25)
        out = save_dir / "C_monthly_ndvi.png"
        plt.savefig(out, dpi=180, bbox_inches="tight")
        plt.close()
        print(f"  [保存] {out}")


def plot_partial_corr(oran, tara, save_dir):
    """生育期内の partial correlation を 2 サイト並置."""
    save_dir = Path(save_dir)
    rows = []
    for site, df in [("Oran", oran), ("Tarazona", tara)]:
        g = df[df.get("phen") == "growing"].dropna(
            subset=["NDVI", "Rn", "LE"]).copy()
        if len(g) < 30:
            continue
        z = lambda s: (s - s.mean()) / s.std()
        z_le = z(g["LE"].values); z_nv = z(g["NDVI"].values); z_rn = z(g["Rn"].values)

        def part(y, x, w):
            ry = y - w * (np.dot(w, y) / np.dot(w, w))
            rx = x - w * (np.dot(w, x) / np.dot(w, w))
            return float(np.corrcoef(ry, rx)[0, 1])
        rows.append((site, part(z_le, z_nv, z_rn), part(z_le, z_rn, z_nv)))

    if not rows:
        return
    with plt.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(7, 4.2))
        sites = [r[0] for r in rows]
        nv_vals = [r[1] for r in rows]
        rn_vals = [r[2] for r in rows]
        x = np.arange(len(sites)); width = 0.36
        ax.bar(x - width / 2, nv_vals, width, label=r"NDVI | Rn",
               color="#6c8ebf", edgecolor="white")
        ax.bar(x + width / 2, rn_vals, width, label=r"Rn | NDVI",
               color="#d6a64a", edgecolor="white")
        for i, (n, r) in enumerate(zip(nv_vals, rn_vals)):
            ax.text(i - width / 2, n + 0.02 * np.sign(n or 1), f"{n:+.2f}",
                    ha="center", fontsize=10)
            ax.text(i + width / 2, r + 0.02 * np.sign(r or 1), f"{r:+.2f}",
                    ha="center", fontsize=10)
        ax.set_xticks(x); ax.set_xticklabels(sites)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_ylabel("Partial correlation with LE")
        ax.set_title("Independent contribution of NDVI vs Rn to LE\n(growing season)",
                     loc="left", fontweight="bold")
        ax.legend(frameon=False, loc="upper left")
        ax.grid(True, alpha=0.25, axis="y")
        out = save_dir / "C_partial_corr.png"
        plt.savefig(out, dpi=180, bbox_inches="tight")
        plt.close()
        print(f"  [保存] {out}")


# ================================================================
# H2: NDVI が高 LAI 域で飽和していないか EVI で検証
#
#   仮説: Tarazona アーモンド樹冠は LAI が高くても NDVI は 0.5–0.6 で
#         頭打ち。生育期内 NDVI 分散が情報を持たない可能性。
#   検証: NDVI と EVI の関係が線形なら飽和なし。NDVI が頭打ちで EVI が
#         伸びるなら飽和あり (= EVI を使う方が植生活性度を捉えられる)。
# ================================================================
def ndvi_saturation_check(merged, site, save_dir):
    if "EVI" not in merged.columns or merged["EVI"].notna().sum() < 30:
        print(f"\n[H2 {site}] EVI 列なし or サンプル不足; skip"); return None

    print(f"\n{'=' * 60}\n[H2 {site}] NDVI vs EVI 飽和検証\n{'=' * 60}")
    sub = merged.dropna(subset=["NDVI", "EVI"]).copy()
    if len(sub) < 30:
        return None

    r_lin, _ = stats.spearmanr(sub["NDVI"], sub["EVI"])
    print(f"  全データ ρ(NDVI, EVI) = {r_lin:+.3f}  n={len(sub)}")

    # 高 NDVI 領域 (>= p67) と低-中 NDVI 領域 (<= p67) で線形関係が変わるか
    p67 = sub["NDVI"].quantile(0.67)
    low = sub[sub["NDVI"] < p67]
    hi  = sub[sub["NDVI"] >= p67]
    if len(hi) >= 10 and len(low) >= 10:
        r_low, _ = stats.spearmanr(low["NDVI"], low["EVI"])
        r_hi,  _ = stats.spearmanr(hi["NDVI"],  hi["EVI"])
        slope_low = np.polyfit(low["NDVI"], low["EVI"], 1)[0]
        slope_hi  = np.polyfit(hi["NDVI"],  hi["EVI"], 1)[0]
        print(f"  NDVI < p67 (low-mid, n={len(low)}): "
              f"ρ={r_low:+.3f}  slope EVI/NDVI={slope_low:+.3f}")
        print(f"  NDVI ≥ p67 (high,    n={len(hi)}):  "
              f"ρ={r_hi:+.3f}  slope EVI/NDVI={slope_hi:+.3f}")
        if slope_hi < slope_low * 0.6:
            print(f"  → 高 NDVI 領域で傾きが落ちている = NDVI 飽和の兆候")
        else:
            print(f"  → 飽和の明確な兆候なし")

    # EVI を使うと NDVI で見えなかった LE 相関が出るか
    if "LE" in merged.columns:
        gg = merged[merged.get("phen") == "growing"].dropna(
            subset=["NDVI", "EVI", "LE"])
        if len(gg) >= 30:
            r_n, _ = stats.spearmanr(gg["NDVI"], gg["LE"])
            r_e, _ = stats.spearmanr(gg["EVI"],  gg["LE"])
            print(f"  生育期内 (n={len(gg)}): ρ(NDVI,LE)={r_n:+.3f}  "
                  f"ρ(EVI,LE)={r_e:+.3f}")

    # Figure
    with plt.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(5.5, 5))
        ax.scatter(sub["NDVI"], sub["EVI"], s=10, alpha=0.45,
                   color=SITE_COLOR[site], edgecolors="white", linewidths=0.3)
        m, b = np.polyfit(sub["NDVI"], sub["EVI"], 1)
        xs = np.linspace(sub["NDVI"].min(), sub["NDVI"].max(), 50)
        ax.plot(xs, m * xs + b, "k--", lw=1.0,
                label=f"linear fit: y={m:+.2f}x{b:+.3f}")
        ax.plot(xs, xs, ":", color="gray", lw=0.8, label="y = x")
        ax.set_xlabel("NDVI"); ax.set_ylabel("EVI")
        ax.set_title(f"{site}: NDVI vs EVI", loc="left", fontweight="bold")
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.25)
        out = Path(save_dir) / f"C_H2_ndvi_evi_{site}.png"
        plt.savefig(out, dpi=180, bbox_inches="tight"); plt.close()
        print(f"  [保存] {out}")


# ================================================================
# H4: Oran の NDVI~H 負相関がアルベド経由か検証
#
#   仮説: NDVI 上昇 → アルベド低下 → 利用可能エネルギー増 → H 減?
#         それとも気孔開放 → LE 増 → エネルギーバランス上 H 減?
#   検証: Oran CSV の ALB 列を使って、
#     1) NDVI ~ ALB が負相関か
#     2) NDVI の H 寄与が ALB を控除しても残るか (partial corr)
# ================================================================
def albedo_feedback_check(oran_raw_path, oran_merged):
    print(f"\n{'=' * 60}\n[H4 Oran] アルベド・フィードバック検証\n{'=' * 60}")
    df = pd.read_csv(oran_raw_path)
    if "ALB" not in df.columns:
        print("  ALB 列なし; skip"); return
    src_col = "TIMESTAMP" if "TIMESTAMP" in df.columns else "DateTime"
    ts = df[src_col].astype(str).str.strip()
    dt = pd.Series(pd.NaT, index=df.index, dtype="datetime64[ns]")
    for fmt in ["%Y/%m/%d %H:%M:%S", "%Y/%m/%d"]:
        miss = dt.isna()
        if miss.any():
            dt = dt.where(~miss,
                          pd.to_datetime(ts[miss], format=fmt, errors="coerce"))
    df["datetime"] = dt
    df = df.dropna(subset=["datetime"]).copy()

    df["ALB"] = pd.to_numeric(df["ALB"], errors="coerce")
    df["SW_IN"] = pd.to_numeric(df["SW_IN"], errors="coerce")
    df = df.where(df > -9000)  # 念のためセンチネル
    # 日中のみ ALB を集計 (SW_IN > 100 W/m²)
    daytime = df[df["SW_IN"] > 100].copy()
    daytime["date"] = daytime["datetime"].dt.normalize()
    alb_daily = daytime.groupby("date", as_index=False).agg(ALB_day=("ALB", "median"))
    alb_daily = alb_daily[(alb_daily["ALB_day"] > 0.05) &
                           (alb_daily["ALB_day"] < 0.5)]
    print(f"  ALB 日中中央値: n={len(alb_daily)} 日, "
          f"med={alb_daily['ALB_day'].median():.3f}  "
          f"range=[{alb_daily['ALB_day'].quantile(0.05):.3f}, "
          f"{alb_daily['ALB_day'].quantile(0.95):.3f}]")

    merged = oran_merged.merge(alb_daily, on="date", how="left")
    sub = merged.dropna(subset=["NDVI", "H", "ALB_day"])
    if len(sub) < 30:
        print("  サンプル不足"); return

    r_na, _ = stats.spearmanr(sub["NDVI"], sub["ALB_day"])
    r_nh, _ = stats.spearmanr(sub["NDVI"], sub["H"])
    r_ah, _ = stats.spearmanr(sub["ALB_day"], sub["H"])
    print(f"  ρ(NDVI, ALB) = {r_na:+.3f}   "
          f"ρ(NDVI, H)   = {r_nh:+.3f}   "
          f"ρ(ALB, H)    = {r_ah:+.3f}")

    # 部分相関: NDVI と H の関係が ALB を控除すると残るか
    z = lambda s: (s.values - s.mean()) / s.std()
    z_n, z_h, z_a = z(sub["NDVI"]), z(sub["H"]), z(sub["ALB_day"])

    def part(y, x, w):
        ry = y - w * (np.dot(w, y) / np.dot(w, w))
        rx = x - w * (np.dot(w, x) / np.dot(w, w))
        return float(np.corrcoef(ry, rx)[0, 1])
    pr = part(z_h, z_n, z_a)
    print(f"  partial r(H, NDVI | ALB) = {pr:+.3f}")
    if abs(pr) < abs(r_nh) * 0.5:
        print("  → ALB を控除すると NDVI~H が大きく減衰 = "
              "アルベド経由が主機構")
    else:
        print("  → ALB 控除後も NDVI~H が残存 = "
              "気孔/エネルギー再分配など別機構が主")


# ================================================================
# H7: 月別 NDVI のセカンドピーク検出 (二期作見落としチェック)
#
#   Oran 月フィルター 11–6 月で生育期と仮定しているが、もし 7–10 月にも
#   小さなピーク (夏作) があれば、生育期定義から漏れている。
# ================================================================
def detect_second_peak(per_site):
    print(f"\n{'=' * 60}\n[H7] 月別 NDVI セカンドピーク検出\n{'=' * 60}")
    for site, (df, _) in per_site.items():
        sub = df.dropna(subset=["NDVI"]).copy()
        if "date" not in sub.columns or len(sub) == 0:
            continue
        sub["month"] = pd.to_datetime(sub["date"]).dt.month
        med = sub.groupby("month")["NDVI"].median()
        med_arr = med.reindex(range(1, 13)).values
        # 単純なピーク検出: 局所最大が 2 つ以上あるか?
        peaks = []
        for i in range(12):
            left  = med_arr[(i - 1) % 12]
            right = med_arr[(i + 1) % 12]
            v = med_arr[i]
            if np.isfinite(v) and np.isfinite(left) and np.isfinite(right):
                if v >= left and v >= right and v > np.nanmean(med_arr):
                    peaks.append((i + 1, v))
        peaks_str = "  ".join(f"M{m:02d}={v:.2f}" for m, v in peaks)
        print(f"  {site:8s} peaks: {peaks_str}  "
              f"(全{len([p for p in peaks])}峰)")
        if len(peaks) >= 2:
            print(f"    → セカンドピーク検出. 二期作の可能性")
        else:
            print(f"    → 単峰. 単作と整合")


# ================================================================
# H1+H8: 灌漑効果と深根効果の分離 (Tarazona のみ)
#
#   仮説: Tarazona の高 EF は深根アクセスではなく、灌漑による
#         安定的な表層水供給だけで説明できるかもしれない。
#   検証: 「灌漑からの経過日数」別に LE/EF を見る。
#         - 灌漑だけが効いているなら、灌漑後 1 週間以内に EF が高く、
#           7 日以降は急減衰
#         - 深根が効いているなら、長期間 (7 日以降も) EF が維持される
# ================================================================
def irrigation_lag_analysis(tara_df, save_dir):
    """Irrig_mm の入った Tarazona DataFrame で「最終灌漑からの日数」別 LE/EF."""
    if "Irrig_mm" not in tara_df.columns:
        print("\n[H1+H8] Irrig_mm 列なし; skip"); return None

    print(f"\n{'=' * 60}\n[H1+H8] 灌漑経過日数別 LE/EF (Tarazona)\n{'=' * 60}")

    df = tara_df.copy().sort_values("date").reset_index(drop=True)
    df["Irrig_mm"] = pd.to_numeric(df["Irrig_mm"], errors="coerce").fillna(0)
    irrig_event = df["Irrig_mm"] > 0.5  # 0.5 mm 以上を灌漑日と定義
    n_events = int(irrig_event.sum())
    print(f"  灌漑日 (>0.5 mm): {n_events} 日")
    if n_events < 5:
        print("  灌漑イベント不足; skip"); return None

    # 各日について、直近の灌漑日からの経過日数
    days_since = np.full(len(df), np.nan)
    last = -10**9
    for i in range(len(df)):
        if irrig_event.iloc[i]:
            last = i
            days_since[i] = 0
        elif last >= 0:
            days_since[i] = i - last
    df["days_since_irrig"] = days_since

    # ビンを区切って集計
    bins = [-0.5, 0.5, 3.5, 7.5, 14.5, 30.5, 1e9]
    labels = ["irrig day", "1-3", "4-7", "8-14", "15-30", ">30"]
    df["lag_bin"] = pd.cut(df["days_since_irrig"], bins=bins, labels=labels)

    # 生育期に絞る
    growing = df[df.get("phen") == "growing"].copy() if "phen" in df.columns else df

    print(f"  生育期 n={len(growing)}")
    rows = []
    for lab in labels:
        sub = growing[growing["lag_bin"] == lab]
        if len(sub) >= 5:
            le_med = sub["LE"].median()
            ef_med = sub["EF"].median()
            et_med = sub["ET"].median()
            rows.append((lab, len(sub), le_med, ef_med, et_med))
            print(f"    lag={lab:>10s}  n={len(sub):3d}  "
                  f"LE_med={le_med:6.1f}  EF_med={ef_med:.3f}  ET_med={et_med:.2f}")

    # Mann-Whitney: 灌漑直後 (lag 1-3) vs 1 週間以上後 (8-14)
    sub_a = growing[growing["lag_bin"] == "1-3"]["EF"].dropna()
    sub_b = growing[growing["lag_bin"] == "8-14"]["EF"].dropna()
    if len(sub_a) >= 5 and len(sub_b) >= 5:
        u, p = stats.mannwhitneyu(sub_a, sub_b, alternative="greater")
        print(f"  EF (lag 1-3) > EF (lag 8-14)?  "
              f"medians: {sub_a.median():.3f} vs {sub_b.median():.3f}  "
              f"p={p:.2e}")
        if p > 0.05:
            print("  → 灌漑後 1 週間以上経っても EF が維持されている "
                  "= 深根アクセスの状況証拠")
        else:
            print("  → 灌漑直後の方が有意に高い "
                  "= 灌漑寄与が支配的")

    # 図: 経過日数 vs EF の box
    plot_data = []
    plot_labels = []
    for lab in labels:
        sub = growing[growing["lag_bin"] == lab]["EF"].dropna()
        if len(sub) >= 3:
            plot_data.append(sub.values); plot_labels.append(f"{lab}\n(n={len(sub)})")
    if plot_data:
        with plt.rc_context(PAPER_RC):
            fig, ax = plt.subplots(figsize=(8, 4.5))
            bp = ax.boxplot(plot_data, labels=plot_labels, patch_artist=True,
                            widths=0.55, showfliers=False)
            for patch in bp["boxes"]:
                patch.set_facecolor(SITE_COLOR["Tarazona"]); patch.set_alpha(0.6)
            ax.set_ylabel("EF")
            ax.set_xlabel("Days since last irrigation")
            ax.set_title("Tarazona: EF vs days since irrigation\n(growing season)",
                         loc="left", fontweight="bold")
            ax.grid(True, alpha=0.25, axis="y")
            ax.set_ylim(-0.05, 1.55)
            out = Path(save_dir) / "C_H1_irrigation_lag.png"
            plt.savefig(out, dpi=180, bbox_inches="tight"); plt.close()
            print(f"  [保存] {out}")

    return df


def deep_dive_gpp(df, site):
    """
    NDVI~GPP_proxy 信号を 生育期限定 / VPD 層別 / 月別ラグ で詳査。
    Tarazona は ρ=+0.67 と強いので、その信号源が本当に植生活性か、
    単なる季節同位相 (Rn と NDVI が同時に上がるだけ) かを切り分ける。
    """
    print(f"\n{'='*60}\n[GPP 深掘り {site}]\n{'='*60}")
    if "phen" not in df.columns:
        print("  phen 列なし; スキップ"); return

    g = df[df["phen"] == "growing"].dropna(subset=["NDVI"])
    print(f"  生育期 n={len(g)}")

    for var in ["LE", "ET", "EF", "GPP_proxy"]:
        if var not in g.columns: continue
        sub = g[["NDVI", var]].dropna()
        if len(sub) < 10: continue
        r, p = stats.spearmanr(sub["NDVI"], sub[var])
        print(f"  growing only: NDVI ~ {var:9s} ρ={r:+.3f}  p={p:.2e}  n={len(sub)}")

    # GPP_proxy の構成要素を分離: NDVI 単独 vs Rn 単独 と LE
    if {"Rn", "NDVI", "LE"}.issubset(g.columns):
        sub = g.dropna(subset=["Rn", "NDVI", "LE"])
        if len(sub) > 20:
            r1, _ = stats.spearmanr(sub["NDVI"], sub["LE"])
            r2, _ = stats.spearmanr(sub["Rn"],   sub["LE"])
            r3, _ = stats.spearmanr(sub["NDVI"] * sub["Rn"], sub["LE"])
            print(f"  growing only [n={len(sub)}]: ρ(NDVI,LE)={r1:+.3f}  "
                  f"ρ(Rn,LE)={r2:+.3f}  ρ(NDVI×Rn,LE)={r3:+.3f}")
            print(f"  → NDVI×Rn が NDVI 単独より強ければ、植生×放射の積が効いている")

    # VPD 三分位での層別
    if "VPD" in g.columns:
        sub = g.dropna(subset=["NDVI", "LE", "VPD"]).copy()
        if len(sub) > 30:
            sub["VPD_q"] = pd.qcut(sub["VPD"], q=3, labels=["low", "mid", "high"])
            print("  VPD 層別 NDVI~LE:")
            for q in ["low", "mid", "high"]:
                ss = sub[sub["VPD_q"] == q]
                if len(ss) >= 10:
                    r, p = stats.spearmanr(ss["NDVI"], ss["LE"])
                    print(f"    VPD={q:4s} (n={len(ss):3d}): ρ={r:+.3f}  p={p:.2e}")


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

    # 生 CSV ヘッダ確認 (列名のズレを発見するため)
    inspect_raw_oran_csv(A_PATHS["oran_ec"])
    inspect_raw_oran_csv(A_PATHS["tara_ec"])

    # 先頭行のタイムスタンプを目視 + 生 NETRAD 検証
    dump_oran_first_rows(A_PATHS["oran_ec"], n=10)
    deep_inspect_oran_raw(A_PATHS["oran_ec"])
    diagnose_oran_parse_failures(A_PATHS["oran_ec"])
    inspect_tarazona_units(A_PATHS["tara_ec"])

    print("\n--- C クリーン版 (TIMESTAMP 修正 + 単位整合) ---")
    oran_ec = normalize_swc(load_oran_ec_clean(A_PATHS["oran_ec"]), "Oran")
    diagnose_ec(oran_ec, "Oran (clean)")

    tara_ec = normalize_swc(load_tarazona_ec_clean(A_PATHS["tara_ec"]), "Tarazona")
    diagnose_ec(tara_ec, "Tarazona (clean)")

    oran_m, oran_sg, oran_meta = run_site_C("Oran",     oran_ec, NDVI_APPEEARS_CSV)
    tara_m, tara_sg, tara_meta = run_site_C("Tarazona", tara_ec, NDVI_APPEEARS_CSV)

    per_site = {"Oran": (oran_m, oran_sg), "Tarazona": (tara_m, tara_sg)}
    plot_ndvi_timeseries(per_site, SAVE_DIR)
    plot_ndvi_vs_flux(per_site, SAVE_DIR)

    # A 連結 — 厳しい / 緩和 / Otsu 由来の3パターンを並べる
    a_res_strict   = cross_site_dry_canopy(oran_m, tara_m, mode="ndvi_p67")
    a_res_relaxed  = cross_site_dry_canopy(oran_m, tara_m, mode="ndvi_p50")
    a_res_growing  = cross_site_dry_canopy(oran_m, tara_m, mode="growing_flag")
    a_results = a_res_relaxed if a_res_relaxed else a_res_strict

    # GPP 信号源の深掘り
    deep_dive_gpp(oran_m, "Oran")
    deep_dive_gpp(tara_m, "Tarazona")

    # 部分相関で NDVI / Rn の独立寄与を切り分け
    partial_correlation_analysis(oran_m, "Oran")
    partial_correlation_analysis(tara_m, "Tarazona")

    # H1+H8: 灌漑経過日数別 LE/EF (Tarazona のみ)
    # tara_m は Irrig_mm を含む load_tarazona_ec_clean 由来であることを期待
    if "Irrig_mm" in tara_m.columns:
        irrigation_lag_analysis(tara_m, SAVE_DIR)
    else:
        # ローダから落ちている場合は再取得
        from data_loaders import load_tarazona_ec_clean as _ltec
        tara_full = _ltec(A_PATHS["tara_ec"], verbose=False)
        if "Irrig_mm" in tara_full.columns:
            tara_with_irr = tara_m.merge(
                tara_full[["date", "Irrig_mm"]], on="date", how="left")
            irrigation_lag_analysis(tara_with_irr, SAVE_DIR)

    # H2: NDVI 飽和を EVI で検証 (両サイト)
    ndvi_saturation_check(oran_m, "Oran",     SAVE_DIR)
    ndvi_saturation_check(tara_m, "Tarazona", SAVE_DIR)

    # H4: アルベド・フィードバック (Oran のみ)
    albedo_feedback_check(A_PATHS["oran_ec"], oran_m)

    # H7: 月別 NDVI のセカンドピーク検出
    detect_second_peak(per_site)

    # paper-quality figures
    plot_final_summary(oran_m, tara_m, SAVE_DIR)
    plot_monthly_ndvi(per_site, SAVE_DIR)
    plot_partial_corr(oran_m, tara_m, SAVE_DIR)

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
