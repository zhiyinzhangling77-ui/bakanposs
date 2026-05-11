"""
data_loaders.py — EC データの正しいローダ (解析 A/B/C 共通)

解析C で発見した v9 ローダのバグ群を全て修正したクリーン版を、
解析 A/B でも import して使えるよう切り出したモジュール。

主な修正:
1. Oran タイムスタンプ: pandas の自動推定で先頭行の '2018/01/01' から
   '%Y/%m/%d' を確定→残り全部を NaT 化していた問題を、
   日時付き → 日付のみ → mixed の三段フォールバックで解消
2. -9999 等のセンチネル値マスク (実測 0 件だが保険として保持)
3. Oran VPD: hPa → kPa 統一 (÷10)、|VPD|>10 kPa の外れ値除去
4. Tarazona ET: ET_avg ではなく ET_sum (日合計 mm) を採用
5. Tarazona VPD: VPD_kPa 列を直接採用 (VPD_mean は Pa)

最終単位:
  LE / H / G / Rn : W/m²
  ET              : mm/day
  VPD             : kPa
  SWC             : %  (m³/m³ で来た場合は normalize_swc で要変換)
"""

import numpy as np
import pandas as pd

EF_DENOM_MIN = 10.0       # Rn-G > これ未満の日は EF を計算しない (W/m²)
SENTINEL_THR = -9000.0    # これ以下は欠損扱い
VPD_MAX_KPA  = 10.0       # 物理的にあり得ない VPD はマスク


def _recover_datetime_from_julian(df, mask):
    """year + Julian + Time_hours 列から datetime を復元.
    TIMESTAMP/DateTime が空欄 ('nan') の行に対するフォールバック."""
    needed = ["year", "Julian", "Time_hours"]
    if not all(c in df.columns for c in needed):
        return None
    sub = df.loc[mask, needed].apply(pd.to_numeric, errors="coerce")
    valid = sub.notna().all(axis=1)
    if not valid.any():
        return None
    sub = sub.loc[valid]
    base = pd.to_datetime(
        sub["year"].astype(int).astype(str)
        + sub["Julian"].astype(int).astype(str).str.zfill(3),
        format="%Y%j", errors="coerce")
    hours_td = pd.to_timedelta(sub["Time_hours"].astype(float), unit="h")
    recovered = base + hours_td
    out = pd.Series(pd.NaT, index=df.index, dtype="datetime64[ns]")
    out.loc[sub.index] = recovered.values
    return out


def load_oran_ec_clean(filepath, verbose=True):
    """Oran 半時間値 CSV → 日次集計.
    Ameriflux の TIMESTAMP 列を二段＋mixed の三段フォールバックでパース.
    最後に year/Julian/Time_hours 列から復元 (TIMESTAMP='nan' 行救済)."""
    df = pd.read_csv(filepath)
    n_raw = len(df)

    src_col = "TIMESTAMP" if "TIMESTAMP" in df.columns else "DateTime"
    ts = df[src_col].astype(str).str.strip()
    dt = pd.Series(pd.NaT, index=df.index, dtype="datetime64[ns]")
    fmts = [
        "%Y/%m/%d %H:%M:%S", "%Y/%m/%d %H:%M",
        "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M",
        "%Y/%m/%d", "%Y-%m-%d",
        "%m/%d/%Y %H:%M:%S", "%m/%d/%Y",
    ]
    for fmt in fmts:
        miss = dt.isna()
        if not miss.any():
            break
        dt2 = pd.to_datetime(ts[miss], format=fmt, errors="coerce")
        dt = dt.where(~miss, dt2)
    if dt.notna().sum() < n_raw * 0.95:
        miss = dt.isna()
        dt2 = pd.to_datetime(ts[miss], format="mixed", errors="coerce")
        dt = dt.where(~miss, dt2)
    n_after_str_parse = int(dt.notna().sum())

    # フォールバック: year/Julian/Time_hours から復元
    miss = dt.isna()
    if miss.any():
        recovered = _recover_datetime_from_julian(df, miss)
        if recovered is not None:
            n_recover = int((recovered.notna() & miss).sum())
            dt = dt.where(~miss, recovered)
            if verbose and n_recover > 0:
                print(f"[Oran loader] year/Julian/Time_hours から "
                      f"{n_recover} 行を復元")

    df["datetime"] = dt
    n_parsed = int(df["datetime"].notna().sum())
    if verbose:
        days_n = max(1, df["datetime"].dt.normalize().nunique())
        recover_gain = n_parsed - n_after_str_parse
        print(f"[Oran loader] {n_parsed}/{n_raw} subdaily records parsed "
              f"(~{n_parsed / days_n:.1f}/day; "
              f"+{recover_gain} from Julian fallback)")

    df = df.dropna(subset=["datetime"]).sort_values("datetime").reset_index(drop=True)

    flux_cols = ["SWC_1_1_1", "LE", "H", "G", "NETRAD", "VPD", "ET"]
    for col in flux_cols:
        s = pd.to_numeric(df[col], errors="coerce")
        df[col] = s.where(s > SENTINEL_THR, np.nan)

    df["date"] = df["datetime"].dt.normalize()
    daily = df.groupby("date", as_index=False).agg(
        SWC=("SWC_1_1_1", "mean"), LE=("LE", "mean"), H=("H", "mean"),
        G=("G", "mean"), Rn=("NETRAD", "mean"),
        VPD=("VPD", "mean"), ET=("ET", "sum"))

    denom = daily["Rn"] - daily["G"]
    valid = denom > EF_DENOM_MIN
    daily["EF"] = np.nan
    daily.loc[valid, "EF"] = (daily.loc[valid, "LE"] / denom[valid]).clip(0, 1.5)
    daily = daily[(daily["SWC"] > 0) & (daily["SWC"] < 100)].reset_index(drop=True)

    daily["VPD"] = daily["VPD"] / 10.0  # hPa → kPa
    daily.loc[(daily["VPD"] > VPD_MAX_KPA) | (daily["VPD"] < 0), "VPD"] = np.nan
    daily["site"] = "Oran"
    if verbose:
        print(f"[Oran loader] {len(daily)} 日分 (VPD: hPa→kPa, 外れ値除去)")
    return daily


def load_tarazona_ec_clean(filepath, verbose=True):
    """Tarazona 日次集計 CSV → ET_sum / VPD_kPa を優先採用."""
    df = pd.read_csv(filepath)
    df["date"] = pd.to_datetime(df["date"], errors="coerce")
    df = df.dropna(subset=["date"]).sort_values("date").reset_index(drop=True)

    rename = {"SWC_avg": "SWC", "LE_avg": "LE", "H_avg": "H",
              "G_avg": "G", "NetRad_avg": "Rn"}
    df = df.rename(columns=rename)

    if "ET_sum" in df.columns:
        df["ET"] = pd.to_numeric(df["ET_sum"], errors="coerce")
        et_src = "ET_sum"
    else:
        df["ET"] = pd.to_numeric(df["ET_avg"], errors="coerce")
        et_src = "ET_avg"

    if "VPD_kPa" in df.columns:
        df["VPD"] = pd.to_numeric(df["VPD_kPa"], errors="coerce")
        vpd_src = "VPD_kPa"
    elif "VPD_mean" in df.columns:
        df["VPD"] = pd.to_numeric(df["VPD_mean"], errors="coerce") / 1000.0
        vpd_src = "VPD_mean/1000"
    else:
        df["VPD"] = np.nan
        vpd_src = "missing"

    for c in ["SWC", "LE", "H", "G", "Rn"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    denom = df["Rn"] - df["G"]
    valid = denom > EF_DENOM_MIN
    df["EF"] = np.nan
    df.loc[valid, "EF"] = (df.loc[valid, "LE"] / denom[valid]).clip(0, 1.5)

    keep = ["date", "SWC", "LE", "H", "G", "Rn", "VPD", "ET", "EF"]
    extra = [c for c in ["Irrig_mm", "Rain_mm", "IrrigRain_mm",
                         "GPP_avg", "NDVI_orig", "NDVI_interp",
                         "Fv_interp"] if c in df.columns]
    df = df[[c for c in keep if c in df.columns] + extra].copy()
    df["site"] = "Tarazona"
    if verbose:
        print(f"[Tarazona loader] {len(df)} 日分 (ET←{et_src}, VPD←{vpd_src})")
    return df


def normalize_swc(df, site_name):
    """SWC が m³/m³ (0–1) で来ていたら % (0–100) に変換."""
    if "SWC" not in df.columns:
        return df
    if df["SWC"].max() <= 1.0:
        df = df.copy()
        df["SWC"] = df["SWC"] * 100.0
        print(f"  SWC [{site_name}]: m³/m³ → %")
    return df
