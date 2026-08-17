"""エネルギー収支の閉合 EBR を site-year 品質指標として計算するコア。

渦相関(EC)は乱流で運ばれる熱・水・炭素を測るが、乱流が足りないと H+LE が正味放射
Rn−G を下回る。閉合率 **EBR = Σ(H+LE) / Σ(Rn−G)** は、その site-year がどれだけ乱流を
捉えられたかの物理的サイン。EBR が低い(<0.7)＝乱流不足＝炭素フラックスも取りこぼしの
疑い＝除外候補、という**物理裏付けのある除外基準**。

RK_VARS(11変数)に Rn/G は含まれないため、生 CSV から H/LE/Rn/G 列を直接読む。
`rank_sites`(品質ランキング列)と `research/energy_closure_step28.py`(合成検証つき CLI)の
両方がこのモジュールを使う（計算を一箇所に集約してズレを防ぐ）。
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import AnalysisConfig
from .sites import SiteSpec
from .preprocess import (
    _read_table_header, _read_table_columns, find_corevars_files)

# H/LE/Rn/G の実カラム候補(japanflux2024 / FLUXNET2015 の表記揺れを吸収)。
# gap-fill 済み(_F_MDS)を優先し、無ければ生列(_1_1_1)へフォールバック。
CAND: dict[str, list[str]] = {
    "H":  ["H_F_MDS", "H_F", "H_1_1_1", "H"],
    "LE": ["LE_F_MDS", "LE_F", "LE_1_1_1", "LE"],
    "Rn": ["NETRAD_F_MDS", "NETRAD_F", "NETRAD_1_1_1", "NETRAD", "RN", "NET_RAD"],
    "G":  ["G_F_MDS", "G_F_MDS_1", "G_F", "G_1_1_1", "G"],
}


def _resolve(header, cands):
    hset = set(header)
    for c in cands:
        if c in hset:
            return c
    return None


def load_energy(site: SiteSpec, months, qc_max):
    """生 CSV から H/LE/Rn/G を直接読み、対象月に絞った 30 分値 DataFrame。

    必須は H/LE/Rn。G(地中熱)は無いサイト(例 JP-Ta2)があるので任意扱い＝G=0 近似。
    --qc-max 指定時は各列の ``_QC`` で低品質補完を NaN 化する。

    Returns
    -------
    (raw, cols, missing, g_approx):
        raw は H/LE/Rn/G 列を持つ DataFrame(必須欠如時は None)。cols は解決した実列名。
        missing は欠けた必須列(空なら成功)。g_approx は G を 0 近似したか。
    """
    cfg = AnalysisConfig(qc_max=qc_max) if qc_max is not None else AnalysisConfig()
    files = find_corevars_files(site)
    header0 = _read_table_header(files[0])
    cols = {k: _resolve(header0, cands) for k, cands in CAND.items()}
    missing = [k for k in ("H", "LE", "Rn") if cols[k] is None]
    if missing:
        return None, cols, missing, False
    g_approx = cols["G"] is None
    present = {k: c for k, c in cols.items() if c is not None}

    qc_of = {}
    if qc_max is not None:
        for k, c in present.items():
            qc = c + "_QC"
            if qc in set(header0):
                qc_of[k] = qc
    want = {"TIMESTAMP_START", *present.values(), *qc_of.values()}
    parts = []
    for f in files:
        df = _read_table_columns(f, want)
        ts = pd.to_datetime(
            pd.to_numeric(df["TIMESTAMP_START"]).astype("int64").astype(str),
            format="%Y%m%d%H%M")
        df = df.drop(columns=["TIMESTAMP_START"]); df.index = ts
        parts.append(df)
    raw = pd.concat(parts)
    raw = raw[~raw.index.duplicated(keep="first")].sort_index()
    raw = raw.replace(cfg.na_sentinel, np.nan)
    for k, qc in qc_of.items():
        raw[present[k]] = raw[present[k]].where(raw[qc] <= qc_max)
    raw = raw.rename(columns={v: k for k, v in present.items()})
    if g_approx:
        raw["G"] = 0.0
    raw = raw[["H", "LE", "Rn", "G"]]
    if months:
        raw = raw[raw.index.month.isin(months)]
    return raw, cols, [], g_approx


def closure(df):
    """EBR=Σ(H+LE)/Σ(Rn−G) と、原点通し回帰の傾き・R²・n を返す。"""
    ok = df[["H", "LE", "Rn", "G"]].notna().all(axis=1)
    d = df.loc[ok]
    turb = (d["H"] + d["LE"]).to_numpy()
    avail = (d["Rn"] - d["G"]).to_numpy()
    n = len(d)
    if n < 200:
        return None
    ebr = float(turb.sum() / avail.sum()) if avail.sum() != 0 else np.nan
    s = float((avail @ turb) / (avail @ avail)) if (avail @ avail) != 0 else np.nan
    ss_res = float(((turb - s * avail) ** 2).sum())
    ss_tot = float(((turb - turb.mean()) ** 2).sum())
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return {"ebr": ebr, "slope": s, "r2": r2, "n": n}


def by_year(df):
    """年ごとの closure() 辞書。"""
    out = {}
    for y, g in df.groupby(df.index.year):
        c = closure(g)
        if c:
            out[int(y)] = c
    return out


def verdict(ebr):
    if not np.isfinite(ebr):
        return "—"
    if ebr < 0.7:
        return "⚠ 不閉合(乱流不足の疑い/除外候補)"
    if ebr <= 1.05:
        return "✅ 良好(0.7–1.05)"
    return "△ 過閉合(Rn/G測器 or 貯留無視)"


def site_ebr(site: SiteSpec, months, qc_max=1):
    """1 サイトの EBR 要約を返す簡便関数（rank_sites 用）。

    Returns dict {ebr, slope, r2, n, g_approx, n_years, n_bad} or {"note": 理由}。
    n_bad は EBR<0.7 の年数。既定 qc_max=1（実測寄り）。
    """
    raw, cols, missing, g_approx = load_energy(site, months, qc_max)
    if missing:
        return {"note": "必須列欠如:" + ",".join(missing)}
    c = closure(raw)
    if c is None:
        return {"note": "データ不足"}
    yr = by_year(raw)
    n_bad = sum(1 for cc in yr.values()
                if np.isfinite(cc["ebr"]) and cc["ebr"] < 0.7)
    return {**c, "g_approx": g_approx, "n_years": len(yr), "n_bad": n_bad}
