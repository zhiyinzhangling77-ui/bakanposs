"""旗39：呼吸残差の記憶は「生物」か「分割窓アーティファクト」か（DT/NT で判別・タワーのみ）。

旗38 で GOSIF SIF が記憶を説明しなかった。残る2説：(A)GOSIFが力不足、(B)記憶が分割アルゴリズムの
移動窓(Reichstein/Lasslop 4-7日)の平滑さ＝生物でない。**(B) は独立衛星でも説明できない**ので旗38の null と整合。
これをタワーだけで安く判別する：

  もし記憶が**分割窓**なら、昼分割(RECO_DT)と夜分割(RECO_NT)は違う仮定・違う窓なので
  **記憶の時間スケール(e-folding)が変わる**はず。もし記憶が**生物(基質等)**なら、同じ呼吸を別法で割っただけ
  なので **DT/NT で同じ**はず。＝TROPOMI に行く前の安い判別。

判定：DT と NT の e-folding/ACF が**大きく違えば→分割窓の寄与が濃い(アーティファクト寄り)**、
**同じなら→生物由来と整合(＝真SIFで探す価値がある)**。非対称な検定（「違う」は情報量大、「同じ」は
「同一窓 or 生物」で確定はしない）。

    python research/memory_partition_step39.py                                 # 合成で検証
    python research/memory_partition_step39.py --sites JP-Tak JP-Mse JP-Ta2 CN-HaM MN-Hst --qc-max 1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from memory_timescale_step37 import DRIVERS, _autocorr, efolding_days, fit_residual


def _nt_of(dt_col):
    return dt_col.replace("_DT_", "_NT_").replace("_DT", "_NT") if "DT" in dt_col else None


def load_dt_nt_daily(site, months, qc_max):
    """生CSVから RECO_DT/RECO_NT＋駆動を読み、夏を日平均に。"""
    import pandas as pd
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import (
        _read_table_header, _read_table_columns, find_corevars_files)
    spec = get_site(site)
    cfg = AnalysisConfig(qc_max=qc_max) if qc_max is not None else AnalysisConfig()
    vmap = spec.var_map()
    dt_col = vmap["GER"]; nt_col = _nt_of(dt_col)
    files = find_corevars_files(spec)
    hset = set(_read_table_header(files[0]))
    if nt_col is None or nt_col not in hset or dt_col not in hset:
        return None
    # 駆動の実カラム
    drv = {v: vmap[v] for v in DRIVERS if v in vmap and vmap[v] in hset and not vmap[v].startswith("@")}
    want = {"TIMESTAMP_START", dt_col, nt_col, *drv.values()}
    parts = []
    for f in files:
        df = _read_table_columns(f, want)
        ts = pd.to_datetime(pd.to_numeric(df["TIMESTAMP_START"]).astype("int64").astype(str),
                            format="%Y%m%d%H%M")
        df = df.drop(columns=["TIMESTAMP_START"]); df.index = ts
        parts.append(df)
    raw = pd.concat(parts)
    raw = raw[~raw.index.duplicated(keep="first")].sort_index().replace(cfg.na_sentinel, np.nan)
    ren = {dt_col: "GER_DT", nt_col: "GER_NT", **{c: v for v, c in drv.items()}}
    raw = raw.rename(columns=ren)
    raw = raw[raw.index.month.isin(sorted(months))]
    keep = ["GER_DT", "GER_NT"] + [v for v in DRIVERS if v in raw.columns]
    return raw[keep].groupby(raw.index.normalize()).mean().dropna()


def memory_of(daily, ger_col):
    d = daily.copy()
    d["GER"] = d[ger_col]
    r2, res = fit_residual(d)
    ef, _ = efolding_days(res)
    return {"r2": r2, "ac": _autocorr(res), "ef": ef}


def make_synth(kind, days=900, seed=0):
    """DT/NT に相当する2系列。same=同じ隠れ記憶、diff=違う時間スケールの記憶。"""
    import pandas as pd
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2005-06-01", periods=days, freq="D")
    Ta = 20 + rng.normal(0, 2, days); Rg = 300 + rng.normal(0, 40, days)
    th = np.clip(0.3 + rng.normal(0, .05, days), .05, .6); VPD = np.clip(.5 + .06 * (Ta - 15), .1, None)

    def ar(phi):
        x = np.zeros(days)
        for i in range(1, days):
            x[i] = phi * x[i - 1] + rng.normal(0, 0.5)
        return (x - x.mean()) / x.std()
    base = 2.0 * np.exp(0.06 * (Ta - 20))
    S = ar(0.8)                                   # 共通の隠れ記憶（≈4-5日）
    df = pd.DataFrame({"Rg": Rg, "Ta": Ta, "VPD": VPD, "Ts": Ta - 1, "th": th,
                       "P": np.clip(rng.normal(0, 1, days), 0, None),
                       "gH": Rg * .2, "gLE": Rg * .3}, index=idx)
    if kind == "same":       # DT/NT とも同じ記憶 S（生物由来を模す）
        df["GER_DT"] = base + 1.5 * S + rng.normal(0, .15, days)
        df["GER_NT"] = base + 1.5 * S + rng.normal(0, .25, days)
    else:                    # DT/NT で違う時間スケール（分割窓の違いを模す）
        df["GER_DT"] = base + 1.5 * S + rng.normal(0, .15, days)
        df["GER_NT"] = base + 1.5 * ar(0.4) + rng.normal(0, .15, days)   # 速い記憶
    return df.dropna()


def _report(daily, tag):
    dt = memory_of(daily, "GER_DT"); nt = memory_of(daily, "GER_NT")
    dEF = abs(dt["ef"] - nt["ef"])
    print(f"  {tag:<8} DT: R²={dt['r2']:.2f} ACF={dt['ac']:+.2f} e-fold={dt['ef']}日"
          f"  ｜ NT: R²={nt['r2']:.2f} ACF={nt['ac']:+.2f} e-fold={nt['ef']}日  ｜Δe-fold={dEF}日")
    return dt, nt, dEF


def main():
    p = argparse.ArgumentParser(description="呼吸残差の記憶は生物か分割窓か(DT/NT判別)")
    p.add_argument("--sites", nargs="+")
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--qc-max", type=int, default=None)
    a = p.parse_args()

    if not a.sites:
        print("=== 旗39 合成検証：DT/NT の記憶が同じか違うかを判別できるか ===")
        _report(make_synth("same"), "同記憶")
        _report(make_synth("diff"), "異記憶")
        print("\n  → 同記憶は Δe-fold≈0（生物由来と整合）、異記憶は Δe-fold 大（分割窓の違いを検出）が期待。")
        return

    print(f"=== 旗39 実データ 記憶は生物か分割窓か（DT/NT 判別, QC≤{a.qc_max}, 月={a.month}）===")
    print("  DT/NT で e-folding が同じ→生物由来と整合(真SIFで探す価値)／大きく違う→分割窓の寄与濃い(アーティファクト寄り)\n")
    from japanflux_pn.sites import get_site  # noqa: F401
    diffs = []
    for s in a.sites:
        try:
            daily = load_dt_nt_daily(s, a.month, a.qc_max)
        except Exception as e:
            print(f"  {s:<8} SKIP {type(e).__name__}: {e}"); continue
        if daily is None or len(daily) < 60:
            print(f"  {s:<8} NT列なし/データ不足"); continue
        _, _, dEF = _report(daily, s)
        diffs.append(dEF)
    if diffs:
        med = np.median(diffs)
        print(f"\n  Δe-fold の中央値 = {med:.0f}日")
        if med <= 1:
            print("  → ○ DT/NT で記憶がほぼ同じ＝生物由来と整合（分割窓だけでは説明しにくい）＝真SIF(TROPOMI)で探す価値あり")
        elif med >= 3:
            print("  → ⚠ DT/NT で記憶が大きく違う＝分割窓の寄与が濃い（アーティファクト寄り）＝記憶の相当部分は分割由来の疑い")
        else:
            print("  → △ 中間＝生物と分割窓が混在の可能性")
    print("  留保：非対称な検定（『違う』は情報量大、『同じ』は同一窓 or 生物で確定せず）。DT/NT は同じNEEを別法で割った派生量。")


if __name__ == "__main__":
    main()
