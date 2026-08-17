"""旗38：SIF は「気象で作れなかった呼吸残差の記憶」を説明するか（交絡を破る検証）。

旗37 で GER 日残差に約4日の記憶があり、気象の後方窓（1〜30日）では作れないと分かった。この4日は
文献の基質供給ラグ(1〜5日)と一致するが、GER 分割の移動窓(4〜7日)とも一致＝生物か分割かフラックス
単独では切れない。**SIF(太陽光誘起蛍光)は分割に依らない独立な光合成シグナル**なので、もし SIF が残差の
記憶を説明したら＝基質供給が本物と証明でき、分割窓アーティファクト説を棄却できる（旗37 の交絡を破る唯一の道）。

やること：各サイトで
  1. GER 日残差＝旗37 の瞬間気象モデルの残差（記憶を持つ）。
  2. SIF 時系列（4-day 等）を日次へ整列し、SIF・遅延SIF(cumSIF 4/8日)を残差に足す。
  3. **SIF で残差の自己相関(記憶)が落ちるか**を測る（気象では落ちなかったものが SIF で落ちれば＝基質供給が正体）。
     ＋残差と SIF の偏相関、R² 上昇。
判定：SIF で ACF が落ちる→★基質供給が記憶の正体（分割アーティファクトでなく生物）。落ちない→SIFでも説明できず更に深い未観測。

SIF は `--sif <csv>`（列: date, sif）で渡す（取得は SIF_PIPELINE.md 参照＝GEE等でローカル抽出）。
--sif 無しなら合成で検証（隠れ基質を SIF が捉える場合/捉えない場合を分離）。

    python research/sif_respiration_step38.py                                  # 合成で検証
    python research/sif_respiration_step38.py --site JP-Tak --sif JP-Tak_sif.csv --qc-max 1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from memory_timescale_step37 import DRIVERS, _z, _autocorr, efolding_days, fit_residual


def load_flux_daily(site, months, qc_max):
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    from japanflux_pn.run_robustness import get_site_years
    cfg = AnalysisConfig(qc_max=qc_max) if qc_max is not None else AnalysisConfig()
    years, mo = get_site_years(site)
    ms = sorted(months or mo)
    raw = load_raw_all(get_site(site), cfg)
    raw = raw[raw.index.month.isin(ms)]
    keep = ["GER"] + [v for v in DRIVERS if v in raw.columns]
    return raw[keep].groupby(raw.index.normalize()).mean().dropna()


def load_sif(path):
    """SIF csv (date, sif) を日次 Series に。4/8-day は日次へ前方補間。"""
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    dc = cols.get("date") or cols.get("time") or list(df.columns)[0]
    sc = cols.get("sif") or [c for c in df.columns if c != dc][0]
    s = pd.Series(pd.to_numeric(df[sc], errors="coerce").to_numpy(),
                  index=pd.to_datetime(df[dc]))
    s = s[~s.index.duplicated()].sort_index()
    daily = s.resample("D").mean().interpolate("time", limit=8)   # 4/8日→日次
    return daily


def add_sif_terms(daily, sif_daily):
    """daily に SIF・cumSIF(4/8日, 年内)を付ける。整列できた行のみ。"""
    d = daily.copy()
    d["SIF"] = sif_daily.reindex(d.index)
    d = d.dropna(subset=["SIF"])
    if len(d) < 40:
        return d, []
    # 年内で後方窓平均（基質の蓄積を模す）
    parts = []
    for _, g in d.groupby(d.index.year):
        g = g.copy()
        g["cumSIF4"] = g["SIF"].rolling(4, min_periods=2).mean()
        g["cumSIF8"] = g["SIF"].rolling(8, min_periods=3).mean()
        parts.append(g)
    d = pd.concat(parts).dropna()
    return d, ["SIF", "cumSIF4", "cumSIF8"]


def _partial_spearman(y, x, ctrl):
    def rk(a):
        return pd.Series(np.asarray(a, float)).rank().to_numpy()
    ok = np.isfinite(y) & np.isfinite(x)
    Z = [np.asarray(c, float) for c in ctrl]
    for c in Z:
        ok &= np.isfinite(c)
    if ok.sum() < 20:
        return np.nan
    ry, rx = rk(y[ok]), rk(x[ok])
    A = np.column_stack([rk(c[ok]) for c in Z] + [np.ones(ok.sum())])
    yr = ry - A @ np.linalg.lstsq(A, ry, rcond=None)[0]
    xr = rx - A @ np.linalg.lstsq(A, rx, rcond=None)[0]
    if xr.std() == 0 or yr.std() == 0:
        return np.nan
    return float(np.corrcoef(xr, yr)[0, 1])


def analyze(daily, sif_daily):
    """気象のみ vs 気象+SIF で残差の記憶(ACF)がどう変わるか。"""
    r2_w, res_w = fit_residual(daily)          # 気象のみ
    ac_w = _autocorr(res_w); ef_w, _ = efolding_days(res_w)
    d2, sif_cols = add_sif_terms(daily, sif_daily)
    if not sif_cols:
        return {"note": "SIF 整列不足"}
    # SIF を足す前(この部分集合での気象のみ)と後
    r2_w2, res_w2 = fit_residual(d2)
    r2_s, res_s = fit_residual(d2, sif_cols)
    ac_s = _autocorr(res_s); ef_s, _ = efolding_days(res_s)
    # 残差 vs SIF の偏相関（気象を差し引いた上で SIF が呼吸残差を説明するか）
    pr = _partial_spearman(d2["GER"].to_numpy(),
                           d2["SIF"].to_numpy(),
                           [d2[v].to_numpy() for v in DRIVERS if v in d2])
    return {"n": len(d2), "r2_weather": r2_w2, "r2_sif": r2_s,
            "ac_weather": _autocorr(res_w2), "ac_sif": ac_s,
            "ef_weather": ef_w, "ef_sif": ef_s, "partial_r": pr}


def make_synth(kind, days=900, seed=0):
    """隠れ基質 S を SIF が捉える/捉えない 2 ケース。気象は S を作れない。"""
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2005-06-01", periods=days, freq="D")
    Ta = 20 + 5 * np.sin(2 * np.pi * idx.dayofyear / 365.25) + rng.normal(0, 1.5, days)
    Rg = 300 + rng.normal(0, 40, days); th = np.clip(0.3 + rng.normal(0, 0.05, days), .05, .6)
    VPD = np.clip(0.5 + 0.06 * (Ta - 15), 0.1, None); P = np.clip(rng.normal(0, 1, days), 0, None)
    # 隠れ基質 S＝気象と独立な AR(0.8)（最近の光合成の蓄積を模す, 約4-5日スケール）
    S = np.zeros(days)
    for i in range(1, days):
        S[i] = 0.8 * S[i - 1] + rng.normal(0, 0.5)
    S = _z(S)
    GER = 2.0 * np.exp(0.06 * (Ta - 20)) + 1.5 * S + rng.normal(0, 0.15, days)
    flux = pd.DataFrame({"GER": np.clip(GER, 1e-3, None), "Rg": Rg, "Ta": Ta, "VPD": VPD,
                         "Ts": Ta - 1, "th": th, "P": P, "gH": Rg * .2, "gLE": Rg * .3}, index=idx)
    if kind == "sif_captures":       # SIF が基質 S を捉える
        sif = pd.Series(S + rng.normal(0, 0.3, days), index=idx)
    else:                             # SIF が無関係（対照）
        sif = pd.Series(rng.normal(0, 1, days), index=idx)
    return flux, sif


def _report(r, tag):
    if "note" in r:
        print(f"  {tag}: {r['note']}"); return
    print(f"\n  === {tag}（N={r['n']}）===")
    print(f"  気象のみ  : R²={r['r2_weather']:.3f}  残差ACF={r['ac_weather']:+.2f}  e-fold={r['ef_weather']}")
    print(f"  ＋SIF     : R²={r['r2_sif']:.3f}  残差ACF={r['ac_sif']:+.2f}  e-fold={r['ef_sif']}"
          f"（ACF低下={r['ac_weather']-r['ac_sif']:+.2f}）")
    pr = r["partial_r"]
    print(f"  残差 vs SIF 偏相関（気象差引後）: r={pr:+.2f}" if np.isfinite(pr) else "  偏相関: —")
    drop = r["ac_weather"] - r["ac_sif"]
    if np.isfinite(drop) and r["ac_sif"] < 0.3 and drop >= 0.2:
        print("  → ★ SIF で残差の記憶が消えた＝**基質供給(最近の光合成)が記憶の正体**")
        print("     ＝独立な衛星シグナルが説明＝分割窓アーティファクトでなく生物（旗37の交絡を破った）")
    elif np.isfinite(drop) and drop >= 0.15:
        print("  → ○ SIF で一部説明＝基質供給が寄与するが残りは別の未観測")
    else:
        print("  → ・ SIF でも記憶が落ちない＝基質供給では説明できず更に深い未観測（or SIF整列の問題）")


def main():
    p = argparse.ArgumentParser(description="SIFは呼吸残差の記憶を説明するか")
    p.add_argument("--site")
    p.add_argument("--sif", help="SIF csv（列: date, sif）")
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--qc-max", type=int, default=None)
    a = p.parse_args()

    if not a.site:
        print("=== 旗38 合成検証：SIF が隠れ基質を捉えれば残差の記憶を消せるか ===")
        for kind, lab in [("sif_captures", "SIF が隠れ基質 S を捉える"),
                          ("noise", "SIF が無関係（対照）")]:
            flux, sif = make_synth(kind)
            _report(analyze(flux, sif), lab)
        print("\n  → 捉える場合は SIF で ACF が急落・偏相関が高い、無関係な対照は落ちないのが期待。")
        return

    if not a.sif:
        print("実データには --sif <csv> が必要（SIF_PIPELINE.md で取得）。"); return
    daily = load_flux_daily(a.site, a.month, a.qc_max)
    sif = load_sif(a.sif)
    print(f"=== 旗38 実データ {a.site}（SIFは呼吸残差の記憶を説明するか, QC≤{a.qc_max}）===")
    _report(analyze(daily, sif), f"{a.site}")
    print("\n  読み方：SIFは分割に依らない独立な光合成シグナル。SIFで記憶が落ちれば＝基質供給が本物で")
    print("    旗37の『生物か分割窓か』の交絡を破れる。落ちなければ更に深い未観測（土壌/微生物/深水分）。")
    print("  留保：SIFピクセルとタワーfootprintの空間不一致・4/8日→日次補間・GERは分割派生量。")


if __name__ == "__main__":
    main()
