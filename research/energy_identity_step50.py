"""旗50（前提監査③）：エネルギー背骨の冗長は「地表エネルギー収支の言い換え」か。

前提の穴：層①の中核は「エネルギー系 {Rg,Ta,γH,γLE} が最強の冗長」だが、地表エネルギー収支
H+LE ≈ Rn−G が制約として存在し、Rn は Rg の関数である以上、**この冗長は物理的にほぼ自明**かもしれない。
旗32 で VPD/GER/GEP の派生は監査したのに、**収支制約そのものは監査していない**。旗46 でギャップフィル
由来でないことは確かめたが、「実測でも立つ」ことと「発見である」ことは別問題である。

決着法：**制約だけを保存したサロゲート**と比べる。各時刻で
    総乱流フラックス T = γH + γLE  … そのまま保存（＝利用可能エネルギーとの結びつきを残す）
    ボーエン配分 β = γH / T        … **時間方向にシャッフル**（＝配分が持つ情報だけを壊す）
から γH* = β_shuffled·T, γLE* = (1−β_shuffled)·T を作り、Ω({Rg,Ta,γH*,γLE*}) を測る。
  ・Ω(サロゲート) ≈ Ω(観測) → 冗長は**「H と LE が同じ利用可能エネルギーを分け合う」ことの言い換え**
    ＝背骨は発見でなく収支の再記述。主張をそう書き直す。
  ・Ω(観測) > Ω(サロゲート) → **配分の仕方それ自体が情報を持つ**＝収支だけでは説明できない構造がある。

併せて座標を張り替えた Ω({Rg,Ta,T,β}) も出す（「どれだけ」と「どう分けるか」に分解＝
③が言う『意味があるのは収支からのズレの方』を直接見る）。

    python research/energy_identity_step50.py                       # 合成で検証
    python research/energy_identity_step50.py --sites JP-Tak JP-Fhk
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

BINS = 8


def _omega(cols_raw, bins=BINS):
    from japanflux_pn import information_theory as it
    cols = [it.digitize_series(np.asarray(c, float), bins) for c in cols_raw]
    return float(it.o_information_indices(cols, bins, correct=True))


def partition_surrogate(gH, gLE, rng):
    """総乱流フラックスは保存し、ボーエン配分だけを時間シャッフルしたサロゲート。"""
    T = gH + gLE
    ok = np.abs(T) > 1e-9
    beta = np.full_like(T, np.nan)
    beta[ok] = gH[ok] / T[ok]
    idx = np.flatnonzero(np.isfinite(beta))
    b_sh = beta.copy()
    b_sh[idx] = beta[idx][rng.permutation(len(idx))]
    return b_sh * T, (1 - b_sh) * T, T, beta


def analyze(Rg, Ta, gH, gLE, nrep=20, seed=0):
    ok = np.isfinite(Rg) & np.isfinite(Ta) & np.isfinite(gH) & np.isfinite(gLE)
    Rg, Ta, gH, gLE = Rg[ok], Ta[ok], gH[ok], gLE[ok]
    if len(Rg) < 500:
        return {"note": f"点不足({len(Rg)})"}
    rng = np.random.default_rng(seed)
    obs = _omega([Rg, Ta, gH, gLE])
    sur = []
    for _ in range(nrep):
        h, e, T, beta = partition_surrogate(gH, gLE, rng)
        m = np.isfinite(h) & np.isfinite(e)
        if m.sum() > 500:
            sur.append(_omega([Rg[m], Ta[m], h[m], e[m]]))
    T = gH + gLE
    with np.errstate(invalid="ignore", divide="ignore"):
        beta = np.where(np.abs(T) > 1e-9, gH / T, np.nan)
    m = np.isfinite(beta)
    reparam = _omega([Rg[m], Ta[m], T[m], beta[m]]) if m.sum() > 500 else np.nan
    return {"n": int(len(Rg)), "obs": obs,
            "sur": float(np.mean(sur)) if sur else np.nan,
            "sur_sd": float(np.std(sur)) if sur else np.nan,
            "reparam": reparam}


def verdict(r):
    if "note" in r:
        return "―" + r["note"]
    if not np.isfinite(r["sur"]):
        return "△サロゲート不能"
    d = r["obs"] - r["sur"]
    rel = abs(d) / abs(r["sur"]) if r["sur"] else np.nan
    z = d / r["sur_sd"] if r["sur_sd"] > 0 else np.nan
    if np.isfinite(rel) and rel < 0.10:
        return f"▲収支の言い換え（差 {d:+.4f}, {rel:.0%}）"
    if d > 0:
        return f"★配分自体が情報を持つ（差 {d:+.4f}, z={z:+.1f}）"
    return f"○観測の方が低い（差 {d:+.4f}, z={z:+.1f}）"


# ---------- 合成 -------------------------------------------------------------------
def _synth(kind, n=6000, seed=0):
    """T(総乱流)は放射に従う。kind で『配分βが情報を持つか』を変える。"""
    rng = np.random.default_rng(seed)
    Rg = np.abs(rng.normal(0, 1, n))
    Ta = 0.7 * Rg + rng.normal(0, 0.6, n)
    T = 2.0 * Rg + 0.3 * Ta + np.abs(rng.normal(0, 0.3, n))   # 利用可能エネルギーに従う総量
    if kind == "beta_informative":
        beta = 1 / (1 + np.exp(-(0.9 * Rg - 0.6 * Ta + rng.normal(0, 0.3, n))))
    else:                                                      # 配分は無情報（収支だけ）
        beta = rng.uniform(0.2, 0.8, n)
    return Rg, Ta, beta * T, (1 - beta) * T


def run_synth():
    print("=== 旗50 合成検証：配分βが情報を持つ場合／持たない場合を見分けられるか ===")
    print("  どちらも総乱流フラックス T は放射に従う（＝収支の制約は共通）。\n")
    print(f"  {'仕込み':<24}{'Ω(観測)':>10}{'Ω(配分シャッフル)':>18}{'Ω(Rg,Ta,T,β)':>14}  判定")
    for kind, lab in [("beta_uninformative", "配分は無情報（収支のみ）"),
                      ("beta_informative", "配分が駆動と結ぶ")]:
        Rg, Ta, gH, gLE = _synth(kind)
        r = analyze(Rg, Ta, gH, gLE)
        print(f"  {lab:<24}{r['obs']:>10.4f}{r['sur']:>18.4f}{r['reparam']:>14.4f}  {verdict(r)}")
    print("\n  → 上が▲（収支の言い換え）、下が★（配分自体が情報）と出れば検出器は妥当。")


# ---------- 実データ ---------------------------------------------------------------
def run_real(sites, months):
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.run_robustness import get_site_years
    from japanflux_pn.preprocess import (load_raw_all, slice_and_anomaly,
                                         slice_span_and_anomaly)
    cfg = AnalysisConfig()
    print("=== 旗50 実データ：エネルギー背骨の冗長は収支の言い換えか ===")
    print("  Ω>0=冗長。観測 と『総量は保存・配分だけシャッフル』を比べる。\n")
    print(f"  {'site-year':<16}{'N':>7}{'Ω(観測)':>10}{'Ω(配分シャッフル)':>18}"
          f"{'Ω(Rg,Ta,T,β)':>14}  判定")
    rows = []
    for site in sites:
        years, mons = get_site_years(site)
        if months:
            mons = months
        try:
            raw = load_raw_all(get_site(site), cfg)
        except Exception as e:
            print(f"  {site}: 読み込み失敗 {type(e).__name__}"); continue
        for y in years:
            try:
                if len(mons) == 1:
                    an, va = slice_and_anomaly(raw, y, mons[0], cfg)
                else:
                    an, va = slice_span_and_anomaly(raw, y, mons, cfg)
            except Exception:
                continue
            d = an.loc[va]
            if len(d) < 500:
                continue
            r = analyze(d["Rg"].to_numpy(), d["Ta"].to_numpy(),
                        d["gH"].to_numpy(), d["gLE"].to_numpy())
            if "note" in r:
                continue
            rows.append(r)
            print(f"  {site+' '+str(y):<16}{r['n']:>7}{r['obs']:>10.4f}"
                  f"{r['sur']:>18.4f}{r['reparam']:>14.4f}  {verdict(r)}", flush=True)
    if not rows:
        print("  有効な site-year なし"); return
    obs = np.mean([r["obs"] for r in rows]); sur = np.mean([r["sur"] for r in rows])
    rep = np.mean([r["reparam"] for r in rows])
    n_id = sum(1 for r in rows if "▲" in verdict(r))
    print(f"\n  === まとめ（{len(rows)} site-year 平均）===")
    print(f"  Ω(観測) {obs:.4f} ／ Ω(配分シャッフル) {sur:.4f} ／ 差 {obs-sur:+.4f}"
          f" ({abs(obs-sur)/abs(sur):.0%})")
    print(f"  『収支の言い換え』判定：{n_id}/{len(rows)} site-year")
    print(f"  座標張替え Ω(Rg,Ta,総量,配分) の平均 {rep:.4f}")
    print("\n  読み：差がほぼ無い＝**背骨の冗長は『H と LE が同じ利用可能エネルギーを分け合う』ことの")
    print("        言い換え**＝発見でなく収支の再記述。主張をそう書き直す（前提③の対処）。")
    print("        差が大きい＝配分の仕方それ自体が放射・温度と結んでいる＝収支だけでは説明できない構造。")


def main():
    p = argparse.ArgumentParser(description="背骨の冗長が収支の言い換えかを判定")
    p.add_argument("--sites", nargs="+")
    p.add_argument("--month", type=int, nargs="+", default=None)
    a = p.parse_args()
    if a.sites:
        run_real(a.sites, a.month)
    else:
        run_synth()


if __name__ == "__main__":
    main()
