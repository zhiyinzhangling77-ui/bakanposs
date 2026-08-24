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


def sd_coords(gH, gLE):
    """アノマリでも安全な直交座標：S=総量(γH+γLE), D=配分コントラスト(γH−γLE)。
    比 β=γH/(γH+γLE) はアノマリだと分母がゼロを跨いで発散するので**使わない**（旗50 第1回の失敗）。"""
    return gH + gLE, gH - gLE


def surrogate(gH, gLE, rng, shuffle="D"):
    """S,D の一方だけを時間シャッフルする。割り算を含まないので発散しない。
    shuffle='D'：総量Sを保存し**配分だけ壊す** ／ 'S'：配分Dを保存し**総量だけ壊す**。"""
    S, D = sd_coords(gH, gLE)
    if shuffle == "D":
        D = D[rng.permutation(len(D))]
    else:
        S = S[rng.permutation(len(S))]
    return (S + D) / 2.0, (S - D) / 2.0


def analyze(Rg, Ta, gH, gLE, nrep=20, seed=0):
    ok = np.isfinite(Rg) & np.isfinite(Ta) & np.isfinite(gH) & np.isfinite(gLE)
    Rg, Ta, gH, gLE = Rg[ok], Ta[ok], gH[ok], gLE[ok]
    if len(Rg) < 500:
        return {"note": f"点不足({len(Rg)})"}
    rng = np.random.default_rng(seed)
    obs = _omega([Rg, Ta, gH, gLE])
    surD, surS = [], []
    for _ in range(nrep):
        h, e = surrogate(gH, gLE, rng, "D")      # 総量は保存・配分だけ壊す
        surD.append(_omega([Rg, Ta, h, e]))
        h, e = surrogate(gH, gLE, rng, "S")      # 配分は保存・総量だけ壊す
        surS.append(_omega([Rg, Ta, h, e]))
    S, D = sd_coords(gH, gLE)
    reparam = _omega([Rg, Ta, S, D])
    return {"n": int(len(Rg)), "obs": obs,
            "sur": float(np.mean(surD)), "sur_sd": float(np.std(surD)),
            "surS": float(np.mean(surS)), "surS_sd": float(np.std(surS)),
            "reparam": reparam}


def rho(r):
    """ρ＝Δ配分／Δ総量。配分を壊した損失が総量を壊した損失に比べてどれだけ大きいか。
    **閾値で切らず、合成の2基準点（配分が無情報／情報あり）と比べて読む**（旗50の較正）。"""
    if "note" in r:
        return np.nan
    dD, dS = r["obs"] - r["sur"], r["obs"] - r["surS"]
    return dD / dS if dS > 0 else np.nan


def verdict(r, ref_lo=None, ref_hi=None):
    if "note" in r:
        return "―" + r["note"]
    q = rho(r)
    if not np.isfinite(q):
        return "△判定不能"
    if ref_lo is None or ref_hi is None:
        return f"ρ={q:.2f}"
    mid = (ref_lo + ref_hi) / 2
    if q < ref_lo * 1.2:
        return f"▲収支の言い換え寄り（ρ={q:.2f} ≈ 無情報基準 {ref_lo:.2f}）"
    if q > mid:
        return f"★配分が情報を持つ（ρ={q:.2f} ≳ 情報あり基準 {ref_hi:.2f}）"
    return f"○中間（ρ={q:.2f}, 基準 {ref_lo:.2f}〜{ref_hi:.2f}）"


def synth_references(nrep=8):
    """合成の2基準点 ρ を返す：(配分が無情報, 配分が情報を持つ)。実データの物差しにする。"""
    out = []
    for kind in ("beta_uninformative", "beta_informative"):
        qs = []
        for k in range(nrep):
            Rg, Ta, gH, gLE = _synth(kind, seed=k)
            qs.append(rho(analyze(Rg, Ta, gH, gLE, nrep=6, seed=k)))
        out.append(float(np.nanmean(qs)))
    return out[0], out[1]


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
    print(f"  {'仕込み':<24}{'Ω(観測)':>10}{'配分壊す':>10}{'総量壊す':>10}{'ρ=Δ配分/Δ総量':>16}")
    for kind, lab in [("beta_uninformative", "配分は無情報（収支のみ）"),
                      ("beta_informative", "配分が駆動と結ぶ")]:
        Rg, Ta, gH, gLE = _synth(kind)
        r = analyze(Rg, Ta, gH, gLE)
        print(f"  {lab:<24}{r['obs']:>10.4f}{r['sur']:>10.4f}{r['surS']:>10.4f}{rho(r):>16.3f}")
    lo, hi = synth_references()
    print(f"\n  基準点（8シード平均）：配分が無情報 ρ={lo:.2f} ／ 配分が情報を持つ ρ={hi:.2f}")
    print("  → 2つが十分に離れていれば、実データの ρ をこの物差しで読める。")
    print("  ※閾値による二値判定は諦めた：どんなサロゲートも乗法結合を壊すため、")
    print("    『配分が無情報』でも Δ配分 は完全にはゼロにならない（第2回の設計判断）。")


# ---------- 実データ ---------------------------------------------------------------
def run_real(sites, months):
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.run_robustness import get_site_years
    from japanflux_pn.preprocess import (load_raw_all, slice_and_anomaly,
                                         slice_span_and_anomaly)
    cfg = AnalysisConfig()
    print("=== 旗50 実データ：エネルギー背骨の冗長は収支の言い換えか ===")
    print("  合成で物差しを作ってから読む（数十秒）…", flush=True)
    ref_lo, ref_hi = synth_references()
    print(f"  基準点：配分が無情報 ρ={ref_lo:.2f} ／ 配分が情報を持つ ρ={ref_hi:.2f}")
    print("  ρ＝Δ配分/Δ総量。ρ が無情報基準に近い＝背骨は収支の言い換え。\n")
    print(f"  {'site-year':<16}{'N':>7}{'Ω(観測)':>10}{'配分壊す':>10}{'総量壊す':>10}"
          f"{'Ω(Rg,Ta,S,D)':>13}  判定")
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
            print(f"  {site+' '+str(y):<16}{r['n']:>7}{r['obs']:>10.4f}{r['sur']:>10.4f}"
                  f"{r['surS']:>10.4f}{r['reparam']:>13.4f}  {verdict(r, ref_lo, ref_hi)}", flush=True)
    if not rows:
        print("  有効な site-year なし"); return
    obs = np.mean([r["obs"] for r in rows]); sur = np.mean([r["sur"] for r in rows])
    surS = np.mean([r["surS"] for r in rows]); rep = np.mean([r["reparam"] for r in rows])
    n_id = sum(1 for r in rows if "▲" in verdict(r, ref_lo, ref_hi))
    qs = np.array([rho(r) for r in rows], float)
    print(f"\n  === まとめ（{len(rows)} site-year 平均）===")
    print(f"  Ω(観測) {obs:.4f}")
    print(f"  配分Dを壊す {sur:.4f}（差 {obs-sur:+.4f}）／ 総量Sを壊す {surS:.4f}（差 {obs-surS:+.4f}）")
    print(f"  ρ 中央値 {np.nanmedian(qs):.2f}（基準：無情報 {ref_lo:.2f} ／ 情報あり {ref_hi:.2f}）")
    print(f"  『収支の言い換え寄り』判定：{n_id}/{len(rows)} site-year")
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
