"""旗44：水分依存Q10(旗26/27/42)の『温度エイリアシング』交絡を叩く。

前提の穴：Q10 は温度そのものの関数（Lloyd & Taylor 1994＝低温ほど見かけQ10が大きい）。
旗26/42 の q10_by_moisture は θビンごとの平均温度を制御していない。温帯・冷温帯林では
湿ったビン＝春秋(低温)／乾いたビン＝真夏(高温)にエイリアスするので、**水分が感度に一切効かなくても
Lloyd-Taylor曲率だけで「湿→高Q10」が出る**（本ファイルの合成でr=+0.90〜1.00を再現＝旗42と同じ値）。
旗42のパターン(温帯冷温帯で正・熱帯泥炭と凍土で逆・Fairbanks Q10=25という当てはめ崩壊)は
この交絡の予測とも一致してしまう＝旗42はまだ温度交絡を除いていない。

決着法：温度の曲率を明示的に吸収した上で「水分が"温度感度"を変えるか」を交互作用係数 d で測る。
    ln R = a + b·Tc + e·Tc² + c·θz + d·(Tc·θz)          (Tc=T−平均T, θz=θの標準化)
  ・e·Tc² が Lloyd-Taylor 曲率（＝見かけQ10の温度依存）を吸収する。
  ・d>0 (CI>0) なら「湿るほど温度感度が大きい」が曲率を除いてなお残る＝**水分依存Q10は本物**。
  ・曲率なしモデルでは d>0 なのに曲率ありで消えるなら＝**温度エイリアシングの産物**。
平均温度での Q10 を θz=∓1SD について exp(10(b∓d)) で示す（旗26/42のQ10と直接比較できる形）。

    python research/q10_confound_step44.py                                    # 合成で検証
    python research/q10_confound_step44.py --cosore-dir /mnt/hdd/cosore-0.7.0 # 実データ(森林)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore

T0_LT = -46.02          # Lloyd-Taylor の下限温度 (°C)


def _design(Tc, thz, quad):
    cols = [Tc, thz, Tc * thz, np.ones_like(Tc)]
    if quad:
        cols.insert(1, Tc ** 2)
    return np.column_stack(cols)


def _fit_d(Tc, thz, lnR, quad):
    """交互作用係数 d（＝水分が温度感度を変える強さ）と温度感度 b を返す。"""
    A = _design(Tc, thz, quad)
    coef = np.linalg.lstsq(A, lnR, rcond=None)[0]
    b = coef[0]
    d = coef[3] if quad else coef[2]
    return float(b), float(d)


def _boot_d(Tc, thz, lnR, quad, nboot=200, seed=0, block=48 * 7):
    """ブロックブート（自己相関を潰さないよう連続ブロックで再標本）。"""
    rng = np.random.default_rng(seed); n = len(Tc); ds = []
    nblk = max(1, n // block)
    for _ in range(nboot):
        starts = rng.integers(0, max(1, n - block), nblk)
        idx = np.concatenate([np.arange(s, min(s + block, n)) for s in starts])
        try:
            ds.append(_fit_d(Tc[idx], thz[idx], lnR[idx], quad)[1])
        except np.linalg.LinAlgError:
            continue
    if len(ds) < 30:
        return None
    return (float(np.percentile(ds, 2.5)), float(np.percentile(ds, 97.5)))


def analyze(T, th, R):
    ok = np.isfinite(T) & np.isfinite(th) & np.isfinite(R) & (R > 0)
    T, th, R = T[ok], th[ok], R[ok]
    if len(T) < 1000 or T.max() - T.min() < 5 or th.std() == 0:
        return {"note": "点不足/温度レンジ不足"}
    Tc = T - T.mean(); thz = (th - th.mean()) / th.std(); lnR = np.log(R)
    b_nq, d_nq = _fit_d(Tc, thz, lnR, quad=False)     # 曲率を入れない（＝旗26/42と同じ土俵）
    b_q, d_q = _fit_d(Tc, thz, lnR, quad=True)        # Lloyd-Taylor曲率を吸収
    ci = _boot_d(Tc, thz, lnR, quad=True)
    return {"n": int(len(T)), "d_noquad": d_nq, "d_quad": d_q, "ci": ci,
            "q10_dry": float(np.exp(10 * (b_q - d_q))), "q10_wet": float(np.exp(10 * (b_q + d_q))),
            "corr_thT": float(np.corrcoef(th, T)[0, 1])}


def verdict(res):
    if "note" in res:
        return "―" + res["note"]
    ci = res["ci"]
    if ci is None:
        return "△CI不定"
    if ci[0] > 0:
        return "★曲率制御後も水分依存Q10(本物)"
    if ci[1] < 0:
        return "×曲率制御後は逆(乾で高感度)"
    if res["d_noquad"] > 0:
        return "▲曲率なしで正→制御で消失=温度エイリアシングの産物"
    return "△制御後CI0跨ぎ"


def _synth(kind, n=60000, seed=0):
    """Lloyd-Taylor(低温ほど見かけQ10大) × θとTのエイリアス(湿=低温)。kind=alias は水分が感度に効かない。"""
    rng = np.random.default_rng(seed)
    doy = rng.uniform(0, 365, n)
    T = 12 + 12 * np.sin(2 * np.pi * (doy - 100) / 365) + rng.normal(0, 2, n)
    th = 0.30 - 0.05 * np.sin(2 * np.pi * (doy - 100) / 365) + rng.normal(0, 0.07, n)
    lnR = 320.0 * (1.0 / (10 - T0_LT) - 1.0 / (T - T0_LT))   # Lloyd-Taylor（Q10が温度の関数）
    lnR += 0.8 * (th - 0.30) / 0.10                           # 水分は"量"を変えるが感度は変えない
    if kind == "true":
        lnR += 0.6 * ((th - 0.30) / 0.10) * (T - 12) / 10     # 真に感度を変える
    return T, th, np.exp(lnR + rng.normal(0, 0.15, n))


def main():
    p = argparse.ArgumentParser(description="水分依存Q10の温度エイリアシング交絡を叩く")
    p.add_argument("--cosore-dir"); p.add_argument("--igbp", default="forest")
    p.add_argument("--month", type=int, nargs="+", default=None)
    a = p.parse_args()

    if not a.cosore_dir:
        print("=== 旗44 合成検証：温度交絡だけで『湿→高Q10』が出るか ===")
        for kind, lab in [("alias", "水分は感度に効かない(温度エイリアスのみ)"),
                          ("true", "真に水分が感度を変える")]:
            T, th, R = _synth(kind)
            r = analyze(T, th, R)
            print(f"  {lab}")
            print(f"    θ-T相関={r['corr_thT']:+.2f}  d(曲率なし)={r['d_noquad']:+.4f}  "
                  f"d(曲率あり)={r['d_quad']:+.4f} CI={r['ci']}")
            print(f"    平均温度でのQ10：乾(−1SD) {r['q10_dry']:.2f} → 湿(+1SD) {r['q10_wet']:.2f}")
            print(f"    → {verdict(r)}\n")
        print("  期待：上=曲率なしで正でも曲率制御で消える(▲)、下=制御後も正(★)。この差が付けば検出器は妥当。")
        return

    import pandas as pd
    root = Path(a.cosore_dir); desc = pd.read_csv(root / "description.csv")
    print(f"=== 旗44 実データ：水分依存Q10は温度エイリアシングか（{a.igbp}）===")
    print(f"  {'dataset':<32} {'θ-T':>6} {'d(曲率なし)':>11} {'d(曲率あり)':>11} {'95%CI':>17} "
          f"{'Q10乾→湿':>11}  判定")
    tally = {}
    for _, dd in desc.iterrows():
        ds = str(dd["CSR_DATASET"]); igbp = str(dd.get("CSR_IGBP", ""))
        if a.igbp and a.igbp.lower() not in igbp.lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            df, st, sm = load_cosore(f, a.month)
            if "Tsoil" not in df or "SM" not in df:
                continue
            res = analyze(df["Tsoil"].to_numpy(), df["SM"].to_numpy(), df["Rs"].to_numpy())
        except Exception as e:
            print(f"  {ds:<32} SKIP {type(e).__name__}"); continue
        if "note" in res:
            continue
        v = verdict(res); key = v.split("(")[0][:14]
        tally[key] = tally.get(key, 0) + 1
        ci = res["ci"]
        s_ci = f"[{ci[0]:+.4f},{ci[1]:+.4f}]" if ci else "—"
        print(f"  {ds:<32} {res['corr_thT']:>+6.2f} {res['d_noquad']:>+11.4f} {res['d_quad']:>+11.4f} "
              f"{s_ci:>17} {res['q10_dry']:>5.2f}→{res['q10_wet']:<5.2f}  {v}")
    print("\n  === まとめ ===")
    for k, v in sorted(tally.items(), key=lambda x: -x[1]):
        print(f"    {k:<18} {v}")
    print("  読み：★が多数＝水分依存Q10は温度曲率でなく本物(旗42を支持)。")
    print("        ▲が多数＝旗42の『湿→高Q10』は温度エイリアシングの産物＝正直に格下げが必要。")
    print("  注：dは ln単位/(℃·θ1SD)。ブロックブート(7日)で自己相関に配慮。")


if __name__ == "__main__":
    main()
