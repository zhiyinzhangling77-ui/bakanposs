"""旗42：水分依存Q10(旗26/27)をチャンバー呼吸(COSORE)で直接検証＝穴④を叩く。

旗26/27で「湿るほど呼吸の温度感度Q10が上がる」は最厳QC・DT/NT分割で格下げされ、頑健なのは
JP-Fhk/Fjyのみ＝分割派生GER依存の脆い発見だった。**チャンバーRsは分割を通さない直接測定**なので、
同じ解析(θビンごとに ln(Rs)=a+b·Tsoil → Q10=exp(10b)、Q10 vs θ の Spearman)を直接測定でやる：
  ・多数の森林で Q10 が土壌水分で上がる(r>0,CI>0) → 水分依存Q10は本物(分割由来でない)。
  ・出ない/割れる → 旗26の水分依存Q10は分割アルゴリズム由来 or 特定サイト固有だった。

旗26 の q10_by_moisture / _boot_trend を再利用。COSORE の Rs=CSR_FLUX_CO2・土壌温度・水分を使う。

    python research/cosore_q10_step42.py                                   # 合成で検証
    python research/cosore_q10_step42.py --cosore-dir /mnt/hdd/cosore-0.7.0 # 全森林で一括
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from respiration_q10_moisture_step26 import q10_by_moisture, _boot_trend, make_synth
from cosore_memory_step40 import load_cosore


def analyze_file(path, months, nbin=5):
    df, st, sm = load_cosore(path, months)
    if "Tsoil" not in df or "SM" not in df:
        return {"note": "土壌温度/水分なし"}
    Ta = df["Tsoil"].to_numpy(); th = df["SM"].to_numpy(); Rs = df["Rs"].to_numpy()
    ok = np.isfinite(Ta) & np.isfinite(th) & np.isfinite(Rs) & (Rs > 0)
    if ok.sum() < 1000:
        return {"note": f"点不足({int(ok.sum())})"}
    r, ci, fr = _boot_trend(Ta[ok], th[ok], Rs[ok], nbin, 200)
    rows = q10_by_moisture(Ta[ok], th[ok], Rs[ok], nbin)
    return {"n": int(ok.sum()), "r": r, "ci": ci, "st": st, "sm": sm,
            "q10lo": rows[0]["Q10"] if rows else np.nan,
            "q10hi": rows[-1]["Q10"] if rows else np.nan}


def _verdict(ci):
    if not isinstance(ci, tuple):
        return "—"
    if ci[0] > 0:
        return "★水分依存Q10(正)"
    if ci[1] < 0:
        return "×逆(乾で高Q10)"
    return "△CI0跨ぎ"


def run_batch(cosore_dir, months, igbp_filter="forest"):
    import pandas as pd
    root = Path(cosore_dir)
    desc = pd.read_csv(root / "description.csv")
    print(f"=== 旗42 バッチ：COSORE チャンバーRsの水分依存Q10（{igbp_filter or '全'}, 月={months or '全'}）===")
    print(f"  {'dataset':<32}{'IGBP':<13} {'N':>7} {'Q10乾→湿':>10} {'r(Q10 vs θ)':>11} {'95%CI':>15}  判定")
    rows = []
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); igbp = str(d.get("CSR_IGBP", ""))
        if igbp_filter and igbp_filter.lower() not in igbp.lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            r = analyze_file(f, months)
        except Exception as e:
            print(f"  {ds:<32}{igbp[:11]:<13} SKIP {type(e).__name__}"); continue
        if "note" in r:
            continue                       # 土壌欠測/点不足は静かに飛ばす
        v = _verdict(r["ci"])
        cistr = f"[{r['ci'][0]:+.2f},{r['ci'][1]:+.2f}]" if isinstance(r["ci"], tuple) else "—"
        rows.append({"ds": ds, "r": r["r"], "ci": r["ci"], "v": v})
        print(f"  {ds:<32}{igbp[:11]:<13} {r['n']:>7} {r['q10lo']:>4.2f}→{r['q10hi']:<4.2f} "
              f"{r['r']:>+11.2f} {cistr:>15}  {v}")
    if not rows:
        print("  判定可能サイトなし（土壌温度・水分が要る）"); return
    pos = sum(1 for x in rows if isinstance(x["ci"], tuple) and x["ci"][0] > 0)
    neg = sum(1 for x in rows if isinstance(x["ci"], tuple) and x["ci"][1] < 0)
    print(f"\n  === まとめ（判定可能 n={len(rows)}）===")
    print(f"  ★水分依存Q10(正,CI>0)：{pos}／×逆(乾で高Q10)：{neg}／△CI0跨ぎ：{len(rows)-pos-neg}")
    print("  読み方：分割を通さない直接測定の多数森林で★が多数なら＝水分依存Q10は本物(旗26の分割由来説を否定)。")
    print("    割れる/少ないなら＝旗26の水分依存Q10は分割アルゴリズム由来 or 特定サイト固有だった。留保：チャンバーは点測定・土壌呼吸のみ。")


def main():
    p = argparse.ArgumentParser(description="水分依存Q10をチャンバー直接測定で検証")
    p.add_argument("--cosore-dir")
    p.add_argument("--igbp", default="forest")
    p.add_argument("--month", type=int, nargs="+", default=None)
    a = p.parse_args()

    if a.cosore_dir:
        run_batch(a.cosore_dir, a.month, igbp_filter=(a.igbp or None)); return

    print("=== 旗42 合成検証：Q10 vs θ を直接測定ロジックで検出できるか ===")
    for kind, lab in [("const", "Q10一定(水分非依存)"), ("wet_up", "湿るほどQ10大")]:
        Ta, th, GER = make_synth(kind)
        r, ci, fr = _boot_trend(Ta, th, GER, 5, 200, nboot=300)
        print(f"  {lab:<20} r={r:+.2f} CI={ci if isinstance(ci,tuple) else '—'}  {_verdict(ci)}")
    print("\n  → 一定は r≈0(CI0跨ぎ)、湿るほど大は r>0・CI>0(★) が期待。")


if __name__ == "__main__":
    main()
