"""旗22：光合成の"掛け算モデル(光利用効率)"は破れているか？呼吸と直接対比。

呼吸は掛け算 R=f(Ta)g(θ) を有意に破った（旗20, JP-Tak）。光合成の標準モデルも掛け算：
  GPP = ε · Rg · f(VPD)（光利用効率 LUE。光×環境スカラーの分離可能形）。
これが破れるか＝**放射の光利用が VPD で非加法に変わるか**を、GEP(Rg,VPD) の応答曲面で測る。
旗13（乾いた年ほど I(Rg;GEP) が下がる＝放射↔光合成の脱結合）の"モデル形での現れ"になりうる。

旗20 の検証済み関数（interaction_fraction / surrogate_pvalue）をそのまま再利用（総称的）。
同じ実行で呼吸 GER(Ta,θ) も出して**光合成 vs 呼吸のモデル破れを対比**する。

    python research/gep_interaction_step22.py                     # 合成で検証
    python research/gep_interaction_step22.py --site JP-Tak --years 1999 ... --month 7 8 --deyear --nperm 500
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from interaction_surface_step20 import interaction_fraction, surrogate_pvalue

# (名前, 駆動1, 駆動2, 目標)。炭素×2＋水（エネルギー）で三輪の標準モデル形を対比。
#  光合成=光利用効率 GPP=ε·Rg·f(VPD)／呼吸=R=f(Ta)g(θ)／蒸発=Penman-Monteith/PT: LE~f(Rg,VPD)。
PAIRS = [("光合成 GEP(Rg,VPD)", "Rg", "VPD", "GEP"),
         ("呼吸   GER(Ta,θ)",   "Ta", "th", "GER"),
         ("蒸発   γLE(Rg,VPD)", "Rg", "VPD", "gLE")]


def _deyear(g, years_idx, do):
    if not do:
        return g
    out = g.copy()
    for idx in years_idx:
        gp = g[idx][np.isfinite(g[idx]) & (g[idx] > 0)]
        gm = np.exp(np.mean(np.log(gp))) if gp.size else 1.0
        out[idx] = g[idx] / gm
    return out


def run_site(site, months, deyear, nbins, min_cell, nperm):
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    from japanflux_pn.run_robustness import get_site_years
    import pandas as pd
    cfg = AnalysisConfig()
    years, mo = get_site_years(site)
    ms = sorted(months or mo)
    raw_all = load_raw_all(get_site(site), cfg)
    need = sorted({d for _, d1, d2, t in PAIRS for d in (d1, d2, t)})
    cols = {v: [] for v in need}
    yr_slices = []
    pos = 0
    used = 0
    for y in years:
        start = pd.Timestamp(year=y, month=ms[0], day=1)
        end = pd.Timestamp(year=y, month=ms[-1], day=1) + pd.offsets.MonthBegin(1)
        r = raw_all[(raw_all.index >= start) & (raw_all.index < end)]
        if r.empty:
            continue
        n = len(r)
        for v in need:
            cols[v].append(r[v].to_numpy(float))
        yr_slices.append(np.arange(pos, pos + n)); pos += n; used += 1
    if not used:
        return None, 0
    data = {v: np.concatenate(cols[v]) for v in need}
    rows = []
    for name, d1, d2, tgt in PAIRS:
        z = _deyear(data[tgt], yr_slices, deyear)
        if nperm > 0:
            frac, null, pv = surrogate_pvalue(data[d1], data[d2], z, nbins, min_cell, nperm)
        else:
            frac = interaction_fraction(data[d1], data[d2], z, nbins, min_cell)[0]
            null = pv = np.nan
        rows.append((name, frac, null, pv))
    return rows, used


def make_synth():
    rng = np.random.default_rng(0); n = 40000
    Rg = rng.uniform(50, 900, n); VPD = rng.uniform(2, 30, n)
    Rgn = (Rg - 475) / 425; VPDn = (VPD - 16) / 14
    # 分離可能な LUE（掛け算, 交互作用ゼロ期待）
    sep = np.exp(0.8 * Rgn - 0.4 * VPDn) * (1 + rng.normal(0, 0.03, n))
    # 非加法：高 VPD で光利用が落ちる（Rg の効きが VPD で変わる）
    inter = np.exp(0.8 * Rgn - 0.4 * VPDn - 0.7 * Rgn * VPDn) * (1 + rng.normal(0, 0.03, n))
    return Rg, VPD, np.clip(sep, 1e-3, None), np.clip(inter, 1e-3, None)


def main():
    p = argparse.ArgumentParser(description="光合成の掛け算モデル(LUE)の破れ＋呼吸対比")
    p.add_argument("--site")
    p.add_argument("--years", type=int, nargs="+")
    p.add_argument("--month", type=int, nargs="+", default=None)
    p.add_argument("--nbins", type=int, default=6)
    p.add_argument("--min-cell", type=int, default=20)
    p.add_argument("--deyear", action="store_true")
    p.add_argument("--nperm", type=int, default=0)
    a = p.parse_args()

    if not a.site:
        Rg, VPD, sep, inter = make_synth()
        print("=== 旗22 合成検証：光合成 GEP(Rg,VPD) の掛け算(LUE)からのズレ ===")
        for lab, g in [("分離可能 LUE（掛け算）", sep), ("非加法（高VPDで光利用低下）", inter)]:
            frac, null, pv = surrogate_pvalue(Rg, VPD, g, a.nbins, a.min_cell, 300)
            print(f"  {lab:<26} 交互作用={frac:.3f}  ヌル={null:.3f}  p={pv:.3f}")
        print("  → 分離可能≈0・非加法が大きく出れば検出成功。実データは --site で。")
        return

    rows, used = run_site(a.site, a.month, a.deyear, a.nbins, a.min_cell, a.nperm)
    if rows is None:
        print("有効年なし"); return
    dy = "・deyear" if a.deyear else ""
    print(f"=== 旗22 実データ {a.site}（生プール {used}年{dy}, {a.nbins}×{a.nbins}）"
          f"光合成 vs 呼吸のモデル破れ ===\n")
    print(f"  {'系（掛け算モデルの形）':<20} {'交互作用':>8} {'ヌル':>7} {'p':>7}  判定")
    for name, frac, null, pv in rows:
        sig = np.isfinite(pv) and pv < 0.05
        mark = ("★有意に破れる（大）" if (sig and frac >= 0.10) else
                "有意・中" if (sig and frac >= 0.05) else
                "有意だが小" if sig else "非有意（掛け算で概ね書ける）")
        print(f"  {name:<20} {frac:>8.3f} {null:>7.3f} {pv:>7.3f}  {mark}")
    print("\n  読み方：光合成(LUE)が破れる＝放射の光利用が VPD で非加法に変わる（旗13 脱結合のモデル形）。")
    print("         呼吸との対比：どちらの標準モデル形が観測でより破れているか。")
    print("  留保：GEP/GER は NEE 分割の派生量。生プール・年レベルは deyear で除去可。単一サイト。")


if __name__ == "__main__":
    main()
