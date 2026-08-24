"""旗46（前提監査⑤）：背骨（放射共通駆動の冗長）はギャップフィルの産物か。

前提の穴：層①「エネルギー系 {Rg,Ta,γH,γLE} が最強の冗長」は**既定の gap-fill 込み（_F_MDS）**で
計算している。MDS（Reichstein 2005）は**±7日窓の気象類似日の平均で欠測を埋める**＝埋めた点は
Rg/Ta/VPD の決定的関数であり、**気象-フラックス間の冗長を機械的に押し上げる**。背骨は本研究の
headline の一つなのに、実測のみ（QC=0）で確認したことが一度もない。

**公平な比較のための落とし穴**：QC=0 にすると点数 N が減る。O-info の z はシャッフルヌルとの距離で
N に依存するので、**「z が下がった＝gap-fill の押し上げだった」とは限らない**（ただ点が減っただけかも）。
そこで3本を並べる：
  (a) gap-fill 込み（既定）        … 現状の主張
  (b) 実測のみ QC=0               … 埋めた点を除く
  (c) gap-fill 込みを (b) と同じ N に無作為間引き … N の効果だけを再現した対照
判定は **(b) vs (c)**：同じ N で (b) が (c) より明確に小さければ **gap-fill が冗長を押し上げていた**。
差が無ければ **背骨は実測だけで立つ＝穴⑤は解消**。

    python research/backbone_gapfill_step46.py                 # 合成で検証（機構と検出器）
    python research/backbone_gapfill_step46.py --sites JP-Tak JP-Fhk JP-Tmd
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

BACKBONE = "エネルギー Rg,Ta,γH,γLE"


# ---------- 合成：MDS 風の穴埋めが冗長を押し上げることの実証 -----------------------
def _synth_omega(fill_frac, n=8000, seed=0, bins=11):
    """Rg を共通駆動とする4変数。fill_frac の割合をドライバの決定的関数で"穴埋め"する。"""
    from japanflux_pn import information_theory as it
    rng = np.random.default_rng(seed)
    Rg = rng.normal(0, 1, n)
    Ta = 0.7 * Rg + rng.normal(0, 0.7, n)
    gH = 0.6 * Rg + 0.3 * Ta + rng.normal(0, 0.9, n)      # 独立なフラックス変動を持つ
    gLE = 0.5 * Rg + 0.2 * Ta + rng.normal(0, 0.9, n)
    n_fill = int(n * fill_frac)
    if n_fill:
        m = rng.choice(n, n_fill, replace=False)
        gH[m] = 0.6 * Rg[m] + 0.3 * Ta[m]                 # MDS＝気象類似日の平均＝独自変動ゼロ
        gLE[m] = 0.5 * Rg[m] + 0.2 * Ta[m]
    cols = [it.digitize_series(v, bins) for v in (Rg, Ta, gH, gLE)]
    return it.o_information_indices(cols, bins, correct=True), (Rg, Ta, gH, gLE)


def _omega_subset(vars4, keep_idx, bins=11):
    from japanflux_pn import information_theory as it
    cols = [it.digitize_series(v[keep_idx], bins) for v in vars4]
    return it.o_information_indices(cols, bins, correct=True)


def run_synth():
    print("=== 旗46 合成検証：MDS風の穴埋めは冗長(Ω>0)を押し上げるか、そして検出できるか ===")
    print("  設定：Rg 共通駆動の4変数。穴埋め点はドライバの決定的関数（＝独自変動ゼロ）で置換。\n")
    print(f"  {'穴埋め率':>8} {'Ω(gap-fill込み)':>16}   ※8シード平均（推定ノイズを均す）")
    base = None
    for f in (0.0, 0.15, 0.30, 0.50):
        om = float(np.mean([_synth_omega(f, seed=k)[0] for k in range(8)]))
        if f == 0.0:
            base = om
        print(f"  {f:>7.0%} {om:>16.4f}" + ("  ← 実測のみに相当" if f == 0 else
                                            f"  (+{om - base:.4f})"))
    print("\n  → 穴埋めが増えるほど Ω(冗長)が上がる＝**gap-fill は冗長を機械的に押し上げる**（機構の実証）。\n")

    # 検出器：(b)実測のみ vs (c)同じNに間引いた gap-fill 込み
    rng = np.random.default_rng(1)
    A, B, C = [], [], []
    for k in range(8):                                 # 同じくシード平均
        om_gf_k, vars4 = _synth_omega(0.30, seed=k)
        n = len(vars4[0]); n_meas = int(n * 0.70)
        _, vars_clean = _synth_omega(0.0, seed=k)
        A.append(om_gf_k)
        B.append(_omega_subset(vars_clean, np.arange(n_meas)))            # (b) 実測のみ
        C.append(_omega_subset(vars4, rng.choice(n, n_meas, replace=False)))  # (c) 同N対照
    om_gf, om_b, om_c = map(lambda x: float(np.mean(x)), (A, B, C))
    print(f"  (a) gap-fill込み 全N      Ω={om_gf:.4f}")
    print(f"  (b) 実測のみ    N={n_meas}  Ω={om_b:.4f}")
    print(f"  (c) 同N間引き対照 N={n_meas}  Ω={om_c:.4f}")
    print(f"  → (b)−(c) = {om_b - om_c:+.4f}：負に大きいほど「冗長は gap-fill 由来」。")
    print("     N を揃えた (c) と比べるので、点数が減ったことによる見かけの低下を除いてある。")


# ---------- 実データ ---------------------------------------------------------------
def run_real(sites, obins):
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.run_robustness import get_site_years
    from japanflux_pn.preprocess import (load_raw_all, slice_and_anomaly,
                                         slice_span_and_anomaly, PreprocessResult)
    from japanflux_pn.oinfo_analysis import o_information_subsystems

    cfg_gf = AnalysisConfig()                       # gap-fill 込み（既定）
    cfg_m0 = AnalysisConfig(qc_max=0)               # 実測のみ
    rng = np.random.default_rng(0)

    def _slice(raw, y, months, cfg):
        if len(months) == 1:
            return slice_and_anomaly(raw, y, months[0], cfg)
        return slice_span_and_anomaly(raw, y, months, cfg)

    print("=== 旗46 実データ：背骨（エネルギー系の冗長）はギャップフィルの産物か ===")
    print("  Ω>0=冗長支配。(b)実測のみ と (c)同N間引きの gap-fill込み を比べるのが判定。\n")
    rows, skipped = [], []
    for site in sites:
        years, months = get_site_years(site)
        spec = get_site(site)
        try:
            raw_gf = load_raw_all(spec, cfg_gf)
            raw_m0 = load_raw_all(spec, cfg_m0)
        except Exception as e:
            print(f"  {site}: 読み込み失敗 {type(e).__name__}"); continue
        for y in years:
            try:
                an_gf, va_gf = _slice(raw_gf, y, months, cfg_gf)
                an_m0, va_m0 = _slice(raw_m0, y, months, cfg_m0)
            except Exception:
                continue
            n_gf, n_m0 = int(va_gf.sum()), int(va_m0.sum())
            if n_m0 < 500 or n_gf < 500:
                skipped.append((site, y, n_gf, n_m0))
                continue
            def _tbl(anom, valid, cfg):
                pre = PreprocessResult(anomaly=anom, valid=valid, site=site, year=y,
                                       month=months[0], config=cfg, months=months)
                return o_information_subsystems(pre, obins).set_index("subsystem")
            t_gf, t_m0 = _tbl(an_gf, va_gf, cfg_gf), _tbl(an_m0, va_m0, cfg_m0)
            # (c) 同N対照：gap-fill込みの有効時刻を無作為に n_m0 点へ間引く
            idx_true = np.flatnonzero(va_gf.to_numpy())
            keep = rng.choice(idx_true, n_m0, replace=False)
            va_c = va_gf.copy(); va_c[:] = False; va_c.iloc[np.sort(keep)] = True
            t_c = _tbl(an_gf, va_c, cfg_gf)
            for sub in t_gf.index:
                rows.append({"site": site, "year": y, "sub": sub,
                             "n_gf": n_gf, "n_m0": n_m0,
                             "gf": t_gf.loc[sub, "Omega"], "m0": t_m0.loc[sub, "Omega"],
                             "ctl": t_c.loc[sub, "Omega"],
                             "z_gf": t_gf.loc[sub, "z"], "z_m0": t_m0.loc[sub, "z"]})
            print(f"  {site} {y}: 実測率 {n_m0/n_gf:.0%} (N {n_gf}→{n_m0})", flush=True)

    if skipped:
        print(f"\n  除外 {len(skipped)} site-year（実測点<500＝O-info推定に不足）：")
        for st, y, a_, b_ in skipped[:12]:
            print(f"    {st} {y}: N {a_}→{b_}（実測率 {b_/a_:.0%}）")
        if len(skipped) > 12:
            print(f"    …他 {len(skipped)-12}")
    if not rows:
        print("  有効な site-year なし"); return
    import pandas as pd
    df = pd.DataFrame(rows)
    print(f"\n  === サブシステム別（{df['site'].nunique()}サイト {len(df.groupby(['site','year']))} site-year 平均）===")
    print(f"  {'subsystem':<24} {'Ω(a)gap込':>10} {'Ω(b)実測':>10} {'Ω(c)同N対照':>12} "
          f"{'(b)−(c)':>9}  判定")
    for sub, g in df.groupby("sub"):
        a, b, c = g["gf"].mean(), g["m0"].mean(), g["ctl"].mean()
        d = b - c
        n_lower = int((g["m0"] < g["ctl"]).sum())            # (b)<(c) の site-year 数
        n_flip = int(((g["gf"] > 0) & (g["m0"] < 0)).sum())  # 冗長→相乗の符号反転
        scale = max(abs(c), 0.05)                            # 0近傍で割らない（相対値の暴走を防ぐ）
        if n_flip >= 0.6 * len(g):
            v = f"▲符号反転 gap込=冗長→実測=相乗 ({n_flip}/{len(g)})"
        elif abs(d) / scale < 0.15:
            v = f"★実測だけで立つ ({n_lower}/{len(g)}が低下)"
        elif d < 0:
            v = f"▲gap-fillが押し上げ ({n_lower}/{len(g)}が低下)"
        else:
            v = f"○実測の方が強い ({len(g)-n_lower}/{len(g)})"
        print(f"  {sub:<24} {a:>10.4f} {b:>10.4f} {c:>12.4f} {d:>+9.4f}  {v}")
    print("\n  注：相対比は |Ω(c)| が 0 近傍だと発散するので下限 0.05 で丸めてある。")
    print("      符号反転（gap込みでは冗長 Ω>0 なのに実測では相乗 Ω<0）が多数なら、その系の"
          "『冗長支配』は gap-fill が作っていた可能性が高い。")
    bb = df[df["sub"] == BACKBONE]
    if not bb.empty:
        print(f"\n  === 背骨 {BACKBONE} ===")
        print(f"  実測率 中央値 {(bb['n_m0']/bb['n_gf']).median():.0%}／"
              f"Ω 実測 {bb['m0'].mean():.4f} vs 同N対照 {bb['ctl'].mean():.4f}／"
              f"z 実測 {bb['z_m0'].mean():+.1f}")
        print("  読み：(b)≈(c) なら背骨は実測だけで立つ＝穴⑤解消。(b)≪(c) なら冗長の相当部分は"
              "MDS が気象から作った写り込み＝主張の格下げが要る。")


def main():
    p = argparse.ArgumentParser(description="背骨の冗長がgap-fill由来かを判定")
    p.add_argument("--sites", nargs="+")
    p.add_argument("--obins", type=int, default=8)
    a = p.parse_args()
    if a.sites:
        run_real(a.sites, a.obins)
    else:
        run_synth()


if __name__ == "__main__":
    main()
