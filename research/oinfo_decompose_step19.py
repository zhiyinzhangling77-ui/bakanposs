"""旗19：呼吸系の O-information を冗長成分(TC)と相乗成分(DTC)に分解して、
水田の"冗長支配"が「冗長が増えた(TC↑)」のか「相乗が減った(DTC↓)」のかを見る。

O-information Ω = TC − DTC。
  TC（全相関）＝変数間で共有された依存の総量＝**冗長寄り**。
  DTC（双対全相関）＝組にしか宿らない依存＝**相乗寄り**。
Ω>0 冗長支配。だが「なぜ冗長支配か」は Ω だけでは分からない：TC が大きいのか、DTC が小さいのか。
旗18 で水田の θ はやや固定と分かった。θ が動かないと、θ が絡む相乗(DTC)が減るのか、
Rg 共通駆動の冗長(TC)が相対的に効くのか——を、TC/DTC を分けて測って切り分ける。

呼吸系 {Rg,Ta,θ,GER}。各サイト健全年ごとに TC・Ω・DTC（Miller-Madow 補正）を出し、年平均で比較。

    python research/oinfo_decompose_step19.py --sites JP-Tak JP-Ta2 CN-HaM JP-Mse
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

SUBSYS = ["Rg", "Ta", "th", "GER"]


def main():
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import (load_raw_all, slice_and_anomaly,
                                         slice_span_and_anomaly)
    from japanflux_pn.run_robustness import get_site_years
    from japanflux_pn import information_theory as it

    p = argparse.ArgumentParser(description="呼吸系O-infoをTC(冗長)とDTC(相乗)に分解")
    p.add_argument("--sites", nargs="+", default=["JP-Tak", "JP-Ta2", "CN-HaM", "JP-Mse"])
    p.add_argument("--obins", type=int, default=6)
    a = p.parse_args()
    cfg = AnalysisConfig()
    m = a.obins
    paddy = {"JP-Mse", "KR-CRK"}

    print(f"=== 呼吸系 {{{','.join(SUBSYS)}}} の O-info 分解（旗19）===")
    print("  Ω=TC−DTC。TC=冗長寄り(共有依存)、DTC=相乗寄り(組にしか宿る)。Ω>0 冗長支配。")
    print("  水田の冗長支配は『TCが大』か『DTCが小』か？を切り分ける\n")
    print("  ★判定は生Ωの絶対符号でなく z（シャッフルヌル比）で（4変数Ωは負バイアスのため）。")
    print("    TC/DTC の絶対値も同じバイアスを含むが、サイト間の相対差は読める。\n")
    print(f"  {'サイト':<8} {'年':>3} {'TC(冗長)':>9} {'DTC(相乗)':>9} {'Ω':>9} {'z(Ω)':>7} {'判定':>8}")
    for s in a.sites:
        try:
            years, months = get_site_years(s)
            raw_all = load_raw_all(get_site(s), cfg)
        except Exception as e:
            print(f"  {s:<8} 読込失敗 {type(e).__name__}: {e}"); continue
        TCs, DTCs, Oms, Zs = [], [], [], []
        for y in years:
            try:
                if len(months) == 1:
                    anom, valid = slice_and_anomaly(raw_all, y, months[0], cfg)
                else:
                    anom, valid = slice_span_and_anomaly(raw_all, y, months, cfg)
                if int(valid.sum()) < 500:
                    continue
                vf = anom.loc[valid]
                cols = [it.digitize_series(vf[v].to_numpy(float), m) for v in SUBSYS]
                tc = it.total_correlation_indices(cols, m, correct=True)
                om = it.o_information_indices(cols, m, correct=True)
                rng = np.random.default_rng(cfg.seed)
                st = it.surrogate_o_information_stats(
                    cols, m, cfg.n_surrogates, cfg.sig_c, rng, correct=True)
                z = (om - st["mu"]) / st["sigma"] if st["sigma"] > 0 else np.nan
                TCs.append(tc); Oms.append(om); DTCs.append(tc - om); Zs.append(z)
            except Exception:
                continue
        if not TCs:
            print(f"  {s:<8} 有効年なし"); continue
        tc, dtc, om, mz = np.mean(TCs), np.mean(DTCs), np.mean(Oms), np.nanmean(Zs)
        verdict = ("冗長支配" if mz >= 2.36 else "相乗支配" if mz <= -2.36 else "曖昧")
        tag = " ←水田" if s in paddy else ""
        print(f"  {s:<8} {len(TCs):>3} {tc:>9.4f} {dtc:>9.4f} {om:>9.4f} {mz:>7.1f} "
              f"{verdict:>8}{tag}")

    print("\n  読み方:")
    print("   ・判定(冗長/相乗)は z の符号で。生Ωが全サイト負でも、それは疎性バイアス（絶対符号は無意味）。")
    print("   ・(1)の答え＝TC/DTC のサイト間『相対差』：水田で DTC(相乗成分)が自然より小さければ")
    print("     『相乗が減って冗長寄りに』、TC(冗長成分)が大きければ『冗長が増えて』。TC はどこも同程度なら前者。")
    print("   ・＝旗18(θやや固定)が系の相乗成分(DTC)を削るのか冗長成分(TC)を増すのかを特定。")
    print("   ※4変数系レベルのDTC。旗17のθ×温度ペア相乗(GER目標)とは別の量（両立しうる）。")


if __name__ == "__main__":
    main()
