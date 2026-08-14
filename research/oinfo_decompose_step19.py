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
    print(f"  {'サイト':<8} {'年':>3} {'TC(冗長)':>9} {'DTC(相乗)':>9} {'Ω=TC-DTC':>10} {'判定':>8}")
    for s in a.sites:
        try:
            years, months = get_site_years(s)
            raw_all = load_raw_all(get_site(s), cfg)
        except Exception as e:
            print(f"  {s:<8} 読込失敗 {type(e).__name__}: {e}"); continue
        TCs, DTCs, Oms = [], [], []
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
                TCs.append(tc); Oms.append(om); DTCs.append(tc - om)
            except Exception:
                continue
        if not TCs:
            print(f"  {s:<8} 有効年なし"); continue
        tc, dtc, om = np.mean(TCs), np.mean(DTCs), np.mean(Oms)
        verdict = "冗長支配" if om > 0 else "相乗支配"
        tag = " ←水田" if s in paddy else ""
        print(f"  {s:<8} {len(TCs):>3} {tc:>9.4f} {dtc:>9.4f} {om:>10.4f} {verdict:>8}{tag}")

    print("\n  読み方（水田 vs 自然の TC・DTC を比べる）:")
    print("   ・水田で DTC(相乗)が自然より小さい → θが動かず組の相乗が減った＝『相乗が減って冗長支配』")
    print("   ・水田で TC(冗長)が自然より大きい → Rg共通駆動などの冗長が増えた＝『冗長が増えて冗長支配』")
    print("   ・両方 → 合わせ技。＝旗18(θやや固定)が系のどこに効くかを、冗長/相乗の内訳で特定。")
    print("   ※これは4変数系レベルのDTC。旗17のθ×温度ペア相乗(GER目標)とは別の量（両立しうる）。")


if __name__ == "__main__":
    main()
