"""旗18：水田で θ は本当に"情報を失っている"のか？を直接測る。

fig6 の機構説明は「湛水→θ が固定→θ が情報を運べない→呼吸の相乗崩壊」。だが旗17 で水田は
θ×温度 の相乗が最強＝θ は情報を持つ、と出て機構と矛盾した。決着は θ のエントロピーを直接見る：
  ・季節内（各年）の θ の分散・エントロピー（＝湛水で"固定"なら小さいはず）
  ・全年プールの θ の分散・エントロピー（＝年々差があれば季節内より大きい）
両者の差で、「θ は季節内は固定だが年間では動く（＝旗17 はプールの年々差由来）」かを判定できる。

    python research/theta_entropy_step18.py --sites JP-Tak JP-Ta2 CN-HaM JP-Mse
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def _entropy_digitized(x, m):
    x = np.asarray(x, float); lo, hi = x.min(), x.max()
    if hi <= lo:
        return 0.0
    idx = np.clip(np.floor((x - lo) / (hi - lo) * m).astype(int), 0, m - 1)
    c = np.bincount(idx, minlength=m); p = c[c > 0] / len(idx)
    return float(-np.sum(p * np.log(p)))


def main():
    from japanflux_pn.run_robustness import get_site_years
    from japanflux_pn.preprocess import load_corevars_hh

    p = argparse.ArgumentParser(description="θ が水田で固定=無情報かを直接測る")
    p.add_argument("--sites", nargs="+", default=["JP-Tak", "JP-Ta2", "CN-HaM", "JP-Mse"])
    p.add_argument("--var", default="th", help="調べる変数（既定 th=土壌水分）")
    p.add_argument("--obins", type=int, default=6)
    a = p.parse_args()
    paddy = {"JP-Mse", "KR-CRK"}

    print(f"=== {a.var}(土壌水分) は水田で本当に情報を失うか（旗18）===")
    print("  季節内=各年の分散/エントロピー（湛水で固定なら小）、プール=全年まとめ（年々差で増える）\n")
    # CV=std/|mean|（スケール不変の"固定"指標）、H=分布形（min/max規格化）
    print(f"  {'サイト':<8} {'年':>3} {'季節内CV':>9} {'プールCV':>9} {'CV比(プ/内)':>10} "
          f"{'季節内H':>8} {'プールH':>8}")
    for s in a.sites:
        try:
            years, months = get_site_years(s)
        except Exception as e:
            print(f"  {s:<8} get_site_years失敗 {e}"); continue
        within_cv, within_H, pooled = [], [], []
        used = 0
        for y in years:
            try:
                v = load_corevars_hh(s, y, months, None).valid_frame[a.var].to_numpy(float)
            except Exception:
                continue
            mu = abs(float(np.mean(v)))
            within_cv.append(float(np.std(v)) / mu if mu > 1e-9 else np.nan)
            within_H.append(_entropy_digitized(v, a.obins))
            pooled.append(v); used += 1
        if not used:
            print(f"  {s:<8} 有効年なし"); continue
        allv = np.concatenate(pooled)
        wcv, wH = float(np.nanmean(within_cv)), float(np.mean(within_H))
        pmu = abs(float(np.mean(allv)))
        pcv = float(np.std(allv)) / pmu if pmu > 1e-9 else np.nan
        pH = _entropy_digitized(allv, a.obins)
        cvr = pcv / wcv if wcv > 1e-9 else np.inf
        tag = " ←水田" if s in paddy else ""
        print(f"  {s:<8} {used:>3} {wcv:>9.3f} {pcv:>9.3f} {cvr:>10.2f} "
              f"{wH:>8.3f} {pH:>8.3f}{tag}")

    print("\n  読み方（CV=std/mean がスケール不変の『固定』指標。小さいほど θ が動かない）:")
    print("   ・水田の『季節内CV』が自然サイトよりずっと小 → θ は季節内で固定＝『無情報』を支持")
    print("     （fig6 の機構は季節内で正しく、旗17 の相乗はプール年々差 θ 由来＝プーリング交絡）。")
    print("   ・水田の『季節内CV』が自然と同程度 → θ は固定でない＝fig6 の機構『θ無情報』は誤り")
    print("     ＝水田 0/8 の冗長化は別機構（要再解釈）。")
    print("   ・『CV比(プ/内)』が水田で特に大 → 年々差が大＝旗17 の相乗はプール由来の可能性大。")


if __name__ == "__main__":
    main()
