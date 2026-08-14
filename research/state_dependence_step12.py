"""旗12：非不変の正体は"状態依存"。ある駆動の効き（係数）が、その年の状態（乾湿）で変わるかを見る。

旗11 で θ→GEP は「強いのに不変でない（CV大）」と出た。その正体は状態依存の可能性：
  水が効くのは乾いた年だけ（湿った年は水が余って効かない）。
そこで各年の θ→GEP の係数 b_年 を、その年の乾燥度（例: 年平均 VPD）と結びつける。
b が状態と相関すれば「効きが状態で変わる＝固定係数で外挿できない＝状態依存」を定量化できる。
＝fig1（乾いた夏に結合が変わる）と同じ物語を、係数の可変性として示す。

    python research/state_dependence_step12.py --site JP-Tak \
        --years 1999 2000 ... 2021 --month 7 8 \
        --target GEP --driver th --state VPD --self-lag 1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from invariance_real_step11 import slope_of


def draw(states, slopes, labels, path, driver, target, state):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    jpf = "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"
    jp = fm.FontProperties(fname=jpf) if Path(jpf).exists() else None
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.scatter(states, slopes, s=55, color="#1f6fb2", zorder=3, edgecolor="white")
    if len(states) >= 2:
        m, b = np.polyfit(states, slopes, 1)
        xs = np.linspace(min(states), max(states), 50)
        ax.plot(xs, m * xs + b, color="#c0392b", lw=2)
    for s, sl, la in zip(states, slopes, labels):
        ax.annotate(la, (s, sl), fontsize=7, color="#666")
    ax.set_xlabel(f"その年の状態: 平均 {state}（乾燥度）", fontproperties=jp)
    ax.set_ylabel(f"{driver}→{target} の効き（年ごとの係数）", fontproperties=jp)
    ax.set_title(f"状態依存: {driver}→{target} の効きは {state} で変わるか",
                 fontproperties=jp)
    fig.savefig(path, bbox_inches="tight", dpi=130); plt.close(fig)


def main() -> None:
    from japanflux_pn.preprocess import load_corevars_hh

    p = argparse.ArgumentParser(description="駆動の効きが状態(乾湿)で変わるか（状態依存）")
    p.add_argument("--site", required=True)
    p.add_argument("--years", type=int, nargs="+", required=True)
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--target", required=True)
    p.add_argument("--driver", required=True, help="効きを見たい駆動（例 th）")
    p.add_argument("--state", default="VPD", help="年の状態を表す変数（例 VPD, th）")
    p.add_argument("--self-lag", type=int, default=1)
    a = p.parse_args()

    states, slopes, labels = [], [], []
    for y in a.years:
        try:
            pre = load_corevars_hh(a.site, y, a.month, None)
            vf = pre.valid_frame
            Y = vf[a.target].to_numpy(dtype=float)
            b = slope_of(Y, vf[a.driver].to_numpy(dtype=float), a.self_lag)
            st = float(vf[a.state].to_numpy(dtype=float).mean())  # 年平均＝その年の状態
            slopes.append(b); states.append(st); labels.append(str(y))
        except Exception as e:
            print(f"  {y}: SKIP {type(e).__name__}: {e}")

    states = np.array(states); slopes = np.array(slopes)
    print(f"\n=== 状態依存: {a.driver}→{a.target} の効き vs 年平均 {a.state} / {a.site} ===")
    print(f"  {'年':>6}  {'年平均'+a.state:>10}  {a.driver+'→'+a.target+'の係数':>16}")
    for la, st, sl in sorted(zip(labels, states, slopes), key=lambda z: z[1]):
        print(f"  {la:>6}  {st:>10.3f}  {sl:>16.3f}")

    if len(states) >= 3:
        r = float(np.corrcoef(states, slopes)[0, 1])
        print(f"\n  相関 r(状態 {a.state} , 係数 {a.driver}→{a.target}) = {r:+.3f}")
        strong = abs(r) >= 0.4
        print("  判定: " + (f"✅ 効きが状態({a.state})で変わる＝状態依存"
                            f"（固定係数では外挿できない）" if strong
                            else "△ 明確な状態依存は見えない（別の状態変数/駆動を試す）"))
        out = Path(__file__).resolve().parent / f"state_dep_{a.site}_{a.driver}_{a.target}.png"
        draw(states, slopes, labels, out, a.driver, a.target, a.state)
        print(f"  [図] {out}")
    print("\n  → 効きが状態で変わる＝その関係は固定関数で将来へ外挿できない（非線形/状態依存）。")
    print("     ＝旗10/11 の『不変な関係だけが転移する』の裏返しを、係数の可変性で定量化。")


if __name__ == "__main__":
    main()
