"""旗8：残差マップを実データで。目的変数を既知予測子でモデル化し、残差の構造を年ごとに地図化。

問2a：既知変数で説明 → 残差の構造を場所・季節ごとに見る → 偏りのある所を特定 →
未観測要因を仮説化 → プロキシで検証。ここは「実データで残差の構造を年ごとに出す」段。

    # GEP を Rg,Ta,VPD で説明した残差の構造を、年ごとに地図化
    python research/residual_map_real.py --site JP-Tak \
        --years 2003 2004 2005 2006 2007 2008 --month 7 8 \
        --target GEP --predictors Rg Ta VPD

正直な前置き：残差の構造＝「選んだ予測子で説明しきれない分」。未観測変数だけでなく、
非線形・入れ忘れた観測変数でも出る。だから『未観測原因の"候補"の在り処』を示すだけで、
断定はしない（→ 領域知識で候補を絞り、プロキシで検証するのが次段）。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from residual_footprint_step5 import _autocorr1
from unobserved_common_cause_step4 import mutual_info


def residual_structure(y: np.ndarray, xmat: np.ndarray, m: int = 8) -> dict:
    """y を xmat（列＝予測子）で線形回帰し、残差の構造を測る。"""
    A = np.column_stack([xmat, np.ones(len(y))])          # 切片つき計画行列
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ coef
    r2 = 1.0 - resid.var() / y.var() if y.var() > 0 else 0.0
    ac = _autocorr1(resid)                                # 残差の lag1 自己相関
    mi_rr = mutual_info(resid[1:], resid[:-1], m)         # I(r_t; r_{t-1})
    return {"r2": float(r2), "ac": float(ac), "mi_rr": float(mi_rr)}


def draw_map(labels, ac_map, path, target, predictors):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    jp_path = "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"
    jp = fm.FontProperties(fname=jp_path) if Path(jp_path).exists() else None
    r = np.arange(len(ac_map))
    fig, ax = plt.subplots(figsize=(1.2 + 0.5 * len(labels), 4.4))
    ax.bar(r, ac_map, color="#2e8b57")
    ax.set_xticks(r); ax.set_xticklabels(labels, rotation=0)
    ax.set_xlabel("年", fontproperties=jp)
    ax.set_ylabel("残差の自己相関（説明しきれない構造）", fontproperties=jp)
    ax.set_title(f"残差マップ: {target} ← {'+'.join(predictors)}（構造が高い年＝未モデル化の候補）",
                 fontproperties=jp, fontsize=11)
    fig.savefig(path, bbox_inches="tight", dpi=130); plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description="実データの残差マップ（年ごと）")
    p.add_argument("--site", required=True)
    p.add_argument("--years", type=int, nargs="+", required=True)
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--target", required=True, help="説明したい変数（RK 名）")
    p.add_argument("--predictors", nargs="+", required=True, help="既知の予測子（RK 名）")
    p.add_argument("--m", type=int, default=8)
    p.add_argument("--self-lag", type=int, default=0,
                   help="目的変数自身の過去 Y_{t-1..t-L} も予測子に入れる（L=1 推奨）。"
                        "GEP 自身の自己相関を除き、"
                        "『本当に説明しきれない構造』だけを残すため。")
    a = p.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from japanflux_pn.preprocess import load_corevars_hh

    print(f"=== 残差マップ（実データ）: {a.target} ← {'+'.join(a.predictors)} / {a.site} ===")
    print(f"  {'年':>6}  {'n':>6}  {'R^2':>6}  {'残差自己相関':>12}  {'I(r_t;r_(t-1))':>13}")
    labels, ac_map = [], []
    for y in a.years:
        try:
            pre = load_corevars_hh(a.site, y, a.month, None)
            vf = pre.valid_frame
            Yall = vf[a.target].to_numpy(dtype=float)
            cols = [vf[c].to_numpy(dtype=float) for c in a.predictors]
            L = a.self_lag
            if L > 0:
                # Y_t を X_t と Y_{t-1..t-L} で予測（Y 自身の自己相関を除く）
                Yt = Yall[L:]
                Xcols = [c[L:] for c in cols]
                for k in range(1, L + 1):
                    Xcols.append(Yall[L - k:len(Yall) - k])
                s = residual_structure(Yt, np.column_stack(Xcols), m=a.m)
            else:
                s = residual_structure(Yall, np.column_stack(cols), m=a.m)
            labels.append(str(y)); ac_map.append(s["ac"])
            print(f"  {y:>6}  {pre.n_points:>6}  {s['r2']:>6.3f}  "
                  f"{s['ac']:>12.3f}  {s['mi_rr']:>13.4f}")
        except Exception as e:
            print(f"  {y:>6}  SKIP: {type(e).__name__}: {e}")

    if ac_map:
        out = Path(__file__).resolve().parent / f"residual_map_real_{a.site}_{a.target}.png"
        draw_map(labels, ac_map, out, a.target, a.predictors)
        hi = labels[int(np.argmax(ac_map))]
        print(f"\n  残差の構造が最大の年 = {hi}（＝選んだ予測子で説明しきれない構造が強い）")
        print(f"  [図] {out}")
    print("\n  ※構造が高い年/場所は『未観測 or 未モデル化の候補』。断定でなく、")
    print("    領域知識で候補（例: 干ばつ年の土壌水分深さ・フェノロジー等）を絞り、プロキシで検証（問2a）。")


if __name__ == "__main__":
    main()
