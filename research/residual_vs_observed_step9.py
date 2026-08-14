"""旗9：残差が"他の観測変数"と相関するかを調べ、未使用の観測駆動と真に未観測を切り分ける。

問2a：既知予測子で目的変数をモデル化した残差が、
  - モデルに入れていない**他の観測変数 V** と相関する → V は「観測してるのに入れ忘れた駆動」
    （＝入れれば説明が改善する。まだ"未観測"ではない）
  - どの観測変数とも相関しない → 残るのは**真に未観測**の候補（領域知識で仮説化→プロキシ検証）

各候補 V について I(残差; V) を測り、大きい順に並べる。＝「次に入れるべき観測変数」の順位表。

    python research/residual_vs_observed_step9.py --site JP-Tak \
        --years 2003 2004 2005 2006 --month 7 8 \
        --target GEP --predictors Rg Ta VPD --self-lag 1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))        # research/
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))    # リポジトリルート（japanflux_pn）
from unobserved_common_cause_step4 import mutual_info


def fit_residual(Yall: np.ndarray, xcols: list[np.ndarray], self_lag: int):
    """Y を xcols(+ Y 自身の過去 self_lag 本) で線形回帰した残差と開始位置 L を返す。"""
    L = self_lag
    Yt = Yall[L:]
    cols = [c[L:] for c in xcols]
    for k in range(1, L + 1):
        cols.append(Yall[L - k:len(Yall) - k])
    A = np.column_stack(cols + [np.ones(len(Yt))])
    coef, *_ = np.linalg.lstsq(A, Yt, rcond=None)
    resid = Yt - A @ coef
    r2 = 1.0 - resid.var() / Yt.var() if Yt.var() > 0 else 0.0
    return resid, float(r2), L


def main() -> None:
    from japanflux_pn.config import RK_VARS
    p = argparse.ArgumentParser(description="残差と他の観測変数の相関で『次に入れるべき変数』を探す")
    p.add_argument("--site", required=True)
    p.add_argument("--years", type=int, nargs="+", required=True)
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--target", required=True)
    p.add_argument("--predictors", nargs="+", required=True)
    p.add_argument("--self-lag", type=int, default=1)
    p.add_argument("--m", type=int, default=8)
    p.add_argument("--include-carbon", action="store_true",
                   help="GEP/GER/NEE も候補に含める（既定は除外：分割の代数的相関を避けるため）")
    a = p.parse_args()

    from japanflux_pn.preprocess import load_corevars_hh

    used = set(a.predictors) | {a.target}
    # GEP/GER/NEE は同じ NEE 分割の親戚で互いに独立でない（代数的相関のアーティファクト）。
    # 目的が炭素変数のときは、残りの炭素変数を候補から除く（--include-carbon で復活）。
    CARBON = {"GEP", "GER", "NEE"}
    drop = set() if a.include_carbon else (CARBON - {a.target})
    candidates = [v for v in RK_VARS if v not in used and v not in drop]
    acc: dict[str, list[float]] = {v: [] for v in candidates}
    r2s = []
    for y in a.years:
        try:
            pre = load_corevars_hh(a.site, y, a.month, None)
            vf = pre.valid_frame
            Yall = vf[a.target].to_numpy(dtype=float)
            xcols = [vf[c].to_numpy(dtype=float) for c in a.predictors]
            resid, r2, L = fit_residual(Yall, xcols, a.self_lag)
            r2s.append(r2)
            for v in candidates:
                vv = vf[v].to_numpy(dtype=float)[L:]           # 残差と同じ時刻に揃える
                acc[v].append(mutual_info(resid, vv, a.m))
        except Exception as e:
            print(f"  {y}: SKIP {type(e).__name__}: {e}")

    print(f"\n=== 残差 vs 他の観測変数: {a.target} ← {'+'.join(a.predictors)}"
          f"(+自己ラグ{a.self_lag})  {a.site} ===")
    print(f"  モデルの平均 R^2 = {np.mean(r2s):.3f}（{len(r2s)}年）")
    print(f"  残差と最も相関する観測変数＝『次に入れるべき駆動の候補』\n")
    ranked = sorted(candidates, key=lambda v: -np.mean(acc[v]))
    print(f"  {'候補変数':>6}  {'I(残差; V) 平均':>14}")
    for v in ranked:
        print(f"  {v:>6}  {np.mean(acc[v]):>14.4f}")
    top = ranked[0]
    print(f"\n  → 最有力: {top}（I={np.mean(acc[top]):.4f}）。")
    print("     大きければ『観測してるのに入れ忘れた駆動』＝入れれば説明改善。")
    print("     全部小さければ、残差は既存の観測変数では説明できない＝真に未観測の候補")
    print("     （領域知識で仮説化→衛星等プロキシで検証）。")


if __name__ == "__main__":
    main()
