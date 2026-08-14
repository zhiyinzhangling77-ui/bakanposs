"""合成 → 実データへの橋渡し：実フラックスの2変数で TE を1回動かす。

research/ の TE（自前実装）を、この研究本体の前処理（load_corevars_hh：5日アノマリ＋
listwise）で作った実データの2変数に当てる。向き・有意性・遅れ（lag）を出す。

    python research/te_real_data.py --site JP-Tak --year 2003 --month 7 8 --x Rg --y GEP

※正直な前置き：これは**2変数だけの見かけの情報の流れ**。共通原因（他の変数）は差し引いて
  いないので、因果の断定はできない（それをやるのが本体の condition_driver / PCMCI）。
  ここでの目的は「合成で確かめた道具を、実データで1回動かす」こと。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# research/ 内の TE 関数
sys.path.insert(0, str(Path(__file__).resolve().parent))
from transfer_entropy_step1 import transfer_entropy
from transfer_entropy_step2 import te_with_significance


def analyze_pair(x: np.ndarray, y: np.ndarray, xname: str, yname: str,
                 m: int = 8, lagmax: int = 8, n_surr: int = 200) -> None:
    print(f"\n=== TE 解析: {xname} と {yname}  (n={len(x)}, m={m}) ===")
    # 1) lag スキャンで両方向の TE を見る
    print(f"  {'lag':>4}  {'TE(%s→%s)'%(xname,yname):>14}  {'TE(%s→%s)'%(yname,xname):>14}")
    best_lag, best = 1, -1.0
    for lag in range(1, lagmax + 1):
        a = transfer_entropy(x, y, m=m, lag=lag)
        b = transfer_entropy(y, x, m=m, lag=lag)
        if a > best:
            best, best_lag = a, lag
        print(f"  {lag:>4}  {a:14.4f}  {b:14.4f}")
    print(f"  → TE({xname}→{yname}) 最大の lag = {best_lag}")

    # 2) その lag で両方向の有意性（サロゲート）
    fwd = te_with_significance(x, y, m=m, lag=best_lag, n_surr=n_surr)
    rev = te_with_significance(y, x, m=m, lag=best_lag, n_surr=n_surr)
    print(f"\n  [lag={best_lag}] 有意性（μ+2.36σ, α≈0.01）")
    for nm, r in [(f"{xname}→{yname}", fwd), (f"{yname}→{xname}", rev)]:
        mk = "✅有意" if r["significant"] else "×非有意"
        print(f"    TE {nm}: {r['te']:.4f} bit  z={r['z']:.1f}σ  p={r['p']:.3f}  → {mk}")
    print("\n  ※2変数のみ。共通原因は未除去なので因果の断定は不可（→ 本体の条件付け/PCMCI）。")


def main() -> None:
    p = argparse.ArgumentParser(description="実フラックス2変数の TE を1回動かす")
    p.add_argument("--site", required=True)
    p.add_argument("--year", type=int, required=True)
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--x", required=True, help="送り手変数（RK 名: Rg, Ta, VPD, Ts, P, th, gH, gLE, GER, NEE, GEP）")
    p.add_argument("--y", required=True, help="受け手変数（RK 名）")
    p.add_argument("--m", type=int, default=8)
    p.add_argument("--lagmax", type=int, default=8)
    p.add_argument("--n-surr", type=int, default=200)
    a = p.parse_args()

    # 本体の前処理を再利用（5日アノマリ＋11変数 listwise）
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from japanflux_pn.preprocess import load_corevars_hh
    pre = load_corevars_hh(a.site, a.year, a.month, None)
    vf = pre.valid_frame
    print(f"[preprocess] {a.site} {a.year}-{pre.month_label}: n_points={pre.n_points}")
    x = vf[a.x].to_numpy(dtype=float)
    y = vf[a.y].to_numpy(dtype=float)
    analyze_pair(x, y, a.x, a.y, m=a.m, lagmax=a.lagmax, n_surr=a.n_surr)


if __name__ == "__main__":
    main()
