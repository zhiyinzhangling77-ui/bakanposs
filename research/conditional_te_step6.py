"""旗6：条件付き Transfer Entropy。共通原因 Z を条件づけると"見かけの流れ"が消える。

実データ（Rg↔GEP）は両方向の TE がほぼ同じ大きさで両方有意だった＝共通駆動（交絡）の署名。
そこで「Z を条件に入れた TE」を作り、共通原因を差し引くと対称な流れが崩れるかを見る。
  条件付き TE:  TE(X→Y|Z) = I(Y_t ; X_{t-lag} | Y_{t-1}, Z_{t-condlag})
これは指針 問2a「I(X;Y|W) で交絡を除く」の TE 版。

まず合成（Z→X, Z→Y の共通原因）で「条件づけ前は両方向とも高い→Z を条件づけると両方消える」
を確認する（実データで見た対称パターンの再現と解決）。

    python research/conditional_te_step6.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from transfer_entropy_step1 import digitize, transfer_entropy
from unobserved_common_cause_step4 import make_common_cause


def conditional_te(x, y, z, m=6, lag=1, condlag=1) -> float:
    """TE(X→Y|Z) = I(Y_t ; X_{t-lag} | Y_{t-1}, Z_{t-condlag})（ビン推定・4次元）。
    ★条件づけで次元が増えるので、ビン数 m は小さめ（6 前後）にしないとセルが疎になる。
    """
    xi, yi, zi = digitize(x, m), digitize(y, m), digitize(z, m)
    n = len(yi)
    t0 = max(1, lag, condlag)
    yf = yi[t0:n]              # Y_t
    yp = yi[t0 - 1:n - 1]      # Y_{t-1}（条件その1）
    xs = xi[t0 - lag:n - lag]  # X_{t-lag}
    zc = zi[t0 - condlag:n - condlag]  # Z_{t-condlag}（条件その2＝共通原因）

    # 軸: (yf, xs, yp, zc)。条件 C=(yp,zc) をまとめて扱う。
    j = np.zeros((m, m, m, m)); np.add.at(j, (yf, xs, yp, zc), 1.0); j /= j.sum()
    p_c = j.sum(axis=(0, 1))          # p(yp,zc)
    p_ac = j.sum(axis=1)              # p(yf,yp,zc)
    p_bc = j.sum(axis=0)              # p(xs,yp,zc)
    s = 0.0
    for a, b, c, d in np.argwhere(j > 0):   # a=yf,b=xs,c=yp,d=zc
        num = j[a, b, c, d] * p_c[c, d]
        den = p_ac[a, c, d] * p_bc[b, c, d]
        if num > 0 and den > 0:
            s += j[a, b, c, d] * np.log2(num / den)
    return float(s)


if __name__ == "__main__":
    M = 6
    print("=== 条件付き TE：共通原因 Z を差し引く（合成 Z→X, Z→Y）===")
    print("  ★条件づけ後の TE(X→Y|Z) が 0 に近づけば『対称な流れは共通原因の影』と分かる。")
    print("  ★ただし4次元ビンは疎になりやすく、データ量 N で結果が激変することを見せる。\n")

    print(f"  {'N':>8}  {'TE(X→Y)':>9}  {'TE(X→Y|Z)':>10}  {'TE(Y→X)':>9}  {'TE(Y→X|Z)':>10}  判定")
    for N in (6000, 40000, 120000):
        x, y, z = make_common_cause(n=N, seed=0)
        te_xy = transfer_entropy(x, y, m=M, lag=1)
        te_yx = transfer_entropy(y, x, m=M, lag=1)
        cte_xy = conditional_te(x, y, z, m=M, lag=1, condlag=1)
        cte_yx = conditional_te(y, x, z, m=M, lag=1, condlag=1)
        ok = cte_xy < 0.5 * te_xy and cte_yx < 0.5 * te_yx
        mk = "✅崩れる" if ok else "⚠バイアス支配"
        print(f"  {N:>8}  {te_xy:9.4f}  {cte_xy:10.4f}  {te_yx:9.4f}  {cte_yx:10.4f}  {mk}")

    print("\n=== 教訓（実データへの含意）===")
    print("  ・概念は正しい: 十分な N では TE(X→Y|Z) が 0 に近づき、対称な流れが共通原因の影だと分かる。")
    print("  ・だが N が小さい(≈6000)と4次元ビンが疎になり、推定バイアスが信号を飲む（＝見かけ減らない）。")
    print("  ・実フラックスの1サイト1夏は n≈3000。だから naïve なビン条件付けは危険で、")
    print("    本体は Miller-Madow 補正・m 調整・CMIknn(KNN) を使う（condition_driver / PCMCI）。")
    print("  → Rg↔GEP の対称 TE を『解く』には、少数標本に強い推定＝本体の条件付け/PCMCI が要る。")
