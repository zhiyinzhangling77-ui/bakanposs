"""次の旗：Transfer Entropy にサロゲート検定を足し、「その TE は偶然か本物か」を判定する。

step1 で TE の数字は出た。だが TE はビン推定で少し正に偏る（逆向きでも 0 にならない）。
そこで「送り手を時間シャッフルして、向きの関係を壊した TE の分布（帰無分布）」と観測 TE を
比べ、観測がその上ずっと外側にあれば "偶然でない（有意）" と言える。
＝指針 問2a の「条件付き独立性で足跡を検出する」の 2 変数版。

    python research/transfer_entropy_step2.py

ゴール（合成データ）:
  本物の向き TE(X→Y) は有意、無いはずの TE(Y→X) は非有意、になることを確認する。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# step1 の関数を再利用（自前実装・numpy のみ）
sys.path.insert(0, str(Path(__file__).resolve().parent))
from transfer_entropy_step1 import transfer_entropy, make_synthetic


# ---------------------------------------------------------------------------
# サロゲート検定：送り手 (source) を時間シャッフルして「向きの結合」を壊し、
# その TE を n_surr 回計算して帰無分布を作る。観測 TE がその分布の
# どれだけ外側かで有意性を判定する。
#   ★n_surr: 多いほど帰無分布が滑らか（100〜500 が定番）。
#   ★閾値 c=2.36: 片側 α≈0.01（μ+2.36σ）。R&K の設定に合わせた。
# ---------------------------------------------------------------------------
def te_with_significance(src: np.ndarray, dst: np.ndarray, m: int = 8, lag: int = 3,
                         n_surr: int = 200, c: float = 2.36, seed: int = 1) -> dict:
    te_obs = transfer_entropy(src, dst, m=m, lag=lag)     # 観測 TE(src→dst)
    rng = np.random.default_rng(seed)
    null = np.empty(n_surr)
    for i in range(n_surr):
        src_shuf = rng.permutation(src)                   # 送り手だけ時間を壊す
        null[i] = transfer_entropy(src_shuf, dst, m=m, lag=lag)
    mu, sd = float(null.mean()), float(null.std())
    thr = mu + c * sd                                     # 有意閾値（μ+2.36σ）
    z = (te_obs - mu) / sd if sd > 0 else np.inf          # 何σ外側か
    p = float((null >= te_obs).mean())                    # 経験 p 値
    return {"te": te_obs, "null_mu": mu, "null_sd": sd,
            "threshold": thr, "z": z, "p": p, "significant": te_obs > thr}


if __name__ == "__main__":
    M, DELAY = 8, 3
    LAG = DELAY
    x, y = make_synthetic(n=4000, delay=DELAY, coupling=0.7, seed=0)

    xy = te_with_significance(x, y, m=M, lag=LAG, n_surr=200)   # X→Y（本物）
    yx = te_with_significance(y, x, m=M, lag=LAG, n_surr=200)   # Y→X（無いはず）

    print("=== Transfer Entropy ＋ サロゲート有意性（合成データ）===")
    print(f"  設定: m={M}, lag={LAG}, 真の遅れ={DELAY}, サロゲート=200, 閾値=μ+2.36σ(α≈0.01)\n")
    for name, r in [("X→Y（本物の向き）", xy), ("Y→X（無いはずの向き）", yx)]:
        mark = "✅有意" if r["significant"] else "×非有意"
        print(f"  TE {name}: {r['te']:.4f} bit   "
              f"帰無 μ={r['null_mu']:.4f} σ={r['null_sd']:.4f}  "
              f"閾値={r['threshold']:.4f}  z={r['z']:.1f}σ  p={r['p']:.3f}  → {mark}")

    ok = xy["significant"] and not yx["significant"]
    print("\n  判定: " + ("✅ 本物の向きだけが有意（偶然と区別できた）"
                          if ok else "⚠ 期待どおりに分離しない（m/lag/サロゲート数を見直す）"))
    # 読み方:
    #  ・帰無 μ ≒ ビン推定の上ずれ（バイアスの床）。観測がこの床より 2.36σ 以上高ければ有意。
    #  ・z は「偶然の分布から何σ外れているか」。p は「偶然でこれ以上の TE が出る割合」。
    #  ・本物の向きは z が大きく p≒0、無い向きは z が小さく p が大きい、が期待。
