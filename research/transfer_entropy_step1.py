"""最初の旗：Transfer Entropy を 2 本の時系列に対して 1 回だけ計算する最小構成。

目的（これだけ）:
  - 2 本の時系列 X, Y から、両方向の Transfer Entropy TE(X→Y) と TE(Y→X) を計算し、
    数字として出す。
  - 既知の因果（X が d ステップ遅れて Y を動かす）を持つ合成データで、
    TE(X→Y) > TE(Y→X) になることを自分の目で確認する。

ライブラリの選択:
  - 自前実装（numpy のみ）。追加インストール不要で依存が最も軽く、ビン分割や lag の
    扱いがコード上で丸見えになるため（JIDT は Java、pyinform は C 拡張の追加が必要）。

TE とは（1 行）:
  TE(X→Y) = 「Y の次の値を予想するとき、Y 自身の過去に加えて X の過去を知ると、
             どれだけ不確かさ（ビット）が減るか」。＝X から Y への情報の流れ。
  形式的には条件付き相互情報量 I(Y_t ; X_{t-d} | Y_{t-1})。
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# 1. ビン分割：連続値を m 個の等幅ビン（0..m-1 の整数ラベル）に離散化する。
#    ★パラメータ m（ビン数）: 結果を左右する。大きいほど関係を細かく捉えるが、
#      1 ビンあたりのデータが減って推定が不安定（＋TE が上振れ）になる。
#      少ないデータなら m は小さめ（8 前後）が無難。
# ---------------------------------------------------------------------------
def digitize(v: np.ndarray, m: int) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    edges = np.linspace(v.min(), v.max(), m + 1)   # 等幅のビン境界
    idx = np.clip(np.digitize(v, edges[1:-1]), 0, m - 1)  # 0..m-1 のラベル
    return idx.astype(int)


# ---------------------------------------------------------------------------
# 2. Transfer Entropy（ビン推定・履歴長 1）
#    TE(X→Y; lag=d) = I(Y_t ; X_{t-d} | Y_{t-1})
#      = Σ p(yf, yp, xs) · log2[ p(yf,yp,xs)·p(yp) / ( p(yf,yp)·p(yp,xs) ) ]
#    yf = Y_t（予想したい未来）, yp = Y_{t-1}（Y 自身の過去＝これで条件付け）,
#    xs = X_{t-d}（送り手 X の過去）。
#    ★パラメータ lag=d（遅れ）: 「X が何ステップ前に Y を動かすか」を指定する。
#      合成データの真の遅れに合わせると向きが最もはっきり出る。ズレると弱く出る。
#    単位は bit（log の底が 2）。
# ---------------------------------------------------------------------------
def transfer_entropy(x: np.ndarray, y: np.ndarray, m: int = 8, lag: int = 1) -> float:
    xi = digitize(x, m)
    yi = digitize(y, m)
    n = len(yi)
    t0 = max(1, lag)                       # y_{t-1} と x_{t-lag} が両方存在する最初の t
    yf = yi[t0:n]                          # Y_t
    yp = yi[t0 - 1:n - 1]                  # Y_{t-1}
    xs = xi[t0 - lag:n - lag]              # X_{t-lag}

    # 3 次元の同時ヒストグラム（頻度）→ 確率
    joint = np.zeros((m, m, m), dtype=float)   # 軸: (yf, yp, xs)
    np.add.at(joint, (yf, yp, xs), 1.0)
    joint /= joint.sum()

    p_yp = joint.sum(axis=(0, 2))              # p(yp)
    p_fy = joint.sum(axis=2)                   # p(yf, yp)
    p_yx = joint.sum(axis=0)                   # p(yp, xs)

    te = 0.0
    nz = np.argwhere(joint > 0)                # 0 のセルは log で無視
    for a, b, c in nz:                         # a=yf, b=yp, c=xs
        num = joint[a, b, c] * p_yp[b]
        den = p_fy[a, b] * p_yx[b, c]
        if num > 0 and den > 0:
            te += joint[a, b, c] * np.log2(num / den)
    return float(te)


# ---------------------------------------------------------------------------
# 3. 合成データ：X が d ステップ遅れて Y を動かす（Y→X の経路は無い）。
#    これで「向き（X→Y）」が正しく出るかを検証する。
# ---------------------------------------------------------------------------
def make_synthetic(n: int = 4000, delay: int = 3, coupling: float = 0.7,
                   seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    x = np.zeros(n)
    y = np.zeros(n)
    ex = rng.normal(0, 1, n)   # X のノイズ
    ey = rng.normal(0, 1, n)   # Y のノイズ
    for t in range(1, n):
        x[t] = 0.5 * x[t - 1] + ex[t]                     # X は自分の過去だけで動く
        src = x[t - delay] if t - delay >= 0 else 0.0      # ★X を delay ステップ遅らせて
        y[t] = coupling * src + 0.3 * y[t - 1] + ey[t]     #   Y に注入（X→Y の一方向）
    return x, y


# ---------------------------------------------------------------------------
# 4. 実行：両方向の TE を計算して表示。
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    M = 8          # ビン数
    DELAY = 3      # 真の遅れ（X が 3 ステップ先の Y を動かす）
    LAG = DELAY    # TE を測る遅れ（真の遅れに合わせる）

    x, y = make_synthetic(n=4000, delay=DELAY, coupling=0.7, seed=0)

    te_xy = transfer_entropy(x, y, m=M, lag=LAG)   # X→Y（本物の向き）
    te_yx = transfer_entropy(y, x, m=M, lag=LAG)   # Y→X（無いはずの向き）

    print("=== 最小構成 Transfer Entropy（合成データ・1 回計算）===")
    print(f"  設定: ビン数 m={M}, 遅れ lag={LAG}, 真の遅れ delay={DELAY}, N={len(x)}")
    print(f"  TE(X→Y) = {te_xy:.4f} bit   ← 本物の因果の向き")
    print(f"  TE(Y→X) = {te_yx:.4f} bit   ← 無いはずの向き")
    print(f"  差 TE(X→Y) − TE(Y→X) = {te_xy - te_yx:+.4f} bit")
    ok = te_xy > te_yx
    print("  判定: " + ("✅ TE(X→Y) > TE(Y→X) ＝ 向きを正しく捉えた"
                        if ok else "⚠ 向きが出ていない（m や lag を見直す）"))

    # 読み方（コメント）:
    #  ・単位は bit。TE(X→Y)=v は「X の過去を知ると Y の次の値の不確かさが v bit 減る」。
    #  ・0 に近い＝その向きの情報の流れは（この遅れでは）ほぼ無い。
    #  ・大きいほど強い一方向の流れ。ただし絶対値は m・結合・ノイズ・データ量で変わるので、
    #    "どちらが大きいか（向き）" を主に読む。厳密な有意性はサロゲート検定が要る（次段）。
