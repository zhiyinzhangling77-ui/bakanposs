"""旗5：残差の構造＝未観測原因の足跡（問2a「残差を主役に」の最小実験）。

既知の予測子で Y をモデル化し、その残差に構造が残るかを見る。
- 隠れた原因 Z（自己相関あり）が Y を動かしているのに Z を測っていない場合、
  観測できる予測子だけで Y を説明しても、**残差に Z 由来の構造（自己相関）が残る**。
  ＝「まだ説明できていない原因がいる」という足跡。
- 逆に、原因を全部観測できていれば残差は白色（構造なし）になる。

detection（Z を知らないふりで使える）:
  残差 r_t の自己相関、I(r_t; r_{t-1})。構造があれば未観測原因シグナル。
reveal（種明かし）:
  I(r_t; Z_{t-1}) が高い＝残差の構造は隠れた Z の足跡だった、を確認。

    python research/residual_footprint_step5.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from unobserved_common_cause_step4 import mutual_info   # ビン MI を再利用


def _linfit_resid(y, x):
    """y を x で線形回帰し（y≈a+bx）、残差 r=y-(a+bx) を返す。"""
    b, a = np.polyfit(x, y, 1)      # 傾き b, 切片 a
    return y - (a + b * x)


def _autocorr1(r):
    """残差の lag1 自己相関（構造の指標）。"""
    return float(np.corrcoef(r[1:], r[:-1])[0, 1])


def _autocorr_null(r, n=200, seed=1):
    """残差をシャッフルした時の |自己相関| の95%点（＝偶然の上限の目安）。"""
    rng = np.random.default_rng(seed)
    vals = [abs(np.corrcoef((s := rng.permutation(r))[1:], s[:-1])[0, 1]) for _ in range(n)]
    return float(np.quantile(vals, 0.95))


def run_scenario(name, x, y, z, m=8):
    r = _linfit_resid(y, x)                 # 既知の x で y をモデル化した残差
    ac = _autocorr1(r)
    ac_null = _autocorr_null(r)
    mi_rr = mutual_info(r[1:], r[:-1], m)   # I(r_t; r_{t-1})：残差自身の構造
    mi_rz = mutual_info(r[1:], z[:-1], m)   # I(r_t; Z_{t-1})：種明かし（隠れ Z との関係）
    structured = abs(ac) > ac_null
    print(f"【{name}】")
    print(f"  残差の自己相関 r_t~r_(t-1) = {ac:+.3f}  (偶然の上限≈{ac_null:.3f}) → "
          + ("構造あり＝未観測原因の足跡" if structured else "構造なし＝白色（説明しきれている）"))
    print(f"  I(r_t; r_(t-1)) = {mi_rr:.4f} bit   （残差自身にまだ情報がある＝説明不足の量）")
    print(f"  I(r_t; Z_(t-1)) = {mi_rz:.4f} bit   （種明かし：残差は隠れ Z の足跡か）\n")
    return structured


def make_hidden_driver(n=4000, seed=0):
    """隠れた自己相関 Z が X,Y を動かす（Z は未観測）。X は Z の雑な代理。"""
    rng = np.random.default_rng(seed)
    z = np.zeros(n); x = np.zeros(n); y = np.zeros(n)
    ez, ex, ey = (rng.normal(0, 1, n) for _ in range(3))
    for t in range(1, n):
        z[t] = 0.7 * z[t - 1] + ez[t]        # 自己相関のある隠れ駆動
        x[t] = 0.6 * z[t - 1] + ex[t]        # X ← Z（雑な代理）
        y[t] = 1.0 * z[t - 1] + ey[t]        # Y ← Z
    return x, y, z


def make_fully_observed(n=4000, seed=0):
    """原因 X を全部観測できている（X→Y、隠れ駆動なし）。残差は白色になるはず。"""
    rng = np.random.default_rng(seed)
    x = np.zeros(n); y = np.zeros(n); z = np.zeros(n)
    ex, ey, ez = (rng.normal(0, 1, n) for _ in range(3))
    for t in range(1, n):
        x[t] = 0.7 * x[t - 1] + ex[t]
        y[t] = 1.0 * x[t - 1] + ey[t]        # Y ← X（原因を観測済み）
        z[t] = 0.7 * z[t - 1] + ez[t]        # 無関係
    return x, y, z


if __name__ == "__main__":
    print("=== 残差の構造で未観測原因を検出（問2a）===\n")

    # A: 隠れ駆動あり → 残差に構造が残る
    xa, ya, za = make_hidden_driver()
    a_struct = run_scenario("隠れ駆動 Z あり（Z は未観測）", xa, ya, za)

    # B: 原因を全部観測 → 残差は白色（対比）
    xb, yb, zb = make_fully_observed()
    #   Y_t を X_{t-1} で回帰（本物の予測子）→ 残差
    b_struct = run_scenario("原因を全部観測（隠れ駆動なし）",
                            xb[:-1], yb[1:], zb[:-1])

    print("=== 判定 ===")
    ok = a_struct and not b_struct
    print("  " + ("✅ 隠れ駆動ありの時だけ残差に構造＝未観測原因の足跡を検出できた"
                  if ok else "⚠ 期待どおり分離しない（結合/ノイズ/ビン数を確認）"))
    print("  → 実データでは: 既知変数で説明→残差を場所・季節で地図化→構造の偏りを探す→")
    print("     その場所に共通する未観測要因を生態学の知識で仮説化→プロキシで検証（指針 問2a）。")
