"""旗4：未観測共通原因（問2a の核心）。

「Z が X と Y を両方動かす」と、X と Y には直接の因果が無いのに関連が出る（見かけのつながり）。
- Z を測っていない → I(X;Y) は高く、X↔Y が本物に見えてしまう（＝未観測交絡の罠）。
- Z を測って条件づける → I(X;Y|Z) が消える → 「Z が共通原因だった」と分かる。
これは指針 問2a の道具「条件付き独立性 I(X;Y|W) で交絡因子を検出」の最小実験。

対比として、本当に X→Y があるケースでは、無関係な Z で条件づけても I(X;Y|Z) は残る、も見る。

    python research/unobserved_common_cause_step4.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from transfer_entropy_step1 import digitize   # ビン分割を再利用


# --- 相互情報量 I(A;B) [bit]（ビン推定）--------------------------------------
def mutual_info(a, b, m=8):
    ai, bi = digitize(a, m), digitize(b, m)
    j = np.zeros((m, m)); np.add.at(j, (ai, bi), 1.0); j /= j.sum()
    pa, pb = j.sum(1), j.sum(0)
    s = 0.0
    for i in range(m):
        for k in range(m):
            if j[i, k] > 0:
                s += j[i, k] * np.log2(j[i, k] / (pa[i] * pb[k]))
    return float(s)


# --- 条件付き相互情報量 I(A;B|C) [bit] = Σ p(abc) log2[ p(abc)p(c)/(p(ac)p(bc)) ] --
def cond_mutual_info(a, b, c, m=8):
    ai, bi, ci = digitize(a, m), digitize(b, m), digitize(c, m)
    j = np.zeros((m, m, m)); np.add.at(j, (ai, bi, ci), 1.0); j /= j.sum()
    pc = j.sum((0, 1)); pac = j.sum(1); pbc = j.sum(0)
    s = 0.0
    for i, k, l in np.argwhere(j > 0):
        num = j[i, k, l] * pc[l]
        den = pac[i, l] * pbc[k, l]
        if num > 0 and den > 0:
            s += j[i, k, l] * np.log2(num / den)
    return float(s)


def make_common_cause(n=4000, seed=0):
    """Z が X, Y を両方動かす（X–Y に直接の因果は無い）。"""
    rng = np.random.default_rng(seed)
    z = np.zeros(n); x = np.zeros(n); y = np.zeros(n)
    ez, ex, ey = (rng.normal(0, 1, n) for _ in range(3))
    for t in range(1, n):
        z[t] = 0.5 * z[t - 1] + ez[t]
        x[t] = 0.9 * z[t - 1] + ex[t]     # X ← Z
        y[t] = 0.9 * z[t - 1] + ey[t]     # Y ← Z（X とは Z 経由でしか繋がらない）
    return x, y, z


def make_direct(n=4000, seed=0):
    """X→Y が本物で、Z は無関係。"""
    rng = np.random.default_rng(seed)
    z = np.zeros(n); x = np.zeros(n); y = np.zeros(n)
    ez, ex, ey = (rng.normal(0, 1, n) for _ in range(3))
    for t in range(1, n):
        x[t] = 0.5 * x[t - 1] + ex[t]
        y[t] = 0.9 * x[t - 1] + ey[t]     # Y ← X（本物）
        z[t] = 0.5 * z[t - 1] + ez[t]     # Z は誰とも繋がらない
    return x, y, z


if __name__ == "__main__":
    M = 8
    print("=== 未観測共通原因の検出（I(X;Y) vs I(X;Y|Z)）===\n")

    # 1) 共通原因ケース：Z→X, Z→Y（X–Y 直接なし）
    x, y, z = make_common_cause()
    # X_t と Y_t は Z_{t-1} を共有 → 揃えて Z_{t-1} で条件づける
    mi = mutual_info(x[1:], y[1:], M)
    cmi = cond_mutual_info(x[1:], y[1:], z[:-1], M)
    print("【共通原因ケース：Z が X と Y を動かす（X–Y 直接なし）】")
    print(f"  I(X;Y)      = {mi:.4f} bit   ← Z を測らないと『繋がって見える』")
    print(f"  I(X;Y | Z)  = {cmi:.4f} bit   ← Z で条件づけると消える")
    print(f"  減少率 = {100*(1-cmi/mi):.0f}%  → 見かけの関連は Z（共通原因）の影\n")

    # 2) 直接ケース：X→Y が本物、Z 無関係
    x2, y2, z2 = make_direct()
    mi2 = mutual_info(x2[:-1], y2[1:], M)                 # X_{t-1} と Y_t（本物の向き）
    cmi2 = cond_mutual_info(x2[:-1], y2[1:], z2[:-1], M)  # 無関係 Z で条件づけ
    print("【対比・直接ケース：X→Y が本物、Z は無関係】")
    print(f"  I(X;Y)      = {mi2:.4f} bit")
    print(f"  I(X;Y | Z)  = {cmi2:.4f} bit   ← 無関係な Z で条件づけても残る（本物）\n")

    print("=== 問2a の教訓 ===")
    print("  ・共通原因ケース: Z を『観測して条件づける』と偽の関連が消える＝Z が交絡だと検出できる。")
    print("  ・もし Z を測っていなければ、I(X;Y) の高さだけ見て X↔Y を本物と誤認する（未観測交絡の罠）。")
    print("  ・直接ケース: 無関係な変数で条件づけても本物の関連は残る＝『何を条件づけるか』が要。")
    print("  → 未観測の原因は『残差/条件付けの構造』＋『領域知識でZの候補を持ち込む』で捕まえる。")
