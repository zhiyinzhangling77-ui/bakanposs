#!/usr/bin/env python3
"""旗131 段3 — Gaumont-Guay 2006 が印字した 3 つの связ（相関）を、**同じ標本の相関行列として成立するか**当たる。

## なぜこれをやるか

**A-2 の差分は「交絡を分離した」である。** 先行が分離していないことを、
**言葉ではなく、先行自身が印字した数で**言えるなら、その方が強い（旗127 で Dörr & Münnich に対してやった作法）。

Gaumont-Guay 2006（AFM 140, 220–235）は、生長期の日別 n=110 について次の 3 つを印字している：

  (i)  Q10Rs = 9.01 θ + 0.62      r² = 0.55   （p.228 本文・Fig. 6b）→ r(Q10, θ) = +√0.55
  (ii) Q10Rs = −0.12 Ts + 3.94    r² = 0.19   （p.228 本文・Fig. 6c）→ r(Q10, Ts) = −√0.19
  (iii)「both climate variables were highly correlated during the growing season (r = 0.82)」（p.229）

**(iii) は符号を付けずに 0.82 と印字されている＝正である。**

## 当てること（**結果を見る前に決める**）

1. **可能性の門**：3 つを同じ標本の相関行列と読んだとき、行列は半正定値か
   （det ≥ 0 か）。**det < 0 なら、3 つは同じ標本から出ていない。**
2. **符号の分岐**：r(θ,Ts) = +0.82 と −0.82 の両方で 1 を回す。
3. **通った枝でだけ**、Q10 を (θ, Ts) に同時回帰したときの標準化偏回帰係数を閉形式で出す。
   **これは「先行の数から、先行がやらなかった同時推定を復元する」ことに当たる。**

## 門①（対照）— **この算術が正しいことを、答えの分かっている場合で示す**

閉形式（β = R⁻¹ r）と、**同じ相関行列から実際に標本を発生させて最小二乗で当てた値**が
一致するかを見る。**一致しなければ、この道具の出す数は信用しない。**
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from runlog import tee_stdout  # noqa: E402

R2_THETA = 0.55     # Q10Rs ~ θ   （正の傾き 9.01）
R2_TS = 0.19        # Q10Rs ~ Ts  （負の傾き −0.12）
N = 110             # 本文「n = 110 for Rs10 and Q10Rs」


def corr_matrix(r_yth, r_yts, r_thts):
    return np.array([[1.0, r_yth, r_yts],
                     [r_yth, 1.0, r_thts],
                     [r_yts, r_thts, 1.0]])


def partials(r_yth, r_yts, r_thts):
    """Q10 を (θ, Ts) に同時回帰したときの標準化偏回帰係数と R²。"""
    R = np.array([[1.0, r_thts], [r_thts, 1.0]])
    r = np.array([r_yth, r_yts])
    beta = np.linalg.solve(R, r)
    return beta, float(beta @ r)


def control_synthetic(r_yth, r_yts, r_thts, seed=0):
    """門① — 同じ相関行列から標本を作り、最小二乗が閉形式に一致するか。"""
    C = corr_matrix(r_yth, r_yts, r_thts)
    w, V = np.linalg.eigh(C)
    if w.min() < 0:
        return None, "この枝は半正定値でないので標本を作れない（＝門①も回せない）"
    L = V @ np.diag(np.sqrt(np.clip(w, 0, None)))
    rng = np.random.default_rng(seed)
    X = (L @ rng.standard_normal((3, 400000))).T
    y, th, ts = X[:, 0], X[:, 1], X[:, 2]
    A = np.column_stack([th, ts])
    A = (A - A.mean(0)) / A.std(0)
    yz = (y - y.mean()) / y.std()
    b, *_ = np.linalg.lstsq(A, yz, rcond=None)
    return b, None


def main():
    tee_stdout("step131_partial")
    print("=== 旗131 段3 — Gaumont-Guay 2006 の 3 つの相関は同じ標本の行列として成立するか ===")
    r_yth = +np.sqrt(R2_THETA)     # 傾き 9.01 が正なので正
    r_yts = -np.sqrt(R2_TS)        # 傾き −0.12 が負なので負
    print(f"  本文から：r(Q10, θ)  = +√{R2_THETA} = {r_yth:+.4f}   （Q10Rs = 9.01θ + 0.62）")
    print(f"            r(Q10, Ts) = −√{R2_TS} = {r_yts:+.4f}   （Q10Rs = −0.12Ts + 3.94）")
    print(f"            n = {N}")

    for label, r_thts in [("印字どおり r(θ,Ts) = +0.82", +0.82),
                          ("符号を逆に読む r(θ,Ts) = −0.82", -0.82)]:
        print("\n" + "-" * 74)
        print(f"### 枝: {label}")
        C = corr_matrix(r_yth, r_yts, r_thts)
        det = float(np.linalg.det(C))
        eig = np.linalg.eigvalsh(C)
        print(f"  det = {det:+.4f}   固有値 = {np.array2string(eig, precision=4)}")
        if det < 0 or eig.min() < 0:
            print("  判定: **NOT_PSD — この 3 つは同じ標本の相関行列になりえない**")
            print("        ＝(iii) の 0.82 は Fig. 6 の n=110 とは別の標本／別の時間刻みか、"
                  "でなければ符号が違う。**どちらかは、この場では決められない。**")
            continue
        print("  判定: PSD — 同じ標本の相関行列として成立しうる")
        beta, r2 = partials(r_yth, r_yts, r_thts)
        print(f"  閉形式の標準化偏回帰係数: β_θ = {beta[0]:+.4f}  β_Ts = {beta[1]:+.4f}   R² = {r2:.4f}")
        print(f"  単変量（先行が印字したもの）:  r_θ = {r_yth:+.4f}  r_Ts = {r_yts:+.4f}")
        flip = np.sign(beta) != np.sign([r_yth, r_yts])
        print(f"  符号が反転する変数: θ={bool(flip[0])}  Ts={bool(flip[1])}")
        b_emp, err = control_synthetic(r_yth, r_yts, r_thts)
        if err:
            print(f"  門①: {err}")
        else:
            d = np.abs(b_emp - beta).max()
            print(f"  門①（合成 40 万標本の最小二乗）: β = {b_emp[0]:+.4f}, {b_emp[1]:+.4f}  "
                  f"閉形式との最大差 = {d:.4f} → {'一致（算術は正しい）' if d < 0.01 else '★不一致：この道具を信用しない'}")
    print("\n" + "=" * 74)
    print("  ※ ここで出した偏回帰係数は **先行の本文には無い**。"
          "先行が印字した相関から復元したものであり、先行の主張ではない。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
