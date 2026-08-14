"""旗3：lag スキャン。TE を lag=1..L で計算し、ピークが"真の遅れ"に立つのを見る。

step1 は真の遅れを知って lag を合わせた。ここでは lag を知らないふりをして総当たりし、
「TE が最大になる lag ＝ X が Y を動かす遅れ」を当てられるかを確認する。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from transfer_entropy_step1 import transfer_entropy, make_synthetic

if __name__ == "__main__":
    M, DELAY, LMAX = 8, 3, 8
    x, y = make_synthetic(n=4000, delay=DELAY, coupling=0.7, seed=0)

    print(f"=== lag スキャン（真の遅れ={DELAY}, m={M}）===")
    print(f"  {'lag':>4}  {'TE(X→Y)':>9}  {'TE(Y→X)':>9}")
    best_lag, best_te = 0, -1.0
    for lag in range(1, LMAX + 1):
        a = transfer_entropy(x, y, m=M, lag=lag)   # X→Y
        b = transfer_entropy(y, x, m=M, lag=lag)   # Y→X
        star = " ←最大" if a > best_te else ""
        if a > best_te:
            best_te, best_lag = a, lag
        print(f"  {lag:>4}  {a:9.4f}  {b:9.4f}{star}")

    print(f"\n  TE(X→Y) が最大の lag = {best_lag}（真の遅れ={DELAY}）")
    print("  判定: " + ("✅ ピークが真の遅れに一致＝遅れを当てられた"
                        if best_lag == DELAY else
                        "⚠ ピークがズレた（結合の強さ/ノイズ/ビン数を確認）"))
    # 読み方: X→Y の曲線が lag=3 で山になり、他の lag では低い。
    #  Y→X はどの lag でもバイアスの床付近（＝向きが無い）。
    #  ＝「向き」と「遅れ」を同時に読める。
