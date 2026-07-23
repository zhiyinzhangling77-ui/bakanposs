"""合成データ生成 (R&K §2.2 の結合ロジスティック写像ほか)。"""

from __future__ import annotations

import numpy as np


def coupled_logistic(
    n: int,
    coupling: float,
    r: float = 4.0,
    burn_in: int = 500,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """単方向結合ロジスティック写像 X → Y (ラグ 1 で結合)。

    x_{n+1} = r x_n (1 - x_n)
    y_{n+1} = r m_n (1 - m_n),  m_n = e·x_n + (1-e)·y_n

    coupling e > 0 で X が Y を 1 step 先行して駆動する。TE(X→Y) > TE(Y→X) が
    成り立つはず (Schreiber 2000; R&K §2.2)。
    """
    rng = np.random.default_rng(seed)
    x = np.empty(n + burn_in)
    y = np.empty(n + burn_in)
    x[0] = rng.uniform(0.1, 0.9)
    y[0] = rng.uniform(0.1, 0.9)
    e = coupling
    for k in range(n + burn_in - 1):
        x[k + 1] = r * x[k] * (1.0 - x[k])
        m = e * x[k] + (1.0 - e) * y[k]
        y[k + 1] = r * m * (1.0 - m)
    return x[burn_in:], y[burn_in:]


def independent_noise(n: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """互いに独立な一様乱数 2 本。"""
    rng = np.random.default_rng(seed)
    return rng.random(n), rng.random(n)
