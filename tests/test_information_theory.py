"""information_theory の数値検証 (合成データ)。"""

from __future__ import annotations

import numpy as np
import pytest

from japanflux_pn import information_theory as it
from synthetic import coupled_logistic, independent_noise


M = 11  # ビン数


# ---------------------------------------------------------------------------
# エントロピー・MI の基本恒等式
# ---------------------------------------------------------------------------
def test_uniform_entropy_near_log_m():
    """m 個のビンに一様分布 → H ≈ log(m) (正規化で ≈ 1)。"""
    rng = np.random.default_rng(0)
    x = rng.random(20000)
    h = it.shannon_entropy(x, M)
    assert h / np.log(M) == pytest.approx(1.0, abs=0.02)


def test_mi_self_equals_entropy():
    """I(X, X) = H(X)。"""
    rng = np.random.default_rng(1)
    x = rng.random(5000)
    xi = it.digitize_series(x, M)
    hx = it._entropy_of_indices([xi], M)
    mi = it.mutual_information_indices(xi, xi, M)
    assert mi == pytest.approx(hx, rel=1e-9)


def test_mi_independent_near_zero():
    """独立系列の MI は正のバイアス程度で ~0 (log m 正規化で小)。"""
    x, y = independent_noise(3000, seed=2)
    mi = it.mutual_information(x, y, M) / np.log(M)
    assert mi < 0.05


def test_mi_symmetric():
    x, y = coupled_logistic(3000, coupling=0.4, seed=3)
    xi, yi = it.digitize_series(x, M), it.digitize_series(y, M)
    assert it.mutual_information_indices(xi, yi, M) == pytest.approx(
        it.mutual_information_indices(yi, xi, M), rel=1e-9
    )


# ---------------------------------------------------------------------------
# Transfer entropy: 方向性・ピークラグ・有意性
# ---------------------------------------------------------------------------
def test_te_directional_asymmetry():
    """単方向結合 X→Y: T(X→Y, 1) > T(Y→X, 1)。"""
    x, y = coupled_logistic(4000, coupling=0.4, seed=4)
    xi, yi = it.digitize_series(x, M), it.digitize_series(y, M)
    t_xy = it.transfer_entropy_indices(xi, yi, 1, M, gap_guard=False)
    t_yx = it.transfer_entropy_indices(yi, xi, 1, M, gap_guard=False)
    assert t_xy > t_yx
    assert t_xy > 0.05  # 実質的な情報流があること


def test_te_peak_at_lag_one():
    """結合はラグ 1 step で最大。"""
    x, y = coupled_logistic(4000, coupling=0.4, seed=5)
    xi, yi = it.digitize_series(x, M), it.digitize_series(y, M)
    lags = list(range(1, 11))
    curve = it.te_lag_curve(xi, yi, lags, M, gap_guard=False)
    assert lags[int(np.argmax(curve))] == 1


def test_te_significant_over_surrogate():
    """結合 X→Y は τ=1 でサロゲートしきい値を超え、逆向きは超えない。"""
    x, y = coupled_logistic(4000, coupling=0.4, seed=6)
    xi, yi = it.digitize_series(x, M), it.digitize_series(y, M)
    rng = np.random.default_rng(100)
    t_xy = it.transfer_entropy_indices(xi, yi, 1, M, gap_guard=False)
    stats = it.surrogate_te_stats(xi, yi, [1], M, n_surrogates=100, c=2.36, rng=rng)
    assert t_xy > stats["threshold"][0]

    # 逆向き: 駆動されていない X の未来は Y の過去から予測できない
    t_yx = it.transfer_entropy_indices(yi, xi, 1, M, gap_guard=False)
    rng2 = np.random.default_rng(101)
    stats_yx = it.surrogate_te_stats(yi, xi, [1], M, n_surrogates=100, c=2.36, rng=rng2)
    assert t_yx <= stats_yx["threshold"][0]


def test_te_independent_not_significant():
    """独立ノイズ 2 本は全ラグで有意にならない。"""
    x, y = independent_noise(2000, seed=7)
    xi, yi = it.digitize_series(x, M), it.digitize_series(y, M)
    lags = list(range(1, 6))
    curve = it.te_lag_curve(xi, yi, lags, M, gap_guard=False)
    rng = np.random.default_rng(200)
    stats = it.surrogate_te_stats(xi, yi, lags, M, n_surrogates=100, c=2.36, rng=rng)
    assert np.all(curve <= stats["threshold"])


def test_shuffle_destroys_te():
    """時間並べ替えで結合 (TE) が消える。

    m=11 の 3 次元ヒストグラムは有限標本で正のバイアス (floor) を持つため、
    シャッフル後の TE は 0 ではなくサロゲート平均 μ_ss 付近に落ちる。よって
    絶対値 0 ではなく「サロゲートしきい値を下回り、実測より十分小さい」で判定する。
    """
    x, y = coupled_logistic(4000, coupling=0.4, seed=8)
    xi, yi = it.digitize_series(x, M), it.digitize_series(y, M)
    rng = np.random.default_rng(9)
    xs = xi[rng.permutation(len(xi))]
    ys = yi[rng.permutation(len(yi))]
    t_shuf = it.transfer_entropy_indices(xs, ys, 1, M, gap_guard=False)
    t_real = it.transfer_entropy_indices(xi, yi, 1, M, gap_guard=False)

    stats = it.surrogate_te_stats(xi, yi, [1], M, n_surrogates=100, c=2.36, rng=rng)
    assert t_shuf < stats["threshold"][0]          # シャッフルは有意でない
    assert t_shuf == pytest.approx(stats["mu"][0], abs=4 * stats["sigma"][0])  # μ_ss 付近
    assert t_real > stats["threshold"][0]           # 実測は有意
    assert t_shuf < 0.5 * t_real                    # 実測より十分小さい


# ---------------------------------------------------------------------------
# gap_guard: ギャップ跨ぎ三つ組の除外
# ---------------------------------------------------------------------------
def test_gap_guard_drops_boundary_triples():
    """step_index が不連続な箇所の三つ組が除外される。"""
    xi = np.array([0, 1, 2, 3, 4, 5, 6, 7], dtype=np.int64)
    yi = xi.copy()
    # k に 1 箇所ギャップ (3→10) を入れる
    step = np.array([0, 1, 2, 3, 10, 11, 12, 13], dtype=np.int64)
    x_tl, y_m1, y_t = it._lag_triples(xi, yi, tau=1, step_index=step, gap_guard=True)
    # 連続する隣接ペアのみ: (0,1)(1,2)(2,3) と (4,5)(5,6)(6,7) → 6 組、境界(3,4)除外
    assert len(y_t) == 6
