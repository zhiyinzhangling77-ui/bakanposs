"""エントロピー・相互情報量・transfer entropy とサロゲート有意性検定。

数値カーネル。変数名やサイトを一切知らず、離散化済みの 1 次元ビンインデックス
配列 (と任意でグリッド位置 ``step_index``) だけを受け取る。

推定方式は Ruddell & Kumar (2009) に従う固定区間ビン (m=11)。多次元同時分布は
各次元のビンインデックスを混合基数で 1 整数へエンコードし ``np.bincount`` で
度数化する。エントロピーは自然対数 ``H = -Σ p log p``。正規化 (``/log m``) は
出力側で行い、ここでは生値 (nats) を返す。
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# 離散化
# ---------------------------------------------------------------------------
def bin_edges(x: np.ndarray, n_bins: int) -> np.ndarray:
    """系列 x の [min, max] を n_bins 等分する境界 (長さ n_bins+1)。

    範囲が退化 (min==max) している場合は微小幅を与えて 1 本のビンに落とす。
    """
    x = np.asarray(x, dtype=float)
    lo, hi = np.nanmin(x), np.nanmax(x)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        hi = lo + 1.0
    return np.linspace(lo, hi, n_bins + 1)


def to_bin_indices(x: np.ndarray, edges: np.ndarray, n_bins: int) -> np.ndarray:
    """x を bin index 0..n_bins-1 に離散化する。

    内側境界 ``edges[1:-1]`` で :func:`np.digitize` し、両端を含めてクリップする
    (最大値は最終ビンに入る)。ラグ版の同一変数は必ず同じ ``edges`` を共有させる
    こと (``X_{t-τ}`` と ``X_t`` が同一格子)。
    """
    idx = np.digitize(np.asarray(x, dtype=float), edges[1:-1])
    return np.clip(idx, 0, n_bins - 1).astype(np.int64)


def digitize_series(x: np.ndarray, n_bins: int) -> np.ndarray:
    """生系列を独自の [min,max] レンジで bin index 化する簡便版。"""
    return to_bin_indices(x, bin_edges(x, n_bins), n_bins)


# ---------------------------------------------------------------------------
# エントロピー (ビンインデックスから)
# ---------------------------------------------------------------------------
def _entropy_of_indices(cols: list[np.ndarray], n_bins: int) -> float:
    """同時エントロピー H(cols...) を混合基数エンコード + bincount で計算。

    cols: それぞれ同じ長さ N の bin index 配列 (値域 0..n_bins-1)。
    自然対数 (nats)。
    """
    n = len(cols[0])
    if n == 0:
        return 0.0
    code = np.zeros(n, dtype=np.int64)
    for c in cols:
        code = code * n_bins + c
    counts = np.bincount(code)
    counts = counts[counts > 0]
    p = counts / n
    return float(-np.sum(p * np.log(p)))


def shannon_entropy(x: np.ndarray, n_bins: int) -> float:
    """生系列の Shannon エントロピー H(X) [nats]。"""
    return _entropy_of_indices([digitize_series(x, n_bins)], n_bins)


# ---------------------------------------------------------------------------
# 相互情報量 (ゼロラグ)
# ---------------------------------------------------------------------------
def mutual_information_indices(xi: np.ndarray, yi: np.ndarray, n_bins: int) -> float:
    """I(X, Y) = H(X) + H(Y) - H(X, Y) [nats]。入力は bin index。"""
    hx = _entropy_of_indices([xi], n_bins)
    hy = _entropy_of_indices([yi], n_bins)
    hxy = _entropy_of_indices([xi, yi], n_bins)
    return hx + hy - hxy


def mutual_information(x: np.ndarray, y: np.ndarray, n_bins: int) -> float:
    """生系列版の I(X, Y)。"""
    return mutual_information_indices(
        digitize_series(x, n_bins), digitize_series(y, n_bins), n_bins
    )


# ---------------------------------------------------------------------------
# Transfer entropy (Knuth 2005 form, R&K 式 (5))
# ---------------------------------------------------------------------------
def _lag_triples(
    xi: np.ndarray,
    yi: np.ndarray,
    tau: int,
    step_index: np.ndarray | None,
    gap_guard: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(X_{t-τ}, Y_{t-1}, Y_t) の整列済み三つ組を返す。

    時間順に詰めた配列上で位置ラグを取る。``gap_guard`` かつ ``step_index`` 有りの
    場合、レギュラ格子上でちょうど連続する三つ組 (k_t-k_{t-1}=1 かつ k_t-k_{t-τ}=τ)
    のみ残す (欠測ギャップ跨ぎを除外)。
    """
    n = len(xi)
    if tau < 1 or tau >= n:
        empty = np.empty(0, dtype=np.int64)
        return empty, empty, empty
    x_tl = xi[: n - tau]         # X_{t-τ},  t = τ..n-1
    y_m1 = yi[tau - 1 : n - 1]   # Y_{t-1}
    y_t = yi[tau:]               # Y_t
    if gap_guard and step_index is not None:
        k = step_index
        kt = k[tau:]
        ktm1 = k[tau - 1 : n - 1]
        ktl = k[: n - tau]
        keep = (kt - ktm1 == 1) & (kt - ktl == tau)
        x_tl, y_m1, y_t = x_tl[keep], y_m1[keep], y_t[keep]
    return x_tl, y_m1, y_t


def transfer_entropy_indices(
    xi: np.ndarray,
    yi: np.ndarray,
    tau: int,
    n_bins: int,
    step_index: np.ndarray | None = None,
    gap_guard: bool = True,
) -> float:
    """T(X→Y, τ) [nats]。入力は bin index 配列。

    T = H(X_{t-τ}, Y_{t-1}) + H(Y_t, Y_{t-1}) - H(Y_{t-1}) - H(X_{t-τ}, Y_t, Y_{t-1})
    """
    x_tl, y_m1, y_t = _lag_triples(xi, yi, tau, step_index, gap_guard)
    if len(y_t) == 0:
        return np.nan
    h_xtl_ym1 = _entropy_of_indices([x_tl, y_m1], n_bins)
    h_yt_ym1 = _entropy_of_indices([y_t, y_m1], n_bins)
    h_ym1 = _entropy_of_indices([y_m1], n_bins)
    h_xtl_yt_ym1 = _entropy_of_indices([x_tl, y_t, y_m1], n_bins)
    return h_xtl_ym1 + h_yt_ym1 - h_ym1 - h_xtl_yt_ym1


def transfer_entropy(
    x: np.ndarray, y: np.ndarray, tau: int, n_bins: int, gap_guard: bool = False
) -> float:
    """生系列版の T(X→Y, τ)。テスト・単発利用向け。"""
    return transfer_entropy_indices(
        digitize_series(x, n_bins),
        digitize_series(y, n_bins),
        tau,
        n_bins,
        gap_guard=gap_guard,
    )


def te_lag_curve(
    xi: np.ndarray,
    yi: np.ndarray,
    lags: list[int],
    n_bins: int,
    step_index: np.ndarray | None = None,
    gap_guard: bool = True,
) -> np.ndarray:
    """各ラグ τ における T(X→Y, τ) の配列。"""
    return np.array(
        [
            transfer_entropy_indices(xi, yi, t, n_bins, step_index, gap_guard)
            for t in lags
        ]
    )


# ---------------------------------------------------------------------------
# サロゲート有意性検定 (shuffled surrogate, R&K §A2)
# ---------------------------------------------------------------------------
def surrogate_te_stats(
    xi: np.ndarray,
    yi: np.ndarray,
    lags: list[int],
    n_bins: int,
    n_surrogates: int,
    c: float,
    rng: np.random.Generator,
) -> dict[str, np.ndarray]:
    """ラグ毎の サロゲート T 分布から (μ_ss, σ_ss, Δ) を返す [nats]。

    X, Y を独立にランダム置換して時間結合を破壊し、同じ推定器で T を計算する
    (サロゲートはギャップ意味を持たないので gap_guard=False)。しきい値は
    Δ(T) = μ_ss + c·σ_ss。戻り値の配列は ``lags`` と同順。
    """
    n = len(xi)
    m = len(lags)
    samples = np.empty((n_surrogates, m), dtype=float)
    for s in range(n_surrogates):
        xs = xi[rng.permutation(n)]
        ys = yi[rng.permutation(n)]
        for j, t in enumerate(lags):
            samples[s, j] = transfer_entropy_indices(
                xs, ys, t, n_bins, step_index=None, gap_guard=False
            )
    mu = np.nanmean(samples, axis=0)
    sigma = np.nanstd(samples, axis=0)
    return {"mu": mu, "sigma": sigma, "threshold": mu + c * sigma}


def surrogate_mi_stats(
    xi: np.ndarray,
    yi: np.ndarray,
    n_bins: int,
    n_surrogates: int,
    c: float,
    rng: np.random.Generator,
) -> dict[str, float]:
    """サロゲート MI 分布から (μ_ss, σ_ss, Δ) を返す [nats]。"""
    n = len(xi)
    samples = np.empty(n_surrogates, dtype=float)
    for s in range(n_surrogates):
        samples[s] = mutual_information_indices(
            xi[rng.permutation(n)], yi[rng.permutation(n)], n_bins
        )
    mu = float(np.mean(samples))
    sigma = float(np.std(samples))
    return {"mu": mu, "sigma": sigma, "threshold": mu + c * sigma}
