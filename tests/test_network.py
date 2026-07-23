"""network の分類ロジックと build_network 統合テスト。"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from japanflux_pn.config import AnalysisConfig, RK_VARS
from japanflux_pn import network as nw
from japanflux_pn.preprocess import PreprocessResult
from synthetic import coupled_logistic


# ---------------------------------------------------------------------------
# 分類ロジック
# ---------------------------------------------------------------------------
def test_classify_types():
    assert nw.classify_coupling(t_significant=False, i_significant=True, tz=np.nan) == 1
    assert nw.classify_coupling(t_significant=True, i_significant=True, tz=0.5) == 2
    assert nw.classify_coupling(t_significant=True, i_significant=True, tz=3.0) == 3
    assert nw.classify_coupling(t_significant=False, i_significant=False, tz=np.nan) == 4


def test_first_significant_peak_basic():
    lags = [1, 2, 3, 4, 5]
    curve = np.array([0.1, 0.9, 0.4, 0.2, 0.1])
    thr = np.full(5, 0.5)
    tau, j = nw.first_significant_peak(curve, thr, lags)
    assert tau == 2 and j == 1


def test_first_significant_peak_none():
    lags = [1, 2, 3]
    curve = np.array([0.1, 0.2, 0.15])
    thr = np.full(3, 0.5)
    tau, j = nw.first_significant_peak(curve, thr, lags)
    assert tau is None and j is None


def test_first_significant_peak_handles_nan():
    lags = [1, 2, 3]
    curve = np.array([np.nan, 0.9, 0.3])
    thr = np.full(3, 0.5)
    tau, j = nw.first_significant_peak(curve, thr, lags)
    assert tau == 2


def test_peak_min_run_rejects_isolated_crossing():
    """min_run=2 で単発クロスは棄却、連続有意帯は採用。"""
    lags = [1, 2, 3, 4, 5]
    thr = np.full(5, 0.5)
    isolated = np.array([0.1, 0.9, 0.1, 0.1, 0.1])   # lag2 のみ有意 (単発)
    assert nw.first_significant_peak(isolated, thr, lags, min_run=1)[0] == 2
    assert nw.first_significant_peak(isolated, thr, lags, min_run=2)[0] is None

    coherent = np.array([0.1, 0.9, 0.8, 0.1, 0.1])   # lag2-3 が連続有意
    assert nw.first_significant_peak(coherent, thr, lags, min_run=2)[0] == 2


# ---------------------------------------------------------------------------
# build_network 統合 (合成データ, 既知の結合)
# ---------------------------------------------------------------------------
def _synthetic_pre(n: int = 1000) -> PreprocessResult:
    """11 変数フレーム。VPD→gLE を lag1 で結合、他は独立乱数。"""
    cfg = AnalysisConfig(lag_max=6, n_surrogates=30, seed=0)
    rng = np.random.default_rng(42)
    driver, target = coupled_logistic(n, coupling=0.45, seed=11)
    grid = pd.date_range("2020-07-01", periods=n,
                         freq=pd.Timedelta(minutes=30))
    data = {v: rng.random(n) for v in RK_VARS}
    data["VPD"] = driver
    data["gLE"] = target
    frame = pd.DataFrame(data, index=grid)[RK_VARS]
    valid = pd.Series(True, index=grid)
    return PreprocessResult(anomaly=frame, valid=valid, site="SYNTH",
                            year=2020, month=7, config=cfg)


def test_build_network_shapes_and_coupling():
    pre = _synthetic_pre()
    net = nw.build_network(pre)

    # 形状と対称性
    for mat in (net.AI, net.ATz, net.Gamma, net.ctype):
        assert mat.shape == (11, 11)
        assert list(mat.index) == RK_VARS and list(mat.columns) == RK_VARS
    ai = net.AI.to_numpy(dtype=float)
    assert np.allclose(ai, ai.T, equal_nan=True)  # AI 対称

    # 既知の結合 VPD→gLE が検出される
    assert int(net.ctype.loc["VPD", "gLE"]) in (2, 3)
    assert np.isfinite(net.Gamma.loc["VPD", "gLE"])
    assert net.Gamma.loc["VPD", "gLE"] == pytest.approx(0.5, abs=1e-9)  # lag1 = 0.5h
    # 逆向き gLE→VPD は駆動されていないので強制/フィードバックにならない
    assert int(net.ctype.loc["gLE", "VPD"]) in (1, 4)

    # 独立ペアはほぼ非結合 (type 4 が多数派)
    off = [int(net.ctype.loc[a, b]) for a in RK_VARS for b in RK_VARS if a != b]
    assert off.count(4) > len(off) // 2
