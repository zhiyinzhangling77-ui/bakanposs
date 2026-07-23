"""preprocess のアノマリ・欠測処理の単体テスト (ファイル不要)。"""

from __future__ import annotations

import numpy as np
import pandas as pd

from japanflux_pn.config import AnalysisConfig, RK_VARS
from japanflux_pn import preprocess as pp


def _make_grid(n_days: int, sp: int = 48) -> pd.DatetimeIndex:
    start = pd.Timestamp("2020-07-01")
    freq = pd.Timedelta(minutes=24 * 60 // sp)
    return pd.date_range(start, periods=n_days * sp, freq=freq)


def test_forward_window_anomaly_constant_diurnal():
    """時刻帯ごとに一定 (日変化のみ) の系列 → アノマリはほぼ 0。"""
    sp, w = 48, 5
    n_days = 10
    diurnal = np.sin(np.linspace(0, 2 * np.pi, sp, endpoint=False))
    values = np.tile(diurnal, n_days)
    anom = pp._forward_window_anomaly(values, sp, w)
    # 窓が満たせる先頭側はほぼ 0
    head = anom[: (n_days - w) * sp]
    assert np.nanmax(np.abs(head)) < 1e-10
    # 末尾 (w-1) 日は窓不足で NaN
    assert np.all(np.isnan(anom[-(w - 1) * sp :]))


def test_forward_window_removes_forward_mean():
    """既知の値で前方窓平均が正しく引かれているか。"""
    sp, w = 48, 5
    n = sp * (w + 1)
    values = np.arange(n, dtype=float)
    anom = pp._forward_window_anomaly(values, sp, w)
    t = 0
    expected = values[t] - np.mean([values[t + d * sp] for d in range(w)])
    assert np.isclose(anom[t], expected)


def test_missing_propagates_to_nan():
    """窓内に 1 点でも欠測があればアノマリ NaN。"""
    sp, w = 48, 5
    n = sp * (w + 2)
    values = np.ones(n, dtype=float)
    values[3 * sp + 7] = np.nan  # ある時刻帯 7 の day3 を欠測に
    anom = pp._forward_window_anomaly(values, sp, w)
    # 時刻帯 7 で、その欠測を窓に含む t (day0..day3) が NaN になる
    assert np.isnan(anom[7])          # day0 の窓は day0..day4 を含む → 欠測含む
    assert np.isnan(anom[3 * sp + 7])


def test_compute_anomaly_listwise_deletion():
    """1 変数でも欠測なら valid=False (listwise)。"""
    cfg = AnalysisConfig()
    grid = _make_grid(12)
    rng = np.random.default_rng(0)
    raw = pd.DataFrame(
        {v: rng.random(len(grid)) for v in RK_VARS}, index=grid
    )
    raw.iloc[100, raw.columns.get_loc("VPD")] = np.nan
    anom, valid = pp.compute_anomaly(raw, cfg)
    assert not valid.iloc[100]
    assert valid.dtype == bool
    assert list(anom.columns) == RK_VARS


def test_step_index_matches_grid_positions():
    """step_index が格子上の整数位置と一致する。"""
    cfg = AnalysisConfig()
    grid = _make_grid(12)
    rng = np.random.default_rng(1)
    raw = pd.DataFrame({v: rng.random(len(grid)) for v in RK_VARS}, index=grid)
    anom, valid = pp.compute_anomaly(raw, cfg)
    res = pp.PreprocessResult(
        anomaly=anom, valid=valid, site="TEST", year=2020, month=7, config=cfg
    )
    k = res.step_index
    # 単調増加、かつ有効行数と一致
    assert len(k) == res.n_points
    assert np.all(np.diff(k) > 0)
