"""O-information の正準ケースとサブシステム解析の検証。"""

from __future__ import annotations

import numpy as np
import pandas as pd

from japanflux_pn.config import AnalysisConfig, RK_VARS
from japanflux_pn import information_theory as it
from japanflux_pn import oinfo_analysis as oi
from japanflux_pn.preprocess import PreprocessResult


def test_o_information_redundancy_positive():
    """3 変数コピー X1=X2=X3 → 冗長支配 Ω>0。"""
    rng = np.random.default_rng(0)
    m = 4
    x = rng.integers(0, m, size=6000)
    omega = it.o_information_indices([x.copy(), x.copy(), x.copy()], m)
    assert omega > 0.3


def test_o_information_synergy_negative():
    """3 ビット・パリティ X3=X1⊕X2 → 相乗支配 Ω<0。"""
    rng = np.random.default_rng(1)
    n = 9000
    x1 = rng.integers(0, 2, size=n)
    x2 = rng.integers(0, 2, size=n)
    x3 = np.bitwise_xor(x1, x2)
    omega = it.o_information_indices([x1, x2, x3], 2)
    assert omega < -0.3


def test_o_information_independent_near_zero_vs_surrogate():
    """独立 3 変数 → Ω はシャッフルヌルと区別できない (|z| 小)。"""
    rng = np.random.default_rng(2)
    m = 6
    n = 4000
    cols = [it.digitize_series(rng.normal(size=n), m) for _ in range(4)]
    omega = it.o_information_indices(cols, m, correct=True)
    stats = it.surrogate_o_information_stats(cols, m, 60, 2.36,
                                             np.random.default_rng(3), correct=True)
    z = (omega - stats["mu"]) / stats["sigma"]
    assert abs(z) < 3.0                     # 独立系は有意な冗長/相乗を示さない


def test_o_equals_tc_minus_dtc():
    """Ω = TC − DTC の恒等式 (DTC = ΣH(X_-i) − (n-1)H(X))。"""
    rng = np.random.default_rng(4)
    m = 5
    n = 3000
    cols = [rng.integers(0, m, size=n) for _ in range(3)]
    # 相関を持たせる
    cols[1] = (cols[0] + rng.integers(0, 2, size=n)) % m
    omega = it.o_information_indices(cols, m)
    tc = it.total_correlation_indices(cols, m)
    h_all = it._entropy_of_indices(cols, m)
    dtc = sum(it._entropy_of_indices([cols[j] for j in range(3) if j != i], m)
              for i in range(3)) - (3 - 1) * h_all
    assert abs(omega - (tc - dtc)) < 1e-9


def _synth_pre(n=4000):
    """土壌–呼吸系を相乗的に、フラックス系を Rg 冗長に作った合成。"""
    cfg = AnalysisConfig(n_surrogates=40, seed=0)
    rng = np.random.default_rng(7)
    rg = rng.normal(size=n)
    data = {v: rng.normal(size=n) for v in RK_VARS}
    data["Rg"] = rg
    data["gH"] = 0.8 * rg + 0.3 * rng.normal(size=n)     # Rg 冗長
    data["gLE"] = 0.8 * rg + 0.3 * rng.normal(size=n)    # Rg 冗長
    data["Ta"] = 0.3 * rg + rng.normal(size=n)
    # GER をパリティ的に Ta と th の相互作用で (相乗を仕込む)
    th = rng.integers(0, 2, size=n)
    ta_hi = (data["Ta"] > np.median(data["Ta"])).astype(int)
    data["th"] = th.astype(float) + 0.1 * rng.normal(size=n)
    data["GER"] = np.bitwise_xor(th, ta_hi).astype(float) + 0.1 * rng.normal(size=n)
    grid = pd.date_range("2020-07-01", periods=n, freq=pd.Timedelta(minutes=30))
    frame = pd.DataFrame(data, index=grid)[RK_VARS]
    valid = pd.Series(True, index=grid)
    return PreprocessResult(anomaly=frame, valid=valid, site="S", year=2020,
                            month=7, config=cfg, months=[7])


def test_subsystem_analysis_runs_and_separates():
    pre = _synth_pre()
    tbl = oi.o_information_subsystems(pre, obins=6)
    assert len(tbl) == len(oi.SUBSYSTEMS)
    row = tbl[tbl["vars"] == "Rg,Ta,gH,gLE"].iloc[0]
    # エネルギー系 (gH,gLE が Rg 冗長) は Ω>0 側
    assert row["Omega"] > row["null_mu"]
