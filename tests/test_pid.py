"""PID (Williams-Beer I_min) の正準 3 ケースと共通駆動分解の検証。"""

from __future__ import annotations

import numpy as np
import pandas as pd

from japanflux_pn.config import AnalysisConfig, RK_VARS
from japanflux_pn import information_theory as it
from japanflux_pn import pid_analysis as pid
from japanflux_pn.preprocess import PreprocessResult


def test_pid_pure_redundancy():
    """T=S1=S2 (完全同一) → 冗長 R が支配、U と S ≈ 0。"""
    rng = np.random.default_rng(0)
    m = 4
    t = rng.integers(0, m, size=6000)
    r = it.pid_williams_beer(t.copy(), t.copy(), t.copy(), m)
    assert r["R"] > 0.5 * r["I_joint"]          # 冗長が大半
    assert abs(r["U1"]) < 1e-6 and abs(r["U2"]) < 1e-6
    assert abs(r["S"]) < 1e-6


def test_pid_pure_synergy_xor():
    """T = S1 XOR S2 (独立ビット) → 相乗 S が支配、R≈0, U≈0。"""
    rng = np.random.default_rng(1)
    n = 8000
    s1 = rng.integers(0, 2, size=n)
    s2 = rng.integers(0, 2, size=n)
    t = np.bitwise_xor(s1, s2)
    r = it.pid_williams_beer(t, s1, s2, 2)
    assert r["S"] > 0.8 * r["I_joint"]          # ほぼ全部が相乗
    assert r["R"] < 0.05
    assert r["U1"] < 0.05 and r["U2"] < 0.05


def test_pid_pure_unique():
    """T=S1, S2 独立 → U1 が支配、R≈0, S≈0。"""
    rng = np.random.default_rng(2)
    m = 5
    n = 8000
    t = rng.integers(0, m, size=n)
    s2 = rng.integers(0, m, size=n)
    r = it.pid_williams_beer(t.copy(), t.copy(), s2, m)
    assert r["U1"] > 0.8 * r["I_joint"]
    assert r["R"] < 0.05
    assert abs(r["S"]) < 0.1


def test_specific_information_sums_to_mi():
    """Σ_s p(s) i(s;A) = I(T;A) の整合。"""
    rng = np.random.default_rng(3)
    m = 6
    n = 5000
    t = rng.integers(0, m, size=n)
    a = (t + rng.integers(0, 2, size=n)) % m       # T と相関
    i_spec, p_t = it.specific_information(t, a, m)
    mi_from_spec = float(np.sum(p_t * i_spec))
    mi_direct = it.mutual_information_indices(t, a, m)
    assert abs(mi_from_spec - mi_direct) < 1e-9


def _synth_pre(n=4000):
    """Rg が Ta・VPD・gLE を共通駆動、GER は Ta 独自の熱情報も持つ合成。"""
    cfg = AnalysisConfig(seed=0)
    rng = np.random.default_rng(7)
    rg = rng.normal(size=n)
    data = {v: rng.normal(size=n) for v in RK_VARS}
    data["Rg"] = rg
    data["VPD"] = 0.7 * rg + 0.5 * rng.normal(size=n)   # 主に Rg
    data["gLE"] = 0.7 * rg + 0.5 * rng.normal(size=n)   # 主に Rg (gLE↔VPD は Rg 冗長)
    data["Ta"] = 0.2 * rg + rng.normal(size=n)          # 大きな非 Rg 成分
    data["GER"] = data["Ta"] + 0.3 * rng.normal(size=n)  # Ta に連動 (非 Rg 熱情報)
    grid = pd.date_range("2020-07-01", periods=n, freq=pd.Timedelta(minutes=30))
    frame = pd.DataFrame(data, index=grid)[RK_VARS]
    valid = pd.Series(True, index=grid)
    return PreprocessResult(anomaly=frame, valid=valid, site="S", year=2020,
                            month=7, config=cfg, months=[7])


def test_pid_with_driver_separates_redundant_and_unique():
    pre = _synth_pre()
    res = pid.pid_with_driver(pre, driver="Rg")
    t = res.table

    def row(y, x):
        return t[(t["target"] == y) & (t["source"] == x)].iloc[0]

    # gLE←VPD: 両方 Rg 共通駆動 → VPD の情報は Rg と冗長 (R/I 高い)
    assert row("gLE", "VPD")["redundancy_frac"] > 0.6
    # GER←Ta: Ta は Rg 外の熱情報を運ぶ → X 固有が残る (R/I 低い)
    assert row("GER", "Ta")["redundancy_frac"] < 0.4
    assert row("GER", "Ta")["U_source_pct"] > row("gLE", "VPD")["U_source_pct"]
