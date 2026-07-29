"""climate_response の順位相関と年別指標の検証 (合成)。"""

from __future__ import annotations

import numpy as np
import pandas as pd

from japanflux_pn.config import AnalysisConfig, RK_VARS
from japanflux_pn import climate_response as cr


def test_spearman_monotonic():
    x = np.array([1, 2, 3, 4, 5], dtype=float)
    assert cr._spearman(x, 2 * x + 1) > 0.99          # 単調増加 → +1
    assert cr._spearman(x, -x) < -0.99                # 単調減少 → -1
    assert abs(cr._spearman(x, np.array([3, 1, 4, 1, 5.0]))) < 0.9


def test_spearman_p_significance():
    """強い単調関係は小さい p、無相関は大きい p。"""
    x = np.arange(12, dtype=float)
    r_mono, p_mono = cr._spearman_p(x, 3 * x - 2, n_perm=2000)
    assert r_mono > 0.99 and p_mono < 0.01
    rng = np.random.default_rng(0)
    r_rand, p_rand = cr._spearman_p(x, rng.permutation(x), n_perm=2000)
    assert p_rand > 0.10                      # 無相関は有意でない


def test_spearman_constant_is_nan():
    x = np.array([1.0, 2, 3])
    assert np.isnan(cr._spearman(x, np.array([5.0, 5, 5])))


def _raw_year(year, vpd_level, n_days=67, seed=0):
    """1 年分の生データ。VPD の平均レベルを year ごとに変え、乾燥年ほど
    Rg→GEP 結合を弱める（気孔閉鎖の模擬）。"""
    idx = pd.date_range(f"{year}-07-01", periods=n_days * 48, freq="30min")
    rng = np.random.default_rng(seed)
    n = len(idx)
    rg = rng.normal(size=n)
    data = {v: rng.normal(size=n) for v in RK_VARS}
    data["Rg"] = rg
    data["VPD"] = vpd_level + 0.5 * rng.normal(size=n)     # 平均レベルが年で違う
    data["th"] = -vpd_level + rng.normal(size=n)           # 乾燥年ほど土壌水分低
    data["Ta"] = 0.5 * rg + rng.normal(size=n)
    # 乾燥(高VPD)ほど Rg→GEP を弱める: 結合係数 ∝ (1/vpd_level)
    coup = max(0.1, 1.5 - 0.5 * vpd_level)
    data["GEP"] = coup * rg + rng.normal(size=n)
    df = pd.DataFrame(data, index=idx)[RK_VARS]
    return df


def test_year_metrics_and_stress_response():
    cfg = AnalysisConfig()
    # 4 年、VPD レベルを 0.5→2.0 に上げる（乾燥化）
    frames = []
    for i, vpd in enumerate([0.5, 1.0, 1.5, 2.0]):
        frames.append(_raw_year(2000 + i, vpd, seed=i))
    raw_all = pd.concat(frames).sort_index()

    rows = []
    for i in range(4):
        r = cr.year_metrics(raw_all, 2000 + i, [7], cfg)
        assert r is not None
        rows.append(r)
    df = pd.DataFrame(rows)

    # VPD 平均が年で単調増加
    assert df["VPD_mean"].is_monotonic_increasing
    # 目玉: 乾燥(高VPD)で I(Rg;GEP) が下がる → 負相関
    r = cr._spearman(df["I_Rg_GEP"].to_numpy(), df["VPD_mean"].to_numpy())
    assert r < -0.5
