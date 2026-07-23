"""条件付き相互情報 I(X;Y|Z) と共通駆動分離の検証。"""

from __future__ import annotations

import numpy as np
import pandas as pd

from japanflux_pn.config import AnalysisConfig, RK_VARS
from japanflux_pn import information_theory as it
from japanflux_pn import condition_driver as cd
from japanflux_pn.preprocess import PreprocessResult


def test_cmi_collapses_for_common_driver():
    """X = Z+ノイズ, Y = Z+ノイズ (共通駆動)。I(X;Y)>0 だが I(X;Y|Z)≈0。"""
    m = 11
    rng = np.random.default_rng(0)
    n = 4000
    z = rng.normal(size=n)
    x = z + 0.3 * rng.normal(size=n)
    y = z + 0.3 * rng.normal(size=n)
    xi = it.digitize_series(x, m)
    yi = it.digitize_series(y, m)
    zi = it.digitize_series(z, m)

    mi = it.mutual_information_indices(xi, yi, m)
    cmi = it.conditional_mutual_information_indices(xi, yi, [zi], m)
    assert mi > 0.3                      # 共通駆動で強く相関
    # Z を与えると大半が消える (ビン化残差ぶんだけ残るので 0 ではない)
    assert cmi < 0.2 * mi                # 80%以上の減少


def test_cmi_survives_for_direct_chain():
    """X→Y 直接結合 + Z は無関係 → I(X;Y|Z) は有意に残る。"""
    m = 11
    rng = np.random.default_rng(2)
    n = 4000
    x = rng.normal(size=n)
    y = 0.8 * x + 0.4 * rng.normal(size=n)   # 直接 X→Y
    z = rng.normal(size=n)                    # 無関係な "駆動"
    xi, yi, zi = (it.digitize_series(v, m) for v in (x, y, z))

    cmi = it.conditional_mutual_information_indices(xi, yi, [zi], m)
    stats = it.surrogate_cmi_stats(xi, yi, [zi], m, 40, 2.36,
                                   np.random.default_rng(3))
    assert cmi > stats["threshold"]      # 無関係な Z を条件付けても直接依存は残る


def test_cmi_equals_mi_when_z_independent():
    """独立な Z で条件付けると I(X;Y|Z) ≈ I(X;Y) (Z が情報を持たない)。"""
    m = 8
    rng = np.random.default_rng(4)
    n = 5000
    x = rng.normal(size=n)
    y = 0.7 * x + 0.5 * rng.normal(size=n)
    z = rng.normal(size=n)
    xi, yi, zi = (it.digitize_series(v, m) for v in (x, y, z))
    mi = it.mutual_information_indices(xi, yi, m)
    cmi = it.conditional_mutual_information_indices(xi, yi, [zi], m)
    # Z は無関係だが次元が増えるぶん推定は上振れしうる。桁として一致することを確認
    assert cmi == mi or abs(cmi - mi) < 0.6 * mi + 0.05


def _synth_pre(n=3000):
    """Rg が Ta と VPD を共通駆動、gLE は VPD に直接依存する合成フレーム。"""
    cfg = AnalysisConfig(n_surrogates=30, seed=0)
    rng = np.random.default_rng(7)
    rg = rng.normal(size=n)
    data = {v: rng.normal(size=n) for v in RK_VARS}
    data["Rg"] = rg
    data["Ta"] = rg + 0.3 * rng.normal(size=n)     # ほぼ Rg 共通駆動
    data["VPD"] = 0.5 * rg + rng.normal(size=n)    # Rg 一部 + 大きな独立成分
    data["gLE"] = data["VPD"] + 0.3 * rng.normal(size=n)  # VPD に直接連動
    grid = pd.date_range("2020-07-01", periods=n, freq=pd.Timedelta(minutes=30))
    frame = pd.DataFrame(data, index=grid)[RK_VARS]
    valid = pd.Series(True, index=grid)
    return PreprocessResult(anomaly=frame, valid=valid, site="S",
                            year=2020, month=7, config=cfg, months=[7])


def test_condition_on_driver_separates_common_and_direct():
    pre = _synth_pre()
    res = cd.condition_on_driver(pre, driver="Rg")

    def drop(a, b):
        return 1 - res.cmi.loc[a, b] / res.mi.loc[a, b]

    # 両ペアとも素の MI は有意
    assert bool(res.mi_sig.loc["Ta", "VPD"]) and bool(res.mi_sig.loc["VPD", "gLE"])
    # Ta-VPD (共通駆動) は Rg 条件付けで大きく減る / VPD-gLE (直接) はよく残る
    assert drop("Ta", "VPD") > 0.6
    assert drop("VPD", "gLE") < 0.5
    assert drop("Ta", "VPD") > drop("VPD", "gLE")
