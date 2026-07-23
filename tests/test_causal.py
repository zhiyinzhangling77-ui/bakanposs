"""causal_network の graph 解析ロジック検証 (tigramite 不要部分)。

PCMCI 本体 (tigramite) が入っていれば統合スモークも走る。無ければ skip。
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from japanflux_pn.config import AnalysisConfig, RK_VARS
from japanflux_pn import causal_network as cn


def _mock_results(n, tau_max):
    """tigramite 風の graph / val_matrix を作る。"""
    graph = np.full((n, n, tau_max + 1), "", dtype="<U3")
    val = np.zeros((n, n, tau_max + 1), dtype=float)
    return graph, val


def test_extract_links_parses_directed_and_contemp():
    cfg = AnalysisConfig()
    n = len(RK_VARS)
    graph, val = _mock_results(n, cfg.lag_max)
    i_rg, i_ta, i_vpd, i_gle = (RK_VARS.index(v) for v in ("Rg", "Ta", "VPD", "gLE"))

    # Rg → Ta を lag1 で有向、強度 0.5
    graph[i_rg, i_ta, 1] = "-->"
    val[i_rg, i_ta, 1] = 0.5
    # VPD → gLE を lag2 で有向、強度 -0.3 (符号あり)
    graph[i_vpd, i_gle, 2] = "-->"
    val[i_vpd, i_gle, 2] = -0.3
    # 同時・向き未確定 Rg o-o VPD
    graph[i_rg, i_vpd, 0] = "o-o"
    graph[i_vpd, i_rg, 0] = "o-o"
    val[i_rg, i_vpd, 0] = 0.4

    links = cn.extract_links({"graph": graph, "val_matrix": val}, RK_VARS, cfg)
    directed = links[links["kind"] == "directed"]
    contemp = links[links["kind"] == "contemp_undirected"]

    assert len(directed) == 2
    assert len(contemp) == 1
    # 強度絶対値で降順 → Rg→Ta (0.5) が先頭
    top = links.iloc[0]
    assert top["src"] == "Rg" and top["dst"] == "Ta"
    assert top["lag_h"] == pytest.approx(0.5)     # lag1 = 0.5h
    vg = directed[(directed["src"] == "VPD") & (directed["dst"] == "gLE")].iloc[0]
    assert vg["lag_h"] == pytest.approx(1.0)       # lag2 = 1.0h
    assert vg["strength"] == pytest.approx(-0.3)


def test_extract_links_empty_graph():
    cfg = AnalysisConfig()
    graph, val = _mock_results(len(RK_VARS), cfg.lag_max)
    links = cn.extract_links({"graph": graph, "val_matrix": val}, RK_VARS, cfg)
    assert links.empty


def test_run_pcmci_requires_tigramite_or_full_coverage():
    """tigramite が無ければ ImportError、欠測ありなら ValueError を出す設計の確認。"""
    pytest.importorskip("tigramite")
    # tigramite があるなら、欠測ありフレームで ValueError になることを確認
    from japanflux_pn.preprocess import PreprocessResult
    cfg = AnalysisConfig()
    grid = pd.date_range("2020-07-01", periods=500, freq=pd.Timedelta(minutes=30))
    frame = pd.DataFrame({v: np.random.default_rng(0).random(500) for v in RK_VARS},
                         index=grid)[RK_VARS]
    valid = pd.Series(True, index=grid)
    valid.iloc[10] = False
    pre = PreprocessResult(anomaly=frame, valid=valid, site="S", year=2020,
                           month=7, config=cfg, months=[7])
    with pytest.raises(ValueError):
        cn.run_pcmci(pre)
