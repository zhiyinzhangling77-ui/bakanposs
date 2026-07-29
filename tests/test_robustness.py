"""run_robustness の集計ロジック検証 (tigramite 不要)。"""

from __future__ import annotations

import pandas as pd

from japanflux_pn import run_robustness as rb


def _links(pairs):
    """(src, dst, lag_h, strength) タプル列から directed リンク DataFrame を作る。"""
    rows = [{"src": s, "dst": d, "lag_h": lag, "strength": st, "kind": "directed"}
            for s, d, lag, st in pairs]
    return pd.DataFrame(rows)


def test_aggregate_counts_frequency_and_core():
    per_year = {
        2001: _links([("Rg", "gLE", 0.5, 0.2), ("GEP", "NEE", 0.5, 0.3),
                      ("Ta", "Rg", 5.0, 0.05)]),          # Ta→Rg は偽陽性
        2002: _links([("Rg", "gLE", 0.5, 0.25), ("GEP", "NEE", 0.5, 0.28)]),
        2003: _links([("Rg", "gLE", 1.0, 0.22), ("GEP", "NEE", 0.5, 0.31),
                      ("VPD", "gLE", 0.5, 0.1)]),         # VPD→gLE は 1 年のみ
        2004: None,                                        # 欠測年 → 除外
    }
    tbl, n_total = rb.aggregate(per_year)
    assert n_total == 3                                    # None は数えない

    def row(s, d):
        return tbl[(tbl["src"] == s) & (tbl["dst"] == d)].iloc[0]

    # Rg→gLE, GEP→NEE は 3 年全部 → 頻度 1.0 (コア)
    assert row("Rg", "gLE")["n_years"] == 3
    assert row("Rg", "gLE")["frequency"] == 1.0
    assert row("GEP", "NEE")["frequency"] == 1.0
    # VPD→gLE は 1 年のみ → 散発
    assert row("VPD", "gLE")["n_years"] == 1
    # Ta→Rg (偽陽性) は 1 年のみ
    assert row("Ta", "Rg")["n_years"] == 1
    # 頻度降順に並ぶ
    assert tbl.iloc[0]["frequency"] >= tbl.iloc[-1]["frequency"]


def test_aggregate_empty():
    tbl, n_total = rb.aggregate({2001: None, 2002: None})
    assert n_total == 0
    assert tbl.empty


def test_aggregate_ignores_nondirected():
    per_year = {
        2001: pd.DataFrame([
            {"src": "Rg", "dst": "gLE", "lag_h": 0.5, "strength": 0.2,
             "kind": "directed"},
            {"src": "Rg", "dst": "Ta", "lag_h": 0.0, "strength": 0.4,
             "kind": "contemp_undirected"},   # 無向は集計しない
        ]),
    }
    tbl, _ = rb.aggregate(per_year)
    assert len(tbl) == 1 and tbl.iloc[0]["dst"] == "gLE"
