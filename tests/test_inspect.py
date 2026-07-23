"""inspect_site の列マッピング検出と year-scan の検証 (合成 CSV)。"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from japanflux_pn.config import AnalysisConfig, RK_VARS
from japanflux_pn.sites import SiteSpec, DEFAULT_VAR_MAP
from japanflux_pn import inspect_site as insp


def _write_fluxnet_csv(path, cols: list[str], year: int, seed: int = 0,
                       missing_frac: float = 0.1):
    """FLUXNET2015 互換の HH CSV を 7/1..8/6 分だけ生成 (7月+前方窓)。"""
    idx = pd.date_range(f"{year}-07-01", f"{year}-08-06 23:30", freq="30min")
    rng = np.random.default_rng(seed)
    df = pd.DataFrame(index=idx)
    df["TIMESTAMP_START"] = idx.strftime("%Y%m%d%H%M")
    for c in cols:
        x = rng.random(len(idx))
        mask = rng.random(len(idx)) < missing_frac
        x[mask] = -9999.0
        df[c] = x
    df.to_csv(path, index=False)


@pytest.fixture
def good_site(tmp_path):
    """全 11 変数が標準名で揃った合成サイト。"""
    d = tmp_path / "JP-TEST"
    d.mkdir()
    cols = [DEFAULT_VAR_MAP[v] for v in RK_VARS]
    for yr in (2003, 2004):
        _write_fluxnet_csv(d / f"JP-TEST_COREVARS_HH_{yr}-{yr}.csv", cols, yr, seed=yr)
    return SiteSpec(code="JP-TEST", data_dir=str(d))


def test_check_mapping_all_present(good_site):
    header = insp._read_header(insp.find_corevars_files(good_site)[0])
    present, missing = insp.check_mapping(header, good_site)
    assert len(present) == 11
    assert missing == {}


def test_year_scan_counts(good_site):
    cfg = AnalysisConfig()
    tbl = insp.year_scan(good_site, cfg, years=[2003, 2004], months=[7])
    assert set(tbl["year"]) == {2003, 2004}
    # 7月 (~1488 格子) から前方窓損失 + ~10% 欠測 × 11 変数の listwise
    for _, r in tbl.iterrows():
        assert 0 < r["n_points"] < r["n_grid"]
        assert all(0.0 <= r[c] <= 1.0 for c in tbl.columns if c.startswith("cov_"))


def test_missing_column_suggests_candidate(tmp_path):
    """TS_F_MDS_1 の代わりに TS_F_MDS_2 → Ts の候補として提示される。"""
    d = tmp_path / "JP-ALT"
    d.mkdir()
    cols = [DEFAULT_VAR_MAP[v] for v in RK_VARS if v != "Ts"] + ["TS_F_MDS_2"]
    _write_fluxnet_csv(d / "JP-ALT_COREVARS_HH_2010-2010.csv", cols, 2010)
    site = SiteSpec(code="JP-ALT", data_dir=str(d))
    header = insp._read_header(insp.find_corevars_files(site)[0])
    present, missing = insp.check_mapping(header, site)
    assert "Ts" in missing
    assert "TS_F_MDS_2" in missing["Ts"]["candidates"]

    # override を当てれば解決する
    site2 = SiteSpec(code="JP-ALT", data_dir=str(d),
                     var_overrides={"Ts": "TS_F_MDS_2"})
    _, missing2 = insp.check_mapping(insp._read_header(
        insp.find_corevars_files(site2)[0]), site2)
    assert missing2 == {}


def test_report_smoke(good_site, capsys, monkeypatch):
    """report() が例外なく走り、健全年サジェストまで出力する。"""
    from japanflux_pn import sites as sites_mod
    monkeypatch.setitem(sites_mod.SITES, "JP-TEST", good_site)
    insp.report("JP-TEST", months=[7], years=[2003, 2004])
    out = capsys.readouterr().out
    assert "mapping" in out and "year-scan" in out
    assert "11/11 variables mapped OK" in out
