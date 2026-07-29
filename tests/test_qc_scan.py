"""qc_scan の QC 列解決と listwise 集計の検証 (合成 CSV)。"""

from __future__ import annotations

import numpy as np
import pandas as pd

from japanflux_pn.config import AnalysisConfig, RK_VARS
from japanflux_pn.sites import SiteSpec, DEFAULT_VAR_MAP
from japanflux_pn import qc_scan as qs


def _write_csv_with_qc(path, year, seed=0):
    """値列 + 一部の QC 列を持つ合成 HH CSV。"""
    idx = pd.date_range(f"{year}-07-01", f"{year}-07-31 23:30", freq="30min")
    rng = np.random.default_rng(seed)
    df = pd.DataFrame(index=idx)
    df["TIMESTAMP_START"] = idx.strftime("%Y%m%d%H%M")
    for v in RK_VARS:
        col = DEFAULT_VAR_MAP[v]
        df[col] = rng.random(len(idx))
        # NEE と TS には QC 列を付ける (半分が実測=0, 半分が補完=2)
        if v in ("NEE", "Ts", "gLE"):
            qc = np.where(rng.random(len(idx)) < 0.5, 0, 2)
            df[col + "_QC"] = qc
    df.to_csv(path, index=False)


def _site(tmp_path):
    d = tmp_path / "JP-QC"
    d.mkdir()
    _write_csv_with_qc(d / "JP-QC_COREVARS_HH_2003-2003.csv", 2003)
    return SiteSpec(code="JP-QC", data_dir=str(d))


def test_resolve_qc_columns_and_derived_inherit(tmp_path):
    site = _site(tmp_path)
    header = list(pd.read_csv(qs.find_corevars_files(site)[0], nrows=0).columns)
    qcmap = qs.resolve_qc_columns(header, site)
    assert qcmap["NEE"] == "NEE_vUT_QC"
    assert qcmap["Ts"] == "TS_F_MDS_QC"
    # 派生炭素 GER/GEP は自前 QC が無いので NEE の QC を継承
    assert qcmap["GER"] == "NEE_vUT_QC"
    assert qcmap["GEP"] == "NEE_vUT_QC"
    # QC 列が無く派生でもない変数は None (常に実測扱い)
    assert qcmap["Rg"] is None


def test_read_corevars_raw_applies_qc_mask(tmp_path):
    """config.qc_max=0 で QC>0 の値が NaN 化される（実測のみ残る）。"""
    from japanflux_pn import preprocess as pp

    site = _site(tmp_path)
    f = pp.find_corevars_files(site)[0]

    # qc_max=None: gap-fill 込み → Ts は全点有限
    raw_full = pp.read_corevars_raw(f, site, AnalysisConfig())
    assert raw_full["Ts"].notna().all()

    # qc_max=0: Ts_QC>0 の点が NaN、実測(QC=0)のみ残る (~50%)
    raw_qc = pp.read_corevars_raw(f, site, AnalysisConfig(qc_max=0))
    frac = raw_qc["Ts"].notna().mean()
    assert 0.3 < frac < 0.7
    # QC 列が無い Rg は qc_max=0 でも全点残る
    assert raw_qc["Rg"].notna().all()


def test_qc_scan_listwise_drops_with_strict_qc(tmp_path):
    site = _site(tmp_path)
    cfg = AnalysisConfig()
    strict = qs.qc_scan(site, cfg, 2003, [7], qc_max=0)
    loose = qs.qc_scan(site, cfg, 2003, [7], qc_max=2)
    # QC 列を持つ変数は実測率 ~50%、listwise は厳しくなる
    assert 0 < strict["n_measured"] < loose["n_measured"]
    assert loose["n_measured"] == loose["n_grid"]   # QC≤2 なら全部残る
    assert strict["per_var"]["NEE"]["measured_frac"] < 0.7
    assert strict["limiting"] in ("NEE", "Ts", "gLE", "GER", "GEP")
