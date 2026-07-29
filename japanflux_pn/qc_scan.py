"""QC 被覆スキャン: 実測のみ (gap-fill 除外) に絞ったとき、11 変数 listwise で
何点残るかを測る。QC 制限再解析が feasible か（n≥500 か）を先に判断するための下調べ。

FLUXNET2015/JapanFlux2024 の QC 慣習: gap-fill 変数 (`_F`, `_MDS`) は QC=0 実測、
1/2/3 が補完品質。閾値以下だけ残す。派生炭素 (RECO/GPP) は自前 QC を持たないため
NEE の QC を代理に使う。

    python -m japanflux_pn.qc_scan --site JP-Tak --year 2003 --month 7 8
    python -m japanflux_pn.qc_scan --site JP-Tak --year 2003 --month 7 8 --qc-max 1

`n_measured` が 500 を大きく割るなら QC=0 は非現実的 → 閾値を緩める(1)か、
「実測+良質補完」で妥協する、等の判断材料になる。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig, RK_VARS
from .sites import SiteSpec, get_site, resolve_qc_columns
from .preprocess import find_corevars_files, _regular_grid


def _load_values_and_qc(files: list[Path], site: SiteSpec, config: AnalysisConfig):
    """値列と QC 列を読み、(values[RK], qc[RK or None]) を返す。"""
    vmap = site.var_map()
    header = list(pd.read_csv(files[0], nrows=0).columns)
    qcmap = resolve_qc_columns(header, site)
    val_cols = [vmap[v] for v in RK_VARS]
    qc_cols = sorted({c for c in qcmap.values() if c})
    want = set(["TIMESTAMP_START"] + val_cols + qc_cols)

    frames = []
    for f in files:
        df = pd.read_csv(f, usecols=lambda c: c in want)
        ts = pd.to_datetime(df["TIMESTAMP_START"].astype("int64").astype(str),
                            format="%Y%m%d%H%M")
        df = df.drop(columns=["TIMESTAMP_START"])
        df.index = ts
        frames.append(df)
    allc = pd.concat(frames)
    allc = allc[~allc.index.duplicated(keep="first")].sort_index()
    allc = allc.replace(config.na_sentinel, np.nan)
    return allc, vmap, qcmap


def qc_scan(site: SiteSpec, config: AnalysisConfig, year: int, months: list[int],
            qc_max: int) -> dict:
    """対象期間で、各変数の実測率と 11 変数 listwise の実測点数を集計。"""
    files = find_corevars_files(site)
    allc, vmap, qcmap = _load_values_and_qc(files, site, config)

    months = sorted(months)
    start = pd.Timestamp(year=year, month=months[0], day=1)
    end = pd.Timestamp(year=year, month=months[-1], day=1) + pd.offsets.MonthBegin(1)
    step = pd.Timedelta(minutes=24 * 60 // config.steps_per_day)
    grid = _regular_grid(start, end - step, config.steps_per_day)
    sub = allc.reindex(grid)

    per_var = {}
    measured_mask = pd.DataFrame(index=grid)
    for v in RK_VARS:
        val = sub[vmap[v]]
        has_val = val.notna()
        qc = qcmap[v]
        if qc is None:
            meas = has_val                       # QC 無し → 値があれば実測扱い
            note = "no-QC"
        else:
            meas = has_val & (sub[qc] <= qc_max)  # 値あり かつ QC≤閾値
            note = qc
        measured_mask[v] = meas
        per_var[v] = {"measured_frac": float(meas.mean()), "qc_col": note}

    listwise = measured_mask.all(axis=1)
    n_grid = len(grid)
    limiting = min(RK_VARS, key=lambda v: per_var[v]["measured_frac"])
    return {"n_grid": n_grid, "n_measured": int(listwise.sum()),
            "per_var": per_var, "limiting": limiting, "qc_max": qc_max}


def report(site_code: str, year: int, months: list[int], qc_max: int,
           config: AnalysisConfig | None = None) -> None:
    config = config or AnalysisConfig()
    site = get_site(site_code)
    print(f"### QC 被覆スキャン {site_code} {year} months={months} (QC≤{qc_max}) ###\n")
    r = qc_scan(site, config, year, months, qc_max)

    print(f"  {'var':>4} {'measured%':>10}  QC列")
    for v in RK_VARS:
        pv = r["per_var"][v]
        print(f"  {v:>4} {pv['measured_frac']:9.0%}   {pv['qc_col']}")

    n, ng = r["n_measured"], r["n_grid"]
    print(f"\n  listwise 全11変数実測: n_measured = {n}/{ng} ({n/ng:.0%})")
    print(f"  律速変数 (最低実測率): {r['limiting']} "
          f"({r['per_var'][r['limiting']]['measured_frac']:.0%})")
    if n >= 500:
        print(f"\n  → feasible (n={n} ≥ 500)。QC≤{qc_max} で再解析可能。")
    else:
        print(f"\n  → 点数不足 (n={n} < 500)。QC 閾値を緩める(--qc-max {qc_max+1})か、"
              f"月をプールする(--month 7 8)、または QC 制限は断念して gap-fill 版を"
              f"本文・QC 版は感度解析、等を検討。")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="QC 被覆スキャン (実測点数の feasibility)")
    p.add_argument("--site", required=True)
    p.add_argument("--year", type=int, required=True)
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--qc-max", type=int, default=0,
                   help="残す QC 上限 (0=実測のみ, 1=実測+良質補完)")
    args = p.parse_args(argv)
    report(args.site, args.year, args.month, args.qc_max)


if __name__ == "__main__":
    main()
