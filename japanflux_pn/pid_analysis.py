"""部分情報分解 (PID) — 同期を冗長/固有/相乗に分解 (Goodwell & Kumar 2017 系)。

条件付け診断は「Rg で説明できる同期は何割か (drop%)」を測った。PID はさらに厳密に、
目標 Y と 2 源 {Rg, X} が Y について持つ情報 I(Y; Rg, X) を分解する:

    R      … Rg と X が共有する冗長情報 (= 共通駆動の定量化)
    U_Rg   … Rg 固有 (X には無い)
    U_X    … X 固有 (Rg では説明できない = 真の直接結合)
    S      … Rg と X の組でしか得られない相乗情報

「フラックス相互同期は Rg との冗長 (R が支配) / 熱結合は X 固有 (U_X が残る)」を定量化。

    python -m japanflux_pn.pid_analysis --site JP-Tak --year 2003 --month 7 8
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig, RK_VARS
from .preprocess import load_corevars_hh, PreprocessResult
from . import information_theory as it


@dataclass
class PIDResult:
    driver: str
    table: pd.DataFrame     # target, source, I2, R, U_source, U_driver, S, redundancy_frac
    config: AnalysisConfig

    def save(self, outdir) -> None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        self.table.to_csv(outdir / f"pid_with_{self.driver}.csv", index=False)


def pid_with_driver(pre: PreprocessResult, driver: str = "Rg") -> PIDResult:
    """各 (目標 Y, 源 X) について I(Y; driver, X) を PID 分解する (X,Y ≠ driver)。"""
    cfg = pre.config
    m = cfg.n_bins
    vf = pre.valid_frame
    idx = {v: it.digitize_series(vf[v].to_numpy(dtype=float), m) for v in RK_VARS}
    dz = idx[driver]
    logm = cfg.log_m

    rows = []
    others = [v for v in RK_VARS if v != driver]
    for y in others:
        for x in others:
            if x == y:
                continue
            p = it.pid_williams_beer(idx[y], dz, idx[x], m)  # 源1=driver, 源2=X
            i2 = p["I2"]                                      # I(Y; X)
            red_frac = p["R"] / i2 if i2 > 1e-12 else np.nan
            rows.append({
                "target": y, "source": x,
                "I_YX_pct": 100.0 * i2 / logm,               # 素の I(Y;X)
                "R_pct": 100.0 * p["R"] / logm,              # Rg と冗長
                "U_source_pct": 100.0 * p["U2"] / logm,      # X 固有
                "U_driver_pct": 100.0 * p["U1"] / logm,      # Rg 固有
                "S_pct": 100.0 * p["S"] / logm,              # 相乗
                "redundancy_frac": red_frac,                 # R / I(Y;X)
            })
    table = pd.DataFrame(rows)
    return PIDResult(driver=driver, table=table, config=cfg)


def calibrate(pre: PreprocessResult, driver: str = "Rg") -> pd.DataFrame:
    """相乗モードの頑健性較正: 測度不変な相互作用情報 II（plugin と Miller-Madow）と
    I_min / MMI の S を並べる。II<0 が正味相乗（測度非依存）。MM で 3D バイアスを除く。"""
    cfg = pre.config
    m = cfg.n_bins
    vf = pre.valid_frame
    idx = {v: it.digitize_series(vf[v].to_numpy(dtype=float), m) for v in RK_VARS}
    dz = idx[driver]
    logm = cfg.log_m

    rows = []
    others = [v for v in RK_VARS if v != driver]
    for y in others:
        for x in others:
            if x == y:
                continue
            wb = it.pid_williams_beer(idx[y], dz, idx[x], m)
            mmi = it.pid_mmi(idx[y], dz, idx[x], m)
            ii_plug = it.interaction_information_indices(idx[y], dz, idx[x], m, False)
            ii_mm = it.interaction_information_indices(idx[y], dz, idx[x], m, True)
            rows.append({
                "target": y, "source": x,
                "I_YX_pct": 100.0 * wb["I2"] / logm,
                "S_imin_pct": 100.0 * wb["S"] / logm,
                "S_mmi_pct": 100.0 * mmi["S"] / logm,
                "II_plugin_pct": 100.0 * ii_plug / logm,   # >0 冗長, <0 相乗
                "II_mm_pct": 100.0 * ii_mm / logm,          # バイアス補正版
            })
    return pd.DataFrame(rows)


def report_calibration(site: str, year: int, months: list[int], driver: str = "Rg",
                       config: AnalysisConfig | None = None) -> pd.DataFrame:
    config = config or AnalysisConfig()
    pre = load_corevars_hh(site, year, months, config)
    print(f"[preprocess] {site} {year}-{pre.month_label}: n_points={pre.n_points}")
    tbl = calibrate(pre, driver)
    sig = tbl[tbl["I_YX_pct"] >= 3.0].copy()

    # I_min で相乗が大きかったペアが、測度不変 II とバイアス補正でも相乗のままか
    syn_imin = sig.sort_values("S_imin_pct", ascending=False).head(12)
    print(f"\n=== 相乗モードの較正 (II<0=正味相乗, 測度非依存) ===")
    print(f"  {'Y←X':<12} {'I(Y;X)':>7} {'S_Imin':>7} {'S_MMI':>7} "
          f"{'II_plug':>8} {'II_mm':>7}  判定")
    for _, r in syn_imin.iterrows():
        lab = f"{r['target']}←{r['source']}"
        verdict = "相乗(頑健)" if r["II_mm_pct"] < -0.5 else (
                  "相乗(弱/消)" if r["II_mm_pct"] < 0.5 else "→冗長に反転")
        print(f"  {lab:<12} {r['I_YX_pct']:6.1f} {r['S_imin_pct']:6.1f} "
              f"{r['S_mmi_pct']:6.1f} {r['II_plugin_pct']:7.1f} "
              f"{r['II_mm_pct']:6.1f}  {verdict}")

    # フラックス冗長ペアが II で正味冗長のままか (対照)
    red = sig.sort_values("II_mm_pct", ascending=False).head(6)
    print(f"\n  [対照: 正味冗長 II_mm>0 上位] "
          + ", ".join(f"{r['target']}←{r['source']}({r['II_mm_pct']:.1f})"
                      for _, r in red.iterrows()))
    n_syn_mm = int((sig["II_mm_pct"] < -0.5).sum())
    n_red_mm = int((sig["II_mm_pct"] > 0.5).sum())
    print(f"\n  MM 補正後: 正味相乗 {n_syn_mm} / 正味冗長 {n_red_mm} "
          f"/ ほぼ0 {len(sig)-n_syn_mm-n_red_mm} (計 {len(sig)} ペア)")
    return tbl


def report(site: str, year: int, months: list[int], driver: str = "Rg",
           config: AnalysisConfig | None = None,
           outroot: str | Path | None = None) -> PIDResult:
    config = config or AnalysisConfig()
    pre = load_corevars_hh(site, year, months, config)
    print(f"[preprocess] {site} {year}-{pre.month_label}: n_points={pre.n_points}")
    res = pid_with_driver(pre, driver)
    tbl = res.table

    # I(Y;X) が大きい上位ペアで、その情報を R/U_X/S へどう割るか
    top = tbl.sort_values("I_YX_pct", ascending=False).head(20)
    print(f"\n=== PID: I(Y; X) を {driver} と冗長(R)/X固有(U_X)/相乗(S) に分解 ===")
    print(f"  (上位 20 ペア, % は /log m 正規化)")
    print(f"  {'Y←X':<12} {'I(Y;X)':>7} {'R':>6} {'U_X':>6} {'S':>6} {'R/I':>6}")
    for _, r in top.iterrows():
        lab = f"{r['target']}←{r['source']}"
        print(f"  {lab:<12} {r['I_YX_pct']:6.1f} {r['R_pct']:6.1f} "
              f"{r['U_source_pct']:6.1f} {r['S_pct']:6.1f} {r['redundancy_frac']:5.0%}")

    # 冗長支配 vs 固有支配で仕分け (I(Y;X) が意味を持つペアのみ)
    sig = tbl[tbl["I_YX_pct"] >= 3.0].copy()
    red_dom = sig[sig["redundancy_frac"] >= 0.6]
    uniq_dom = sig[sig["redundancy_frac"] < 0.4]
    print(f"\n  [Rg 冗長支配] R/I≥60% ({len(red_dom)}): "
          + ", ".join(f"{r['target']}←{r['source']}" for _, r in
                      red_dom.sort_values('redundancy_frac', ascending=False)
                      .head(12).iterrows()))
    print(f"  [X 固有支配] R/I<40% ({len(uniq_dom)}): "
          + ", ".join(f"{r['target']}←{r['source']}" for _, r in
                      uniq_dom.sort_values('redundancy_frac').head(12).iterrows()))
    syn = sig.sort_values("S_pct", ascending=False).head(5)
    print(f"  [相乗 上位] "
          + ", ".join(f"{r['target']}←{r['source']}({r['S_pct']:.1f})"
                      for _, r in syn.iterrows()))

    if outroot is not None:
        outdir = Path(outroot) / f"{site}_{year}{pre.month_label}_pid{driver}"
        res.save(outdir)
        print(f"\n[output] {outdir}")
    return res


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="部分情報分解 (PID) with 共通駆動")
    p.add_argument("--site", required=True)
    p.add_argument("--year", type=int, required=True)
    p.add_argument("--month", type=int, nargs="+", default=[7])
    p.add_argument("--driver", default="Rg")
    p.add_argument("--bins", type=int, default=None)
    p.add_argument("--qc-max", type=int, default=None)
    p.add_argument("--calibrate", action="store_true",
                   help="相乗モードの較正 (測度不変な II + Miller-Madow 補正)")
    p.add_argument("--outroot", default=None)
    args = p.parse_args(argv)
    kw = {}
    if args.bins is not None:
        kw["n_bins"] = args.bins
    if args.qc_max is not None:
        kw["qc_max"] = args.qc_max
    config = AnalysisConfig(**kw)
    if args.calibrate:
        report_calibration(args.site, args.year, args.month, args.driver, config)
    else:
        report(args.site, args.year, args.month, args.driver, config, args.outroot)


if __name__ == "__main__":
    main()
