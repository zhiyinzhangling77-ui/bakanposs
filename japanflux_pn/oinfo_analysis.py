"""高次情報構造 — O-information (Rosas et al. 2019) で系レベルの冗長支配 vs 相乗支配。

PID は 2 源→目標の分解。O-information はサブシステム全体を 1 数値で特徴づける:
Ω>0 冗長支配（共通駆動で情報重複）、Ω<0 相乗支配（情報が組にしか宿らない創発構造）。
物理的に意味のある小集合（呼吸制御・蒸発制御・光合成制御…）を 3 生態系で測り、
「土壌・呼吸系は相乗支配か、水田では湛水でそれが壊れるか」という高次の署名を問う。

次元の呪い対策: n=4 変数の同時分布は疎（m^4 セル）なので、粗ビン（既定 m=6）
＋ Miller-Madow 補正 ＋ シャッフルサロゲート（同じ疎性のヌルでバイアスを吸収）。
z=(Ω−μ_ss)/σ_ss の符号と大きさで冗長/相乗を判定する。

    python -m japanflux_pn.oinfo_analysis --site JP-Tak --year 2003 --month 7 8
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig, RK_VARS
from .preprocess import load_corevars_hh, PreprocessResult
from . import information_theory as it


# 高次同時分布の疎性を抑える粗ビン（PID/TE の m=11 とは別）
OINFO_BINS_DEFAULT = 6

# 物理的に意味のあるサブシステム（4 変数）
SUBSYSTEMS: dict[str, list[str]] = {
    "呼吸制御 Rg,Ta,θ,GER":    ["Rg", "Ta", "th", "GER"],
    "蒸発制御 Rg,VPD,θ,γLE":   ["Rg", "VPD", "th", "gLE"],
    "光合成制御 Rg,Ta,VPD,GEP": ["Rg", "Ta", "VPD", "GEP"],
    "炭素コア Rg,GEP,GER,NEE":  ["Rg", "GEP", "GER", "NEE"],
    "土壌熱 Ta,Ts,θ,GER":       ["Ta", "Ts", "th", "GER"],
    "エネルギー Rg,Ta,γH,γLE":  ["Rg", "Ta", "gH", "gLE"],
}


def o_information_subsystems(pre: PreprocessResult, obins: int = OINFO_BINS_DEFAULT
                             ) -> pd.DataFrame:
    """各サブシステムの O-information とシャッフルヌルからの z を返す。"""
    cfg = pre.config
    vf = pre.valid_frame
    idx = {v: it.digitize_series(vf[v].to_numpy(dtype=float), obins) for v in RK_VARS}
    rng = np.random.default_rng(cfg.seed)

    rows = []
    for name, vs in SUBSYSTEMS.items():
        cols = [idx[v] for v in vs]
        omega = it.o_information_indices(cols, obins, correct=True)
        tc = it.total_correlation_indices(cols, obins, correct=True)
        stats = it.surrogate_o_information_stats(
            cols, obins, cfg.n_surrogates, cfg.sig_c, rng, correct=True)
        z = (omega - stats["mu"]) / stats["sigma"] if stats["sigma"] > 0 else np.nan
        rows.append({
            "subsystem": name, "vars": ",".join(vs),
            "Omega": omega, "TC": tc,
            "null_mu": stats["mu"], "z": z,
        })
    return pd.DataFrame(rows)


def scan_years(site: str, obins: int = OINFO_BINS_DEFAULT,
               config: AnalysisConfig | None = None) -> pd.DataFrame:
    """健全全年で各サブシステムの O-info z を計算 (相乗崩壊の年々安定性の検証)。"""
    from .run_robustness import SITE_YEARS
    from .sites import get_site
    from .preprocess import (load_raw_all, slice_and_anomaly,
                             slice_span_and_anomaly, PreprocessResult)

    config = config or AnalysisConfig()
    years, months = SITE_YEARS[site]
    raw_all = load_raw_all(get_site(site), config)
    rows = []
    for y in years:
        if len(months) == 1:
            anom, valid = slice_and_anomaly(raw_all, y, months[0], config)
        else:
            anom, valid = slice_span_and_anomaly(raw_all, y, months, config)
        if int(valid.sum()) < 500:
            continue
        pre = PreprocessResult(anomaly=anom, valid=valid, site=site, year=y,
                               month=months[0], config=config, months=months)
        tbl = o_information_subsystems(pre, obins)
        for _, r in tbl.iterrows():
            rows.append({"year": y, "subsystem": r["subsystem"],
                         "Omega": r["Omega"], "z": r["z"]})
        print(f"  {site} {y}: done", flush=True)
    return pd.DataFrame(rows)


def report_years(site: str, obins: int = OINFO_BINS_DEFAULT,
                 config: AnalysisConfig | None = None,
                 outroot: str | Path | None = None) -> pd.DataFrame:
    print(f"===== {site} O-information 複数年安定性 (obins={obins}) =====")
    long = scan_years(site, obins, config)
    if long.empty:
        print("  (有効年なし)")
        return long
    n_years = long["year"].nunique()
    print(f"\n=== 各サブシステムの相乗/冗長の年々一貫性 (有効 {n_years} 年) ===")
    print(f"  {'subsystem':<22} {'mean z':>7} {'相乗年':>7} {'冗長年':>7}  傾向")
    summary = []
    for name, g in long.groupby("subsystem", sort=False):
        mz = float(g["z"].mean())
        n_syn = int((g["z"] <= -2.36).sum())
        n_red = int((g["z"] >= 2.36).sum())
        tot = len(g)
        if n_syn >= 0.6 * tot:
            trend = "★相乗が安定"
        elif n_red >= 0.6 * tot:
            trend = "冗長が安定"
        else:
            trend = "混在/不安定"
        summary.append({"subsystem": name, "mean_z": mz,
                        "n_synergy": n_syn, "n_redundant": n_red, "n_years": tot})
        print(f"  {name:<22} {mz:7.1f} {n_syn:4d}/{tot:<2d} {n_red:4d}/{tot:<2d}  {trend}")

    if outroot is not None:
        outdir = Path(outroot)
        outdir.mkdir(parents=True, exist_ok=True)
        long.to_csv(outdir / f"{site}_oinfo_multiyear.csv", index=False)
        pd.DataFrame(summary).to_csv(
            outdir / f"{site}_oinfo_multiyear_summary.csv", index=False)
        print(f"\n[output] {outdir}/{site}_oinfo_multiyear*.csv")
    return long


def report(site: str, year: int, months: list[int], obins: int = OINFO_BINS_DEFAULT,
           config: AnalysisConfig | None = None,
           outroot: str | Path | None = None) -> pd.DataFrame:
    config = config or AnalysisConfig()
    pre = load_corevars_hh(site, year, months, config)
    print(f"[preprocess] {site} {year}-{pre.month_label}: n_points={pre.n_points} "
          f"(O-info bins={obins})")
    tbl = o_information_subsystems(pre, obins)

    print(f"\n=== O-information: 系レベルの冗長 vs 相乗 (シャッフルヌル比) ===")
    print(f"  ※ 4変数ヒストは疎で Ω が負にバイアスするため、絶対符号でなく z で判定。")
    print(f"    z>0(有意)=独立より冗長, z<0(有意)=独立より相乗。|z|≥2.36 で有意。")
    print(f"  {'subsystem':<22} {'Omega':>8} {'null_mu':>8} {'z':>7}  判定")
    for _, r in tbl.iterrows():
        z = r["z"]
        if not np.isfinite(z) or abs(z) < 2.36:
            verdict = "有意でない"
        elif z > 0:
            verdict = "冗長的 (共通駆動)"
        else:
            verdict = "★相乗的 (創発)"
        name = r["subsystem"]
        print(f"  {name:<22} {r['Omega']:8.3f} {r['null_mu']:8.3f} "
              f"{z:7.1f}  {verdict}")

    if outroot is not None:
        outdir = Path(outroot)
        outdir.mkdir(parents=True, exist_ok=True)
        tbl.to_csv(outdir / f"{site}_{year}{pre.month_label}_oinfo.csv", index=False)
        print(f"\n[output] {outdir}/{site}_{year}{pre.month_label}_oinfo.csv")
    return tbl


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="O-information 高次情報構造")
    p.add_argument("--site", required=True)
    p.add_argument("--year", type=int, default=None,
                   help="対象年 (--multiyear 時は不要)")
    p.add_argument("--month", type=int, nargs="+", default=[7])
    p.add_argument("--obins", type=int, default=OINFO_BINS_DEFAULT,
                   help=f"O-info 用ビン数 (既定 {OINFO_BINS_DEFAULT}, 疎性対策で粗め)")
    p.add_argument("--multiyear", action="store_true",
                   help="健全全年で相乗/冗長の年々安定性を集計 (--year 無視)")
    p.add_argument("--qc-max", type=int, default=None)
    p.add_argument("--outroot", default=None)
    args = p.parse_args(argv)
    kw = {}
    if args.qc_max is not None:
        kw["qc_max"] = args.qc_max
    config = AnalysisConfig(**kw)
    if args.multiyear:
        report_years(args.site, args.obins, config, args.outroot)
    else:
        if args.year is None:
            p.error("--year は必須です (--multiyear なら不要)")
        report(args.site, args.year, args.month, args.obins, config, args.outroot)


if __name__ == "__main__":
    main()
