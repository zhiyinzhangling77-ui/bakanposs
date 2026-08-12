"""感度スイープ (Phase 1): ビン数 m を振って fig2b の結論が不変かを確認する。

fig2b の主張は「共通駆動 Rg で条件付けると、エネルギー・炭素フラックス結合は崩れ、
気温・土壌・呼吸の結合は残る」。この定性的な二分が **ビン数 m の選び方に依存しない**
ことを示せれば、「m=11 を都合よく選んだ」という批判を封じられる。

各 m について I(X;Y) と I(X;Y|Rg) を計算し、2 群の平均 drop% を集計する:
  - flux 群 (崩れるはず): gH-GEP, gLE-GEP, gH-gLE, NEE-GEP
  - thermal 群 (残るはず): Ta-Ts, Ta-GER, Ts-GER, th-GER
複数年を指定すると年で平均し、頑健性を上げる。

    python -m japanflux_pn.run_sensitivity --site JP-Tak --years 2003 2004 2005 --month 7 8

夜間・切断耐性:
    nohup python -m japanflux_pn.run_sensitivity --site JP-Tak \\
        --years 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 \\
        --bins 7 9 11 13 15 --month 7 8 \\
        --outroot ~/bakanposs/japanflux_pn/outputs_sensitivity \\
        > sens_Tak.log 2>&1 &
"""

from __future__ import annotations

import argparse
import time
from dataclasses import replace
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig, RK_VARS
from .preprocess import load_corevars_hh
from .condition_driver import condition_on_driver

# fig2b と同じ 2 群 (RK 変数名)
FLUX_PAIRS = [("gH", "GEP"), ("gLE", "GEP"), ("gH", "gLE"), ("NEE", "GEP")]
THERMAL_PAIRS = [("Ta", "Ts"), ("Ta", "GER"), ("Ts", "GER"), ("th", "GER")]


def _drop_pct(res, a: str, b: str) -> float:
    """Rg 条件付けによる MI の減少率 [%] (fig2b と同じ定義)。"""
    i, c = res.mi.loc[a, b], res.cmi.loc[a, b]
    return 100.0 * (1 - c / i) if (i is not None and i > 0) else np.nan


def _group_mean_drop(res, pairs) -> float:
    vals = [_drop_pct(res, a, b) for a, b in pairs]
    vals = [v for v in vals if v == v]  # drop NaN
    return float(np.mean(vals)) if vals else np.nan


def sweep(site: str, years: list[int], months: list[int], bins: list[int],
          driver: str = "Rg", outroot: str | Path | None = None) -> pd.DataFrame:
    base = AnalysisConfig()
    t0 = time.time()
    n_jobs = len(bins) * len(years)
    done = 0
    print(f"===== 感度スイープ {site} | bins={bins} | years={years} | months={months} =====")
    print(f"  開始 {time.strftime('%Y-%m-%d %H:%M:%S')}  (全 {n_jobs} ジョブ)\n", flush=True)

    records: list[dict] = []
    for m in bins:
        cfg = replace(base, n_bins=m)
        per_year_flux: list[float] = []
        per_year_therm: list[float] = []
        for y in years:
            done += 1
            try:
                pre = load_corevars_hh(site, y, months, cfg)
                res = condition_on_driver(pre, driver)
                fd = _group_mean_drop(res, FLUX_PAIRS)
                td = _group_mean_drop(res, THERMAL_PAIRS)
                per_year_flux.append(fd)
                per_year_therm.append(td)
                elapsed = time.time() - t0
                eta = elapsed / done * (n_jobs - done)
                print(f"  [m={m:2d} y={y}] n={pre.n_points:5d} "
                      f"flux群 drop={fd:5.1f}%  thermal群 drop={td:5.1f}%  "
                      f"({done}/{n_jobs}, ETA {eta/60:.1f}min)", flush=True)
            except Exception as e:  # 1 年の欠測等で落ちても続行
                print(f"  [m={m:2d} y={y}] SKIP: {type(e).__name__}: {e}", flush=True)
        fmean = float(np.nanmean(per_year_flux)) if per_year_flux else np.nan
        tmean = float(np.nanmean(per_year_therm)) if per_year_therm else np.nan
        records.append({
            "n_bins": m,
            "flux_mean_drop_pct": fmean,
            "thermal_mean_drop_pct": tmean,
            "separation": tmean - fmean,  # 正なら「flux崩れ・thermal残る」の二分が成立
            "n_years_used": len(per_year_flux),
        })
        print(f"  => m={m}: flux群平均 drop={fmean:.1f}%  thermal群平均 drop={tmean:.1f}%\n",
              flush=True)

    df = pd.DataFrame.from_records(records)
    print("===== まとめ =====")
    print(df.to_string(index=False))
    # 定性的結論: 全 m で flux>thermal の drop (=崩れる) かつ分離>0 か
    ok = bool((df["flux_mean_drop_pct"] > df["thermal_mean_drop_pct"]).all()) \
        and bool((df["separation"] > 0).all())
    verdict = ("✅ 全ビンで結論不変: flux群は崩れ、thermal群は残る（m非依存）"
               if ok else
               "⚠ 一部ビンで二分が崩れる: 原因を確認（データ量/推定バイアス）")
    print("\n" + verdict, flush=True)

    if outroot:
        outdir = Path(outroot)
        outdir.mkdir(parents=True, exist_ok=True)
        path = outdir / f"sensitivity_bins_{site}.csv"
        df.to_csv(path, index=False)
        print(f"[saved] {path}", flush=True)
    return df


def main() -> None:
    p = argparse.ArgumentParser(description="fig2b のビン数感度スイープ (Phase 1)")
    p.add_argument("--site", required=True)
    p.add_argument("--years", type=int, nargs="+", required=True)
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--bins", type=int, nargs="+", default=[7, 9, 11, 13, 15])
    p.add_argument("--driver", default="Rg")
    p.add_argument("--outroot", default=None)
    a = p.parse_args()
    sweep(a.site, a.years, a.month, a.bins, a.driver, a.outroot)


if __name__ == "__main__":
    main()
