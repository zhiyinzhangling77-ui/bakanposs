"""複数年ロバスト性: 各サイトを複数の健全年で PCMCI にかけ、因果リンクの年々一貫性を集計。

因果骨格 (Rg→γLE, GEP→NEE, Ta→Ts …) が単年の偶然でなく安定特性かを確認する。各年で
PCMCI+ の有向リンクを取り、(src→dst) ごとに「何年出現したか」を数える。コアリンクは
高頻度、偽陽性 (例: Rg への逆向き) は散発的になるはず。

breadth を稼ぐため既定は高速な parcorr、tau_max=6 (有向リンクは ≤1h に集中)。完全被覆
でない年 (センサ故障等) は自動でスキップする。

    python -m japanflux_pn.run_robustness --site JP-Tak
    python -m japanflux_pn.run_robustness --site JP-Mse --test cmiknn --tau-max 6 --sig-samples 200

夜間・切断耐性:
    nohup python -m japanflux_pn.run_robustness --site JP-Tak \\
        --outroot ~/bakanposs/japanflux_pn/outputs_robust > robust_Tak.log 2>&1 &
"""

from __future__ import annotations

import argparse
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig
from . import causal_network as cn


# inspect_site の year-scan で 7・8 月とも健全だった候補年 (欠測年は driver が自動スキップ)
SITE_YEARS: dict[str, tuple[list[int], list[int]]] = {
    "JP-Tak": (list(range(1999, 2010)) + list(range(2012, 2022)), [7, 8]),
    "JP-Mse": (list(range(2002, 2010)), [7, 8]),
    "JP-BBY": ([2015, 2016, 2017, 2019, 2020], [7, 8]),
}


def collect_links_for_years(
    site: str, years: list[int], months: list[int], test: str,
    tau_max: int, pc_alpha: float, config: AnalysisConfig,
    sig_samples: int, knn: float, max_conds_dim: int | None,
) -> dict[int, pd.DataFrame | None]:
    """年ごとに PCMCI の有向リンクを取る。欠測年やエラー年は None。"""
    from .preprocess import load_corevars_hh

    out: dict[int, pd.DataFrame | None] = {}
    for y in years:
        t0 = time.time()
        try:
            pre = load_corevars_hh(site, y, months, config)
            if not bool(pre.valid.all()):
                print(f"  {site} {y}: 欠測 {int((~pre.valid).sum())} 点 → skip", flush=True)
                out[y] = None
                continue
            results, _, var_names = cn.run_pcmci(
                pre, tau_max=tau_max, pc_alpha=pc_alpha, test=test, knn=knn,
                sig_samples=sig_samples, max_conds_dim=max_conds_dim, verbosity=0)
            links = cn.extract_links(results, var_names, config)
            n_dir = int((links["kind"] == "directed").sum()) if not links.empty else 0
            out[y] = links
            print(f"  {site} {y}: {n_dir} 有向リンク ({time.time()-t0:.0f}s)", flush=True)
        except Exception as e:  # noqa: BLE001
            print(f"  {site} {y}: ERROR {type(e).__name__}: {e} → skip", flush=True)
            out[y] = None
    return out


def aggregate(per_year: dict[int, pd.DataFrame | None]) -> tuple[pd.DataFrame, int]:
    """(src→dst) ごとに出現年数・頻度・平均強度・代表ラグを集計する。"""
    valid_years = [y for y, v in per_year.items() if v is not None]
    n_total = len(valid_years)
    rec: dict[tuple[str, str], dict] = defaultdict(
        lambda: {"years": set(), "strengths": [], "lags": []})
    for y, links in per_year.items():
        if links is None or links.empty:
            continue
        d = links[links["kind"] == "directed"]
        for (s, t), grp in d.groupby(["src", "dst"]):
            best = grp["strength"].abs().to_numpy().argmax()
            rec[(s, t)]["years"].add(y)
            rec[(s, t)]["strengths"].append(float(grp["strength"].abs().iloc[best]))
            rec[(s, t)]["lags"].append(float(grp["lag_h"].iloc[best]))

    rows = []
    for (s, t), v in rec.items():
        rows.append({
            "src": s, "dst": t,
            "n_years": len(v["years"]),
            "frequency": len(v["years"]) / max(n_total, 1),
            "mean_strength": float(np.mean(v["strengths"])),
            "median_lag_h": float(np.median(v["lags"])),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["frequency", "mean_strength"], ascending=False)
        df = df.reset_index(drop=True)
    return df, n_total


def report(site: str, test: str = "parcorr", tau_max: int = 6, pc_alpha: float = 0.01,
           sig_samples: int = 200, knn: float = 0.1, max_conds_dim: int | None = 3,
           config: AnalysisConfig | None = None,
           outroot: str | Path | None = None) -> pd.DataFrame:
    config = config or AnalysisConfig()
    years, months = SITE_YEARS[site]
    print(f"===== {site} 複数年ロバスト性 (test={test}, tau_max={tau_max}) =====")
    print(f"  候補年 {years[0]}–{years[-1]} ({len(years)} 年), months={months}\n", flush=True)

    per_year = collect_links_for_years(
        site, years, months, test, tau_max, pc_alpha, config,
        sig_samples, knn, max_conds_dim)
    tbl, n_total = aggregate(per_year)

    print(f"\n=== {site}: 有向因果リンクの年々一貫性 (有効 {n_total} 年) ===")
    if tbl.empty:
        print("  (リンク無し)")
        return tbl
    print(f"  {'link':<12} {'出現':>8} {'頻度':>6} {'平均|強度|':>10} {'代表lag':>8}")
    for _, r in tbl.iterrows():
        arrow = f"{r['src']}→{r['dst']}"
        flag = "  ⚠→Rg(外生)" if r["dst"] == "Rg" else (
               "  ★コア" if r["frequency"] >= 0.7 else "")
        print(f"  {arrow:<12} {int(r['n_years']):3d}/{n_total:<3d} "
              f"{r['frequency']:5.0%} {r['mean_strength']:10.3f} "
              f"{r['median_lag_h']:6.1f}h{flag}")

    robust = tbl[(tbl["frequency"] >= 0.7) & (tbl["dst"] != "Rg")]
    print(f"\n  [ロバストなコア骨格] 頻度≥70%: "
          + (", ".join(f"{r['src']}→{r['dst']}({r['frequency']:.0%})"
                       for _, r in robust.iterrows()) or "(該当なし)"))

    if outroot is not None:
        outdir = Path(outroot)
        outdir.mkdir(parents=True, exist_ok=True)
        tbl.to_csv(outdir / f"{site}_link_consistency_{test}.csv", index=False)
        print(f"\n[output] {outdir}/{site}_link_consistency_{test}.csv")
    return tbl


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="複数年ロバスト性 (因果リンクの年々一貫性)")
    p.add_argument("--site", required=True, choices=list(SITE_YEARS))
    p.add_argument("--test", default="parcorr", choices=["parcorr", "cmiknn"])
    p.add_argument("--tau-max", type=int, default=6)
    p.add_argument("--pc-alpha", type=float, default=0.01)
    p.add_argument("--sig-samples", type=int, default=200)
    p.add_argument("--knn", type=float, default=0.1)
    p.add_argument("--max-conds-dim", type=int, default=3)
    p.add_argument("--outroot", default=None)
    args = p.parse_args(argv)
    report(args.site, args.test, args.tau_max, args.pc_alpha, args.sig_samples,
           args.knn, args.max_conds_dim, outroot=args.outroot)


if __name__ == "__main__":
    main()
