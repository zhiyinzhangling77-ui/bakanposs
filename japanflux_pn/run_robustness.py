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
    # 2 つ目の森林（常緑針葉樹林）。実データ範囲が不明なので広めに取り、
    # 有効点数 <500 の年は各 driver が自動スキップする。inspect で健全年を確認後、
    # ここを精緻化するとよい。
    "JP-Ta2": (list(range(1994, 2016)), [7, 8]),
    # --- A3 気候勾配拡張: 乾燥・大陸性の草原 (JapanFlux2024 内, 11/11 マッピング) ---
    # 「湿潤日本 ↔ 乾燥モンゴル/青海」の気候軸。inspect_site で 7・8 月とも健全だった年。
    "CN-HaM": ([2002, 2003, 2004, 2005, 2006, 2008, 2009, 2012, 2013, 2014], [7, 8]),
    "MN-Kbu": ([2003, 2004, 2005, 2006, 2007], [7, 8]),
}

_AUTO_YEARS_CACHE: dict[tuple[str, tuple[int, ...]], tuple[list[int], list[int]]] = {}


def get_site_years(site: str, months: list[int] | None = None
                   ) -> tuple[list[int], list[int]]:
    """解析対象年を返す。手登録 :data:`SITE_YEARS` を優先し、無ければ自動検出。

    自動検出はデータの観測期間（ファイルの TIMESTAMP 範囲）から全年を採り、months は
    既定 [7,8]。有効点数 <500 の疎な年は各ドライバが自動スキップするので、範囲を広めに
    採ってよい（JP-Ta2 で確立した方針を全 11/11 サイトへ一般化）。これで rank_sites が
    示した 31 サイトを手登録なしで robustness / oinfo / climate に流せる。
    """
    if site in SITE_YEARS:
        return SITE_YEARS[site]
    months = list(months) if months else [7, 8]
    key = (site, tuple(months))
    if key in _AUTO_YEARS_CACHE:
        return _AUTO_YEARS_CACHE[key]
    import re
    from .sites import get_site
    from .preprocess import find_corevars_files
    from . import inspect_site as insp
    spec = get_site(site)
    files = find_corevars_files(spec)
    if files[0].suffix.lower() in (".xlsx", ".xls"):
        # xlsx (ChinaFlux 年別 / KoFlux は年範囲入りの 1 ファイル)。TIMESTAMP を
        # 軽く読めないので、ファイル名の 4 桁年から範囲を採る。疎な年は driver が skip。
        yrs = set()
        for f in files:
            for m in re.findall(r"((?:19|20)\d{2})", f.name):
                y = int(m)
                if 1990 <= y <= 2035:
                    yrs.add(y)
        res = (list(range(min(yrs), max(yrs) + 1)) if yrs
               else list(range(1990, 2025)), months)
    else:
        lo, _ = insp._timestamp_span(files[0])
        _, hi = insp._timestamp_span(files[-1])
        res = (list(range(lo.year, hi.year + 1)), months)
    _AUTO_YEARS_CACHE[key] = res
    return res


def collect_links_for_years(
    site: str, years: list[int], months: list[int], test: str,
    tau_max: int, pc_alpha: float, config: AnalysisConfig,
    sig_samples: int, knn: float, max_conds_dim: int | None,
) -> dict[int, pd.DataFrame | None]:
    """年ごとに PCMCI の有向リンクを取る。欠測年やエラー年は None。

    巨大な COREVARS CSV は 1 回だけ読み (load_raw_all)、各年はメモリ上で切り出す。
    """
    from .preprocess import (load_raw_all, slice_and_anomaly,
                             slice_span_and_anomaly, PreprocessResult)
    from .sites import get_site

    site_spec = get_site(site)
    print(f"  [読込] {site} の COREVARS を 1 回だけロード中...", flush=True)
    t_load = time.time()
    raw_all = load_raw_all(site_spec, config)          # ← 1 回だけ
    print(f"  [読込] 完了 ({time.time()-t_load:.0f}s)\n", flush=True)

    out: dict[int, pd.DataFrame | None] = {}
    for y in years:
        t0 = time.time()
        try:
            if len(months) == 1:
                anom, valid = slice_and_anomaly(raw_all, y, months[0], config)
            else:
                anom, valid = slice_span_and_anomaly(raw_all, y, months, config)
            pre = PreprocessResult(anomaly=anom, valid=valid, site=site, year=y,
                                   month=months[0], config=config, months=months)
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
           outroot: str | Path | None = None,
           years: list[int] | None = None) -> pd.DataFrame:
    config = config or AnalysisConfig()
    auto_years, months = get_site_years(site)
    years = years if years else auto_years   # --years で代表数年に絞れる (CMIknn 用)
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
    p.add_argument("--site", required=True,
                   help="サイトコード (手登録外でも 11/11 なら自動で健全年検出)")
    p.add_argument("--test", default="parcorr", choices=["parcorr", "cmiknn"])
    p.add_argument("--tau-max", type=int, default=6)
    p.add_argument("--pc-alpha", type=float, default=0.01)
    p.add_argument("--sig-samples", type=int, default=200)
    p.add_argument("--knn", type=float, default=0.1)
    p.add_argument("--max-conds-dim", type=int, default=3)
    p.add_argument("--outroot", default=None)
    p.add_argument("--years", type=int, nargs="+", default=None,
                   help="対象年を明示 (省略時は自動の健全年全部)。CMIknn を数年に絞る用。")
    args = p.parse_args(argv)
    report(args.site, args.test, args.tau_max, args.pc_alpha, args.sig_samples,
           args.knn, args.max_conds_dim, outroot=args.outroot, years=args.years)


if __name__ == "__main__":
    main()
