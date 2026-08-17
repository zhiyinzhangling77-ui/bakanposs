"""全サイトを走査し「使える健全 site-year の実数」でランキングする。

名目の観測期間 (span) ではなく、**11 変数が listwise で揃い、目標点数 (既定 500–1500)
に収まる (年×月) の数**で並べる。検定力 (climate_response の年々相関, O-info の複数年
安定性) は site-year 数で決まるので、「どのサイトにデータが溜まっているか」を厳密に測る。

ローカル (/mnt/hdd が見える環境) で:

    python -m japanflux_pn.rank_sites
    python -m japanflux_pn.rank_sites --months 7 8 --top 25
    python -m japanflux_pn.rank_sites --months 7 8 --min 800 --max 3000 --csv rank.csv

各サイトについて: マッピング n/11, 健全 site-year 数, 期間, 律速変数を出す。マッピングが
11 未満のサイト (土壌欠測など) は健全数を測れないので理由付きで別掲する。走査は全 CSV を
読むため時間がかかる (サイトごとに進捗を表示)。
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd

from .config import AnalysisConfig, RK_VARS
from .sites import SITES, JAPANFLUX_ROOT, SiteSpec, discover_japanflux_sites
from .preprocess import find_corevars_files
from . import inspect_site as insp


def _resolve_sites(root: str) -> dict[str, SiteSpec]:
    """指定 root で発見したサイトに、手登録サイト (curated) を上書きした解決表。"""
    resolved = dict(discover_japanflux_sites(root))
    resolved.update(SITES)
    return resolved


def scan_all(
    root: str = JAPANFLUX_ROOT,
    months: list[int] | None = None,
    n_min: int = 500,
    n_max: int = 1500,
    config: AnalysisConfig | None = None,
    ebr: bool = False,
    ebr_qc: int = 1,
) -> pd.DataFrame:
    """全サイトを走査し、1 サイト 1 行の集計 DataFrame を返す。

    ``ebr=True`` で各サイトのエネルギー収支閉合 EBR も計算する（H/LE/Rn/G 列を追加で
    読むため遅くなる）。``ebr_qc`` は EBR 算出時の QC 閾値（既定 1＝実測寄り）。
    """
    months = months or [7, 8]
    config = config or AnalysisConfig()
    rows = []
    resolved = _resolve_sites(root)
    codes = sorted(resolved)
    for i, code in enumerate(codes, 1):
        t0 = time.time()
        rec = {"site": code, "n_mapped": 0, "healthy_site_years": 0,
               "years_span": "", "n_files": 0, "note": ""}
        if ebr:
            rec.update({"ebr": float("nan"), "ebr_bad_years": 0, "ebr_note": ""})
        try:
            site = resolved[code]
            files = find_corevars_files(site)
            rec["n_files"] = len(files)
            header = insp._read_header(files[0])
            present, missing = insp.check_mapping(header, site)
            rec["n_mapped"] = len(present)
            if missing:
                rec["note"] = "欠測: " + ",".join(sorted(missing))
            else:
                tbl = insp.year_scan(site, config, years=None, months=months)
                if not tbl.empty:
                    healthy = tbl[(tbl["n_points"] >= n_min) & (tbl["n_points"] <= n_max)]
                    rec["healthy_site_years"] = int(len(healthy))
                    yrs = sorted(set(int(y) for y in tbl["year"]))
                    if yrs:
                        rec["years_span"] = f"{yrs[0]}-{yrs[-1]}"
            if ebr:
                _fill_ebr(rec, site, months, ebr_qc)
        except Exception as e:  # noqa: BLE001
            rec["note"] = f"error: {type(e).__name__}: {e}"
        rows.append(rec)
        etag = (f" EBR={rec['ebr']:.2f}" if ebr and rec.get("ebr") == rec.get("ebr")
                else (f" EBR:{rec['ebr_note']}" if ebr else ""))
        print(f"  [{i:>3}/{len(codes)}] {code:<10} "
              f"map={rec['n_mapped']}/11 healthy={rec['healthy_site_years']:>3}"
              f"{etag} ({time.time()-t0:.0f}s) {rec['note']}", flush=True)
    df = pd.DataFrame(rows)
    return df.sort_values(["healthy_site_years", "n_mapped"],
                          ascending=False).reset_index(drop=True)


def _fill_ebr(rec: dict, site: SiteSpec, months: list[int], ebr_qc: int) -> None:
    """rec に EBR 情報を書き込む（閉合コアは energy_closure に集約）。"""
    from .energy_closure import site_ebr
    r = site_ebr(site, months, qc_max=ebr_qc)
    if "note" in r:
        rec["ebr_note"] = r["note"]
        return
    rec["ebr"] = round(r["ebr"], 3)
    rec["ebr_bad_years"] = int(r["n_bad"])
    rec["ebr_note"] = ("G=0近似" if r["g_approx"] else "") + \
        (f" 不閉合{r['n_bad']}/{r['n_years']}年" if r["n_bad"] else "")


def report(root: str = JAPANFLUX_ROOT, months: list[int] | None = None,
           n_min: int = 500, n_max: int = 1500, top: int | None = None,
           csv: str | None = None, ebr: bool = False, ebr_qc: int = 1) -> pd.DataFrame:
    months = months or [7, 8]
    print(f"### 全サイト走査 (months={months}, 健全域 {n_min}-{n_max} 点"
          f"{'・EBR併記' if ebr else ''})\n", flush=True)
    df = scan_all(root, months, n_min, n_max, ebr=ebr, ebr_qc=ebr_qc)

    full = df[df["n_mapped"] == 11].copy()
    partial = df[df["n_mapped"] < 11].copy()

    print(f"\n=== 11/11 マッピング済み × 健全 site-year 数 ランキング "
          f"(months={months}) ===")
    if ebr:
        print(f"  {'site':<10} {'健全年':>6} {'EBR':>6} {'閉合':>4} "
              f"{'期間':>12} {'files':>6}   EBR注")
    else:
        print(f"  {'site':<10} {'健全年':>6} {'期間':>12} {'files':>6}")
    shown = full if top is None else full.head(top)
    for _, r in shown.iterrows():
        if ebr:
            e = r.get("ebr")
            estr = f"{e:>6.2f}" if e == e else f"{'—':>6}"    # NaN 判定
            ok = "✅" if (e == e and e >= 0.7) else ("⚠" if e == e else "—")
            print(f"  {r['site']:<10} {int(r['healthy_site_years']):>6} {estr} "
                  f"{ok:>4} {r['years_span']:>12} {int(r['n_files']):>6}   "
                  f"{r.get('ebr_note', '')}")
        else:
            print(f"  {r['site']:<10} {int(r['healthy_site_years']):>6} "
                  f"{r['years_span']:>12} {int(r['n_files']):>6}")
    print(f"\n  合計 (11/11 サイト): {len(full)} サイト, "
          f"健全 site-year 総数 = {int(full['healthy_site_years'].sum())} "
          f"(×{len(months)} 月/年)")

    if not partial.empty:
        print(f"\n=== 変数欠測などで full 解析不可 ({len(partial)} サイト) ===")
        for _, r in partial.iterrows():
            print(f"  {r['site']:<10} map={int(r['n_mapped'])}/11  {r['note']}")

    if csv:
        Path(csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(csv, index=False)
        print(f"\n[output] {csv}")
    return df


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="rank sites by healthy site-years")
    p.add_argument("--root", default=JAPANFLUX_ROOT)
    p.add_argument("--months", type=int, nargs="+", default=[7, 8])
    p.add_argument("--min", dest="n_min", type=int, default=500)
    p.add_argument("--max", dest="n_max", type=int, default=1500)
    p.add_argument("--top", type=int, default=None)
    p.add_argument("--csv", default=None)
    p.add_argument("--ebr", action="store_true",
                   help="エネルギー収支閉合 EBR も計算（H/LE/Rn/G を追加で読むため遅い）")
    p.add_argument("--ebr-qc", type=int, default=1,
                   help="EBR 算出時の QC 閾値（既定 1＝実測寄り）")
    args = p.parse_args(argv)
    report(args.root, args.months, args.n_min, args.n_max, args.top, args.csv,
           ebr=args.ebr, ebr_qc=args.ebr_qc)


if __name__ == "__main__":
    main()
