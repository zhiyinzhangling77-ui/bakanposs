"""実データ投入前の下調べ: 列名マッピングの実地確認と、健全年の探索。

ローカル (/mnt/hdd が見える環境) で:

    python -m japanflux_pn.inspect_site --site JP-Tak
    python -m japanflux_pn.inspect_site --site JP-Tak --months 7 8 --years 2000-2023

やること:
1. COREVARS HH ファイル一覧と観測期間を表示
2. 11 変数 (RK_VARS) の既定マッピングが実ヘッダに存在するか確認。欠けている
   ものには、FLUXNET2015 の同系プレフィックスから override 候補を提示
3. 全変数が揃っていれば、年 × 月ごとに前処理後の有効点数 (listwise) を出し、
   目標 500-1500 に収まる健全年を探す。各変数の生欠測率も併記

この段階では :mod:`sites` を書き換えず、``var_overrides`` に貼るべき提案だけ出す。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig, RK_VARS
from .sites import BASE_RH_COL, SiteSpec, DEFAULT_VAR_MAP, VPD_FROM_TA_RH, get_site
from .preprocess import (
    find_corevars_files, compute_anomaly, _regular_grid, _vpd_from_ta_rh,
)


# 欠けた変数の override 候補を探すためのプレフィックス (FLUXNET2015 命名)
CANDIDATE_PREFIX: dict[str, tuple[str, ...]] = {
    "Rg":  ("SW_IN",),
    "Ta":  ("TA",),
    "VPD": ("VPD",),
    "Ts":  ("TS",),
    "P":   ("P_", "P"),
    "th":  ("SWC",),
    "gH":  ("H_", "H"),
    "gLE": ("LE",),
    "GER": ("RECO",),
    "NEE": ("NEE",),
    "GEP": ("GPP",),
}


def _read_header(path: Path) -> list[str]:
    return list(pd.read_csv(path, nrows=0).columns)


def _timestamp_span(path: Path) -> tuple[pd.Timestamp, pd.Timestamp]:
    df = pd.read_csv(path, usecols=["TIMESTAMP_START"])
    ts = pd.to_datetime(df["TIMESTAMP_START"].astype("int64").astype(str),
                        format="%Y%m%d%H%M")
    return ts.min(), ts.max()


def check_mapping(header: list[str], site: SiteSpec) -> tuple[dict, dict]:
    """R&K 変数ごとに (present) と、欠けたものの候補列を返す。"""
    hset = set(header)
    vmap = site.var_map()
    # 導出変数 (VPD@BASE) は TA+RH が揃えば「present」扱い。
    present = {v: vmap[v] for v in RK_VARS
               if vmap[v] in hset or vmap[v].startswith("@")}
    missing = {}
    for v in RK_VARS:
        if vmap[v] in hset or vmap[v].startswith("@"):
            continue
        cands = sorted(
            c for c in header
            if any(c.upper().startswith(p.upper()) for p in CANDIDATE_PREFIX[v])
        )
        missing[v] = {"expected": vmap[v], "candidates": cands}
    return present, missing


def _load_tolerant(files: list[Path], columns: list[str], na: float) -> pd.DataFrame:
    """存在する列のみ寛容に読み込み、RK 表記 index の生データを返す。"""
    frames = []
    for f in files:
        have = [c for c in columns if c in set(_read_header(f))]
        usecols = ["TIMESTAMP_START"] + have
        df = pd.read_csv(f, usecols=lambda c: c in usecols)
        ts = pd.to_datetime(df["TIMESTAMP_START"].astype("int64").astype(str),
                            format="%Y%m%d%H%M")
        df = df.drop(columns=["TIMESTAMP_START"])
        df.index = ts
        frames.append(df)
    out = pd.concat(frames)
    out = out[~out.index.duplicated(keep="first")].sort_index()
    return out.replace(na, np.nan)


def year_scan(
    site: SiteSpec,
    config: AnalysisConfig,
    years: list[int] | None,
    months: list[int],
) -> pd.DataFrame:
    """年 × 月ごとに前処理後の有効点数と各変数の生欠測率を集計する。"""
    vmap = site.var_map()
    derived_vpd = vmap["VPD"] == VPD_FROM_TA_RH
    real_cols = [vmap[v] for v in RK_VARS if not vmap[v].startswith("@")]
    if derived_vpd:                     # VPD 導出のため RH も読む
        real_cols = real_cols + [BASE_RH_COL]
    files = find_corevars_files(site)
    raw_real = _load_tolerant(files, real_cols, config.na_sentinel)
    if derived_vpd and vmap["Ta"] in raw_real.columns and BASE_RH_COL in raw_real.columns:
        raw_real[vmap["VPD"]] = _vpd_from_ta_rh(
            raw_real[vmap["Ta"]].to_numpy(dtype=float),
            raw_real[BASE_RH_COL].to_numpy(dtype=float),
        )
    # 実カラム名 → R&K 表記へ (存在する列のみ; VPD マーカも含む)
    inv = {vmap[v]: v for v in RK_VARS if vmap[v] in raw_real.columns}
    raw = raw_real.rename(columns=inv)
    for v in RK_VARS:  # 欠けた変数は全 NaN 列として補い、listwise を厳密化
        if v not in raw.columns:
            raw[v] = np.nan
    raw = raw[RK_VARS]

    if years is None:
        years = sorted({t.year for t in raw.index})
    step = pd.Timedelta(minutes=24 * 60 // config.steps_per_day)

    rows = []
    for y in years:
        for mo in months:
            try:
                m_start = pd.Timestamp(year=y, month=mo, day=1)
            except ValueError:
                continue
            m_end = m_start + pd.offsets.MonthBegin(1)
            buf_end = m_end + pd.Timedelta(days=config.anomaly_window_days)
            grid = _regular_grid(m_start, buf_end - step, config.steps_per_day)
            sub = raw.reindex(grid)
            if sub.notna().any(axis=None) is False or sub.dropna(how="all").empty:
                continue
            keep = _regular_grid(m_start, m_end - step, config.steps_per_day)
            anom, valid = compute_anomaly(sub, config, keep_index=keep)
            raw_month = raw.reindex(keep)
            cover = {f"cov_{v}": float(raw_month[v].notna().mean()) for v in RK_VARS}
            rows.append({"year": y, "month": mo, "n_points": int(valid.sum()),
                         "n_grid": len(keep), **cover})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# レポート出力
# ---------------------------------------------------------------------------
def report(site_code: str, months: list[int], years: list[int] | None,
           config: AnalysisConfig | None = None) -> None:
    config = config or AnalysisConfig()
    site = get_site(site_code)
    print(f"### site {site_code}  ({site.description})")
    print(f"### data_dir: {site.data_dir}\n")

    files = find_corevars_files(site)
    print(f"[files] {len(files)} {site.source} HH file(s):")
    for f in files:
        try:
            lo, hi = _timestamp_span(f)
            span = f"{lo:%Y-%m-%d} .. {hi:%Y-%m-%d}"
        except Exception as e:  # noqa: BLE001
            span = f"(span unreadable: {e})"
        print(f"   {f.name}   [{span}]")

    header = _read_header(files[0])
    present, missing = check_mapping(header, site)
    print(f"\n[mapping] {len(present)}/{len(RK_VARS)} variables mapped OK")
    for v in RK_VARS:
        if v in present:
            print(f"   OK  {v:>3} -> {present[v]}")
        else:
            exp = missing[v]["expected"]
            cands = missing[v]["candidates"]
            print(f"   !!  {v:>3} -> {exp}  MISSING. candidates: {cands or '(none)'}")

    if missing:
        print("\n[action] sites.py の SiteSpec に var_overrides を追加してください、例:")
        ov = {v: (missing[v]["candidates"][0] if missing[v]["candidates"] else "???")
              for v in missing}
        print(f"   var_overrides={ov!r}")
        print("   ↑ 候補を確認して確定後、再度 inspect を実行して year-scan に進む。")
        return

    print(f"\n[year-scan] months={months}  (target n_points 500-1500)")
    tbl = year_scan(site, config, years, months)
    if tbl.empty:
        print("   (no data in requested years/months)")
        return
    cov_cols = [c for c in tbl.columns if c.startswith("cov_")]
    tbl = tbl.sort_values(["month", "year"]).reset_index(drop=True)
    for _, r in tbl.iterrows():
        n = int(r["n_points"])
        flag = "  <-- OK (500-1500)" if 500 <= n <= 1500 else (
               "  (<500 sparse)" if n < 500 else "  (>1500 dense)")
        worst = min(cov_cols, key=lambda c: r[c])
        print(f"   {int(r['year'])}-{int(r['month']):02d}: "
              f"n_points={n:4d}/{int(r['n_grid'])}{flag}   "
              f"limiting cov: {worst[4:]}={r[worst]:.0%}")

    ok = tbl[(tbl["n_points"] >= 500) & (tbl["n_points"] <= 1500)]
    if not ok.empty:
        best = ok.loc[ok["n_points"].idxmax()]
        print(f"\n[suggest] 健全年候補 (最大有効点数): "
              f"{int(best['year'])}-{int(best['month']):02d} "
              f"n_points={int(best['n_points'])}")
        print(f"   → python -m japanflux_pn.run_site --site {site_code} "
              f"--year {int(best['year'])} --month {int(best['month'])} --peak-min-run 2")


def _parse_years(s: str | None) -> list[int] | None:
    if not s:
        return None
    out: list[int] = []
    for part in s.split(","):
        part = part.strip()
        if "-" in part:
            a, b = part.split("-")
            out.extend(range(int(a), int(b) + 1))
        elif part:
            out.append(int(part))
    return out or None


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="JapanFlux2024 site inspection")
    p.add_argument("--site", required=True)
    p.add_argument("--months", type=int, nargs="+", default=[7, 8])
    p.add_argument("--years", type=str, default=None,
                   help="例 2003 / 2000-2010 / 2003,2005,2008 (既定=全年)")
    args = p.parse_args(argv)
    report(args.site, args.months, _parse_years(args.years))


if __name__ == "__main__":
    main()
