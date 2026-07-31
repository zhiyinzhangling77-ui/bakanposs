"""気候ストレスに対する情報構造の応答 — 情報ネットワークによる機能ストレス検出。

仮説: 水ストレス年（高 VPD・低土壌水分）には気孔閉鎖で光合成が光に律速されなくなり、
放射→光合成の情報結合 I(Rg; GPP) が弱まる。呼吸は水分律速が強まり I(θ; GER) が変化する。

各サイト・各健全年で「ストレス指標（生の平均 VPD/土壌水分/気温・降水量）」と「結合指標
（アノマリ空間の相互情報 I(Rg;GPP) 等）」を出し、年々の Spearman 相関を取る。単変数の
平均フラックスでは見えない、情報構造の気候応答を捉える。

    python -m japanflux_pn.climate_response --site JP-Tak

平均フラックスと違い、ここで見るのは「結合の強さが気候でどう組み変わるか」＝機能的可塑性。
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig, RK_VARS
from .sites import get_site
from .preprocess import (load_raw_all, slice_and_anomaly, slice_span_and_anomaly)
from . import information_theory as it
from .run_robustness import get_site_years


# ストレス指標と結合指標
STRESS_VARS = ["Ta", "VPD", "th", "P"]           # 生の平均/合計 (実レベル)
COUPLINGS = [("Rg", "GEP"), ("Rg", "gLE"), ("Rg", "NEE"),
             ("Ta", "GER"), ("th", "GER"), ("VPD", "gLE")]


def _spearman(x: np.ndarray, y: np.ndarray) -> float:
    """順位相関 (タイは平均順位, pandas.rank)。"""
    rx = pd.Series(x).rank().to_numpy()
    ry = pd.Series(y).rank().to_numpy()
    if np.std(rx) == 0 or np.std(ry) == 0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def _spearman_p(x: np.ndarray, y: np.ndarray, n_perm: int = 5000,
                seed: int = 0) -> tuple[float, float]:
    """Spearman r と、順列検定による両側 p 値 (年ラベルをシャッフル)。"""
    r_obs = _spearman(x, y)
    if not np.isfinite(r_obs):
        return r_obs, np.nan
    rng = np.random.default_rng(seed)
    y = np.asarray(y, dtype=float)
    count = 0
    for _ in range(n_perm):
        if abs(_spearman(x, rng.permutation(y))) >= abs(r_obs) - 1e-12:
            count += 1
    return r_obs, (count + 1) / (n_perm + 1)


def year_metrics(raw_all: pd.DataFrame, year: int, months: list[int],
                 config: AnalysisConfig) -> dict | None:
    """1 年のストレス指標（生平均）と結合指標（アノマリ MI）をまとめる。"""
    months = sorted(months)
    start = pd.Timestamp(year=year, month=months[0], day=1)
    end = pd.Timestamp(year=year, month=months[-1], day=1) + pd.offsets.MonthBegin(1)
    raw = raw_all[(raw_all.index >= start) & (raw_all.index < end)]
    if raw.empty:
        return None
    stress = {
        "Ta_mean": float(raw["Ta"].mean()),
        "VPD_mean": float(raw["VPD"].mean()),
        "SWC_mean": float(raw["th"].mean()),        # 土壌水分 θ
        "P_sum": float(raw["P"].sum()),
    }

    if len(months) == 1:
        anom, valid = slice_and_anomaly(raw_all, year, months[0], config)
    else:
        anom, valid = slice_span_and_anomaly(raw_all, year, months, config)
    vf = anom.loc[valid]
    if len(vf) < 200:                                # 情報量推定に最低限
        return None
    m = config.n_bins
    logm = config.log_m
    idx = {v: it.digitize_series(vf[v].to_numpy(dtype=float), m) for v in RK_VARS}

    def MI(a: str, b: str) -> float:
        return 100.0 * it.mutual_information_indices(idx[a], idx[b], m, True) / logm

    coup = {f"I_{a}_{b}": MI(a, b) for a, b in COUPLINGS}
    return {"year": year, "n": len(vf), **stress, **coup}


def scan(site: str, config: AnalysisConfig | None = None) -> pd.DataFrame:
    config = config or AnalysisConfig()
    years, months = get_site_years(site)
    raw_all = load_raw_all(get_site(site), config)
    rows = []
    for y in years:
        t0 = time.time()
        r = year_metrics(raw_all, y, months, config)
        if r is None:
            print(f"  {site} {y}: skip (欠測/点数不足)", flush=True)
            continue
        rows.append(r)
        print(f"  {site} {y}: n={r['n']} VPD={r['VPD_mean']:.2f} "
              f"SWC={r['SWC_mean']:.1f} I(Rg;GEP)={r['I_Rg_GEP']:.1f} "
              f"({time.time()-t0:.0f}s)", flush=True)
    return pd.DataFrame(rows)


def report(site: str, config: AnalysisConfig | None = None,
           outroot: str | Path | None = None) -> pd.DataFrame:
    config = config or AnalysisConfig()
    print(f"===== {site} 気候ストレスへの情報構造応答 =====")
    df = scan(site, config)
    if len(df) < 4:
        print(f"\n  有効年 {len(df)} < 4。相関評価に不足。")
        return df

    print(f"\n=== 年々 Spearman 相関 (結合 vs ストレス指標, 有効 {len(df)} 年) ===")
    print(f"  {'結合':<12} {'vs VPD':>8} {'vs SWC':>8} {'vs Ta':>8}")
    for a, b in COUPLINGS:
        col = f"I_{a}_{b}"
        rV = _spearman(df[col].to_numpy(), df["VPD_mean"].to_numpy())
        rS = _spearman(df[col].to_numpy(), df["SWC_mean"].to_numpy())
        rT = _spearman(df[col].to_numpy(), df["Ta_mean"].to_numpy())
        print(f"  {a+'→'+b:<12} {rV:8.2f} {rS:8.2f} {rT:8.2f}")

    # 目玉: 光利用結合 I(Rg;GPP) が乾燥で下がるか (順列検定 p 値つき)
    rV, pV = _spearman_p(df["I_Rg_GEP"].to_numpy(), df["VPD_mean"].to_numpy())
    rS, pS = _spearman_p(df["I_Rg_GEP"].to_numpy(), df["SWC_mean"].to_numpy())
    print(f"\n  [仮説検定] I(Rg;GPP) vs 乾燥 (順列検定 5000):")
    print(f"    VPD↑ で結合↓ (負相関を期待): r={rV:+.2f}  p={pV:.3f}")
    print(f"    土壌水分↑ で結合↑ (正相関を期待): r={rS:+.2f}  p={pS:.3f}")
    sig = (rV <= -0.4 and pV < 0.05) or (rS >= 0.4 and pS < 0.05)
    verdict = ("支持・有意 (乾燥で光合成が放射から脱結合)" if sig
               else ("傾向あり・非有意" if (rV <= -0.4 or rS >= 0.4)
                     else "不明瞭 (このサイト/期間では弱い)"))
    print(f"    → {verdict}")

    if outroot is not None:
        outdir = Path(outroot)
        outdir.mkdir(parents=True, exist_ok=True)
        df.to_csv(outdir / f"{site}_climate_response.csv", index=False)
        print(f"\n[output] {outdir}/{site}_climate_response.csv")
    return df


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="気候ストレスへの情報構造応答")
    p.add_argument("--site", required=True,
                   help="サイトコード (手登録外でも 11/11 なら自動で健全年検出)")
    p.add_argument("--outroot", default=None)
    args = p.parse_args(argv)
    report(args.site, outroot=args.outroot)


if __name__ == "__main__":
    main()
