"""複数サイト一括の気候応答: 光利用脱結合 I(Rg;GPP)〜乾燥 をサイト横断で検定。

JP-Tef (n=21) で有意 (VPD r=-0.44 p=0.049, SWC r=+0.49 p=0.027) だったように、
脱結合の検出には長期記録 (検定力) が要る。全 11/11 サイトで年々 Spearman 相関と
順列 p 値を出し、「どの生態系・どの記録長で脱結合が有意か」を 1 表にまとめる。

    python -m japanflux_pn.batch_climate --all --csv batch_climate.csv
    python -m japanflux_pn.batch_climate --sites JP-Tak JP-Tef JP-Tmd JP-SMF JP-Mse

各サイト: 有効年数, I(Rg;GPP) vs VPD の (r,p), vs SWC の (r,p), 判定。記録が長いほど
検定力が高い点に注意 (n<4 は相関不可)。生態系ラベルは ecosystem.py で別途付与できる。
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig
from .sites import JAPANFLUX_ROOT
from .climate_response import scan, _spearman_p


def run(sites: list[str], config: AnalysisConfig | None = None) -> pd.DataFrame:
    config = config or AnalysisConfig()
    rows = []
    for i, site in enumerate(sites, 1):
        t0 = time.time()
        try:
            df = scan(site, config)
        except Exception as e:  # noqa: BLE001
            print(f"  [{i}/{len(sites)}] {site}: skip ({type(e).__name__}: {e})",
                  flush=True)
            continue
        n = len(df)
        if n < 4:
            print(f"  [{i}/{len(sites)}] {site}: 有効年 {n} < 4 (相関不可)", flush=True)
            continue
        rV, pV = _spearman_p(df["I_Rg_GEP"].to_numpy(), df["VPD_mean"].to_numpy())
        rS, pS = _spearman_p(df["I_Rg_GEP"].to_numpy(), df["SWC_mean"].to_numpy())
        sig = (rV <= -0.4 and pV < 0.05) or (rS >= 0.4 and pS < 0.05)
        rows.append({"site": site, "n_years": n,
                     "r_VPD": rV, "p_VPD": pV, "r_SWC": rS, "p_SWC": pS,
                     "decoupling_sig": bool(sig)})
        print(f"  [{i}/{len(sites)}] {site}: n={n:>2} "
              f"VPD r={rV:+.2f} p={pV:.3f} | SWC r={rS:+.2f} p={pS:.3f} "
              f"{'★有意' if sig else ''} ({time.time()-t0:.0f}s)", flush=True)
    return pd.DataFrame(rows)


def report(sites: list[str], config: AnalysisConfig | None = None,
           csv: str | None = None) -> pd.DataFrame:
    print(f"### 一括 気候応答 (光利用脱結合) ({len(sites)} サイト)\n", flush=True)
    df = run(sites, config)
    if df.empty:
        print("  (結果なし)")
        return df

    df = df.sort_values(["decoupling_sig", "n_years"], ascending=False)
    print("\n=== I(Rg;GPP) 〜 乾燥 の脱結合 (記録長順・有意優先) ===")
    print(f"  {'site':<10} {'年数':>4} {'r(VPD)':>7} {'p(VPD)':>7} "
          f"{'r(SWC)':>7} {'p(SWC)':>7}  判定")
    for _, r in df.iterrows():
        mark = "★有意" if r["decoupling_sig"] else (
            "傾向" if (r["r_VPD"] <= -0.4 or r["r_SWC"] >= 0.4) else "")
        print(f"  {r['site']:<10} {int(r['n_years']):>4} "
              f"{r['r_VPD']:>7.2f} {r['p_VPD']:>7.3f} "
              f"{r['r_SWC']:>7.2f} {r['p_SWC']:>7.3f}  {mark}")

    n_sig = int(df["decoupling_sig"].sum())
    long = df[df["n_years"] >= 15]
    n_sig_long = int(long["decoupling_sig"].sum())
    print(f"\n  脱結合が有意: {n_sig}/{len(df)} サイト "
          f"(長期 n≥15 に限ると {n_sig_long}/{len(long)})")
    print("  ※ 記録が長いほど検定力が高い。生態系タイプ別集計は ecosystem.py と結合を。")

    if csv:
        Path(csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(csv, index=False)
        print(f"\n[output] {csv}")
    return df


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="batch multi-site climate decoupling")
    p.add_argument("--sites", nargs="+", default=None)
    p.add_argument("--all", action="store_true",
                   help="11/11 マッピングの全サイトを対象")
    p.add_argument("--root", default=JAPANFLUX_ROOT)
    p.add_argument("--csv", default=None)
    args = p.parse_args(argv)
    if args.all:
        from .batch_oinfo import _eleven_mapped_sites
        sites = _eleven_mapped_sites(args.root)
        print(f"[all] 11/11 サイト {len(sites)} 件\n")
    elif args.sites:
        sites = args.sites
    else:
        p.error("--sites か --all を指定してください")
    report(sites, csv=args.csv)


if __name__ == "__main__":
    main()
