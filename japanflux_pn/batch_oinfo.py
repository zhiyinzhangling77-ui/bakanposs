"""複数サイト一括 O-information: 高次相乗の"生態系タイプ越え反復"を1表にまとめる。

「自然生態系は地下の高次相乗を持ち、管理水田は失う」という主張を、n=2 森林でなく
11/11 マッピングの全サイト (rank_sites 上位) で検定する。各サイトについて全健全年で
O-info z を計算し、サブシステムごとに **相乗年/全年 (z≤−2.36) の割合** を出す。

    # 手元の 11/11 全サイトを一括 (時間がかかる。進捗表示あり)
    python -m japanflux_pn.batch_oinfo --all --csv batch_oinfo.csv
    # サイトを明示
    python -m japanflux_pn.batch_oinfo --sites JP-Tak JP-Tef JP-Fhk JP-Mse CN-HaM MN-Kbu

出力: サイト×サブシステムの相乗割合表と、鍵サブシステム (呼吸制御/土壌熱) のまとめ。
site-year 数が検定力を決めるので、健全年の多いサイト (JP-Tak 44, JP-Tef 40…) ほど強い。
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig
from .sites import JAPANFLUX_ROOT
from .oinfo_analysis import scan_years, SUBSYSTEMS, OINFO_BINS_DEFAULT

SYN_Z = -2.36   # z ≤ この値で「相乗」(片側 α=0.01)
RED_Z = 2.36    # z ≥ この値で「冗長」
# 地下クラスタ (主張の核) を最初に置く
KEY_SUBSYSTEMS = ["呼吸制御 Rg,Ta,θ,GER", "土壌熱 Ta,Ts,θ,GER"]


def _eleven_mapped_sites(root: str) -> list[str]:
    """11/11 マッピングのサイトのみ返す (欠測サイトは O-info を計算できない)。"""
    from .rank_sites import _resolve_sites
    from .preprocess import find_corevars_files
    from . import inspect_site as insp
    out = []
    for code, site in sorted(_resolve_sites(root).items()):
        try:
            header = insp._read_header(find_corevars_files(site)[0])
            present, _ = insp.check_mapping(header, site)
            if len(present) == 11:
                out.append(code)
        except Exception:  # noqa: BLE001
            continue
    return out


def _summarize(long: pd.DataFrame) -> dict[str, dict]:
    """1 サイトの scan_years 出力 → サブシステムごとの相乗割合など。"""
    res = {}
    for name, g in long.groupby("subsystem", sort=False):
        tot = len(g)
        n_syn = int((g["z"] <= SYN_Z).sum())
        n_red = int((g["z"] >= RED_Z).sum())
        res[name] = {
            "n_years": tot, "n_syn": n_syn, "n_red": n_red,
            "syn_frac": n_syn / tot if tot else np.nan,
            "mean_z": float(g["z"].mean()) if tot else np.nan,
        }
    return res


def run(sites: list[str], obins: int = OINFO_BINS_DEFAULT,
        config: AnalysisConfig | None = None) -> pd.DataFrame:
    config = config or AnalysisConfig()
    rows = []
    for i, site in enumerate(sites, 1):
        t0 = time.time()
        try:
            long = scan_years(site, obins, config)
        except Exception as e:  # noqa: BLE001
            print(f"  [{i}/{len(sites)}] {site}: skip ({type(e).__name__}: {e})",
                  flush=True)
            continue
        if long.empty:
            print(f"  [{i}/{len(sites)}] {site}: 有効年なし", flush=True)
            continue
        summ = _summarize(long)
        for name, s in summ.items():
            rows.append({"site": site, "subsystem": name, **s})
        key = summ.get(KEY_SUBSYSTEMS[0], {})
        print(f"  [{i}/{len(sites)}] {site}: 呼吸制御 相乗 "
              f"{key.get('n_syn',0)}/{key.get('n_years',0)} "
              f"(mean z={key.get('mean_z',float('nan')):+.1f})  "
              f"({time.time()-t0:.0f}s)", flush=True)
    return pd.DataFrame(rows)


def report(sites: list[str], obins: int = OINFO_BINS_DEFAULT,
           config: AnalysisConfig | None = None,
           csv: str | None = None) -> pd.DataFrame:
    print(f"### 一括 O-information ({len(sites)} サイト, obins={obins})\n", flush=True)
    df = run(sites, obins, config)
    if df.empty:
        print("  (結果なし)")
        return df

    for key in KEY_SUBSYSTEMS:
        sub = df[df["subsystem"] == key].sort_values("syn_frac", ascending=False)
        if sub.empty:
            continue
        print(f"\n=== {key}: サイト別の相乗安定性 ===")
        print(f"  {'site':<10} {'相乗年/全年':>10} {'相乗割合':>8} {'mean z':>8}  判定")
        for _, r in sub.iterrows():
            frac = r["syn_frac"]
            verdict = ("★相乗が安定" if frac >= 0.6 else
                       "冗長寄り" if r["n_red"] >= 0.6 * r["n_years"] else "混在")
            print(f"  {r['site']:<10} {int(r['n_syn']):>4}/{int(r['n_years']):<4} "
                  f"{frac:>7.0%} {r['mean_z']:>8.1f}  {verdict}")
        n_syn_sites = int((sub["syn_frac"] >= 0.6).sum())
        print(f"  → 相乗が安定なサイト: {n_syn_sites}/{len(sub)}")

    if csv:
        Path(csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(csv, index=False)
        print(f"\n[output] {csv}")
    return df


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="batch multi-site O-information")
    p.add_argument("--sites", nargs="+", default=None, help="サイトコード列")
    p.add_argument("--all", action="store_true",
                   help="11/11 マッピングの全サイトを対象")
    p.add_argument("--root", default=JAPANFLUX_ROOT)
    p.add_argument("--obins", type=int, default=OINFO_BINS_DEFAULT)
    p.add_argument("--csv", default=None)
    args = p.parse_args(argv)
    if args.all:
        sites = _eleven_mapped_sites(args.root)
        print(f"[all] 11/11 サイト {len(sites)} 件: {', '.join(sites)}\n")
    elif args.sites:
        sites = args.sites
    else:
        p.error("--sites か --all を指定してください")
    report(sites, args.obins, csv=args.csv)


if __name__ == "__main__":
    main()
