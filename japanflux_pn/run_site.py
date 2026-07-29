"""エントリポイント: サイト名を引数に取り、一連の解析を回す。

    python -m japanflux_pn.run_site --site JP-Tak --year 2003 --month 7

前処理 → ネットワーク構築 → 隣接行列 CSV とネットワーク図・診断プロットを
``japanflux_pn/outputs/<site>_<year><month>/`` に出力する。R&K Bondville との
定性比較 (sanity check) をコンソールに表示する。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from .config import AnalysisConfig
from .preprocess import load_corevars_hh
from .network import build_network, NetworkResult


# sanity check で見る主要ペア (R&K Bondville の期待挙動)
SANITY_PAIRS: list[tuple[str, str, str]] = [
    ("Ta", "Rg", "Ta↔Rg: type 1 想定 (強同期, 情報流ほぼ無し)"),
    ("Ts", "GER", "Ts→GER: 土壌温度→呼吸の強い共変 (type 1-3)"),
    ("VPD", "gLE", "VPD→LE: 有意な情報流を想定"),
    ("GEP", "NEE", "GEP↔NEE: type 2 (双方向フィードバック) 想定"),
    ("P", "th", "P→θ: 降水→土壌水分"),
]

# 診断プロットに出すペア
DIAG_PAIRS: list[tuple[str, str]] = [
    ("VPD", "gLE"), ("NEE", "GEP"), ("GEP", "NEE"), ("P", "th"),
    ("Ts", "GER"), ("Rg", "Ta"),
]


def run(site: str, year: int, month: int | list[int],
        config: AnalysisConfig | None = None,
        outroot: str | Path | None = None) -> NetworkResult:
    config = config or AnalysisConfig()
    pre = load_corevars_hh(site, year, month, config)
    if pre.n_points < 100:
        print(f"[warn] 有効データ点数 {pre.n_points} が少ない (推奨 500-1500)")
    print(f"[preprocess] {site} {year}-{pre.month_label}: n_points={pre.n_points}")

    net = build_network(pre)

    outroot = Path(outroot) if outroot else Path(__file__).parent / "outputs"
    outdir = outroot / f"{site}_{year}{pre.month_label}"
    outdir.mkdir(parents=True, exist_ok=True)
    net.save(outdir)
    print(f"[output] {outdir}  (adjacency CSV + meta)")
    _draw_figures(net, outdir)

    _print_sanity(net)
    return net


def _draw_figures(net: NetworkResult, outdir: Path) -> None:
    """描画ライブラリがあれば図を出す。無ければ警告してスキップ (解析本体は済)。"""
    try:
        from . import viz
    except ImportError as e:  # matplotlib / networkx 未導入
        print(f"[warn] 図の生成をスキップ: {e.name} が未導入 "
              f"(`pip install networkx matplotlib` で有効化)")
        return
    viz.draw_network(net, outdir / "network.png")
    viz.plot_type_matrix(net, outdir / "coupling_type.png")
    viz.plot_lag_diagnostics(net, DIAG_PAIRS, outdir / "lag_diagnostics.png")
    print(f"[output] {outdir}  (figures: network / coupling_type / lag_diagnostics .png)")


def _print_sanity(net: NetworkResult) -> None:
    tnames = {1: "1 sync", 2: "2 feedback", 3: "3 forcing", 4: "4 uncoupled"}
    print("\n=== sanity check (R&K Bondville との定性比較) ===")
    for src, dst, note in SANITY_PAIRS:
        t = int(net.ctype.loc[src, dst])
        tz = net.ATz.loc[src, dst]
        lag = net.Gamma.loc[src, dst]
        ai = net.AI.loc[src, dst]
        rev = int(net.ctype.loc[dst, src])
        tz_s = f"{tz:.2f}" if np.isfinite(tz) else "  - "
        lag_s = f"{lag:.1f}h" if np.isfinite(lag) else "  -  "
        print(f"  {src:>3}→{dst:<3} type={tnames[t]:<11} "
              f"Tz={tz_s} τ'={lag_s} I'={ai:4.1f}%  (逆向 type={tnames[rev]})")
        print(f"        └ {note}")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="JapanFlux2024 process-network analysis")
    p.add_argument("--site", required=True, help="サイトコード (例 JP-Tak)")
    p.add_argument("--year", type=int, required=True)
    p.add_argument("--month", type=int, nargs="+", default=[7],
                   help="対象月。単月 (既定 7) か連続複数月 (例: --month 7 8) で "
                        "プールし n を稼ぐ")
    p.add_argument("--bins", type=int, default=None, help="ビン数 m (既定 11)")
    p.add_argument("--surrogates", type=int, default=None, help="サロゲート数 M (既定 100)")
    p.add_argument("--peak-min-run", type=int, default=None,
                   help="連続有意ラグ帯の最小長 (既定 1=忠実R&K, 実データ推奨 2)")
    p.add_argument("--qc-max", type=int, default=None,
                   help="QC≤この値のみ残す (0=実測, 1=実測+良質補完)。既定 None=gap-fill込み")
    p.add_argument("--outroot", default=None)
    args = p.parse_args(argv)

    kw = {}
    if args.bins is not None:
        kw["n_bins"] = args.bins
    if args.surrogates is not None:
        kw["n_surrogates"] = args.surrogates
    if args.peak_min_run is not None:
        kw["peak_min_run"] = args.peak_min_run
    if args.qc_max is not None:
        kw["qc_max"] = args.qc_max
    config = AnalysisConfig(**kw)
    run(args.site, args.year, args.month, config, args.outroot)


if __name__ == "__main__":
    main()
