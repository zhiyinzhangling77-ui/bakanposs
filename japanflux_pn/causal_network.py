"""多変量因果ネットワーク: PCMCI+ (Tigramite) による有向グラフ推定。

ペアワイズ TE/MI は共通駆動 (放射 Rg) を分離できず、条件付き MI も固定ビン 3D の
疎性で弱結合が読めない。PCMCI+ は各候補リンクを「他の全変数の過去」で条件付けて
偽の共通駆動・間接経路を落とし、自己相関に頑健な有向グラフを与える
(Runge et al. 2019, Sci. Adv.; 地球科学時系列の因果探索の標準)。

依存 (任意): ``pip install tigramite``

    python -m japanflux_pn.causal_network --site JP-Tak --year 2003 --month 7 8
    python -m japanflux_pn.causal_network --site JP-Tak --year 2003 --month 7 8 --test cmiknn

7+8 月プールの健全サイト年は 30 分格子が完全被覆 (n=2976, 欠測無し) なので、アノマリ
系列をそのままレギュラ時系列として渡せる。欠測を含む年は Tigramite の mask 対応が
必要 (将来課題) なので、ここでは完全被覆を要求する。

test:
  parcorr … 線形部分相関 (高速, 既定)。線形近似の有向グラフ。
  cmiknn  … KNN ベース条件付き MI (非線形, 低速)。固定ビンの限界を脱する本命。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig, RK_VARS, RK_LABELS
from .preprocess import load_corevars_hh, PreprocessResult


def _require_tigramite():
    try:
        from tigramite.pcmci import PCMCI
        from tigramite import data_processing as tp
        from tigramite.independence_tests.parcorr import ParCorr
    except ImportError as e:  # noqa: BLE001
        raise ImportError(
            "tigramite が未導入です。`pip install tigramite` を実行してください。"
        ) from e
    try:
        from tigramite.independence_tests.cmiknn import CMIknn
    except Exception:  # noqa: BLE001
        CMIknn = None
    return PCMCI, tp, ParCorr, CMIknn


def run_pcmci(
    pre: PreprocessResult,
    tau_max: int | None = None,
    pc_alpha: float = 0.01,
    test: str = "parcorr",
    knn: float = 0.1,
    verbosity: int = 0,
):
    """PCMCI+ を回して (results, pcmci, var_names) を返す。

    完全被覆 (欠測無し) のレギュラ系列を前提。健全年 7+8 月プールがこれに該当。
    """
    PCMCI, tp, ParCorr, CMIknn = _require_tigramite()
    cfg = pre.config
    tau_max = int(tau_max if tau_max is not None else cfg.lag_max)

    if not bool(pre.valid.all()):
        n_bad = int((~pre.valid).sum())
        raise ValueError(
            f"欠測 {n_bad} 点あり。PCMCI アダプタは完全被覆のみ対応。健全年の "
            f"7+8 月プール (n_points == n_grid) を使うか、mask 対応を実装のこと。"
        )

    data = pre.anomaly[RK_VARS].to_numpy(dtype=float)
    dataframe = tp.DataFrame(data, var_names=list(RK_VARS))

    if test == "cmiknn":
        if CMIknn is None:
            raise ImportError("CMIknn が使えません (tigramite のバージョン確認)。")
        cond = CMIknn(significance="shuffle_test", knn=knn)
    elif test == "parcorr":
        cond = ParCorr(significance="analytic")
    else:
        raise ValueError(f"unknown test {test!r}; 'parcorr' or 'cmiknn'")

    pcmci = PCMCI(dataframe=dataframe, cond_ind_test=cond, verbosity=verbosity)
    results = pcmci.run_pcmciplus(tau_min=1, tau_max=tau_max, pc_alpha=pc_alpha)
    return results, pcmci, list(RK_VARS)


# ---------------------------------------------------------------------------
# 結果の解釈
# ---------------------------------------------------------------------------
def extract_links(results, var_names, config: AnalysisConfig) -> pd.DataFrame:
    """PCMCI+ の graph から有向リンク一覧 (src, dst, lag_h, strength) を作る。

    graph[i, j, tau] の '-->' を i(過去 tau) → j の有向因果とみなす。tau=0 の
    '-->' は向き付き同時因果、'o-o' は向き未確定の同時結合として別扱い。
    """
    graph = results["graph"]
    val = results["val_matrix"]
    n = len(var_names)
    tau_max = graph.shape[2] - 1

    rows = []
    for i in range(n):
        for j in range(n):
            for tau in range(0, tau_max + 1):
                mark = graph[i, j, tau]
                if mark == "":
                    continue
                if tau == 0 and i >= j:
                    continue  # 同時リンクの重複を避ける
                directed = mark == "-->"
                contemp_undirected = (tau == 0 and mark == "o-o")
                if not (directed or contemp_undirected):
                    continue
                rows.append({
                    "src": var_names[i],
                    "dst": var_names[j],
                    "lag_h": config.lag_hours(tau),
                    "strength": float(val[i, j, tau]),
                    "kind": "directed" if directed else "contemp_undirected",
                })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.reindex(df["strength"].abs().sort_values(ascending=False).index)
        df = df.reset_index(drop=True)
    return df


def report(site: str, year: int, months: list[int], test: str = "parcorr",
           pc_alpha: float = 0.01, config: AnalysisConfig | None = None,
           outroot: str | Path | None = None) -> pd.DataFrame:
    config = config or AnalysisConfig()
    pre = load_corevars_hh(site, year, months, config)
    print(f"[preprocess] {site} {year}-{pre.month_label}: n_points={pre.n_points}")
    print(f"[pcmci] test={test} pc_alpha={pc_alpha} tau_max={config.lag_max} "
          f"({config.lag_hours(config.lag_max):.0f} h)")

    results, pcmci, var_names = run_pcmci(pre, pc_alpha=pc_alpha, test=test)
    links = extract_links(results, var_names, config)

    directed = links[links["kind"] == "directed"] if not links.empty else links
    contemp = links[links["kind"] == "contemp_undirected"] if not links.empty else links
    print(f"\n=== PCMCI+ 有向因果リンク ({len(directed)} 本) ===")
    print(f"  {'link':<14} {'lag':>6} {'|strength|':>10}")
    for _, r in directed.iterrows():
        arrow = f"{RK_LABELS[r['src']]}→{RK_LABELS[r['dst']]}"
        print(f"  {arrow:<14} {r['lag_h']:5.1f}h {abs(r['strength']):10.3f}")
    if len(contemp):
        print(f"\n  同時 (向き未確定) {len(contemp)} 本: "
              + ", ".join(f"{RK_LABELS[r['src']]}–{RK_LABELS[r['dst']]}"
                          for _, r in contemp.iterrows()))

    if outroot is not None:
        outdir = Path(outroot) / f"{site}_{year}{pre.month_label}_pcmci_{test}"
        outdir.mkdir(parents=True, exist_ok=True)
        links.to_csv(outdir / "causal_links.csv", index=False)
        print(f"\n[output] {outdir}")
    return links


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="PCMCI+ 多変量因果ネットワーク")
    p.add_argument("--site", required=True)
    p.add_argument("--year", type=int, required=True)
    p.add_argument("--month", type=int, nargs="+", default=[7])
    p.add_argument("--test", default="parcorr", choices=["parcorr", "cmiknn"])
    p.add_argument("--pc-alpha", type=float, default=0.01)
    p.add_argument("--outroot", default=None)
    args = p.parse_args(argv)
    report(args.site, args.year, args.month, args.test, args.pc_alpha,
           outroot=args.outroot)


if __name__ == "__main__":
    main()
