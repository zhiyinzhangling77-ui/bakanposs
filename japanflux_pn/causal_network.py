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
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig, RK_VARS, RK_LABELS
from .preprocess import load_corevars_hh, PreprocessResult


def _patch_numpy_corrcoef() -> None:
    """NumPy 2.x で削除された corrcoef の ddof/bias 引数を握り潰す互換シム。

    tigramite 5.2 のブロックシャッフル有意化 (_get_acf) が ``np.corrcoef(..., ddof=0)``
    を呼ぶが、NumPy 2.0 で ddof/bias は削除された (元々 no-op)。cmiknn を通すため、
    これらのキーワードを落として本来の corrcoef に委譲する。
    """
    if getattr(np.corrcoef, "_ddof_shim", False):
        return
    _orig = np.corrcoef

    def _compat(*args, **kwargs):
        kwargs.pop("ddof", None)
        kwargs.pop("bias", None)
        return _orig(*args, **kwargs)

    _compat._ddof_shim = True
    np.corrcoef = _compat


def _make_progress_pcmci(PCMCI):
    """PC-stable の変数ごとに進捗 (k/N・経過・ETA) を出す PCMCI サブクラス。

    ``_run_pc_stable_single`` は各ターゲット変数の条件選択を担い、cmiknn ではここが
    支配的コスト。1 変数終わるごとに経過時間から ETA を外挿して表示する。位相が
    向き付け (orientation) に移ると N を超えるので、その旨を出す。
    """
    class _ProgressPCMCI(PCMCI):
        def _run_pc_stable_single(self, j, *args, **kwargs):
            if not hasattr(self, "_prog_t0"):
                self._prog_t0 = time.time()
                self._prog_done = 0
            res = super()._run_pc_stable_single(j, *args, **kwargs)
            self._prog_done += 1
            n = self.N
            el = time.time() - self._prog_t0
            k = self._prog_done
            if k <= n:
                rate = el / k
                eta = rate * (n - k)
                print(f"  [PC {100*k/n:3.0f}%] {k}/{n} vars | "
                      f"elapsed {el/60:5.1f}min | ETA(PC) ~{eta/60:5.1f}min "
                      f"({rate/60:.1f}min/var)", flush=True)
            else:
                print(f"  [orient] test {k-n} | elapsed {el/60:5.1f}min", flush=True)
            return res
    return _ProgressPCMCI


def _require_tigramite():
    try:
        from tigramite.pcmci import PCMCI
        from tigramite import data_processing as tp
        from tigramite.independence_tests.parcorr import ParCorr
    except ImportError as e:  # noqa: BLE001
        raise ImportError(
            "tigramite が未導入です。`pip install tigramite` を実行してください。"
        ) from e
    cmiknn_err = None
    try:
        from tigramite.independence_tests.cmiknn import CMIknn
    except Exception as e:  # noqa: BLE001
        CMIknn = None
        cmiknn_err = e
    return PCMCI, tp, ParCorr, CMIknn, cmiknn_err


def run_pcmci(
    pre: PreprocessResult,
    tau_max: int | None = None,
    pc_alpha: float = 0.01,
    test: str = "parcorr",
    knn: float = 0.1,
    sig_samples: int = 250,
    max_conds_dim: int | None = None,
    verbosity: int = 0,
):
    """PCMCI+ を回して (results, pcmci, var_names) を返す。

    完全被覆 (欠測無し) のレギュラ系列を前提。健全年 7+8 月プールがこれに該当。
    """
    PCMCI, tp, ParCorr, CMIknn, cmiknn_err = _require_tigramite()
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
            raise ImportError(
                "CMIknn を import できません。多くは numba/scikit-learn 未導入が原因: "
                "`pip install numba scikit-learn`。 元エラー: "
                f"{type(cmiknn_err).__name__}: {cmiknn_err}"
            )
        _patch_numpy_corrcoef()   # NumPy 2.x とのブロックシャッフル互換
        cond = CMIknn(significance="shuffle_test", knn=knn, sig_samples=sig_samples)
    elif test == "parcorr":
        cond = ParCorr(significance="analytic")
    else:
        raise ValueError(f"unknown test {test!r}; 'parcorr' or 'cmiknn'")

    PCMCIcls = _make_progress_pcmci(PCMCI) if verbosity >= 1 else PCMCI
    pcmci = PCMCIcls(dataframe=dataframe, cond_ind_test=cond, verbosity=0)
    results = pcmci.run_pcmciplus(tau_min=1, tau_max=tau_max, pc_alpha=pc_alpha,
                                  max_conds_dim=max_conds_dim)
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
                if i == j:
                    kind = "auto"          # 自己相関 (自己回帰項)
                elif directed:
                    kind = "directed"      # 変数間の有向因果
                else:
                    kind = "contemp_undirected"
                rows.append({
                    "src": var_names[i],
                    "dst": var_names[j],
                    "lag_h": config.lag_hours(tau),
                    "strength": float(val[i, j, tau]),
                    "kind": kind,
                })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.reindex(df["strength"].abs().sort_values(ascending=False).index)
        df = df.reset_index(drop=True)
    return df


def report(site: str, year: int, months: list[int], test: str = "parcorr",
           pc_alpha: float = 0.01, tau_max: int | None = None,
           sig_samples: int = 250, max_conds_dim: int | None = None,
           knn: float = 0.1, config: AnalysisConfig | None = None,
           outroot: str | Path | None = None) -> pd.DataFrame:
    config = config or AnalysisConfig()
    pre = load_corevars_hh(site, year, months, config)
    tau_max = int(tau_max if tau_max is not None else config.lag_max)
    print(f"[preprocess] {site} {year}-{pre.month_label}: n_points={pre.n_points}")
    print(f"[pcmci] test={test} pc_alpha={pc_alpha} tau_max={tau_max} "
          f"({config.lag_hours(tau_max):.0f} h)"
          + (f" sig_samples={sig_samples} knn={knn} "
             f"max_conds_dim={max_conds_dim}" if test == "cmiknn" else ""))

    # cmiknn は重いので進捗を出す (verbosity=1)
    results, pcmci, var_names = run_pcmci(
        pre, tau_max=tau_max, pc_alpha=pc_alpha, test=test, knn=knn,
        sig_samples=sig_samples, max_conds_dim=max_conds_dim,
        verbosity=1 if test == "cmiknn" else 0)
    links = extract_links(results, var_names, config)

    if links.empty:
        print("\n(有意なリンク無し)")
        return links
    directed = links[links["kind"] == "directed"]
    contemp = links[links["kind"] == "contemp_undirected"]
    auto = links[links["kind"] == "auto"]

    # 変数間の有向因果 (自己相関を除いた本命)
    print(f"\n=== PCMCI+ 変数間 有向因果リンク ({len(directed)} 本) ===")
    print(f"  (自己相関 X→X {len(auto)} 本は除外)")
    print(f"  {'link':<14} {'lag':>6} {'|strength|':>10}")
    for _, r in directed.iterrows():
        arrow = f"{RK_LABELS[r['src']]}→{RK_LABELS[r['dst']]}"
        flag = "  ⚠逆向き(Rg外生)" if r["dst"] == "Rg" else ""
        print(f"  {arrow:<14} {r['lag_h']:5.1f}h {abs(r['strength']):10.3f}{flag}")

    # ハブ構造: 出次数 (何変数を駆動するか) / 入次数
    out_deg = directed.groupby("src")["dst"].nunique().sort_values(ascending=False)
    in_deg = directed.groupby("dst")["src"].nunique().sort_values(ascending=False)
    print(f"\n  [ソースハブ] 出次数: "
          + ", ".join(f"{RK_LABELS[v]}={d}" for v, d in out_deg.head(4).items()))
    print(f"  [シンクハブ] 入次数: "
          + ", ".join(f"{RK_LABELS[v]}={d}" for v, d in in_deg.head(4).items()))

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
    p.add_argument("--tau-max", type=int, default=None,
                   help="最大ラグ (step)。既定は config の lag_max=36。cmiknn は "
                        "重いので 6 など短縮推奨")
    p.add_argument("--sig-samples", type=int, default=250,
                   help="cmiknn シャッフル数 (既定 250)。速度と精度のトレードオフ")
    p.add_argument("--knn", type=float, default=0.1,
                   help="cmiknn 近傍数 (割合 or 整数, 既定 0.1)。小さいほど速い")
    p.add_argument("--max-conds-dim", type=int, default=None,
                   help="PC-stable の条件次元上限 (小さいほど速い、既定 無制限)")
    p.add_argument("--outroot", default=None)
    args = p.parse_args(argv)
    report(args.site, args.year, args.month, args.test, args.pc_alpha,
           tau_max=args.tau_max, sig_samples=args.sig_samples,
           max_conds_dim=args.max_conds_dim, knn=args.knn, outroot=args.outroot)


if __name__ == "__main__":
    main()
