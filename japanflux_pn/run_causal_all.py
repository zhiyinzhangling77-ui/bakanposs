"""3 サイトの PCMCI+ 因果ネットワークを一括・完全版で回す (夜間実行向け)。

- 完全版 (簡略化なし): tau_max=36 (18h), sig_samples=500, knn=0.1, 条件次元無制限
- サイトごとに CSV 保存 (チェックポイント)。途中で落ちても済んだサイトは残る
- 1 サイトが失敗しても残りを続行 (クラッシュ隔離)
- PC 位相の変数ごとに進捗 %・経過・ETA を表示、サイト間の残り時間も外挿

夜間・切断耐性のある実行 (推奨):

    nohup python -m japanflux_pn.run_causal_all --test cmiknn \\
        --outroot ~/bakanposs/japanflux_pn/outputs_pcmci_full \\
        > pcmci_overnight.log 2>&1 &

    tail -f pcmci_overnight.log      # 進捗を眺める
"""

from __future__ import annotations

import argparse
import time
import traceback
from pathlib import Path

from .config import AnalysisConfig
from . import causal_network as cn


# 既定バッチ: 森林 / 水田 / 湿原 の健全年 7+8 月プール
DEFAULT_SPECS: list[tuple[str, int, list[int]]] = [
    ("JP-Tak", 2003, [7, 8]),
    ("JP-Mse", 2003, [7, 8]),
    ("JP-BBY", 2015, [7, 8]),
]


def run_batch(specs, test, tau_max, pc_alpha, sig_samples, knn, max_conds_dim,
              outroot) -> None:
    config = AnalysisConfig()
    t_start = time.time()
    site_times: list[float] = []
    n = len(specs)
    print(f"===== PCMCI+ バッチ ({n} サイト, test={test}) =====")
    print(f"  完全版設定: tau_max={tau_max or config.lag_max} "
          f"sig_samples={sig_samples} knn={knn} max_conds_dim={max_conds_dim}")
    print(f"  開始 {time.strftime('%Y-%m-%d %H:%M:%S')}\n", flush=True)

    for i, (site, year, months) in enumerate(specs):
        head = f"[{i+1}/{n}] {site} {year} months={months}"
        print(f"\n########## {head} ##########", flush=True)
        t0 = time.time()
        try:
            cn.report(site, year, months, test=test, pc_alpha=pc_alpha,
                      tau_max=tau_max, sig_samples=sig_samples, knn=knn,
                      max_conds_dim=max_conds_dim, config=config, outroot=outroot)
        except Exception as e:  # noqa: BLE001  クラッシュ隔離
            print(f"\n[ERROR] {site} 失敗: {type(e).__name__}: {e}", flush=True)
            traceback.print_exc()
            continue
        dt = time.time() - t0
        site_times.append(dt)
        avg = sum(site_times) / len(site_times)
        remaining = n - (i + 1)
        print(f"\n[site done] {site} 所要 {dt/60:.1f}min "
              f"({time.strftime('%H:%M:%S')})", flush=True)
        if remaining:
            print(f"[cross-site ETA] 残り {remaining} サイト ~{avg*remaining/60:.1f}min",
                  flush=True)

    print(f"\n===== 全 {n} サイト完了。総計 {(time.time()-t_start)/60:.1f}min "
          f"({time.strftime('%Y-%m-%d %H:%M:%S')}) =====", flush=True)


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="PCMCI+ 3 サイト一括 (完全版)")
    p.add_argument("--test", default="cmiknn", choices=["parcorr", "cmiknn"])
    p.add_argument("--tau-max", type=int, default=None,
                   help="最大ラグ (step)。既定は config の lag_max=36 (完全版)")
    p.add_argument("--pc-alpha", type=float, default=0.01)
    p.add_argument("--sig-samples", type=int, default=500,
                   help="cmiknn シャッフル数 (完全版 500)")
    p.add_argument("--knn", type=float, default=0.1)
    p.add_argument("--max-conds-dim", type=int, default=None,
                   help="条件次元上限 (既定 無制限 = 完全版)")
    p.add_argument("--outroot", default=None, required=True,
                   help="サイトごとの CSV 保存先 (チェックポイント)")
    args = p.parse_args(argv)
    run_batch(DEFAULT_SPECS, args.test, args.tau_max, args.pc_alpha,
              args.sig_samples, args.knn, args.max_conds_dim, args.outroot)


if __name__ == "__main__":
    main()
