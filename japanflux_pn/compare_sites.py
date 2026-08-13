"""多サイト横断比較：各サイトの robustness CSV を集め、因果リンクの生態系間の異同を見る。

run_robustness が各サイトに保存する `<site>_link_consistency_<test>.csv`（列 src,dst,frequency…）
を複数集め、リンク×サイトの出現頻度マトリクスにする。全サイトでコア（≥70%）に入るリンクは
**普遍的な骨格**、一部だけのものは**生態系固有**の候補。

    # ディレクトリ内の CSV を全部集める
    python -m japanflux_pn.compare_sites --dir outputs_robust_parcorr --heatmap cross_sites.png

    # CSV を明示
    python -m japanflux_pn.compare_sites --csvs A_link_consistency_parcorr.csv B_...csv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


def _site_from_name(path) -> str:
    name = Path(path).stem
    m = re.match(r"(.+?)_link_consistency", name)
    return m.group(1) if m else name


def build_matrix(csvs: list[str], core: float = 0.7) -> pd.DataFrame:
    """リンク(src→dst) × サイト の頻度マトリクス。欠損は 0（そのサイトで出なかった）。"""
    cols: dict[str, dict[str, float]] = {}
    for c in csvs:
        df = pd.read_csv(c)
        site = _site_from_name(c)
        cols[site] = {f"{r.src}→{r.dst}": float(r.frequency)
                      for r in df.itertuples(index=False)}
    mat = pd.DataFrame(cols).fillna(0.0)
    # 「何サイトでコアか」→「平均頻度」で降順ソート（普遍リンクを上に）
    n_core = (mat >= core).sum(axis=1)
    mat = (mat.assign(_n=n_core, _m=mat.mean(axis=1))
              .sort_values(["_n", "_m"], ascending=False)
              .drop(columns=["_n", "_m"]))
    return mat


def draw_heatmap(mat: pd.DataFrame, path, core: float, top: int = 22) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    jp_path = "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"
    jp = fm.FontProperties(fname=jp_path) if Path(jp_path).exists() else None
    m = mat.head(top)
    fig, ax = plt.subplots(figsize=(1.4 + 0.9 * m.shape[1], 0.42 * m.shape[0] + 1.5))
    im = ax.imshow(m.to_numpy() * 100, cmap="YlGnBu", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=30, ha="right")
    ax.set_yticks(range(m.shape[0]))
    ax.set_yticklabels(m.index, fontproperties=jp)   # リンク名（→は IPAGothic にある）
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = int(round(m.iat[i, j] * 100))
            ax.text(j, i, v, ha="center", va="center", fontsize=8,
                    color="white" if v >= 55 else "#333")
    ax.set_title(f"因果リンクの生態系間比較（出現頻度 %／コア={core:.0%}）",
                 fontproperties=jp)
    cb = fig.colorbar(im, ax=ax, fraction=0.03)
    cb.set_label("出現頻度 %", fontproperties=jp)
    fig.savefig(path, bbox_inches="tight", dpi=130); plt.close(fig)


def report(csvs: list[str], core: float = 0.7, out_csv=None, heatmap=None) -> pd.DataFrame:
    mat = build_matrix(csvs, core)
    sites = list(mat.columns)
    n_core = (mat >= core).sum(axis=1)
    universal = n_core[n_core == len(sites)].index.tolist()
    specific = n_core[(n_core >= 1) & (n_core < len(sites))].index.tolist()

    print(f"=== 多サイト横断比較（{len(sites)} サイト: {', '.join(sites)}）===\n")
    print(f"■ 全サイトでコア（≥{core:.0%}）＝普遍的な骨格: {len(universal)} 本")
    for l in universal:
        print(f"    {l:16s}  " + "  ".join(f"{s}:{mat.loc[l, s]*100:3.0f}%" for s in sites))
    print(f"\n■ 一部サイトのみコア＝生態系固有の候補: {len(specific)} 本")
    for l in specific[:20]:
        hi = [s for s in sites if mat.loc[l, s] >= core]
        print(f"    {l:16s}  コア: {', '.join(hi)}")

    print("\n--- 出現頻度マトリクス（%、上位）---")
    print((mat.head(22) * 100).round(0).astype(int).to_string())

    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        mat.to_csv(out_csv)
        print(f"\n[saved] {out_csv}")
    if heatmap:
        draw_heatmap(mat, heatmap, core)
        print(f"[saved] {heatmap}")
    return mat


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="多サイト横断: 因果リンクの生態系間比較")
    p.add_argument("--csvs", nargs="+", default=None, help="robustness CSV を列挙")
    p.add_argument("--dir", default=None, help="このディレクトリ内の *_link_consistency_*.csv を集める")
    p.add_argument("--core", type=float, default=0.7)
    p.add_argument("--out-csv", default=None)
    p.add_argument("--heatmap", default=None)
    a = p.parse_args(argv)
    csvs = list(a.csvs or [])
    if a.dir:
        csvs += [str(x) for x in sorted(Path(a.dir).glob("*_link_consistency_*.csv"))]
    if not csvs:
        p.error("--csvs か --dir を指定してください")
    report(csvs, a.core, a.out_csv, a.heatmap)


if __name__ == "__main__":
    main()
