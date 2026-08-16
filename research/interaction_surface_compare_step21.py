"""旗21：呼吸の"掛け算モデルの破れ"（Ta×θ 非加法）は生態系をまたいで一般的か。

旗20 で JP-Tak の呼吸は掛け算モデル R=f(Ta)g(θ) を有意に破る（deyear 0.139, p=0.002）と分かった。
これが自然生態系で一般的か（JP-Ta2 常緑林・CN-HaM 高山草原）／水田でどうか（JP-Mse）を、
各サイト健全年プール＋年レベル除去(deyear)＋順列検定で並べる。旗20 の関数を再利用。

    python research/interaction_surface_compare_step21.py --nperm 500 --fig interaction_cross.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from interaction_surface_step20 import interaction_fraction, surrogate_pvalue

DEFAULT_SITES = ["JP-Tak", "JP-Ta2", "CN-HaM", "JP-Mse"]


def load_pooled(site, months, deyear=True):
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    from japanflux_pn.run_robustness import get_site_years
    import pandas as pd
    cfg = AnalysisConfig()
    years, mo = get_site_years(site)
    months = months or mo
    ms = sorted(months)
    raw_all = load_raw_all(get_site(site), cfg)
    Ta, th, GER = [], [], []
    used = 0
    for y in years:
        start = pd.Timestamp(year=y, month=ms[0], day=1)
        end = pd.Timestamp(year=y, month=ms[-1], day=1) + pd.offsets.MonthBegin(1)
        r = raw_all[(raw_all.index >= start) & (raw_all.index < end)]
        if r.empty:
            continue
        g = r["GER"].to_numpy(float)
        if deyear:
            gp = g[np.isfinite(g) & (g > 0)]
            gm = np.exp(np.mean(np.log(gp))) if gp.size else 1.0
            g = g / gm
        GER.append(g); Ta.append(r["Ta"].to_numpy(float)); th.append(r["th"].to_numpy(float))
        used += 1
    if not used:
        return None
    return (np.concatenate(Ta), np.concatenate(th), np.concatenate(GER), used)


def draw(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    jpf = "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"
    jp = fm.FontProperties(fname=jpf) if Path(jpf).exists() else None
    sites = [r["site"] for r in rows]
    fr = [r["frac"] for r in rows]
    nul = [r["null"] for r in rows]
    cols = ["#c0392b" if r["paddy"] else "#1f7a3d" for r in rows]
    x = np.arange(len(sites))
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.bar(x, fr, color=cols, width=0.6, label="観測")
    ax.plot(x, nul, "k_", ms=18, mew=2, label="順列ヌル平均")
    for i, r in enumerate(rows):
        star = "*" if (np.isfinite(r["p"]) and r["p"] < 0.05) else ""
        ax.text(i, fr[i] + 0.005, f"{fr[i]:.2f}{star}", ha="center", fontproperties=jp, fontsize=10)
    ax.axhline(0.05, color="#999", ls=":", lw=1)
    ax.set_xticks(x); ax.set_xticklabels(sites, fontproperties=jp)
    ax.set_ylabel("掛け算モデルからの交互作用割合", fontproperties=jp)
    ax.set_title("呼吸の掛け算モデルの破れ（Ta×θ 非加法）は生態系をまたぐか\n"
                 "（* = 順列検定 p<0.05, 年レベル除去済み）", fontproperties=jp)
    ax.legend(prop=jp, frameon=False)
    fig.savefig(path, bbox_inches="tight", dpi=140); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description="掛け算モデルの破れの生態系間比較")
    p.add_argument("--sites", nargs="+", default=DEFAULT_SITES)
    p.add_argument("--month", type=int, nargs="+", default=None)
    p.add_argument("--nbins", type=int, default=6)
    p.add_argument("--min-cell", type=int, default=20)
    p.add_argument("--nperm", type=int, default=500)
    p.add_argument("--no-deyear", action="store_true")
    p.add_argument("--fig", default=None)
    a = p.parse_args()
    paddy = {"JP-Mse", "KR-CRK"}

    print("=== 呼吸の掛け算モデルの破れ（Ta×θ 非加法）の生態系間比較（旗21）===")
    print("  各サイト健全年プール・年レベル除去(deyear)・順列検定。frac=掛け算からのズレ\n")
    print(f"  {'サイト':<8} {'年':>3} {'交互作用':>8} {'ヌル平均':>8} {'p':>7}  判定")
    rows = []
    for s in a.sites:
        d = load_pooled(s, a.month, deyear=not a.no_deyear)
        if d is None:
            print(f"  {s:<8} 有効年なし"); continue
        Ta, th, GER, used = d
        if a.nperm > 0:
            frac, null, pv = surrogate_pvalue(Ta, th, GER, a.nbins, a.min_cell, a.nperm)
        else:
            frac = interaction_fraction(Ta, th, GER, a.nbins, a.min_cell)[0]; null = pv = np.nan
        sig = (np.isfinite(pv) and pv < 0.05)
        mark = ("★有意に破れる" if (sig and frac >= 0.05) else
                "有意だが小" if sig else "非有意（掛け算で概ね書ける）")
        rows.append({"site": s, "frac": frac, "null": null, "p": pv,
                     "paddy": s in paddy, "used": used})
        print(f"  {s:<8} {used:>3} {frac:>8.3f} {null:>8.3f} {pv:>7.3f}  {mark}")

    nat = [r for r in rows if not r["paddy"]]
    nsig = sum(1 for r in nat if np.isfinite(r["p"]) and r["p"] < 0.05 and r["frac"] >= 0.05)
    print(f"\n  自然生態系で有意に破れる：{nsig}/{len(nat)}")
    if nat and nsig == len(nat):
        print("    → ✅ 掛け算モデルの破れは自然生態系で一般的（旗20 が転移）")
    elif nsig >= 1:
        print("    → ○ 一部で有意（サイト/年を増やして確認）")
    print("  留保：GER は NEE 分割の派生量（分割は Ta 主効果を押し付けうるが Ta×θ 交互作用は押し付けない）。"
          "生プール・セル平均応答曲面。")

    if a.fig and rows:
        draw(rows, a.fig); print(f"\n  [図] {a.fig}")


if __name__ == "__main__":
    main()
