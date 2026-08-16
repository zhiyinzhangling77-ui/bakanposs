"""旗27：呼吸の"水分依存 Q10"は長期森林で一般的か（JP-Tak 固有か検証）。

旗26 で JP-Tak は Q10 が土壌水分で増える(r=+0.90, QC でも +1.0)と分かった。これが
長期森林で一般的かを、rank_sites で健全年数の多い森林を並べて再現テストする。旗26 の
q10_by_moisture / _boot_trend を再利用。各サイト健全年プールで Q10 vs θ の Spearman＋CI を出す。

既定サイト＝健全年数≥20の11/11森林（rank_sites より）＋対照に草原/高山。

    python research/respiration_q10_compare_step27.py --fig q10_cross.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from respiration_q10_moisture_step26 import q10_by_moisture, _boot_trend

# 健全年数の多い 11/11 森林（rank_sites）＋対照（高山草原・草原ステップ）
DEFAULT_SITES = ["JP-Tak", "JP-Tef", "JP-Fhk", "JP-Tmd", "JP-SMF", "JP-Spp",
                 "JP-Ta2", "JP-Fjy", "CN-HaM", "MN-Kbu"]
FOREST = {"JP-Tak", "JP-Tef", "JP-Fhk", "JP-Tmd", "JP-SMF", "JP-Spp", "JP-Ta2", "JP-Fjy"}


def load_site(site, months, qc_max, cfg_cls, get_site, load_raw_all, get_site_years):
    import pandas as pd
    cfg = cfg_cls(qc_max=qc_max) if qc_max is not None else cfg_cls()
    years, mo = get_site_years(site)
    ms = sorted(months or mo)
    raw_all = load_raw_all(get_site(site), cfg)
    parts, used = [], 0
    for y in years:
        start = pd.Timestamp(year=y, month=ms[0], day=1)
        end = pd.Timestamp(year=y, month=ms[-1], day=1) + pd.offsets.MonthBegin(1)
        r = raw_all[(raw_all.index >= start) & (raw_all.index < end)]
        if not r.empty:
            parts.append(r[["GER", "Ta", "th"]]); used += 1
    if not parts:
        return None, 0
    return pd.concat(parts).dropna(), used


def draw(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    jpf = "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"
    jp = fm.FontProperties(fname=jpf) if Path(jpf).exists() else None
    sites = [r["site"] for r in rows]; rs = [r["r"] for r in rows]
    los = [r["ci"][0] if isinstance(r["ci"], tuple) else np.nan for r in rows]
    his = [r["ci"][1] if isinstance(r["ci"], tuple) else np.nan for r in rows]
    cols = ["#1f7a3d" if r["forest"] else "#b8860b" for r in rows]
    x = np.arange(len(sites))
    fig, ax = plt.subplots(figsize=(1.4 + 0.85 * len(sites), 4.6))
    for i, r in enumerate(rows):
        if isinstance(r["ci"], tuple):
            ax.plot([i, i], [los[i], his[i]], color=cols[i], lw=2, zorder=1)
    ax.scatter(x, rs, c=cols, s=70, zorder=3)
    ax.axhline(0, color="#999", lw=1)
    ax.set_xticks(x); ax.set_xticklabels(sites, rotation=30, ha="right", fontproperties=jp)
    ax.set_ylabel("Q10 vs 土壌水分 の Spearman r（＋95%CI）", fontproperties=jp)
    ax.set_title("呼吸の『水分依存 Q10』は長期森林で一般的か\n"
                 "（緑=森林 / 金=草原・高山, r>0 かつ CI が 0 を跨がない＝水分依存 Q10）",
                 fontproperties=jp)
    fig.savefig(path, bbox_inches="tight", dpi=140); plt.close(fig)


def main():
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    from japanflux_pn.run_robustness import get_site_years

    p = argparse.ArgumentParser(description="水分依存Q10は長期森林で一般的か")
    p.add_argument("--sites", nargs="+", default=DEFAULT_SITES)
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--nbin", type=int, default=5)
    p.add_argument("--qc-max", type=int, default=None)
    p.add_argument("--fig", default=None)
    a = p.parse_args()

    print("=== 呼吸の『水分依存 Q10』の生態系間比較（旗27）===")
    print("  各サイト健全年プールで Q10 vs θ の Spearman r＋ブートCI。r>0 かつ CI が 0 跨がず＝水分依存Q10\n")
    print(f"  {'サイト':<8} {'年':>3} {'N':>7} {'r(Q10 vs θ)':>11} {'95%CI':>18}  判定")
    rows = []
    for s in a.sites:
        try:
            df, used = load_site(s, a.month, a.qc_max, AnalysisConfig, get_site,
                                 load_raw_all, get_site_years)
        except Exception as e:
            print(f"  {s:<8} SKIP {type(e).__name__}: {e}"); continue
        if df is None or len(df) < 3000:
            print(f"  {s:<8} データ不足"); continue
        Ta = df["Ta"].to_numpy(); th = df["th"].to_numpy(); GER = df["GER"].to_numpy()
        r, ci, fr = _boot_trend(Ta, th, GER, a.nbin, 200)
        cistr = f"[{ci[0]:+.2f},{ci[1]:+.2f}]" if isinstance(ci, tuple) else "—"
        sig = isinstance(ci, tuple) and ci[0] > 0            # 水分依存 Q10（正）
        neg = isinstance(ci, tuple) and ci[1] < 0
        mark = ("✅水分依存Q10(正)" if sig else "×逆(乾で高Q10)" if neg else "△CI が0跨ぎ")
        rows.append({"site": s, "r": r, "ci": ci, "forest": s in FOREST, "n": used})
        print(f"  {s:<8} {used:>3} {len(df):>7} {r:>+11.2f} {cistr:>18}  {mark}")

    forests = [r for r in rows if r["forest"] and isinstance(r["ci"], tuple)]
    npos = sum(1 for r in forests if r["ci"][0] > 0)
    print(f"\n  森林で『水分依存Q10(r>0, CI>0)』：{npos}/{len(forests)}")
    if forests and npos == len(forests):
        print("    → ✅ 長期森林で一般的＝JP-Tak 固有でない（水分依存 Q10 は森林の普遍署名）")
    elif npos >= max(2, len(forests) - 1):
        print("    → ○ 大半の森林で一般的（例外は要精査）")
    else:
        print("    → △ 森林で割れる＝JP-Tak 突出の可能性（旗21 と同じ構図か要確認）")
    print("  留保：生プール(年間差・季節)、GER 分割派生、Q10 絶対値でなく傾向を見る。--qc-max 1 で穴埋め依存確認。")

    if a.fig and rows:
        draw(rows, a.fig); print(f"\n  [図] {a.fig}")


if __name__ == "__main__":
    main()
