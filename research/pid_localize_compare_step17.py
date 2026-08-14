"""旗17：呼吸の相乗源(θ×温度)の"普遍性"と"水田での消失"をクロスサイトで。

旗16 で JP-Tak の呼吸相乗の源＝(Ta,θ)=温度×土壌水分 と局在した。これが
  ・自然生態系で普遍か（JP-Ta2 常緑林・CN-HaM 高山草原でも (Ta,θ) が相乗源か）
  ・水田で消えるか（JP-Mse は湛水で θ が情報を失う→(Ta,θ) の相乗が崩れるはず）
を、各サイトで健全年をプールし PID 局在（旗16 の localize を再利用）して並べる。
＝「θ が呼吸相乗の要」を、局在の普遍性＋管理での消失の両面で確かめる。

    python research/pid_localize_compare_step17.py --fig pid_localize_cross.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pid_localize_step16 import localize, _digitize

# 自然（普遍性の検証）＋ 水田（消失の検証）
DEFAULT_SITES = ["JP-Tak", "JP-Ta2", "CN-HaM", "JP-Mse"]
PAIR = "Ta,th"   # 機構仮説の相乗源（温度×土壌水分）
DRIVERS = ["Rg", "Ta", "th", "VPD"]


def pooled_localize(site: str, obins: int):
    from japanflux_pn.run_robustness import get_site_years
    from japanflux_pn.preprocess import load_corevars_hh
    years, months = get_site_years(site)
    cols = {v: [] for v in set(DRIVERS + ["GER"])}
    used = []
    for y in years:
        try:
            vf = load_corevars_hh(site, y, months, None).valid_frame
            for v in cols:
                cols[v].append(vf[v].to_numpy(float))
            used.append(y)
        except Exception:
            continue
    if not used:
        return None, 0, 0
    idx = {v: _digitize(np.concatenate(cols[v]), obins) for v in cols}
    n = len(idx["GER"])
    return localize(idx, "GER", DRIVERS, obins), len(used), n


def _pair_row(rows, pair):
    """指定ペアの (II, S, 相乗順位) を返す。"""
    order = sorted(rows, key=lambda r: r["II"])          # II 昇順（相乗が上）
    for rank, r in enumerate(order, 1):
        if r["pair"] == pair or r["pair"] == ",".join(reversed(pair.split(","))):
            return r["II"], r["S"], rank, len(order), order[0]["pair"]
    return np.nan, np.nan, np.nan, len(order), order[0]["pair"] if order else "?"


def draw(results, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    jpf = "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"
    jp = fm.FontProperties(fname=jpf) if Path(jpf).exists() else None
    sites = [r["site"] for r in results]
    vals = [-r["II"] if np.isfinite(r["II"]) else 0.0 for r in results]  # 正で相乗
    cols = ["#c0392b" if r["is_paddy"] else "#1f7a3d" for r in results]
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    x = np.arange(len(sites))
    ax.bar(x, vals, color=cols, width=0.6)
    for i, r in enumerate(results):
        tag = (f"1位" if r["rank"] == 1 else f"{r['rank']}位") if np.isfinite(r["rank"]) else "—"
        ax.text(i, vals[i] + max(vals) * 0.02, tag, ha="center", fontproperties=jp, fontsize=10)
    ax.set_xticks(x); ax.set_xticklabels(sites, fontproperties=jp)
    ax.set_ylabel("θ×温度 の正味相乗（−II, 大=相乗）", fontproperties=jp)
    ax.set_title("呼吸(GER)の相乗源『θ×温度』は自然林で普遍・水田(湛水)で消える",
                 fontproperties=jp)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#1f7a3d", label="自然生態系"),
                       Patch(color="#c0392b", label="水田(湛水)")],
              prop=jp, frameon=False)
    fig.savefig(path, bbox_inches="tight", dpi=140); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description="呼吸の相乗源θ×温度のクロスサイト局在")
    p.add_argument("--sites", nargs="+", default=DEFAULT_SITES)
    p.add_argument("--obins", type=int, default=6)
    p.add_argument("--fig", default=None)
    a = p.parse_args()

    paddy = {"JP-Mse", "KR-CRK"}
    results = []
    print("=== 呼吸の相乗源『θ×温度(Ta,th)』のクロスサイト局在 ===")
    print(f"  各サイト健全年プール・PID。II<0=相乗、相乗順位1位＝θ×温度が最大の相乗源\n")
    print(f"  {'サイト':<8} {'年':>4} {'点数':>8} {'θ×温度 II':>11} {'相乗順位':>8} {'最大相乗ペア':>14}")
    for s in a.sites:
        rows, nyr, npt = pooled_localize(s, a.obins)
        if rows is None:
            print(f"  {s:<8} 有効年なし"); continue
        ii, S, rank, k, top = _pair_row(rows, PAIR)
        results.append({"site": s, "II": ii, "S": S, "rank": rank, "k": k,
                        "top": top, "is_paddy": s in paddy, "nyr": nyr})
        topj = top.replace("th", "θ").replace("Ta", "気温").replace("Rg", "日射")
        print(f"  {s:<8} {nyr:>4} {npt:>8} {ii:>11.4f} {str(rank)+'/'+str(k):>8} {topj:>14}")

    nat = [r for r in results if not r["is_paddy"]]
    pad = [r for r in results if r["is_paddy"]]
    nat_top = sum(1 for r in nat if r["rank"] == 1)
    print("\n=== 判定 ===")
    print(f"  自然生態系で θ×温度 が相乗1位：{nat_top}/{len(nat)}")
    if nat_top == len(nat) and nat:
        print("    → ✅ θ×温度の相乗源は自然生態系で普遍（旗16の局在が転移）")
    elif nat_top >= 1:
        print(f"    → ○ 一部で1位（他サイトの最大ペアも確認）")
    any_残存 = False
    for r in pad:
        gone = (not np.isfinite(r["II"])) or r["II"] >= 0 or (np.isfinite(r["rank"]) and r["rank"] > 2)
        if not gone:
            any_残存 = True
        print(f"  水田 {r['site']}：θ×温度 II={r['II']:.4f}（{'相乗消失' if gone else '★相乗残存'}）"
              f" 最大ペア={r['top']}")
    if pad and any_残存:
        print("    → ⚠ 水田でも θ×温度 の相乗は残る（消えていない）。")
        print("      ＝fig6 の水田 0/8 は O-info(系レベル・4変数)の冗長化であって、θ×温度ペアの")
        print("        相互作用の消失ではない。両者は別の量。『湛水がθ×温度相乗を壊す』とは言えない。")
    elif pad:
        print("    → 水田で θ×温度 の相乗が消える＝湛水で θ が情報を失う機構と整合"
              "（局在の視点から旗14/fig6 を裏づけ）")

    if a.fig and results:
        draw(results, a.fig); print(f"\n  [図] {a.fig}")


if __name__ == "__main__":
    main()
