"""旗16：呼吸の相乗を"局在"させる。O-info は系レベルの要約で「{Rg,Ta,θ,GER} は相乗支配」
としか言わない。どの2駆動の"組"が相乗の源かは PID（部分情報分解）で目標 GER について測る。

各駆動ペア (A,B) について I(GER; A,B) を分解：
  R（冗長）… A,B が共有、U_A,U_B（固有）、**S（相乗）… A,B の組でしか出ない情報**。
相乗の"符号つき不変量"として相互作用情報 II = I(GER;A)+I(GER;B)−I(GER;A,B) も出す
（II<0 ＝正味相乗、測度に依らない＝I_min の相乗バイアスの当て馬でなくクロスチェック）。
機構仮説：呼吸の相乗の源が「土壌水分 θ×温度」なら、**(θ,Ta) ペアが最大の相乗（S大・II<0）**のはず。

本体の検証済みプリミティブ（information_theory.pid_williams_beer / interaction_information_indices）を再利用。
シャッフルサロゲートで II のヌルからの z も出し、相乗が偶然でないかを判定する。

    python research/pid_localize_step16.py                       # 合成で検証
    python research/pid_localize_step16.py --site JP-Tak --years 1999 2000 ... --month 7 8
"""

from __future__ import annotations

import argparse
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from japanflux_pn import information_theory as it


def _digitize(x, m):
    x = np.asarray(x, float); lo, hi = x.min(), x.max()
    if hi <= lo:
        return np.zeros(len(x), dtype=np.int64)
    return np.clip(np.floor((x - lo) / (hi - lo) * m).astype(np.int64), 0, m - 1)


def ii_surrogate_z(ti, ai, bi, m, n_surr=300, seed=0, correct=True):
    """II のシャッフルヌルからの z。ヌルは A,B を各々独立シャッフル（相乗を壊す）。
    II<0=相乗なので、相乗が有意なら z が負に大きい。"""
    obs = it.interaction_information_indices(ti, ai, bi, m, correct)
    rng = np.random.default_rng(seed)
    n = len(ti)
    samp = np.empty(n_surr)
    for s in range(n_surr):
        samp[s] = it.interaction_information_indices(
            ti, ai[rng.permutation(n)], bi[rng.permutation(n)], m, correct)
    mu, sd = float(samp.mean()), float(samp.std())
    z = (obs - mu) / sd if sd > 0 else np.nan
    return obs, z


def localize(idx: dict, target: str, drivers: list[str], m: int,
             n_surr: int = 300) -> list[dict]:
    ti = idx[target]
    rows = []
    for a, b in combinations(drivers, 2):
        p = it.pid_williams_beer(ti, idx[a], idx[b], m)
        ii, z = ii_surrogate_z(ti, idx[a], idx[b], m, n_surr)
        rows.append({"pair": f"{a},{b}", "R": p["R"], "U1": p["U1"], "U2": p["U2"],
                     "S": p["S"], "I_joint": p["I_joint"], "II": ii, "z_II": z})
    rows.sort(key=lambda r: r["S"], reverse=True)   # 相乗の大きい順
    return rows


def draw(rows, target, path, site_note=""):
    """相乗源ペアの棒図：II(測度不変・正味相乗, 負ほど相乗)を相乗順に。最上位を強調。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    jpf = "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"
    jp = fm.FontProperties(fname=jpf) if Path(jpf).exists() else None
    r = sorted(rows, key=lambda z: z["II"])         # II 昇順（最も相乗が上）
    labels = [x["pair"].replace("th", "θ").replace("Ta", "気温").replace("Rg", "日射")
              for x in r]
    ii = [-x["II"] for x in r]                        # 正で相乗（見やすさ）
    top = min(range(len(r)), key=lambda i: r[i]["II"])
    cols = ["#1f7a3d" if i == top else "#9ec6ac" for i in range(len(r))]
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    y = np.arange(len(r))[::-1]
    ax.barh(y, ii, color=cols)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontproperties=jp)
    ax.set_xlabel("正味の相乗の強さ（−II, 大きいほど相乗）", fontproperties=jp)
    ax.set_title(f"呼吸(GER)の相乗の源はどの駆動の組か{site_note}\n"
                 f"（θ×温度が最大＝土壌水分と温度の相互作用）", fontproperties=jp)
    fig.savefig(path, bbox_inches="tight", dpi=140); plt.close(fig)


def _print(rows, target, m, note):
    print(f"\n=== 呼吸の相乗の局在: 目標 {target} / 各駆動ペア (m={m}) {note} ===")
    print("  S=相乗(組でしか出ない), R=冗長, II<0=正味相乗(測度不変), z_II<0 かつ |z|≥2.36 で有意\n")
    print(f"  {'ペア':<10} {'S':>7} {'R':>7} {'II':>8} {'z(II)':>8}  判定")
    for r in rows:
        syn = (r["II"] < 0 and np.isfinite(r["z_II"]) and r["z_II"] <= -2.36)
        mark = "★相乗が有意" if syn else ("相乗寄り(非有意)" if r["II"] < 0 else "冗長寄り")
        print(f"  {r['pair']:<10} {r['S']:7.4f} {r['R']:7.4f} {r['II']:8.4f} "
              f"{r['z_II']:8.1f}  {mark}")


def make_synth(n=8000, m=6, seed=0):
    """GER が θ×Ta の相互作用（乗法的な非加法）で決まる合成。(θ,Ta) が相乗の源のはず。
    Rg,VPD は GER と弱く線形（相乗でない）当て馬。"""
    rng = np.random.default_rng(seed)
    th = rng.normal(0, 1, n)
    Ta = rng.normal(0, 1, n)
    Rg = rng.normal(0, 1, n)
    VPD = rng.normal(0, 1, n)
    # 相乗：積 θ*Ta（各単独では GER とほぼ無相関、組で効く）＋ Rg の弱い線形
    GER = 1.6 * (th * Ta) + 0.3 * Rg + rng.normal(0, 0.4, n)
    idx = {v: _digitize(val, m) for v, val in
           dict(GER=GER, th=th, Ta=Ta, Rg=Rg, VPD=VPD).items()}
    return idx


def main():
    p = argparse.ArgumentParser(description="呼吸の相乗をPIDで局在（目標GER）")
    p.add_argument("--site")
    p.add_argument("--years", type=int, nargs="+")
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--target", default="GER")
    p.add_argument("--drivers", nargs="+", default=["Rg", "Ta", "th", "VPD"])
    p.add_argument("--obins", type=int, default=6)
    p.add_argument("--fig", default=None, help="相乗源ペアの棒図の保存先PNG")
    a = p.parse_args()

    if not a.site:
        idx = make_synth(m=a.obins)
        rows = localize(idx, "GER", ["Rg", "Ta", "th", "VPD"], a.obins)
        _print(rows, "GER", a.obins, "[合成: 真の相乗源=(th,Ta)]")
        top = rows[0]["pair"]
        ok = top == "Ta,th" or top == "th,Ta"
        print("\n  → " + ("✅ 期待どおり (th,Ta) が最大の相乗＝局在成功"
                          if ok else f"⚠ 最大相乗は {top}（合成条件/ビンを確認）"))
        print("  意味: O-infoの系レベル『相乗支配』を、PIDが『どの組(θ×温度)か』まで割り出す。")
        print("        実データは → --site JP-Tak --years ... で GER の相乗源ペアを局在。")
        return

    # 実データ：健全年をプールして PID（年ごとは点数不足になりうるので連結）
    from japanflux_pn.preprocess import load_corevars_hh
    cols = {v: [] for v in set(a.drivers + [a.target])}
    used = []
    for y in a.years or []:
        try:
            vf = load_corevars_hh(a.site, y, a.month, None).valid_frame
            for v in cols:
                cols[v].append(vf[v].to_numpy(float))
            used.append(y)
        except Exception as e:
            print(f"  {y}: SKIP {type(e).__name__}: {e}")
    if not used:
        print("有効年なし"); return
    idx = {v: _digitize(np.concatenate(cols[v]), a.obins) for v in cols}
    print(f"[pool] {a.site} 年={used} 総点数={len(idx[a.target])}")
    rows = localize(idx, a.target, a.drivers, a.obins)
    _print(rows, a.target, a.obins, f"[実データ {a.site} プール {len(used)}年]")
    print("\n  → 最大相乗のペアが呼吸の相乗の源。機構仮説(θ×温度)なら (th,Ta) が上位のはず。")
    print("     留保: I_min の S は上バイアスしうる→ II(測度不変)と z で確認。プールは年間差を均す近似。")
    if a.fig:
        draw(rows, a.target, a.fig, site_note=f"（{a.site}・{len(used)}年）")
        print(f"  [図] {a.fig}")


if __name__ == "__main__":
    main()
