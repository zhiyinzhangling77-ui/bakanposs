"""旗20：標準的な呼吸モデルの"掛け算仮定"は破れているか？を直接テストする。

多くの Reco モデルは R = R_ref · f(Ta) · g(θ)（温度関数×水分関数の**分離可能＝掛け算**）。
この研究の主張は「モデルは高次（非加法）相互作用を取りこぼす」。旗16-19 で情報の相乗は示したが、
「情報の相乗」と「掛け算モデルで再現できない」は別物。ここでは**GER の応答曲面 E[GER|Ta,θ] が
掛け算（分離可能）からどれだけズレるか**を直接測る＝モデル形の妥当性の観測的検証。

やること：生の GER を Ta×θ の格子（既定 6×6）に入れ、各セルの平均 GER で応答曲面 M[i,j] を作る。
log 空間で加法モデル（=元空間で掛け算）L_hat[i,j]=μ+r[i]+c[j] を最小二乗で当て、
**交互作用（残差）の割合 = SS_interaction / SS_total** を出す。大きいほど掛け算では書けない
＝温度と水分が非加法に絡む（旗16 の θ×温度 相乗の"モデル形での現れ"）。

    python research/interaction_surface_step20.py                        # 合成で検証
    python research/interaction_surface_step20.py --site JP-Tak --years 1999 ... --month 7 8
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def _bin(x, nb):
    x = np.asarray(x, float); lo, hi = np.nanmin(x), np.nanmax(x)
    if hi <= lo:
        return np.zeros(len(x), dtype=int)
    return np.clip(((x - lo) / (hi - lo) * nb).astype(int), 0, nb - 1)


def interaction_fraction(Ta, th, GER, nb=6, min_cell=20):
    """応答曲面 E[GER|Ta,θ] の掛け算（分離可能）からのズレ＝交互作用の割合を返す。

    log 空間で 2元加法モデル μ+r[i]+c[j] を当て（=元空間の掛け算 a[i]*b[j]）、
    交互作用 SS / 全 SS を返す。0=完全に掛け算で書ける、大=非加法な絡みがある。
    """
    Ta, th, GER = map(lambda a: np.asarray(a, float), (Ta, th, GER))
    ok = np.isfinite(Ta) & np.isfinite(th) & np.isfinite(GER) & (GER > 0)
    Ta, th, GER = Ta[ok], th[ok], GER[ok]
    bi, bj = _bin(Ta, nb), _bin(th, nb)
    M = np.full((nb, nb), np.nan); N = np.zeros((nb, nb))
    for i in range(nb):
        for j in range(nb):
            m = (bi == i) & (bj == j)
            if m.sum() >= min_cell:
                M[i, j] = GER[m].mean(); N[i, j] = m.sum()
    mask = np.isfinite(M)
    if mask.sum() < nb + 2:            # セルが少なすぎ
        return np.nan, M, N, mask
    L = np.log(M)
    # 重み付き 2元加法当てはめ（欠セルは無視）。反復で行・列効果を推定。
    w = np.where(mask, N, 0.0)
    mu = np.sum(w * np.where(mask, L, 0)) / w.sum()
    r = np.zeros(nb); c = np.zeros(nb)
    for _ in range(200):
        for i in range(nb):
            wi = w[i]; s = wi.sum()
            r[i] = (np.sum(wi * (np.where(mask[i], L[i], 0) - mu - c)) / s) if s > 0 else 0
        for j in range(nb):
            wj = w[:, j]; s = wj.sum()
            c[j] = (np.sum(wj * (np.where(mask[:, j], L[:, j], 0) - mu - r)) / s) if s > 0 else 0
    L_hat = mu + r[:, None] + c[None, :]
    resid = np.where(mask, L - L_hat, 0.0)
    tot = np.where(mask, L - mu, 0.0)
    ss_int = float(np.sum(w * resid ** 2))
    ss_tot = float(np.sum(w * tot ** 2))
    frac = ss_int / ss_tot if ss_tot > 0 else np.nan
    return frac, M, N, mask


def surrogate_pvalue(Ta, th, GER, nb, min_cell, nperm=500, seed=0):
    """Taビン内でθをシャッフル（Ta主効果・θ周辺分布は保存、Ta×θ依存だけ壊す）帰無から
    交互作用割合の p 値を出す。観測がヌルを超えれば非加法は有意。"""
    Ta, th, GER = map(lambda a: np.asarray(a, float), (Ta, th, GER))
    ok = np.isfinite(Ta) & np.isfinite(th) & np.isfinite(GER) & (GER > 0)
    Ta, th, GER = Ta[ok], th[ok], GER[ok]
    obs = interaction_fraction(Ta, th, GER, nb, min_cell)[0]
    if not np.isfinite(obs):
        return obs, np.nan, np.nan
    bi = _bin(Ta, nb)
    groups = [np.where(bi == i)[0] for i in range(nb)]
    rng = np.random.default_rng(seed)
    null = np.empty(nperm)
    for s in range(nperm):
        th_s = th.copy()
        for g in groups:
            if g.size > 1:
                th_s[g] = th[g[rng.permutation(g.size)]]
        null[s] = interaction_fraction(Ta, th_s, GER, nb, min_cell)[0]
    null = null[np.isfinite(null)]
    p = (np.sum(null >= obs) + 1) / (null.size + 1)
    return obs, float(np.mean(null)), float(p)


def make_synth(kind, n=40000, seed=0):
    rng = np.random.default_rng(seed)
    Ta = rng.uniform(10, 30, n); th = rng.uniform(0.1, 0.5, n)
    fTa = np.exp(0.06 * (Ta - 20))          # Q10 的
    gth = th / (0.2 + th)                    # 水分飽和関数
    if kind == "mult":                        # 純・掛け算（交互作用ゼロを期待）
        GER = 2.0 * fTa * gth * (1 + rng.normal(0, 0.05, n))
    elif kind == "interact":                  # 非加法：水分が温度感度を変える
        GER = 2.0 * fTa ** (gth) * gth * (1 + rng.normal(0, 0.05, n))
    GER = np.clip(GER, 1e-3, None)
    return Ta, th, GER


def main():
    p = argparse.ArgumentParser(description="呼吸の掛け算モデル仮定の破れを測る")
    p.add_argument("--site")
    p.add_argument("--years", type=int, nargs="+")
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--nbins", type=int, default=6)
    p.add_argument("--min-cell", type=int, default=20)
    p.add_argument("--deyear", action="store_true",
                   help="各年の幾何平均でGERを割り年間レベル差を除く（プーリング交絡の制御）")
    p.add_argument("--nperm", type=int, default=0,
                   help="Taビン内θシャッフルの順列検定回数（0=しない, 例500）")
    a = p.parse_args()

    if not a.site:
        print("=== 旗20 合成検証：応答曲面の掛け算からのズレ（交互作用割合）===")
        for kind, label in [("mult", "純・掛け算 R=f(Ta)g(θ)"),
                            ("interact", "非加法 水分が温度感度を変える")]:
            Ta, th, GER = make_synth(kind)
            frac, *_ = interaction_fraction(Ta, th, GER, a.nbins, a.min_cell)
            print(f"  {label:<28} 交互作用割合 = {frac:6.3f}")
        print("  → 掛け算合成≈0・非加法合成が大きく出れば検出成功（掛け算モデルの破れを測れる）。")
        print("     実データは → --site JP-Tak --years <全健全年> --month 7 8")
        return

    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    import pandas as pd
    cfg = AnalysisConfig()
    raw_all = load_raw_all(get_site(a.site), cfg)
    ms = sorted(a.month)
    cols = {"GER": [], "Ta": [], "th": []}
    used = []
    for y in a.years or []:
        start = pd.Timestamp(year=y, month=ms[0], day=1)
        end = pd.Timestamp(year=y, month=ms[-1], day=1) + pd.offsets.MonthBegin(1)
        r = raw_all[(raw_all.index >= start) & (raw_all.index < end)]
        if r.empty:
            continue
        g = r["GER"].to_numpy(float)
        if a.deyear:                      # 年ごとに幾何平均で割る（log 空間で年中心化）
            gp = g[np.isfinite(g) & (g > 0)]
            gm = np.exp(np.mean(np.log(gp))) if gp.size else 1.0
            g = g / gm
        cols["GER"].append(g)
        cols["Ta"].append(r["Ta"].to_numpy(float))
        cols["th"].append(r["th"].to_numpy(float))
        used.append(y)
    if not used:
        print("有効年なし"); return
    Ta = np.concatenate(cols["Ta"]); th = np.concatenate(cols["th"]); GER = np.concatenate(cols["GER"])
    frac, M, N, mask = interaction_fraction(Ta, th, GER, a.nbins, a.min_cell)
    dy = "・年レベル除去(deyear)" if a.deyear else ""
    print(f"=== 旗20 実データ {a.site}（生 GER・Ta・θ プール {len(used)}年{dy}, {a.nbins}×{a.nbins}格子）===")
    print(f"  有効セル {int(mask.sum())}/{a.nbins*a.nbins}、総点数 {np.isfinite(GER).sum()}")
    print(f"\n  ★掛け算モデルからの交互作用割合 = {frac:.3f}")
    if a.nperm > 0:
        _, null_mu, pval = surrogate_pvalue(Ta, th, GER, a.nbins, a.min_cell, a.nperm)
        sig = "有意（非加法は偶然でない）" if (np.isfinite(pval) and pval < 0.05) else "非有意"
        print(f"  順列検定（Taビン内θシャッフル {a.nperm}回）: ヌル平均={null_mu:.3f}  p={pval:.3f}  → {sig}")
    if np.isfinite(frac):
        if frac >= 0.15:
            v = "大＝温度と水分が非加法に絡む＝掛け算モデル R=f(Ta)g(θ) では書けない（旗16 相乗のモデル形での現れ）"
        elif frac >= 0.05:
            v = "中＝掛け算からの無視できないズレ（要・年/サイト積み増し）"
        else:
            v = "小＝ほぼ掛け算で書ける（この格子・期間では非加法性は弱い）"
        print(f"  → {v}")
    print("  留保：生データのプールは年間差・季節を含む。セル平均の応答曲面なので同時分布の高次とは別角度。")
    print("       これは『モデル形（掛け算）の妥当性』の直接検証＝旗16(情報の相乗)と相補的。")


if __name__ == "__main__":
    main()
