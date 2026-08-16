"""旗26：呼吸の非加法(Ta×θ)の"中身"＝温度感度Q10は土壌水分で変わるか。

旗20/22 で呼吸は掛け算モデル R=f(Ta)g(θ) を破る(Ta×θ 非加法)と分かった。その"中身"を開ける：
生態学の定番の問い「**呼吸の温度感度 Q10 は土壌水分で変わるか**」。もし変わるなら、それが
非加法の正体（温度と水分が独立に効くのでなく、水分が温度応答を組み替える）。

やること：土壌水分 θ を分位ビンに分け、各ビン内で ln(GER)=a+b·Ta を回帰→傾き b、Q10=exp(10b)。
Q10 が θ で系統的に変われば「水分依存の Q10」＝非加法の中身。ブートストラップで傾向の有意性、
Q10 vs 水分の Spearman も出す。＝PRESENTATION §1 背景「土壌水分で Q10 が変わる」の実データ検証。

    python research/respiration_q10_moisture_step26.py                    # 合成で検証
    python research/respiration_q10_moisture_step26.py --site JP-Tak --years 1999 ... --month 7 8 --fig q10.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def q10_by_moisture(Ta, th, GER, nbin=5, min_n=200):
    """θ 分位ビンごとに Q10=exp(10·d ln GER/d Ta) を出す。"""
    Ta, th, GER = map(lambda a: np.asarray(a, float), (Ta, th, GER))
    ok = np.isfinite(Ta) & np.isfinite(th) & np.isfinite(GER) & (GER > 0)
    Ta, th, GER = Ta[ok], th[ok], GER[ok]
    edges = np.quantile(th, np.linspace(0, 1, nbin + 1))
    rows = []
    for i in range(nbin):
        lo, hi = edges[i], edges[i + 1]
        m = (th >= lo) & (th <= hi) if i == nbin - 1 else (th >= lo) & (th < hi)
        if m.sum() < min_n:
            continue
        A = np.column_stack([Ta[m], np.ones(m.sum())])
        b, a = np.linalg.lstsq(A, np.log(GER[m]), rcond=None)[0]
        rows.append({"th_center": float(np.median(th[m])), "b": float(b),
                     "Q10": float(np.exp(10 * b)), "n": int(m.sum()),
                     "GER_mean": float(GER[m].mean())})
    return rows


def _boot_trend(Ta, th, GER, nbin, min_n, nboot=500, seed=0):
    """Q10 vs θ の Spearman r と、ブートストラップ CI・p(符号の安定性)。"""
    def spear(rows):
        if len(rows) < 3:
            return np.nan
        x = np.array([r["th_center"] for r in rows]); y = np.array([r["Q10"] for r in rows])
        import pandas as pd
        return float(np.corrcoef(pd.Series(x).rank(), pd.Series(y).rank())[0, 1])
    r_obs = spear(q10_by_moisture(Ta, th, GER, nbin, min_n))
    rng = np.random.default_rng(seed)
    n = len(Ta); rs = []
    for _ in range(nboot):
        idx = rng.integers(0, n, n)
        rs.append(spear(q10_by_moisture(Ta[idx], th[idx], GER[idx], nbin, min_n)))
    rs = np.array([x for x in rs if np.isfinite(x)])
    if rs.size < 10:
        return r_obs, np.nan, np.nan
    lo, hi = np.percentile(rs, [2.5, 97.5])
    # 符号の安定性（0 を跨がないか）
    frac_same = np.mean(np.sign(rs) == np.sign(r_obs)) if np.isfinite(r_obs) else np.nan
    return r_obs, (lo, hi), frac_same


def make_synth(kind, n=60000, seed=0):
    rng = np.random.default_rng(seed)
    Ta = rng.uniform(8, 30, n); th = rng.uniform(0.1, 0.5, n)
    thn = (th - 0.3) / 0.2
    if kind == "const":       # Q10 一定（水分非依存）
        b = 0.07 + 0 * thn
    elif kind == "wet_up":    # 湿るほど Q10 大（水分が温度感度を上げる）
        b = 0.05 + 0.04 * thn
    GER = np.exp(b * (Ta - 20) + 0.5 * thn) * (1 + rng.normal(0, 0.05, n))
    return Ta, th, np.clip(GER, 1e-3, None)


def draw(rows, path, site, r, ci):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    jpf = "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"
    jp = fm.FontProperties(fname=jpf) if Path(jpf).exists() else None
    x = [r_["th_center"] for r_ in rows]; y = [r_["Q10"] for r_ in rows]
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.plot(x, y, "o-", color="#1f7a3d", ms=9, lw=2)
    for r_ in rows:
        ax.annotate(f"n={r_['n']}", (r_["th_center"], r_["Q10"]), fontsize=8,
                    color="#666", xytext=(0, 8), textcoords="offset points", ha="center")
    ax.axhline(1, color="#bbb", ls=":", lw=1)
    ax.set_xlabel("土壌水分 θ（ビン中央値）", fontproperties=jp)
    ax.set_ylabel("呼吸の温度感度 Q10", fontproperties=jp)
    cis = f"（Spearman r={r:+.2f}）" if np.isfinite(r) else ""
    ax.set_title(f"呼吸の温度感度 Q10 は土壌水分で変わるか（{site}）{cis}", fontproperties=jp)
    fig.savefig(path, bbox_inches="tight", dpi=140); plt.close(fig)


def _report(rows, r, ci, frac, tag):
    print(f"\n  === {tag} ===")
    print(f"  {'θ中央':>7} {'Q10':>7} {'GER平均':>8} {'n':>7}")
    for r_ in rows:
        print(f"  {r_['th_center']:>7.3f} {r_['Q10']:>7.2f} {r_['GER_mean']:>8.2f} {r_['n']:>7}")
    if np.isfinite(r):
        cis = f"95%CI[{ci[0]:+.2f},{ci[1]:+.2f}]" if isinstance(ci, tuple) else ""
        stable = "✅ 符号安定＝水分依存の Q10 は本物" if (isinstance(ci, tuple) and (ci[0] > 0 or ci[1] < 0)) \
            else "△ CI が 0 を跨ぐ＝傾向はあるが確証弱い"
        print(f"  → Q10 vs 土壌水分: Spearman r={r:+.2f} {cis}  {stable}")


def main():
    p = argparse.ArgumentParser(description="呼吸の温度感度Q10は土壌水分で変わるか")
    p.add_argument("--site")
    p.add_argument("--years", type=int, nargs="+")
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--nbin", type=int, default=5)
    p.add_argument("--qc-max", type=int, default=None)
    p.add_argument("--fig", default=None)
    a = p.parse_args()

    if not a.site:
        print("=== 旗26 合成検証：Q10 は土壌水分で変わるか ===")
        for kind, lab in [("const", "Q10 一定（水分非依存）"), ("wet_up", "湿るほど Q10 大")]:
            Ta, th, GER = make_synth(kind)
            rows = q10_by_moisture(Ta, th, GER, a.nbin)
            r, ci, fr = _boot_trend(Ta, th, GER, a.nbin, 200, nboot=300)
            _report(rows, r, ci, fr, lab)
        print("\n  → 一定は r≈0、湿るほど大は r>0 で CI が 0 を跨がない、が期待。")
        return

    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    import pandas as pd
    cfg = AnalysisConfig(qc_max=a.qc_max) if a.qc_max is not None else AnalysisConfig()
    raw_all = load_raw_all(get_site(a.site), cfg)
    ms = sorted(a.month)
    parts = []
    used = []
    for y in a.years or []:
        start = pd.Timestamp(year=y, month=ms[0], day=1)
        end = pd.Timestamp(year=y, month=ms[-1], day=1) + pd.offsets.MonthBegin(1)
        r = raw_all[(raw_all.index >= start) & (raw_all.index < end)]
        if not r.empty:
            parts.append(r[["GER", "Ta", "th"]]); used.append(y)
    if not parts:
        print("有効年なし"); return
    df = pd.concat(parts).dropna()
    qtag = f"・QC≤{a.qc_max}" if a.qc_max is not None else "・gap-fill込み"
    Ta = df["Ta"].to_numpy(); th = df["th"].to_numpy(); GER = df["GER"].to_numpy()
    rows = q10_by_moisture(Ta, th, GER, a.nbin)
    r, ci, fr = _boot_trend(Ta, th, GER, a.nbin, 200)
    print(f"=== 旗26 実データ {a.site}（呼吸の Q10 vs 土壌水分, {len(used)}年{qtag}, N={len(df)}）===")
    _report(rows, r, ci, fr, f"{a.site} 呼吸 Q10")
    print("\n  読み方：Q10 が θ で系統的に変われば＝非加法(Ta×θ)の中身は『水分依存の温度感度』。")
    print("     ＝掛け算モデル R=f(Ta)g(θ) が不足する具体的理由（温度応答自体が水分で組み変わる）。")
    print("  留保：生プール（年間差・季節含む）、GER は分割派生量。")
    if a.fig and rows:
        draw(rows, a.fig, a.site, r, ci); print(f"  [図] {a.fig}")


if __name__ == "__main__":
    main()
