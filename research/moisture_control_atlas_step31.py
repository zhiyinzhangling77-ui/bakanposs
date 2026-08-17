"""旗31：生態系ごとの「土壌水分θの炭素制御」アトラス（他生態系で比較）。

軸にする変数の関係＝**θ→炭素**。ただし θ と気温・放射は共変動する（乾＝暑い・晴れ）ので、
生の相関は温度/光の効きを拾ってしまう。so **主要駆動を差し引いた偏相関**で「θ独自の制御」を測る：
  ・**θ→GER | Ta**：温度を差し引いた上で、水分が呼吸をどれだけ動かすか
  ・**θ→GEP | Rg,VPD**：光と乾燥を差し引いた上で、水分が光合成をどれだけ動かすか
偏 Spearman（順位化してから制御変数を線形除去し残差相関）で、外れ値・非線形単調にも頑健に。

各サイト健全年の夏を日平均に集約（日周期を除き、水分が動く数日〜季節スケールを残す）→
IGBP 生態系でグループ化→「どの生態系で水分制御が強い/符号が違うか」を比較。EBR（旗29）を
品質重みとして併記。仮説：乾燥草原=強い正（水分律速）／湿潤林=弱い／湛水水田=特殊（θ不動・
嫌気で符号反転もありうる）。

    python research/moisture_control_atlas_step31.py                         # 合成で検証
    python research/moisture_control_atlas_step31.py --fig moisture_atlas.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

# 生態系を張る既定サイト（森林・草原・湿地・水田）。--sites で上書き可。
DEFAULT_SITES = ["JP-Tak", "JP-Ta2", "JP-Tef", "JP-Fhk", "JP-Fjy",   # 森林
                 "CN-HaM", "MN-Kbu",                                  # 草原
                 "JP-BBY",                                            # 湿地
                 "JP-Mse"]                                            # 水田


def _rank(a):
    import pandas as pd
    return pd.Series(np.asarray(a, float)).rank().to_numpy()


def partial_spearman(y, x, controls):
    """偏 Spearman：y,x,controls を順位化→controls を線形除去した残差同士の相関。"""
    y = np.asarray(y, float); x = np.asarray(x, float)
    ok = np.isfinite(y) & np.isfinite(x)
    Zc = []
    for c in controls:
        c = np.asarray(c, float); ok = ok & np.isfinite(c); Zc.append(c)
    if ok.sum() < 20:
        return np.nan, int(ok.sum())
    ry = _rank(y[ok]); rx = _rank(x[ok])
    if not Zc:
        return float(np.corrcoef(ry, rx)[0, 1]), int(ok.sum())
    Z = np.column_stack([_rank(c[ok]) for c in Zc] + [np.ones(ok.sum())])
    bx = np.linalg.lstsq(Z, rx, rcond=None)[0]
    by = np.linalg.lstsq(Z, ry, rcond=None)[0]
    xr = rx - Z @ bx; yr = ry - Z @ by
    if xr.std() == 0 or yr.std() == 0:
        return np.nan, int(ok.sum())
    return float(np.corrcoef(xr, yr)[0, 1]), int(ok.sum())


def _boot_ci(y, x, controls, nboot=400, seed=0):
    """日ブートで偏 Spearman の 95%CI。CI が 0 を跨がなければ符号有意。"""
    r0, n = partial_spearman(y, x, controls)
    if not np.isfinite(r0):
        return r0, None, n
    rng = np.random.default_rng(seed)
    ys = np.asarray(y, float)
    rs = []
    for _ in range(nboot):
        idx = rng.integers(0, len(ys), len(ys))
        r, _ = partial_spearman(ys[idx], np.asarray(x)[idx],
                                [np.asarray(c)[idx] for c in controls])
        if np.isfinite(r):
            rs.append(r)
    if len(rs) < 20:
        return r0, None, n
    lo, hi = np.percentile(rs, [2.5, 97.5])
    return r0, (float(lo), float(hi)), n


def daily_summer(site, months, qc_max):
    """健全年の夏を日平均に集約した DataFrame（th,Ta,Rg,VPD,GER,GEP）。"""
    import pandas as pd
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    from japanflux_pn.run_robustness import get_site_years
    cfg = AnalysisConfig(qc_max=qc_max) if qc_max is not None else AnalysisConfig()
    years, mo = get_site_years(site)
    ms = sorted(months or mo)
    raw = load_raw_all(get_site(site), cfg)
    raw = raw[raw.index.month.isin(ms)]
    keep = [c for c in ["th", "Ta", "Rg", "VPD", "GER", "GEP"] if c in raw.columns]
    raw = raw[keep]
    daily = raw.groupby(raw.index.normalize()).mean()
    return daily.dropna(), len(years)


def analyze_site(site, months, qc_max):
    d, nyr = daily_summer(site, months, qc_max)
    if len(d) < 60:
        return {"note": f"日数不足({len(d)})"}
    out = {"n_days": len(d), "n_years": nyr}
    # θ→GER | Ta（温度を差し引いた水分の呼吸制御）
    out["ger"] = _boot_ci(d["GER"].to_numpy(), d["th"].to_numpy(),
                          [d["Ta"].to_numpy()])
    # θ→GEP | Rg,VPD（光・乾燥を差し引いた水分の光合成制御）
    ctrl = [d[c].to_numpy() for c in ("Rg", "VPD") if c in d]
    out["gep"] = _boot_ci(d["GEP"].to_numpy(), d["th"].to_numpy(), ctrl)
    return out


def make_synth(kind, days=1000, seed=0):
    """水分律速 vs 水分飽和 の2生態系を合成（θ→炭素の偏相関の符号を検証）。"""
    import pandas as pd
    rng = np.random.default_rng(seed)
    Ta = 20 + rng.normal(0, 4, days)
    th = np.clip(0.3 + rng.normal(0, 0.08, days), 0.05, 0.6)
    Rg = 300 + rng.normal(0, 60, days); VPD = np.clip(1 + 0.1 * (Ta - 20), 0.1, None)
    if kind == "water_limited":       # 水分律速：θが炭素を独自に上げる
        GER = np.exp(0.06 * (Ta - 20)) * (0.5 + 1.5 * th) + rng.normal(0, 0.1, days)
        GEP = (0.01 * Rg) * (0.4 + 1.6 * th) + rng.normal(0, 0.3, days)
    else:                              # 水分飽和：θは効かない（温度/光のみ）
        GER = np.exp(0.06 * (Ta - 20)) + rng.normal(0, 0.1, days)
        GEP = 0.01 * Rg - 0.3 * VPD + rng.normal(0, 0.3, days)
    idx = pd.date_range("2001-07-01", periods=days, freq="D")
    return pd.DataFrame({"th": th, "Ta": Ta, "Rg": Rg, "VPD": VPD,
                         "GER": GER, "GEP": GEP}, index=idx)


def _fmt(res):
    if res is None or not isinstance(res, tuple):
        r, ci, n = (res if isinstance(res, tuple) else (np.nan, None, 0))
    r, ci, n = res
    if not np.isfinite(r):
        return f"{'—':>6} {'':>16}"
    cistr = f"[{ci[0]:+.2f},{ci[1]:+.2f}]" if ci else "[—]"
    sig = "✓" if (ci and (ci[0] > 0 or ci[1] < 0)) else "·"
    return f"{r:>+6.2f} {cistr:>16}{sig}"


def draw(rows, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    jpf = "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"
    jp = fm.FontProperties(fname=jpf) if Path(jpf).exists() else None
    rows = [r for r in rows if isinstance(r.get("ger"), tuple) and np.isfinite(r["ger"][0])]
    if not rows:
        return
    grp_color = {"森林": "#1f7a3d", "草原": "#b8860b", "湿地": "#3a6ea5",
                 "農地": "#a63d40", "?": "#888"}
    labels = [f"{r['site']}\n{r['type']}" for r in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(1, 2, figsize=(2.0 + 1.0 * len(rows), 4.8), sharey=True)
    for ax, key, title in [(axes[0], "ger", "θ→呼吸GER | 温度Ta"),
                           (axes[1], "gep", "θ→光合成GEP | 光Rg,乾燥VPD")]:
        for i, r in enumerate(rows):
            res = r.get(key)
            if not (isinstance(res, tuple) and np.isfinite(res[0])):
                continue
            rr, ci, _ = res
            col = grp_color.get(r["type"], "#888")
            if ci:
                ax.plot([i, i], [ci[0], ci[1]], color=col, lw=2, zorder=1)
            ax.scatter([i], [rr], c=col, s=70, zorder=3)
        ax.axhline(0, color="#999", lw=1)
        ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right",
                                             fontproperties=jp, fontsize=8)
        ax.set_title(title, fontproperties=jp)
    axes[0].set_ylabel("偏 Spearman r（水分の独自制御, ＋95%CI）", fontproperties=jp)
    fig.suptitle("生態系ごとの土壌水分θの炭素制御（旗31）\n"
                 "緑=森林 金=草原 青=湿地 赤=農地／CI が 0 を跨がねば水分制御が有意",
                 fontproperties=jp)
    fig.savefig(path, bbox_inches="tight", dpi=140); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description="生態系ごとの土壌水分の炭素制御アトラス")
    p.add_argument("--sites", nargs="+", default=DEFAULT_SITES)
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--qc-max", type=int, default=None)
    p.add_argument("--ebr", action="store_true", help="EBR(旗29)を品質列に併記")
    p.add_argument("--fig", default=None)
    a = p.parse_args()

    if a.sites == ["SYNTH"] or a.sites == ["synth"]:
        a.sites = []
    if not a.sites:
        print("=== 旗31 合成検証：θ→炭素の偏相関で水分律速 vs 飽和を分けられるか ===")
        for kind, lab in [("water_limited", "水分律速（θが炭素を独自制御）"),
                          ("saturated", "水分飽和（θは効かない）")]:
            d = make_synth(kind)
            rg = _boot_ci(d["GER"].to_numpy(), d["th"].to_numpy(), [d["Ta"].to_numpy()])
            gp = _boot_ci(d["GEP"].to_numpy(), d["th"].to_numpy(),
                          [d["Rg"].to_numpy(), d["VPD"].to_numpy()])
            print(f"  {lab:<26} θ→GER|Ta={_fmt(rg)}   θ→GEP|Rg,VPD={_fmt(gp)}")
        print("\n  → 水分律速は両方 r>0 で CI が 0 を跨がず、飽和は r≈0 が期待。")
        return

    from japanflux_pn.ecosystem import classify
    cls = classify(sites=a.sites).set_index("site")
    ebr_of = {}
    if a.ebr:
        from japanflux_pn.energy_closure import site_ebr
        from japanflux_pn.sites import get_site

    qtag = f"QC≤{a.qc_max}" if a.qc_max is not None else "gap-fill込み"
    print(f"=== 旗31 実データ 生態系ごとの土壌水分θの炭素制御（{qtag}, 月={a.month}）===")
    print("  偏 Spearman：θ→GER|Ta（温度差引後の水分の呼吸制御）／θ→GEP|Rg,VPD（光・乾燥差引後）。")
    print("  ✓=CI が 0 を跨がず符号有意。EBR は品質重み（<0.7 は不閉合＝割引）。\n")
    ehdr = f" {'EBR':>5}" if a.ebr else ""
    print(f"  {'サイト':<8} {'生態系':>6} {'年':>3} {'日数':>5}{ehdr}  "
          f"{'θ→GER|Ta':>10}{'95%CI':>16}  {'θ→GEP|Rg,VPD':>12}{'95%CI':>16}")
    rows = []
    for s in a.sites:
        typ = cls.loc[s, "type"] if s in cls.index else "?"
        try:
            res = analyze_site(s, a.month, a.qc_max)
        except Exception as e:
            print(f"  {s:<8} {typ:>6} SKIP {type(e).__name__}: {e}"); continue
        if "note" in res:
            print(f"  {s:<8} {typ:>6} {res['note']}"); continue
        rec = {"site": s, "type": typ, **res}
        estr = ""
        if a.ebr:
            try:
                e = site_ebr(get_site(s), a.month, qc_max=1)
                ev = e.get("ebr", float("nan"))
                rec["ebr"] = ev
                estr = f" {ev:>5.2f}" if ev == ev else f" {'—':>5}"
            except Exception:
                estr = f" {'—':>5}"
        rows.append(rec)
        print(f"  {s:<8} {typ:>6} {res['n_years']:>3} {res['n_days']:>5}{estr}  "
              f"{_fmt(res['ger'])}  {_fmt(res['gep'])}")

    # 生態系グループごとに「水分制御が有意なサイト数」を集計
    print("\n  === 生態系グループ別まとめ（θ→GER|Ta の符号有意サイト）===")
    from collections import defaultdict
    grp = defaultdict(list)
    for r in rows:
        if isinstance(r.get("ger"), tuple) and np.isfinite(r["ger"][0]):
            grp[r["type"]].append(r)
    for typ, rs in sorted(grp.items()):
        pos = sum(1 for r in rs if r["ger"][1] and r["ger"][1][0] > 0)
        neg = sum(1 for r in rs if r["ger"][1] and r["ger"][1][1] < 0)
        med = np.median([r["ger"][0] for r in rs])
        print(f"  {typ:<6} n={len(rs)}  正(水分↑呼吸)={pos} 負={neg} 中央r={med:+.2f}")
    print("\n  読み方：生態系で符号がそろえば『水分制御の向きは生態系署名』、割れれば固有。")
    print("  留保：夏内の季節トレンド（θ低下と炭素低下の共変動）は主要駆動の偏相関で概ね除くが")
    print("    完全でない。生プール（年間差）。GER/GEP は分割派生量。EBR<0.7 のサイトは割引。")
    if a.fig and rows:
        draw(rows, a.fig); print(f"\n  [図] {a.fig}")


if __name__ == "__main__":
    main()
