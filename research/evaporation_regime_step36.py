"""旗36：蒸発レジームの生態系アトラス（水制限 vs エネルギー制限）。

旗35 で「独立測定間のリンクが本物」と分かった。潜熱 γLE(蒸発) は**独立測定＝派生でない**ので、
θ→γLE は派生量問題ゼロの綺麗な軸。旗31（θ→炭素）と同じ枠組みで標的を蒸発にする。

問い：土壌水分 θ は蒸発 γLE を、放射 Rg を差し引いた上でどれだけ制御するか。
  ・**θ→γLE | Rg**（偏 Spearman）：放射を差し引いた水分の蒸発制御
     正で強い＝**水制限（water-limited ET）**（水が蒸発を律速, 乾燥系）
     ≈0＝**エネルギー制限（energy-limited ET）**（蒸発は放射に従う, 湿潤系）
  ・対照 corr(Rg, γLE)：蒸発がどれだけ放射に従うか（エネルギー制限の指標）
＝Budyko の水/エネルギー制限を情報理論で。旗34 の「乾燥ステップで潜熱が放射から脱結合」を生態系全体で検証。
使う変数 Rg, θ, γLE(gLE) は全部独立測定＝旗35 の"岩盤"（γH/γLE の共通 w' は留保）。

    python research/evaporation_regime_step36.py                       # 合成で検証
    python research/evaporation_regime_step36.py --qc-max 1 --ebr --block --fig evap_regime.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from moisture_control_atlas_step31 import partial_spearman, _boot_ci, DEFAULT_SITES


def daily_energy(site, months, qc_max):
    """健全年の夏を日平均に集約（th,Rg,gLE,gH,Ta）。"""
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
    keep = [c for c in ["th", "Rg", "gLE", "gH", "Ta"] if c in raw.columns]
    daily = raw[keep].groupby(raw.index.normalize()).mean()
    return daily.dropna(), len(years)


def analyze_site(site, months, qc_max, block=False):
    d, nyr = daily_energy(site, months, qc_max)
    if len(d) < 60 or "gLE" not in d:
        return {"note": f"日数不足/欠測({len(d)})"}
    yr = d.index.year.to_numpy() if block else None
    out = {"n_days": len(d), "n_years": nyr}
    # θ→γLE | Rg（放射を差し引いた水分の蒸発制御）＝水制限の指標
    out["le"] = _boot_ci(d["gLE"].to_numpy(), d["th"].to_numpy(),
                         [d["Rg"].to_numpy()], blocks=yr)
    # 対照：蒸発が放射に従う度合い（エネルギー制限）
    out["rg_le"] = float(np.corrcoef(d["Rg"], d["gLE"])[0, 1])
    # 参考：θ→γH | Rg（顕熱・逆を期待＝乾くと顕熱↑）
    if "gH" in d:
        out["h"] = _boot_ci(d["gH"].to_numpy(), d["th"].to_numpy(),
                            [d["Rg"].to_numpy()], blocks=yr)
    return out


def make_synth(kind, days=1000, seed=0):
    import pandas as pd
    rng = np.random.default_rng(seed)
    Rg = np.clip(300 + rng.normal(0, 60, days), 0, None)
    th = np.clip(0.3 + rng.normal(0, 0.08, days), 0.05, 0.6)
    if kind == "water_limited":       # 蒸発は水分律速：γLE = Rg·g(θ)
        gLE = Rg * (0.2 + 1.4 * th) + rng.normal(0, 15, days)
        gH = Rg * (0.9 - 1.2 * th) + rng.normal(0, 15, days)     # 乾くと顕熱↑
    else:                              # エネルギー制限：γLE = Rg のみ（θ無関係）
        gLE = 0.6 * Rg + rng.normal(0, 15, days)
        gH = 0.3 * Rg + rng.normal(0, 15, days)
    idx = pd.date_range("2001-07-01", periods=days, freq="D")
    return pd.DataFrame({"th": th, "Rg": Rg, "gLE": np.clip(gLE, 0, None),
                         "gH": np.clip(gH, 0, None), "Ta": 20 + rng.normal(0, 3, days)},
                        index=idx)


def _fmt(res):
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
    rows = [r for r in rows if isinstance(r.get("le"), tuple) and np.isfinite(r["le"][0])]
    if not rows:
        return
    gc = {"森林": "#1f7a3d", "草原": "#b8860b", "湿地": "#3a6ea5", "農地": "#a63d40", "?": "#888"}
    labels = [f"{r['site']}\n{r['type']}" for r in rows]
    x = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(2 + 0.9 * len(rows), 4.6))
    for i, r in enumerate(rows):
        rr, ci, _ = r["le"]; col = gc.get(r["type"], "#888")
        if ci:
            ax.plot([i, i], [ci[0], ci[1]], color=col, lw=2, zorder=1)
        ax.scatter([i], [rr], c=col, s=70, zorder=3)
    ax.axhline(0, color="#999", lw=1)
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=45, ha="right",
                                         fontproperties=jp, fontsize=8)
    ax.set_ylabel("θ→蒸発γLE | 放射Rg の偏Spearman（＋95%CI）", fontproperties=jp)
    ax.set_title("蒸発レジームの生態系アトラス（旗36）\n"
                 "正で強い＝水制限ET（乾燥系）／≈0＝エネルギー制限ET（湿潤系）",
                 fontproperties=jp)
    fig.savefig(path, bbox_inches="tight", dpi=140); plt.close(fig)


def main():
    p = argparse.ArgumentParser(description="蒸発レジーム(水制限vsエネルギー制限)の生態系アトラス")
    p.add_argument("--sites", nargs="+", default=DEFAULT_SITES)
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--qc-max", type=int, default=None)
    p.add_argument("--ebr", action="store_true")
    p.add_argument("--block", action="store_true", help="年ブロック・ブート")
    p.add_argument("--fig", default=None)
    a = p.parse_args()

    if a.sites == ["SYNTH"] or a.sites == ["synth"]:
        a.sites = []
    if not a.sites:
        print("=== 旗36 合成検証：θ→蒸発の偏相関で 水制限 vs エネルギー制限 を分けられるか ===")
        for kind, lab in [("water_limited", "水制限ET（θが蒸発を制御）"),
                          ("energy_limited", "エネルギー制限ET（蒸発は放射のみ）")]:
            d = make_synth(kind)
            le = _boot_ci(d["gLE"].to_numpy(), d["th"].to_numpy(), [d["Rg"].to_numpy()])
            rgle = np.corrcoef(d["Rg"], d["gLE"])[0, 1]
            print(f"  {lab:<28} θ→γLE|Rg={_fmt(le)}  corr(Rg,γLE)={rgle:+.2f}")
        print("\n  → 水制限は θ→γLE|Rg が正で CI が 0 を跨がず、エネルギー制限は ≈0 が期待。")
        return

    from japanflux_pn.ecosystem import classify
    cls = classify(sites=a.sites).set_index("site")
    if a.ebr:
        from japanflux_pn.energy_closure import site_ebr
        from japanflux_pn.sites import get_site

    qtag = f"QC≤{a.qc_max}" if a.qc_max is not None else "gap-fill込み"
    rtag = "・年ブロックブート" if a.block else ""
    print(f"=== 旗36 実データ 蒸発レジーム 水制限 vs エネルギー制限（{qtag}{rtag}, 月={a.month}）===")
    print("  θ→γLE|Rg（放射差引後の水分の蒸発制御）。正で強い=水制限ET / ≈0=エネルギー制限ET。")
    print("  corr(Rg,γLE)=蒸発が放射に従う度合い(エネルギー制限の指標)。✓=CI が 0 を跨がず有意。\n")
    eh = f" {'EBR':>5}" if a.ebr else ""
    print(f"  {'サイト':<8} {'生態系':>6} {'年':>3} {'日数':>5}{eh} {'Rg~γLE':>6}  "
          f"{'θ→γLE|Rg':>10}{'95%CI':>16}  判定")
    rows = []
    for s in a.sites:
        typ = cls.loc[s, "type"] if s in cls.index else "?"
        try:
            res = analyze_site(s, a.month, a.qc_max, block=a.block)
        except Exception as e:
            print(f"  {s:<8} {typ:>6} SKIP {type(e).__name__}: {e}"); continue
        if "note" in res:
            print(f"  {s:<8} {typ:>6} {res['note']}"); continue
        le = res["le"]
        wl = isinstance(le[1], tuple) and le[1][0] > 0
        verdict = "水制限ET(水↑蒸発)" if wl else \
                  ("逆" if isinstance(le[1], tuple) and le[1][1] < 0 else "エネルギー制限ET")
        estr = ""
        if a.ebr:
            try:
                ev = site_ebr(get_site(s), a.month, qc_max=1).get("ebr", float("nan"))
                estr = f" {ev:>5.2f}" if ev == ev else f" {'—':>5}"
            except Exception:
                estr = f" {'—':>5}"
        rows.append({"site": s, "type": typ, **res})
        print(f"  {s:<8} {typ:>6} {res['n_years']:>3} {res['n_days']:>5}{estr} "
              f"{res['rg_le']:>+6.2f}  {_fmt(le)}  {verdict}")

    print("\n  === 生態系グループ別（θ→γLE|Rg の符号有意）===")
    from collections import defaultdict
    grp = defaultdict(list)
    for r in rows:
        if isinstance(r.get("le"), tuple) and np.isfinite(r["le"][0]):
            grp[r["type"]].append(r)
    for typ, rs in sorted(grp.items()):
        pos = sum(1 for r in rs if r["le"][1] and r["le"][1][0] > 0)
        med = np.median([r["le"][0] for r in rs])
        print(f"  {typ:<6} n={len(rs)}  水制限(正,CI>0)={pos}  中央r={med:+.2f}")
    print("\n  読み方：乾燥草原=水制限(正)・湿潤森林=エネルギー制限(≈0) が Budyko の予想。")
    print("    旗31(呼吸の水分律速)と対になる第2の生態系署名（蒸発の水分律速）。全変数独立測定＝派生量問題なし。")
    print("  留保：γH/γLE は共通 w'（旗35）。日平均・生プール・θ深度不統一（旗33）。")
    if a.fig and rows:
        draw(rows, a.fig); print(f"\n  [図] {a.fig}")


if __name__ == "__main__":
    main()
