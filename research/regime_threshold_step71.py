"""旗71：**何が水制限とエネルギー制限を分けるのか**——旗70 が残した問い。（v3）

旗70 で、**同じモンゴル半乾燥ステップ・同じ草原・113km 離れた MN-Skt** で Bowen 反転が
再現しなかった。＝問いは「乾燥地で一般に成り立つか」から**何が両者を分けるのか**に変わった。

## 仮説（**事後に立てた＝探索的**）

**θ 閾値仮説**：蒸発の水分律速は **θ がある水準を下回るときだけ**現れる。MN-Skt は
律速域に入る日が少ないため、夏全体で均すとエネルギー制限に見える。
＝**Seneviratne et al. (2010, Earth-Sci Rev) の soil moisture regime と同型で、新説ではない**。
確かめるのは**自分のデータがその枠組みと整合するか**である。

## 設計（サイト間比較を避ける）

**θ の深度がサイト間で不統一**（旗33）＝**絶対値を横に比べない**。
  ・サイト水準の乾湿は**降水量 P** で比べる（単位共通・深度問題なし）。
  ・閾値の検定は**各サイト内部**で、θ を**そのサイトの百分位**で三分位に切って行う。
  ・副次的利点：検定が**サイト内の日**で行われるため、**独立地点 n の少なさに依存しない**。

## 自分の道具の欠陥（v1・v2 で見つけた2件）と、その修正

**【10件目】プラセボが帰無になっていなかった。**
v1 は `np.roll(th, len//2)` の単一シフト。本解析は**夏（7–8月）だけを連結した系列**なので、
年あたり日数の整数倍だけずらすと**7月が7月に写り季節進行が壊れない**。実際 MN-Hst
（248日/6年≒62日/年、シフト124＝ちょうど2年分）で**プラセボが実測とほぼ同じ反転**を出した。
v2 は複数シフト化したが**不十分**だった——`n·k/(N+1)` は k が大きいと **n に近いシフト＝実質わずかな
ずれ**になる（MN-Nkh：310日でシフト286＝**実質 −24日**）。**θ は自分自身と強く相関したまま**である。
→ **v3 は循環シフトを捨て、「年の入れ替え」に変える**：ある年の θ を**別の年の**フラックスに当てる。
   同じ**日付**（夏の中の位置）どうしを対応させるので、**θ の分布・夏内構造・自己相関は保たれ**、
   **日々の対応と年の対応だけが壊れる**。＝標準的な帰無の作り方。

**【11件目】ブロック・ブートを渡し忘れていた。**
`_boot_ci` は `blocks=` を渡さないと**日ブート**になり、docstring 自身が
「**自己相関で CI が過小になる**」と警告している。旗36 は `blocks=年` を渡していたのに、
**旗71 v1/v2 は渡していなかった**＝**すべての ✓（有意）が過大**だった。
→ v3 は**層の中の年ラベルで年ブロック・ブート**を行う。

## 合成検証（3通り）

  ・`threshold`：θ<θ* でだけ水律速 → **乾の層でだけ反転**が出るはず。
  ・`none`：閾値なし → **どの層も似た弱さ**のはず。
  ・`season`：**θ とフラックスが共通の夏内トレンドだけで結ばれる**（因果なし）
    → **日付を制御すると消え、プラセボにも同じ模様が出る**はず＝交絡を捕まえられるかの試験。

    python research/regime_threshold_step71.py --synth
    python research/regime_threshold_step71.py --sites MN-Skt MN-Hst MN-Nkh MN-Kbu --month 7 8
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from moisture_control_atlas_step31 import _boot_ci
from evaporation_regime_step36 import daily_energy

LABELS = ["乾(下位1/3)", "中(中位1/3)", "湿(上位1/3)"]


def _controls(sub, with_doy):
    """偏相関の制御変数。`with_doy` なら**夏の中の位置**も制御して季節トレンドを除く。"""
    ctrl = [sub["Rg"].to_numpy()]
    if with_doy:
        ctrl.append(sub.index.dayofyear.to_numpy().astype(float))
    return ctrl


def strata_stats(d, thcol="th", with_doy=False):
    """θ の**サイト内百分位**で三分位に切り、各層で偏相関を測る（**年ブロック・ブート**）。"""
    th = d[thcol].to_numpy()
    ok = np.isfinite(th)
    if ok.sum() < 30:
        return None
    q1, q2 = np.nanpercentile(th[ok], [100 / 3, 200 / 3])
    out = []
    for lo, hi in [(-np.inf, q1), (q1, q2), (q2, np.inf)]:
        m = np.isfinite(th) & (th > lo) & (th <= hi)
        sub = d.loc[m]
        if len(sub) < 20 or "gLE" not in sub:
            out.append(None); continue
        ctrl = _controls(sub, with_doy)
        # **年ブロック・ブート**（v1/v2 は渡し忘れ＝日ブートで CI が過小だった＝欠陥11件目）
        yr = sub.index.year.to_numpy()
        le = _boot_ci(sub["gLE"].to_numpy(), sub[thcol].to_numpy(), ctrl, blocks=yr)
        h = (_boot_ci(sub["gH"].to_numpy(), sub[thcol].to_numpy(), ctrl, blocks=yr)
             if "gH" in sub else (np.nan, None, 0))
        out.append({"n": int(len(sub)), "sd": float(np.nanstd(sub[thcol])),
                    "mean": float(np.nanmean(sub[thcol])),
                    "nyr": int(len(np.unique(yr))), "le": le, "h": h})
    return out


def year_swap(d, k, thcol="th"):
    """**年の入れ替え**プラセボ：ある年の θ を、**k 年ずらした別の年**の同じ日付に当てる。

    θ の分布・夏内構造・自己相関は保たれ、**日々の対応と年の対応だけが壊れる**。
    循環シフト（v1/v2）と違い、**実質わずかなずれ**になってしまう危険がない。
    """
    years = np.unique(d.index.year)
    if len(years) < 2:
        return None
    pos = {y: i for i, y in enumerate(years)}
    key = list(zip(d.index.year.to_numpy(), d.index.dayofyear.to_numpy()))
    lookup = dict(zip(key, d[thcol].to_numpy()))
    new = [lookup.get((years[(pos[y] + k) % len(years)], doy), np.nan) for y, doy in key]
    d2 = d.copy()
    d2[thcol] = np.asarray(new, float)
    return d2


def placebo_strata(d, thcol="th", with_doy=False):
    """全ての年ずらし（1..nyear-1）を試し、|r| が**最大**になった回を基準線に採る（最悪値）。"""
    years = np.unique(d.index.year)
    best = None
    for k in range(1, len(years)):
        d2 = year_swap(d, k, thcol)
        if d2 is None:
            continue
        rows = strata_stats(d2, thcol, with_doy)
        if rows is None:
            continue
        worst = max([abs(r["le"][0]) for r in rows if r and np.isfinite(r["le"][0])],
                    default=np.nan)
        if best is None or (np.isfinite(worst) and worst > best[0]):
            best = (worst, k, rows)
    return (best[2], best[1]) if best else (None, None)


def _f(ci):
    r, c, _ = ci
    if not np.isfinite(r):
        return f"{'—':>21}"
    if not isinstance(c, tuple):
        return f"{r:+.2f}{'':>15}"
    sig = "✓" if (c[0] > 0 or c[1] < 0) else "·"
    return f"{r:+.2f} [{c[0]:+.2f},{c[1]:+.2f}]{sig}"


def show(rows):
    print(f"    {'層':<12}{'n':>5}{'年':>3}{'θ平均':>8}{'θのsd':>7}  "
          f"{'θ→γLE|Rg':<23}{'θ→γH|Rg':<23}判定")
    for lab, r in zip(LABELS, rows or []):
        if r is None:
            print(f"    {lab:<12}{'—':>5}{'—':>3}{'—':>8}{'—':>7}  日数不足"); continue
        le_pos = isinstance(r["le"][1], tuple) and r["le"][1][0] > 0
        h_neg = isinstance(r["h"][1], tuple) and r["h"][1][1] < 0
        v = ("**Bowen反転**" if (le_pos and h_neg) else
             "水制限(蒸発のみ)" if le_pos else "反転なし")
        print(f"    {lab:<12}{r['n']:>5}{r['nyr']:>3}{r['mean']:>8.3f}{r['sd']:>7.3f}  "
              f"{_f(r['le']):<23}{_f(r['h']):<23}{v}")


def make_synth(kind, nyear=6, seed=0):
    """**夏だけ・複数年**の合成（実データと同じ構造）。

    θ には**夏内の乾燥トレンド・年ごとの水準差・自己相関**を入れる。
    `season` は θ とフラックスを**共通の季節トレンドだけ**で結ぶ（因果なし）＝交絡の試験。
    """
    rng = np.random.default_rng(seed)
    idx = pd.DatetimeIndex(np.concatenate(
        [pd.date_range(f"{2010+i}-07-01", f"{2010+i}-08-31").values for i in range(nyear)]))
    s = (idx.dayofyear.to_numpy() - 182) / 61.0            # 夏の中の位置 0→1
    yeff = rng.normal(0, 1, nyear)
    ye = yeff[idx.year.to_numpy() - 2010]
    ar = np.zeros(len(idx))
    for i in range(1, len(idx)):
        ar[i] = 0.8 * ar[i - 1] + rng.normal(0, 1)          # 自己相関する擾乱
    th = np.clip(0.20 - 0.10 * s + 0.04 * ye + 0.012 * ar, 0.02, None)
    Rg = 220 - 40 * s + 40 * rng.standard_normal(len(idx))
    if kind == "threshold":
        f = np.clip(th / 0.15, 0, 1)                        # θ*=0.15
        gLE = Rg * (0.15 + 0.5 * f) + rng.normal(0, 12, len(idx))
        gH = Rg * (0.55 - 0.4 * f) + rng.normal(0, 12, len(idx))
    elif kind == "season":
        # **θ は使わない**。季節トレンドだけでフラックスが決まる＝θ との相関は見せかけ。
        gLE = Rg * (0.20 + 0.40 * (1 - s)) + rng.normal(0, 12, len(idx))
        gH = Rg * (0.50 - 0.30 * (1 - s)) + rng.normal(0, 12, len(idx))
    else:
        gLE = Rg * (0.35 + 0.5 * th) + rng.normal(0, 12, len(idx))
        gH = Rg * (0.35 - 0.2 * th) + rng.normal(0, 12, len(idx))
    return pd.DataFrame({"th": th, "Rg": Rg, "gLE": np.clip(gLE, 0, None),
                         "gH": np.clip(gH, 0, None)}, index=idx)


def report(d):
    print("    【Rg のみ制御】")
    show(strata_stats(d))
    rows, k = placebo_strata(d)
    print(f"    【↑のプラセボ＝年の入れ替え・最悪値（採用 {k} 年ずらし）】")
    show(rows)
    print("    【Rg＋日付を制御＝共通の季節トレンドを除く】")
    show(strata_stats(d, with_doy=True))
    rows2, k2 = placebo_strata(d, with_doy=True)
    print(f"    【↑のプラセボ（採用 {k2} 年ずらし）】")
    show(rows2)


def main():
    p = argparse.ArgumentParser(description="蒸発レジームを分けるのは θ の水準か（探索的）")
    p.add_argument("--sites", nargs="+", default=[])
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--qc-max", type=int, default=None)
    p.add_argument("--synth", action="store_true")
    a = p.parse_args()

    print("=== 旗71 v3：蒸発レジームを分けるのは θ の水準か（**探索的**）===")
    print("  v3 の修正：**年ブロック・ブート**（v1/v2 は日ブート＝CI が過小＝欠陥11件目）と、")
    print("  **年の入れ替えプラセボ**（循環シフトは実質わずかなずれになりうる＝欠陥10件目の再修正）。")
    print("  θ の深度はサイト間で不統一（旗33）＝**絶対値を横に比べない**。検定は各サイト内部。\n")

    if a.synth or not a.sites:
        print("  ── 合成検証 ──")
        for kind, lab in [
                ("threshold", "**閾値あり**（θ*=0.15 未満でだけ水律速）"),
                ("none", "閾値なし（全域で同じ依存）"),
                ("season", "**季節交絡だけ**（θ は無関係・共通の夏内トレンドのみ）")]:
            print(f"  ━ {lab} ━")
            report(make_synth(kind))
            print()
        print("  → 期待：閾値ありは**乾の層でだけ反転**しプラセボは反転なし。")
        print("     閾値なしはどの層も似た弱さ。**季節交絡は日付制御で消え、プラセボにも同じ模様が出る**。")
        if not a.sites:
            return
        print()

    print("  ── 実データ ──")
    for s in a.sites:
        try:
            d, nyr = daily_energy(s, a.month, a.qc_max)
        except Exception as e:
            print(f"  ━ {s} ━ 読み込み失敗 {type(e).__name__}: {str(e)[:120]}\n"); continue
        if "th" not in d or "gLE" not in d or len(d) < 60:
            print(f"  ━ {s} ━ 日数/変数が足りない（{len(d)}日）\n"); continue
        print(f"  ━ {s} ━（{len(d)}日・{nyr}年・月={a.month}）")
        try:
            from japanflux_pn.config import AnalysisConfig
            from japanflux_pn.sites import get_site
            from japanflux_pn.preprocess import load_raw_all
            cfg = AnalysisConfig(qc_max=a.qc_max) if a.qc_max is not None else AnalysisConfig()
            raw = load_raw_all(get_site(s), cfg)
            raw = raw[raw.index.month.isin(a.month)]
            if "P" in raw.columns:
                pyr = raw["P"].groupby(raw.index.year).sum()
                print(f"    対象月の降水量：年あたり中央 {np.nanmedian(pyr):.0f} mm"
                      f"（{np.nanmin(pyr):.0f}–{np.nanmax(pyr):.0f}）← **乾湿はこれで比べる**")
        except Exception as e:
            print(f"    降水量の取得に失敗 {type(e).__name__}: {str(e)[:80]}")
        report(d)
        print()

    print("  === 読み方（v3）===")
    print("  ★とみなせるのは、**実測が同じ層のプラセボを明確に上回る**場合だけである。")
    print("  **プラセボと同程度なら、それは θ の効果ではない**（v1/v2 はここで誤りかけた）。")
    print("  **日付制御で消える結果は、θ ではなく季節進行で説明できる**。")
    print("  留保：")
    print("   ・**仮説は旗70 の結果を見た後に立てた＝探索的**。確認には別データか事前登録が要る。")
    print("   ・**範囲制限**：層に切ると θ の幅が狭まり相関は機械的に縮む。各層の sd を併記した。")
    print("     sd が小さい層で r が小さくても、それだけでは『依存が無い』とは言えない。")
    print("   ・年ブロック・ブートは**年が少ないと CI が広くなる**（4–6年）。")
    print("     ＝v1/v2 より ✓ が減るのは**正しい方向の変化**であって、検出力が落ちたのではない。")
    print("   ・Seneviratne et al. (2010) の soil moisture regime と同型＝**新発見ではない**。")
    print("   ・γH/γLE は共通 w'（旗35）・θ 深度不統一（旗33）は解消していない。")


if __name__ == "__main__":
    main()
