"""旗71：**何が水制限とエネルギー制限を分けるのか**——旗70 が残した問い。

旗70 で、**同じモンゴル半乾燥ステップ・同じ草原・113km 離れた MN-Skt** で Bowen 反転が
再現しなかった。しかも corr(Rg,γLE)=+0.51（4サイト中最大）＝**本当にエネルギー制限寄り**らしい。

＝問いは「乾燥地で一般に成り立つか」ではなく、**何が両者を分けるのか**に変わった。

## 仮説（**事後に立てた＝探索的**であることを明記する）

**θ 閾値仮説**：蒸発の水分律速は**θ がある水準を下回るときだけ**現れる。
MN-Skt は単に**湿っている**ため、θ が律速域に入らず、エネルギー制限に見えている。
＝地点固有の謎ではなく、**θ の分布の違い**で説明できる。

これは新説ではない。**Seneviratne et al. (2010, Earth-Sci Rev) の
soil moisture–climate regime（wet / transitional / dry を θ_crit が分ける）そのもの**である。
＝本ツールが確かめるのは新発見ではなく、**自分のデータがその枠組みと整合するか**である
（Q10 が DAMM の再発見かもしれないのと同型の位置づけ）。

## 設計上の壁と、その回避

**θ の深度がサイト間で不統一**（旗33）＝**絶対値の θ をサイト間で比べてはいけない**。
そこで：

  ・**サイト水準の乾湿比較には降水量 P を使う**（単位が同じで深度問題が無い）。
  ・**閾値の検定は各サイトの内部で行い、θ は「そのサイトの中での百分位」で切る**
    ＝深度が違っても**サイト内の相対的な乾湿**は意味を持つ。

これにより「サイトが4つしかない」問題も回避できる——**検定はサイト内の日で行う**ため、
**独立地点 n の少なさに依存しない**。

## 何を測るか

各サイト・各 θ 三分位（乾／中／湿）について：
  ・偏 Spearman θ→γLE|Rg と θ→γH|Rg（旗31/36 と同一の統計量）
  ・**その層の θ の標準偏差と n**＝**範囲制限（range restriction）で相関が縮む効果を読者が見える**ように
  ・**プラセボ**：θ を位相シフトした系列で同じ三分位を切る（旗54 と同じ device）。
    層別という操作そのものが模様を作るなら、プラセボにも同じ模様が出るはずである。

## 予測（探索的だが、先に書く）

**乾いた層でだけ Bowen 反転が出る**——**MN-Skt を含めて**。
出なければ θ 閾値仮説は棄却され、「地点差は θ の分布では説明できない」と記す。

    # 合成で先に検出器を確かめる（閾値あり／なし）
    python research/regime_threshold_step71.py --synth
    # 実データ
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
from moisture_control_atlas_step31 import partial_spearman, _boot_ci
from evaporation_regime_step36 import daily_energy

LABELS = ["乾(下位1/3)", "中(中位1/3)", "湿(上位1/3)"]

# ── 旗71 第1版の欠陥（自分の道具の欠陥・10件目）と、その修正 ──
# 第1版のプラセボは `np.roll(th, len//2)` の**単一シフト**だった。しかし本解析は
# **夏（7–8月）だけを連結した系列**なので、年あたりの日数のちょうど整数倍だけずらすと
# **7月が7月に、8月が8月に写る**＝**夏の中の季節進行が壊れない**。
# 実際 MN-Hst は 248日/6年 ≒ 62日/年 で 124 = ちょうど2年分となり、
# **プラセボが実測とほぼ同じ反転（+0.38/−0.25 対 実測 +0.40/−0.42）を出した**。
# ＝プラセボが機能していなかった。
#
# これは同時に**実質的な交絡の指摘**でもある：θ と γLE が**共通の夏内トレンド**を持てば、
# Rg を制御しても偏相関は出る。
#
# 修正は二つ：
#   (1) プラセボを**複数シフト**にし、**年境界の整数倍付近を除く**。基準線は |r| の**最大値**。
#   (2) 偏相関の制御に**日付（夏の中の位置）を加える**＝共通の季節トレンドを直接除去する。
N_PLACEBO = 12          # 試すシフトの本数
DOY_GUARD = 4           # 年境界の整数倍から±この日数は使わない


def _controls(sub, with_doy):
    """偏相関の制御変数。`with_doy` なら**夏の中の位置**も制御して季節トレンドを除く。"""
    ctrl = [sub["Rg"].to_numpy()]
    if with_doy and isinstance(sub.index, pd.DatetimeIndex):
        ctrl.append(sub.index.dayofyear.to_numpy().astype(float))
    return ctrl


def strata_stats(d, thcol="th", with_doy=False):
    """θ の**サイト内百分位**で三分位に切り、各層で偏相関を測る。"""
    th = d[thcol].to_numpy()
    ok = np.isfinite(th)
    if ok.sum() < 30:
        return None
    q1, q2 = np.nanpercentile(th[ok], [100 / 3, 200 / 3])
    edges = [(-np.inf, q1), (q1, q2), (q2, np.inf)]
    out = []
    for lo, hi in edges:
        m = np.isfinite(th) & (th > lo) & (th <= hi)
        sub = d.loc[m]
        if len(sub) < 20 or "gLE" not in sub:
            out.append(None); continue
        ctrl = _controls(sub, with_doy)
        le = _boot_ci(sub["gLE"].to_numpy(), sub[thcol].to_numpy(), ctrl)
        h = (_boot_ci(sub["gH"].to_numpy(), sub[thcol].to_numpy(), ctrl)
             if "gH" in sub else (np.nan, None, 0))
        out.append({"n": int(len(sub)), "sd": float(np.nanstd(sub[thcol])),
                    "mean": float(np.nanmean(sub[thcol])), "le": le, "h": h})
    return out


def _shifts(n, nyear):
    """使うシフト量。**年あたり日数の整数倍付近を除く**（第1版の欠陥＝夏内位相が壊れない）。"""
    per = max(int(round(n / max(nyear, 1))), 1)
    cand = []
    for k in range(1, N_PLACEBO + 1):
        s = int(round(n * k / (N_PLACEBO + 1)))
        if s <= 0 or s >= n:
            continue
        if per > 1 and min(s % per, per - s % per) <= DOY_GUARD:
            continue                                   # 年境界の整数倍付近は捨てる
        cand.append(s)
    return cand or [max(n // 3, 1)]


def placebo_strata(d, nyear=1, thcol="th", with_doy=False):
    """**複数の位相シフト**で同じ層別を行い、|r| が**最大**になった回を基準線として返す。

    θ の分布も自己相関も保ったまま γLE/γH との対応を壊す。単一シフトだと
    **たまたま夏内位相が保たれる**ことがあるため（旗71 第1版の欠陥）、複数試して**最悪値**を採る。
    """
    best = None
    for sh in _shifts(len(d), nyear):
        d2 = d.copy()
        d2[thcol] = np.roll(d2[thcol].to_numpy(), sh)
        rows = strata_stats(d2, thcol, with_doy)
        if rows is None:
            continue
        worst = max([abs(r["le"][0]) for r in rows if r and np.isfinite(r["le"][0])],
                    default=np.nan)
        if best is None or (np.isfinite(worst) and worst > best[0]):
            best = (worst, sh, rows)
    return (best[2], best[1]) if best else (None, None)


def _f(ci):
    """(r, CI, n) を『r [lo,hi]✓』の形に。CI が 0 を跨がねば ✓。"""
    r, c, _ = ci
    if not np.isfinite(r):
        return f"{'—':>20}"
    if not isinstance(c, tuple):
        return f"{r:+.2f}{'':>14}"
    sig = "✓" if (c[0] > 0 or c[1] < 0) else "·"
    return f"{r:+.2f} [{c[0]:+.2f},{c[1]:+.2f}]{sig}"


def show(name, rows):
    print(f"    {'層':<12}{'n':>5}{'θ平均':>8}{'θのsd':>8}  {'θ→γLE|Rg':<22}{'θ→γH|Rg':<22}判定")
    for lab, r in zip(LABELS, rows or []):
        if r is None:
            print(f"    {lab:<12}{'—':>5}{'—':>8}{'—':>8}  日数不足"); continue
        le_pos = isinstance(r["le"][1], tuple) and r["le"][1][0] > 0
        h_neg = isinstance(r["h"][1], tuple) and r["h"][1][1] < 0
        v = ("**Bowen反転**" if (le_pos and h_neg) else
             "水制限(蒸発のみ)" if le_pos else "反転なし")
        print(f"    {lab:<12}{r['n']:>5}{r['mean']:>8.3f}{r['sd']:>8.3f}  "
              f"{_f(r['le']):<22}{_f(r['h']):<22}{v}")


def make_synth(kind, days=900, seed=0):
    """閾値**あり**（θ<θ*でだけ蒸発が水に律速）と**なし**（どこでも同じ弱い依存）。"""
    rng = np.random.default_rng(seed)
    Rg = 200 + 60 * rng.standard_normal(days)
    th = np.clip(0.05 + 0.25 * rng.random(days), 0, 1)
    if kind == "threshold":
        f = np.clip(th / 0.15, 0, 1)              # θ*=0.15 で頭打ち
        gLE = Rg * (0.15 + 0.5 * f) + rng.normal(0, 12, days)
        gH = Rg * (0.55 - 0.4 * f) + rng.normal(0, 12, days)
    else:                                          # 閾値なし（全域で同じ弱い線形）
        gLE = Rg * (0.35 + 0.15 * th) + rng.normal(0, 12, days)
        gH = Rg * (0.35 - 0.05 * th) + rng.normal(0, 12, days)
    return pd.DataFrame({"th": th, "Rg": Rg,
                         "gLE": np.clip(gLE, 0, None), "gH": np.clip(gH, 0, None)})


def main():
    p = argparse.ArgumentParser(description="蒸発レジームを分けるのは θ の水準か（探索的）")
    p.add_argument("--sites", nargs="+", default=[])
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--qc-max", type=int, default=None)
    p.add_argument("--synth", action="store_true")
    a = p.parse_args()

    print("=== 旗71：蒸発レジームを分けるのは θ の水準か（**探索的**・旗70 が残した問い）===")
    print("  仮説は**結果を見た後に立てた**＝確認的検定ではない。確認には別データが要る。")
    print("  θ の深度はサイト間で不統一（旗33）＝**絶対値を横に比べない**。")
    print("  検定は**各サイトの内部**で、θ を**そのサイトの百分位**で三分位に切って行う。\n")

    if a.synth or not a.sites:
        print("  ── 合成検証：検出器は閾値を見つけるか、無いのに作らないか ──")
        for kind, lab in [("threshold", "**閾値あり**（θ*=0.15 未満でだけ水律速）"),
                          ("none", "閾値なし（全域で同じ弱い依存）")]:
            print(f"  {lab}")
            show(kind, strata_stats(make_synth(kind)))
            rows, sh = placebo_strata(make_synth(kind), nyear=1)
            print(f"    プラセボ（**複数シフトの最悪値**・採用シフト={sh}）：")
            show(kind, rows)
            print()
        print("  → 期待：閾値ありは**乾の層でだけ反転**、閾値なしは**どの層も似た弱さ**。")
        print("     プラセボはどちらも反転なし＝層別そのものは模様を作らない。")
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
        # サイト水準の乾湿は**深度問題の無い P** で見る（θ の絶対値では比べない）
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
                      f"（範囲 {np.nanmin(pyr):.0f}–{np.nanmax(pyr):.0f}）"
                      f"  ← **サイト間の乾湿はこれで比べる**")
            else:
                print("    降水量 P が無い＝サイト間の乾湿を P で比べられない")
        except Exception as e:
            print(f"    降水量の取得に失敗 {type(e).__name__}: {str(e)[:80]}")
        show(s, strata_stats(d))
        rows, sh = placebo_strata(d, nyear=nyr)
        print(f"    プラセボ（**複数シフトの最悪値**・採用シフト={sh}日）：")
        show(s, rows)
        print(f"    **季節も制御**（Rg に加えて夏の中の位置＝日付も偏らせる）：")
        show(s, strata_stats(d, with_doy=True))
        rows2, sh2 = placebo_strata(d, nyear=nyr, with_doy=True)
        print(f"    ↑のプラセボ（複数シフトの最悪値・採用シフト={sh2}日）：")
        show(s, rows2)
        print()

    print("  === 読み方 ===")
    print("  **どのサイトでも乾の層でだけ反転**が出れば、θ 閾値仮説と整合＝")
    print("  地点差は『θ の分布がどこにあるか』で説明でき、**地点固有の謎ではなくなる**。")
    print("  **MN-Skt の乾の層でも反転が出るか**が核心（出なければ仮説は棄却）。")
    print("  留保：")
    print("   ・**仮説は旗70 の結果を見た後に立てた＝探索的**。確認には別データか事前登録が要る。")
    print("   ・**範囲制限**：層に切ると θ の幅が狭まり相関は機械的に縮む。各層の sd を併記したのは")
    print("     読者がその効果を見えるようにするため。**sd が小さい層で r が小さくても、それだけでは")
    print("     『依存が無い』とは言えない**。")
    print("   ・プラセボに同じ模様が出るなら、それは層別操作そのものが作った模様である。")
    print("   ・**旗71 第1版の欠陥（自分の道具の欠陥10件目）**：プラセボが単一シフトだったため、")
    print("     夏だけを連結した系列で**年あたり日数の整数倍**ずれると**7月が7月に写り季節進行が")
    print("     壊れず**、MN-Hst でプラセボが実測とほぼ同じ反転を出した。→ **複数シフトの最悪値**に")
    print("     変更し、さらに**日付を制御した版**を併記して**共通の夏内トレンド**を直接除く。")
    print("   ・**日付を制御した版で消える結果は、θ ではなく季節進行で説明できる**＝そう読むこと。")
    print("   ・Seneviratne et al. (2010) の soil moisture regime の枠組みと同型＝**新発見ではない**。")
    print("     本ツールが確かめるのは『自分のデータがその枠組みと整合するか』である。")
    print("   ・γH/γLE は共通 w'（旗35）・θ 深度不統一（旗33）は解消していない。")


if __name__ == "__main__":
    main()
