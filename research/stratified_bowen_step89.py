"""旗89：**Bowen 反転の季節依存は、θ と Rg の水準に還元されるか**（事前登録 step89）。

旗88 で **A-3 の半分が季節に依る**と分かった——**蒸発側（θ→γLE|Rg>0）は北米 9/9 で成立**するが、
**Bowen 反転（＋θ→γH|Rg<0）は 4/9**（秋 3・冬 1）。
**「季節」は札であって物理量ではない。** 秋と春で何が違うのかを**量で**言えるかを確かめる。

**候補**：反転が見えるには「水がある」だけでなく**「使えるエネルギーもある」**必要があるのでは。

**層別の切り方は `PREREGISTRATION_step89.md` で固定済み**（実データを見る前に commit した）：
  ・層別変数は **θ と Rg の日平均生値**（アノマリではない＝**水準**を問うている）
  ・しきい値は**サイトごとの中央値**を**全月プール**で計算（季節を層の定義に埋め込まない）
  ・**2×2 の 4 セル**（三分位にしない＝検出力を優先）
  ・統計量は**旗36 と同一のものを呼ぶ**（新しい検出器を作らない）
  ・**判定は符号と CI のみ。セル間で r の大きさを比べない**（範囲制限で必ず減衰するため）

**安全弁（事前登録済み）**：Rg は季節と強く相関するので、**θ高×Rg高 が秋ばかり**になりうる。
**単一季節が 70% 以上を占めたセルは「季節と分離できていない」**とし、**還元の証拠に使わない**。

    python research/stratified_bowen_step89.py                 # 合成で検証（既定）
    python research/stratified_bowen_step89.py --real \\
        --sites US-Wkg US-Whs US-SRM MN-Hst MN-Nkh MN-Kbu
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
from evaporation_regime_step36 import daily_energy, _fmt

CELLS = ("θ高×Rg高", "θ高×Rg低", "θ低×Rg高", "θ低×Rg低")
MIN_DAYS, MIN_YEARS, SEASON_CAP = 60, 3, 0.70   # 事前登録で固定した下限・上限

# **季節の定義は「セルの中身を記述する」ためだけに使う**（層の定義には使わない）。
# 全 12 月を覆う必要があるので**気象学的四季**を採る。**検定そのものには影響しない**。
SEASON = {3: "春", 4: "春", 5: "春", 6: "夏", 7: "夏", 8: "夏",
          9: "秋", 10: "秋", 11: "秋", 12: "冬", 1: "冬", 2: "冬"}


def cell_of(d):
    """サイトごとの中央値で 2×2 に割る。**全月プールで計算する**。"""
    tmed, rmed = d["th"].median(), d["Rg"].median()
    hi_t, hi_r = d["th"] >= tmed, d["Rg"] >= rmed
    lab = pd.Series(index=d.index, dtype=object)
    lab[hi_t & hi_r] = "θ高×Rg高"
    lab[hi_t & ~hi_r] = "θ高×Rg低"
    lab[~hi_t & hi_r] = "θ低×Rg高"
    lab[~hi_t & ~hi_r] = "θ低×Rg低"
    return lab, float(tmed), float(rmed)


def test_cell(d):
    """1 セルの Bowen 反転を測る。**統計量は旗36 と同一**（`_boot_ci`・年ブロック）。"""
    if len(d) < MIN_DAYS or d.index.year.nunique() < MIN_YEARS:
        return None
    yr = d.index.year.to_numpy()
    le = _boot_ci(d["gLE"].to_numpy(), d["th"].to_numpy(), [d["Rg"].to_numpy()], blocks=yr)
    h = _boot_ci(d["gH"].to_numpy(), d["th"].to_numpy(), [d["Rg"].to_numpy()], blocks=yr)
    return {"n": len(d), "yrs": int(d.index.year.nunique()), "le": le, "h": h}


def reversed_(res):
    """Bowen 反転＝θ→γLE>0 と θ→γH<0 が**両方とも CI で 0 を跨がない**。"""
    if res is None:
        return None
    (rl, cl, _), (rh, ch, _) = res["le"], res["h"]
    ok_l = np.isfinite(rl) and cl and cl[0] > 0
    ok_h = np.isfinite(rh) and ch and ch[1] < 0
    return bool(ok_l and ok_h)


def season_mix(idx):
    s = pd.Series([SEASON[m] for m in idx.month])
    frac = s.value_counts(normalize=True)
    return frac, float(frac.iloc[0]), str(frac.index[0])


def run_site(site, d):
    print(f"\n  ━━ {site} ━━")
    lab, tmed, rmed = cell_of(d)
    print(f"    しきい値（全月プールの中央値）：θ={tmed:.3f}／Rg={rmed:.1f}"
          f"／全日数 {len(d):,}／年数 {d.index.year.nunique()}")
    print(f"    {'セル':<10}{'日数':>7}{'年':>4}  {'季節内訳':<26}"
          f"{'θ→γLE|Rg':>10}{'CI':>16}  {'θ→γH|Rg':>10}{'CI':>16}  判定")
    out = {}
    for c in CELLS:
        sub = d[lab == c]
        if sub.empty:
            print(f"    {c:<10}{0:>7}"); out[c] = None; continue
        frac, top, topname = season_mix(sub.index)
        mix = " ".join(f"{k}{v:.0%}" for k, v in frac.items())
        res = test_cell(sub)
        if res is None:
            print(f"    {c:<10}{len(sub):>7}{sub.index.year.nunique():>4}  {mix:<26}"
                  f"{'—':>10}{'':>16}  {'—':>10}{'':>16}  **判定しない**（下限未満）")
            out[c] = None; continue
        rev = reversed_(res)
        # **事前登録の安全弁**：単一季節が 70% 以上なら証拠に使わない
        blocked = top >= SEASON_CAP
        mark = ("**季節と分離できていない**" if blocked
                else ("**Bowen反転**" if rev else "反転せず"))
        print(f"    {c:<10}{res['n']:>7}{res['yrs']:>4}  {mix:<26}"
              f"{_fmt(res['le'])}  {_fmt(res['h'])}  {mark}")
        if blocked:
            print(f"       ↑ {topname} が {top:.0%}＝**層別が季節を分離していない**"
                  f"（事前登録の 70% 規則）")
        out[c] = None if blocked else rev
    # ── **本当の判別**：θ高×Rg高 の**中を季節で割る**（修正版・下の注記を参照）──
    # **合成検証で分かった欠陥**：θ の季節周期が「反転する季節」とそろっていると、
    # **θ で切ることが季節で切ることと同じになり**、`season`（秋だけ反転）と
    # `th_only`（θ だけで決まる）が**同じ絵になる**。**セルの季節内訳を見ても分からない**
    # ——反転を駆動する季節がセル内で少数派でも、セル全体としては反転して見えるため。
    # ＝**還元されたかどうかは、セルの中を季節で割って初めて言える。**
    hh = d[lab == "θ高×Rg高"]
    print(f"    ── θ高×Rg高 の**中を季節で割る**（ここが本当の判別）──")
    per = {}
    for s in ("春", "夏", "秋", "冬"):
        sub = hh[[SEASON[m] == s for m in hh.index.month]]
        res = test_cell(sub)
        if res is None:
            print(f"      {s}：日数 {len(sub):>5}＝**判定しない**（下限未満）")
            per[s] = None; continue
        rev = reversed_(res)
        per[s] = rev
        print(f"      {s}：日数 {res['n']:>5}／年 {res['yrs']:>3}  "
              f"{_fmt(res['le'])}  {_fmt(res['h'])}  "
              f"{'**Bowen反転**' if rev else '反転せず'}")
    ok = [s for s, v in per.items() if v is True]
    ng = [s for s, v in per.items() if v is False]
    if len(ok) + len(ng) < 2:
        print(f"      → **判定できる季節が 2 未満**＝**このサイトでは還元を判定できない**")
        out["_within"] = None
    elif not ng:
        print(f"      → **判定できた季節すべて（{'・'.join(ok)}）で反転**"
              f"＝**このサイトでは季節が効いていない**")
        out["_within"] = True
    else:
        print(f"      → **季節で分かれる**（反転：{'・'.join(ok) or 'なし'}／"
              f"反転せず：{'・'.join(ng)}）＝**θ と Rg では説明しきれていない**")
        out["_within"] = False
    return out


def synth(kind, years=8, seed=0):
    """**道具が「θ と Rg がともに高いときだけ反転する」系列を検出できるか**を確かめる。

    ``kind``：
      ・``both``   —— **θ高 かつ Rg高のときだけ**顕熱が潜熱へ振り替わる（**H1/H2 が真**）
      ・``th_only``—— **θ高なら Rg に依らず**振り替わる（**Rg は要らない**）
      ・``season``—— **秋の日だけ**振り替わる（**季節そのものが要因＝還元されない**）
    **季節構造（Rg と θ の年周期）を必ず入れる**——入れなければ層別の試験にならない。
    """
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2010-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy()
    Rg = 150 + 120 * np.sin(2 * np.pi * (doy - 80) / 365) + rng.normal(0, 25, len(idx))
    Rg = np.clip(Rg, 5, None)
    # θ は秋に高く春に低い（モンスーン後に湿る＝北米南西部を模す）
    th = np.clip(0.18 + 0.10 * np.sin(2 * np.pi * (doy - 200) / 365)
                 + rng.normal(0, 0.03, len(idx)), 0.02, 0.6)
    hi_t, hi_r = th >= np.median(th), Rg >= np.median(Rg)
    if kind == "both":
        on = hi_t & hi_r
    elif kind == "th_only":
        on = hi_t
    else:                                   # season
        on = pd.Series(idx.month).isin([9, 10]).to_numpy()
    # 反転が「効いている」日だけ、θ が顕熱を潜熱へ振り替える
    beta = np.where(on, 1.6, 0.0)
    gLE = Rg * (0.25 + beta * (th - th.mean())) + rng.normal(0, 8, len(idx))
    gH = Rg * (0.45 - beta * (th - th.mean())) + rng.normal(0, 8, len(idx))
    return pd.DataFrame({"th": th, "Rg": Rg,
                         "gLE": np.clip(gLE, 0, None), "gH": np.clip(gH, 0, None)},
                        index=idx)


def main():
    p = argparse.ArgumentParser(description="旗89：季節を θ と Rg に置き換える")
    p.add_argument("--real", action="store_true", help="実データを使う")
    p.add_argument("--sites", nargs="+",
                   default=["US-Wkg", "US-Whs", "US-SRM", "MN-Hst", "MN-Nkh", "MN-Kbu"])
    p.add_argument("--qc-max", type=int, default=None)
    a = p.parse_args()

    print("=== 旗89：Bowen 反転の季節依存は θ と Rg の水準に還元されるか ===")
    print("  **事前登録 step89 で層別の切り方・判定規則・安全弁を固定済み**。")
    print("  判定は**符号と CI のみ**。**セル間で r の大きさは比べない**"
          "（θ で切ると必ず減衰する＝保守的な検定）。")
    print(f"  下限：日数≥{MIN_DAYS}・年数≥{MIN_YEARS}／"
          f"単一季節 ≥{SEASON_CAP:.0%} のセルは**証拠に使わない**。")

    tally = {c: {"yes": 0, "no": 0, "skip": 0} for c in CELLS}
    within = {"yes": 0, "no": 0, "skip": 0}
    if not a.real:
        print("\n  【合成データで道具を検証する】**実データに触れる前に**、")
        print("  『θ高かつRg高のときだけ反転する』『θだけで決まる』『秋だけ反転する』の")
        print("  三つを作り、**道具がそれぞれを区別できるか**を見る。")
        for kind, want in (("both", "θ高×Rg高 だけで反転すべき"),
                           ("th_only", "θ高の 2 セルで反転すべき"),
                           ("season", "季節が要因＝70% 規則に触れるか反転が偏るべき")):
            print(f"\n  ===== 合成 `{kind}` —— 期待：{want} =====")
            run_site(f"合成-{kind}", synth(kind))
        print("\n  → **三つを区別できていれば、道具は使える。**")
        print("     区別できなければ**実データに進まない**（旗52 の作法）。")
        return

    for s in a.sites:
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max)
        except Exception as e:
            print(f"\n  ━━ {s} ━━\n    読み込み失敗 {type(e).__name__}: {str(e)[:120]}")
            continue
        if len(d) < MIN_DAYS or not {"th", "Rg", "gLE", "gH"} <= set(d.columns):
            print(f"\n  ━━ {s} ━━\n    変数不足（{len(d)} 日・列 {list(d.columns)}）")
            continue
        res = run_site(s, d)
        for c in CELLS:
            key = "skip" if res[c] is None else ("yes" if res[c] else "no")
            tally[c][key] += 1
        w = res.get("_within")
        within["skip" if w is None else ("yes" if w else "no")] += 1

    print(f"\n  === 集計（事前登録の判定規則に当てる）===")
    print(f"    {'セル':<10}{'反転あり':>8}{'反転なし':>8}{'使えない':>9}")
    for c in CELLS:
        t = tally[c]
        print(f"    {c:<10}{t['yes']:>8}{t['no']:>8}{t['skip']:>9}")
    print(f"\n    **θ高×Rg高 の中を季節で割った結果**（**修正版の主判別**）：")
    print(f"      季節に依らず反転 {within['yes']} サイト／"
          f"季節で分かれる {within['no']} サイト／判定できない {within['skip']} サイト")
    n_judged = {c: tally[c]["yes"] + tally[c]["no"] for c in CELLS}
    n_w = within["yes"] + within["no"]
    hh, hl, lh = "θ高×Rg高", "θ高×Rg低", "θ低×Rg高"
    print(f"\n  === 結論 ===")
    if n_w >= 3 and within["no"] > within["yes"]:
        print("  **▲還元されない**——**θ高×Rg高 の中でも季節によって反転したりしなかったりする**。")
        print("  ＝**季節そのものに別の要因**（フェノロジー等）が在る。**季節依存は未解明のまま残す。**")
        print("  （**この判別は合成検証で設計を直して初めて可能になった**＝修正前の 2×2 だけでは")
        print("    『θ だけで決まる』と見分けがつかなかった。）")
    elif n_w >= 3 and within["yes"] > within["no"]:
        print("  **★季節は θ と Rg の水準に還元された**——**θ高×Rg高 の中では、"
              "どの季節でも反転する**。")
        print("  ＝A-3 を『**θ と Rg がともに高いとき Bowen 反転が起きる**』と書き換える。")
    elif n_judged[hh] < 4:
        print(f"  **判定しない**——θ高×Rg高 で判定できたサイトが {n_judged[hh]} で 4 未満、")
        print(f"  かつ季節割りで判定できたサイトも {n_w} で 3 未満。")
    else:
        maj = lambda c: tally[c]["yes"] > n_judged[c] / 2 if n_judged[c] else False
        if maj(hh) and not maj(hl) and not maj(lh):
            print("  **★季節は θ と Rg の水準に還元された**"
                  "＝A-3 を『θ と Rg がともに高いとき Bowen 反転が起きる』と書き換える。")
        elif maj(hh) and maj(hl):
            print("  **○θ だけで決まる**＝**Rg は要らない**。候補『エネルギーも要る』を取り下げる。")
        elif maj(hh) and maj(lh):
            print("  **○Rg だけで決まる**＝**『水マスター』という読み方を見直す**"
                  "（A-3 に不利な結果である）。")
        else:
            print("  **▲還元されない**＝**季節そのものに別の要因**（フェノロジー等）が在る。"
                  "季節依存は未解明のまま残す。")
    print("\n  留保（事前登録に書いた通り）：")
    print("   ・**独立クラスタは 3 つ**（モンゴル1・Walnut Gulch・Santa Rita）＝6サイト≠6反復。")
    print("   ・**θ で切ると相関は必ず減衰する**＝**出なかったことは「無い」を意味しない**。")
    print("   ・**『θ と Rg がともに高い』はサイトごとの中央値**＝『同じ物理条件』ではなく")
    print("     **『そのサイトにとって湿って明るい日』**である。")
    print("   ・**反転が出ても機構は分からない**（旗71）。**いつ出るかを精密にするだけ**。")


if __name__ == "__main__":
    main()
