"""旗90：**春と秋の差は、土壌温度 Ts で説明できるか**（事前登録 step90）。

旗89 で**季節依存は θ と Rg に還元されなかった**——**θ高×Rg高 に限っても春は 0/6**。
残る候補で**測定量で書けるもの**が `Ts`（地温は放射に遅れるので、**同じ Rg でも
春の土壌は冷たく秋は暖かい**）。**`Ts` は分割派生量ではない**（旗32/35 の問題が無い）。

**この検定は共通支持に懸かっている**：**春と秋の Ts 分布が重ならなければ、
どんな層別でも答えは出ない**。**だから重なりを先に測り、足りなければ「判定できない」で終える**
——**それは失敗ではなく「地温で分離している」という結果**である。

**事前登録 step90 で固定済み**：
  ・θ・Rg のしきい値は**旗89 と同一**（サイトごとの中央値・全月プール）＝**作り直さない**
  ・帯 ＝ **[春秋の 10 パーセンタイルの大きい方, 90 パーセンタイルの小さい方]**
  ・下限 ＝ 帯の中で**春・秋それぞれ 60 日以上かつ 3 年以上**
  ・判定は**符号と CI のみ**。**帯の内外で r の大きさを比べない**

**旗89 の欠陥26（層別変数が季節と周期を共有すると分離できない）を構造的に回避している**——
**セルの中を季節で割るのではなく、Ts を揃えたうえで季節ごとに別々に検定する**。

    python research/soiltemp_match_step90.py                   # 合成で検証（既定）
    python research/soiltemp_match_step90.py --real \\
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
from evaporation_regime_step36 import daily_energy, _fmt
# **旗89 の関数をそのまま呼ぶ**——しきい値も検定も**作り直さない**
from stratified_bowen_step89 import cell_of, test_cell, reversed_, MIN_DAYS, MIN_YEARS

SPRING, AUTUMN = (3, 4, 5), (9, 10, 11)
PLO, PHI = 10, 90                      # 帯の定義（事前登録で固定）


def band(ts_sp, ts_au):
    """春・秋の Ts 分布の**重なり帯**。事前登録どおり 10/90 パーセンタイルで作る。"""
    if len(ts_sp) == 0 or len(ts_au) == 0:
        return None
    lo = max(np.percentile(ts_sp, PLO), np.percentile(ts_au, PLO))
    hi = min(np.percentile(ts_sp, PHI), np.percentile(ts_au, PHI))
    return (float(lo), float(hi)) if hi > lo else (float(lo), float(hi))


def run_site(site, d):
    print(f"\n  ━━ {site} ━━")
    if "Ts" not in d.columns or d["Ts"].notna().sum() < MIN_DAYS:
        print("    **Ts が無い／足りない**＝このサイトでは判定できない")
        return None
    lab, tmed, rmed = cell_of(d)
    hh = d[(lab == "θ高×Rg高") & d["Ts"].notna()]
    sp = hh[[m in SPRING for m in hh.index.month]]
    au = hh[[m in AUTUMN for m in hh.index.month]]
    print(f"    θ・Rg のしきい値（旗89 と同一）：θ={tmed:.3f}／Rg={rmed:.1f}"
          f"／θ高×Rg高 は {len(hh):,} 日")
    if len(sp) == 0 or len(au) == 0:
        print(f"    春 {len(sp)} 日／秋 {len(au)} 日＝**片方が空**＝判定できない")
        return None
    print(f"    Ts 分布（θ高×Rg高 の中）："
          f"春 n={len(sp)} 中央{sp['Ts'].median():.1f} "
          f"[{np.percentile(sp['Ts'],PLO):.1f}–{np.percentile(sp['Ts'],PHI):.1f}]"
          f"／秋 n={len(au)} 中央{au['Ts'].median():.1f} "
          f"[{np.percentile(au['Ts'],PLO):.1f}–{np.percentile(au['Ts'],PHI):.1f}]")
    lo, hi = band(sp["Ts"].to_numpy(), au["Ts"].to_numpy())
    if hi <= lo:
        print(f"    → **帯幅が 0 以下**（[{lo:.1f}, {hi:.1f}]）"
              f"＝**春と秋は地温で完全に分離している**＝**この設計では答えられない**")
        return None
    spb = sp[(sp["Ts"] >= lo) & (sp["Ts"] <= hi)]
    aub = au[(au["Ts"] >= lo) & (au["Ts"] <= hi)]
    print(f"    **重なり帯 Ts ∈ [{lo:.1f}, {hi:.1f}]**（幅 {hi-lo:.1f}）"
          f"／残存：春 {len(spb)}日({len(spb)/len(sp):.0%})・"
          f"秋 {len(aub)}日({len(aub)/len(au):.0%})")
    out = {}
    for name, sub in (("春", spb), ("秋", aub)):
        res = test_cell(sub)
        if res is None:
            print(f"      {name}（帯の中）：日数 {len(sub)}／年 "
                  f"{sub.index.year.nunique() if len(sub) else 0}"
                  f"＝**判定しない**（下限 {MIN_DAYS}日・{MIN_YEARS}年 未満）")
            out[name] = None; continue
        rev = reversed_(res)
        out[name] = rev
        print(f"      {name}（帯の中）：日数 {res['n']:>5}／年 {res['yrs']:>3}  "
              f"{_fmt(res['le'])}  {_fmt(res['h'])}  "
              f"{'**Bowen反転**' if rev else '反転せず'}")
    if out.get("春") is None or out.get("秋") is None:
        print("      → **このサイトは判定しない**（帯の中の日数・年数が下限未満）")
        return None
    if out["春"] and out["秋"]:
        print("      → **帯の中では春も秋も反転**＝**Ts を揃えると差が消える**")
        return "both"
    if out["秋"] and not out["春"]:
        print("      → **帯の中でも秋は反転・春は反転せず**＝**Ts では説明されない**")
        return "autumn_only"
    if not out["秋"]:
        print("      → **帯の中では秋も反転しない**＝**帯に絞って検出力が落ちた公算**"
              "＝**春について何も言えない**")
        return "neither"
    return "spring_only"


def synth(kind, years=10, seed=0):
    """**Ts で説明される系列**と**されない系列**を作り、道具が区別できるか確かめる。

      ・``ts`` / ``season`` —— **現実的な地温**（年振幅 12℃・日々のばらつき 3.5℃）。
        **帯が 1℃ ほどしか作れず「判定できない」で終わる**——
        **これは実データへの予測でもある**：地温の季節振幅が日々のばらつきより
        大きければ、**春と秋には共通の土俵がほとんど無い**。
      ・``nooverlap`` —— ばらつきを 0.3℃ にした極端例。**帯幅が 0 以下**になるべき。
      ・``ts_overlap`` / ``season_overlap`` —— **作為的にばらつきを 9℃ に広げた検証用**。
        **肯定側の分岐（春も反転／秋だけ反転）を通すためだけ**に用意した。
        **実データがこうであることを主張するものではない。**
    """
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2010-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy()
    Rg = np.clip(150 + 120 * np.sin(2 * np.pi * (doy - 80) / 365)
                 + rng.normal(0, 25, len(idx)), 5, None)
    # **θ は二峰性にする**——第1版は秋に単峰の正弦にしたため、
    # **春が θ高×Rg高 に 1 日も入らず、検定したい状況を再現できていなかった**。
    # 実データでは春も 242–596 日入っている（旗89）。北米南西部は
    # **冬春の前線性降雨と夏のモンスーン**の二峰性なので、それを模す。
    wet = (np.exp(-0.5 * ((doy - 90) / 35.0) ** 2)       # 冬春の雨
           + np.exp(-0.5 * ((doy - 230) / 35.0) ** 2))   # モンスーン
    th = np.clip(0.14 + 0.10 * wet + rng.normal(0, 0.035, len(idx)), 0.02, 0.6)
    # **地温は放射に遅れる**（位相を 40 日ずらす）＝春は冷たく秋は暖かい。
    # **重なりを作るのは日々のばらつき**（曇天・降雨の連なり）である。
    # `*_overlap` は**作為的にばらつきを大きくした検証用**——
    # **道具の肯定側の分岐（春も反転する／秋だけ反転する）を通すためだけに用意した**。
    # **実データがこうであることを主張するものではない。**
    noise = {"nooverlap": 0.3, "ts_overlap": 9.0, "season_overlap": 9.0}.get(kind, 3.5)
    Ts = 15 + 12 * np.sin(2 * np.pi * (doy - 80 - 40) / 365) + rng.normal(0, noise, len(idx))
    if kind in ("season", "season_overlap"):
        on = pd.Series(idx.month).isin(AUTUMN).to_numpy()
    else:                                   # ts / nooverlap
        on = Ts >= np.median(Ts)
    beta = np.where(on, 1.6, 0.0)
    gLE = Rg * (0.25 + beta * (th - th.mean())) + rng.normal(0, 8, len(idx))
    gH = Rg * (0.45 - beta * (th - th.mean())) + rng.normal(0, 8, len(idx))
    return pd.DataFrame({"th": th, "Rg": Rg, "Ts": Ts,
                         "gLE": np.clip(gLE, 0, None), "gH": np.clip(gH, 0, None)},
                        index=idx)


def main():
    p = argparse.ArgumentParser(description="旗90：春と秋を Ts で揃える")
    p.add_argument("--real", action="store_true")
    p.add_argument("--sites", nargs="+",
                   default=["US-Wkg", "US-Whs", "US-SRM", "MN-Hst", "MN-Nkh", "MN-Kbu"])
    p.add_argument("--qc-max", type=int, default=None)
    a = p.parse_args()

    print("=== 旗90：春と秋の差は土壌温度 Ts で説明できるか ===")
    print("  **事前登録 step90 で帯の定義・下限・判定規則を固定済み**。")
    print("  **共通支持が無ければ「判定できない」で終える**——それは失敗ではなく、")
    print("  **『春と秋は地温で分離している』という結果**である。")

    if not a.real:
        print("\n  【合成データで道具を検証する】**実データに触れる前に**、")
        print("  **五つの場合**を作り、**道具がそれぞれ正しい枝に入るか**を見る。")
        for kind, want in (
                ("ts", "地温の年振幅が日々のばらつきより大きい＝**帯が細く判定できない**はず"),
                ("season", "同上（**現実的な地温では重なりが足りない**という予測でもある）"),
                ("nooverlap", "帯が作れず判定できないべき"),
                ("ts_overlap", "**重なりを作為的に広げた検証用**：春も秋も反転すべき"),
                ("season_overlap", "**同上**：秋だけ反転すべき")):
            print(f"\n  ===== 合成 `{kind}` —— 期待：{want} =====")
            print(f"  【判定】{run_site(f'合成-{kind}', synth(kind))}")
        print("\n  → **五つとも正しい枝に入れば道具は使える**。"
              "できなければ**実データに進まない**（旗52/89 の作法）。")
        print("  **`ts`/`season` が「判定できない」で終わるのは正しい挙動**であり、")
        print("  **同時に実データへの予測でもある**——**地温では春と秋の土俵が重ならないかもしれない**。")
        return

    tally = {}
    for s in a.sites:
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max, extra=("Ts",))
        except Exception as e:
            print(f"\n  ━━ {s} ━━\n    読み込み失敗 {type(e).__name__}: {str(e)[:120]}")
            continue
        tally[s] = run_site(s, d)

    print(f"\n  === 集計（事前登録の判定規則に当てる）===")
    for s, v in tally.items():
        label = {"both": "帯の中で春も秋も反転（Ts で説明される）",
                 "autumn_only": "帯の中でも秋だけ反転（Ts では説明されない）",
                 "neither": "帯の中では秋も反転せず（検出力低下＝何も言えない）",
                 "spring_only": "帯の中で春だけ反転（想定外）",
                 None: "判定できない"}.get(v, str(v))
        print(f"    {s:<9}{label}")
    n_both = sum(v == "both" for v in tally.values())
    n_au = sum(v == "autumn_only" for v in tally.values())
    n_nei = sum(v == "neither" for v in tally.values())
    n_judged = n_both + n_au + n_nei
    print(f"\n  === 結論 ===")
    if n_judged < 3:
        print(f"  **判定しない**——判定できたサイトが {n_judged} で 3 未満。")
        print(f"  **重なりが無かったのか日数が足りなかったのか**は、上の各サイトの行に書いてある。")
    elif n_both > n_au and n_both > n_nei:
        print("  **★Ts で説明された**——**地温を揃えると春と秋の差が消える**。")
        print("  ＝A-3 を『**θ が高く、かつ土壌が暖かいとき Bowen 反転が起きる**』と書き換える。")
        print("  **ただし「Ts が原因」とは言わない**——**Ts を揃えると差が消える**までである。")
    elif n_au >= n_both and n_au > n_nei:
        print("  **▲Ts でも説明されない**——**地温を揃えても春だけ反転しない**。")
        print("  ＝**水・エネルギー・熱のどれでもない要因**（フェノロジー等）が残る。")
        print("  **測定量での説明は打ち切る**と確定して記録する。")
    else:
        print("  **○判定保留**——**帯の中では秋も反転しなくなった**サイトが多い。")
        print("  ＝**帯に絞って検出力が落ちただけ**の公算が高く、**春について何も言えない**。")
    print("\n  留保（事前登録どおり）：")
    print("   ・**独立クラスタは 3 つ**＝6 サイト≠6 反復。")
    print("   ・**帯に絞れば n は必ず減り相関は減衰する**＝**符号と CI のみで判定**した。")
    print("   ・**Ts の深度は不明**（旗33/80：`TS_F_MDS_1` は最浅の慣例に従っただけ）。")
    print("   ・**Ts と一緒に動く量（フェノロジー・根の活性・大気の乾燥度）は揃っていない**")
    print("     ＝**『Ts で説明された』と出ても『Ts が原因』ではない**。")


if __name__ == "__main__":
    main()
