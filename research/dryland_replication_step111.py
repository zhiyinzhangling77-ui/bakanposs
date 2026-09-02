"""旗111：**同じ母集団で、「雨または秋」の型は再現するか**（CN-Du2・事前登録 step111）。

## 旗109 との違い

**旗109 は「別の母集団への外挿」だった**（湿潤な森林・年 688–1356 mm）。
**本検定は違う**——**CN-Du2 は年 468 mm で、既知の乾燥地（293–354 mm）と同じ側にある**（旗110）。
**＝これは「再現」の検定である。**

## 二段の問い

  1. **CN-Du2 で Bowen 反転が起きるか**（＝**A-3 が 4 つ目の独立クラスタで再現するか**）
     **中国のサイトで反転を測るのは、この研究で初めてである。**
  2. **起きるなら、その型は何か**（`rain_only`／`season_only`／`both`／**`or`**）

## **道具は旗109 のものをそのまま使う**

**`run_site` を import して呼ぶ**——**統計量・しきい値・下限・`DELTA_FLOOR` を一つも変えない。**
**US-SRM も同じ道具で一緒に走らせる**（**旗106/107 は別々の道具で測った**）。
**一本の経路で両者を出すので、「型が同じか」を比べるときに実装差が入らない。**
**US-SRM の結果が旗106/107 と食い違ったら、それ自体を欠陥として記録する。**

    python research/dryland_replication_step111.py            # 合成で検証（既定）
    python research/dryland_replication_step111.py --real     # 実データ（/mnt/hdd）
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from humid_forest_pattern_step109 import run_site, synth, DELTA_FLOOR, PATTERN
from rain_history_probe_step103 import daily_precip
from evaporation_regime_step36 import daily_energy
from runlog import tee_stdout

NEW = "CN-Du2"                 # 旗110 で唯一、目安（531 mm 以下）を満たした新サイト
REF = "US-SRM"                 # **同じ道具で走らせ直す**（旗106/107 との突き合わせ用）
# 旗106/107 が US-SRM について出した値（**食い違いを検出するために書き留める**）
SRM_KNOWN = {"秋全体 反転": True, "秋直後 反転": True, "秋遠い 反転": True,
             "遠い Δ_H": -0.38, "直後 Δ_H": -0.25, "旗107 の型": "or（雨または秋）"}


def main():
    ap = argparse.ArgumentParser(description="旗111：同じ母集団で型は再現するか")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()
    tee_stdout("step111")

    print("=== 旗111：同じ母集団で、「雨または秋」の型は再現するか（CN-Du2）===")
    print("  **旗109 と違い、これは『再現』の検定である**——"
          "**CN-Du2 は年 468 mm で既知の乾燥地と同じ側**（旗110）。")
    print("  **道具は旗109 のものをそのまま使う**（統計量・しきい値・下限を一つも変えない）。")
    print(f"  **`DELTA_FLOOR` = {DELTA_FLOOR}**／**US-SRM も同じ道具で走らせて突き合わせる。**")
    print("  **二段**：①CN-Du2 で反転するか（A-3 の再現）②型は何か。")

    if not a.real:
        print("\n  【合成データで検証する】**旗109 と同じ 5 通りを走らせ、"
              "5 枝すべてが到達することを確かめる**。")
        print("  **道具を変えていないので結果は同じはず**——**同じであることの確認が目的。**")
        print("  **違っていたら、私がどこかを壊した証拠なので実データに進まない。**")
        want = {"rain_only": "rain_only（雨だけ）", "season_only": "season_only（季節だけ）",
                "both": "both（雨かつ秋）", "or": "or（雨または秋）",
                "no_reversal": None}
        got = {}
        for k, w in want.items():
            print(f"\n  ===== 合成 `{k}` =====")
            d, P = synth(k)
            got[k] = run_site(f"合成/{k}", d, P)
        print("\n  === 合成のまとめ ===")
        ok = True
        for k, w in want.items():
            hit = (got[k] == w)
            ok &= hit
            print(f"    {k:<12}期待 {str(w):<24}実際 {got[k]}  {'✔' if hit else '**✘**'}")
        print(f"\n  → **5 枝すべて一致：{ok}**"
              f"{'' if ok else '  ← **一致しない枝がある＝実データに進まない**'}")
        return

    res = {}
    for s in (NEW, REF):
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max)
            P = daily_precip(s, a.qc_max)
        except Exception as e:
            print(f"\n  ━━ {s} ━━\n    読めない {type(e).__name__}: {str(e)[:90]}")
            res[s] = None; continue
        if P is None or P.dropna().empty:
            print(f"\n  ━━ {s} ━━\n    **降水 P が無い**"); res[s] = None; continue
        res[s] = run_site(s, d, P)

    print("\n  === 第一段：A-3 は 4 つ目の独立クラスタで再現したか ===")
    if res.get(NEW) is None:
        print(f"  **{NEW} は門①-a を通らなかった**——**秋全体で Bowen 反転しない。**")
        print("  → **▲A-3 は 4 つ目の乾燥地で再現しなかった。**")
        print("  **「乾燥地で水がマスター」の一般性は下がる。** **そう書く。**")
        print("  **これは失敗ではなく結果である**——"
              "**A-3 は旗70 で一度 n=1 まで縮み、旗82 で 3 に回復した。また縮みうる。**")
    else:
        print(f"  **★{NEW} の秋全体で Bowen 反転した**"
              f"——**A-3 が 4 つ目の独立クラスタで再現**（**中国では初めて**）。")

    print("\n  === 第二段：型 ===")
    for s in (NEW, REF):
        print(f"    {s:<10}{res.get(s) or '判定しない（門①-a を通らない）'}")
    if res.get(NEW) and res.get(REF):
        if res[NEW] == res[REF]:
            print(f"\n  **○同じ型が 2 つの独立クラスタで出た：{res[NEW]}**")
            if res[NEW].startswith("or"):
                print("  **旗107 で事後に見つけた『雨または秋』が、事前登録した再現を得た。**")
        else:
            print(f"\n  **○型はクラスタごとに違う**（{NEW}: {res[NEW]}／{REF}: {res[REF]}）。")
            print("  **「雨または秋」は Santa Rita 固有のまま。**")
        print("  **どちらにせよ「2 クラスタで見た」までしか言えない。** **「一般に」とは書かない。**")

    print(f"\n  === 突き合わせ：{REF} は旗106/107 と同じ答えを返したか ===")
    print(f"    旗107 が出した型：**{SRM_KNOWN['旗107 の型']}**"
          f"／今回：**{res.get(REF)}**")
    if res.get(REF) and res[REF] != SRM_KNOWN["旗107 の型"]:
        print("    → **食い違った。** **これは実装差か、私がどこかを壊した証拠である。**")
        print("    **食い違いを解くまで、本検定の結論は保留する。**")
    elif res.get(REF):
        print("    → **一致した。** **一本の経路で両クラスタを比べられている。**")
    print(f"    （旗107 の実測：`遠い` Δ_H {SRM_KNOWN['遠い Δ_H']:+.2f}／"
          f"`直後` Δ_H {SRM_KNOWN['直後 Δ_H']:+.2f}）")

    print("\n  留保（事前登録どおり）：")
    print("   ・**新しいクラスタは 1 つだけ。** **型の一致は「2 クラスタで」までしか言えない。**")
    print("   ・**CN-Du2 は 10 年**で US-SRM（19–22 年）より短い。**CI は広くなる。**")
    print("   ・**CN-Du2 は冬の雨がほとんど無い**（0%／北米の乾燥地は 14–20%）。")
    print("     **「同じ乾燥地」と一括りにできない。**")
    print("   ・**θ の単位が違う**（ChinaFlux は分数・旗110 の欠陥 #41）。")
    print("     **`cell_of` も偏 Spearman も単位に依らないので検定は無事**だが、")
    print("     **θ の値を書くときは％に直す。**")
    print("   ・**ChinaFlux をこの研究で本格的に使うのは初めて**である。")
    print("     **警告が出たら、結果より先に警告を読む。**")


if __name__ == "__main__":
    main()
