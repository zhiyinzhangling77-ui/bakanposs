"""旗113：**気候の札をやめ、門に決めさせる**（残り 3 サイト・事前登録 step113）。

## なぜ札をやめるのか

**旗110 の 531 mm も、旗112 の指標 28.1 も、私が引いた線である。**
**旗112 の線は、答えの分かっている CN-Du2 を外した**（指標 36.4 ＝湿潤側と判定したが、
**旗111 で反転が出たサイト**）。**原因は較正集合の交絡**——**手元の乾燥地はすべて暑く、
湿潤側は涼しい**ので、**指標は気温の差をなぞって分離していた。**

**本当に効いているのは「そのサイトで反転が起きるか」であり、気候の札ではない。**
**それは門①-a が直接測っている。** **札で前もって選ばず、門に決めさせる。**

## 対象（**これで 4 群が揃うサイトの列挙は尽きる**）

**CN-Ha2（570 mm・−0.9 °C）・CN-Aro（618 mm・0.6 °C）・CN-Xsh（1356 mm・20.6 °C・熱帯）。**
**CN-HbC は除く**——**`description` が別サイト（`CN-Erg`）を指しており**（旗110）、
**フォルダの中身を確かめるまで使えない。** **気候の理由ではなく、来歴の理由である。**

## **釣りにならないための縛り**（事前登録 step113 で固定）

**型は 4 つあるので、当てずっぽうでも 1 サイトが `or` になる確率は約 1/4。**
**3 サイトなら少なくとも 1 つが `or` になる確率は約 58%。**
  ・**`or` が 1 サイトだけ出ても「再現した」とは書かない。**
  ・**`or` が 2 クラスタ以上で初めて「型が繰り返した」と書く。**
  ・**門を通ったサイトは、型が何であれ全部書く。拾わない。**

    python research/gate_decides_step113.py            # 合成で検証（既定）
    python research/gate_decides_step113.py --real     # 実データ（/mnt/hdd）
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from humid_forest_pattern_step109 import run_site, synth, DELTA_FLOOR
from rain_history_probe_step103 import daily_precip
from evaporation_regime_step36 import daily_energy
from runlog import tee_stdout

# 旗110 の線で落ちていた 3 件（**気候で選ばず、門に決めさせる**）
SITES = ("CN-Ha2", "CN-Aro", "CN-Xsh")
CLIMATE = {"CN-Ha2": "570 mm・−0.9 °C（寒冷）",
           "CN-Aro": "618 mm・0.6 °C（寒冷）",
           "CN-Xsh": "1356 mm・20.6 °C（**熱帯**）"}
# これまでに門を通り、型が付いたクラスタ（**旗106/107/111 の実測**）
KNOWN = {"Santa Rita（US-SRM）": "or（雨または秋）",
         "Duolun（CN-Du2）": "rain_only（雨だけ）"}


def main():
    ap = argparse.ArgumentParser(description="旗113：門に決めさせる")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()
    tee_stdout("step113")

    print("=== 旗113：気候の札をやめ、門に決めさせる（残り 3 サイト）===")
    print("  **旗112 の指標は、答えの分かっている CN-Du2 を外した**——**札は使わない。**")
    print("  **3 サイトすべてを走らせ、結果に関わらず全部報告する。**")
    print(f"  **道具は旗109/111 と同一**（`DELTA_FLOOR` = {DELTA_FLOOR}）。")
    print("  **`or` が 1 サイトだけ出ても「再現した」とは書かない**（当てずっぽうで約 58%）。")
    print("  **これで 4 群が揃うサイトの列挙は尽きる。**")

    if not a.real:
        print("\n  【合成データで検証する】**旗109/111 と同じ 5 通り**。")
        print("  **道具を変えていないので結果は同じはず**——**同じであることの確認が目的。**")
        want = {"rain_only": "rain_only（雨だけ）", "season_only": "season_only（季節だけ）",
                "both": "both（雨かつ秋）", "or": "or（雨または秋）", "no_reversal": None}
        got, ok = {}, True
        for k, w in want.items():
            print(f"\n  ===== 合成 `{k}` =====")
            d, P = synth(k)
            got[k] = run_site(f"合成/{k}", d, P)
            ok &= (got[k] == w)
        print("\n  === 合成のまとめ ===")
        for k, w in want.items():
            print(f"    {k:<12}期待 {str(w):<24}実際 {got[k]}  "
                  f"{'✔' if got[k] == w else '**✘**'}")
        print(f"\n  → **5 枝すべて一致：{ok}**"
              f"{'' if ok else '  ← **一致しない枝がある＝実データに進まない**'}")
        return

    res = {}
    for s in SITES:
        print(f"\n  ＜{s}：{CLIMATE[s]}＞")
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max)
            P = daily_precip(s, a.qc_max)
        except Exception as e:
            print(f"    読めない {type(e).__name__}: {str(e)[:90]}")
            res[s] = None; continue
        if P is None or P.dropna().empty:
            print("    **降水 P が無い**"); res[s] = None; continue
        res[s] = run_site(s, d, P)

    passed = {s: v for s, v in res.items() if v}
    print("\n  === 第一段：門①-a を通ったか ===")
    for s in SITES:
        print(f"    {s:<10}{CLIMATE[s]:<28}{res.get(s) or '**通らない**'}")
    if not passed:
        print("\n  **○A-3 の反転は、これ以上手元では見つからない**——**列挙は尽きた。**")
        print("  **手元でクラスタを増やす道は、ここで終わりである。**")
    else:
        print(f"\n  **★門を通った：{list(passed)}**——**独立クラスタがさらに増える。**")
        for s in passed:
            print(f"    **{s} の気候をそのまま書く：{CLIMATE[s]}**")
        if "CN-Xsh" in passed:
            print("    **CN-Xsh（熱帯・1356 mm）が通った**——"
                  "**A-3 の「乾燥地で」という限定が揺らぐ。そう書く。**")

    print("\n  === 第二段：型 ===")
    print("  これまでに型が付いたクラスタ（旗106/107/111）：")
    for k, v in KNOWN.items():
        print(f"    {k:<24}{v}")
    if passed:
        print("  今回：")
        for s, v in passed.items():
            print(f"    {s:<24}{v}")
    allt = Counter(list(KNOWN.values()) + list(passed.values()))
    n_or = sum(v for k, v in allt.items() if k.startswith("or"))
    print(f"\n  **`or` が出たクラスタ数：{n_or}**")
    if n_or >= 2:
        print("  → **○「雨または秋」が 2 クラスタ以上で出た**——**型が繰り返した。**")
        print("  **ただし「2 クラスタで」と書き、「一般に」とは書かない。**")
    else:
        print("  → **○`or` は Santa Rita 固有のまま**（**1 クラスタ**）。")
        print("  **事前登録どおり、1 件では「再現した」と書かない。**")

    print("\n  留保（事前登録どおり）：")
    print("   ・**これで 4 群が揃うサイトの列挙は尽きた。**")
    print("     **手元でクラスタを増やす道は、これが最後である。**")
    print("   ・**CN-Ha2・CN-Aro は寒冷**（−0.9・0.6 °C）。")
    print("     **凍結・融雪が θ の測定と意味に影響しうるが、本研究は扱っていない。**")
    print("   ・**CN-Xsh は熱帯**で、**春・秋という季節区分そのものが北半球中緯度の枠**である。")
    print("     **`SPRING=(3,4,5)`・`AUTUMN=(9,10,11)` は熱帯では意味が違う。**")
    print("   ・**CN-HbC を除いたのは来歴の理由**（`description` が `CN-Erg` を指す・旗110）")
    print("     **であり、結果を見て外したのではない。**")


if __name__ == "__main__":
    main()
