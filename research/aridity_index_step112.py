"""旗112：**乾燥度を、恣意的な線ではなく測定量から決める**（下調べ・**検定はしない**）。

## なぜこれが要るのか

**旗110 で私が引いた線は「既知の乾燥地の最大 354 mm × 1.5 ＝ 531 mm」だった。**
**1.5 に根拠は無い。** **私が決めただけである。**
**その線一本で CN-Ha2（570 mm）と CN-Aro（618 mm）が落ち、CN-Du2（468 mm）だけが残った。**

**いま「線を 620 mm に緩める」ことはしない。**
**旗111 の結果を見た後で線を動かすのは禁じ手である**（旗93 で下限を緩めなかったのと同じ）。

**代わりに、線そのものを作り直す。**

## **どう作り直すか**

**降水量だけでは乾燥度は決まらない**——**同じ 500 mm でも、暑ければ乾く。**
**蒸発要求を入れた指標を、手元の測定量だけで作る**：

    乾燥指標 = 年降水量 P [mm] ／（年平均気温 Ta [°C] + 10）

**`P` も `Ta` も、独立に測っている 8 変数のうちの 2 つ**である（旗32）。
**どちらも本検定の結果変数（γLE・γH）ではない**——**循環しない。**
**`Rg` は使わない**（**セルの定義と偏相関の統制に入っているので、近すぎる**）。

## **しきい値は、文献からではなく「手元のサイト」から決める**

**この形の指標には既知の分類（乾燥・半乾燥…）があるが、私はその一次を確認していない。**
**未確認の閾値を持ち込まない**（旗49/92 の作法）。

**代わりに、手元で答えの分かっているサイトで較正する**：
  ・**乾燥側**：US-Wkg・US-Whs・US-SRM（**旗82/106 で反転が出ている**）
  ・**湿潤側**：US-Ho1・US-NC2・US-WCr（**旗109 で θ→γLE ≈0＝水支配が無い**）

**この 6 つが指標で分離しなければ、指標を使わない。**
**分離したら、境目（乾燥側の最大と湿潤側の最小の中点）を線とする。**

**＝線は「私の勘」ではなく「手元のデータが示す分かれ目」になる。**

## **この道具は検定をしない**

**CN-Ha2・CN-Aro で Bowen 反転を測ったことは一度も無い。**
**本道具でも測らない。** **出すのは乾燥指標だけである。**
**どのサイトを次の検定に入れるかを、結果を見る前に決めるための道具である。**

    python research/aridity_index_step112.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rain_history_probe_step103 import daily_precip
from runlog import tee_stdout

DRY = ("US-Wkg", "US-Whs", "US-SRM")        # 反転が出ている（旗82/106）
WET = ("US-Ho1", "US-NC2", "US-WCr")        # θ→γLE ≈0＝水支配が無い（旗109）
TEST = ("CN-Du2", "CN-Ha2", "CN-Aro", "CN-HbC", "CN-Xsh")   # 判定したい 5 件


def annual_ta(site, qc_max):
    """**年平均気温**。`daily_energy` の `dropna` を通さずに読む（**日が落ちると偏る**）。"""
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    cfg = AnalysisConfig(qc_max=qc_max) if qc_max is not None else AnalysisConfig()
    raw = load_raw_all(get_site(site), cfg)
    if "Ta" not in raw.columns:
        return None, 0
    t = raw["Ta"].dropna()
    if t.empty:
        return None, 0
    return float(t.mean()), int(t.index.year.nunique())


def index_of(site, qc_max):
    P = daily_precip(site, qc_max)
    if P is None or P.dropna().empty:
        return None
    yrs_p = P.index.year.nunique()
    ann_p = float(P.sum()) / max(yrs_p, 1)
    ta, yrs_t = annual_ta(site, qc_max)
    if ta is None:
        return {"P": ann_p, "Ta": None, "idx": None, "yrs": (yrs_p, 0)}
    denom = ta + 10.0
    idx = ann_p / denom if denom > 0 else np.nan     # **Ta ≤ −10 °C では定義できない**
    return {"P": ann_p, "Ta": ta, "idx": idx, "yrs": (yrs_p, yrs_t)}


def main():
    ap = argparse.ArgumentParser(description="旗112：乾燥度を測定量から決める（検定はしない）")
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()
    tee_stdout("step112")

    print("=== 旗112：乾燥度を、恣意的な線ではなく測定量から決める（下調べ・検定はしない）===")
    print("  **旗110 の線（531 mm）は私が決めただけで根拠が無い。**")
    print("  **旗111 の結果を見た後で線を緩めるのは禁じ手**（旗93）。**線そのものを作り直す。**")
    print("  **指標 ＝ 年降水量 P ／（年平均気温 Ta + 10）**")
    print("  **`P` も `Ta` も独立測定の 8 変数のうち**で、**本検定の結果変数ではない**（循環しない）。")
    print("  **`Rg` は使わない**——**セルの定義と統制に入っていて近すぎる。**")
    print("  **しきい値は文献ではなく手元のサイトから決める**（未確認の閾値を持ち込まない）。")
    print("  **CN-Ha2・CN-Aro で反転を測ったことは一度も無い。本道具でも測らない。**\n")

    rows = {}
    for s in DRY + WET + TEST:
        try:
            r = index_of(s, a.qc_max)
        except Exception as e:
            print(f"  {s:<10}**読めない**（{type(e).__name__}: {str(e)[:60]}）"); continue
        if r is None:
            print(f"  {s:<10}**降水 P が無い**"); continue
        rows[s] = r

    print(f"  {'サイト':<10}{'年降水':>9}{'年平均気温':>11}{'指標 P/(Ta+10)':>16}   群")
    for grp, names in (("**乾燥側（反転あり）**", DRY), ("**湿潤側（水支配なし）**", WET),
                       ("判定したい", TEST)):
        for s in names:
            if s not in rows:
                continue
            r = rows[s]
            ta = f"{r['Ta']:>9.1f}°C" if r["Ta"] is not None else "     **無い**"
            ix = f"{r['idx']:>14.1f}" if r["idx"] is not None and np.isfinite(r["idx"]) \
                else "        **出せない**"
            print(f"  {s:<10}{r['P']:>7.0f}mm{ta}{ix}   {grp}")

    dry_i = [rows[s]["idx"] for s in DRY if s in rows and rows[s]["idx"] is not None
             and np.isfinite(rows[s]["idx"])]
    wet_i = [rows[s]["idx"] for s in WET if s in rows and rows[s]["idx"] is not None
             and np.isfinite(rows[s]["idx"])]

    print("\n  === 較正：指標は、答えの分かっている 6 サイトを分離するか ===")
    if len(dry_i) < 2 or len(wet_i) < 2:
        print("  **較正できない**（どちらかの群が 2 未満）。**指標は使わない。**")
        return
    print(f"  乾燥側 {len(dry_i)} 件：{min(dry_i):.1f}–{max(dry_i):.1f}")
    print(f"  湿潤側 {len(wet_i)} 件：{min(wet_i):.1f}–{max(wet_i):.1f}")
    if max(dry_i) >= min(wet_i):
        print("  → **分離しない**（重なっている）。")
        print("  **この指標では乾燥度を切れない。** **使わない。**")
        print("  **旗110 の線（531 mm）のまま、CN-Du2 だけで進む**と記す。")
        return
    cut = (max(dry_i) + min(wet_i)) / 2
    print(f"  → **分離した。** **境目（中点）＝ {cut:.1f}**")
    print("  **この線は「私の勘」ではなく「手元のデータが示す分かれ目」である。**")

    print("\n  === 判定したい 5 件を、この線に当てる ===")
    add = []
    for s in TEST:
        if s not in rows or rows[s]["idx"] is None or not np.isfinite(rows[s]["idx"]):
            print(f"    {s:<10}**指標を出せない**"); continue
        v = rows[s]["idx"]
        side = "**乾燥側**" if v < cut else "湿潤側"
        print(f"    {s:<10}{v:>6.1f}  → {side}")
        if v < cut:
            add.append(s)
    print(f"\n  **乾燥側に入る：{add}**")

    print("\n  === 次の判断（**事前登録の前に決める**）===")
    print("  ・**CN-Du2 以外にも乾燥側の新サイトがある** → **旗113 を事前登録**して、")
    print("    **そのサイトも同じ道具（旗109/111）で走らせる**")
    print("    （**乾燥地クラスタが増える＝型の比較が n=2 より強くなる**）")
    print("  ・**CN-Du2 だけ** → **旗111 の結論のまま**。**新規観測へ渡す。**")
    print("\n  留保：")
    print("   ・**この指標も一つの指標にすぎない。** **蒸発要求を Ta だけで代表している。**")
    print("   ・**しきい値は手元の 6 サイトから決めた**——**別のサイト群なら別の線になる。**")
    print("   ・**指標が乾燥側と言っても、そのサイトで反転が起きるとは限らない。**")
    print("     **それは旗113 の門①-a で受ける。**")
    print("   ・**CN-HbC は説明文が別サイト（CN-Erg）を指している**（旗110）。")
    print("     **指標が乾燥側でも、フォルダの中身を確かめるまで検定に入れない。**")
    print("   ・**反転も Δ も一度も計算していない**（事前登録の前に答えを見ない）。")


if __name__ == "__main__":
    main()
