"""旗101：**COSORE に何層の土壌温度・水分があるか**（下調べ・**検定はしない**）。

## なぜこれを見るのか

**旗100 で「タワーの気象では説明されない」と出た**（判定4組中1組のみ説明された）。
**残る候補は「地点規模の生物過程」か「土壌水文」。**
**土壌水文なら、深い層の水分がチャンバー駆動から抜けている可能性がある。**

**そして——本研究は 40 旗にわたって、一貫して浅い 1 層しか使っていない。**
`cosore_memory_step40._pick_sm` は **`CSR_SM<深さ>` のうち 5 cm に最も近い層**を選び、
`_pick_soil_temp` も **5 cm に最も近い層**を選ぶ。**それ以外の層は一度も読んでいない。**

**旗97 で FLUXNET の `SWC_F_MDS_1..8` について同じ見落としをしていた**
（「多深度の水分は新規観測が要る」と書いたが、既に手元にあった）。
**COSORE 側でも同じことが起きていないかを確かめる。**

## **この道具は実行可能性だけを出す。記憶も相関も計算しない。**

**それは検定の答えそのものであり、事前登録の前に見てはいけない**（旗94/98 と同じ作法）。
出すのは**層の数・深さ・有効日数・年数**だけである。

    python research/cosore_depths_probe_step101.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from same_site_arc_step66 import PAIRS, SENSITIVITY_ONLY
from runlog import tee_stdout      # **出力を最初からファイルに残す**（旗110 の反省）

MIN_DAYS, MIN_YEARS = 60, 3


def depth_cols(cols, prefix):
    """`CSR_SM12.5` のような列を **(深さ, 列名)** で拾う。"""
    out = []
    for c in cols:
        m = re.fullmatch(rf"CSR_{prefix}(\d+\.?\d*)", c)
        if m:
            out.append((float(m.group(1)), c))
    return sorted(out)


def summarize(path):
    df = pd.read_csv(path, low_memory=False)
    cols = list(df.columns)
    if "CSR_FLUX_CO2" not in cols:
        return None
    tcol = "CSR_TIMESTAMP_BEGIN" if "CSR_TIMESTAMP_BEGIN" in cols else "CSR_TIMESTAMP_END"
    ts = pd.to_datetime(df[tcol], errors="coerce")
    day = ts.dt.normalize()
    rec = {"sm": [], "t": []}
    for key, pref in (("sm", "SM"), ("t", "T")):
        for d, c in depth_cols(cols, pref):
            v = pd.to_numeric(df[c], errors="coerce")
            ok = v.notna() & day.notna()
            if not ok.any():
                rec[key].append((d, c, 0, 0)); continue
            dd = day[ok]
            rec[key].append((d, c, dd.nunique(), dd.dt.year.nunique()))
    return rec


def main():
    ap = argparse.ArgumentParser(description="COSORE の層数を数える（検定はしない）")
    ap.add_argument("--cosore-dir", default="/mnt/hdd/cosore-0.7.0")
    ap.add_argument("--all", action="store_true",
                    help="旗86 の 44 組だけでなく全データセットを見る")
    a = ap.parse_args()

    tee_stdout("step101")
    print("=== 旗101：COSORE に何層の土壌温度・水分があるか（下調べ・検定はしない）===")
    print("  **本研究は 40 旗にわたって、5 cm に最も近い 1 層しか使っていない。**")
    print("  **旗97 で FLUXNET の多深度水分を見落としていたのと同じことが、")
    print("  COSORE 側でも起きていないかを確かめる。**")
    print("  **記憶も相関も計算しない**——**それは検定の答えである。**")
    print(f"  下限は**有効 {MIN_DAYS} 日・{MIN_YEARS} 暦年**。\n")

    root = Path(a.cosore_dir)
    if a.all:
        targets = sorted(p.stem[5:] for p in (root / "datasets").glob("data_*.csv"))
    else:
        targets = sorted({ds for site, ds, km in PAIRS
                          if (site, ds) not in SENSITIVITY_ONLY})

    n_multi_sm = n_multi_t = n_read = 0
    multi_list = []
    for ds in targets:
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            rec = summarize(f)
        except Exception as e:
            print(f"  {ds:<34}**読めない**（{type(e).__name__}: {str(e)[:50]}）")
            continue
        if rec is None:
            print(f"  {ds:<34}CSR_FLUX_CO2 が無い"); continue
        n_read += 1
        sm_ok = [r for r in rec["sm"] if r[2] >= MIN_DAYS and r[3] >= MIN_YEARS]
        t_ok = [r for r in rec["t"] if r[2] >= MIN_DAYS and r[3] >= MIN_YEARS]
        tag = ""
        if len(sm_ok) >= 2:
            n_multi_sm += 1; tag += " **水分が多層**"
        if len(t_ok) >= 2:
            n_multi_t += 1; tag += " **温度が多層**"
        print(f"  ━━ {ds} ━━{tag}")
        if rec["sm"]:
            print("    水分：" + " ／ ".join(
                f"{d:g}cm {n:,}日/{y}年" + ("" if n >= MIN_DAYS and y >= MIN_YEARS else "（下限未満）")
                for d, c, n, y in rec["sm"]))
        else:
            print("    水分：**無し**")
        if rec["t"]:
            print("    温度：" + " ／ ".join(
                f"{d:g}cm {n:,}日/{y}年" + ("" if n >= MIN_DAYS and y >= MIN_YEARS else "（下限未満）")
                for d, c, n, y in rec["t"]))
        else:
            print("    温度：**無し**")
        if len(sm_ok) >= 2 or len(t_ok) >= 2:
            multi_list.append((ds, [f"{d:g}" for d, c, n, y in sm_ok],
                               [f"{d:g}" for d, c, n, y in t_ok]))

    print("\n  === まとめ ===")
    print(f"  読めたデータセット：{n_read}")
    print(f"  **下限を満たす水分層が 2 層以上：{n_multi_sm} 本**")
    print(f"  **下限を満たす温度層が 2 層以上：{n_multi_t} 本**")
    if multi_list:
        print("\n  多層のデータセット（**使える層の深さ**）：")
        for ds, sm, t in multi_list:
            print(f"    {ds:<34}水分 [{', '.join(sm) or '—'}] cm／温度 [{', '.join(t) or '—'}] cm")

    print("\n  === 次の判断（**事前登録の前に決める**）===")
    print("  ・**多層の水分が 3 本以上** → **旗102 を事前登録**して、")
    print("    **深い層で駆動を組み直すと ~4 日メモリが消えるか**を検定する")
    print("  ・**多層の水分が 3 本未満** → **手 E は実行できない**と記して打ち切り、")
    print("    **「土壌水文か生物過程か」は手元では割れない**と結論して新規観測へ渡す")
    print("  ・**多層の温度が 3 本以上** → **別枠の一手**（**熱の浸透深度**）が立つ。")
    print("    **本研究は温度も 5 cm 1 層しか使っていない。**")
    print("\n  留保：")
    print("   ・**深さの表記は列名に依存する**（`CSR_SM<深さcm>`）。")
    print("     **実際の設置深度が列名どおりかは、我々には確かめられない**（旗97 と同じ限界）。")
    print("   ・**層が多いことと、層が独立なことは別**——")
    print("     **深い層が浅い層とほぼ同じ動きなら、足しても何も増えない。**")
    print("     **その確認は旗102 の設計で受ける**（**ここでは相関を見ない**）。")


if __name__ == "__main__":
    main()
