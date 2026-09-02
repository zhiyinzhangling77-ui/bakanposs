"""旗110：**4 群が揃った 9 サイトは、乾燥地か**（下調べ・**検定はしない**）。

## なぜこれが要るのか

**旗108 を直した道具で走らせ直したら、4 群すべてが揃うサイトが 4 → 9 件に増えた**
（**CN-Aro・CN-Du2・CN-Ha2・CN-HbC・CN-Xsh・US-Ho1・US-NC2・US-SRM・US-WCr**）。
**うち 5 件は ChinaFlux で、旗108 の第1版では数えてすらいなかった。**

**だが「その 5 件が乾燥地か」を、私は知らない。**
**サイト名から推測してはいけない**——**それは研究不正の一歩手前である。**
**A-3 は「乾燥地で水がマスター変数」という主張なので、
母集団が乾燥地かどうかで、検定の意味がまったく変わる。**

## **この道具は、乾燥度を「データから」出す。名前で決めない。**

**出すのは測定量だけ**：
  ・**年降水量**（`P` の合計 ÷ 年数）——**測定量そのもの**
  ・**θ の中央値と四分位**——**測定量そのもの**
  ・**`P` の季節分布**（春・夏・秋・冬の割合）——**「秋に雨が降る気候か」が分かる**
  ・**`SiteSpec.description`**（ChinaFlux 発見が拾った生態系フォルダ名など）
    ——**メタデータであって測定ではない。参考として出し、判定には使わない。**

**既知の乾燥地 3 サイト（US-Wkg・US-Whs・US-SRM）を必ず並べる**
——**絶対値ではなく、手元の乾燥地と比べて判断するため。**

**相関も Δ も計算しない**（旗94/98/101/103/108 と同じ作法）。

    python research/aridity_probe_step110.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stratified_bowen_step89 import cell_of
from soiltemp_match_step90 import SPRING, AUTUMN
from rain_history_probe_step103 import daily_precip
from evaporation_regime_step36 import daily_energy

# 旗108（直した版）で 4 群すべてが揃った 9 件
NEW = ("CN-Aro", "CN-Du2", "CN-Ha2", "CN-HbC", "CN-Xsh",
       "US-Ho1", "US-NC2", "US-WCr")
REF = ("US-Wkg", "US-Whs", "US-SRM")     # **手元の既知の乾燥地**（物差しとして並べる）
SUMMER, WINTER = (6, 7, 8), (12, 1, 2)


def profile(site, qc_max):
    d, _ = daily_energy(site, list(range(1, 13)), qc_max)
    P = daily_precip(site, qc_max)
    if P is None or P.dropna().empty:
        return None
    yrs = P.index.year.nunique()
    ann = float(P.sum()) / max(yrs, 1)
    seas = {}
    tot = float(P.sum())
    for nm, mon in (("春", SPRING), ("夏", SUMMER), ("秋", AUTUMN), ("冬", WINTER)):
        s = float(P[[m in mon for m in P.index.month]].sum())
        seas[nm] = s / tot if tot > 0 else np.nan
    lab, tmed, rmed = cell_of(d)
    th = d["th"].to_numpy()
    from japanflux_pn.sites import get_site
    desc = getattr(get_site(site), "description", "") or ""
    return {"ann": ann, "yrs": yrs, "seas": seas, "tmed": tmed,
            "th_q": (np.nanpercentile(th, 25), np.nanmedian(th),
                     np.nanpercentile(th, 75)),
            "cell": int((lab == "θ高×Rg高").sum()), "desc": desc[:60]}


def main():
    ap = argparse.ArgumentParser(description="旗110：9 サイトは乾燥地か（検定はしない）")
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()

    print("=== 旗110：4 群が揃った 9 サイトは、乾燥地か（下調べ・検定はしない）===")
    print("  **サイト名から推測しない。** **年降水量と θ をデータから出して比べる。**")
    print("  **既知の乾燥地 3 件を物差しとして並べる**（絶対値ではなく相対で見る）。")
    print("  **相関も Δ も計算しない**——**それは検定の答えである。**\n")

    rows = {}
    for s in list(REF) + list(NEW):
        try:
            p = profile(s, a.qc_max)
        except Exception as e:
            print(f"  {s:<10}**読めない**（{type(e).__name__}: {str(e)[:60]}）")
            continue
        if p is None:
            print(f"  {s:<10}**降水 P が無い**"); continue
        rows[s] = p

    print(f"\n  {'サイト':<10}{'年降水':>9}{'年数':>5}   "
          f"{'θ 中央 [25–75%]':<24}{'春':>6}{'夏':>6}{'秋':>6}{'冬':>6}  セル")
    for s in list(REF) + list(NEW):
        if s not in rows:
            continue
        p = rows[s]
        q = p["th_q"]
        tag = "  ← **既知の乾燥地**" if s in REF else ""
        print(f"  {s:<10}{p['ann']:>7.0f}mm{p['yrs']:>5}   "
              f"{q[1]:>6.2f} [{q[0]:.2f}–{q[2]:.2f}]      "
              f"{p['seas']['春']:>6.0%}{p['seas']['夏']:>6.0%}"
              f"{p['seas']['秋']:>6.0%}{p['seas']['冬']:>6.0%}{p['cell']:>6}{tag}")

    print("\n  参考（**メタデータであって測定ではない。判定には使わない**）：")
    for s in list(REF) + list(NEW):
        if s in rows and rows[s]["desc"]:
            print(f"    {s:<10}{rows[s]['desc']}")

    ref_ann = [rows[s]["ann"] for s in REF if s in rows]
    if ref_ann:
        hi = max(ref_ann)
        print(f"\n  === 目安 ===")
        print(f"  **既知の乾燥地 3 件の年降水は {min(ref_ann):.0f}–{hi:.0f} mm。**")
        near = [s for s in NEW if s in rows and rows[s]["ann"] <= hi * 1.5]
        print(f"  **その 1.5 倍（{hi*1.5:.0f} mm）以下の新サイト：{near}**")
        print("  **これは目安であって定義ではない**——**乾燥度は降水だけでは決まらない**")
        print("  （**蒸発要求＝放射・気温・VPD が要る**）。**次の事前登録で扱いを決める。**")

    print("\n  === 次の判断（**事前登録の前に決める**）===")
    print("  ・**乾燥地に近い新サイトが 1 件以上**（**既存 2 クラスタから 50 km 以上**）")
    print("    → **旗111 を事前登録**して、**「雨または秋」の型を同じ母集団で検定する**")
    print("      （**旗109 と違い、これは『再現』の検定になる**）")
    print("  ・**乾燥地に近い新サイトが 0 件** → **旗109 の結論のまま**")
    print("    （**手元の乾燥地クラスタは 2**）**と記して新規観測へ渡す**")
    print("\n  留保：")
    print("   ・**年降水量は乾燥度の一面でしかない。** **同じ 300 mm でも、")
    print("     蒸発要求が違えば水制限の強さは違う。**")
    print("   ・**`description` はメタデータ**——**提供者が書いた文字列であり、測定ではない。**")
    print("     **判定には使わない。** **食い違ったら測定量を採る。**")
    print("   ・**相関も Δ も一度も計算していない**（事前登録の前に答えを見ない）。")


if __name__ == "__main__":
    main()
