"""旗94：**PhenoCam 取得物の下調べ**（登録の前・検定はしない）。

旗93 で「**この観測構成では測れない**」と確定したので観測を足した。
**下調べを通してから事前登録（旗95）を書く。順序を変えない**（旗68/76/79/83/87 と同じ作法）。

**この道具は四つだけ答える。判定はしない。**

  1. **何が置かれているか**——サイト・ROI・列名・日付の範囲
  2. **フラックスと何年重なるか**（**3 年未満なら、そのサイトは判定しない**）
  3. **θ高×Rg高 の春に GCC が何日あるか**（**下限 60 日を超えるか**）
  4. **緑と水が分離できるか**——**緑の日と枯れた日で θ の重なり帯が作れるか**
     （**作れなければ旗93 と同じ壁**＝「水と緑は切り離せない」と確定して終える）

**4 がいちばん大事である。** 旗93 は**検出力**で落ちたが、**この壁は検出力では解けない**
——**GCC を日次で使えても、緑の日と枯れた日で θ が重ならなければ、比べる土俵が無い。**

    python research/phenocam_probe_step94.py --phenocam-dir /mnt/hdd/PhenoCam
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
from evaporation_regime_step36 import daily_energy
from stratified_bowen_step89 import cell_of, MIN_DAYS, MIN_YEARS
from soiltemp_match_step90 import band, SPRING, PLO, PHI

# PhenoCam のサイト名 → 対応する AmeriFlux タワー。**推測を混ぜない**——
# 対応は**座標で確かめるべき**だが、PhenoCam 側の座標が要約ファイルに在るとは限らない。
# **在れば確かめ、無ければ「名前で対応させた・未確認」と明記する**（旗51/79 と同じ注意）。
NAME_MAP = {"kendall": "US-Wkg", "luckyhills": "US-Whs", "santarita": "US-SRM"}
# GCC らしき列（PhenoCam の要約ファイルは版で列名が変わりうる）
GCC_PAT = re.compile(r"^(gcc|smooth_gcc|midday_gcc)", re.IGNORECASE)
DATE_PAT = re.compile(r"^(date|midday_date)$", re.IGNORECASE)


def read_summary(path):
    """PhenoCam 要約 csv を読む。**`#` で始まるコメント行を飛ばす**。"""
    try:
        df = pd.read_csv(path, comment="#", low_memory=False)
    except Exception as e:
        return None, f"{type(e).__name__}: {str(e)[:70]}"
    dcol = next((c for c in df.columns if DATE_PAT.match(str(c))), None)
    if dcol is None:
        return None, f"日付列が無い（列：{list(df.columns)[:10]}）"
    gcols = [c for c in df.columns if GCC_PAT.match(str(c))]
    if not gcols:
        return None, f"GCC らしき列が無い（列：{list(df.columns)[:14]}）"
    out = pd.DataFrame(index=pd.to_datetime(df[dcol], errors="coerce"))
    for c in gcols:
        out[c] = pd.to_numeric(df[c], errors="coerce").to_numpy()
    return out[out.index.notna()].sort_index(), None


def scan(root):
    """サイト名 → [(ROI ラベル, パス)]。**ファイル名から推測せず、中身も見る。**"""
    per = {}
    for p in sorted(Path(root).rglob("*.csv")):
        low = p.name.lower()
        site = next((k for k in NAME_MAP if k in low or k in str(p.parent).lower()), None)
        per.setdefault(site or f"**不明（{p.parent.name}）**", []).append(p)
    return per


def main():
    ap = argparse.ArgumentParser(description="PhenoCam 取得物の下調べ")
    ap.add_argument("--phenocam-dir", required=True)
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()

    print("=== 旗94：PhenoCam 取得物の下調べ（検定はしない）===")
    print("  **下調べを通してから事前登録（旗95）を書く。順序を変えない。**")
    print("  **いちばん大事なのは 4 番**——**緑の日と枯れた日で θ の重なり帯が作れるか**。")
    print("  **作れなければ、日次 GCC を得ても旗93 と同じ壁に当たる。**\n")

    per = scan(a.phenocam_dir)
    if not per:
        print(f"  {a.phenocam_dir} に csv が無い"); return
    print(f"  見つかったサイト：{len(per)} 件 {sorted(per)}\n")

    for site in sorted(per):
        files = per[site]
        amf = NAME_MAP.get(site)
        print(f"  ━━ {site}"
              f"{f' → {amf}（**名前で対応させた・座標未確認**）' if amf else ''} ━━")
        print(f"    csv {len(files)} 件")
        rois = {}
        for f in files:
            g, err = read_summary(f)
            if g is None:
                print(f"    ※{f.name}：**読めない**（{err}）"); continue
            best = max(g.columns, key=lambda c: g[c].notna().sum())
            n_ok = int(g[best].notna().sum())
            print(f"    {f.name}")
            print(f"       列 {list(g.columns)}／採る列 **{best}**（有効 {n_ok:,}）"
                  f"／期間 {g.index.min():%Y-%m-%d}〜{g.index.max():%Y-%m-%d}"
                  f"／年数 {g.index.year.nunique()}")
            rois[f.name] = (g, best)
        if not rois or amf is None:
            if amf is None:
                print("    **対応する AmeriFlux タワーが分からない**＝以降の照合はできない\n")
            else:
                print("    **読める ROI が無い**\n")
            continue

        # ── フラックス側と突き合わせる ──
        try:
            d, _ = daily_energy(amf, list(range(1, 13)), a.qc_max, extra=("Ts",))
        except Exception as e:
            print(f"    タワー側を読めない {type(e).__name__}: {str(e)[:90]}\n"); continue
        lab, tmed, rmed = cell_of(d)
        hh = d[lab == "θ高×Rg高"]
        sp = hh[[m in SPRING for m in hh.index.month]]
        print(f"    タワー {amf}：{len(d):,} 日／{d.index.year.min()}–{d.index.year.max()}"
              f"／θ高×Rg高 の春 {len(sp)} 日")

        for name, (g, col) in rois.items():
            s = g[col].dropna()
            yrs_ov = sorted(set(s.index.year) & set(sp.index.year))
            j = sp.join(s.rename("gcc"), how="inner").dropna(subset=["gcc"])
            print(f"    ── {name}（{col}）──")
            print(f"       **重なる年 {len(yrs_ov)}**"
                  f"{'（' + str(yrs_ov[0]) + '–' + str(yrs_ov[-1]) + '）' if yrs_ov else ''}"
                  f"／**θ高×Rg高 の春で GCC も在る日 {len(j)}**")
            if len(yrs_ov) < MIN_YEARS:
                print(f"       → **重なる年が {MIN_YEARS} 未満＝判定しない**"); continue
            if len(j) < 2 * MIN_DAYS:
                print(f"       → **日数が {2*MIN_DAYS} 未満＝2 群に割ると下限を割る**"); continue
            # ── 4番：**緑と水を分離できるか** ──
            gm = float(j["gcc"].median())
            green = j[j["gcc"] >= gm]; brown = j[j["gcc"] < gm]
            lo, hi = band(green["th"].to_numpy(), brown["th"].to_numpy())
            print(f"       GCC 中央値 {gm:.4f}／緑 {len(green)} 日・枯 {len(brown)} 日")
            print(f"       θ 中央値：緑 {green['th'].median():.3f}／枯 {brown['th'].median():.3f}")
            if hi <= lo:
                print(f"       → **θ の重なり帯が作れない**（[{lo:.3f}, {hi:.3f}]）"
                      f"＝**緑と水が完全に分離＝旗93 と同じ壁**")
                continue
            gb = green[(green["th"] >= lo) & (green["th"] <= hi)]
            bb = brown[(brown["th"] >= lo) & (brown["th"] <= hi)]
            ok = len(gb) >= MIN_DAYS and len(bb) >= MIN_DAYS
            print(f"       **θ の重なり帯 [{lo:.3f}, {hi:.3f}]**（幅 {hi-lo:.3f}）"
                  f"／帯の中：緑 {len(gb)} 日・枯 {len(bb)} 日")
            print(f"       → {'**検定できる見込み**（両群とも下限 ' + str(MIN_DAYS) + ' 日以上）' if ok else '**帯に絞ると下限を割る＝検定できない**'}")
        print()

    print("  === 次の判断 ===")
    print("  ・**帯が作れて両群とも 60 日以上**のサイトが **2 つ以上** → **旗95 を事前登録して検定する**")
    print("  ・**1 つ以下** → **独立クラスタが足りない**と記して、**検定しない**")
    print("  ・**どのサイトでも帯が作れない** → **緑と水は切り離せないと確定**して打ち切る")
    print("  留保：")
    print("   ・**PhenoCam とタワーの対応は名前で付けた・座標未確認**（旗51/79 と同じ注意）。")
    print("     **要約ファイルに座標が在れば確かめること。**")
    print("   ・**カメラの視野とフラックスのフェッチは同じではない**（旗81 と同型の問題）。")
    print("   ・**GCC は色であって光合成ではない。**")
    print("   ・**ROI が複数あるとき、どれを使うかは人が決める**——**この道具は選ばない**。")


if __name__ == "__main__":
    main()
