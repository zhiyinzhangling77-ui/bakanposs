"""旗76：同一地点ペアで **降水と気圧が使えるか**を確かめる（検定はしない・下調べ）。

旗45 と旗75 は、結論のたびに同じ限定を書いてきた——**「降水は COSORE に無い」**。
だが**同一地点のペアが4組ある**（JP-Fhk 0.00km／JP-Tef 0.01km／JP-Yms 0.03・0.69km、旗66/67）。
タワー側は `P_F`（降水）を持つ。＝**この4組については、手元のデータで限定を外せるはず**である。

同時に、**物理側の対抗仮説**も安く潰せるかもしれない：
**気圧の変動はチャンバー測定に既知のアーティファクト（pressure pumping）を生む**。
気圧列がタワーに在れば、**新規観測なしで物理仮説の一部を検定できる**。

本ツールは**下調べだけ**を行う（旗68/69 と同じ作法＝**推測で進めず、実ファイルを見てから決める**）：

  1. 各タワーで **P（降水）が読めるか**、読めるなら**量と頻度**はどうか
  2. **気圧に相当する列が実ファイルに在るか**（PA / PRESS / ATM 等を総当たり）
  3. **チャンバーの観測期間と重なる日数**——ここが足りなければ検定できない

**重なりが十分にあり、降水が実際に降っていて、気圧列が在る**——この3つが揃った組だけが、
次の検定（旗77）に進む価値を持つ。揃わなければ「**手元では無理**」と確定して記録する。

    python research/precip_pressure_probe_step76.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore
from same_site_arc_step66 import PAIRS

# 気圧らしき列を探すためのトークン。**"P" 単独は降水と紛れる**ので使わない
# （`sites.py` にも「"P" は気圧なので使わない」という別形式向けの注意がある＝規約が形式ごとに違う）。
PRESS_TOKENS = ("PA_", "PA)", "PRESS", "ATM", "BARO", "HPA", "KPA")


def tower_columns(code, data_dir):
    """そのサイトの最大 csv の列名（気圧列を探すため）。"""
    root = Path(data_dir)
    csvs = sorted([p for p in root.rglob(f"*{code}*")
                   if p.is_file() and p.suffix.lower() == ".csv"
                   and "__MACOSX" not in p.parts and not p.name.startswith(("._", "~$"))],
                  key=lambda p: p.stat().st_size, reverse=True)
    for f in csvs[:3]:
        try:
            return list(pd.read_csv(f, nrows=2).columns), f.name
        except Exception:
            continue
    return None, None


def tower_daily_p(code):
    """タワーの日降水量（mm/日）。P が無ければ None。"""
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    raw = load_raw_all(get_site(code), AnalysisConfig())
    if "P" not in raw.columns:
        return None
    # 30分値の合計＝日降水量
    return raw["P"].groupby(raw.index.normalize()).sum()


def main():
    p = argparse.ArgumentParser(description="同一地点ペアで降水と気圧が使えるか")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--data-dir", default="/mnt/hdd/JAPANFLUX")
    a = p.parse_args()
    root = Path(a.cosore_dir)

    print("=== 旗76：同一地点ペアで降水と気圧が使えるか（下調べ・検定はしない）===")
    print("  旗45/75 は毎回『降水は COSORE に無い』と限定してきた。")
    print("  だが**同一地点ペアのタワー側は降水を持つ**＝この4組なら限定を外せるはず。")
    print("  併せて**気圧**（pressure pumping ＝物理側の対抗仮説）が在るかも見る。\n")

    done = set()
    for code, ds, km in PAIRS:
        f = root / "datasets" / f"data_{ds}.csv"
        print(f"  ━ {code} ↔ {ds}（{km:.2f} km）━")
        if not f.exists():
            print("    チャンバーのデータが無い\n"); continue
        try:
            df, st, sm = load_cosore(f, None)
        except Exception as e:
            print(f"    チャンバー読み込み失敗 {type(e).__name__}: {str(e)[:90]}\n"); continue
        if df is None or "Rs" not in df:
            print("    チャンバー側に Rs が無い\n"); continue
        ch = df["Rs"].groupby(df.index.normalize()).mean().dropna()
        if ch.empty:
            print("    チャンバー側の日次が空\n"); continue
        span = (ch.index.min(), ch.index.max())
        print(f"    チャンバー観測期間：{span[0]:%Y-%m-%d}〜{span[1]:%Y-%m-%d}（{len(ch)} 日）")

        # ① 降水
        try:
            pser = tower_daily_p(code)
        except Exception as e:
            print(f"    タワー読み込み失敗 {type(e).__name__}: {str(e)[:90]}"); pser = None
        if pser is None:
            print("    **タワーに降水 P が無い**")
        else:
            ov = pser.loc[(pser.index >= span[0]) & (pser.index <= span[1])].dropna()
            both = ch.index.intersection(ov.index)
            if len(ov) == 0:
                print("    **チャンバー期間に重なる降水データが無い**")
            else:
                wet = float((ov > 1.0).mean())
                print(f"    降水 P：重なり **{len(both)} 日**／期間内 {len(ov)} 日"
                      f"・年あたり {ov.sum()/max(len(ov)/365.25,1e-9):.0f} mm"
                      f"・**>1mm の日が {wet:.0%}**・最大 {ov.max():.0f} mm/日")
                if len(both) < 120:
                    print("      ※重なりが 120 日未満＝**検定には足りない可能性が高い**")

        # ② 気圧（同じサイトを2回見ない）
        if code not in done:
            done.add(code)
            cols, fname = tower_columns(code, a.data_dir)
            if cols is None:
                print("    気圧：列を読めるファイルが無い")
            else:
                cand = [c for c in cols if any(t in c.upper() for t in PRESS_TOKENS)]
                if cand:
                    print(f"    **気圧の候補列**：{cand[:6]}（{fname}）")
                else:
                    print(f"    **気圧に相当する列は無い**（{len(cols)} 列を走査・{fname}）")
        print()

    print("  === 次の判断 ===")
    print("  **降水の重なりが十分な組**があれば、旗77 として")
    print("  『降水を駆動に加えてもメモリが残るか』を検定できる＝**旗45/75 の最大の限定が外れる**。")
    print("  **気圧の列がある**なら、pressure pumping（物理側の対抗仮説）も")
    print("  **新規観測なしで**検定できる。")
    print("  どちらも無ければ『**手元では無理**』と確定し、そう記録する。")
    print("  留保：")
    print("   ・タワーの降水は**タワー位置での観測**であり、チャンバーの直上ではない。")
    print("     0.00〜0.69 km なので対流性の強雨では差が出うる（この限定は残る）。")
    print("   ・気圧列があっても、**単位（hPa/kPa）と欠測の扱い**は使う前に確かめること。")
    print("   ・降水は 30 分値の合計を日量とした。**元が強度(mm/h)なら合計は誤り**——")
    print("     値域（年 mm）が常識的かで気づけるので、上の年あたり mm を必ず見ること。")


if __name__ == "__main__":
    main()
