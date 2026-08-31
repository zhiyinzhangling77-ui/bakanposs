"""旗80：旗79 で登録した AmeriFlux サイトが**本当に読めるか**を確かめる（検定はしない）。

登録は済んだが、**二つの仮定が入っている**——

  1. **分割法**：`RECO_DT_VUT_REF`（昼分割）を使った。日本のサイトが `RECO_DT_vUT` なので
     **揃える**ためだが、**AmeriFlux の FLUXNET 版に DT があるかは確認していない**
     （旗79 の候補列表示は 8 件で切れており、NT しか見えていなかった）。
     無ければ **NT に切り替える**——だがそれは**旗29/30/37/39 が扱った「分割法による違い」を
     サイト間に持ち込む**ので、そう明記する必要がある。
  2. **深度**：`TS_F_MDS_1` / `SWC_F_MDS_1`（**_1 が最浅**という慣例）を使った。
     **列名から実際の深度は分からない**＝旗33 の「θ 深度不統一」がここでも残る。

本ツールは各サイトについて、**写像した 11 列が実ヘッダに在るか**を1つずつ照合し、
**無い列を名指しする**。そのうえで実際に読み込んで**行数と期間**を出す。
**ここを通らないものは旗66 に進めない。**

    python research/verify_registration_step80.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

NEW = ["US-SSH", "US-Ha1", "US-MMS", "CA-TP4", "CA-TP3", "CA-TPD",
       "US-Wkg", "US-Whs", "US-SRM",
       # --- 旗86：★で選ばずに取得した 17 本 ---------------------------------
       "BR-Sa3", "CA-Ca1", "CL-SDF", "US-Bi1", "US-Bi2", "US-Me6",
       "US-NC2", "US-NC4", "US-SRS", "US-Ton", "US-Tw3", "US-Uaf",
       "US-Var", "US-WCr", "US-xSE",
       # 旧 FULLSET 配布（個別 BADM 無し）＝座標は旗79 修正版で確認する
       "US-Ho1", "US-UMB"]


def main():
    from japanflux_pn.sites import get_site, RK_VARS
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.preprocess import find_corevars_files, load_raw_all

    print("=== 旗80：登録した AmeriFlux サイトが本当に読めるか（検定はしない）===")
    print("  仮定は二つ：**分割法 DT**（日本側と揃えるため）と**深度 _1**（最浅の慣例）。")
    print("  どちらも**列名から確認できない**ので、**在るかどうかだけ**をここで確かめる。\n")

    ok_sites = []
    for code in NEW:
        print(f"  ━━ {code} ━━")
        try:
            sp = get_site(code)
        except Exception as e:
            print(f"    登録が引けない {type(e).__name__}: {e}\n"); continue
        try:
            files = find_corevars_files(sp)
        except Exception as e:
            print(f"    HH ファイルが見つからない {type(e).__name__}: {str(e)[:120]}\n"); continue
        print(f"    HH ファイル {len(files)} 件：{files[0].name}")
        try:
            head = pd.read_csv(files[0], nrows=2)
        except Exception as e:
            print(f"    ヘッダを読めない {type(e).__name__}: {str(e)[:120]}\n"); continue
        cols = set(head.columns)
        vm = sp.var_map()
        miss = {k: v for k, v in vm.items() if v not in cols}
        if miss:
            print(f"    **写像した列が無い**：{miss}")
            # 代案を探す（DT が無ければ NT、_1 が無ければ他の層）
            for k, v in miss.items():
                alt = [c for c in head.columns
                       if c.split("_")[0] == v.split("_")[0] and not c.endswith("_QC")]
                print(f"       {k}（{v}）の代案：{alt[:6]}")
        else:
            print(f"    **11 列すべて在る**")
        try:
            raw = load_raw_all(sp, AnalysisConfig())
        except Exception as e:
            print(f"    読み込み失敗 {type(e).__name__}: {str(e)[:140]}\n"); continue
        have = [v for v in RK_VARS if v in raw.columns]
        n_ok = {v: int(raw[v].notna().sum()) for v in ("GER", "Ts", "th") if v in raw}
        print(f"    読み込み成功：{len(raw):,} 行／{raw.index.min():%Y-%m}〜{raw.index.max():%Y-%m}")
        print(f"    変数 {len(have)}/{len(RK_VARS)}／有効数 {n_ok}")
        if all(n_ok.get(v, 0) > 1000 for v in ("GER", "Ts", "th")):
            ok_sites.append(code)
            print(f"    → **旗66 に進める**")
        else:
            print(f"    → **GER/Ts/th のどれかが足りない＝旗66 では使えない**")
        print()

    print("  === まとめ ===")
    print(f"  旗66 に進めるサイト：{len(ok_sites)}/{len(NEW)} {ok_sites}")
    print("\n  留保：")
    print("   ・**DT が無く NT に切り替えた場合**、日本のサイト（DT）との比較に")
    print("     **分割法の違いが混ざる**（旗39：NT の記憶は 8 中 7 で DT より長い）。")
    print("     ＝その場合はサイト間の比較ではなく**同一地点内のタワー対チャンバー**に限って読むこと。")
    print("   ・**深度は列名から分からない**。チャンバー側（COSORE の CSR_T<depth>）と")
    print("     **合っている保証は無い**＝旗33 の指摘がそのまま残る。")


if __name__ == "__main__":
    main()
