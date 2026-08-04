#!/usr/bin/env python3
"""FLUXNET形式の半時間(HH)CSVを読み、欠測(-9999)処理・QCフィルタ・主要変数抽出をする。

対応: FLUXNET2015 / ICOS(Warm Winter 2020) / AmeriFlux-FLUXNET / FLUXNET-CH4
(いずれも同じ FLUXNET 形式なので同じコードで読める)。

使い方:
  pip install pandas
  python read_fluxnet.py FLX_XX-Yyy_..._FULLSET_HH_....csv [out_clean.csv]

QC: FLUXNET のフラグは 0=実測, 1=良質ギャップ充填, 2=中, 3=不良。
    既定は QC<=1(実測+良質)以外を欠測(NaN)にする。
"""

from __future__ import annotations

import sys

try:
    import pandas as pd
except Exception:
    sys.exit("pandas が必要です: pip install pandas")

# 使う主要変数(存在するものだけ抽出)。VUT_REF = 変動 u* しきい値の代表値。
KEEP = [
    "TIMESTAMP_START", "TIMESTAMP_END",
    "NEE_VUT_REF", "NEE_VUT_REF_QC",
    "GPP_NT_VUT_REF", "GPP_DT_VUT_REF",
    "RECO_NT_VUT_REF", "RECO_DT_VUT_REF",
    "LE_F_MDS", "LE_F_MDS_QC", "H_F_MDS", "H_F_MDS_QC",
    "NETRAD", "SW_IN_F", "SW_IN_F_QC",
    "VPD_F", "VPD_F_QC", "TA_F", "TA_F_QC",
    "SWC_F_MDS_1", "SWC_F_MDS_1_QC", "TS_F_MDS_1",
    "P_F", "P_F_QC",
    "FCH4", "FCH4_F", "FCH4_F_QC",   # FLUXNET-CH4
]

# 値列 → 対応する QC 列(この QC がしきい値超なら値を NaN にする)
QC_PAIRS = {
    "NEE_VUT_REF": "NEE_VUT_REF_QC",
    "LE_F_MDS": "LE_F_MDS_QC",
    "H_F_MDS": "H_F_MDS_QC",
    "SW_IN_F": "SW_IN_F_QC",
    "VPD_F": "VPD_F_QC",
    "TA_F": "TA_F_QC",
    "SWC_F_MDS_1": "SWC_F_MDS_1_QC",
    "P_F": "P_F_QC",
    "FCH4_F": "FCH4_F_QC",
}


def load_fluxnet(path: str, qc_max: int = 1) -> "pd.DataFrame":
    df = pd.read_csv(path, na_values=["-9999", "-9999.0", -9999, -9999.0])
    cols = [c for c in KEEP if c in df.columns]
    df = df[cols].copy()

    # タイムスタンプ(YYYYMMDDHHMM)を datetime に
    for t in ("TIMESTAMP_START", "TIMESTAMP_END"):
        if t in df.columns:
            df[t] = pd.to_datetime(
                df[t].astype("Int64").astype("string"),
                format="%Y%m%d%H%M", errors="coerce",
            )
    if "TIMESTAMP_START" in df.columns:
        df = df.set_index("TIMESTAMP_START")

    # QC フィルタ(qc_max を超える品質は欠測に)
    for val, qc in QC_PAIRS.items():
        if val in df.columns and qc in df.columns:
            df.loc[df[qc] > qc_max, val] = pd.NA

    return df


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    df = load_fluxnet(sys.argv[1])
    print(f"行数: {len(df)}  期間: {df.index.min()} 〜 {df.index.max()}")
    print("列:", list(df.columns))
    # 各変数の有効データ率(QC後)
    print("\n有効データ率(QC後):")
    for c in df.columns:
        if c.endswith("_QC"):
            continue
        frac = df[c].notna().mean() if len(df) else 0
        print(f"  {c:20s} {frac*100:5.1f}%")
    if len(sys.argv) > 2:
        df.to_csv(sys.argv[2])
        print(f"\n保存: {sys.argv[2]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
