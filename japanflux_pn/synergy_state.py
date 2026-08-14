"""旗15：高次×状態依存の接続。各サブ系の相乗/冗長（O-info z）が、その年の気候状態
（生の平均 VPD・土壌水分）で組み変わるかを、年を環境にした Spearman 相関で見る。

旗14：呼吸系(Rg,Ta,θ,GER)は相乗支配（z<0）が広く見られる。もし相乗の源が「温度×水分の
相互作用」なら、乾湿でその強さが変わるはず＝相乗z が気候状態の関数（旗13の状態依存の高次版）。

本体 `oinfo_analysis.scan_years`（年ごと各サブ系の z）と `climate_response.scan`
（年ごとの生 VPD/土壌水分/気温）を年で突き合わせ、サブ系ごとに
  Spearman r( z , 状態 )  ＋ 順列検定 p
を出す。相乗（z<0）が乾燥（高VPD/低θ）で強まるなら r(z,VPD)<0 / r(z,θ)>0。

    python -m japanflux_pn.synergy_state --site JP-Tak
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig
from .oinfo_analysis import scan_years, SUBSYSTEMS, OINFO_BINS_DEFAULT
from .climate_response import scan as climate_scan, _spearman, _spearman_p

STATES = [("VPD_mean", "VPD", "乾燥↑"), ("SWC_mean", "土壌水分", "湿潤↑"),
          ("Ta_mean", "気温", "高温↑")]


def build(site: str, obins: int, config: AnalysisConfig) -> pd.DataFrame:
    """年 × (各サブ系 z, 状態指標) の突き合わせ表を返す。"""
    long = scan_years(site, obins, config)            # year, subsystem, Omega, z
    if long.empty:
        return pd.DataFrame()
    zwide = long.pivot_table(index="year", columns="subsystem", values="z")
    clim = climate_scan(site, config)                 # year, VPD_mean, SWC_mean, Ta_mean...
    if clim.empty:
        return pd.DataFrame()
    clim = clim.set_index("year")[["VPD_mean", "SWC_mean", "Ta_mean"]]
    df = zwide.join(clim, how="inner")
    return df.reset_index()


def report(site: str, obins: int = OINFO_BINS_DEFAULT,
           config: AnalysisConfig | None = None, out_csv=None) -> pd.DataFrame:
    config = config or AnalysisConfig()
    print(f"===== {site} 高次構造の状態依存（相乗/冗長 z vs 気候状態）=====")
    df = build(site, obins, config)
    if len(df) < 4:
        print(f"  有効年 {len(df)} < 4 → 相関評価に不足。")
        return df
    print(f"  有効年 {len(df)}。z<0=相乗 / z>0=冗長。乾燥で相乗が強まるなら r(z,VPD)<0。\n")
    print(f"  {'サブシステム':<22} {'平均z':>7} {'vs VPD':>8} {'vs 土壌水分':>10} {'vs 気温':>8}")
    subs = [s for s in SUBSYSTEMS if s in df.columns]
    rows = []
    for sub in subs:
        z = df[sub].to_numpy(dtype=float)
        rV = _spearman(z, df["VPD_mean"].to_numpy())
        rS = _spearman(z, df["SWC_mean"].to_numpy())
        rT = _spearman(z, df["Ta_mean"].to_numpy())
        print(f"  {sub:<22} {np.nanmean(z):7.1f} {rV:8.2f} {rS:10.2f} {rT:8.2f}")
        rows.append({"subsystem": sub, "mean_z": float(np.nanmean(z)),
                     "r_vs_VPD": rV, "r_vs_SWC": rS, "r_vs_Ta": rT})

    # 目玉：呼吸系の相乗が乾燥で強まるか（順列検定 p つき）
    resp = next((s for s in subs if s.startswith("呼吸")), None)
    if resp:
        z = df[resp].to_numpy(dtype=float)
        rV, pV = _spearman_p(z, df["VPD_mean"].to_numpy())
        rS, pS = _spearman_p(z, df["SWC_mean"].to_numpy())
        print(f"\n  [仮説検定] 呼吸系の相乗 z vs 乾燥（順列 5000, 有効 {len(df)} 年）")
        print(f"    VPD↑ で z↓(相乗↑) を期待: r={rV:+.2f}  p={pV:.3f}")
        print(f"    土壌水分↑ で z↑(相乗↓) を期待: r={rS:+.2f}  p={pS:.3f}")
        strong = (rV <= -0.4 and pV < 0.05) or (rS >= 0.4 and pS < 0.05)
        trend = (rV <= -0.4 or rS >= 0.4)
        verd = ("✅ 支持・有意＝呼吸の相乗は乾燥で強まる（高次構造が状態依存）"
                if strong else ("△ 傾向あり・非有意（年数/効果が不足）" if trend
                                else "× 明確な状態依存は見えない（相乗の強さは気候状態に依らず）"))
        print(f"    → {verd}")

    print("\n  意味: 相乗z が状態で動く＝『どの高次構造が現れるか』自体が気候で組み変わる")
    print("        ＝旗13（結合『強度』の状態依存）を高次（相乗/冗長）へ拡張。固定モデルでは書けない層。")

    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_csv, index=False)
        print(f"\n[saved] {out_csv}")
    return df


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="高次構造(相乗/冗長)の状態依存")
    p.add_argument("--site", required=True)
    p.add_argument("--obins", type=int, default=OINFO_BINS_DEFAULT)
    p.add_argument("--out-csv", default=None)
    a = p.parse_args(argv)
    report(a.site, a.obins, out_csv=a.out_csv)


if __name__ == "__main__":
    main()
