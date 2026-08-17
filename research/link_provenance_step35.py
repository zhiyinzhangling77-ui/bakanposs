"""旗35：因果骨格の各リンクを「本物の因果」か「定義由来」かに仕分ける（変数の独立性で骨格の意味が変わる）。

問い：データの全変数は独立か？→ NO。旗32 で確定した派生構造を、因果骨格の全リンクに適用する。
派生変数(VPD, GER, GEP)を端に持つリンクは、生態系の因果でなく **定義・分割アルゴリズム由来** の
可能性がある。骨格の各リンクを素性で分類し「どれが本物の物理因果でどれが定義の写り込みか」を明示する。

変数の素性（旗32＋文献 Pastorello2020/Lasslop2010/Reichstein2005）:
  独立測定 : Rg, Ta, Ts, θ(th), P, γH(gH), γLE(gLE), NEE
  派生     : VPD = f(Ta,RH)              … Ta の関数
             GER = NEE を Ta,Rg で分割     … NEE/Ta/Rg の関数（分割モデル）
             GEP = NEE を Rg,VPD で分割    … NEE/Rg/VPD の関数（分割モデル, β(VPD)含む）
  測定共有 : γH,γLE,NEE は同じ渦相関系（超音波 w'）＝フラックス間リンクは共通 w' の写り込みがありうる

分類：
  INDEP     … 両端が独立測定＝本物の物理因果の候補（骨格の"岩盤"）
  DEF_IDENT … NEE↔GER↔GEP の炭素恒等式＝定義（因果でない）
  DEF_INPUT … 派生変数と、その"作成に使った入力"の間＝分割/定義の写り込み（Ta→GER, Rg/VPD→GEP, Ta→VPD 等）
  DERIV_PROXY … 派生変数が入力(Ta)を代理して別の独立変数へ＝元の独立リンクの化け（VPD→Ts は Ta→Ts の化け 等）
  EC_SHARED … γH/γLE/NEE 間＝共通 w' の測定共有がありうる

    python research/link_provenance_step35.py                       # 既定=文書化済みJP-Takコアリンク
    python research/link_provenance_step35.py --csv JP-Tak_link_consistency_parcorr.csv --min-freq 0.7
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

INDEP = {"Rg", "Ta", "Ts", "th", "P", "gH", "gLE", "NEE"}
DERIVED_INPUTS = {                    # 派生変数 → その作成に使った入力
    "VPD": {"Ta"},                    # VPD=f(Ta,RH)
    "GER": {"NEE", "Ta", "Rg"},       # RECO 分割（温度応答モデル）
    "GEP": {"NEE", "Rg", "VPD"},      # GPP 分割（光応答＋β(VPD)）
}
CARBON = {"NEE", "GER", "GEP"}        # 炭素恒等式 NEE=GER−GEP で結ばれる三角
FLUXES = {"gH", "gLE", "NEE"}         # 同じ渦相関系（共通 w'）

# ラベル正規化（γH/θ 等の別表記を吸収）
ALIAS = {"γH": "gH", "γLE": "gLE", "θ": "th", "H": "gH", "LE": "gLE", "SWC": "th"}


def _norm(v):
    return ALIAS.get(v, v)


def classify(src, dst):
    """リンク src→dst の素性を返す (category, 理由)。"""
    s, d = _norm(src), _norm(dst)
    # 炭素三角（定義）
    if s in CARBON and d in CARBON:
        return "DEF_IDENT", "炭素恒等式 NEE=GER−GEP の定義（因果でない）"
    # 派生変数と、その作成入力の間
    if d in DERIVED_INPUTS and s in DERIVED_INPUTS[d]:
        return "DEF_INPUT", f"{d} は {s} を入力に作られる（分割/定義の写り込み）"
    if s in DERIVED_INPUTS and d in DERIVED_INPUTS[s]:
        return "DEF_INPUT", f"{s} は {d} を入力に作られる（分割/定義の写り込み）"
    # 派生変数が入力を代理して独立変数へ（VPD→Ts 等 = 元の Ta→Ts の化け）
    if s in DERIVED_INPUTS and d in INDEP:
        proxied = DERIVED_INPUTS[s] & INDEP
        return "DERIV_PROXY", f"{s} は入力{sorted(proxied)}を代理＝{sorted(proxied)}→{d} の化けの疑い"
    if d in DERIVED_INPUTS and s in INDEP and s not in DERIVED_INPUTS[d]:
        return "DERIV_PROXY", f"{d} は入力を含む派生量＝{s}→{d} は入力経由の写り込みの疑い"
    # フラックス間（共通 w'）
    if s in FLUXES and d in FLUXES:
        return "EC_SHARED", "同じ渦相関系（共通 w'）＝測定共有の写り込みがありうる"
    # 両端が独立測定
    if s in INDEP and d in INDEP:
        return "INDEP", "両端が独立測定＝本物の物理因果の候補（骨格の岩盤）"
    return "OTHER", "未分類"


# RESULTS_3LAYERS に文書化された JP-Tak コアリンク（頻度≥70%）
DEFAULT_LINKS = [
    ("GEP", "NEE", 1.00), ("Rg", "gH", 0.95), ("Rg", "gLE", 0.95),
    ("Rg", "VPD", 0.95), ("Ta", "Ts", 0.90), ("Rg", "Ta", 0.86),
    ("Rg", "Ts", 0.76), ("VPD", "Ts", 0.76), ("Ta", "GER", 0.71),
]

CAT_ORDER = ["INDEP", "EC_SHARED", "DERIV_PROXY", "DEF_INPUT", "DEF_IDENT", "OTHER"]
CAT_LABEL = {
    "INDEP": "✅本物の因果候補", "EC_SHARED": "△共通w'の測定共有",
    "DERIV_PROXY": "⚠入力の化け", "DEF_INPUT": "❌分割/定義の写り込み",
    "DEF_IDENT": "❌炭素恒等式(定義)", "OTHER": "?未分類",
}


def load_csv_links(path, min_freq):
    import pandas as pd
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    src = cols.get("src") or cols.get("source") or cols.get("from")
    dst = cols.get("dst") or cols.get("target") or cols.get("to")
    fr = cols.get("frequency") or cols.get("freq") or cols.get("consistency")
    out = []
    for _, r in df.iterrows():
        f = float(r[fr]) if fr else 1.0
        if f >= min_freq:
            out.append((str(r[src]), str(r[dst]), f))
    return sorted(out, key=lambda x: -x[2])


def main():
    p = argparse.ArgumentParser(description="因果骨格リンクを素性で仕分け")
    p.add_argument("--csv", default=None, help="run_robustness の link_consistency CSV")
    p.add_argument("--min-freq", type=float, default=0.70)
    a = p.parse_args()

    if a.csv:
        links = load_csv_links(a.csv, a.min_freq)
        title = f"{Path(a.csv).name}（頻度≥{a.min_freq}）"
    else:
        links = DEFAULT_LINKS
        title = "文書化済み JP-Tak コアリンク（RESULTS_3LAYERS, 頻度≥0.70）"

    print(f"=== 旗35 因果骨格リンクの素性仕分け：{title} ===")
    print("  全変数は独立でない→派生変数を端に持つリンクは定義/分割の写り込みかも。\n")
    print(f"  {'リンク':<16}{'頻度':>6}  {'分類':<12} 理由")
    from collections import Counter
    cnt = Counter()
    for s, d, f in links:
        cat, why = classify(s, d)
        cnt[cat] += 1
        print(f"  {s+'→'+d:<16}{f:>6.2f}  {CAT_LABEL[cat]:<12} {why}")

    print("\n  === まとめ（骨格の意味の再読）===")
    for cat in CAT_ORDER:
        if cnt[cat]:
            print(f"  {CAT_LABEL[cat]:<16} {cnt[cat]} 本")
    ngenuine = cnt["INDEP"]
    print(f"\n  → 本物の物理因果の候補（両端独立）は {ngenuine}/{len(links)} 本。")
    print("    残りは定義・分割・入力の写り込みで、因果骨格の『リンク』として額面通り読めない。")
    print("    ＝『GEP→NEE が普遍』は恒等式で当たり前、意味ある普遍背骨は Rg→γH 等の独立測定間リンク。")


if __name__ == "__main__":
    main()
