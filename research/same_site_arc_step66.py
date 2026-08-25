"""旗66：**同一地点でタワーとチャンバーを同じ物差しで並べる**＝穴②を閉じる。

旗64/65 で分かったこと：**JP-Fhk ↔ `d20200328_UEYAMA_HOKUROKU` が 0.00 km**で、
チャンバー側にデータがあり **★短メモリ**（R²=0.90, ACF1=+0.75, e-fold=7日）。
JP-Fhk は**閉合最良（EBR 0.85）でSIFも取得済み**＝**タワー・衛星・チャンバーが同一地点に揃う唯一の組**。

これまで穴②（同一地点で3観測系を突き合わせていない）は本研究**最大の弱点**として残っていた。
本ツールは、**同じ検出器（旗53/54 の較正済み・非線形基底）・同じ判定基準**で両側を測り、
**できれば同じ期間**（チャンバーの観測期間にタワーを合わせる）で並べる。

  ・両側で「メモリが在る」と出る → **同一地点で分割の有無に依らず再現**＝弧が閉じる。
  ・片側だけ → **分割（タワー側）か点測定（チャンバー側）に固有**の可能性＝そう記す。

対照として **JP-Tef ↔ `UEYAMA_TESHIO`（0.01km）**も並べる（チャンバー側は季節メモリ e-fold=30日）。
**同じ場所で両側とも「短メモリなし」なら、それも整合の証拠**になる。

    python research/same_site_arc_step66.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore, _acf_gap, _efold_gap
from memory_attribution_flex_step54 import flex_basis, _fit, ACF_THR, EFOLD_MAX

# 旗64/65 で確定した同一地点の対（距離は旗51 の再実行値）
PAIRS = [("JP-Fhk", "d20200328_UEYAMA_HOKUROKU", 0.00),
         ("JP-Tef", "d20200328_UEYAMA_TESHIO", 0.01)]
# 旗38 の実測（記録済み）：季節制御後の SIF 偏相関と、記憶ACFが落ちたか
SIF_NOTE = {"JP-Fhk": "偏相関 +0.03・記憶ACFは落ちず",
            "JP-Tef": "偏相関 +0.09・記憶ACFは落ちず"}


def memory_from_daily(daily, ycol, tcol, wcol):
    """旗53/54 と同じ非線形基底で残差メモリを測る（タワー・チャンバー共通の物差し）。"""
    y = np.log(daily[ycol].where(daily[ycol] > 0)).to_numpy()
    T = daily[tcol].to_numpy()
    W = daily[wcol].to_numpy() if wcol and wcol in daily else None
    res, r2 = _fit(y, flex_basis(T, W))
    if res is None or not np.isfinite(r2):
        return None
    return {"r2": r2, "acf1": _acf_gap(res, 1), "efold": _efold_gap(res),
            "n": int(np.isfinite(res).sum())}


def tower_daily(site, span=None):
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    cfg = AnalysisConfig()
    raw = load_raw_all(get_site(site), cfg)
    d = raw[["GER", "Ts", "th"]].copy()
    if span:
        d = d.loc[(d.index >= span[0]) & (d.index <= span[1])]
    if d.empty:
        return None
    daily = d.groupby(d.index.normalize()).mean()
    return daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D"))


def chamber_daily(path, span=None):
    df, st, sm = load_cosore(path, None)
    if "Tsoil" not in df:
        return None, (st, sm)
    cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
    d = df[cols].copy()
    if span:
        d = d.loc[(d.index >= span[0]) & (d.index <= span[1])]
    if d.empty:
        return None, (st, sm)
    daily = d.groupby(d.index.normalize()).mean()
    return daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D")), (st, sm)


def verdict(m):
    if m is None:
        return "判定不能"
    if not np.isfinite(m["r2"]) or m["r2"] < 0.3:
        return f"駆動弱(R²={m['r2']:.2f})＝判定不能"
    if not np.isfinite(m["acf1"]) or not np.isfinite(m["efold"]):
        return "推定不能"
    if m["acf1"] >= ACF_THR and m["efold"] <= EFOLD_MAX:
        return "★短メモリ"
    if m["acf1"] >= ACF_THR:
        return f"季節メモリ(e-fold={m['efold']:.0f}日)"
    return f"メモリ弱(ACF1={m['acf1']:.2f})"


def main():
    p = argparse.ArgumentParser(description="同一地点でタワーとチャンバーを並べる")
    p.add_argument("--cosore-dir", required=True)
    a = p.parse_args()
    root = Path(a.cosore_dir)

    print("=== 旗66：同一地点でタワー×チャンバーを同じ物差しで並べる（穴②）===")
    print(f"  検出器は旗53/54 の較正済み（非線形基底・R²≥0.3・ACF1≥{ACF_THR}・e-fold≤{EFOLD_MAX}日）。")
    print("  タワー側は GER（分割派生）、チャンバー側は Rs（分割を通さない直接測定）。\n")

    for site, ds, km in PAIRS:
        f = root / "datasets" / f"data_{ds}.csv"
        print(f"  ━━ {site} ↔ {ds}（{km:.2f} km）━━")
        if not f.exists():
            print("    チャンバーのデータファイルが無い\n"); continue
        ch_all, (st, sm) = chamber_daily(f)
        if ch_all is None:
            print("    チャンバー側に土壌温度が無い\n"); continue
        span = (ch_all.index.min(), ch_all.index.max())
        print(f"    チャンバー観測期間：{span[0]:%Y-%m}〜{span[1]:%Y-%m}（列 T={st} SM={sm}）")

        try:
            tw_all = tower_daily(site)
            tw_ov = tower_daily(site, span)
        except Exception as e:
            print(f"    タワー側の読み込み失敗 {type(e).__name__}\n"); continue

        rows = []
        m = memory_from_daily(ch_all, "Rs", "Tsoil", "SM")
        rows.append(("チャンバー Rs（全期間）", m))
        if tw_all is not None:
            rows.append(("タワー GER（全期間）", memory_from_daily(tw_all, "GER", "Ts", "th")))
        if tw_ov is not None and len(tw_ov.dropna(subset=["GER"])) >= 60:
            rows.append(("タワー GER（**同一期間**）", memory_from_daily(tw_ov, "GER", "Ts", "th")))
            rows.append(("チャンバー Rs（同一期間）", m))   # チャンバーは元から同一期間
        else:
            print("    ※タワー側に同一期間のデータが足りず、期間を揃えた比較はできない")

        print(f"    {'系列':<26}{'N':>6}{'R²':>7}{'ACF1':>8}{'e-fold':>8}  判定")
        for lab, mm in rows:
            if mm is None:
                print(f"    {lab:<26}{'—':>6}{'—':>7}{'—':>8}{'—':>8}  判定不能"); continue
            print(f"    {lab:<26}{mm['n']:>6}{mm['r2']:>7.2f}{mm['acf1']:>8.2f}"
                  f"{mm['efold']:>8.0f}  {verdict(mm)}")
        print(f"    衛星SIF（旗38 の記録）：{SIF_NOTE.get(site, '—')}")
        print()

    print("  === 読み方 ===")
    print("  同一地点で**両側とも★短メモリ**＝分割の有無に依らず再現＝**穴②が閉じる**。")
    print("  片側だけ★＝分割（タワー）か点測定（チャンバー）に固有の可能性。")
    print("  両側とも★でない＝その地点にはメモリが無い（それも整合の証拠）。")
    print("  留保：タワーは生態系呼吸（GER, 分割派生）・チャンバーは土壌呼吸（Rs, 点測定）で")
    print("        **測っている対象が違う**。同一地点であっても『同じ量』ではない。")
    print("        チャンバー期間が短い場合、同一期間のタワー側は年数不足で不安定になりうる。")


if __name__ == "__main__":
    main()
