"""旗56（策B）：タワー側のメモリを「夜間の実測NEE」だけで測り直す＝派生量を使わない版。

策B の骨子：**派生量（VPD・GEP・GER）と、収支で縛られた組み合わせを主張から外し、
測定量だけで骨格を組み直す**。旗50（背骨は収支の恒等式）・旗46（炭素の冗長は gap-fill 由来）を受けて、
これは選択ではなく必然になった。

だが「メモリ」だけは、タワー側の主張（旗25/37）が **GER＝NEE を分割した派生量**に依存したままだ。
チャンバー（旗40/53/54）は分割を通さないので測定量だが、**タワーで同じことが言えるか**は未確認。

そこで：**夜間の NEE はそれ自体が生態系呼吸**（光合成がゼロなので NEE ≈ RECO）。
分割アルゴリズムも昼のモデル外挿も通さない。さらに **QC=0（実測のみ）**に絞れば穴埋めも通さない。
＝**タワー側のメモリを、測定量だけで検定できる**。

  ・夜間実測NEE でもメモリが出る → タワー側の主張が**派生量依存から解放される**（策Bの核心）。
  ・出ない → タワーのメモリは**分割の産物**だった可能性＝旗39/52 の増幅説と整合し、主張はチャンバーに限定する。

検出器は旗53 の較正済みのもの（非線形基底・ACF1 ≥ 0.64・短メモリ限定）を使う。

    python research/nightly_nee_memory_step56.py --sites JP-Tak JP-Fhk JP-Tmd JP-Fjy JP-Tef
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import _acf_gap, _efold_gap
from memory_attribution_flex_step54 import flex_basis, _fit, ACF_THR, EFOLD_MAX

NIGHT_RG = 20.0        # 夜間判定 (W/m2)：これ以下なら光合成ほぼゼロ＝NEE≈呼吸


def nightly_daily(raw, night_rg=NIGHT_RG):
    """夜間だけを取り出して日次に均す。**その日の夜が実測で埋まっている日のみ**残す。"""
    d = raw[["NEE", "Ts", "th", "Rg"]].copy()
    d = d[d["Rg"] < night_rg]
    d = d[np.isfinite(d["NEE"]) & (d["NEE"] > 0)]        # 夜のNEE>0＝呼吸（負値は乱流不良の疑い）
    if d.empty:
        return None
    g = d.groupby(d.index.normalize())
    daily = g.mean()
    daily["n_night"] = g.size()
    daily = daily[daily["n_night"] >= 6]                 # 1晩あたり最低6点（3時間相当）
    if len(daily) < 60:
        return None
    return daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D"))


def detect(daily, form="flex"):
    y = np.log(daily["NEE"].where(daily["NEE"] > 0)).to_numpy()
    T = daily["Ts"].to_numpy()
    W = daily["th"].to_numpy() if "th" in daily else None
    cols = flex_basis(T, W) if form == "flex" else ([T] + ([W] if W is not None else []))
    res, r2 = _fit(y, cols)
    if res is None:
        return None
    return {"r2": r2, "acf1": _acf_gap(res, 1), "efold": _efold_gap(res),
            "n": int(np.isfinite(res).sum())}


def main():
    p = argparse.ArgumentParser(description="夜間実測NEEだけでタワーのメモリを測る")
    p.add_argument("--sites", nargs="+", required=True)
    p.add_argument("--month", type=int, nargs="+", default=None)
    p.add_argument("--qc-max", type=int, default=0, help="0=実測のみ（既定）")
    a = p.parse_args()

    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.run_robustness import get_site_years
    from japanflux_pn.preprocess import load_raw_all

    cfg = AnalysisConfig(qc_max=a.qc_max)
    print("=== 旗56：夜間の実測NEEだけでタワーのメモリを測る（派生量なし）===")
    print(f"  夜間 Rg<{NIGHT_RG:.0f} W/m² の NEE>0 のみ、QC≤{a.qc_max}（実測のみ）、"
          f"1晩6点以上の日だけ。分割も昼の外挿も通さない。")
    print(f"  判定は旗53 の較正値：非線形基底で ACF1 ≥ {ACF_THR}、短メモリ限定 e-fold ≤ {EFOLD_MAX}日\n")
    print(f"  {'site':<10}{'年数':>5}{'日数':>6}{'線形ACF1':>10}{'線形ef':>8}"
          f"{'柔軟ACF1':>10}{'柔軟ef':>8}{'柔軟R²':>8}  判定")
    n_star = n_ok = n_weak = 0
    for site in a.sites:
        years, mons = get_site_years(site)
        if a.month:
            mons = a.month
        try:
            raw = load_raw_all(get_site(site), cfg)
        except Exception as e:
            print(f"  {site:<10} 読み込み失敗 {type(e).__name__}"); continue
        raw = raw[raw.index.month.isin(mons)]
        daily = nightly_daily(raw)
        if daily is None:
            print(f"  {site:<10} 夜間の実測点が不足（QC≤{a.qc_max} で残らない）"); continue
        lin, flx = detect(daily, "linear"), detect(daily, "flex")
        if not lin or not flx or not np.isfinite(flx["acf1"]):
            print(f"  {site:<10} 推定不能"); continue
        # **適格条件を先に当てる**：R²<0.3 は旗40以来「駆動弱＝判定不能」であって「メモリ無し」ではない。
        # （第1回はこれを怠り、ノイズに支配されたサイトを「·なし」と表示していた＝誤り）
        if not (np.isfinite(flx["r2"]) and flx["r2"] >= 0.3):
            n_weak += 1
            print(f"  {site:<10}{len(years):>5}{len(daily.dropna(subset=['NEE'])):>6}"
                  f"{lin['acf1']:>10.2f}{lin['efold']:>8.0f}{flx['acf1']:>10.2f}"
                  f"{flx['efold']:>8.0f}{flx['r2']:>8.2f}  駆動弱＝**判定不能**", flush=True)
            continue
        n_ok += 1
        star = (flx["acf1"] >= ACF_THR) and (flx["efold"] <= EFOLD_MAX)
        n_star += int(star)
        note = "★メモリあり" if star else (
            "季節(長)" if flx["acf1"] >= ACF_THR else "·短メモリなし")
        print(f"  {site:<10}{len(years):>5}{len(daily.dropna(subset=['NEE'])):>6}"
              f"{lin['acf1']:>10.2f}{lin['efold']:>8.0f}{flx['acf1']:>10.2f}"
              f"{flx['efold']:>8.0f}{flx['r2']:>8.2f}  {note}", flush=True)

    print(f"\n  === まとめ ===")
    print(f"  **駆動弱（R²<0.3）＝判定不能：{n_weak} サイト**"
          "  ← 夜間NEEはノイズが大きく、駆動モデルが説明できない")
    print(f"  判定できた {n_ok} サイトのうち ★短メモリ：{n_star}")
    if n_ok == 0:
        print("\n  → **この検定は成立しなかった**。『メモリが無い』ではなく"
              "『測定量だけでは信号がノイズに埋もれて判定できない』が正しい読み。")
    print("\n  読み：★が出る＝**タワー側のメモリも派生量（分割）なしで確認できる**")
    print("        ＝策Bの核心（測定量だけで骨格を組み直す）が呼吸メモリについて成立する。")
    print("        ★が出ない＝タワーのメモリは分割の産物だった可能性＝主張をチャンバーに限定する。")
    print("  留保：夜間NEEは乱流不良（u*小）で系統的に過小評価されることが知られる。QC=0 で")
    print("        実測に絞ってもこの偏りは残る。またNEE≈呼吸は光合成ゼロの近似であって厳密ではない。")


if __name__ == "__main__":
    main()
