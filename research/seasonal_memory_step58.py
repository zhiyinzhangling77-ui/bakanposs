"""旗58：~4日メモリは休眠期にも残るか＝残った候補（生物系 vs 物理系）を分ける最初の判別子。

**予測・手続き・判定規則は `research/PREREGISTRATION_step58.md` に実行前確定・commit 済み（bf49f5b）。**
本ファイルはその手続きを実装するだけで、結果を見てから基準を動かさない。

  H1：生育期に短メモリ（ACF1≥0.64 かつ e-fold≤7日）を示すサイトのうち、
      休眠期にも ACF1≥0.64 を示す割合 > 0.5
  成立 → メモリは現在の植物活動に依存しない（物理過程 or 季節を越える遅いプール）
  不成立 → 生育期の生物活動に結びつく（遅い基質・微生物・根の動態）
  両季節そろって判定可能なサイトが 5 未満なら**判定しない**。

    python research/seasonal_memory_step58.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore, _acf_gap, _efold_gap
from memory_attribution_flex_step54 import flex_basis, _fit, ACF_THR, EFOLD_MAX

MIN_DAYS = 60          # 事前登録：各季節の最低日数
MIN_R2 = 0.3           # 事前登録：各季節の最低当てはまり
MIN_AMP = 8.0          # 事前登録：季節振幅(日次Tsoilの90%点−10%点)がこれ未満なら休眠期なしとして除外
N_MIN = 5              # 事前登録：両季節そろって判定可能がこれ未満なら判定しない


def season_split(daily):
    """日次土壌温度の上位40%＝生育期／下位40%＝休眠期（中間20%は捨てる）。"""
    T = daily["Tsoil"]
    lo, hi = T.quantile(0.40), T.quantile(0.60)
    return daily[T >= hi], daily[T <= lo]


def memory_of(sub, full_index):
    """季節部分集合の残差メモリ。**連続グリッドに戻して**ギャップ対応ACFで測る。"""
    if len(sub.dropna(subset=["Rs", "Tsoil"])) < MIN_DAYS:
        return None
    y = np.log(sub["Rs"].where(sub["Rs"] > 0)).to_numpy()
    T = sub["Tsoil"].to_numpy()
    W = sub["SM"].to_numpy() if "SM" in sub else None
    res, r2 = _fit(y, flex_basis(T, W))
    if res is None or not np.isfinite(r2):
        return None
    # 元の日付グリッド上に残差を戻す（季節抽出で日が飛ぶため、ギャップ対応ACFで扱う）
    full = pd.Series(np.nan, index=full_index)
    full.loc[sub.index] = res
    arr = full.to_numpy()
    return {"r2": float(r2), "acf1": _acf_gap(arr, 1), "efold": _efold_gap(arr),
            "n": int(np.isfinite(res).sum())}


def main():
    p = argparse.ArgumentParser(description="休眠期にもメモリが残るか（事前登録 H1）")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--igbp", default="forest")
    a = p.parse_args()

    root = Path(a.cosore_dir); desc = pd.read_csv(root / "description.csv")
    print("=== 旗58：~4日メモリは休眠期にも残るか（事前登録 bf49f5b）===")
    print(f"  季節：日次Tsoil の上位40%＝生育期／下位40%＝休眠期。振幅<{MIN_AMP}℃ のサイトは除外。")
    print(f"  各季節 {MIN_DAYS}日以上・R²≥{MIN_R2} のみ判定。H1：生育期に短メモリのうち休眠期にも"
          f" ACF1≥{ACF_THR} が >0.5\n")
    print(f"  {'dataset':<32}{'振幅℃':>7}{'生ACF1':>8}{'生ef':>6}{'生R²':>7}"
          f"{'休ACF1':>8}{'休ef':>6}{'休R²':>7}  判定")

    n_both = n_grow_short = n_dorm_keep = 0
    skipped = {"振幅小": 0, "季節の点不足": 0}
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        if a.igbp and a.igbp.lower() not in ig.lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            df, st, sm = load_cosore(f, None)          # 通年を使う
            if "Tsoil" not in df:
                continue
            cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
            dd = df[cols].copy()
            daily = dd.groupby(dd.index.normalize()).mean()
            grid = pd.date_range(daily.index.min(), daily.index.max(), freq="D")
            daily = daily.reindex(grid)
            amp = float(daily["Tsoil"].quantile(0.90) - daily["Tsoil"].quantile(0.10))
            if not np.isfinite(amp) or amp < MIN_AMP:
                skipped["振幅小"] += 1
                continue
            grow, dorm = season_split(daily.dropna(subset=["Tsoil"]))
            g, h = memory_of(grow, grid), memory_of(dorm, grid)
        except Exception:
            continue
        if g is None or h is None or g["r2"] < MIN_R2 or h["r2"] < MIN_R2 \
           or not np.isfinite(g["acf1"]) or not np.isfinite(h["acf1"]):
            skipped["季節の点不足"] += 1
            continue
        n_both += 1
        grow_short = (g["acf1"] >= ACF_THR) and (g["efold"] <= EFOLD_MAX)
        dorm_keep = h["acf1"] >= ACF_THR
        n_grow_short += int(grow_short)
        if grow_short:
            n_dorm_keep += int(dorm_keep)
        note = ("生育期★ → 休眠期" + ("も★" if dorm_keep else "は消える")) if grow_short else "生育期に短メモリなし"
        print(f"  {ds:<32}{amp:>7.1f}{g['acf1']:>8.2f}{g['efold']:>6.0f}{g['r2']:>7.2f}"
              f"{h['acf1']:>8.2f}{h['efold']:>6.0f}{h['r2']:>7.2f}  {note}", flush=True)

    print(f"\n  === 結果 ===")
    print(f"  両季節そろって判定可能：{n_both} サイト"
          f"（除外：振幅小 {skipped['振幅小']}／季節の点不足・当てはまり不足 {skipped['季節の点不足']}）")
    if n_both < N_MIN:
        print(f"  → **事前登録の通り判定しない**（{n_both} < {N_MIN}＝検出力不足）。")
        print("     休眠期は呼吸が小さくS/Nが悪化するという、事前に認めた危険がそのまま出た形。")
        return
    print(f"  生育期に短メモリ：{n_grow_short}／{n_both}")
    if n_grow_short == 0:
        print("  → 生育期に短メモリのサイトが無く、H1 を評価できない。")
        return
    frac = n_dorm_keep / n_grow_short
    print(f"  そのうち休眠期にも ACF1≥{ACF_THR}：{n_dorm_keep}／{n_grow_short} = {frac:.2f}"
          f"  （事前登録の閾値 >0.5）")
    print("\n  === 事前登録の判定規則に照らす ===")
    if frac > 0.5:
        print("  → **H1 成立**：メモリは**現在の植物活動に依存しない**＝物理過程、または")
        print("     季節を越える遅いプールを示唆。Stoy 2007 の落葉期観察と整合する向き。")
    else:
        print("  → **H1 不成立**：メモリは**生育期の生物活動に結びつく**＝遅い基質・微生物・")
        print("     根の動態を支持。物理拡散だけでは説明しにくい向き。")
    print("\n  留保（事前登録で先に認めた通り）：休眠期は呼吸が小さくS/Nが悪化する。判定不能が多ければ")
    print("  それは『メモリが無い』ではなく『見えない』。温度分位の季節定義はフェノロジーとの一致を仮定。")
    print("  季節で駆動の分布が変わるため、季節間のACF1比較にはその分の不確かさが乗る。")
    print("  本検定は機構を**絞る**もので**特定する**ものではない。")


if __name__ == "__main__":
    main()
