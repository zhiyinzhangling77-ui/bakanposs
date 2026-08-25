"""旗57：中心的主張を、まだ触っていないデータ（非森林サイト）で再現するか検定する。

**予測・判定規則は `research/PREREGISTRATION_step57.md` に実行前確定・commit 済み（1d676ce）。**
本ファイルはその手続きを実装するだけで、結果を見てから基準を動かさない。

  H1（普及率）：非森林の判定可能サイトのうち短メモリ判定の割合 > 0.25（森林の実測 22/45=0.49）
  H2（説明されなさ）：適格かつプラセボ有りのうち候補がどれも説明しない割合 > 0.6（森林 20/22=0.91）
  判定可能サイトが 6 未満なら**判定しない**（検出力不足）。

森林は除外し、それ以外の IGBP をすべて対象にする（草原・低木地・湿地・サバンナ・農地…）。

    python research/replicate_nonforest_step57.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore, _acf_gap, _efold_gap
from chamber_memory_recount_step53 import detect
from memory_attribution_flex_step54 import analyze as attribute, verdict as attr_verdict, ACF_THR, EFOLD_MAX
from memory_attribution_step45 import load_daily

H1_MIN = 0.25          # 事前登録：普及率の下限
H2_MIN = 0.60          # 事前登録：説明されなさの下限
N_MIN = 6              # 事前登録：これ未満なら判定しない


def main():
    p = argparse.ArgumentParser(description="非森林サイトでの事前登録レプリケーション")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--month", type=int, nargs="+", default=None)
    a = p.parse_args()

    root = Path(a.cosore_dir); desc = pd.read_csv(root / "description.csv")
    print("=== 旗57：中心的主張の再現検定（未使用データ＝非森林サイト）===")
    print(f"  事前登録：H1 普及率>{H1_MIN}／H2 説明されなさ>{H2_MIN}／判定可能<{N_MIN} なら判定しない")
    print("  （森林の実測値：普及率 22/45=0.49、説明されなさ 20/22=0.91）\n")
    print(f"  {'dataset':<32}{'IGBP':<16}{'ACF1':>7}{'ef':>5}{'R²':>7}  メモリ判定 / 帰属")

    n_judge = n_short = 0
    n_attr = n_unexplained = 0
    igbp_seen = {}
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        if "forest" in ig.lower():                 # 森林は既使用なので除外
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        # --- メモリ判定（旗53 と同じ手続き） ---
        try:
            df, st, sm = load_cosore(f, a.month)
            if "Tsoil" not in df:
                continue
            cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
            dd = df[cols].copy()
            daily = dd.groupby(dd.index.normalize()).mean()
            daily = daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D"))
            flx = detect(daily, "Rs", "Tsoil", "SM" if "SM" in df else None, "flex")
        except Exception:
            continue
        if not flx or not np.isfinite(flx.get("efold", np.nan)) or not np.isfinite(flx["r2"]):
            continue
        if flx["r2"] < 0.3:                        # 駆動弱＝判定不能（旗56 の反省）
            continue
        n_judge += 1
        igbp_seen[ig] = igbp_seen.get(ig, 0) + 1
        short = (flx["acf1"] >= ACF_THR) and (flx["efold"] <= EFOLD_MAX)
        n_short += int(short)

        # --- 帰属（旗54 と同じ手続き）：短メモリのサイトだけ ---
        attr = "—"
        if short:
            try:
                dly, meta = load_daily(f, a.month)
                r = attribute(dly)
                v, best = attr_verdict(r)
                attr = v
                if "note" not in r and "プラセボ無し" not in v:
                    n_attr += 1
                    n_unexplained += int(best is None)
            except Exception:
                attr = "（帰属の推定に失敗）"
        print(f"  {ds:<32}{ig[:14]:<16}{flx['acf1']:>7.2f}{flx['efold']:>5.0f}{flx['r2']:>7.2f}"
              f"  {'★短メモリ' if short else '·'} / {attr}", flush=True)

    print(f"\n  === 結果 ===")
    print(f"  判定可能（R²≥0.3）：{n_judge} サイト  内訳 {igbp_seen}")
    if n_judge < N_MIN:
        print(f"  → **事前登録の通り判定しない**（判定可能 {n_judge} < {N_MIN}＝検出力不足）。")
        print("     非森林は土壌センサが疎という、事前に認めた限界がそのまま出た形。")
        return
    p1 = n_short / n_judge
    print(f"  H1 普及率：{n_short}/{n_judge} = {p1:.2f}  （事前登録の閾値 >{H1_MIN}／森林 0.49）")
    if n_attr:
        p2 = n_unexplained / n_attr
        print(f"  H2 説明されなさ：{n_unexplained}/{n_attr} = {p2:.2f}"
              f"  （事前登録の閾値 >{H2_MIN}／森林 0.91）")
    else:
        p2 = np.nan
        print("  H2：帰属を評価できたサイトなし（短メモリ判定かつプラセボ有りが無い）")

    ok1, ok2 = p1 > H1_MIN, (np.isfinite(p2) and p2 > H2_MIN)
    print("\n  === 事前登録の判定規則に照らす ===")
    if ok1 and ok2:
        print("  → ★**再現**：中心的主張は森林に限らない＝**生態系タイプを越える現象**として述べてよい。")
    elif ok1 or ok2:
        which = "H1（普及率）" if ok1 else "H2（説明されなさ）"
        print(f"  → ○**部分再現**：{which} のみ満たす。満たした側だけ一般化し、他方は**森林に限定**して述べる。")
    else:
        print("  → ▲**再現せず**：中心的主張を**「森林で」と明示的に限定**する。")
        print("     「観測の隙間」も森林の話に絞る。")
    print("\n  留保（事前登録で先に認めた通り）：非森林はセンサが疎／撹乱の影響で同じ機構を期待する")
    print("  根拠は強くない＝満たされなくても『主張が誤り』ではなく『適用範囲が森林』を意味しうる。")
    print("  森林の値そのものが予測の出所であり、当たっても機構の証明ではない。")


if __name__ == "__main__":
    main()
