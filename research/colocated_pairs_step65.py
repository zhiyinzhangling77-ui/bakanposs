"""旗65：同一地点のタワー×チャンバー対を точ検し、穴②を閉じられるか確かめる。

旗64 で座標抽出の誤りを直したところ、**同一地点の対が5組**見つかった：
  JP-Tak ↔ `d20200109_KISHIMOTO-MO`（**0.03 km**）／JP-Fhk ↔ UEYAMA_HOKUROKU（0.00 km）／
  JP-Tef ↔ UEYAMA_TESHIO（0.01 km）ほか。

とくに **KISHIMOTO-MO は、我々が「COSORE v0.7 に未同梱だから連絡が要る」と判断していた
まさにその高山サイトのチャンバー**である（`CONTACT_DRAFT_TKY.md` の前提が崩れた）。

だが **旗40/53/54 の解析にこのサイトは一度も現れていない**。＝どこかで落ちている。
本ツールはその理由を突き止め、使えるなら同一地点解析に進めるかを判定する：
  1. データが何を持っているか（Rs・土壌温度の深度・土壌水分・期間・点数）
  2. 較正済み検出器（非線形基底・R²≥0.3・ACF1≥0.64・e-fold≤7日）にかけると何が出るか
  3. 落ちている場合、**どの条件で落ちたか**を明示する（欠測／深度なし／点数不足／当てはまり不足）

    python research/colocated_pairs_step65.py --cosore-dir /mnt/hdd/cosore-0.7.0 \
        --datasets d20200109_KISHIMOTO-MO d20200328_UEYAMA_HOKUROKU d20200328_UEYAMA_TESHIO
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore
from chamber_memory_recount_step53 import detect
from memory_attribution_flex_step54 import ACF_THR, EFOLD_MAX


def inspect(path):
    """生CSVの中身を点検し、なぜ解析に乗らなかったかを説明できる形で返す。"""
    raw = pd.read_csv(path, low_memory=False)
    cols = list(raw.columns)
    tdep = sorted(float(m.group(1)) for c in cols
                  if (m := re.fullmatch(r"CSR_T(\d+\.?\d*)", c)))
    sdep = sorted(float(m.group(1)) for c in cols
                  if (m := re.fullmatch(r"CSR_SM(\d+\.?\d*)", c)))
    info = {"n_rows": len(raw), "T_depths": tdep, "SM_depths": sdep,
            "has_flux": "CSR_FLUX_CO2" in cols,
            "has_tair": any(c in cols for c in ("CSR_TAIR", "CSR_TAIR_AMB"))}
    return info, cols


def main():
    p = argparse.ArgumentParser(description="同一地点対の点検")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--datasets", nargs="+", required=True)
    a = p.parse_args()
    root = Path(a.cosore_dir)

    print("=== 旗65：同一地点のタワー×チャンバー対を点検する ===")
    print("  目的：これらが使えるなら、**穴②（同一地点で3観測系を突き合わせていない）が閉じる**。\n")
    for ds in a.datasets:
        f = root / "datasets" / f"data_{ds}.csv"
        print(f"  --- {ds} ---")
        if not f.exists():
            print("    ファイルが無い\n"); continue
        try:
            info, cols = inspect(f)
        except Exception as e:
            print(f"    読み込み失敗 {type(e).__name__}: {e}\n"); continue
        print(f"    行数 {info['n_rows']}／Rs列 {'あり' if info['has_flux'] else '**なし**'}"
              f"／土壌温度 深度 {info['T_depths'] or '**なし**'}"
              f"／土壌水分 深度 {info['SM_depths'] or 'なし'}"
              f"／気温 {'あり' if info['has_tair'] else 'なし'}")
        if not info["has_flux"]:
            print("    → **Rs が無いので解析対象外**\n"); continue

        try:
            df, st, sm = load_cosore(f, None)
        except Exception as e:
            print(f"    load_cosore 失敗 {type(e).__name__}\n"); continue
        if "Tsoil" not in df:
            print("    → **土壌温度も気温も拾えない＝旗40/53 の対象から外れる**（これが不在の理由）\n")
            continue
        span = f"{df.index.min():%Y-%m}〜{df.index.max():%Y-%m}"
        cols2 = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
        dd = df[cols2].copy()
        daily = dd.groupby(dd.index.normalize()).mean()
        daily = daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D"))
        n_day = int(daily["Rs"].notna().sum())
        flx = detect(daily, "Rs", "Tsoil", "SM" if "SM" in df else None, "flex")
        print(f"    期間 {span}／日数(Rs有効) {n_day}／使用列 T={st} SM={sm}")
        if not flx:
            print("    → **点数不足で当てはめ不能**（旗53 の 40日条件）\n"); continue
        star = (flx["r2"] >= 0.3 and np.isfinite(flx["acf1"]) and np.isfinite(flx["efold"])
                and flx["acf1"] >= ACF_THR and flx["efold"] <= EFOLD_MAX)
        why = []
        if flx["r2"] < 0.3:
            why.append(f"R²={flx['r2']:.2f}<0.3（駆動弱＝判定不能）")
        elif flx["acf1"] < ACF_THR:
            why.append(f"ACF1={flx['acf1']:.2f}<{ACF_THR}")
        elif flx["efold"] > EFOLD_MAX:
            why.append(f"e-fold={flx['efold']:.0f}>{EFOLD_MAX}日（季節メモリ）")
        print(f"    較正済み検出器：R²={flx['r2']:.2f} ACF1={flx['acf1']:+.2f} "
              f"e-fold={flx['efold']:.0f}日 → {'★短メモリ' if star else '／'.join(why) or '·'}")
        print(f"    → {'**同一地点解析に使える**' if star else '短メモリ判定は付かない（理由は上記）'}\n")

    print("  === 次にやること ===")
    print("  ★が付いた対については、**同じ場所で** タワー側（旗37 メモリ／旗38 SIF）と")
    print("  チャンバー側（旗40 メモリ／旗45 帰属）を並べる＝**弧が閉じる**。")
    print("  ★が付かない対でも、**タワー側とチャンバー側で結論が一致するか**は見る価値がある")
    print("  （例：どちらも『短メモリなし』なら、それはそれで整合の証拠）。")
    print("  留保：距離が近くても林分・処理・設置年が違えば同一地点とは言えない。")
    print("        COSORE の site 記述とタワーのメタデータを必ず突き合わせること。")


if __name__ == "__main__":
    main()
