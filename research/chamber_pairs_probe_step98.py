"""旗98：**同一地点に複数チャンバーがある組を数える**（下調べ・検定はしない）。

**手 D の実行可能性だけを確かめる**（`OPEN_QUESTIONS_OPTIONS.md`）。

## 何を測ろうとしているのか（**この道具ではまだ測らない**）

**~4 日メモリの正体**は、候補（基質の遅い転流・微生物動態・根の動態）まで絞れているが、
**どれも「チャンバー＋土壌センサの外側」で区別できない**と書いてきた。
**だが「どの空間スケールの現象か」は手元で測れる**——
**同じ日の残差が、チャンバー間でどれだけ揃うか**を見ればよい。

| 揃い方 | 意味 |
|---|---|
| **同一地点のチャンバー間で強く揃う** | 駆動は**地点スケール**（気象・土壌水文）＝**未観測だが局所生物ではない** |
| **同一地点でも揃わない** | 駆動は**チャンバー・スケール**（微生物群集・根の分布などの微小生息場所） |
| **地点間でも揃う** | **広域の気象**——**我々が測っていない気象量**（大きな示唆） |

**機構を特定せずに、空間スケールで絞れる。**

## **この道具は、実行可能性だけを出す。相関は計算しない。**

**相関は検定の答えそのもの**であり、**事前登録の前に見てはいけない**（旗94 と同じ作法）。
出すのは**組の数・重なる日数・年数・土壌温度の有無**だけである。

  1. **座標で地点をまとめる**（`--km` 以内を同一地点とする。既定 1 km）
  2. **同一地点に 2 本以上あるチャンバー**を列挙し、**対ごとの重なる日数・年数**を出す
  3. **下限（重なり 60 日・3 年）を満たす対が何組あるか**
  4. **地点間の対**（遠く離れた組）についても、**重なる日数**だけ数える

    python research/chamber_pairs_probe_step98.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore
from colocate_step51 import haversine

MIN_DAYS, MIN_YEARS = 60, 3        # 本研究の下限（旗58 以来）


def load_daily(path):
    """チャンバーの日次（Rs・Tsoil・SM）。**残差も相関も作らない。**"""
    try:
        df, st, sm = load_cosore(path, None)
    except Exception as e:
        return None, f"{type(e).__name__}: {str(e)[:60]}", None, None
    if df is None or "Rs" not in df:
        return None, "Rs が無い", st, sm
    cols = ["Rs"] + [c for c in ("Tsoil", "SM") if c in df]
    d = df[cols].groupby(df.index.normalize()).mean()
    return d, None, st, sm


def main():
    ap = argparse.ArgumentParser(description="同一地点の複数チャンバーを数える")
    ap.add_argument("--cosore-dir", required=True)
    ap.add_argument("--km", type=float, default=1.0)
    a = ap.parse_args()

    print("=== 旗98：同一地点に複数チャンバーがある組を数える（下調べ・検定はしない）===")
    print("  **相関は計算しない**——**それは検定の答えであり、事前登録の前に見てはいけない**。")
    print(f"  **{a.km:.1f} km 以内を同一地点**とする。下限は**重なり {MIN_DAYS} 日・{MIN_YEARS} 年**。\n")

    root = Path(a.cosore_dir)
    desc = pd.read_csv(root / "description.csv")
    recs = []
    for _, r in desc.iterrows():
        ds = str(r["CSR_DATASET"])
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            la, lo = float(r["CSR_LATITUDE"]), float(r["CSR_LONGITUDE"])
        except (TypeError, ValueError, KeyError):
            continue
        if not (np.isfinite(la) and np.isfinite(lo)):
            continue
        d, err, st, sm = load_daily(f)
        recs.append({"ds": ds, "lat": la, "lon": lo, "igbp": str(r.get("CSR_IGBP", "")),
                     "d": d, "err": err, "st": st, "sm": sm})
    print(f"  読めたデータセット：{sum(x['d'] is not None for x in recs)}"
          f" / 座標のある {len(recs)}\n")

    # ── 座標で地点をまとめる（単連結）──
    n = len(recs)
    parent = list(range(n))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for i, j in itertools.combinations(range(n), 2):
        if haversine(recs[i]["lat"], recs[i]["lon"], recs[j]["lat"], recs[j]["lon"]) <= a.km:
            parent[find(i)] = find(j)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    multi = {k: v for k, v in groups.items() if len(v) >= 2}
    print(f"  **地点の数（{a.km:.1f} km 単連結）：{len(groups)}**"
          f"／うち**チャンバーが 2 本以上：{len(multi)} 地点**\n")

    ok_pairs = same_pairs = 0
    for k, idxs in sorted(multi.items(), key=lambda kv: -len(kv[1])):
        members = [recs[i] for i in idxs]
        print(f"  ━━ 地点（{len(members)} 本）"
              f"{members[0]['lat']:.4f}, {members[0]['lon']:.4f}"
              f"／{members[0]['igbp']} ━━")
        for m in members:
            if m["d"] is None:
                print(f"    {m['ds']:<32}**読めない**（{m['err']}）"); continue
            has = "T=" + (m["st"] or "無") + " SM=" + (m["sm"] or "無")
            print(f"    {m['ds']:<32}{len(m['d']):>6} 日／"
                  f"{m['d'].index.min():%Y-%m}〜{m['d'].index.max():%Y-%m}／{has}")
        usable = [m for m in members if m["d"] is not None and "Tsoil" in m["d"]]
        if len(usable) < 2:
            print(f"    → **土壌温度のあるチャンバーが 2 本未満**＝この地点は使えない\n")
            continue
        print(f"    ── 対ごとの重なり ──")
        for x, y in itertools.combinations(usable, 2):
            common = x["d"].index.intersection(y["d"].index)
            yrs = pd.Index(common).year.nunique() if len(common) else 0
            same_pairs += 1
            ok = len(common) >= MIN_DAYS and yrs >= MIN_YEARS
            ok_pairs += int(ok)
            print(f"      {x['ds'][:26]:<28}× {y['ds'][:26]:<28}"
                  f"重なり {len(common):>5} 日／{yrs:>2} 年  "
                  f"{'**使える**' if ok else '下限未満'}")
        print()

    # ── 地点間（遠い組）の重なりも数える ──
    usable_all = [r for r in recs if r["d"] is not None and "Tsoil" in r["d"]]
    far_ok = far_all = 0
    for x, y in itertools.combinations(usable_all, 2):
        if haversine(x["lat"], x["lon"], y["lat"], y["lon"]) <= a.km:
            continue
        common = x["d"].index.intersection(y["d"].index)
        if len(common) == 0:
            continue
        far_all += 1
        if len(common) >= MIN_DAYS and pd.Index(common).year.nunique() >= MIN_YEARS:
            far_ok += 1

    print("  === まとめ ===")
    print(f"  **同一地点の対**：{same_pairs} 組（うち**下限を満たす {ok_pairs} 組**）")
    print(f"  **地点間の対**：期間が重なる {far_all} 組（うち**下限を満たす {far_ok} 組**）")
    print("\n  === 次の判断（**事前登録の前に決める**）===")
    print("  ・**同一地点の対が 5 組以上**あり、**地点間の対も数組ある** → **旗99 を事前登録して検定する**")
    print("  ・**同一地点の対が 5 組未満** → **空間スケールの比較はできない**と記して打ち切る")
    print("  ・**地点間の対が 0** → **「広域の気象か」の枝は検定できない**と明記して、")
    print("    **同一地点の中だけ**（地点スケール vs チャンバー・スケール）で組み直す")
    print("  留保：")
    print("   ・**相関はまだ一度も計算していない**（事前登録の前に答えを見ない）。")
    print("   ・**同じ地点でも処理・林分が違えば別の生態系**（旗51 と同じ注意）。")
    print("     **座標が近いことは同一条件の必要条件であって十分条件ではない。**")
    print("   ・**残差の同期を見るには、各チャンバーで駆動を十分に引く必要がある**——")
    print("     **共通の気象が残差に残っていれば自明に揃う**。**そこは旗99 の設計で受ける。**")


if __name__ == "__main__":
    main()
