"""旗51（穴②）：外部連絡の前に確かめる——手元のチャンバーと手元のタワーは同一地点ではないか。

穴②：メモリの弧は フラックス(旗25/37)＋衛星SIF(旗38) が JP-Tak、チャンバー確認(旗40/45) が
別サイト群＝**同一地点で3観測系を突き合わせていない**。従来の対処案は「TKY-Takayama のチャンバーを
Kishimoto-Mo 氏へ依頼する」だったが、その前に確かめるべきことがある：

  **COSORE には日本のチャンバーサイトが入っている**（UEYAMA_TESHIO / YAMASHIRO / HOKUROKU 等）。
  これらが我々の JapanFlux タワーと同一地点なら、**外部連絡なしで弧が閉じる**。

本ツールは COSORE の description.csv の座標と、BADM から抽出したタワー座標を突き合わせ、
**距離が近い組（既定 10km 以内）を列挙**する。同一地点の組が見つかれば、その site で
  ・タワー側：フラックス由来の GER 残差メモリ（旗37）と SIF テスト（旗38）
  ・チャンバー側：分割を通さない Rs 残差メモリ（旗40）と候補帰属（旗45）
を**同じ場所で**並べられる＝弧が完全に閉じる。

    python research/colocate_step51.py --cosore-dir /mnt/hdd/cosore-0.7.0 \
        --sites JP-Tak JP-Fhk JP-Tmd JP-Fjy JP-Tef JP-Ta2 JP-Mse JP-BBY JP-SMF JP-Spp
（タワー座標は research/sif_coords.py と同じ BADM 抽出を使う。--coords で既存CSVを渡してもよい）
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))


def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp, dl = p2 - p1, np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(a))


def tower_coords(sites, data_dir, coords_csv=None):
    import pandas as pd
    if coords_csv and Path(coords_csv).exists():
        df = pd.read_csv(coords_csv)
        col = {c.lower(): c for c in df.columns}
        return {str(r[col.get("site", "site")]): (float(r[col.get("lat", "lat")]),
                                                  float(r[col.get("lon", "lon")]))
                for _, r in df.iterrows()}
    from sif_coords import coords_from_badm                 # BADM から抽出（旗38 と同じ）
    out = {}
    for s in sites:
        try:
            lat, lon = coords_from_badm(data_dir, s)
            if lat is not None and lon is not None:
                out[s] = (float(lat), float(lon))
        except Exception as e:
            print(f"  {s}: 座標取得できず（{type(e).__name__}）")
    return out


def main():
    p = argparse.ArgumentParser(description="COSOREチャンバーとフラックスタワーの同一地点探索")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--sites", nargs="+", required=True)
    p.add_argument("--data-dir", default="/mnt/hdd/JAPANFLUX")
    p.add_argument("--coords", default=None, help="タワー座標CSV（site,lat,lon）があれば")
    p.add_argument("--km", type=float, default=10.0)
    a = p.parse_args()

    import pandas as pd
    desc = pd.read_csv(Path(a.cosore_dir) / "description.csv")
    need = {"CSR_LATITUDE", "CSR_LONGITUDE"}
    if not need <= set(desc.columns):
        print("  description.csv に座標列が無い"); return

    tw = tower_coords(a.sites, a.data_dir, a.coords)
    print(f"=== 旗51：チャンバー×タワーの同一地点探索（{a.km:.0f}km 以内）===")
    print(f"  タワー座標を取得できたサイト：{len(tw)}／{len(a.sites)}")
    if not tw:
        print("  タワー座標が無いので判定できない。research/sif_coords.py の出力CSVを --coords で渡すこと。")
        return

    hits = []
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"])
        try:
            clat = float(d["CSR_LATITUDE"]); clon = float(d["CSR_LONGITUDE"])
        except (TypeError, ValueError):
            continue
        if not (np.isfinite(clat) and np.isfinite(clon)):
            continue
        for s, (tlat, tlon) in tw.items():
            km = haversine(tlat, tlon, clat, clon)
            if km <= a.km:
                hits.append((km, s, ds, str(d.get("CSR_IGBP", "")), clat, clon))

    if not hits:
        print(f"\n  **{a.km:.0f}km 以内の組み合わせなし**＝手元のデータでは同一地点にならない。")
        print("  ＝穴②は外部依存のまま（TKY-Takayama のチャンバーは要連絡）。")
        print("  参考：--km を広げると近傍の組が見えるが、同一地点でなければ弧は閉じない。")
        return

    print(f"\n  **同一地点候補 {len(hits)} 組**")
    print(f"  {'距離km':>7}  {'タワー':<10}{'COSOREチャンバー':<34}{'IGBP':<14}")
    for km, s, ds, ig, clat, clon in sorted(hits):
        print(f"  {km:>7.2f}  {s:<10}{ds:<34}{ig[:12]:<14}")
    print("\n  → この組があれば、**同じ場所で** タワー側(旗37メモリ/旗38 SIF)と")
    print("     チャンバー側(旗40メモリ/旗45帰属)を並べられる＝**穴②が外部連絡なしで閉じる**。")
    print("  次にやること：その site を旗37/38/40/45 に通し、同一地点で結論が一致するかを見る。")
    print("  留保：距離が近くても林分・処理が違えば同一地点とは言えない。COSORE の site 記述と")
    print("        タワーのメタデータで**土地被覆・処理・設置年**を必ず突き合わせること。")


if __name__ == "__main__":
    main()
