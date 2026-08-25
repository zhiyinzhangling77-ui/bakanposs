"""旗64：タワー座標の健全性を点検する——旗51/63 の前提が壊れていた疑いを確かめる。

旗63 の出力で **UEYAMA_YAMASHIRO（34.795°N, 135.846°E＝京都付近）の最近傍が「JP-Tak 3095km」**
と出た。JP-Tak（高山）から京都は**約200km**であり、**3095km はありえない**。
＝旗51/63 が使っている**タワー座標が壊れている**疑いが濃い。

もしそうなら、**旗51 の結論「手元のタワーと COSORE は10km以内に一組も無い」は無効**であり、
「穴②は外部依存で確定」という判断もやり直しになる。

本ツールは BADM から抽出した座標をそのまま並べ、次を点検する：
  1. 緯度が [-90,90]・経度が [-180,180] に収まっているか
  2. **緯度と経度が入れ替わっていないか**（日本のサイトなら lat≈24–46, lon≈123–146）
  3. 既知の距離と矛盾しないか（サイト間距離を出して、明らかに変な値が無いか）

    python research/check_tower_coords_step64.py --sites JP-Tak JP-Fhk JP-Tmd JP-Fjy JP-Tef \
        JP-Ta2 JP-Mse JP-BBY JP-SMF JP-Spp
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from colocate_step51 import haversine, tower_coords

# 日本のサイトが収まるべき範囲（点検用の緩い枠）
JP_LAT = (24.0, 46.0)
JP_LON = (123.0, 146.0)


def main():
    p = argparse.ArgumentParser(description="タワー座標の健全性点検")
    p.add_argument("--sites", nargs="+", required=True)
    p.add_argument("--data-dir", default="/mnt/hdd/JAPANFLUX")
    p.add_argument("--coords", default=None)
    a = p.parse_args()

    tw = tower_coords(a.sites, a.data_dir, a.coords)
    print("=== 旗64：タワー座標の健全性点検 ===")
    print(f"  取得できたサイト：{len(tw)}／{len(a.sites)}\n")
    print(f"  {'site':<10}{'lat':>10}{'lon':>11}  判定")
    bad = []
    for s, (lat, lon) in sorted(tw.items()):
        notes = []
        if not (-90 <= lat <= 90):
            notes.append("**緯度が範囲外**")
        if not (-180 <= lon <= 180):
            notes.append("**経度が範囲外**")
        # 日本のサイト（JP-）なら日本の枠に入るはず。入らず、かつ入れ替えると入るなら swap を疑う
        if s.upper().startswith("JP-"):
            in_jp = JP_LAT[0] <= lat <= JP_LAT[1] and JP_LON[0] <= lon <= JP_LON[1]
            swapped_ok = JP_LAT[0] <= lon <= JP_LAT[1] and JP_LON[0] <= lat <= JP_LON[1]
            if not in_jp:
                notes.append("**日本の範囲外**" + ("＝**緯度経度の入れ替わりを疑う**" if swapped_ok else ""))
        v = "／".join(notes) if notes else "ok"
        if notes:
            bad.append(s)
        print(f"  {s:<10}{lat:>10.4f}{lon:>11.4f}  {v}")

    print(f"\n  === サイト間距離（km）——明らかに変な値が無いか ===")
    names = sorted(tw)
    print("  " + " " * 10 + "".join(f"{n:>10}" for n in names))
    for i in names:
        row = "".join(f"{haversine(*tw[i], *tw[j]):>10.0f}" for j in names)
        print(f"  {i:<10}{row}")

    print("\n  === 判定 ===")
    if bad:
        print(f"  **{len(bad)} サイトの座標がおかしい**：{bad}")
        print("  → 旗51（10km以内に一組も無い）と 旗63（最近傍距離）は**この座標に依存している**ため、")
        print("     **どちらの結論もやり直しが要る**。抽出元（research/sif_coords.py の BADM 解析）を直すこと。")
    else:
        print("  座標に明らかな異常は無い。旗51/63 の距離計算はこの点では信用してよい。")
        print("  ※ただし『範囲内』は正しさを保証しない。可能なら公開サイト情報と突き合わせること。")


if __name__ == "__main__":
    main()
