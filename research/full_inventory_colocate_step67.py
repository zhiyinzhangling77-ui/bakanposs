"""旗67：手元の**全サイト**を棚卸しして、COSORE との同一地点をもう一度探す。

旗51/63/66 は**登録済みの10サイト**としか照合していない。だが `/mnt/hdd/JAPANFLUX` には
**31サイト**あるはずで、**未登録のサイトが COSORE と同一地点かもしれない**。
実際、旗63 の★リストにある `d20190830_LIANG`（31.85°N, 131.30°E）は**宮崎＝日本**である。

穴② は JP-Fhk 1組で閉じたが**1勝1敗**（JP-Tef は食い違い）。**3組目以降が見つかれば強くなる**。
外部データを取りに行く前に、**手元にあるものを数え直す**のが順序。

本ツールは：
  1. データディレクトリを走査して**サイトらしきものを全部**見つける（BADM のフォルダ名・ファイル名から）
  2. 各サイトの座標を**修正済みの抽出器**（旗64）で引く
  3. COSORE の全チャンバーと突き合わせ、**近い組を距離順に**出す（10km で切らない）

    python research/full_inventory_colocate_step67.py --cosore-dir /mnt/hdd/cosore-0.7.0 \
        --data-dir /mnt/hdd/JAPANFLUX
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from sif_coords import coords_from_badm
from colocate_step51 import haversine

# FLUXNET 流のサイトコード（例 JP-Tak, CN-HaM, MN-Hst, US-Wkg）。
# **区切りは `_` や `.` でもよい**——`\b` だと `_JP-Tak_` に効かない（旗67 第1版の誤り）。
# 誤検出を減らすため、ハイフンの後に**英字を1つ以上**含むことを要求する（"HH-20" 等を弾く）。
CODE_RE = re.compile(r"(?<![A-Za-z0-9])([A-Z]{2}-(?=[A-Za-z0-9]{2,6}(?![A-Za-z0-9]))"
                     r"[A-Za-z0-9]*[A-Za-z][A-Za-z0-9]*)(?![A-Za-z0-9])")


def discover_codes(data_dir):
    """データディレクトリ配下のパス名から、サイトコードらしきものを全部拾う。"""
    root = Path(data_dir)
    codes = {}
    for p in root.rglob("*"):
        if "__MACOSX" in p.parts:
            continue
        for m in CODE_RE.finditer(p.name):
            codes.setdefault(m.group(1), 0)
            codes[m.group(1)] += 1
        # ディレクトリ名も見る（BADM_JP-Tak_... のような形）
        for part in p.parts[len(root.parts):]:
            for m in CODE_RE.finditer(part):
                codes.setdefault(m.group(1), 0)
                codes[m.group(1)] += 1
    return sorted(codes)


def main():
    p = argparse.ArgumentParser(description="全サイト棚卸し＋同一地点探索")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--data-dir", default="/mnt/hdd/JAPANFLUX")
    p.add_argument("--km", type=float, default=30.0, help="この距離まで出す（切らずに近い順）")
    p.add_argument("--top", type=int, default=30)
    a = p.parse_args()

    print("=== 旗67：手元の全サイトを棚卸しして同一地点を探し直す ===")
    codes = discover_codes(a.data_dir)
    print(f"  データ配下で見つかったサイトコード：{len(codes)} 件")
    print(f"  {codes}\n")

    tw = {}
    failed = []
    for c in codes:
        try:
            lat, lon = coords_from_badm(a.data_dir, c)
        except Exception:
            lat = lon = None
        if lat is None or lon is None or not (np.isfinite(lat) and np.isfinite(lon)):
            failed.append(c); continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            failed.append(c + "(範囲外)"); continue
        tw[c] = (float(lat), float(lon))
    print(f"  座標を引けたサイト：{len(tw)} 件")
    if failed:
        print(f"  引けなかった：{failed}")
    if not tw:
        print("  → 座標が一つも引けない。--data-dir を確認すること。"); return

    desc = pd.read_csv(Path(a.cosore_dir) / "description.csv")
    pairs = []
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        f = Path(a.cosore_dir) / "datasets" / f"data_{ds}.csv"
        try:
            clat, clon = float(d["CSR_LATITUDE"]), float(d["CSR_LONGITUDE"])
        except (TypeError, ValueError, KeyError):
            continue
        if not (np.isfinite(clat) and np.isfinite(clon)):
            continue
        for s, (tlat, tlon) in tw.items():
            km = haversine(tlat, tlon, clat, clon)
            if km <= a.km:
                pairs.append({"km": km, "site": s, "ds": ds, "igbp": ig,
                              "has_data": f.exists()})
    if not pairs:
        print(f"\n  **{a.km:.0f}km 以内の組は無い**。")
        return
    df = pd.DataFrame(pairs).sort_values("km")
    print(f"\n  **{a.km:.0f}km 以内の組：{len(df)} 件**（近い順）")
    print(f"  {'距離km':>8}  {'タワー':<10}{'COSOREチャンバー':<34}{'IGBP':<15}データ")
    for _, r in df.head(a.top).iterrows():
        print(f"  {r['km']:>8.2f}  {r['site']:<10}{r['ds']:<34}{r['igbp'][:13]:<15}"
              f"{'あり' if r['has_data'] else '**無し**'}")

    usable = df[(df["km"] <= 1.0) & df["has_data"]]
    print(f"\n  === まとめ ===")
    print(f"  1km 以内かつデータあり＝**同一地点解析に使える候補**：{len(usable)} 件")
    if len(usable):
        for _, r in usable.iterrows():
            print(f"    {r['site']} ↔ {r['ds']}（{r['km']:.2f} km）")
    print("\n  → 既知の組（JP-Fhk↔UEYAMA_HOKUROKU, JP-Tef↔UEYAMA_TESHIO）以外が出れば、")
    print("     旗66 の同一地点比較を**3組目以降に広げられる**＝穴②の証拠が強くなる。")
    print("  留保：コードの自動検出はパス名頼りで、取りこぼし・誤検出がありうる。")
    print("        距離が近くても林分・処理・設置年が違えば同一地点とは言えない（旗51 と同じ注意）。")


if __name__ == "__main__":
    main()
