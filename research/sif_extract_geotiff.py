"""SIF GeoTIFF 群から、各タワー座標の画素時系列を抜いて `<site>_sif.csv` を作る（旗38 の入力）。

取得元に依らない抽出器：GOSIF/CSIF を公式配布から直接DL、または GEE でエクスポートした GeoTIFF、
どれでも「日付つき GeoTIFF のフォルダ」さえあれば動く（GEE 認証も実在アセットIDも不要）。
rasterio でタワー(lat,lon)の画素値を読み、ファイル名から日付を解析して時系列にする。

    # 例: GOSIF 8-day (ファイル名 GOSIF_2018001.tif = 年+通日)
    python research/sif_extract_geotiff.py --coords site_coords.csv --tifdir /path/to/GOSIF \
        --date-regex "(\\d{4})(\\d{3})" --date-fmt yyyyddd --scale 0.0001 --nodata 32767

ファイル名の日付形式に合わせて --date-regex（キャプチャ順）と --date-fmt を指定する：
  yyyyddd  … 年(4桁)+通日(3桁)   例 GOSIF_2018001.tif
  yyyymmdd … 年月日(8桁)          例 CSIF_20180101.tif
  ymd      … 年,月,日 の3キャプチャ
--scale/--nodata は各プロダクトのスケール係数・欠測値（GOSIF: scale 0.0001, nodata 32767 等, 要確認）。
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def parse_date(name, rx, fmt):
    m = re.search(rx, name)
    if not m:
        return None
    g = m.groups()
    try:
        if fmt == "yyyyddd":
            return pd.Timestamp(year=int(g[0]), month=1, day=1) + pd.Timedelta(days=int(g[1]) - 1)
        if fmt == "yyyymmdd":
            s = g[0]
            return pd.Timestamp(year=int(s[:4]), month=int(s[4:6]), day=int(s[6:8]))
        if fmt == "ymd":
            return pd.Timestamp(year=int(g[0]), month=int(g[1]), day=int(g[2]))
    except (ValueError, IndexError):
        return None
    return None


def main():
    p = argparse.ArgumentParser(description="SIF GeoTIFF からタワー画素時系列を抽出")
    p.add_argument("--coords", default="site_coords.csv", help="site,lat,lon の csv")
    p.add_argument("--tifdir", required=True, help="SIF GeoTIFF のフォルダ")
    p.add_argument("--glob", default="*.tif*", help="GeoTIFF の glob（既定 *.tif*）")
    p.add_argument("--date-regex", default=r"(\d{4})(\d{3})")
    p.add_argument("--date-fmt", default="yyyyddd", choices=["yyyyddd", "yyyymmdd", "ymd"])
    p.add_argument("--scale", type=float, default=1.0, help="スケール係数（生値×scale）")
    p.add_argument("--nodata", type=float, default=None, help="欠測センチネル（これは NaN 化）")
    p.add_argument("--sites", nargs="+", help="対象サイト（既定=coords 全部で lat/lon 有）")
    a = p.parse_args()

    try:
        import rasterio
    except ImportError:
        print("rasterio が要ります： pip install rasterio"); return

    coords = pd.read_csv(a.coords)
    coords = coords.dropna(subset=["lat", "lon"])
    if a.sites:
        coords = coords[coords["site"].isin(a.sites)]
    if coords.empty:
        print("有効な座標が無い（site_coords.csv の lat/lon を確認）"); return

    tifs = sorted(Path(a.tifdir).glob(a.glob))
    if not tifs:
        print(f"GeoTIFF が無い： {a.tifdir}/{a.glob}"); return
    print(f"  GeoTIFF {len(tifs)} 枚 × サイト {len(coords)} 件を抽出")

    # 各サイトの (lat,lon) を保持、ファイルごとに全サイトの画素を読む（I/O 1 回/ファイル）
    series = {r.site: [] for r in coords.itertuples()}
    n_ok = 0
    for f in tifs:
        d = parse_date(f.name, a.date_regex, a.date_fmt)
        if d is None:
            continue
        try:
            with rasterio.open(f) as ds:
                for r in coords.itertuples():
                    try:
                        val = list(ds.sample([(float(r.lon), float(r.lat))]))[0][0]
                    except Exception:
                        val = None
                    v = float(val) if val is not None else float("nan")
                    if a.nodata is not None and v == a.nodata:
                        v = float("nan")
                    series[r.site].append((d, v * a.scale))
            n_ok += 1
        except Exception as e:
            print(f"  {f.name} 読めず: {type(e).__name__}")
    print(f"  日付解析できた GeoTIFF: {n_ok}/{len(tifs)}")

    for site, rows in series.items():
        df = pd.DataFrame(rows, columns=["date", "sif"]).dropna().sort_values("date")
        out = f"{site}_sif.csv"
        df.to_csv(out, index=False)
        print(f"  [出力] {out}（{len(df)} 点）")
    print("\n  → 次： python research/sif_respiration_step38.py --site JP-Tak --sif JP-Tak_sif.csv --qc-max 1")
    print("  留保：--scale/--nodata はプロダクト仕様を要確認。画素はタワー1点(最近傍)＝footprint不一致は旗38の留保通り。")


if __name__ == "__main__":
    main()
