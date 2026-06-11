"""
load_hsaf.py
============
H-SAF SM products (H26 / H141 / H142) から
Oran・Tarazona のピクセル値を抽出して CSV に出力する。

出力 CSV 列:
  date, site, source,
  h26_swi_L1, h26_swi_L2, h26_swi_L3, h26_swi_L4,   ← H26 4 layer
  h141_swi_L1, h141_swi_L2, h141_swi_L3, h141_swi_L4, ← H141 (〜2021)
  h142_swi_L1, h142_swi_L2, h142_swi_L3, h142_swi_L4  ← H142 (〜2021)

Usage
-----
python load_hsaf.py \
    --h26-dir   /mnt/hdd/Dataset/HSAF/H26 \
    --h141-dir  /mnt/hdd/Dataset/HSAF/H141 \
    --h142-dir  /mnt/hdd/Dataset/HSAF/H142 \
    --output    data/hsaf_daily.csv

サイト座標 (必要なら --sites-json で上書き):
  Tarazona de la Mancha (TzM) : 39.2660°N, -1.9397°W, 727 m.a.s.l.
  Oran                        : 38.8200°N, -1.8600°W, 720 m.a.s.l.

期間ルール (研究設計):
  Oran  2018        → H141 のみ
  Tarazona 2020–2021 → H142 のみ
  Tarazona 2022–2024 → H26 NRT (H141/H142 はこの期間をカバーしない)
"""

import argparse
import os
import json
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import netCDF4 as nc

# -----------------------------------------------------------------------
# サイト定義
# -----------------------------------------------------------------------
DEFAULT_SITES = {
    # Tarazona de la Mancha (TzM): 39.2660°N, -1.9397°W, 727 m.a.s.l.
    "Tarazona": {"lat": 39.2660, "lon": -1.9397},
    # Oran: 38.82°N, -1.86°W, 720 m.a.s.l.
    "Oran":     {"lat": 38.8200, "lon": -1.8600},
}

# 研究期間ごとの SM ソースルール
# H141 カバー: 〜2018-12-31
# H142 カバー: 2019-01-01 〜 2021-12-31
# H26  カバー: 2019-〜 現在 (NRT)
H141_END  = pd.Timestamp("2018-12-31")
H142_START = pd.Timestamp("2019-01-01")
H142_END  = pd.Timestamp("2021-12-31")
H26_START = pd.Timestamp("2019-01-01")   # NRT archive start (実際は〜2019)


def source_for_date(date: pd.Timestamp) -> list:
    """その日付に利用できる SM source リストを返す"""
    sources = []
    if date <= H141_END:
        sources.append("H141")
    if H142_START <= date <= H142_END:
        sources.append("H142")
    if date >= H26_START:
        sources.append("H26")
    return sources


# -----------------------------------------------------------------------
# 座標→最近傍インデックス
# -----------------------------------------------------------------------
def nearest_idx(lat_arr: np.ndarray, lon_arr: np.ndarray,
                target_lat: float, target_lon: float):
    """2D または 1D lat/lon 配列から最近傍インデックスを返す"""
    if lat_arr.ndim == 2:
        dist = (lat_arr - target_lat) ** 2 + (lon_arr - target_lon) ** 2
        idx = np.unravel_index(np.argmin(dist), dist.shape)
        return idx
    else:
        # 1D: H141/H142 は irregular grid の可能性あり
        dist = (lat_arr - target_lat) ** 2 + (lon_arr - target_lon) ** 2
        return (np.argmin(dist),)


# -----------------------------------------------------------------------
# H26 NetCDF loader
# -----------------------------------------------------------------------
def load_h26_file(nc_path: str, sites: dict) -> dict:
    """
    H26 .nc から各サイトの SWI 4 layer を抽出。
    変数名: SWI_002, SWI_005, SWI_010, SWI_020  (or similar T-value names)
    Returns: {site_name: {layer: value, ...}}
    """
    result = {}
    try:
        ds = nc.Dataset(nc_path, "r")
        vars_all = list(ds.variables.keys())

        # lat/lon 変数を探す (名前が製品バージョンで異なる場合あり)
        lat_candidates = [v for v in vars_all if "lat" in v.lower()]
        lon_candidates = [v for v in vars_all if "lon" in v.lower()]
        if not lat_candidates or not lon_candidates:
            ds.close()
            return result

        lat_arr = ds[lat_candidates[0]][:].data
        lon_arr = ds[lon_candidates[0]][:].data

        # SWI 変数を探す
        swi_vars = sorted([v for v in vars_all if v.upper().startswith("SWI")])
        if not swi_vars:
            # フォールバック: layer 変数を探す
            swi_vars = sorted([v for v in vars_all
                               if any(k in v.lower() for k in ["swi", "sm", "wetness"])])

        for site_name, coords in sites.items():
            idx = nearest_idx(lat_arr, lon_arr, coords["lat"], coords["lon"])
            # 最初のファイルだけグリッドセル中心座標をログ (ズレ確認用)
            grid_lat = float(lat_arr[idx[0]]) if lat_arr.ndim == 1 else float(lat_arr[idx[0], idx[1]])
            grid_lon = float(lon_arr[idx[0]]) if lon_arr.ndim == 1 else float(lon_arr[idx[0], idx[1]])
            dist_km  = ((grid_lat - coords["lat"])**2 + (grid_lon - coords["lon"])**2)**0.5 * 111
            key = f"_h26_logged_{site_name}"
            if not getattr(load_h26_file, key, False):
                print(f"  [Grid/H26] {site_name}: target=({coords['lat']:.4f},{coords['lon']:.4f})"
                      f" → grid=({grid_lat:.4f},{grid_lon:.4f})  dist≈{dist_km:.1f} km")
                setattr(load_h26_file, key, True)
            row = {}
            for i, vname in enumerate(swi_vars[:4]):   # 最大 4 layer
                try:
                    val = float(ds[vname][0][idx] if ds[vname].ndim == 3
                                else ds[vname][idx])
                    # スケール・フィル値処理
                    fill = getattr(ds[vname], "_FillValue", None)
                    scale = getattr(ds[vname], "scale_factor", 1.0)
                    offset = getattr(ds[vname], "add_offset", 0.0)
                    if fill is not None and val == fill:
                        val = np.nan
                    else:
                        val = val * scale + offset
                    row[f"h26_swi_L{i+1}"] = round(float(val), 4)
                except Exception:
                    row[f"h26_swi_L{i+1}"] = np.nan
            result[site_name] = row

        ds.close()
    except Exception as e:
        print(f"  [WARN] H26 load error {nc_path}: {e}")
    return result


# -----------------------------------------------------------------------
# H141 / H142 NetCDF loader
# -----------------------------------------------------------------------
def load_h141_h142_file(nc_path: str, product: str, sites: dict) -> dict:
    """
    H141 / H142 .nc から各サイトの SWI 4 layer を抽出。
    変数名: SWI_L1, SWI_L2, SWI_L3, SWI_L4 など (PUM 参照)
    """
    prefix = product.lower()  # "h141" or "h142"
    result = {}
    try:
        ds = nc.Dataset(nc_path, "r")
        vars_all = list(ds.variables.keys())

        lat_candidates = [v for v in vars_all if "lat" in v.lower()]
        lon_candidates = [v for v in vars_all if "lon" in v.lower()]
        if not lat_candidates or not lon_candidates:
            ds.close()
            return result

        lat_arr = ds[lat_candidates[0]][:].data
        lon_arr = ds[lon_candidates[0]][:].data

        # SWI / layer 変数
        layer_vars = sorted([v for v in vars_all if v.upper().startswith("SWI")])
        if not layer_vars:
            layer_vars = sorted([v for v in vars_all
                                 if any(k in v.lower() for k in ["swi", "layer", "wetness"])])

        for site_name, coords in sites.items():
            idx = nearest_idx(lat_arr, lon_arr, coords["lat"], coords["lon"])
            row = {}
            for i, vname in enumerate(layer_vars[:4]):
                try:
                    if ds[vname].ndim == 3:
                        val = float(ds[vname][0][idx])
                    elif ds[vname].ndim == 2:
                        val = float(ds[vname][idx])
                    else:
                        val = float(ds[vname][idx[0]])
                    fill = getattr(ds[vname], "_FillValue", None)
                    scale = getattr(ds[vname], "scale_factor", 1.0)
                    offset = getattr(ds[vname], "add_offset", 0.0)
                    if fill is not None and val == fill:
                        val = np.nan
                    else:
                        val = val * scale + offset
                    row[f"{prefix}_swi_L{i+1}"] = round(float(val), 4)
                except Exception:
                    row[f"{prefix}_swi_L{i+1}"] = np.nan
            result[site_name] = row

        ds.close()
    except Exception as e:
        print(f"  [WARN] {product} load error {nc_path}: {e}")
    return result


# -----------------------------------------------------------------------
# ファイル日付解析
# -----------------------------------------------------------------------
def date_from_filename(nc_path: str) -> pd.Timestamp | None:
    """ファイル名から YYYYMMDD を抽出"""
    stem = Path(nc_path).stem  # 例: h26_2022030100_R01
    for part in stem.split("_"):
        if len(part) >= 8 and part[:8].isdigit():
            try:
                return pd.Timestamp(part[:8])
            except Exception:
                continue
    return None


# -----------------------------------------------------------------------
# メイン処理
# -----------------------------------------------------------------------
def build_index(data_dir: str, pattern: str) -> dict:
    """date → filepath の辞書を作る"""
    idx = {}
    if not data_dir or not os.path.isdir(data_dir):
        return idx
    for fp in glob.glob(os.path.join(data_dir, pattern)):
        date = date_from_filename(fp)
        if date is not None:
            idx[date] = fp
    print(f"  Indexed {len(idx)} files in {data_dir}")
    return idx


def main():
    parser = argparse.ArgumentParser(description="H-SAF SM pixel extractor")
    parser.add_argument("--h26-dir",  default=None)
    parser.add_argument("--h141-dir", default=None)
    parser.add_argument("--h142-dir", default=None)
    parser.add_argument("--output",   default="data/hsaf_daily.csv")
    parser.add_argument("--sites-json", default=None,
                        help="JSON file: {site: {lat, lon}} でデフォルトを上書き")
    parser.add_argument("--start", default="2018-01-01")
    parser.add_argument("--end",   default="2024-12-31")
    args = parser.parse_args()

    # サイト定義
    sites = DEFAULT_SITES.copy()
    if args.sites_json:
        with open(args.sites_json) as f:
            sites.update(json.load(f))
    print(f"Sites: {list(sites.keys())}")

    # ファイルインデックス構築
    print("\nBuilding file indices...")
    h26_idx  = build_index(args.h26_dir,  "h26_*.nc")
    h141_idx = build_index(args.h141_dir, "h141_*.nc")
    h142_idx = build_index(args.h142_dir, "h142_*.nc")

    # 出力ディレクトリ
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    # 全日付ループ
    datelist = pd.date_range(args.start, args.end, freq="D")
    records = []

    print(f"\nExtracting pixels: {args.start} → {args.end}")
    for date in datelist:
        sources = source_for_date(date)

        for site_name in sites:
            row = {"date": date.strftime("%Y-%m-%d"), "site": site_name,
                   "available_sources": ",".join(sources)}

            # H26
            h26_data = {}
            if "H26" in sources and date in h26_idx:
                h26_data = load_h26_file(h26_idx[date], sites).get(site_name, {})
            for k in ["h26_swi_L1","h26_swi_L2","h26_swi_L3","h26_swi_L4"]:
                row[k] = h26_data.get(k, np.nan)

            # H141
            h141_data = {}
            if "H141" in sources and date in h141_idx:
                h141_data = load_h141_h142_file(h141_idx[date], "H141", sites).get(site_name, {})
            for k in ["h141_swi_L1","h141_swi_L2","h141_swi_L3","h141_swi_L4"]:
                row[k] = h141_data.get(k, np.nan)

            # H142
            h142_data = {}
            if "H142" in sources and date in h142_idx:
                h142_data = load_h141_h142_file(h142_idx[date], "H142", sites).get(site_name, {})
            for k in ["h142_swi_L1","h142_swi_L2","h142_swi_L3","h142_swi_L4"]:
                row[k] = h142_data.get(k, np.nan)

            records.append(row)

        if date.day == 1:
            print(f"  {date.strftime('%Y-%m')}  ({len(records)} rows so far)")

    df = pd.DataFrame(records)
    df.to_csv(args.output, index=False)
    print(f"\n[Done] {len(df)} rows → {args.output}")
    print(df.head(10).to_string())


if __name__ == "__main__":
    main()
