"""TROPOMI 真SIF(Caltech/TROPOSIF 等の NetCDF)から、タワー画素の時系列を抜いて `<site>_sif.csv` を作る。

旗38 で GOSIF(再構成8日) は記憶を説明せず。旗39 で記憶は「生物×分割窓の混在」と判明。
**真のSIF(TROPOMI, 実測・ほぼ日次)** は GOSIF の2弱点(再構成・8日)を外した、より公平なテスト
＝「生物成分がどれだけ残るか」を測る。TROPOMI は NetCDF 配布なので、この専用抽出器で読む。

配布形態が多様(変数名・座標名・単一多時刻ファイル/時刻別ファイル)なので、**変数と座標を自動検出**し、
見つからなければ**中身を一覧表示**して終わる(捏造せず、実物に合わせて --var 等を指定してもらう)。

    # まず中身を見る（変数・座標名を確認）
    python research/sif_extract_netcdf.py --ncdir /path/to/tropomi --list
    # 抽出（変数名が分かれば --var で指定, 既定は名前に sif を含む変数を自動選択）
    python research/sif_extract_netcdf.py --coords site_coords.csv --ncdir /path/to/tropomi \
        --var SIF_743 --lat-name lat --lon-name lon

入手先(ローカルで): Caltech TROPOMI SIF  ftp://fluo.gps.caltech.edu/data/tropomi/ (2018-03〜2021-07, ほぼ日次)。
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

LAT_NAMES = ["lat", "latitude", "y", "Latitude", "LAT"]
LON_NAMES = ["lon", "longitude", "x", "Longitude", "LON"]
TIME_NAMES = ["time", "Time", "date", "t"]


def _pick(cands, names):
    for n in names:
        if n in cands:
            return n
    return None


def main():
    p = argparse.ArgumentParser(description="TROPOMI NetCDF SIF からタワー画素時系列を抽出")
    p.add_argument("--coords", default="site_coords.csv")
    p.add_argument("--ncdir", required=True, help="NetCDF フォルダ")
    p.add_argument("--glob", default="*.nc")
    p.add_argument("--var", default=None, help="SIF 変数名（既定=名前に sif を含む変数を自動選択）")
    p.add_argument("--lat-name", default=None)
    p.add_argument("--lon-name", default=None)
    p.add_argument("--time-name", default=None)
    p.add_argument("--date-regex", default=None, help="時刻座標が無いとき、ファイル名から日付を解析")
    p.add_argument("--date-fmt", default="yyyymmdd", choices=["yyyymmdd", "yyyyddd", "ymd"])
    p.add_argument("--sites", nargs="+")
    p.add_argument("--list", action="store_true", help="中身(変数・座標)を表示して終了")
    a = p.parse_args()

    try:
        import xarray as xr
    except ImportError:
        print("xarray が要ります： pip install xarray netcdf4"); return

    files = sorted(Path(a.ncdir).glob(a.glob))
    if not files:
        print(f"NetCDF が無い： {a.ncdir}/{a.glob}"); return

    ds0 = xr.open_dataset(files[0])
    if a.list:
        print(f"=== {files[0].name} の中身（全{len(files)}ファイル）===")
        print("データ変数:")
        for v in ds0.data_vars:
            print(f"  {v}  dims={ds0[v].dims}  shape={ds0[v].shape}")
        print("座標:", list(ds0.coords))
        print("次元:", dict(ds0.sizes))
        print("\n→ SIF 値の変数名を --var に、緯度/経度座標名を --lat-name/--lon-name に指定して再実行。")
        return

    var = a.var or _pick(list(ds0.data_vars), [v for v in ds0.data_vars if "sif" in v.lower()])
    if var is None or var not in ds0.data_vars:
        print(f"SIF 変数が特定できない。--var で指定。候補: {list(ds0.data_vars)}"); return
    latn = a.lat_name or _pick(list(ds0.coords) + list(ds0.dims), LAT_NAMES)
    lonn = a.lon_name or _pick(list(ds0.coords) + list(ds0.dims), LON_NAMES)
    if latn is None or lonn is None:
        print(f"緯度/経度座標が特定できない。--lat-name/--lon-name で指定。座標: {list(ds0.coords)}"); return
    timen = a.time_name or _pick(list(ds0.coords) + list(ds0.dims), TIME_NAMES)
    ds0.close()

    coords = pd.read_csv(a.coords).dropna(subset=["lat", "lon"])
    if a.sites:
        coords = coords[coords["site"].isin(a.sites)]
    if coords.empty:
        print("有効な座標が無い"); return
    print(f"  SIF変数={var}  緯度={latn} 経度={lonn} 時刻={timen or 'ファイル名から'}")
    print(f"  NetCDF {len(files)} 枚 × サイト {len(coords)} 件を抽出")

    if timen is not None:
        # 時刻座標つき（単一/複数ファイルを連結）
        ds = xr.open_mfdataset([str(f) for f in files], combine="by_coords") \
            if len(files) > 1 else xr.open_dataset(files[0])
        for r in coords.itertuples():
            try:
                s = ds[var].sel({latn: float(r.lat), lonn: float(r.lon)}, method="nearest")
                df = s.to_dataframe().reset_index()[[timen, var]].dropna()
                df.columns = ["date", "sif"]
                df["date"] = pd.to_datetime(df["date"])
                df.sort_values("date").to_csv(f"{r.site}_sif.csv", index=False)
                print(f"  [出力] {r.site}_sif.csv（{len(df)} 点）")
            except Exception as e:
                print(f"  {r.site} 抽出失敗: {type(e).__name__}: {e}")
        ds.close()
    else:
        # 時刻座標なし＝ファイル名から日付。各ファイルで全サイトの画素を1回読む
        if not a.date_regex:
            print("時刻座標が無い。--date-regex でファイル名の日付形式を指定。"); return
        rows = {r.site: [] for r in coords.itertuples()}
        for f in files:
            m = re.search(a.date_regex, f.name)
            if not m:
                continue
            g = m.groups()
            try:
                if a.date_fmt == "yyyymmdd":
                    d = pd.Timestamp(g[0][:4] + "-" + g[0][4:6] + "-" + g[0][6:8])
                elif a.date_fmt == "yyyyddd":
                    d = pd.Timestamp(year=int(g[0]), month=1, day=1) + pd.Timedelta(days=int(g[1]) - 1)
                else:
                    d = pd.Timestamp(year=int(g[0]), month=int(g[1]), day=int(g[2]))
            except (ValueError, IndexError):
                continue
            ds = xr.open_dataset(f)
            for r in coords.itertuples():
                try:
                    v = float(ds[var].sel({latn: float(r.lat), lonn: float(r.lon)},
                                          method="nearest").values)
                    rows[r.site].append((d, v))
                except Exception:
                    pass
            ds.close()
        for site, rr in rows.items():
            df = pd.DataFrame(rr, columns=["date", "sif"]).dropna().sort_values("date")
            df.to_csv(f"{site}_sif.csv", index=False)
            print(f"  [出力] {site}_sif.csv（{len(df)} 点）")

    print("\n  → 次： python research/sif_respiration_step38.py --site JP-Tak --sif JP-Tak_sif.csv --qc-max 1")
    print("  留保：TROPOMIは2018+(夏の重なり~4年)＝短いがほぼ日次(4日記憶を分解可)。画素~3.5×5.5kmでfootprint不一致は残る。")


if __name__ == "__main__":
    main()
