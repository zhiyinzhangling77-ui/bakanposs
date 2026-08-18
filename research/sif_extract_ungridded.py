"""TROPOMI ungridded L2 SIF（Caltech 等）から、タワー近傍 sounding を集約して `<site>_sif.csv` を作る。

真SIF(実測)を良い footprint で得るには、Caltech の無料公開＝**ungridded L2（軌道ごとの sounding 点群）**を使う。
格子でないので「タワー画素」を切れない→**タワー(lat,lon)から半径 N km 内の sounding を拾って日次平均**する。
gridded 公開品(TROPOSIF 0.2°≒22km)は GOSIF より粗いので、良い footprint はこの ungridded 経路のみ。

配布の変数名は多様なので、**変数・座標を自動検出、不明なら --list で一覧**（捏造せず実物に合わせる）。
時刻は変数(CF units)から、無ければファイル名から。品質フラグがあれば --qc-var/--qc-max で絞る。

    python research/sif_extract_ungridded.py --ncdir /path/to/L2 --list          # 中身を見る
    python research/sif_extract_ungridded.py --coords site_coords.csv --ncdir /path/to/L2 \
        --var sif --lat-var lat --lon-var lon --time-var TIME --radius 20

入手(ローカル): Caltech ungridded TROPOMI SIF  ftp://fluo.gps.caltech.edu/data/tropomi/ （L2, 2018-03〜2021-07）。
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

LATV = ["lat", "latitude", "sounding_latitude", "Latitude", "LAT", "clat"]
LONV = ["lon", "longitude", "sounding_longitude", "Longitude", "LON", "clon"]
TIMEV = ["TIME", "time", "Time", "delta_time", "datetime"]
SIFV = ["sif", "SIF", "SIF_743", "sif_743", "SIF_740", "SIF_Corr", "sif_dc"]


def haversine_km(lat0, lon0, lats, lons):
    """1点(lat0,lon0)から配列(lats,lons)への大圏距離[km]。"""
    R = 6371.0
    p0 = np.radians(lat0); l0 = np.radians(lon0)
    p1 = np.radians(np.asarray(lats, float)); l1 = np.radians(np.asarray(lons, float))
    dp = p1 - p0; dl = l1 - l0
    a = np.sin(dp / 2) ** 2 + np.cos(p0) * np.cos(p1) * np.sin(dl / 2) ** 2
    return 2 * R * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def _pick(names, cands):
    low = {n.lower(): n for n in names}
    for c in cands:
        if c in names:
            return c
        if c.lower() in low:
            return low[c.lower()]
    return None


def _file_date(name, rx, fmt):
    if not rx:
        return None
    m = re.search(rx, name)
    if not m:
        return None
    g = m.groups()
    try:
        if fmt == "yyyymmdd":
            s = g[0]; return pd.Timestamp(f"{s[:4]}-{s[4:6]}-{s[6:8]}")
        if fmt == "yyyyddd":
            return pd.Timestamp(year=int(g[0]), month=1, day=1) + pd.Timedelta(days=int(g[1]) - 1)
        if fmt == "ymd":
            return pd.Timestamp(year=int(g[0]), month=int(g[1]), day=int(g[2]))
    except (ValueError, IndexError):
        return None


def main():
    p = argparse.ArgumentParser(description="ungridded TROPOMI L2 SIF をタワー近傍で集約")
    p.add_argument("--coords", default="site_coords.csv")
    p.add_argument("--ncdir", required=True)
    p.add_argument("--glob", default="*.nc")
    p.add_argument("--var", default=None, help="SIF 変数名（既定=自動）")
    p.add_argument("--lat-var", default=None)
    p.add_argument("--lon-var", default=None)
    p.add_argument("--time-var", default=None, help="時刻変数（CF units 対応）。無ければ --date-regex")
    p.add_argument("--date-regex", default=None, help="時刻変数が無いとき、ファイル名から日付")
    p.add_argument("--date-fmt", default="yyyymmdd", choices=["yyyymmdd", "yyyyddd", "ymd"])
    p.add_argument("--radius", type=float, default=20.0, help="タワーからの許容半径[km]（既定20）")
    p.add_argument("--qc-var", default=None, help="品質フラグ変数（あれば）")
    p.add_argument("--qc-max", type=float, default=None, help="qc-var ≤ この値のみ採用")
    p.add_argument("--sites", nargs="+")
    p.add_argument("--list", action="store_true")
    a = p.parse_args()

    try:
        import netCDF4
    except ImportError:
        print("netCDF4 が要ります： pip install netCDF4"); return

    files = sorted(Path(a.ncdir).glob(a.glob))
    if not files:
        print(f"NetCDF が無い： {a.ncdir}/{a.glob}"); return

    ds0 = netCDF4.Dataset(files[0])
    allvars = list(ds0.variables)
    if a.list:
        print(f"=== {files[0].name} の変数（全{len(files)}ファイル）===")
        for v in allvars:
            var = ds0.variables[v]
            u = getattr(var, "units", "")
            print(f"  {v}  dims={var.dimensions}  shape={var.shape}  units={u}")
        print("\n→ SIF値/緯度/経度/時刻 の変数名を --var/--lat-var/--lon-var/--time-var に指定して再実行。")
        ds0.close(); return

    var = a.var or _pick(allvars, SIFV)
    latv = a.lat_var or _pick(allvars, LATV)
    lonv = a.lon_var or _pick(allvars, LONV)
    timev = a.time_var or _pick(allvars, TIMEV)
    ds0.close()
    if not all([var, latv, lonv]):
        print(f"変数が特定できない。--var/--lat-var/--lon-var を指定。候補: {allvars}"); return

    coords = pd.read_csv(a.coords).dropna(subset=["lat", "lon"])
    if a.sites:
        coords = coords[coords["site"].isin(a.sites)]
    if coords.empty:
        print("有効な座標が無い"); return
    print(f"  SIF={var} lat={latv} lon={lonv} time={timev or 'ファイル名'} 半径={a.radius}km")
    print(f"  L2 {len(files)} 枚 × サイト {len(coords)} 件を集約")

    rows = {r.site: [] for r in coords.itertuples()}
    for i, f in enumerate(files, 1):
        try:
            ds = netCDF4.Dataset(f)
            lat = np.asarray(ds.variables[latv][:]).ravel()
            lon = np.asarray(ds.variables[lonv][:]).ravel()
            sif = np.asarray(ds.variables[var][:], float).ravel()
            n = min(len(lat), len(lon), len(sif))
            lat, lon, sif = lat[:n], lon[:n], sif[:n]
            qc = None
            if a.qc_var and a.qc_var in ds.variables and a.qc_max is not None:
                qc = np.asarray(ds.variables[a.qc_var][:], float).ravel()[:n]
            # 時刻
            times = None
            if timev and timev in ds.variables:
                tv = ds.variables[timev]
                try:
                    times = pd.to_datetime(netCDF4.num2date(
                        np.asarray(tv[:]).ravel()[:n], getattr(tv, "units", ""),
                        getattr(tv, "calendar", "standard"),
                        only_use_cftime_datetimes=False).astype("datetime64[ns]"))
                except Exception:
                    times = None
            fdate = _file_date(f.name, a.date_regex, a.date_fmt)
            for r in coords.itertuples():
                d = haversine_km(float(r.lat), float(r.lon), lat, lon)
                m = np.isfinite(sif) & np.isfinite(d) & (d <= a.radius)
                if qc is not None:
                    m &= (qc <= a.qc_max)
                if not m.any():
                    continue
                if times is not None:
                    for t, s in zip(np.asarray(times)[m], sif[m]):
                        rows[r.site].append((pd.Timestamp(t).normalize(), float(s)))
                elif fdate is not None:
                    for s in sif[m]:
                        rows[r.site].append((fdate, float(s)))
            ds.close()
        except Exception as e:
            print(f"  {f.name} 読めず: {type(e).__name__}: {e}")
        if i % 20 == 0:
            print(f"    …{i}/{len(files)}")

    for site, rr in rows.items():
        if not rr:
            print(f"  {site}: 半径内 sounding なし（--radius を広げるか座標確認）"); continue
        df = pd.DataFrame(rr, columns=["date", "sif"])
        daily = df.groupby("date")["sif"].mean().reset_index().sort_values("date")
        daily.to_csv(f"{site}_sif.csv", index=False)
        print(f"  [出力] {site}_sif.csv（{len(daily)}日, sounding {len(df)}）")

    print("\n  → 次： python research/sif_respiration_step38.py --site JP-Tak --sif JP-Tak_sif.csv --qc-max 1")
    print("  留保：2018+(~4夏)・半径集約で空間平均・雲/品質フラグは --qc-var で。footprint不一致は残るが真SIF・ほぼ日次。")


if __name__ == "__main__":
    main()
