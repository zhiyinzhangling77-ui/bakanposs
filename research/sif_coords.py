"""SIF 抽出の前提：各サイトのタワー座標(lat/lon)を BADM から取り出す（無ければ手入力テンプレート）。

SIF ピクセルをタワー位置で切り出すには lat/lon が要る。我々のコードには座標が無いので、
配布メタデータ BADM（LOCATION_LAT / LOCATION_LONG）から抽出する。BADM がローカルに無いサイトは
**空欄のテンプレート行**を出す（AsiaFlux/FLUXNET のサイト情報から手入力する。座標の捏造は不可）。

    python research/sif_coords.py --sites JP-Tak JP-Ta2 JP-Mse JP-BBY CN-HaM MN-Hst --out site_coords.csv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

LAT_RE = re.compile(r"LOCATION_LAT", re.IGNORECASE)
LON_RE = re.compile(r"LOCATION_LONG", re.IGNORECASE)
_NUM = re.compile(r"[-+]?\d+\.?\d*")


def _num_in(cells):
    for c in cells:
        m = _NUM.search(str(c))
        if m:
            try:
                return float(m.group())
            except ValueError:
                continue
    return None


def coords_from_badm(data_dir, code):
    """BADM 群から (lat, lon) を抽出。見つからねば (None, None)。"""
    import pandas as pd
    from japanflux_pn.ecosystem import _iter_badm_files
    lat = lon = None
    for f in _iter_badm_files(data_dir, code):
        try:
            df = pd.read_csv(f, header=None, dtype=str, on_bad_lines="skip")
        except Exception:
            continue
        for _, row in df.iterrows():
            cells = [str(x) for x in row.tolist() if x is not None]
            joined = " ".join(cells)
            if lat is None and LAT_RE.search(joined):
                lat = _num_in([c for c in cells if not LAT_RE.search(c)])
            if lon is None and LON_RE.search(joined):
                lon = _num_in([c for c in cells if not LON_RE.search(c)])
        if lat is not None and lon is not None:
            break
    return lat, lon


def main():
    from japanflux_pn.sites import get_site
    p = argparse.ArgumentParser(description="タワー座標を BADM から抽出（無ければ空欄テンプレート）")
    p.add_argument("--sites", nargs="+", required=True)
    p.add_argument("--out", default="site_coords.csv")
    a = p.parse_args()

    rows = []
    print(f"  {'サイト':<8} {'lat':>10} {'lon':>10}  出所")
    for s in a.sites:
        lat = lon = None
        try:
            spec = get_site(s)
            lat, lon = coords_from_badm(spec.data_dir, s)
        except Exception as e:
            print(f"  {s:<8} SKIP {type(e).__name__}: {e}")
        src = "BADM" if (lat is not None and lon is not None) else "★要手入力"
        rows.append((s, lat, lon, src))
        ls = f"{lat:.5f}" if lat is not None else ""
        os_ = f"{lon:.5f}" if lon is not None else ""
        print(f"  {s:<8} {ls:>10} {os_:>10}  {src}")

    with open(a.out, "w") as f:
        f.write("site,lat,lon,source\n")
        for s, lat, lon, src in rows:
            f.write(f"{s},{'' if lat is None else lat},{'' if lon is None else lon},{src}\n")
    print(f"\n  [出力] {a.out}")
    miss = [r[0] for r in rows if r[1] is None or r[2] is None]
    if miss:
        print(f"  ★手入力が要るサイト（AsiaFlux/FLUXNET サイト情報から lat,lon を埋める）: {miss}")
        print("    ※座標は一次情報から。推測・捏造は不可。")


if __name__ == "__main__":
    main()
