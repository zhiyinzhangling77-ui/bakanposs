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
    """セル全体が数値のものだけを採る（"decimal deg ref WGS84" の 84 等の埋め込み数字を弾く）。"""
    for c in cells:
        s = str(c).strip()
        try:
            return float(s)
        except ValueError:
            continue
    return None


def _find_badm(data_dir):
    """BADM ファイル(xlsx優先=JapanFlux2024, csv も)。Site_General を最優先。"""
    root = Path(data_dir)
    pats = ["**/*BADM*Site_General*.xlsx", "**/*BADM*.xlsx",
            "**/*BADM*Site_General*.csv", "**/*BADM*.csv",
            "**/*BIF*.xlsx", "**/*BIF*.csv"]
    out = []
    for p in pats:
        for f in sorted(root.glob(p)):
            if f not in out and not f.name.startswith("~$"):
                out.append(f)
    return out


def _read_any(f):
    """xlsx/csv を header なし・全文字列で読む（BADM の総当たり走査用）。"""
    import pandas as pd
    if str(f).lower().endswith((".xlsx", ".xls")):
        return pd.read_excel(f, header=None, dtype=str)
    return pd.read_csv(f, header=None, dtype=str, on_bad_lines="skip")


def coords_from_badm(data_dir, code):
    """BADM 群から (lat, lon) を抽出。見つからねば (None, None)。"""
    lat = lon = None
    for f in _find_badm(data_dir):
        try:
            df = _read_any(f)
        except Exception:
            continue
        for _, row in df.iterrows():
            cells = [str(x) for x in row.tolist() if x is not None and str(x) != "nan"]
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
        # 妥当域チェック（緯度[-90,90]・経度[-180,180]外は疑い→手入力扱い）
        if lat is not None and not (-90 <= lat <= 90):
            lat = None
        if lon is not None and not (-180 <= lon <= 180):
            lon = None
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
