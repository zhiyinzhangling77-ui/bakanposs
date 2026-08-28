"""旗79：AmeriFlux FLUXNET 取得物の**下調べと座標照合**（登録の前・検定はしない）。

旗78 の一覧に従って AmeriFlux FLUXNET 版を取得した。だが**推測で `sites.py` に登録しない**
（旗68/69/76 と同じ作法）。特に **US-SSH が KAYE のチャンバーと同一地点か**は**未確認の推定**であり、
**座標で確かめてから**でなければ登録してはいけない——**旗64 で、座標が黙って壊れていた前例がある**。

本ツールは三つだけやる：

  1. **何が置かれているか**——サイトごとのファイル、**zip のままか展開済みか**、HH の実体
  2. **列名と変数マップの適合**——既存の2つのマップ（`japanflux` / `base`）でどれだけ埋まるか、
     **埋まらない変数の候補列**は何か（FLUXNET 版は `TS_F_MDS_1` のように**添字が付く**ことがある）
  3. **座標の照合**——BADM から緯度経度を取り、**COSORE の各チャンバーとの距離**を出す
     ＝**同一地点と呼べるか**を 10km 基準（旗51）で判定する

**この3つが揃って初めて登録する。** 揃わないものは「**手元では組めない**」と記録して終わりにする。

    python research/ameriflux_probe_step79.py --amf-dir /mnt/hdd/AmeriFlux_FLUXNET \
        --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from colocate_step51 import haversine

CODE_RE = re.compile(r"(?<![A-Za-z0-9])([A-Z]{2}-[A-Za-z0-9]{2,6})(?![A-Za-z0-9])")
NEED = ("Rg", "Ta", "VPD", "Ts", "P", "th", "gH", "gLE", "GER", "NEE", "GEP")
TOKENS = {"Rg": ("SW_IN",), "Ta": ("TA_",), "VPD": ("VPD",), "Ts": ("TS_",),
          "P": ("P_F", "P_ERA", "PRECIP"), "th": ("SWC",), "gH": ("H_F", "H_CORR"),
          "gLE": ("LE_F", "LE_CORR"), "GER": ("RECO",), "NEE": ("NEE",), "GEP": ("GPP",)}


def scan(amf_dir):
    """サイトコードごとに、ファイル・zip・HH の実体を集める。"""
    root = Path(amf_dir)
    per = defaultdict(lambda: {"files": [], "zips": [], "csv": [], "badm": []})
    for p in root.rglob("*"):
        if not p.is_file() or "__MACOSX" in p.parts or p.name.startswith("._"):
            continue
        m = CODE_RE.search(p.name)
        code = m.group(1) if m else None
        if code is None:
            continue
        d = per[code]
        d["files"].append(p)
        if p.suffix.lower() == ".zip":
            d["zips"].append(p)
        elif p.suffix.lower() == ".csv":
            d["csv"].append(p)
        if "BIF" in p.name.upper() or "BADM" in p.name.upper():
            d["badm"].append(p)
    return per


def read_header(path, inner=None):
    """CSV の列名を読む。``inner`` を渡すと zip の中から読む。"""
    try:
        if inner is None:
            return list(pd.read_csv(path, nrows=2).columns)
        with zipfile.ZipFile(path) as z, z.open(inner) as f:
            return list(pd.read_csv(f, nrows=2).columns)
    except Exception:
        return None


def find_hh(d):
    """HH（30分値）らしき実体を返す：(path, inner_name or None, 列名)。"""
    cands = sorted([p for p in d["csv"] if "_HH" in p.name.upper()],
                   key=lambda p: p.stat().st_size, reverse=True)
    cands += sorted([p for p in d["csv"] if p not in cands],
                    key=lambda p: p.stat().st_size, reverse=True)
    for p in cands[:3]:
        cols = read_header(p)
        if cols:
            return p, None, cols
    # zip のまま置かれている場合は**展開せずに中を読む**
    for z in sorted(d["zips"], key=lambda p: p.stat().st_size, reverse=True)[:3]:
        try:
            with zipfile.ZipFile(z) as zf:
                names = [n for n in zf.namelist() if n.lower().endswith(".csv")]
        except Exception:
            continue
        hh = [n for n in names if "_HH" in n.upper()] or names
        for n in sorted(hh, key=len)[:3]:
            cols = read_header(z, n)
            if cols:
                return z, n, cols
    return None, None, None


def badm_coords(paths, code):
    """BADM（BIF）から LOCATION_LAT / LOCATION_LONG を取る。**そのコードの行だけ**を見る。"""
    for p in paths:
        try:
            if p.suffix.lower() in (".xlsx", ".xls"):
                df = pd.read_excel(p)
            else:
                df = pd.read_csv(p)
        except Exception:
            continue
        cols = {c.upper(): c for c in df.columns}
        sid = next((cols[k] for k in ("SITE_ID", "SITEID") if k in cols), None)
        var = next((cols[k] for k in ("VARIABLE", "VARIABLE_NAME") if k in cols), None)
        val = next((cols[k] for k in ("DATAVALUE", "VALUE") if k in cols), None)
        if not (sid and var and val):
            continue
        sub = df[df[sid].astype(str).str.upper() == code.upper()]
        if sub.empty:
            continue
        def get(name):
            r = sub[sub[var].astype(str).str.upper() == name]
            if r.empty:
                return np.nan
            return pd.to_numeric(r[val], errors="coerce").dropna().astype(float).median()
        lat, lon = get("LOCATION_LAT"), get("LOCATION_LONG")
        if np.isfinite(lat) and np.isfinite(lon):
            return float(lat), float(lon), p.name
    return None, None, None


def main():
    p = argparse.ArgumentParser(description="AmeriFlux 取得物の下調べと座標照合")
    p.add_argument("--amf-dir", default="/mnt/hdd/AmeriFlux_FLUXNET")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--km", type=float, default=25.0)
    a = p.parse_args()

    print("=== 旗79：AmeriFlux 取得物の下調べと座標照合（登録の前・検定はしない）===")
    print("  **推測で sites.py に登録しない**（旗68/69/76 と同じ作法）。")
    print("  特に **US-SSH が KAYE と同一地点か**は未確認の推定＝**座標で確かめてから**でなければ")
    print("  登録してはいけない（旗64 で座標が黙って壊れていた前例がある）。\n")

    per = scan(a.amf_dir)
    if not per:
        print(f"  {a.amf_dir} にサイトコードを含むファイルが無い"); return
    print(f"  見つかったサイト：{len(per)} 件 {sorted(per)}\n")

    desc = pd.read_csv(Path(a.cosore_dir) / "description.csv")
    ch = []
    for _, d in desc.iterrows():
        try:
            clat, clon = float(d["CSR_LATITUDE"]), float(d["CSR_LONGITUDE"])
        except (TypeError, ValueError, KeyError):
            continue
        if np.isfinite(clat) and np.isfinite(clon):
            f = Path(a.cosore_dir) / "datasets" / f"data_{d['CSR_DATASET']}.csv"
            ch.append((str(d["CSR_DATASET"]), clat, clon, f.exists()))

    from japanflux_pn.sites import DEFAULT_VAR_MAP, DEFAULT_VAR_MAP_BASE
    for code in sorted(per):
        d = per[code]
        print(f"  ━━ {code} ━━")
        print(f"    ファイル {len(d['files'])} 件"
              f"（csv {len(d['csv'])}／zip {len(d['zips'])}／BADM {len(d['badm'])}）")
        # ① 中身
        path, inner, cols = find_hh(d)
        if cols is None:
            print("    **HH の実体を読めない**（zip の中も見た）\n"); continue
        where = f"{path.name}" + (f" 内 {inner}" if inner else "")
        print(f"    HH：{where}（{len(cols)} 列）")
        # ② 変数マップの適合
        best = None
        for name, vm in (("japanflux", DEFAULT_VAR_MAP), ("base", DEFAULT_VAR_MAP_BASE)):
            hit = [k for k, v in vm.items() if v in cols]
            miss = [k for k in vm if vm[k] not in cols]
            print(f"    マップ '{name}'：{len(hit)}/{len(vm)} 一致"
                  f"{'／欠け ' + ','.join(miss) if miss else ''}")
            if best is None or len(hit) > best[1]:
                best = (name, len(hit), miss)
        if best and best[2]:
            print(f"    → **欠けている変数の候補列**（FLUXNET 版は添字が付くことがある）：")
            for k in best[2]:
                cand = [c for c in cols if any(t in c.upper() for t in TOKENS.get(k, (k.upper(),)))]
                print(f"       {k:<4}: {cand[:8] if cand else '**候補なし**'}")
        # ③ 座標照合
        lat, lon, src = badm_coords(d["badm"], code)
        if lat is None:
            print("    **BADM から座標を取れない**＝同一地点かどうか判定できない")
        else:
            print(f"    座標：{lat:.5f}, {lon:.5f}（{src}）")
            near = sorted(((haversine(lat, lon, cl, cn), ds, ok) for ds, cl, cn, ok in ch))
            print(f"    **最も近い COSORE チャンバー**（{a.km:.0f}km 以内）：")
            shown = 0
            for km, ds, ok in near:
                if km > a.km:
                    break
                tag = ("**同一地点（10km 以内）**" if km <= 10 else "近いが 10km 超")
                print(f"       {km:>8.2f} km  {ds:<32}{'データあり' if ok else '**データ無し**'}  {tag}")
                shown += 1
                if shown >= 6:
                    break
            if shown == 0:
                print(f"       {a.km:.0f}km 以内に無し（最近傍 {near[0][0]:.1f} km ＝ {near[0][1]}）")
                print(f"       ＝**このタワーでは対を作れない**")
        print()

    print("  === 次の判断 ===")
    print("  **10km 以内にデータのあるチャンバーがある**サイトだけを `sites.py` に登録する。")
    print("  変数の欠けは **var_overrides** で埋める（上の候補列から**人が選ぶ**）。")
    print("  登録後は旗66（同一地点でタワー×チャンバー）を実行する。")
    print("  留保：")
    print("   ・**距離が近いことは同一地点の必要条件であって十分条件ではない**（旗51 と同じ注意）。")
    print("     林分・処理・設置年が違えば、同じ座標でも別の生態系である。")
    print("   ・FLUXNET 版は**ギャップフィル済み**＝旗46 が示した『穴埋めは冗長を押し上げる』が効く。")
    print("     メモリの検定（旗66/74）では駆動が穴埋めでも大きな問題にならないが、**記しておく**。")
    print("   ・**AmeriFlux のデータ利用方針**（CC-BY-4.0）に従い、公表時は各サイトPIへの")
    print("     クレジットと AmeriFlux の引用が要る。COSORE 側の条件（Bond-Lamberty 2020）も併存する。")


if __name__ == "__main__":
    main()
