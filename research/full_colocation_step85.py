"""旗85：**★で選ばない総当たり**——私が作った選択バイアスの大きさを確定する。

旗84 で、同一地点の対を**取得動機**で分けると一致率が
**★で選んだ 82%／付いてきた 25%** と分かれた。**方向は明確だが大きさは決められない**——
「付いてきた」群が乾燥地に偏り、**生態系型と交絡**しているためである（森林だけだと n=2）。

**BADM Only（全 AmeriFlux サイトの座標）があれば、これを確定できる。**

## なぜ座標だけで足りるのか

**タワー側はほぼ常に★だった**（旗66 の 23 組でタワーが★でなかったのは CA-TPD と US-SRM だけ）。
＝**一致率 ≒ チャンバーが★である確率**。
＝**「対を作れた全組のチャンバー★率」が、偏りのない一致率の推定になる**。
**タワーの流束データを落とさなくても測れる。**

## やること

  1. BADM から **AmeriFlux 全サイトの座標**を読む
  2. **COSORE 全チャンバー × AmeriFlux 全サイト**の距離を総当たり
  3. **10km 以内の組を全部列挙**＝**対を作れた可能性のある組**
  4. その**チャンバー★率**（＝偏りのない一致率の推定）と、
     **私が実際に取得した対の★率**を比べる
  5. **取得しなかった組が何組あるか**も出す

## 先に決めておく読み方

  ・**全組の★率と、取得した組の★率がほぼ同じ** → **選択バイアスは小さかった**。旗66 の結論は概ね保つ。
  ・**取得した組の方が明確に高い** → **その差が押し上げ分**。旗66 の「両側★ 9 組」を**その分割り引く**。

    python research/full_colocation_step85.py --cosore-dir /mnt/hdd/cosore-0.7.0 \
        --badm-dir /mnt/hdd/AmeriFlux_BADM
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore
from model_richness_step74 import measure, star
from colocate_step51 import haversine
from same_site_arc_step66 import PAIRS

ACQUIRED = {ds for _, ds, _ in PAIRS}          # 私が実際に対にしたチャンバー


def read_all_sites(badm_dir, max_files=2000):
    """BADM から**全サイトの座標**を読む。

    取得の仕方で**2 通りの形**がありうる：
      ・**multi-site BADM が 1 ファイル**（全サイトが 1 つに入っている）
      ・**サイトごとに 1 ファイル**（837 個に分かれている）
    **どちらでも動くよう、見つかった BADM を全部読んで統合する**
    （第1版は「一番大きい 1 ファイル」しか読まず、**分割配布だと 1 サイトしか取れなかった**）。
    """
    root = Path(badm_dir)
    cands = sorted([p for p in root.rglob("*")
                    if p.is_file() and p.suffix.lower() in (".csv", ".xlsx", ".xls")
                    and ("BIF" in p.name.upper() or "BADM" in p.name.upper())])
    if not cands:
        return 0, pd.DataFrame()
    parts = []
    n_read = 0
    for f in cands[:max_files]:
        try:
            if f.suffix.lower() in (".xlsx", ".xls"):
                df = pd.read_excel(f)
            else:
                df = pd.read_csv(f, low_memory=False, encoding="latin-1",
                                 encoding_errors="replace")
        except Exception:
            continue
        cols = {c.upper(): c for c in df.columns}
        sid = next((cols[k] for k in ("SITE_ID", "SITEID") if k in cols), None)
        var = next((cols[k] for k in ("VARIABLE", "VARIABLE_NAME") if k in cols), None)
        val = next((cols[k] for k in ("DATAVALUE", "VALUE") if k in cols), None)
        if not (sid and var and val):
            continue
        v = df[var].astype(str).str.upper()
        keep = df[v.isin(["LOCATION_LAT", "LOCATION_LONG", "IGBP"])]
        if keep.empty:
            continue
        parts.append(pd.DataFrame({"site": keep[sid].astype(str),
                                   "var": v[keep.index], "val": keep[val]}))
        n_read += 1
    if not parts:
        return n_read, pd.DataFrame()
    allp = pd.concat(parts, ignore_index=True)
    # **同じサイトが複数ファイルに出ることがある**（版違い等）＝最初の値を採る
    piv = (allp.drop_duplicates(subset=["site", "var"], keep="first")
              .pivot(index="site", columns="var", values="val"))
    out = pd.DataFrame(index=piv.index)
    out["lat"] = pd.to_numeric(piv.get("LOCATION_LAT"), errors="coerce")
    out["lon"] = pd.to_numeric(piv.get("LOCATION_LONG"), errors="coerce")
    out["igbp"] = piv.get("IGBP")
    return n_read, out.dropna(subset=["lat", "lon"])


def main():
    p = argparse.ArgumentParser(description="★で選ばない総当たり")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--badm-dir", required=True)
    p.add_argument("--km", type=float, default=10.0)
    a = p.parse_args()

    print("=== 旗85：★で選ばない総当たり（選択バイアスの大きさを確定する）===")
    print("  **タワー側はほぼ常に★**だった（旗66 の 23 組で例外は 2 件）。")
    print("  ＝**一致率 ≒ チャンバーが★である確率**")
    print("  ＝**対を作れた全組のチャンバー★率が、偏りのない一致率の推定**になる。\n")

    n_files, sites = read_all_sites(a.badm_dir)
    if sites.empty:
        print(f"  {a.badm_dir} に座標を含む BADM が無い（読めたファイル {n_files} 件）"); return
    print(f"  BADM ファイル {n_files} 件を統合 → **AmeriFlux サイト {len(sites)} 件**")
    if len(sites) < 100:
        print(f"  ※**100 サイト未満**＝分割配布の一部しか読めていない可能性がある。")
        print(f"    置き場所（--badm-dir）に全ファイルが在るか確認すること。")
    print()

    desc = pd.read_csv(Path(a.cosore_dir) / "description.csv")
    root = Path(a.cosore_dir)
    rows = []
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            clat, clon = float(d["CSR_LATITUDE"]), float(d["CSR_LONGITUDE"])
        except (TypeError, ValueError, KeyError):
            continue
        if not (np.isfinite(clat) and np.isfinite(clon)):
            continue
        dist = np.array([haversine(clat, clon, la, lo)
                         for la, lo in zip(sites["lat"], sites["lon"])])
        k = int(np.argmin(dist))
        near_km = float(dist[k]); near = sites.index[k]
        if near_km > a.km:
            rows.append({"ds": ds, "igbp": ig, "km": near_km, "tower": near,
                         "pairable": False, "star": None, "acq": ds in ACQUIRED})
            continue
        # **対が作れる**＝チャンバー側の判定を計算する（タワーは要らない）
        try:
            df, st, sm = load_cosore(f, None)
        except Exception:
            continue
        s = None
        if df is not None and "Tsoil" in df and "Rs" in df:
            cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
            daily = df[cols].groupby(df.index.normalize()).mean()
            if len(daily) >= 60:
                daily = daily.reindex(pd.date_range(daily.index.min(),
                                                    daily.index.max(), freq="D"))
                m = measure(daily["Rs"].to_numpy(), daily["Tsoil"].to_numpy(),
                            daily["SM"].to_numpy() if "SM" in daily else None,
                            "テンソルビン", True)
                s = star(m)
        rows.append({"ds": ds, "igbp": ig, "km": near_km, "tower": near,
                     "pairable": True, "star": s, "acq": ds in ACQUIRED})
    t = pd.DataFrame(rows)
    if t.empty:
        print("  チャンバーが読めない"); return

    pa = t[t["pairable"]]
    print(f"  {'chamber':<32}{'最近傍タワー':<10}{'km':>7}{'★':>4}  取得")
    for _, r in pa.sort_values("km").iterrows():
        mk = "★" if r["star"] else ("·" if r["star"] is False else "?")
        print(f"  {r['ds']:<32}{str(r['tower']):<10}{r['km']:>7.2f}{mk:>4}  "
              f"{'済' if r['acq'] else ''}")

    print(f"\n  === まとめ ===")
    print(f"  COSORE チャンバー {len(t)} 件のうち、**{a.km:.0f}km 以内にタワーがある**："
          f"**{len(pa)} 件**（＝対を作れた可能性のある組）")
    print(f"  そのうち**私が実際に取得したのは {int(pa['acq'].sum())} 件**"
          f"／**取得しなかったのは {int((~pa['acq']).sum())} 件**")
    for lab, sub in (("**全組**（偏りのない推定）", pa),
                     ("うち取得した組", pa[pa["acq"]]),
                     ("うち取得しなかった組", pa[~pa["acq"]])):
        j = sub[sub["star"].notna()]
        if len(j):
            print(f"    {lab:<28} 判定できた {len(j):>3} 件／★ {int(j['star'].sum()):>3} 件"
                  f"＝**{j['star'].mean():.0%}**")
    # 森林だけでも同じ比較（生態系型の交絡を外す）
    fo = pa[pa["igbp"].astype(str).str.lower().str.contains("forest|plantation")]
    print(f"\n  **森林だけ**（生態系型の交絡を外す）：{len(fo)} 件")
    for lab, sub in (("全組", fo), ("取得した組", fo[fo["acq"]]),
                     ("取得しなかった組", fo[~fo["acq"]])):
        j = sub[sub["star"].notna()]
        if len(j):
            print(f"    {lab:<20} 判定できた {len(j):>3} 件／★ {int(j['star'].sum()):>3} 件"
                  f"＝**{j['star'].mean():.0%}**")

    print("\n  === 読み方（先に決めてある）===")
    print("  ・**全組の★率と取得した組の★率がほぼ同じ** → **選択バイアスは小さかった**。")
    print("    旗66 の「両側★ 9 組」は概ね保てる。")
    print("  ・**取得した組の方が明確に高い** → **その差が押し上げ分**。**その分割り引く**。")
    print("  留保：")
    print("   ・**AmeriFlux は南北アメリカだけ**。日本（JapanFlux）・欧州（ICOS）のタワーは含まれず、")
    print("     **総当たりはアメリカ大陸について完全、それ以外は不完全**である。")
    print("   ・「タワーが 10km 以内にある」は**対が作れる必要条件**であって十分条件ではない")
    print("     （観測期間が重ならなければ検定できない＝旗66 で US-Ha1 が実際にそうだった）。")
    print("   ・**一致率 ≒ チャンバー★率**という近似は、**タワーがほぼ常に★**という観察に基づく。")
    print("     タワーが★でない場所（旗66 では CA-TPD・US-SRM）では成り立たない。")


if __name__ == "__main__":
    main()
