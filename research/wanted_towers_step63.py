"""旗63：穴②を閉じるための「探すべきタワー」一覧を、COSORE 側から逆引きして作る。

旗51 は**手元のタワー（日本のみ）から** COSORE を探し「10km 以内に一組も無い」と結論した。
だがそれは**タワー側の品揃えの問題**であって、「同一地点が存在しない」ではない。
`DRY_SITE_EXPANSION.md` で気づいた通り、**COSORE のチャンバーの多くは、その場所にフラックスタワーがある**
（例：SCOTT_WKG は AmeriFlux US-Wkg と同一研究者・同一名称）。

そこで向きを逆にする：**COSORE 側から「この座標にタワーがあれば弧が閉じる」一覧を作る**。
優先順は本研究の中心的主張に効く順＝**★短メモリを示すサイトが最優先**（旗53 の較正済み基準）。

併せて **旗51 の積み残し**も出す：手元のタワーとの**最近傍距離**（10km で切らずに全部見る）。
JP-Tef と UEYAMA_TESHIO は名前が一致するのに組にならなかったが、**何km離れていたのか**を見ていなかった。

    python research/wanted_towers_step63.py --cosore-dir /mnt/hdd/cosore-0.7.0 \
        --sites JP-Tak JP-Fhk JP-Tmd JP-Fjy JP-Tef JP-Ta2 JP-Mse JP-BBY JP-SMF JP-Spp
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
from chamber_memory_recount_step53 import detect
from memory_attribution_flex_step54 import ACF_THR, EFOLD_MAX
from colocate_step51 import haversine, tower_coords


def main():
    p = argparse.ArgumentParser(description="穴②を閉じるための探すべきタワー一覧")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--sites", nargs="+", default=[])
    p.add_argument("--data-dir", default="/mnt/hdd/JAPANFLUX")
    p.add_argument("--coords", default=None)
    p.add_argument("--top", type=int, default=25)
    a = p.parse_args()

    root = Path(a.cosore_dir); desc = pd.read_csv(root / "description.csv")
    tw = tower_coords(a.sites, a.data_dir, a.coords) if a.sites else {}
    print("=== 旗63：穴②を閉じるための『探すべきタワー』一覧（COSORE 側から逆引き）===")
    print(f"  手元のタワー座標：{len(tw)} 件（最近傍距離の算出に使う）")
    print(f"  優先順＝★短メモリ（非線形基底・R²≥0.3・ACF1≥{ACF_THR}・e-fold≤{EFOLD_MAX}日）を上に。\n")

    rows = []
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            lat, lon = float(d["CSR_LATITUDE"]), float(d["CSR_LONGITUDE"])
        except (TypeError, ValueError, KeyError):
            continue
        if not (np.isfinite(lat) and np.isfinite(lon)):
            continue
        try:
            df, st, sm = load_cosore(f, None)
            if "Tsoil" not in df:
                continue
            cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
            dd = df[cols].copy()
            daily = dd.groupby(dd.index.normalize()).mean()
            daily = daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D"))
            flx = detect(daily, "Rs", "Tsoil", "SM" if "SM" in df else None, "flex")
        except Exception:
            continue
        if not flx or not np.isfinite(flx.get("r2", np.nan)):
            continue
        judged = flx["r2"] >= 0.3 and np.isfinite(flx["acf1"]) and np.isfinite(flx["efold"])
        star = judged and flx["acf1"] >= ACF_THR and flx["efold"] <= EFOLD_MAX
        near_km, near_site = np.inf, "—"
        for s, (tlat, tlon) in tw.items():
            km = haversine(tlat, tlon, lat, lon)
            if km < near_km:
                near_km, near_site = km, s
        rows.append({"ds": ds, "igbp": ig, "lat": lat, "lon": lon,
                     "star": int(star), "judged": int(judged),
                     "acf1": flx["acf1"], "efold": flx["efold"],
                     "near_km": near_km, "near_site": near_site})

    if not rows:
        print("  対象なし"); return
    df = pd.DataFrame(rows).sort_values(["star", "judged", "acf1"], ascending=False)

    print(f"  {'優先':<4}{'dataset':<32}{'IGBP':<15}{'緯度':>8}{'経度':>9}"
          f"{'ACF1':>7}{'ef':>4}  最近傍の手元タワー")
    for _, r in df.head(a.top).iterrows():
        mark = "★" if r["star"] else ("·" if r["judged"] else "△")
        near = f"{r['near_site']} {r['near_km']:.0f}km" if np.isfinite(r["near_km"]) else "—"
        print(f"  {mark:<4}{r['ds']:<32}{r['igbp'][:13]:<15}{r['lat']:>8.3f}{r['lon']:>9.3f}"
              f"{r['acf1']:>7.2f}{r['efold']:>4.0f}  {near}")

    n_star = int(df["star"].sum())
    print(f"\n  === まとめ ===")
    print(f"  ★短メモリのサイト：{n_star} 件＝**この座標にタワーがあれば、弧が閉じる候補**。")
    if np.isfinite(df["near_km"]).any():
        print(f"  手元タワーとの最近傍距離の最小値：{df['near_km'].min():.0f} km"
              f"（{df.loc[df['near_km'].idxmin(), 'ds']}）")
        print("  ＝旗51 が「10km 以内なし」とした中身。**近い組が一つも無い**のか、"
              "**惜しい組があった**のかがこれで分かる。")
    print("\n  次にやること：★の座標を FLUXNET/AmeriFlux/ICOS のサイト一覧と突き合わせ、")
    print("  同一地点のタワーがあるものを特定する。BASE 形式なら本リポジトリは対応済み（改修不要）。")
    print("  留保：座標が近くても林分・処理が違えば同一地点ではない（旗51 と同じ注意）。")
    print("        チャンバーの座標は代表点であり、タワーのフットプリントと一致する保証もない。")


if __name__ == "__main__":
    main()
