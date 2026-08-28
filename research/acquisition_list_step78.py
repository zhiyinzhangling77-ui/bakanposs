"""旗78：**いま取りに行くべきタワーの優先順位表**——旗73〜77 を踏まえて作り直す。

旗63 は同じ逆引きをしたが、**旗53 の内挿基準**で優先順位を付けていた。
その後：

  ・**旗70/71**：水マスター（A-3）は**同じモンゴル内の 113km 先ですら再現しなかった**
    ＝乾燥サイトを増やす目的は「一般化の確認」から「**何がレジームを分けるのか**」に変わった。
  ・**旗73〜77**：A-1（チャンバー呼吸の多日メモリ）は**八つの対抗仮説を跳ね返した**
    ＝**本研究の主軸は A-1** であり、**同一地点の対を増やすことが最も効く**。
  ・**旗74**：判定は**外挿残差**で行うべきと分かった（内挿 22/45 → 外挿 15/44）。

＝**優先順位を付け直す必要がある**。本ツールは COSORE 側から逆引きして、
**外挿基準で★が出るチャンバー**を、**取得の価値が高い順**に並べる。

## 何を出すか

各チャンバーについて：
  ・**外挿残差での判定**（旗74 と同じ：テンソルビン・時間ブロックCV・R²≥0.3・ACF1≥0.64・e-fold≤7日）
  ・記録の長さ（日数・年数）＝**長いほど検定が安定する**
  ・座標と、そこから推定される**地域**
  ・**探すべきネットワーク**（座標からの推定＝**確認が要る**）
  ・**既に対が在るか**（旗66 の4組）

## 優先順位の考え方（先に書いておく）

  1. **外挿で★が出る森林**——A-1 の同一地点検証を増やす。**最優先**。
  2. **記録が長い**——旗74/75/77 はどれも日数不足で判定不能になったサイトがあった。
  3. **既存の対と地理的に離れている**——4組はすべて日本＝**擬似反復の懸念がある**。
     北米・欧州が1組でも入れば、**「同一地点で再現する」が日本に固有でない**と言える。

**乾燥サイト（US-Wkg 等）は優先度を下げる**：旗70/71 で水マスターの一般化そのものが
不確かになったため、**同じ枠組みでの追試より、レジームを分ける要因の解明が先**である。
ただし**取得の手間は同じ**なので、ついでに取る価値はある（そう記す）。

    python research/acquisition_list_step78.py --cosore-dir /mnt/hdd/cosore-0.7.0
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
from same_site_arc_step66 import PAIRS

HAVE = {ds for _, ds, _ in PAIRS}          # 既に対が在るチャンバー


def region_of(lat, lon):
    """座標から地域と、**探すべきネットワークの推定**を返す（推定＝確認が要る）。"""
    if -170 <= lon <= -50 and 10 <= lat <= 75:
        return "北米", "AmeriFlux"
    if -12 <= lon <= 40 and 35 <= lat <= 72:
        return "欧州", "ICOS / European Fluxes"
    if 120 <= lon <= 150 and 24 <= lat <= 46:
        return "日本・東アジア", "JapanFlux / AsiaFlux"
    if 60 <= lon <= 150 and -12 <= lat <= 55:
        return "アジア", "AsiaFlux"
    if -85 <= lon <= -30 and -56 <= lat <= 13:
        return "中南米", "AmeriFlux（LBA 等）"
    if 110 <= lon <= 155 and -45 <= lat <= -10:
        return "豪州", "OzFlux"
    if -20 <= lon <= 52 and -35 <= lat <= 20:
        return "アフリカ", "（該当網は要確認）"
    return "その他", "（要確認）"


def main():
    p = argparse.ArgumentParser(description="取得すべきタワーの優先順位表")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--igbp", default="forest")
    p.add_argument("--top", type=int, default=20)
    a = p.parse_args()
    root = Path(a.cosore_dir)
    desc = pd.read_csv(root / "description.csv")

    print("=== 旗78：いま取りに行くべきタワーの優先順位表（旗73〜77 を踏まえて）===")
    print("  判定は**旗74 と同じ外挿基準**（テンソルビン・時間ブロックCV）。")
    print("  旗63 は旗53 の**内挿**基準で並べていたので、**付け直す**。\n")

    rows = []
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        if a.igbp and a.igbp.lower() not in ig.lower():
            continue
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
        except Exception:
            continue
        if df is None or "Tsoil" not in df or "Rs" not in df:
            continue
        cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
        daily = df[cols].groupby(df.index.normalize()).mean()
        daily = daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D"))
        y = daily["Rs"].to_numpy(); T = daily["Tsoil"].to_numpy()
        W = daily["SM"].to_numpy() if "SM" in daily else None
        m = measure(y, T, W, "テンソルビン", True)
        s = star(m)
        nd = int(np.isfinite(y).sum())
        yrs = (daily.index.max() - daily.index.min()).days / 365.25
        reg, net = region_of(lat, lon)
        rows.append({"ds": ds, "star": s, "acf1": (m or {}).get("acf1", np.nan),
                     "n": nd, "yrs": yrs, "lat": lat, "lon": lon,
                     "reg": reg, "net": net, "have": ds in HAVE})
    if not rows:
        print("  対象が無い"); return
    df2 = pd.DataFrame(rows)
    # 優先順位：★（外挿）→ 既存の対でない → 日数が多い
    df2["rank"] = (df2["star"].fillna(False).astype(int) * 100
                   + (~df2["have"]).astype(int) * 10
                   + np.clip(df2["n"] / 500.0, 0, 9))
    df2 = df2.sort_values("rank", ascending=False)

    print(f"  {'chamber':<30}{'外挿★':>6}{'ACF1':>7}{'日数':>6}{'年':>5}"
          f"  {'緯度':>7}{'経度':>9}  {'地域':<12}{'探すネットワーク（推定）':<24}既存")
    for _, r in df2.head(a.top).iterrows():
        mark = "★" if r["star"] else ("·" if r["star"] is False else "?")
        print(f"  {r['ds']:<30}{mark:>6}"
              f"{(f'{r.acf1:.2f}' if np.isfinite(r.acf1) else '—'):>7}"
              f"{r['n']:>6}{r['yrs']:>5.1f}  {r['lat']:>7.2f}{r['lon']:>9.2f}  "
              f"{r['reg']:<12}{r['net']:<24}{'済' if r['have'] else ''}")

    got = df2[df2["star"] == True]
    print(f"\n  === まとめ ===")
    print(f"  外挿基準で★のチャンバー：{len(got)} 件／うち**まだ対が無い**：{int((~got['have']).sum())} 件")
    if len(got):
        by = got[~got["have"]].groupby("net").size().sort_values(ascending=False)
        print("  **まだ対が無い★を、探すネットワーク別に数えると**：")
        for net, k in by.items():
            print(f"    {net:<26} {k} 件")
    print("\n  === 取得の順序（提案）===")
    print("  1. **★かつ記録が長く、日本以外**——A-1 の同一地点検証を、**日本に固有でない形**にできる。")
    print("     現在の4組はすべて日本＝**擬似反復の懸念がある**（旗43 が名指しした穴と同型）。")
    print("  2. **★かつ記録が長い**（地域は問わない）——判定が安定する。")
    print("     旗74/75/77 では**日数不足で判定不能**になったサイトが複数あった。")
    print("  3. 乾燥サイト（US-Wkg / US-Whs / US-SRM）——**優先度は下げる**。")
    print("     旗70/71 で水マスターの一般化そのものが不確かになったため、")
    print("     **同じ枠組みでの追試より、レジームを分ける要因の解明が先**。")
    print("     ただし取得の手間は同じなので、ついでに取る価値はある。")
    print("\n  留保：")
    print("   ・**ネットワークは座標からの推定**であり、**そのサイトにタワーが在るとは限らない**。")
    print("     チャンバーの近くにタワーが無い場所は当然ある＝**取得先で実在を確認すること**。")
    print("   ・COSORE の座標そのものが正しい保証はない（旗64 で自分の座標抽出が壊れていた）。")
    print("     取得後は**必ず距離を計算して同一地点を確認**する（旗51 と同じ 10km 基準）。")
    print("   ・BASE 形式なら `sites.py` は対応済み＝**登録だけで既存の旗が走る**。")
    print("   ・**データ利用条件**：COSORE は Bond-Lamberty et al. 2020 GCB の引用と")
    print("     各データセット提供者へのクレジットが要る。取得先の網にも各々の規約がある。")


if __name__ == "__main__":
    main()
