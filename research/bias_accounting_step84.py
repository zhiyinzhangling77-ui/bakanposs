"""旗84：**選定バイアスの三つを数字にする**（旗83 で材料があると確認できた分）。

旗83 の結果：
  ・**AmeriFlux 全サイトの座標は手元に無い**（multi-site BADM 0 件）＝**総当たりはできない**。
    → 取得ページの「**BADM Only**」（837 サイト）を落とせば可能になる。
  ・**COSORE に `CSR_PRIMARY_PUB`（55 通り/88）がある**＝**研究群の代理として使える**。
  ・**ふるいの脱落を数える材料は揃っている**。

そこで、**いま測れる三つ**を数字にする。

## ① **★選択バイアスの大きさ**（私が持ち込んだもの）

旗78 で**★が出るチャンバー**を上位に並べ、**その上位を取得するよう勧めた**。
＝同一地点の対のうち、**取得の動機が★だったもの**と、**★でないまま付いてきたもの**がある。

  ・**★で選んだ**：US-SSH（KAYE）／US-Ha1（SAVAGE）／US-MMS（ZHANG）／CA-TP3（TP74）
  ・**付いてきた**：CA-TP4・CA-TPD（Turkey Point 一式）／乾燥地3（優先度を下げたが手間が同じで取得）
  ・**★以前に見つけた**：日本の4組（旗64/67 の座標探索で見つけ、★では選んでいない）

**両群で一致率を比べる**＝**選択の効果が数字になる**。**差が大きければ、9 組という数字は過大**である。

## ② **研究者単位のクラスタ縮約**

旗43 の縮約は**地理だけ**で、**同じ研究者・同じ機材・同じ処理**を独立と数えていた。
`CSR_PRIMARY_PUB` で括り直し、**A-1 の件数が「研究群のうち何群」になるか**を出す。
＝**旗82 の「独立クラスタ 1→3」も、研究群では何になるか**を併記する。

## ③ **ふるいの脱落**

`Tsoil あり → 日数≥60 → R²≥0.3 → 判定` の各段階で**何件落ちたか**を、**森林／非森林別**に数える。
`R²≥0.3` は「駆動でよく説明できる」ふるいなので、**メモリの出方と相関しうる**。

    python research/bias_accounting_step84.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore
from model_richness_step74 import measure, star
from same_site_arc_step66 import PAIRS, memory_from_daily, tower_daily, chamber_daily, verdict
from colocate_step51 import haversine

# **取得の動機**でタワーを分ける（旗78 の順位に従って選んだか、付いてきたか）
STAR_PICKED = {"US-SSH", "US-Ha1", "US-MMS", "CA-TP3"}
CAME_ALONG = {"CA-TP4", "CA-TPD", "US-Wkg", "US-Whs", "US-SRM"}
BEFORE_STAR = {"JP-Fhk", "JP-Tef", "JP-Yms"}      # 旗64/67 の座標探索で見つけた


def pair_agreement(cosore_dir):
    """各対で**同一期間**の判定を取り、両側が一致するかを返す。"""
    root = Path(cosore_dir)
    out = []
    for site, ds, km in PAIRS:
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            ch_all, _ = chamber_daily(f)
        except Exception:
            continue
        if ch_all is None:
            continue
        span = (ch_all.index.min(), ch_all.index.max())
        try:
            tw = tower_daily(site, span)
        except Exception:
            tw = None
        mc = memory_from_daily(ch_all, "Rs", "Tsoil", "SM")
        mt = memory_from_daily(tw, "GER", "Ts", "th") if tw is not None else None
        vc = verdict(mc); vt = verdict(mt) if mt is not None else "判定不能"
        grp = ("★で選んだ" if site in STAR_PICKED else
               "付いてきた" if site in CAME_ALONG else
               "★以前に発見" if site in BEFORE_STAR else "その他")
        judged = ("判定不能" not in vc and "判定不能" not in vt
                  and "推定不能" not in vc and "推定不能" not in vt)
        agree = judged and (("★" in vc) == ("★" in vt))
        out.append({"site": site, "ds": ds, "grp": grp, "vc": vc, "vt": vt,
                    "judged": judged, "agree": agree})
    return pd.DataFrame(out)


def main():
    p = argparse.ArgumentParser(description="選定バイアスを数字にする")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--igbp", default="forest")
    a = p.parse_args()
    root = Path(a.cosore_dir)
    desc = pd.read_csv(root / "description.csv")

    print("=== 旗84：選定バイアスの三つを数字にする ===")
    print("  旗83 で**全サイト座標は手元に無い**と確定＝**総当たりはできない**。")
    print("  代わりに**取得の動機で群を分けて比べる**（弱い形だが、方向と大きさは見える）。\n")

    # ── ① ★選択バイアス ──
    print("  ── ① ★選択バイアスの大きさ ──")
    pa = pair_agreement(a.cosore_dir)
    if pa.empty:
        print("    対が作れない")
    else:
        print(f"  {'タワー':<9}{'チャンバー':<32}{'群':<12}{'一致':>5}  チャンバー判定 / タワー判定")
        for _, r in pa.iterrows():
            mark = "○" if r["agree"] else ("×" if r["judged"] else "—")
            print(f"  {r['site']:<9}{r['ds']:<32}{r['grp']:<12}{mark:>5}  "
                  f"{r['vc'][:18]} / {r['vt'][:18]}")
        print(f"\n    {'群':<14}{'対の数':>6}{'判定できた':>10}{'一致':>6}{'一致率':>8}"
              f"{'タワー単位':>12}")
        for g, sub in pa.groupby("grp"):
            j = sub[sub["judged"]]
            rate = j["agree"].mean() if len(j) else np.nan
            tw = j.groupby("site")["agree"].mean()
            tw_rate = (tw >= 0.5).mean() if len(tw) else np.nan
            print(f"    {g:<14}{len(sub):>6}{len(j):>10}{int(j['agree'].sum()):>6}"
                  f"{(f'{rate:.0%}' if np.isfinite(rate) else '—'):>8}"
                  f"{(f'{tw_rate:.0%}（{len(tw)}本）' if np.isfinite(tw_rate) else '—'):>12}")
        print("    → **★で選んだ群の一致率が明確に高ければ、9 組という数字は過大**である。")
        print("      **タワー単位**も見ること（US-SSH が 8 組を占める＝擬似反復）。")

    # ── ②③ 全データセットの判定と、群・ふるい ──
    print("\n  ── ②③ 研究者単位の縮約と、ふるいの脱落 ──")
    stage = defaultdict(lambda: defaultdict(int))
    rows = []
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        # **欠損を1つの群に潰さない**（自分の道具の欠陥20件目）。
        # `str(NaN)` は "nan" という**中身のある文字列**になるため、
        # `or` のフォールバックが効かず、**出典不明の 12 件が1群に潰れていた**
        # ＝**研究群の数を過小評価**していた（＝擬似反復を実際より軽く見せる方向）。
        # 出典が無いものは**それぞれ別群**として扱う（**同じ研究の可能性は残る**＝そう明記）。
        _pub = str(d.get("CSR_PRIMARY_PUB", "")).strip()
        pub = _pub if _pub and _pub.lower() != "nan" else f"（出典不明:{ds}）"
        kind = "森林" if "forest" in ig.lower() or "plantation" in ig.lower() else "非森林"
        stage[kind]["① description に在る"] += 1
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        stage[kind]["② データファイルが在る"] += 1
        try:
            df, st, sm = load_cosore(f, None)
        except Exception:
            continue
        if df is None or "Tsoil" not in df or "Rs" not in df:
            continue
        stage[kind]["③ Rs と Tsoil が在る"] += 1
        cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
        daily = df[cols].groupby(df.index.normalize()).mean()
        if len(daily) < 60:
            continue
        stage[kind]["④ 日数 ≥ 60"] += 1
        daily = daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D"))
        m = measure(daily["Rs"].to_numpy(), daily["Tsoil"].to_numpy(),
                    daily["SM"].to_numpy() if "SM" in daily else None, "テンソルビン", True)
        if m is None or not np.isfinite(m.get("r2", np.nan)):
            continue
        stage[kind]["⑤ 当てはめができた"] += 1
        if m["r2"] < 0.3:
            continue
        stage[kind]["⑥ R² ≥ 0.3（判定可能）"] += 1
        s = star(m)
        if s:
            stage[kind]["⑦ ★短メモリ"] += 1
        try:
            la, lo = float(d["CSR_LATITUDE"]), float(d["CSR_LONGITUDE"])
        except (TypeError, ValueError, KeyError):
            la = lo = np.nan
        rows.append({"ds": ds, "kind": kind, "pub": pub, "star": bool(s),
                     "acf1": m["acf1"], "r2": m["r2"], "lat": la, "lon": lo})

    print(f"    {'段階':<24}{'森林':>7}{'非森林':>8}")
    keys = ["① description に在る", "② データファイルが在る", "③ Rs と Tsoil が在る",
            "④ 日数 ≥ 60", "⑤ 当てはめができた", "⑥ R² ≥ 0.3（判定可能）", "⑦ ★短メモリ"]
    for k in keys:
        print(f"    {k:<24}{stage['森林'][k]:>7}{stage['非森林'][k]:>8}")
    print("    → **どの段階で落ちているか**。⑥で大きく落ちるなら、")
    print("      **『駆動でよく説明できる』ふるいが母集団を作り替えている**。")

    t = pd.DataFrame(rows)
    if len(t):
        fo = t[t["kind"] == "森林"]
        print(f"\n    **研究者単位（CSR_PRIMARY_PUB）の縮約・森林**")
        g = fo.groupby("pub")["star"].agg(["size", "sum"])
        print(f"      データセット {len(fo)} 件／**研究群 {len(g)} 群**")
        n_unknown = int(fo["pub"].str.startswith("（出典不明").sum())
        print(f"      ★を1つ以上含む群：**{int((g['sum'] > 0).sum())}/{len(g)}**"
              f"（うち**出典不明で1件ずつ別群にした分 {n_unknown} 件**）")
        known = fo[~fo["pub"].str.startswith("（出典不明")]
        if len(known):
            gk = known.groupby("pub")["star"].agg(["size", "sum"])
            print(f"      **出典が分かる分だけ**：{len(known)} 件／{len(gk)} 群／"
                  f"★を含む群 {int((gk['sum'] > 0).sum())}/{len(gk)}")
        print(f"      （データセット単位では ★ {int(fo['star'].sum())}/{len(fo)}）")
        multi = g[g["size"] >= 3].sort_values("size", ascending=False)
        if len(multi):
            print(f"      **3 件以上を出している群**（＝擬似反復の源）：")
            for pub, r in multi.head(8).iterrows():
                print(f"        {int(r['size']):>2} 件（★{int(r['sum'])}）  {pub[:70]}")
        print("      → **群で数えると件数がどう変わるか**が、擬似反復の効き方そのもの。")

        # **論文ではなく「場所」で括る**（これが本来の独立単位）。
        # 論文単位は**出典が欠損すると 1 件ずつ別群になり、独立性を過大評価**する——
        # 実際 KAYE の 8 データセットは **1 つの観測所（Susquehanna Shale Hills CZO）**である。
        print(f"\n    **場所で括った縮約・森林**（座標の単連結。旗43 は 50km を使った）")
        fo2 = fo.dropna(subset=["lat", "lon"])
        for km in (1.0, 50.0):
            lab = np.arange(len(fo2))
            pts = fo2[["lat", "lon"]].to_numpy()
            for i in range(len(pts)):          # 単連結（総当たりで十分な規模）
                for j in range(i + 1, len(pts)):
                    if haversine(*pts[i], *pts[j]) <= km:
                        old, new = lab[j], lab[i]
                        if old != new:
                            lab[lab == old] = new
            gg = fo2.assign(cl=lab).groupby("cl")["star"].agg(["size", "sum"])
            print(f"      {km:>5.0f} km で括ると：**{len(gg)} 箇所**／"
                  f"★を含む箇所 **{int((gg['sum'] > 0).sum())}/{len(gg)}**"
                  f"（{int((gg['sum'] > 0).sum()) / len(gg):.0%}）")
            big = gg[gg["size"] >= 3].sort_values("size", ascending=False)
            if len(big):
                print(f"        3 件以上が同じ箇所：{list(big['size'].astype(int))}")
        print("      → **これが A-1 の本当の分母**：論文でも データセットでもなく、**場所の数**。")

    print("\n  === 読み方 ===")
    print("  ①の差が大きい＝**私の取得選択が結論を押し上げていた**。差が小さければ影響は軽い。")
    print("  ②で**群の数がデータセット数よりずっと少ない**＝**15/44 は独立な 44 ではない**。")
    print("  ③で⑥の脱落が大きい＝**判定できたサイトは「駆動が効く場所」に偏っている**。")
    print("  留保：")
    print("   ・`CSR_PRIMARY_PUB` は**論文**であり研究者そのものではない。")
    print("     **同じ人の別論文は別群になり、共著は1群に潰れる**＝**近似**である。")
    print("   ・①は**総当たりではない**（全サイト座標が無いため）。**方向と大きさの見当**であって、")
    print("     **選択の効果を厳密に推定したものではない**。")
    print("   ・タワー設置バイアス（平坦・均質・フェッチ良好）は**どの方法でも消えない**。")


if __name__ == "__main__":
    main()
