"""旗69：**水マスターの独立地点 n=1 を、手元のデータだけで解消できるか**を確かめる。

主張 A-3（θ が蒸発を支配し Bowen 反転が起きる）は、旗62 で**季節を越えて成立**したが、
旗43 の指摘——**独立地点 n=1**——は残ったままである。使ったのは
**MN-Hst / MN-Nkh / MN-Kbu の3サイト**だが、これらは 50km 単連結で**1クラスタ**に潰れる。

`DRY_SITE_EXPANSION.md` は AmeriFlux（US-Wkg/US-Whs/US-SRM）の取得を計画したが、
**外部データを取りに行く前に手元を数え直す**のが順序である（旗67 と同じ理由）。
旗67 の棚卸しで **MN-Skt / MN-Udg** が手元にあると分かった。**別クラスタなら n が増える**。

本ツールは**事前登録の前**に、二つの前提条件だけを確かめる（**検定はしない**）：

  1. **独立性**：MN-Skt / MN-Udg は既存3サイトから **50km 超**か（＝別クラスタか）。
     旗64 の教訓に従い、**座標が引けなければ「引けない」と出す**（他サイトの値を返さない）。
  2. **実行可能性**：旗36 が要る **th, Rg, gLE, gH, Ta** が実ファイルにあるか。
     **本体HH に無ければ補助ファイルも探す**（旗68 の教訓）。

両方を満たしたサイトだけが、事前登録して検定する価値を持つ。
**満たさなければ「手元では解消できない」と確定し、外部データ取得の判断材料になる**。

    python research/dry_site_feasibility_step69.py --data-dir /mnt/hdd/JAPANFLUX
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from sif_coords import coords_from_badm
from colocate_step51 import haversine

EXISTING = ["MN-Hst", "MN-Nkh", "MN-Kbu"]      # 旗31/36/62 で使った3サイト
CANDIDATES = ["MN-Skt", "MN-Udg"]              # 旗67 の棚卸しで見つかった追加
CLUSTER_KM = 50.0                              # 旗43 と同じ単連結の閾値

# 旗36 が要る変数と、その既定列名（japanflux 形式）＋候補を探すトークン
NEED = {"th": ("SWC_F_MDS", ("SWC", "SM_", "SOIL_W", "VWC")),
        "Rg": ("SW_IN_F", ("SW_IN", "RG")),
        "gLE": ("LE_F_MDS", ("LE_",)),
        "gH": ("H_F_MDS", ("H_",)),
        "Ta": ("TA_F", ("TA_", "TAIR"))}


def get_coords(data_dir, codes):
    out, failed = {}, []
    for c in codes:
        try:
            lat, lon = coords_from_badm(data_dir, c)
        except Exception:
            lat = lon = None
        if lat is None or lon is None or not (np.isfinite(lat) and np.isfinite(lon)):
            failed.append(c); continue
        out[c] = (float(lat), float(lon))
    return out, failed


def check_vars(code, data_dir):
    """旗36 が要る5変数が実ファイルにあるか。**本体に無ければ補助ファイルも探す**。"""
    import pandas as pd
    root = Path(data_dir)
    csvs = sorted([p for p in root.rglob(f"*{code}*")
                   if p.is_file() and p.suffix.lower() == ".csv"
                   and "__MACOSX" not in p.parts and not p.name.startswith(("._", "~$"))],
                  key=lambda p: p.stat().st_size, reverse=True)
    if not csvs:
        return None, "このコードの csv が無い"
    status = {}
    scanned = []
    for f in csvs[:12]:
        try:
            cols = list(pd.read_csv(f, nrows=2).columns)
        except Exception:
            continue
        scanned.append(f.name)
        for k, (exact, toks) in NEED.items():
            if status.get(k, {}).get("exact"):
                continue
            if exact in cols:
                status[k] = {"exact": True, "col": exact, "file": f.name}
            elif k not in status:
                cand = [c for c in cols if any(t in c.upper() for t in toks)]
                if cand:
                    status[k] = {"exact": False, "col": cand[:4], "file": f.name}
    for k in NEED:
        status.setdefault(k, None)
    return status, f"{len(scanned)} ファイルを走査"


def main():
    p = argparse.ArgumentParser(description="乾燥サイト拡張の実行可能性（検定はしない）")
    p.add_argument("--data-dir", default="/mnt/hdd/JAPANFLUX")
    a = p.parse_args()

    print("=== 旗69：水マスターの独立地点 n=1 を手元で解消できるか（事前確認・検定はしない）===")
    print(f"  既存の3サイト {EXISTING} は 50km 単連結で**1クラスタ**＝旗43 の n=1。")
    print(f"  候補 {CANDIDATES} が**別クラスタ**かつ**変数が揃う**なら、事前登録して検定する価値がある。\n")

    codes = EXISTING + CANDIDATES
    coords, failed = get_coords(a.data_dir, codes)
    print("  ── ① 座標と距離（独立性） ──")
    for c in codes:
        if c in coords:
            print(f"    {c:<8} {coords[c][0]:>8.4f}, {coords[c][1]:>9.4f}")
        else:
            print(f"    {c:<8} **座標を引けない**")
    if failed:
        print(f"    引けなかった：{failed}（旗64 の修正により、"
              f"**他サイトの座標を代わりに返すことはしない**）")

    print(f"\n    既存クラスタからの距離（最小値・{CLUSTER_KM:.0f}km 超なら別クラスタ）：")
    verdict_indep = {}
    for c in CANDIDATES:
        if c not in coords:
            verdict_indep[c] = None
            print(f"    {c:<8} 座標が無く判定不能"); continue
        ds = [(haversine(*coords[c], *coords[e]), e) for e in EXISTING if e in coords]
        if not ds:
            verdict_indep[c] = None
            print(f"    {c:<8} 既存側の座標が無く判定不能"); continue
        dmin, near = min(ds)
        indep = dmin > CLUSTER_KM
        verdict_indep[c] = indep
        print(f"    {c:<8} 最近接 {near} まで {dmin:>8.1f} km"
              f"  → {'**別クラスタ＝独立地点が増える**' if indep else '同一クラスタ＝n は増えない'}")
    # 候補どうしも見る（2つが互いに近ければ、増えるのは1地点）
    if all(c in coords for c in CANDIDATES):
        d2 = haversine(*coords[CANDIDATES[0]], *coords[CANDIDATES[1]])
        print(f"    {CANDIDATES[0]} ↔ {CANDIDATES[1]}：{d2:.1f} km"
              f"  → {'互いも別クラスタ' if d2 > CLUSTER_KM else '**互いは同一クラスタ＝増えるのは1地点**'}")

    print("\n  ── ② 変数の有無（実行可能性） ──")
    print("    旗36 が要るのは th, Rg, gLE, gH, Ta。**本体HH に無ければ補助ファイルも探す**。")
    verdict_vars = {}
    for c in CANDIDATES:
        st, note = check_vars(c, a.data_dir)
        print(f"    ━ {c} ━（{note}）")
        if st is None:
            verdict_vars[c] = False
            print("      **ファイルが無い**"); continue
        ok = True
        for k, (exact, _) in NEED.items():
            s = st[k]
            if s is None:
                ok = False
                print(f"      {k:<4}: **どのファイルにも無い**")
            elif s["exact"]:
                print(f"      {k:<4}: {s['col']}（既定名で一致）")
            else:
                ok = False        # 候補はあるが**既定名でない**＝人が確認して決める
                print(f"      {k:<4}: 既定名 {exact} は無し。候補 {s['col']} @ {s['file']}"
                      f"  ← **人が確認して var_overrides を決める**")
        verdict_vars[c] = ok

    print("\n  === まとめ（この結果で次を決める）===")
    usable = [c for c in CANDIDATES
              if verdict_indep.get(c) and verdict_vars.get(c)]
    partial = [c for c in CANDIDATES
               if verdict_indep.get(c) and verdict_vars.get(c) is False]
    for c in CANDIDATES:
        i, v = verdict_indep.get(c), verdict_vars.get(c)
        tag = ("**そのまま使える**" if (i and v) else
               "独立だが変数の確認が要る" if (i and v is False) else
               "同一クラスタ＝n は増えない" if i is False else "判定不能")
        print(f"    {c:<8} 独立={i}  変数={v}  → {tag}")
    if usable:
        print(f"\n  → **{usable} で独立地点が増やせる**。次は**事前登録**（旗62 と同じ形式で、"
              f"\n     予測・統計量・有意性・最低検出力・限界を先に書いて commit）してから検定する。")
    elif partial:
        print(f"\n  → **{partial} は独立だが列名が既定と違う**。**推測で書かず**、"
              f"\n     上の候補列を見て `var_overrides` を決めてから登録する（旗68 と同じ作法）。")
    else:
        print("\n  → **手元では独立地点を増やせない**と確定する。"
              "\n     ＝`DRY_SITE_EXPANSION.md` の外部データ取得（AmeriFlux）が"
              "\n     **唯一の道であることが、推測でなく確認として言える**。")
    print("\n  留保：50km 単連結は旗43 で採った便宜的な閾値であり、**生態学的な独立性の保証ではない**。")
    print("        同じモンゴル半乾燥ステップである以上、距離が離れても**気候レジームは共有**しうる。")
    print("        ＝独立地点が増えても『別の気候帯で再現した』ことにはならない（そこは AmeriFlux が要る）。")


if __name__ == "__main__":
    main()
