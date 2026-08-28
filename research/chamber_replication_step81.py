"""旗81：**チャンバー1本の代表性**を測る——そして「平均が ACF1 を上げていないか」を確かめる。

旗66 拡張版で分かったこと：**US-SSH の 0.1km 以内 8 本のチャンバーが
★6本／メモリ弱1本／季節1本に割れた**。＝**「メモリが在る/無い」を 1 本で判定してきた前提が揺らぐ**。

さらに、道具を読み直して気づいたことがある——
**`load_cosore` はチャンバー個体の識別子を一切見ず、全行をそのまま日平均している**。
COSORE の 1 データセットが複数チャンバーを含むなら、**気づかないまま平均していた**ことになる。

**これは中立ではない**：
**独立なチャンバーを N 本平均すると、白色ノイズは 1/√N に減り、共通信号は残る**
＝**残差の自己相関（ACF1）は機械的に上がる**。
＝**チャンバー数の多いデータセットほど★が出やすい**という交絡がありうる。
**A-1 の件数（外挿 15/44）に直接効く**ので、確かめる。

## 本ツールがやること

  1. **識別子の有無を調べる**——`CSR_PORT` など、チャンバーを区別する列があるか。
     **無ければ「1 データセット＝1 系列」として扱ってきたのは正しかった**と分かる。
  2. **チャンバー数と ACF1 の関係**——数が多いほど ACF1 が高い、という関係があるか。
     **あれば、★の件数は「メモリの普及率」ではなく「平均本数の分布」を映している**ことになる。
  3. **1 本ずつの判定のばらつき**——3 本以上あるデータセットで、
     **チャンバーごとに判定し、多数派と一致する確率**を出す。
     ＝**「1 本で判定したときの当たる確率」**そのもの。

**検定ではなく測定**である。相関が出ても因果とは限らない（本数の多い研究は他の点でも違いうる）。

    python research/chamber_replication_step81.py --cosore-dir /mnt/hdd/cosore-0.7.0 --igbp forest
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import _pick_soil_temp, _pick_sm
from model_richness_step74 import measure, star

# チャンバー個体を区別しうる列（COSORE の規約と、実データで見かけうる別名）
ID_COLS = ("CSR_PORT", "CSR_CHAMBER_ID", "CSR_CHAMBER", "CSR_PLOT", "CSR_COLLAR",
           "CSR_TREATMENT", "CSR_REP")


def load_with_id(path, months=None):
    """`load_cosore` と同じ読み方に、**チャンバー識別子の列を足す**。"""
    df = pd.read_csv(path, low_memory=False)
    cols = list(df.columns)
    if "CSR_FLUX_CO2" not in cols:
        return None, None, None
    tcol = "CSR_TIMESTAMP_BEGIN" if "CSR_TIMESTAMP_BEGIN" in cols else "CSR_TIMESTAMP_END"
    ts = pd.to_datetime(df[tcol], errors="coerce")
    out = pd.DataFrame({"Rs": pd.to_numeric(df["CSR_FLUX_CO2"], errors="coerce").to_numpy()},
                       index=ts)
    st, sm = _pick_soil_temp(cols), _pick_sm(cols)
    if st:
        out["Tsoil"] = pd.to_numeric(df[st], errors="coerce").to_numpy()
    if sm:
        out["SM"] = pd.to_numeric(df[sm], errors="coerce").to_numpy()
    idc = next((c for c in ID_COLS if c in cols), None)
    if idc:
        out["_id"] = df[idc].astype(str).to_numpy()
    out = out[out.index.notna()]
    if months:
        out = out[out.index.month.isin(months)]
    return out.dropna(subset=["Rs"]), idc, (st, sm)


def daily_measure(sub):
    """日平均にして旗74 の外挿基準で測る。"""
    if "Tsoil" not in sub:
        return None
    cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in sub else [])
    daily = sub[cols].groupby(sub.index.normalize()).mean()
    if len(daily) < 60:
        return None
    daily = daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D"))
    y = daily["Rs"].to_numpy(); T = daily["Tsoil"].to_numpy()
    W = daily["SM"].to_numpy() if "SM" in daily else None
    return measure(y, T, W, "テンソルビン", True)


def main():
    p = argparse.ArgumentParser(description="チャンバー1本の代表性を測る")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--igbp", default="forest")
    p.add_argument("--month", type=int, nargs="+", default=None)
    p.add_argument("--min-ch", type=int, default=3, help="1本ずつ判定する最小チャンバー数")
    a = p.parse_args()
    root = Path(a.cosore_dir)
    desc = pd.read_csv(root / "description.csv")

    print("=== 旗81：チャンバー1本の代表性と、平均が ACF1 を上げていないか ===")
    print("  `load_cosore` は**識別子を見ずに全行を日平均**している。")
    print("  **N 本平均するとノイズは 1/√N・共通信号は残る＝ACF1 は機械的に上がる**。")
    print("  ＝**チャンバー数の多いデータセットほど★が出やすい**という交絡がありうる。\n")

    rows = []
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        if a.igbp and a.igbp.lower() not in ig.lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            df, idc, soil = load_with_id(f, a.month)
        except Exception:
            continue
        if df is None or "Tsoil" not in df:
            continue
        nch = int(df["_id"].nunique()) if "_id" in df else 1
        m = daily_measure(df)
        s = star(m)
        rows.append({"ds": ds, "idc": idc or "—", "nch": nch,
                     "acf1": (m or {}).get("acf1", np.nan),
                     "r2": (m or {}).get("r2", np.nan), "star": s,
                     "n": int(df["Rs"].notna().sum()), "df": f})
    if not rows:
        print("  対象が無い"); return
    t = pd.DataFrame(rows)

    # ① 識別子の有無
    print("  ── ① チャンバー識別子はあるか ──")
    print(f"    識別子の列がある：{int((t['idc'] != '—').sum())}/{len(t)} データセット")
    for idc, k in t["idc"].value_counts().items():
        print(f"      {idc:<20} {k} 件")
    multi = t[t["nch"] >= 2]
    print(f"    **2 本以上を含むデータセット：{len(multi)}/{len(t)}**"
          f"（最大 {int(t['nch'].max())} 本）")
    if len(multi) == 0:
        print("    → **1 データセット＝1 系列**だった＝これまでの扱いは正しい。")

    # ② 本数と ACF1 の関係
    print("\n  ── ② チャンバー数と ACF1 の関係（**平均がACF1を上げていないか**）──")
    ok = t[np.isfinite(t["acf1"])]
    if len(ok) >= 8 and ok["nch"].nunique() >= 2:
        from scipy.stats import spearmanr
        r, pv = spearmanr(ok["nch"], ok["acf1"])
        print(f"    Spearman r(本数, ACF1) = {r:+.3f}（p={pv:.3f}, n={len(ok)}）")
        for lo, hi in ((1, 1), (2, 4), (5, 9), (10, 10 ** 6)):
            g = ok[(ok["nch"] >= lo) & (ok["nch"] <= hi)]
            if len(g):
                sr = g["star"].map({True: 1, False: 0}).mean()
                print(f"    本数 {lo}–{hi if hi < 1000 else '∞'}：n={len(g):>3}"
                      f"  ACF1 中央 {g['acf1'].median():+.2f}"
                      f"  ★率 {sr:.0%}" if np.isfinite(sr) else "")
        print("    → **正の相関が強ければ、★の件数は「平均本数の分布」を映している**疑いが立つ。")
        print("      **相関が無ければ、この交絡は効いていない**＝15/44 はそのまま読める。")
    else:
        print("    本数に幅が無く、関係を測れない")

    # ③ 1 本ずつの判定のばらつき
    print(f"\n  ── ③ {a.min_ch} 本以上あるデータセットで、**1 本ずつ判定する** ──")
    tgt = t[t["nch"] >= a.min_ch]
    if len(tgt) == 0:
        print(f"    {a.min_ch} 本以上を含むデータセットが無い＝この測定はできない")
    else:
        agree_all, n_all = 0, 0
        for _, r in tgt.iterrows():
            df, idc, soil = load_with_id(r["df"], a.month)
            per = {}
            for cid, sub in df.groupby("_id"):
                m = daily_measure(sub)
                s = star(m)
                if s is not None:
                    per[cid] = (s, m["acf1"])
            if len(per) < 2:
                continue
            stars = sum(1 for s, _ in per.values() if s)
            maj = stars * 2 >= len(per)
            agr = sum(1 for s, _ in per.values() if s == maj)
            agree_all += agr; n_all += len(per)
            acfs = [a1 for _, a1 in per.values()]
            print(f"    {r['ds']:<30} {len(per)} 本："
                  f"★{stars}／·{len(per)-stars}"
                  f"  ACF1 {min(acfs):+.2f}〜{max(acfs):+.2f}（幅 {max(acfs)-min(acfs):.2f}）"
                  f"  プール判定 {'★' if r['star'] else '·'}")
        if n_all:
            print(f"\n    **1 本を無作為に採ったとき多数派と一致する確率：{agree_all/n_all:.0%}**"
                  f"（{agree_all}/{n_all} 本）")
            print("    ＝**1 本での判定がどれだけ当てにならないか**の直接の測定。")

    print("\n  === 読み方 ===")
    print("  ・②で**正の相関**が出れば、**★の件数は本数の交絡を含む**＝A-1 の 15/44 に注釈が要る。")
    print("  ・③の一致率が**低い**ほど、**1 本で「メモリが在る」と言うことの不確かさが大きい**。")
    print("  ・**プール判定と多数派がずれる**データセットがあれば、")
    print("    **平均が個々のチャンバーには無い性質を作っている**可能性を示す。")
    print("  留保：")
    print("   ・**検定ではなく測定**。本数と ACF1 の相関が出ても因果とは限らない")
    print("     （本数の多い研究は、期間・機材・生態系でも違いうる）。")
    print("   ・1 本ずつにすると**日数が減る**ので、判定不能が増える方向に偏る。")
    print("   ・識別子が**処理区**（例 CSR_TREATMENT）の場合、本数の違いは")
    print("     **空間のばらつきではなく実験処理の違い**である＝そう読むこと。")


if __name__ == "__main__":
    main()
