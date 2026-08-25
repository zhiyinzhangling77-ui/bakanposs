"""旗61：「~4日メモリは森林でなく木本性に伴うか」——**探索的**分析（確認的検定ではない）。

**登録は `research/PREREGISTRATION_step61.md`（72f0afb, 本実行の前に commit）。**
そこに書いた通り、**この問いに対する未使用データは COSORE にもう無い**（森林は旗53/54/58/59、
非森林は旗57/59 で使い、サイト別 ACF1 まで出力で見ている）。
＝**いま同じデータで計算しても確認的検定にはならない**。本ファイルが出すのは
**「登録する価値があるかの下見」＝探索的な数値**であり、証拠として扱わない。

登録した手続き（将来の独立データにそのまま適用する）：
  木本性階級 2=高木林／1=低木・疎林・サバンナ／0=草本・農地・湿地（IGBP から機械的に決める）
  統計量：階級 と **連続量 ACF1**（二値化しない）の Spearman
  検定：階級ラベルの順列（両側・10000回）／主判定は **50km クラスタ縮約**／クラスタ15未満なら判定しない

    python research/woodiness_step61.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore
from chamber_memory_recount_step53 import detect
from out_of_sample_step48 import geo_clusters
from memory_attribution_flex_step54 import ACF_THR, EFOLD_MAX

MIN_CLUSTERS = 15          # 登録：クラスタがこれ未満なら判定しない


def woodiness(igbp: str) -> int | None:
    """IGBP から木本性階級を機械的に決める（判断の余地を残さない）。"""
    g = igbp.lower()
    if "forest" in g or "needleleaf" in g or "broadleaf" in g:
        return 2
    if any(k in g for k in ("shrub", "savanna", "woodland")):
        return 1
    if any(k in g for k in ("grass", "crop", "wetland", "barren", "tundra")):
        return 0
    return None


def _spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 4:
        return np.nan
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    if rx.std() == 0 or ry.std() == 0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def perm_p(cls, val, nperm=10000, seed=0):
    r = _spearman(cls, val)
    if not np.isfinite(r):
        return r, np.nan
    rng = np.random.default_rng(seed)
    null = np.array([_spearman(rng.permutation(cls), val) for _ in range(nperm)])
    null = null[np.isfinite(null)]
    return r, float((np.abs(null) >= abs(r)).mean()) if null.size else np.nan


def main():
    p = argparse.ArgumentParser(description="木本性とメモリ（探索的）")
    p.add_argument("--cosore-dir", required=True)
    a = p.parse_args()
    root = Path(a.cosore_dir); desc = pd.read_csv(root / "description.csv")

    print("=== 旗61：木本性と ~4日メモリ（**探索的**・登録 72f0afb）===")
    print("  ※この問いに対する未使用データは無い。以下は証拠ではなく『下見』である。\n")
    rows = []
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        w = woodiness(ig)
        if w is None:
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
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
        if not flx or not np.isfinite(flx["r2"]) or flx["r2"] < 0.3 or not np.isfinite(flx["acf1"]):
            continue
        try:
            lat, lon = float(d["CSR_LATITUDE"]), float(d["CSR_LONGITUDE"])
        except (TypeError, ValueError, KeyError):
            lat = lon = np.nan
        rows.append({"ds": ds, "igbp": ig, "w": w, "acf1": flx["acf1"],
                     "efold": flx["efold"], "lat": lat, "lon": lon,
                     "short": int(flx["acf1"] >= ACF_THR and flx["efold"] <= EFOLD_MAX)})
    if len(rows) < 10:
        print(f"  判定可能 {len(rows)} 件＝少なすぎる。"); return
    df = pd.DataFrame(rows)

    print(f"  {'木本性':<22}{'n':>4}{'ACF1 中央値':>12}{'短メモリ率':>11}")
    for w, lab in [(2, "2 高木林"), (1, "1 低木・疎林"), (0, "0 草本・農地・湿地")]:
        g = df[df["w"] == w]
        if g.empty:
            continue
        print(f"  {lab:<22}{len(g):>4}{g['acf1'].median():>12.2f}{g['short'].mean():>11.2f}")

    r1, p1 = perm_p(df["w"].to_numpy(), df["acf1"].to_numpy())
    print(f"\n  --- サイト単位（参考）--- n={len(df)}  Spearman r={r1:+.2f}  両側 p={p1:.4f}")

    ok = np.isfinite(df["lat"]) & np.isfinite(df["lon"])
    if ok.all():
        g = geo_clusters(df["lat"].to_numpy(), df["lon"].to_numpy(), km=50.0)
        cl = df.assign(cl=g).groupby("cl").agg(w=("w", "mean"), acf1=("acf1", "mean"))
        r2, p2 = perm_p(cl["w"].to_numpy(), cl["acf1"].to_numpy())
        print(f"  --- **クラスタ単位（主判定）** --- {len(cl)} クラスタ  "
              f"Spearman r={r2:+.2f}  両側 p={p2:.4f}")
        if len(cl) < MIN_CLUSTERS:
            print(f"  → 登録の通り**判定しない**（クラスタ {len(cl)} < {MIN_CLUSTERS}）。")
            r2 = p2 = np.nan
    else:
        r2 = p2 = np.nan
        print("  --- 座標欠損でクラスタ縮約できず ---")

    print("\n  === 下見としての読み（**証拠ではない**）===")
    if np.isfinite(p2) and p2 < 0.05 and r2 > 0:
        print("  → 木本性と ACF1 に正の関連が**見える**。だが同じデータから仮説を作っているので、")
        print("     これは**登録する価値がある**という以上のことを意味しない。将来の独立データで検定する。")
    elif np.isfinite(p2):
        print("  → 木本性との関連は**見えない**。将来データで検定する価値は低いと判断してよい。")
    print("  ＝いずれにせよ現在のスコープ（旗57：普及率は『森林で』）は変えない。")
    print("  留保：IGBP は木本性の粗い代理。木本性は根の深さ・リター量・撹乱頻度と交絡する。")


if __name__ == "__main__":
    main()
