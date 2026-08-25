"""旗53：旗52 の発見を受けて、チャンバーの ~4日メモリ（旗40/45）を較正し直す。

旗52 で判明したこと：**呼吸は駆動の非線形関数なのに日平均へ線形回帰しているため、
未観測駆動がゼロでも自己相関残差が出る**（帰無で線形 e-fold=6.2日 → 非線形基底で 2.6日）。
さらに **旗40 の判定基準（R²≥0.3 × ACF1>0.4 × e-fold 2〜7日）は帰無データを通してしまう**。
＝旗40 の「★短メモリ 12/35」には無視できない偽陽性率がある。

本ツールはチャンバー（＝分割も穴埋めも通さない直接測定）専用に較正し直す：
  1. **合成チャンバー**で、仕込みの強さ hid_sd を 0（帰無）から変えながら、
     線形／非線形基底の検出器が何を出すかを測る＝**帰無分布と検出力曲線**を作る。
  2. その較正から**判別統計と閾値**を決める。較正の結論（実行済み）：
     ・**e-fold は判別に使えない**（帰無 2.5日 vs 4日を仕込んで 3.7日＝ほとんど重なる）。
     ・**非線形基底の ACF1 が分離する**（帰無 0.49 vs 仕込み 0.80、仕込み強度 0.01〜0.08 で不変）。
     ・**旗40 の「ACF1>0.4」は帰無(0.49)をそのまま通す**＝閾値が低すぎた。較正値は **0.64**。
  3. 実データ（COSORE）を**非線形基底＋ACF1≥0.64**で再集計し、旗40 の集計と並べる。

    python research/chamber_memory_recount_step53.py                                   # 較正のみ
    python research/chamber_memory_recount_step53.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import _acf_gap, _efold_gap, load_cosore
from synthetic_tower_step52 import simulate_truth, T0_LT


def detect(daily, col, Tcol, Wcol, form="flex"):
    """日次系列の残差メモリ。form='flex' で非線形基底（Lloyd-Taylor項・二次・交互作用）を使う。"""
    y = np.log(daily[col].where(daily[col] > 0)).to_numpy()
    T = daily[Tcol].to_numpy(); W = daily[Wcol].to_numpy() if Wcol else None
    cols = [T] + ([W] if W is not None else [])
    if form == "flex":
        with np.errstate(divide="ignore", invalid="ignore"):
            lt = 320.0 * (1.0 / (10 - T0_LT) - 1.0 / (T - T0_LT))
        cols += [lt, T ** 2]
        if W is not None:
            cols += [W ** 2, T * W, np.log(np.clip(W, 1e-3, None))]
    X = np.column_stack(cols + [np.ones(len(daily))])
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    if ok.sum() < 40:
        return None
    coef = np.linalg.lstsq(X[ok], y[ok], rcond=None)[0]
    res = np.full(len(y), np.nan); res[ok] = y[ok] - X[ok] @ coef
    ss = np.sum((y[ok] - y[ok].mean()) ** 2)
    return {"r2": float(1 - np.sum(res[ok] ** 2) / ss) if ss > 0 else np.nan,
            "acf1": _acf_gap(res, 1), "efold": _efold_gap(res), "n": int(ok.sum())}


def _daily_from_truth(df):
    d = df[["RECO_true", "Ts", "th"]].copy()
    daily = d.groupby(d.index.normalize()).mean()
    return daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D"))


def calibrate(nrep=6, years=3, sds=(0.0, 0.01, 0.02, 0.04, 0.08), efold_days=4.0):
    """合成チャンバー（分割も穴埋めも通さない）で帰無分布と検出力曲線を作る。"""
    print("=== 旗53 較正：合成チャンバー（直接測定＝処理を通さない）===")
    print("  仕込み強度を変えながら、線形／非線形基底の検出器が何を出すかを見る。\n")
    print(f"  {'仕込み hid_sd':>13}{'線形 ACF1':>11}{'線形 e-fold':>12}"
          f"{'柔軟 ACF1':>11}{'柔軟 e-fold':>12}{'柔軟 R²':>9}")
    out = {}
    for sd in sds:
        rows = {"linear": [], "flex": []}
        for k in range(nrep):
            tru = simulate_truth(years=years, plant_memory=(sd > 0), seed=k,
                                 hid_sd=max(sd, 1e-9), hid_efold_days=efold_days)
            daily = _daily_from_truth(tru)
            for f in ("linear", "flex"):
                r = detect(daily, "RECO_true", "Ts", "th", f)
                if r:
                    rows[f].append(r)
        m = {f: {k2: float(np.nanmean([x[k2] for x in rows[f]])) for k2 in ("r2", "acf1", "efold")}
             for f in rows}
        out[sd] = m
        tag = "（帰無）" if sd == 0 else ""
        print(f"  {sd:>10.3f}{tag:>3}{m['linear']['acf1']:>11.2f}{m['linear']['efold']:>12.1f}"
              f"{m['flex']['acf1']:>11.2f}{m['flex']['efold']:>12.1f}{m['flex']['r2']:>9.2f}")
    null = out[0.0]
    pl = [out[sd]["flex"]["acf1"] for sd in sds if sd > 0]
    print(f"\n  **帰無**：線形 ACF1={null['linear']['acf1']:.2f}/e-fold={null['linear']['efold']:.1f}日、"
          f"柔軟 ACF1={null['flex']['acf1']:.2f}/e-fold={null['flex']['efold']:.1f}日")
    print(f"  **仕込み**（柔軟）：ACF1 {min(pl):.2f}〜{max(pl):.2f}")
    print("  → **e-fold は判別に使えない**：帰無 {:.1f}日 vs 仕込み {:.1f}日 でほとんど重なる。".format(
        null["flex"]["efold"], np.mean([out[sd]["flex"]["efold"] for sd in sds if sd > 0])))
    print("    旗40 が e-fold を主判定に使っていたのは**不適切**だった（帰無をそのまま通す）。")
    thr = float((null["flex"]["acf1"] + np.mean(pl)) / 2)
    print(f"  → 判別統計を **非線形基底の ACF1** に変え、閾値は帰無({null['flex']['acf1']:.2f})と"
          f"仕込み({np.mean(pl):.2f})の中点 **{thr:.2f}** を採る。")
    return thr, null


def run_real(cosore_dir, igbp, months, thr):
    root = Path(cosore_dir)
    desc = pd.read_csv(root / "description.csv")
    print(f"\n=== 旗53 実データ：チャンバーRs を非線形基底で再集計（{igbp or '全'}）===")
    print(f"  {'dataset':<32}{'線形 ACF1':>10}{'線形ef':>8}{'柔軟 ACF1':>10}{'柔軟ef':>8}"
          f"{'柔軟R²':>8}  旗40基準 / 較正後")
    tally = {"旗40★": 0, "較正後★": 0, "判定可能": 0}
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        if igbp and igbp.lower() not in ig.lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            df, st, sm = load_cosore(f, months)
            if "Tsoil" not in df:
                continue
            cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
            dd = df[cols].copy()
            daily = dd.groupby(dd.index.normalize()).mean()
            daily = daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D"))
            lin = detect(daily, "Rs", "Tsoil", "SM" if "SM" in df else None, "linear")
            flx = detect(daily, "Rs", "Tsoil", "SM" if "SM" in df else None, "flex")
        except Exception:
            continue
        if not lin or not flx or not np.isfinite(lin["efold"]) or not np.isfinite(flx["efold"]):
            continue
        tally["判定可能"] += 1
        old = (lin["r2"] >= 0.3) and (lin["acf1"] > 0.4) and (2 <= lin["efold"] <= 7)   # 旗40 の基準
        new = (flx["r2"] >= 0.3) and (flx["acf1"] >= thr)                               # 較正後
        tally["旗40★"] += int(old); tally["較正後★"] += int(new)
        print(f"  {ds:<32}{lin['acf1']:>10.2f}{lin['efold']:>8.1f}{flx['acf1']:>10.2f}"
              f"{flx['efold']:>8.1f}{flx['r2']:>8.2f}  {'★' if old else '·'} / {'★' if new else '·'}")
    print(f"\n  === まとめ（判定可能 {tally['判定可能']}）===")
    print(f"  旗40 の基準（線形・e-fold 2〜7日）で ★：{tally['旗40★']}")
    print(f"  較正後（非線形基底・ACF1 ≥ {thr:.2f}）で ★：{tally['較正後★']}")
    print("  読み：★が大きく減る＝旗40 の集計には**検出器由来の偽陽性**が含まれていた。")
    print("        減らない＝メモリは非線形の取り残しでは説明できない＝旗40 の主張は較正後も立つ。")
    print("  留保：較正は合成に依存する。閾値は『帰無＋余裕』であって物理的な境界ではない。")


def main():
    p = argparse.ArgumentParser(description="チャンバーのメモリを非線形基底で再集計")
    p.add_argument("--cosore-dir"); p.add_argument("--igbp", default="forest")
    p.add_argument("--month", type=int, nargs="+", default=None)
    p.add_argument("--nrep", type=int, default=6)
    a = p.parse_args()
    thr, _ = calibrate(nrep=a.nrep)
    if a.cosore_dir:
        run_real(a.cosore_dir, a.igbp, a.month, thr)


if __name__ == "__main__":
    main()
