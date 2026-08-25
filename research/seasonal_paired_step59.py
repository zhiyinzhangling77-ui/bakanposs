"""旗59：季節でメモリが弱まるかを、二値化しない対比較で検定する（旗58 の設計ミスの修正）。

**手続き・判定規則は `research/PREREGISTRATION_step59.md` に実行前確定・commit 済み（0a2a563）。**

  対統計量：d_i = ACF1(休眠期)_i − ACF1(生育期)_i の平均（二値化しない）
  検定：符号反転による対順列検定（両側・10000回）＋ブートストラップ95%CI
  **森林＝探索的**（旗58 でサイト別数値を既に見たため確認的検定として扱わない）
  **非森林＝確認的**（季節の問いには未使用）。判定規則は非森林に適用する。
  対になるサイトが 8 未満なら判定しない。

    python research/seasonal_paired_step59.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore
from seasonal_memory_step58 import season_split, memory_of, MIN_R2, MIN_AMP

N_MIN_PAIR = 8          # 事前登録：対がこれ未満なら判定しない


def sign_flip_p(d, nperm=10000, seed=0):
    """符号反転による対順列検定（両側）。帰無：季節による系統的な差は無い。"""
    d = np.asarray(d, float)
    obs = float(d.mean())
    rng = np.random.default_rng(seed)
    flips = rng.choice([-1.0, 1.0], size=(nperm, len(d)))
    null = (flips * d).mean(axis=1)
    return obs, float((np.abs(null) >= abs(obs)).mean())


def boot_ci(d, nboot=2000, seed=0):
    d = np.asarray(d, float)
    rng = np.random.default_rng(seed)
    ms = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(nboot)]
    return float(np.percentile(ms, 2.5)), float(np.percentile(ms, 97.5))


def collect(root, desc, forest: bool):
    """森林/非森林それぞれについて、両季節そろったサイトの (生育期ACF1, 休眠期ACF1) を集める。"""
    rows = []
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        is_forest = "forest" in ig.lower()
        if is_forest != forest:
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
            grid = pd.date_range(daily.index.min(), daily.index.max(), freq="D")
            daily = daily.reindex(grid)
            amp = float(daily["Tsoil"].quantile(0.90) - daily["Tsoil"].quantile(0.10))
            if not np.isfinite(amp) or amp < MIN_AMP:
                continue
            g, h = season_split(daily.dropna(subset=["Tsoil"]))
            rg, rh = memory_of(g, grid), memory_of(h, grid)
        except Exception:
            continue
        if rg is None or rh is None or rg["r2"] < MIN_R2 or rh["r2"] < MIN_R2:
            continue
        if not (np.isfinite(rg["acf1"]) and np.isfinite(rh["acf1"])):
            continue
        rows.append({"ds": ds, "igbp": ig, "grow": rg["acf1"], "dorm": rh["acf1"],
                     "d": rh["acf1"] - rg["acf1"]})
    return pd.DataFrame(rows)


def report(df, label, confirmatory):
    print(f"\n  === {label}（{'**確認的**' if confirmatory else '探索的＝データ既見, 証拠として扱わない'}）===")
    if df.empty:
        print("    対になるサイトなし"); return
    print(f"  {'dataset':<32}{'生育期':>8}{'休眠期':>8}{'差 d':>8}")
    for _, r in df.sort_values("d").iterrows():
        print(f"  {r['ds']:<32}{r['grow']:>8.2f}{r['dorm']:>8.2f}{r['d']:>+8.2f}")
    n = len(df)
    print(f"\n  対サイト数 n={n}")
    if n < N_MIN_PAIR:
        print(f"  → **事前登録の通り判定しない**（n={n} < {N_MIN_PAIR}＝検出力不足）。")
        return
    obs, p = sign_flip_p(df["d"].to_numpy())
    lo, hi = boot_ci(df["d"].to_numpy())
    n_neg = int((df["d"] < 0).sum())
    print(f"  mean(d) = {obs:+.3f}  95%CI [{lo:+.3f}, {hi:+.3f}]  "
          f"符号反転検定 両側 p={p:.4f}  （d<0 のサイト {n_neg}/{n}）")
    if not confirmatory:
        print("  ※探索的：この p 値は事前登録された証拠ではない（同じデータで統計量を選び直しているため）。")
        return
    print("\n  === 事前登録の判定規則に照らす（確認的）===")
    if p < 0.05 and obs < 0:
        print("  → **メモリは休眠期に弱まる**＝生育期の**生物活動に結びつく**")
        print("     （遅い基質・微生物・根の動態）。物理拡散だけでは説明しにくい。")
    elif p < 0.05 and obs > 0:
        print("  → **休眠期にむしろ強まる**＝植物活動と逆向き＝物理過程・遅いプールを示唆。")
    else:
        print("  → **季節差の証拠なし**。ただし事前登録で認めた通り、これは『差が無い』ではなく")
        print("     **『この標本では検出できない』**。効果量とCIで解釈すること。")


def main():
    p = argparse.ArgumentParser(description="季節差の対比較（旗59）")
    p.add_argument("--cosore-dir", required=True)
    a = p.parse_args()
    root = Path(a.cosore_dir); desc = pd.read_csv(root / "description.csv")
    print("=== 旗59：季節でメモリは弱まるか（対比較・事前登録 0a2a563）===")
    print("  d = ACF1(休眠期) − ACF1(生育期)。二値化せず、符号反転の対順列検定で判定する。")
    fo = collect(root, desc, forest=True)
    nf = collect(root, desc, forest=False)
    report(fo, "森林", confirmatory=False)
    report(nf, "非森林", confirmatory=True)
    print("\n  留保（事前登録通り）：休眠期は呼吸が小さく S/N が悪化。温度分位の季節定義は")
    print("  フェノロジーとの一致を仮定。対8サイトでは小さな差は原理的に検出できない＝")
    print("  p≥0.05 は帰無の支持ではない。森林側は探索的で、非森林と食い違っても")
    print("  『再現しなかった』とは言えない（標本の性格が違う）。")


if __name__ == "__main__":
    main()
