"""旗126 の宿題（GATE-22）：**実データの残差と回帰子は、どれくらい自己相関しているか。**

## なぜ要るのか

**旗126 は合成で「残差が AR(1) のとき旗44 のブロックブート CI は被覆を保つか」を測った。**
**だが実サイトがどの φ に当たるかは合成からは言えない。**
**旗126 で分かったのは、被覆の挙動が二つの量に強く依存することである**——
**(a) 残差の自己相関 φ、(b) 回帰子 `Tc·θz` の自己相関（＝実質は θ の自己相関）。**
**(b) が白いと、(a) をいくら強くしても推定量の分散は膨らまない**（旗126 の追補）。
**＝実データのこの 2 つを測らないと、旗126 の水準表をどこに当てればよいかが決まらない。**

**これは下調べであって検定ではない。** **`analyze` の判定も旗44 の規則も一切呼ばない。**

## 何をするか

**旗44 と同じ入り口（`load_cosore`）・同じ設計行列（`_design(quad=True)`）で当てはめ、**

  ・**当てはめ残差の lag-1 自己相関**（**刻みが揃った隣接ペアだけで測る**——
    **欠測をまたいだペアを混ぜると、自己相関は実際より低く出る**）
  ・**回帰子 `Tc·θz` の lag-1 自己相関**（旗126 の追補が要点にした量）
  ・**θ 自身の lag-1 自己相関**
  ・**`corr(θ,T)`**（**GATE-21 が欲しがっている量。同じ 1 回の実行で出る**）

をサイトごとに印字し、**最後に分布（中央値・四分位・旗126 の水準への振り分け）**を出す。

    .venv/bin/python research/resid_acf_probe_step126.py                          # ★自己検証（合成・データ不要）
    .venv/bin/python research/resid_acf_probe_step126.py --cosore-dir /mnt/hdd/cosore-0.7.0

**自己検証を先に走らせること**——**既知の φ を入れて、この道具が φ を取り戻せるかを確かめる。**
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ar1_coverage_step126 import synth_ar1  # noqa: E402
from q10_confound_step44 import _design  # noqa: E402
from runlog import tee_stdout  # noqa: E402

PHI_BANDS = ((0.0, 0.25, "ほぼ白(φ<0.25)"), (0.25, 0.7, "中(0.25–0.7)"),
             (0.7, 0.95, "強(0.7–0.95)"), (0.95, 1.01, "極強(≥0.95)"))


def acf1_adjacent(x, gap_ok=None):
    """**刻みが揃った隣接ペアだけ**で測る lag-1 自己相関。

    `gap_ok[i]` は「`x[i]` と `x[i+1]` が刻み 1 つぶんだけ離れている」の真偽。
    **`None` なら全ペアを使う**（合成＝等間隔のとき）。
    戻り値 `(r, 使ったペア数)`。**使えるペアが 100 未満なら `(nan, 数)`。**
    """
    x = np.asarray(x, float)
    a, b = x[:-1], x[1:]
    m = np.isfinite(a) & np.isfinite(b)
    if gap_ok is not None:
        m &= gap_ok
    if m.sum() < 100:
        return float("nan"), int(m.sum())
    a, b = a[m], b[m]
    if a.std() == 0 or b.std() == 0:
        return float("nan"), int(m.sum())
    return float(np.corrcoef(a, b)[0, 1]), int(m.sum())


def measure(T, th, R, gap_ok=None):
    """旗44 の設計行列で当てはめ、残差・回帰子・θ の自己相関と corr(θ,T) を返す。"""
    ok = np.isfinite(T) & np.isfinite(th) & np.isfinite(R) & (R > 0)
    T, th, R = T[ok], th[ok], R[ok]
    g = None if gap_ok is None else gap_ok[ok[:-1] & ok[1:]] if False else gap_ok
    if len(T) < 1000 or th.std() == 0:
        return None
    Tc = T - T.mean()
    thz = (th - th.mean()) / th.std()
    A = _design(Tc, thz, quad=True)                 # [Tc, Tc², θz, Tc·θz, 1]
    resid = np.log(R) - A @ np.linalg.lstsq(A, np.log(R), rcond=None)[0]
    r_res, n_res = acf1_adjacent(resid, g)
    r_x, _ = acf1_adjacent(A[:, 3], g)
    r_th, _ = acf1_adjacent(th, g)
    return {"n": int(len(T)), "acf_resid": r_res, "acf_x": r_x, "acf_th": r_th,
            "n_pairs": n_res, "corr_thT": float(np.corrcoef(th, T)[0, 1])}


def selftest(n=20000):
    """**既知の φ を入れて取り戻せるかを確かめる**（データ不要・これを先に走らせる）。"""
    print("=== 自己検証：既知の φ / φ_θ を、この道具が取り戻せるか ===")
    print(f"  {'入れた φ(残差)':>14} {'入れた φ_θ':>10} | {'測れた 残差ACF1':>15} "
          f"{'測れた 回帰子ACF1':>16} {'測れた θのACF1':>14}")
    ok = True
    for phi, phi_th in ((0.0, 0.0), (0.9, 0.0), (0.98, 0.0), (0.0, 0.99), (0.9, 0.99)):
        s = synth_ar1(-0.5, phi, n, 7, 0.0, phi_th)
        m = measure(s[0], s[1], s[2])
        print(f"  {phi:>14.2f} {phi_th:>10.2f} | {m['acf_resid']:>15.3f} "
              f"{m['acf_x']:>16.3f} {m['acf_th']:>14.3f}")
        if abs(m["acf_resid"] - phi) > 0.05:
            ok = False
    print("\n  " + ("○残差の φ を ±0.05 で取り戻せている＝この道具は使ってよい。"
                    if ok else "**▲取り戻せていない。実データに当てる前に直すこと。**"))
    print("  **注：回帰子の ACF1 は φ_θ より低く出る**——`Tc·θz` は θ と Tc の積で、")
    print("  Tc 側の白い雑音が混ざるため。**φ_θ の推定値ではなく、回帰子そのものの赤さの指標である。**")
    return ok


def main():
    p = argparse.ArgumentParser(description="旗126 の宿題：実データの残差・回帰子の自己相関（下調べ）")
    p.add_argument("--cosore-dir")
    p.add_argument("--igbp", default="forest")
    a = p.parse_args()
    tee_stdout("step126_acf")

    if not a.cosore_dir:
        selftest()
        print("\n  実データを測るには --cosore-dir /mnt/hdd/cosore-0.7.0 を付ける（GATE-22）。")
        return

    import pandas as pd
    from cosore_memory_step40 import load_cosore

    root = Path(a.cosore_dir)
    desc = pd.read_csv(root / "description.csv")
    print(f"=== 旗126 の宿題：実データの自己相関（{a.igbp}）===")
    print("  **下調べであって検定ではない。旗44 の判定は一切呼んでいない。**")
    print(f"  {'dataset':<32} {'n':>7} {'刻み':>8} {'残差ACF1':>9} {'回帰子ACF1':>11} "
          f"{'θのACF1':>9} {'corr(θ,T)':>10}")
    rows = []
    for _, dd in desc.iterrows():
        ds = str(dd["CSR_DATASET"])
        if a.igbp and a.igbp.lower() not in str(dd.get("CSR_IGBP", "")).lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            df, st, sm = load_cosore(f)
            if "Tsoil" not in df or "SM" not in df:
                continue
            idx = df.index.to_series()
            dt = idx.diff().dt.total_seconds().to_numpy()[1:]
            step = float(np.nanmedian(dt))
            gap_ok = np.isclose(dt, step, rtol=0.02)       # **隣接＝刻み 1 つぶんのペアだけ**
            m = measure(df["Tsoil"].to_numpy(), df["SM"].to_numpy(), df["Rs"].to_numpy(), gap_ok)
        except Exception as e:                             # **黙って飛ばさない**（旗85 の作法）
            print(f"  {ds:<32} SKIP {type(e).__name__}: {str(e)[:40]}")
            continue
        if m is None:
            print(f"  {ds:<32} SKIP 点不足/θ が一定")
            continue
        rows.append((ds, m))
        print(f"  {ds:<32} {m['n']:>7} {step / 60:>6.0f}分 {m['acf_resid']:>9.3f} "
              f"{m['acf_x']:>11.3f} {m['acf_th']:>9.3f} {m['corr_thT']:>10.2f}")

    if not rows:
        print("\n  **1 本も測れなかった。** 入力を疑うこと（旗64 の作法）。")
        return

    print(f"\n  === 分布（n={len(rows)} 本）===")
    for key, lab in (("acf_resid", "残差の ACF1"), ("acf_x", "回帰子 Tc·θz の ACF1"),
                     ("acf_th", "θ の ACF1"), ("corr_thT", "corr(θ,T)")):
        v = np.array([m[key] for _, m in rows], float)
        v = v[np.isfinite(v)]
        if not len(v):
            print(f"    {lab:<22} **測れた本数 0**")
            continue
        print(f"    {lab:<22} 中央値 {np.median(v):+.3f}  四分位 [{np.percentile(v, 25):+.3f}, "
              f"{np.percentile(v, 75):+.3f}]  最小 {v.min():+.3f}  最大 {v.max():+.3f}  n={len(v)}")

    print("\n  === 旗126 の水準への振り分け（残差の ACF1）===")
    v = np.array([m["acf_resid"] for _, m in rows], float)
    for lo, hi, lab in PHI_BANDS:
        k = int(np.sum(np.isfinite(v) & (v >= lo) & (v < hi)))
        print(f"    {lab:<18} {k:>3} 本")
    print("    **旗126 が試したのは φ=0/0.5/0.9/0.98。** 上の分布がその範囲の外なら、")
    print("    **旗126 の結論はそのままでは当てはまらない**（水準を足して測り直すこと）。")
    print("\n  === GATE-21 の分（corr(θ,T) が −0.85 より強い負相関の側にある本数）===")
    c = np.array([m["corr_thT"] for _, m in rows], float)
    print(f"    corr(θ,T) ≤ −0.85 のサイト：{int(np.sum(c <= -0.85))} 本 / {len(c)} 本")
    for ds, m in rows:
        if m["corr_thT"] <= -0.85:
            print(f"      {ds}  corr(θ,T)={m['corr_thT']:+.2f}")


if __name__ == "__main__":
    main()
