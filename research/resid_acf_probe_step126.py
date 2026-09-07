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


def measure(T, th, R, t=None, port=None):
    """旗44 の設計行列で当てはめ、残差・回帰子・θ の自己相関と corr(θ,T) を返す。

    `t` は各点の**時刻（秒）**、`port` は**チャンバー識別子**（整数コード）。
    与えると、**欠測行を落とした後に「チャンバー内で時刻順」に並べ直してから**隣接判定を作る
    ——**「同じチャンバーの、刻み 1 つぶんだけ離れたペア」だけで自己相関を測る**ため。

    **道具の欠陥 #59（旗135）**：以前は「落とす前の長さ N−1 の `gap_ok`」をそのまま
    受け取っており、**欠測が 1 つも無いデータセット以外は全部 broadcast エラーで落ちていた**
    （39 本中 36 本が SKIP）。
    **道具の欠陥 #60（旗135）**：`load_cosore` は**チャンバーを分けない**。
    回帰（旗44）は行順に依らないので影響を受けないが、**自己相関は行順に依る**。
    分けずに測ると、**同時刻の別チャンバー同士の相関を「lag-1 自己相関」と取り違える**
    （刻みの中央値が 0 秒になる JASSAL・SIHI で顕著）。
    """
    ok = np.isfinite(T) & np.isfinite(th) & np.isfinite(R) & (R > 0)
    T, th, R = T[ok], th[ok], R[ok]
    if len(T) < 1000 or th.std() == 0:
        return None
    g, step, order, n_port = None, float("nan"), np.arange(len(T)), 1
    if t is not None:
        tt = np.asarray(t, float)[ok]
        pp = (np.zeros(len(tt), int) if port is None
              else np.asarray(port)[ok].astype(np.int64))
        n_port = int(len(np.unique(pp)))
        order = np.lexsort((tt, pp))                # **チャンバー主・時刻従**で並べ替え
        tt, pp = tt[order], pp[order]
        d = np.diff(tt)
        same = pp[:-1] == pp[1:]                    # チャンバーをまたぐペアは使わない
        pos = same & (d > 0)
        if pos.sum():
            step = float(np.median(d[pos]))
        g = same & (d > 0) & np.isclose(d, step, rtol=0.02)
    Tc = T - T.mean()
    thz = (th - th.mean()) / th.std()
    A = _design(Tc, thz, quad=True)                 # [Tc, Tc², θz, Tc·θz, 1]
    resid = np.log(R) - A @ np.linalg.lstsq(A, np.log(R), rcond=None)[0]
    r_res, n_res = acf1_adjacent(resid[order], g)
    r_x, _ = acf1_adjacent(A[order, 3], g)
    r_th, _ = acf1_adjacent(th[order], g)
    return {"n": int(len(T)), "acf_resid": r_res, "acf_x": r_x, "acf_th": r_th,
            "n_pairs": n_res, "step": step, "n_port": n_port,
            "corr_thT": float(np.corrcoef(th, T)[0, 1])}


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
    return ok and selftest_chambers(n)


def selftest_chambers(n=20000, phi=0.9):
    """**対照（道具の欠陥 #60）：チャンバーを混ぜたら壊れ、分けたら戻るか。**

    **同じ時刻に 2 本のチャンバーが独立に測っている**データを合成し、行を交互に並べる
    （COSORE の `CSR_PORT` があるデータセットの実際の並び）。
    **`port` を渡さなければ「同時刻の別チャンバー同士の相関」を lag-1 と取り違えるはず**——
    **その誤りが実際に起きることを見てから、直った版で φ を取り戻せることを確かめる。**
    """
    print("\n=== 対照：チャンバーを混ぜたら壊れるか／分ければ戻るか（欠陥 #60）===")
    a = synth_ar1(-0.5, phi, n, 11, 0.0, 0.99)
    b = synth_ar1(-0.5, phi, n, 12, 0.0, 0.99)
    T = np.empty(2 * n); th = np.empty(2 * n); R = np.empty(2 * n)
    for k, s in ((0, a), (1, b)):
        T[k::2], th[k::2], R[k::2] = s[0], s[1], s[2]
    t = np.repeat(np.arange(n, dtype=float) * 1800.0, 2)   # 同時刻に 2 本
    port = np.tile(np.array([0, 1]), n)
    mix = measure(T, th, R, t)                              # port を渡さない＝混ざったまま
    sep = measure(T, th, R, t, port)                        # 分けた版
    print(f"  入れた φ(残差) = {phi:.2f}")
    print(f"  {'port を渡さない（混在）':<26} 残差ACF1 = {mix['acf_resid']:+.3f}  "
          f"刻み = {mix['step']:.0f}秒  隣接ペア = {mix['n_pairs']}")
    print(f"  {'port を渡す（分離）':<28} 残差ACF1 = {sep['acf_resid']:+.3f}  "
          f"刻み = {sep['step']:.0f}秒  隣接ペア = {sep['n_pairs']}")
    broke = not (abs(mix["acf_resid"] - phi) <= 0.05)
    fixed = abs(sep["acf_resid"] - phi) <= 0.05
    print("  " + ("○混在では φ を外し、分離すれば ±0.05 で取り戻せる＝欠陥 #60 は実在し、直っている。"
                 if (broke and fixed) else
                 f"**▲対照が期待どおりでない（混在で外れた={broke} / 分離で戻った={fixed}）。**"))
    return broke and fixed


def _ports(path, n_expect):
    """`CSR_PORT` を、`load_cosore` が残した行と**同じ順・同じ本数**で取り出す。

    `load_cosore` は「時刻が読めない行」と「`Rs` が欠測の行」を落とす。**同じ条件をここで再現する。**
    **長さが合わなければ `None` を返す**——**チャンバーを取り違えるくらいなら、混在のまま測って
    その旨を印字するほうがまし**（**黙って別の行に貼り付けない**）。
    """
    import pandas as pd
    raw = pd.read_csv(path)
    if "CSR_PORT" not in raw.columns:
        return None
    tcol = ("CSR_TIMESTAMP_BEGIN" if "CSR_TIMESTAMP_BEGIN" in raw.columns
            else "CSR_TIMESTAMP_END")
    keep = (pd.to_datetime(raw[tcol], errors="coerce").notna()
            & pd.to_numeric(raw["CSR_FLUX_CO2"], errors="coerce").notna())
    port = pd.factorize(raw.loc[keep, "CSR_PORT"])[0]
    return port if len(port) == n_expect else None


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
          f"{'θのACF1':>9} {'corr(θ,T)':>10} {'隣接ペア':>9} {'室数':>5}")
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
            t = df.index.to_series().astype("int64").to_numpy() / 1e9   # 秒
            port = _ports(f, len(df))                      # **チャンバー識別子**（欠陥 #60）
            m = measure(df["Tsoil"].to_numpy(), df["SM"].to_numpy(), df["Rs"].to_numpy(),
                        t, port)
        except Exception as e:                             # **黙って飛ばさない**（旗85 の作法）
            print(f"  {ds:<32} SKIP {type(e).__name__}: {str(e)[:40]}")
            continue
        if m is None:
            print(f"  {ds:<32} SKIP 点不足/θ が一定")
            continue
        rows.append((ds, m))
        print(f"  {ds:<32} {m['n']:>7} {m['step'] / 60:>6.0f}分 {m['acf_resid']:>9.3f} "
              f"{m['acf_x']:>11.3f} {m['acf_th']:>9.3f} {m['corr_thT']:>10.2f} "
              f"{m['n_pairs']:>9} {m['n_port']:>5}{'' if port is not None else ' ★port不整合'}")

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
