"""旗72：**帰無の作り方の監査**——シャッフルが自己相関を壊していないか。

旗71 で「**自己相関があるのに日ブートで CI を出していた**」という欠陥（11件目）を見つけた。
**同じ種類の誤りが他にないか**を系統的に見る。

## 監査の結果（コード読みで確定した分）

- **旗31／旗36**：実データ経路は `blocks=年` を渡している＝**無傷**。
  `blocks` 無しの呼び出しは**合成検証（i.i.d.）**なので日ブートで正しい。
  ＝**欠陥11件目は旗71 に限局**。旗62/70 も `--block` で実行済み。
- **O-information（旗14／19／60）**：帰無は
  `japanflux_pn.information_theory.surrogate_o_information_stats` の
  **`col[rng.permutation(n)]`＝素の並べ替え**である。
  これは依存を壊すと同時に**自己相関も壊す**。実データは自己相関を持つので
  **実効標本数は N より小さく、Ω のばらつきは大きい**。なのに帰無は i.i.d. のばらつきしか持たない。
  → **σ が過小 → z が過大 → 有意が出やすい**という疑いがある。
  本研究は**呼吸の多日メモリ**そのものを主題にしているので、この疑いは軽くない。

## 本ツールがやること（**実データを要さない・コンテナ内で決着する**）

**自己相関はあるが、変数間の依存はゼロ**の系列を作る（＝**真の Ω は 0**）。
現行の帰無で z を何度も計算し、**その分布を見る**：

  ・帰無が正しければ **z ≈ N(0,1)**、**|z|>2 は約 5%**。
  ・**|z| が系統的に大きければ、帰無が壊れている**＝これまでの z は過大。

φ（自己相関の強さ）を 0（i.i.d.）と 0.8（実データ相当）で比べれば、
**自己相関が原因であることまで特定できる**。

## 修正案も同時に試す

**ブロック並べ替え**（長さ L の連続塊ごと順序を入れ替える）なら、
**塊の中の自己相関は保たれたまま**、変数間の対応だけが壊れる。
これで校正が戻るかを同じ試験で確かめる。

    python research/surrogate_audit_step72.py
    python research/surrogate_audit_step72.py --reps 60 --n 800   # 精度を上げる
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from japanflux_pn import information_theory as it


def ar1(n, phi, rng):
    """AR(1) 系列（φ=0 なら i.i.d.）。"""
    x = np.zeros(n)
    e = rng.standard_normal(n)
    if phi == 0:
        return e
    x[0] = e[0]
    for i in range(1, n):
        x[i] = phi * x[i - 1] + e[i]
    return x


def _digit(cols, bins):
    return [it.digitize_series(np.asarray(c, float), bins) for c in cols]


def z_plain(idx, bins, nsur, rng):
    """**現行の帰無**：各変数を独立に素の並べ替え（自己相関が壊れる）。"""
    om = it.o_information_indices(idx, bins, correct=True)
    st = it.surrogate_o_information_stats(idx, bins, nsur, 0.0, rng, correct=True)
    s = st["sigma"]
    return om, (om - st["mu"]) / s if s > 0 else np.nan


def _block_shuffle(a, L, rng):
    """長さ L の連続塊ごと順序を入れ替える（塊の中の自己相関は保たれる）。"""
    n = len(a)
    nb = int(np.ceil(n / L))
    blocks = [a[i * L:(i + 1) * L] for i in range(nb)]
    order = rng.permutation(nb)
    return np.concatenate([blocks[i] for i in order])[:n]


def z_block(idx, bins, nsur, rng, L):
    """**修正案の帰無**：ブロック並べ替え（自己相関を保つ）。"""
    om = it.o_information_indices(idx, bins, correct=True)
    samples = np.empty(nsur)
    for s in range(nsur):
        shuf = [_block_shuffle(c, L, rng) for c in idx]
        samples[s] = it.o_information_indices(shuf, bins, correct=True)
    mu, sd = float(np.mean(samples)), float(np.std(samples))
    return om, (om - mu) / sd if sd > 0 else np.nan


def run(phi, n, nvar, bins, nsur, reps, L, seed=0):
    """**真の Ω が 0**（変数間は独立）の系列で、z の分布を測る。"""
    out = {"plain": [], "block": []}
    for r in range(reps):
        rng = np.random.default_rng(seed + r)
        cols = [ar1(n, phi, rng) for _ in range(nvar)]     # **互いに独立**＝真の Ω=0
        idx = _digit(cols, bins)
        out["plain"].append(z_plain(idx, bins, nsur, np.random.default_rng(seed + 1000 + r))[1])
        out["block"].append(z_block(idx, bins, nsur, np.random.default_rng(seed + 2000 + r), L)[1])
    return {k: np.asarray(v, float) for k, v in out.items()}


def summarize(name, z):
    z = z[np.isfinite(z)]
    if len(z) == 0:
        return f"    {name:<22} 計算できず"
    frac = float(np.mean(np.abs(z) > 2))
    return (f"    {name:<22} z の平均 {np.mean(z):+6.2f}  sd {np.std(z):5.2f}  "
            f"**|z|>2 の割合 {frac:5.1%}**（正しければ約 5%）")


def main():
    p = argparse.ArgumentParser(description="O-information の帰無が自己相関を壊していないかの監査")
    p.add_argument("--n", type=int, default=500, help="系列長（日数相当）")
    p.add_argument("--nvar", type=int, default=4)
    p.add_argument("--bins", type=int, default=8)
    p.add_argument("--nsur", type=int, default=100)
    p.add_argument("--reps", type=int, default=30)
    p.add_argument("--block-len", type=int, default=10, help="ブロック長（メモリの時間尺度より長く）")
    a = p.parse_args()

    print("=== 旗72：O-information の帰無が自己相関を壊していないかの監査 ===")
    print("  **変数間の依存がゼロ**の系列を使う＝**真の Ω は 0**。")
    print("  帰無が正しければ z≈N(0,1) で **|z|>2 は約 5%** に収まるはず。")
    print(f"  設定：系列長 {a.n}・変数 {a.nvar}・ビン {a.bins}・サロゲート {a.nsur}・"
          f"反復 {a.reps}・ブロック長 {a.block_len}\n")

    for phi, lab in [(0.0, "φ=0.0（自己相関なし＝合成検証と同じ条件）"),
                     (0.8, "**φ=0.8（実データ相当＝多日メモリがある）**")]:
        print(f"  ━ {lab} ━")
        res = run(phi, a.n, a.nvar, a.bins, a.nsur, a.reps, a.block_len)
        print(summarize("現行（素の並べ替え）", res["plain"]))
        print(summarize("修正案（ブロック）", res["block"]))
        print()

    print("  === 読み方 ===")
    print("  φ=0 で現行が約5%なら、**実装そのものは正しい**（合成検証が通っていた理由）。")
    print("  **φ=0.8 で現行の |z|>2 が 5% を大きく超えるなら、帰無が壊れている**")
    print("  ＝**旗14/19/60 の z は過大**であり、有意判定をやり直す必要がある。")
    print("  修正案（ブロック）が φ=0.8 でも約5%に戻るなら、それが正しい帰無である。")
    print("  留保：")
    print("   ・AR(1) は実データの記憶構造の**単純化**である（実際は非線形・季節つき）。")
    print("     ここで測れるのは『**自己相関があると帰無が壊れるか**』という定性的な問いまで。")
    print("   ・ブロック長は**メモリの時間尺度より長く**採る必要がある（既定10日＝呼吸の4日より長い）。")
    print("     短すぎると自己相関を壊してしまい、現行と同じ問題が残る。")
    print("   ・反復回数が少ないと『割合』自体の誤差が大きい（30回なら ±約4%）。")


if __name__ == "__main__":
    main()
