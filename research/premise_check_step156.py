"""**旗156 の前提の確認**（事前登録を書く前に、前提にする事実を実測する・旗97 の作法）。

**測るのは 4 つだけ。判定は一つも含まない。**

  ・**(a) 決定論か** —— `phase_asymmetry_a3_step154.py` の `_synth_rate` の中身
    （`synth(seed=7000+i)` → `perm_delta_p(seed=i)`）を同じ `i` で二度回して完全一致するか。
    **一致するなら「先頭 200 本が旗154 の 13/200・22/200 に一致する」を門①の対照に使える。**
  ・**(b) 1 本あたりの所要**（20 年・間引きなし・2000 置換）。**reps をいくつにするかを決めるため。**
  ・**(c) 旗154 の 13/200・22/200 の Clopper-Pearson 区間**と、
    **0.07 と区別するのに要る `n`**（穴 #94 の作法：閾値を区間で示す `n` を同時に見積もる）。
  ・**(d) McNemar（discordant の exact binomial）の検出力**——
    対応の強さが分からないので、潜在相関 ρ を変えて模擬する。

    .venv/bin/python research/premise_check_step156.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import phase_asymmetry_step139 as M          # noqa: E402  **一行も書き換えない**
from runlog import tee_stdout                # noqa: E402


def cp(k: int, n: int) -> tuple[float, float]:
    lo = 0.0 if k == 0 else float(stats.beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(stats.beta.ppf(0.975, k + 1, n - k))
    return lo, hi


def power(n: int, ph: float, ple: float, rho: float, sims: int, rng) -> float:
    """**ガウシアン・コピュラで対応のある 2 値を作り、McNemar exact の検出力を出す。**"""
    zh, zle = stats.norm.ppf(ph), stats.norm.ppf(ple)
    L = np.linalg.cholesky(np.array([[1.0, rho], [rho, 1.0]]))
    hit = 0
    z = rng.standard_normal((sims, n, 2)) @ L.T
    sh = z[:, :, 0] < zh
    sle = z[:, :, 1] < zle
    b = (sh & ~sle).sum(axis=1)
    c = (sle & ~sh).sum(axis=1)
    for bb, cc in zip(b, c):
        nd = int(bb + cc)
        if nd == 0:
            continue
        hit += int(stats.binomtest(int(bb), nd, 0.5).pvalue < 0.05)
    return hit / sims


def main() -> int:
    tee_stdout("step156_premise")
    print("=== 旗156 前提の確認（判定は含まない） ===")

    print("\n[a][b] 決定論と所要（20 年・間引きなし・2000 置換）")
    out = {}
    for run in (1, 2):
        rows = []
        t0 = time.time()
        for i in range(4):
            d = M.synth("none", years=20, seed=7000 + i, thin=False)
            r = M.perm_delta_p(d, nperm=2000, seed=i)
            rows.append((i, round(float(r["h"]["delta"]), 12), round(float(r["h"]["p"]), 12),
                         round(float(r["le"]["delta"]), 12), round(float(r["le"]["p"]), 12)))
        out[run] = (rows, time.time() - t0)
    same = out[1][0] == out[2][0]
    print(f"  **二度回して完全一致：{'○ 決定論' if same else '**× 決定論ではない**'}**")
    for row in out[1][0]:
        print("    i=%d  Δh=%+.6f p_h=%.6f  Δle=%+.6f p_le=%.6f" % row)
    per = out[1][1] / 4
    print("  経過：run1 %.1f s ／ run2 %.1f s → **1 本あたり %.2f s**" % (out[1][1], out[2][1], per))
    for n in (200, 400, 600, 800, 1000):
        print("    reps=%4d → %5.1f 分" % (n, per * n / 60))

    print("\n[c] 旗154 の 200 本の Clopper-Pearson（両側 95%）")
    for nm, k in (("θ→γH ", 13), ("θ→γLE", 22)):
        lo, hi = cp(k, 200)
        print("  %s %2d/200 = %.3f  [%.3f, %.3f]" % (nm, k, k / 200, lo, hi))
    print("  真値 0.110 のとき CP 下端 > 0.07 になる n（k = round(0.110 n)）")
    for n in (200, 300, 400, 500, 600, 800, 1000):
        k = int(round(0.110 * n))
        lo, hi = cp(k, n)
        print("    n=%4d k=%3d → [%.3f, %.3f] %s" % (n, k, lo, hi, "○" if lo > 0.07 else "×"))
    print("  真値 0.065 のとき CP 上端 ≤ 0.07 になる n（k = round(0.065 n)）")
    for n in (200, 600, 1000, 2000, 3000, 5000):
        k = int(round(0.065 * n))
        lo, hi = cp(k, n)
        print("    n=%4d k=%4d → [%.3f, %.3f] %s" % (n, k, lo, hi, "○" if hi <= 0.07 else "×"))

    print("\n[d] McNemar（discordant の exact binomial・両側 α=0.05）の検出力")
    print("  真の周辺率を 0.065（H）/ 0.110（LE）に固定し、潜在相関 ρ を変えて模擬（各 2000 回）")
    rng = np.random.default_rng(0)
    for rho in (0.0, 0.5, 0.8, 0.95):
        line = []
        for n in (200, 600, 1000):
            line.append("n=%4d %.3f" % (n, power(n, 0.065, 0.110, rho, 2000, rng)))
        print("    ρ=%.2f  " % rho + " ／ ".join(line))
    print("\n  **注意：これは「真の周辺率が旗154 の点推定に等しい」と仮定した見積りである。**")
    print("  **点推定自体が標本誤差を含む（[c] の区間）ので、検出力の数字も点推定である。**")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
