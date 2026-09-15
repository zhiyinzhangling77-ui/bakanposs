"""**旗157 の前提の確認**（旗97 の作法・事前登録を書く**前**に走らせる）。

## 何を確かめるか

**旗156 が最優先に置いた問い**：**20 年規模の `none` で `Δ` の帰無が緩い**
（`θ→γH` 0.073・`θ→γLE` 0.097・どちらも名目 0.05 より上）**のはなぜか。**
**旗144 の `vshare`（θ の分散に年成分が占める取り分）が効いているか。**

**その問いを検定にする前に、設計が乗る前提を 4 つ実測する**：

  ・**P-1**  `year_contrast_step144.decompose` の `r_between + r_within` が、
             `phase_asymmetry_step139` の検定が使っている `partial_spearman` の `r` と**同じ量か**。
             **違えば、分解しているのは検定と別の量になる**（旗144 の欠陥 #78 と同じ形の事故）。
  ・**P-2**  `synth("none", years=20, seed=7000+i, thin=False)` を**作り直すと
             旗156 の CSV の `delta_h`・`delta_le` が再現するか**（種は `i` だけで決まるはず）。
             **再現しなければ、CSV の `sig` と本周の分解を同じ標本の上で突き合わせられない。**
  ・**P-3**  600 本ぶんの**作り直し＋分解に何分かかるか**（前景 600 s の上限に当たるか）。
  ・**P-4**  `none`(20 年) の MAM・SON で **`vshare` が実際いくつか**。
             **旗144 の表（0.019〜0.462）のどのあたりか**を先に知らないと、予測を書けない。

    .venv/bin/python research/premise_check_step157.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import phase_asymmetry_step139 as M                      # noqa: E402
from year_contrast_step144 import decompose              # noqa: E402
from moisture_control_atlas_step31 import partial_spearman  # noqa: E402

CSV = Path(__file__).resolve().parent / "logs" / "step156_pairs_none_20260914_200040.csv"
SEED0 = 7000
YEARS = 20


def _seasons(d):
    sp = d[np.isin(d.index.month, M.SPRING)]
    au = d[np.isin(d.index.month, M.AUTUMN)]
    return sp, au


def _dec(sub, col):
    return decompose(sub[col].to_numpy(float), sub["th"].to_numpy(float),
                     sub["Rg"].to_numpy(float), sub.index.year.to_numpy())


def p1(nrep=5):
    print("\n【P-1】`decompose` の和は検定の `partial_spearman` と同じ量か")
    worst = 0.0
    for i in range(nrep):
        d = M.synth("none", years=YEARS, seed=SEED0 + i, thin=False)
        for tag, sub in zip(("MAM", "SON"), _seasons(d)):
            for k, col in (("h", "gH"), ("le", "gLE")):
                b, w, _ = _dec(sub, col)
                r_ps = float(partial_spearman(sub[col].to_numpy(float),
                                              sub["th"].to_numpy(float),
                                              [sub["Rg"].to_numpy(float)])[0])
                e = abs((b + w) - r_ps)
                worst = max(worst, e)
                if i == 0:
                    print(f"    i=0 {tag} {k}: r_ps {r_ps:+.6f}"
                          f" ／ between {b:+.6f} ＋ within {w:+.6f} = {b + w:+.6f}"
                          f" ／差 {e:.2e}")
    print(f"    → **{nrep} 本 × 2 季 × 2 系列の最大差 = {worst:.3e}**")
    print(f"       {'**同じ量とみなせる**' if worst < 1e-9 else '**別の量＝設計を変える**'}")
    return worst


def p2(nrep=8):
    print("\n【P-2】作り直しで旗156 の CSV の Δ が再現するか")
    if not CSV.exists():
        print(f"    **CSV が無い**（{CSV}）→ 再現の確認ができない")
        return None
    ref = np.genfromtxt(CSV, delimiter=",", names=True)
    print(f"    CSV：{CSV.name}（{len(ref)} 行）")
    worst = 0.0
    for i in range(nrep):
        d = M.synth("none", years=YEARS, seed=SEED0 + i, thin=False)
        got = M._delta_hat(d)
        row = ref[ref["i"] == i]
        if not len(row):
            print(f"    i={i}: CSV に無い")
            continue
        for k, cn in (("h", "delta_h"), ("le", "delta_le")):
            e = abs(float(got[k][2]) - float(row[cn][0]))
            worst = max(worst, e)
            if i < 2:
                print(f"    i={i} {k}: 作り直し {got[k][2]:+.6f}"
                      f" ／ CSV {float(row[cn][0]):+.6f} ／差 {e:.2e}")
    print(f"    → **{nrep} 本の最大差 = {worst:.3e}**"
          f"（CSV は小数 6 桁で書かれているので 5e-7 までは丸め）")
    print(f"       {'**同じ標本＝突き合わせてよい**' if worst < 1e-6 else '**別の標本＝突き合わせない**'}")
    return worst


def p3(nrep=10):
    print("\n【P-3】600 本ぶんの所要時間")
    t0 = time.time()
    for i in range(nrep):
        d = M.synth("none", years=YEARS, seed=SEED0 + i, thin=False)
        M._delta_hat(d)
        for sub in _seasons(d):
            for col in ("gH", "gLE"):
                _dec(sub, col)
    per = (time.time() - t0) / nrep
    print(f"    1 本 {per:.3f} s → 600 本 {per * 600 / 60:.1f} 分"
          f"（`calendar_driven` 100 本を足して {per * 700 / 60:.1f} 分）")
    print(f"    {'**前景で走らせてよい**' if per * 700 < 500 else '**背景タスクに投げる**'}")
    return per


def p4(nrep=30):
    print("\n【P-4】`none`(20 年) の `vshare`（θ の分散に年成分が占める取り分）")
    acc = {("MAM", "h"): [], ("MAM", "le"): [], ("SON", "h"): [], ("SON", "le"): []}
    for i in range(nrep):
        d = M.synth("none", years=YEARS, seed=SEED0 + i, thin=False)
        for tag, sub in zip(("MAM", "SON"), _seasons(d)):
            for k, col in (("h", "gH"), ("le", "gLE")):
                acc[(tag, k)].append(_dec(sub, col)[2])
    for key, v in acc.items():
        a = np.asarray(v, float)
        print(f"    {key[0]} {key[1]}：vshare 中央 {np.median(a):.3f}"
              f"（四分位 {np.percentile(a, 25):.3f}–{np.percentile(a, 75):.3f}・{nrep} 本）")
    print("    **旗144 の表（`|r_between|` 中央）との対応**：")
    print("      vshare 0.019→0.006 ／ 0.070→0.021 ／ 0.205→0.058 ／ 0.321→0.099（3 年）")
    print("      vshare 0.029→0.006 ／ 0.107→0.012 ／ 0.304→0.039 ／ 0.462→0.059（20 年）")
    return {k: float(np.median(v)) for k, v in acc.items()}


def main() -> int:
    print("**旗157 の前提の確認**（事前登録を書く前・旗97 の作法）")
    print(f"**季節：MAM {M.SPRING} ／ SON {M.AUTUMN}／年数 {YEARS}／種 {SEED0}+i**")
    p1()
    p2()
    p3()
    p4()
    print("\n**ここまでが前提。ここから先（事前登録）は、この 4 つを踏まえて書く。**")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
