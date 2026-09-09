"""**追補 G（探索）**：A-3 側の `Δ(θ→γLE)` を、`ES-FcO` と同じ道具・同じ月・同じ帰無で測る。

## なぜ書くか

**旗142 が残し、旗153 が最優先に置いた宿題である**——
`ES-FcO` の `Δ(θ→γLE) = −0.218`（両側 p=0.0145）と**向きを比べる相手が、A-3 側に無かった**。
**旗88 は春と秋の `r` を別々に印字しているが、その差 `Δ` と、`Δ` の帰無を一度も出していない。**

## ★探索である（`PREREGISTRATION_step154.md` の冒頭）

**結論には使わない。** `NOVELTY_ASSESSMENT.md` の A-3 の主張にも `MEASURED_ONLY_SPINE.md` の骨格にも
足さない。**判定表（★／▲）は当てない**——判定表は `ES-FcO` の識別検定のために書かれたものだからである。

## 道具を書き換えない

**`phase_asymmetry_step139.py` は import するだけで一行も変えない**（旗142 の作法：
判定を見る周に道具をいじらない）。**帰無も較正済みの実装（追補 D-2）をそのまま呼ぶ。**

**ただし旗88 の月（春 4–5・秋 9–10）で出すときだけ、`step139` の module 変数
`SPRING`／`AUTUMN` を一時的に差し替える**（`_months` の with 節）。
**ファイルは書き換えていないが、これは実行時の書き換えである**——**必ず有効な月を印字する。**

## 門①（`PREREGISTRATION_step154.md` 3 節・**実行前に固定**）

  ・**G-A（検出器）**  旗88 の 12 個（3 サイト × 春 4–5／秋 9–10 × `γLE`／`γH`）を ±0.02 で再現
  ・**G-B（水準）**    `synth("none", years=20, thin=False)` 200 本で `p_delta<0.05` ≤ 0.07
  ・**G-C（到達）**    `synth("calendar_driven", years=20, thin=False)` 100 本で
                       `p_delta<0.05 かつ Δ(h)<0` ≥ 0.80

    .venv/bin/python research/phase_asymmetry_a3_step154.py --gates
    .venv/bin/python research/phase_asymmetry_a3_step154.py --real
"""
from __future__ import annotations

import argparse
import contextlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import phase_asymmetry_step139 as M          # noqa: E402  **一行も書き換えない**
from evaporation_regime_step36 import daily_energy   # noqa: E402
from runlog import tee_stdout                        # noqa: E402

SITES = ("US-Wkg", "US-Whs", "US-SRM")

# **旗88 の記録**（`FLAGS_LOG.md` 旗88 の H1 判定表）。**実行前に写した。**
G88 = {
    "US-Wkg": {"sp": {"le": +0.46, "h": +0.01}, "au": {"le": +0.72, "h": -0.47}},
    "US-Whs": {"sp": {"le": +0.64, "h": -0.13}, "au": {"le": +0.84, "h": -0.55}},
    "US-SRM": {"sp": {"le": +0.56, "h": -0.05}, "au": {"le": +0.81, "h": -0.56}},
}
G88_MONTHS = {"sp": (4, 5), "au": (9, 10)}

TOL = 0.02            # G-A の許容（旗139 の G1_TOL と同じ・**実行前に固定**）
LEVEL_MAX = 0.07      # G-B の上限（追補 D-3 と同じ）
REACH_MIN = 0.80      # G-C の下限（追補 D-3 と同じ）
SYNTH_YEARS = 20      # A-3 の年数（18–22）に合わせる
ESFCO_D_LE = -0.218   # 旗142 の `ES-FcO` の Δ(θ→γLE)


@contextlib.contextmanager
def _months(spring, autumn):
    """**`step139` の季節定義を一時的に差し替える**（旗88 の月で出すときだけ）。"""
    old = (M.SPRING, M.AUTUMN)
    M.SPRING, M.AUTUMN = tuple(spring), tuple(autumn)
    try:
        yield
    finally:
        M.SPRING, M.AUTUMN = old


def _load(site, months):
    d, nyr = daily_energy(site, list(months), None)
    return d, nyr


# ------------------------------------------------------------------ 門①
def gate_a() -> bool:
    """**G-A（検出器の健全性）**：旗88 の 12 個を ±0.02・符号一致で再現する。"""
    print(f"\n  【G-A】旗88 の 12 個を再現する（許容 ±{TOL:.2f}・符号一致・**実行前に固定**）")
    hits = 0
    for site in SITES:
        for tag, months in G88_MONTHS.items():
            d, nyr = _load(site, months)
            got = M.season_r(d)
            line = [f"    {site} {tag}({months[0]}–{months[1]} 月) "
                    f"{len(d)} 日／{d.index.year.nunique()} 年（健全年 {nyr}）"]
            for k, nm in (("le", "θ→γLE"), ("h", "θ→γH")):
                r = float(got[k][0])
                ref = G88[site][tag][k]
                hit = bool(np.isfinite(r) and np.sign(r) == np.sign(ref) and abs(r - ref) <= TOL)
                hits += int(hit)
                line.append(f"      {nm} = {r:+.3f}  旗88 = {ref:+.2f}"
                            f"／差 {r - ref:+.3f} → {'○' if hit else '**×**'}")
            print("\n".join(line))
    ok = hits == 12
    print(f"    → G-A は {'○合格' if ok else '**×不合格**'}（{hits}/12）")
    return ok


def _synth_rate(kind: str, reps: int, nperm: int, need_neg_h: bool) -> dict:
    """**合成 `reps` 本で `p_delta < 0.05` の割合を出す**（G-B・G-C 共通の中身）。"""
    cnt = {"h": 0, "le": 0}
    cnt_dir = 0
    deltas = {"h": [], "le": []}
    n_ok = 0
    for i in range(reps):
        d = M.synth(kind, years=SYNTH_YEARS, seed=7000 + i, thin=False)
        r = M.perm_delta_p(d, nperm=nperm, seed=i)
        if not all(np.isfinite(r[k]["delta"]) for k in ("h", "le")):
            continue
        n_ok += 1
        for k in ("h", "le"):
            cnt[k] += int(r[k]["p"] < M.ALPHA)
            deltas[k].append(r[k]["delta"])
        if need_neg_h:
            cnt_dir += int(r["h"]["p"] < M.ALPHA and r["h"]["delta"] < 0)
        if (i + 1) % 25 == 0:
            print(f"      … {i + 1}/{reps} 本")
    out = {k: cnt[k] / max(n_ok, 1) for k in ("h", "le")}
    out["dir"] = cnt_dir / max(n_ok, 1)
    out["n"] = n_ok
    for k in ("h", "le"):
        a = np.asarray(deltas[k], float)
        out[f"med_{k}"] = float(np.median(a)) if len(a) else np.nan
    return out


def gate_b(reps: int, nperm: int) -> bool:
    """**G-B（水準）**：A-3 の標本規模（20 年・間引きなし）で `Δ` の偽陽性が ≤ 0.07 か。"""
    print(f"\n  【G-B】水準：`none` を {reps} 本（{SYNTH_YEARS} 年・間引きなし・各 {nperm} 置換）")
    r = _synth_rate("none", reps, nperm, need_neg_h=False)
    ok = r["h"] <= LEVEL_MAX and r["le"] <= LEVEL_MAX
    for k, nm in (("h", "θ→γH"), ("le", "θ→γLE")):
        print(f"    {nm}：p<0.05 の割合 {r[k]:.3f}（要求 ≤ {LEVEL_MAX:.2f}）"
              f" {'○' if r[k] <= LEVEL_MAX else '**×**'}／Δ の中央値 {r[f'med_{k}']:+.3f}")
    print(f"    → G-B は {'○合格' if ok else '**×不合格**'}（有効 {r['n']}/{reps} 本）")
    return ok


def gate_c(reps: int, nperm: int) -> bool:
    """**G-C（到達）**：同じ規模で本物の季節差（`calendar_driven`）を捕まえるか。"""
    print(f"\n  【G-C】到達：`calendar_driven` を {reps} 本"
          f"（{SYNTH_YEARS} 年・間引きなし・各 {nperm} 置換）")
    r = _synth_rate("calendar_driven", reps, nperm, need_neg_h=True)
    ok = r["dir"] >= REACH_MIN
    print(f"    θ→γH：p<0.05 かつ Δ<0 の割合 {r['dir']:.3f}（要求 ≥ {REACH_MIN:.2f}）"
          f" {'○' if ok else '**×**'}／Δ の中央値 {r['med_h']:+.3f}")
    print(f"    （参考）θ→γLE：p<0.05 の割合 {r['le']:.3f}／Δ の中央値 {r['med_le']:+.3f}")
    print(f"    → G-C は {'○合格' if ok else '**×不合格**'}（有効 {r['n']}/{reps} 本）")
    return ok


def run_gates(reps_b: int, reps_c: int, nperm: int) -> dict:
    print("\n  ===== 門①（対照）——**3 本とも合格しなければ実データの p を読まない** =====")
    a = gate_a()
    b = gate_b(reps_b, nperm)
    c = gate_c(reps_c, nperm)
    print(f"\n  **門①：G-A {'○' if a else '×'}／G-B {'○' if b else '×'}／G-C {'○' if c else '×'}**")
    return {"a": a, "b": b, "c": c}


# ------------------------------------------------------------------ 実データ
def _one_site(site, spring, autumn, nperm: int) -> dict | None:
    with _months(spring, autumn):
        d, nyr = _load(site, tuple(spring) + tuple(autumn))
        res = M.perm_result(d, nperm=nperm, seed=0)
    if res is None:
        print(f"    {site}: **下限未満＝出さない**（{len(d)} 日）")
        return None
    sp_n = int(np.isin(d.index.month, tuple(spring)).sum())
    au_n = int(np.isin(d.index.month, tuple(autumn)).sum())
    print(f"    {site}: 春 {sp_n} 日／秋 {au_n} 日／{d.index.year.nunique()} 年（健全年 {nyr}）")
    for k, nm in (("le", "θ→γLE（主）"), ("h", "θ→γH（併記）")):
        r = res[k]
        print(f"      {nm}：春 {r['r_sp']:+.3f}（片側 p={r['p_sp']:.4f}）"
              f" / 秋 {r['r_au']:+.3f}（片側 p={r['p_au']:.4f}）")
        print(f"        **Δ = 秋 − 春 = {r['delta']:+.3f}（両側 p={r['p_delta']:.4f}）**"
              f"／巡回シフトの通り数 春 {r['n_shift_sp']}・秋 {r['n_shift_au']}")
    return res


def run_real(nperm: int, gates_ok: dict) -> None:
    print("\n  ===== 実データ（**探索・結論には使わない**） =====")
    print("  **判定表（★／▲）は当てない**（`PREREGISTRATION_step154.md` 冒頭）。")
    if not gates_ok["a"]:
        print("  **G-A が落ちた＝読み込み経路が旗88 と違う。何も出さずに止める。**")
        return
    if not gates_ok["b"]:
        print("  **G-B が落ちた＝`p` は読まない。Δ の点推定だけを記述として残す。**")
    if not gates_ok["c"]:
        print("  **G-C が落ちた＝`p ≥ 0.05` を「差が無い」と読まない**（検出力不足と区別できない）。")

    print(f"\n  ----- 本体：`ES-FcO` と同じ月（MAM {(3, 4, 5)} ／ SON {(9, 10, 11)}）"
          f"・追補 D-2 の置換 {nperm} 本 -----")
    main = {s: _one_site(s, (3, 4, 5), (9, 10, 11), nperm) for s in SITES}

    print(f"\n  ----- 併記（**判定に使わない**）：旗88 の月（春 {G88_MONTHS['sp']}"
          f" ／ 秋 {G88_MONTHS['au']}） -----")
    alt = {s: _one_site(s, G88_MONTHS["sp"], G88_MONTHS["au"], nperm) for s in SITES}

    print("\n  ----- **`ES-FcO` との向きの比較**（記述） -----")
    print(f"    `ES-FcO`（CRO・3 年・春に緑）：Δ(θ→γLE) = {ESFCO_D_LE:+.3f}（両側 p=0.0145・旗142）")
    for s in SITES:
        if main[s] is None:
            continue
        dl = main[s]["le"]["delta"]
        same = "同じ向き" if np.sign(dl) == np.sign(ESFCO_D_LE) else "**逆向き**"
        print(f"    {s}：Δ(θ→γLE) = {dl:+.3f}（p={main[s]['le']['p_delta']:.4f}）→ {same}")

    _score(main, alt)


def _score(main: dict, alt: dict) -> None:
    """**★事前予測 H20・H21・H22 の勝敗**（`PREREGISTRATION_step154.md` 5 節・数字を見る前に固定）。"""
    print("\n  ----- ★事前予測の勝敗（**書いたとおりの形で数える**） -----")
    ref = {"US-Wkg": +0.26, "US-Whs": +0.20, "US-SRM": +0.25}   # 旗88 の月での差（点推定）

    n20 = 0
    for s in SITES:
        if main[s] is None:
            continue
        dl = main[s]["le"]["delta"]
        hit = abs(dl - ref[s]) <= 0.10
        n20 += int(hit)
        print(f"    H20 {s}：Δ(θ→γLE) = {dl:+.3f} 対 旗88 月の差 {ref[s]:+.2f}"
              f"／|差| {abs(dl - ref[s]):.3f} ≤ 0.10 → {'○' if hit else '**×**'}")
    print(f"    **H20 は {'当たり' if n20 == 3 else '**外れ**'}（{n20}/3・要求 3/3）**")

    n21 = sum(1 for s in SITES if main[s] is not None and main[s]["le"]["p_delta"] < M.ALPHA)
    print(f"    **H21 は {'当たり' if n21 >= 2 else '**外れ**'}"
          f"（Δ(θ→γLE) の p<0.05 が {n21}/3・要求 2 以上）**")

    n22 = 0
    for s in SITES:
        if main[s] is None:
            continue
        dh, dl = main[s]["h"]["delta"], main[s]["le"]["delta"]
        hit = abs(dh) > abs(dl)
        n22 += int(hit)
        print(f"    H22 {s}：|Δ(θ→γH)| {abs(dh):.3f} > |Δ(θ→γLE)| {abs(dl):.3f}"
              f" → {'○' if hit else '**×**'}")
    print(f"    **H22 は {'当たり' if n22 == 3 else '**外れ**'}（{n22}/3・要求 3/3）**")

    print("\n    **併記（多重性・実行前に決めた形）**：サイト間の族は 3 検定＝Bonferroni なら 0.0167。")
    for s in SITES:
        if main[s] is None:
            continue
        p = main[s]["le"]["p_delta"]
        print(f"      {s}：p={p:.4f} → 0.0167 {'を下回る' if p < 0.05 / 3 else '**を下回らない**'}")
    print("    **補正は当てていない**（穴 #83 と同じ形の穴を本周も塞がない・探索だから）。")

    print("\n    **旗88 の月での Δ（併記・判定に使わない）**：")
    for s in SITES:
        if alt[s] is None:
            continue
        print(f"      {s}：Δ(θ→γLE) = {alt[s]['le']['delta']:+.3f}"
              f"（p={alt[s]['le']['p_delta']:.4f}）／"
              f"Δ(θ→γH) = {alt[s]['h']['delta']:+.3f}（p={alt[s]['h']['p_delta']:.4f}）")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gates", action="store_true", help="門① G-A/G-B/G-C だけ走らせる")
    ap.add_argument("--real", action="store_true", help="門①を走らせてから実データへ")
    ap.add_argument("--nperm", type=int, default=M.NPERM)
    ap.add_argument("--reps-b", type=int, default=200)
    ap.add_argument("--reps-c", type=int, default=100)
    a = ap.parse_args()
    tee_stdout("step154")
    print("**追補 G（探索）：A-3 側の Δ(θ→γLE) を `ES-FcO` と同じ道具で測る**")
    print("**事前登録 `research/PREREGISTRATION_step154.md`（実データの Δ を見る前に確定）**")
    print(f"**季節（既定）：MAM {M.SPRING} ／ SON {M.AUTUMN}**")
    if not (a.gates or a.real):
        ap.error("--gates か --real を指定する")
    g = run_gates(a.reps_b, a.reps_c, a.nperm)
    if a.real:
        run_real(a.nperm, g)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
