"""旗151：旗44 の「生の陽性」を、点推定ではなく CI で数え直す（本文に書く分数の言い方を直すため）。

旗44 の `verdict()` は ▲（温度エイリアシングの産物）を `d_noquad > 0` という**点推定の符号だけ**で
決めている（`q10_confound_step44.py` 94行）。★は曲率ありモデルの CI が 0 を跨がないことを要求する。
＝分子と分母で基準が非対称であり、「生の陽性 N のうち M 本が交絡の産物」という分数は
そのままでは定義が揃っていない。ここでは曲率なしモデルの d にも同じブロックブート CI を付け、
「曲率なしでも CI が 0 を跨がない」＝**厳密な意味での生の陽性**を数える。

★旗44 の判定（★/▲/×/△）は作り直さない（旗102 の作法）。数え直すのは本文の言い方のためだけ。

    .venv/bin/python research/q10_rawpositive_step151.py                          # 門①（合成 G7/G8）
    .venv/bin/python research/q10_rawpositive_step151.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore
from q10_confound_step44 import T0_LT, _boot_d, _fit_d, _synth, verdict


def analyze2(T, th, R):
    """旗44 の `analyze` と同じ手続きに、曲率なしモデルの CI を 1 本足しただけのもの。"""
    ok = np.isfinite(T) & np.isfinite(th) & np.isfinite(R) & (R > 0)
    T, th, R = T[ok], th[ok], R[ok]
    if len(T) < 1000 or T.max() - T.min() < 5 or th.std() == 0:
        return {"note": "点不足/温度レンジ不足"}
    Tc = T - T.mean(); thz = (th - th.mean()) / th.std(); lnR = np.log(R)
    b_nq, d_nq = _fit_d(Tc, thz, lnR, quad=False)
    b_q, d_q = _fit_d(Tc, thz, lnR, quad=True)
    q10_dry, q10_wet = np.exp(10 * (b_q - d_q)), np.exp(10 * (b_q + d_q))
    if not (0.3 < q10_dry < 30 and 0.3 < q10_wet < 30) or abs(d_q) > 0.3:
        return {"note": f"当てはめ崩壊(Q10={q10_dry:.2g}→{q10_wet:.2g})=判定不能"}
    return {"n": int(len(T)), "d_noquad": d_nq, "d_quad": d_q,
            "ci": _boot_d(Tc, thz, lnR, quad=True),          # 旗44 と同じ（seed も同じ）
            "ci_noquad": _boot_d(Tc, thz, lnR, quad=False),  # ★本周が足した 1 本
            "q10_dry": float(q10_dry), "q10_wet": float(q10_wet),
            "corr_thT": float(np.corrcoef(th, T)[0, 1])}


def raw_positive(res):
    """厳密な生の陽性＝曲率なしモデルでも CI の下端が 0 より上。"""
    ci = res.get("ci_noquad")
    if ci is None:
        return None
    return ci[0] > 0


def _gate():
    print("=== 旗151 門①（対照）：厳密な生の陽性の基準が働くか ===")
    ok = {}
    for kind, lab, want in [("alias", "G7 陰性＝水分は感度に効かない（温度エイリアスのみ）", "生の陽性かつ▲"),
                            ("true", "G8 陽性＝真に水分が感度を変える", "生の陽性かつ★")]:
        T, th, R = _synth(kind)
        r = analyze2(T, th, R)
        rp = raw_positive(r)
        v = verdict(r)
        print(f"  {lab}")
        print(f"    d(曲率なし)={r['d_noquad']:+.4f} CI={r['ci_noquad']}  → 生の陽性={rp}")
        print(f"    d(曲率あり)={r['d_quad']:+.4f} CI={r['ci']}  → 旗44 判定={v}")
        got = bool(rp) and (v.startswith("▲") if kind == "alias" else v.startswith("★"))
        ok[kind] = got
        print(f"    要求＝{want} → {'○合格' if got else '×不合格'}\n")
    print(f"  門①：G7={'○' if ok['alias'] else '×'} G8={'○' if ok['true'] else '×'}"
          f" → {'2 本とも合格。実データを判定として読んでよい。' if all(ok.values()) else '不合格。実データを判定として読まない。'}")
    return all(ok.values())


def main():
    p = argparse.ArgumentParser(description="旗44 の生の陽性を CI で数え直す")
    p.add_argument("--cosore-dir"); p.add_argument("--igbp", default="forest")
    p.add_argument("--month", type=int, nargs="+", default=None)
    a = p.parse_args()

    if not a.cosore_dir:
        _gate()
        return

    import pandas as pd
    print("※ 実データの前に門①を走らせる（事前登録どおり・結果を読む前）\n")
    gate_ok = _gate()

    root = Path(a.cosore_dir); desc = pd.read_csv(root / "description.csv")
    print(f"\n=== 旗151 実データ：厳密な生の陽性（{a.igbp}）===")
    print(f"  {'dataset':<32} {'d(曲率なし)':>11} {'CI(曲率なし)':>19} {'生の陽性':>8}  旗44 判定")
    rows = []
    for _, dd in desc.iterrows():
        ds = str(dd["CSR_DATASET"]); igbp = str(dd.get("CSR_IGBP", ""))
        if a.igbp and a.igbp.lower() not in igbp.lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            df, st, sm = load_cosore(f, a.month)
            if "Tsoil" not in df or "SM" not in df:
                continue
            res = analyze2(df["Tsoil"].to_numpy(), df["SM"].to_numpy(), df["Rs"].to_numpy())
        except Exception as e:
            print(f"  {ds:<32} SKIP {type(e).__name__}"); continue
        if "note" in res:
            continue
        v = verdict(res); rp = raw_positive(res); ci = res["ci_noquad"]
        s_ci = f"[{ci[0]:+.4f},{ci[1]:+.4f}]" if ci else "—"
        print(f"  {ds:<32} {res['d_noquad']:>+11.4f} {s_ci:>19} "
              f"{('○' if rp else '×') if rp is not None else '?':>8}  {v}")
        rows.append((ds, res["d_noquad"], rp, v[0]))

    print("\n  === まとめ ===")
    n = len(rows)
    strict = [r for r in rows if r[2]]
    point = [r for r in rows if r[1] > 0]
    print(f"    データセット数                                  {n}")
    print(f"    点推定だけの『生の陽性』(d_noquad>0)             {len(point)}")
    print(f"    ★厳密な生の陽性(曲率なしでも CI が 0 を跨がない) {len(strict)}")
    for mark, lab in [("★", "曲率制御後も残る"), ("▲", "制御で消失"), ("×", "制御後は逆"), ("△", "制御後CI0跨ぎ")]:
        sub = [r for r in rows if r[3] == mark]
        print(f"      うち旗44 {mark}（{lab}）: {len(sub)} 本中 厳密な生の陽性 {sum(1 for r in sub if r[2])} 本")
    tri = [r for r in rows if r[3] == "▲"]
    print(f"\n    ★H18（▲ 10 本のうち 3 本以上が厳密には生の陽性でない）："
          f"該当 {sum(1 for r in tri if not r[2])} 本 / {len(tri)} 本")
    star = [r for r in rows if r[3] == "★"]
    print(f"    ★H19（★ は全本が厳密にも生の陽性）：{sum(1 for r in star if r[2])} 本 / {len(star)} 本")
    print(f"\n  門①は {'合格' if gate_ok else '不合格'}。旗44 の判定（★/▲/×/△）は作り直していない。")


if __name__ == "__main__":
    main()
