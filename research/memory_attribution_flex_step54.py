"""旗54：旗45（メモリの正体を候補で当てる）を、旗52/53 の較正を通してやり直す。

旗45 は「先行湿潤・先行水分・深層水分・非線形Birch・熱慣性のどれも ~4日メモリを説明せず（12/12）」
と結論し、それが本研究の中心（「観測の隙間」）の根拠になっている。だが旗45 は
**基本モデルが線形**（ln Rs ~ 表層Tsoil + 表層SM）で、旗52 が示した通り
**非線形系への線形当てはめは、それ自体が自己相関残差を作る**。すると：

  ・説明すべき「メモリ」の相当部分が**検出器の産物**だった可能性
  ・その大きな産物に埋もれて、**本物の水履歴シグナルが見えていなかった**可能性（＝偽の帰無）

の両方がある。よって基本モデルを**非線形基底**（Lloyd-Taylor項・二次・交互作用・log θ）に替え、
適格条件も旗53 の較正値（柔軟 ACF1 ≥ 0.64 かつ e-fold ≤ 7日）に替えて、同じ候補比較をやり直す。
プラセボ（同次元・位相ずらし）併走は旗45 のまま維持する。

  ・やはりどの候補も説明しない → **「観測の隙間」は最も強い交絡制御を通しても立つ**。
  ・今度は候補が説明する → 旗45 の帰無は**線形当てはめの取り残しに埋もれていただけ**＝結論を差し替える。

    python research/memory_attribution_flex_step54.py                                   # 合成で検証
    python research/memory_attribution_flex_step54.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import _acf_gap, _efold_gap
from memory_attribution_step45 import load_daily, build_blocks, _synth
from synthetic_tower_step52 import T0_LT

ACF_THR = 0.64          # 旗53 の較正値（帰無0.49 と 4日仕込み0.80 の中点）
EFOLD_MAX = 7           # 短メモリに限定（季節メモリを除く）
BEAT_PLACEBO = 0.10     # プラセボをこれだけ上回った候補のみ採用（旗45 と同じ）


def flex_basis(T, W):
    """非線形基底：Lloyd-Taylor 項・二次・交互作用・log θ（旗53 と同じ形）。"""
    cols = [T]
    with np.errstate(divide="ignore", invalid="ignore"):
        cols.append(320.0 * (1.0 / (10 - T0_LT) - 1.0 / (T - T0_LT)))
    cols.append(T ** 2)
    if W is not None:
        cols += [W, W ** 2, T * W, np.log(np.clip(W, 1e-3, None))]
    return cols


def _fit(y, cols):
    X = np.column_stack(cols + [np.ones(len(y))])
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    if ok.sum() < max(60, 5 * X.shape[1]):
        return None, np.nan
    coef = np.linalg.lstsq(X[ok], y[ok], rcond=None)[0]
    res = np.full(len(y), np.nan); res[ok] = y[ok] - X[ok] @ coef
    ss = np.sum((y[ok] - y[ok].mean()) ** 2)
    return res, float(1 - np.sum(res[ok] ** 2) / ss) if ss > 0 else np.nan


def analyze(daily):
    if "T_sh" not in daily:
        return {"note": "土壌温度なし"}
    y = np.log(daily["Rs"].where(daily["Rs"] > 0)).to_numpy()
    T = daily["T_sh"].to_numpy()
    W = daily["SM_sh"].to_numpy() if "SM_sh" in daily else None
    base = flex_basis(T, W)
    res0, r2_0 = _fit(y, base)
    if res0 is None:
        return {"note": "点不足"}
    ac0, ef0 = _acf_gap(res0, 1), _efold_gap(res0)
    if not (np.isfinite(ac0) and np.isfinite(r2_0)):
        return {"note": "推定不能"}
    if r2_0 < 0.3 or ac0 < ACF_THR or ef0 > EFOLD_MAX:
        return {"note": f"対象外(R2={r2_0:.2f}, ACF1={ac0:.2f}, e-fold={ef0})"}

    out = {"r2_0": r2_0, "ac0": ac0, "ef0": ef0, "cands": {}}
    for name, blk in build_blocks(daily).items():
        res, r2 = _fit(y, base + [blk[c].to_numpy() for c in blk.columns])
        if res is None:
            continue
        ac, ef = _acf_gap(res, 1), _efold_gap(res)
        out["cands"][name] = {"r2": r2, "ac": ac, "ef": ef,
                              "dac": ac0 - ac if np.isfinite(ac) else np.nan}
    return out


def verdict(res):
    if "note" in res:
        return "―" + res["note"], None
    cs = res["cands"]
    pl_raw = next((v["dac"] for k, v in cs.items() if "プラセボ" in k), None)
    real = {k: v for k, v in cs.items() if "プラセボ" not in k and np.isfinite(v["dac"])}
    if not real:
        return "―候補なし", None
    # **プラセボが作れなかったサイトで★を出さない**（旗54 第1回の欠陥：GUTIERREZ は
    # プラセボ行が無いまま ★ が付いていた＝過剰適合の対照が無い判定だった）
    if pl_raw is None or not np.isfinite(pl_raw):
        best0 = max(real, key=lambda k: real[k]["dac"])
        return f"△プラセボ無し＝判定保留（最大は {best0} Δ{real[best0]['dac']:+.2f}）", None
    pl = float(pl_raw)
    best = max(real, key=lambda k: real[k]["dac"])
    d = real[best]["dac"]
    if d > pl + BEAT_PLACEBO and real[best]["ac"] < ACF_THR:
        return f"★{best}が記憶を説明", best
    if d > pl + BEAT_PLACEBO:
        return f"○{best}が部分的に説明", best
    return "―どの候補も説明せず（観測の外側）", None


def run_synth():
    print("=== 旗54 合成検証：非線形基底にしても、仕込んだ正体を当てられるか ===")
    print("  旗45 と同じ仕込み（線形Birch／非線形Birch／深層水分／観測外の隠れAR）を、")
    print(f"  非線形基底の基本モデル＋較正値（ACF1≥{ACF_THR}・e-fold≤{EFOLD_MAX}日）で判定する。\n")
    for kind, lab in [("birch", "正体=湿潤パルス(線形Birch)"),
                      ("birch_nl", "正体=非線形Birch"), ("deep", "正体=深層水分"),
                      ("unknown", "正体=観測外の隠れAR（当ててはいけない）")]:
        d, _ = _synth(kind)
        r = analyze(d)
        v, _ = verdict(r)
        if "note" in r:
            print(f"  {lab:<30} {v}"); continue
        print(f"  {lab:<30} 基本 R2={r['r2_0']:.2f} ACF1={r['ac0']:.2f} e-fold={r['ef0']}日")
        for k, c in r["cands"].items():
            print(f"      {k:<24} ΔACF1={c['dac']:+.2f}  → {c['ac']:+.2f}")
        print(f"      → {v}\n")
    print("  期待：Birch/深層は該当候補が★or○、隠れARは『どの候補も説明せず』。")
    print("  ※非線形基底は候補ブロックの一部（θの二次・logθ）と重なるので、旗45 より")
    print("    候補の見かけの寄与は小さくなる＝より保守的な検定になる。")


def run_real(cosore_dir, igbp, months):
    root = Path(cosore_dir)
    desc = pd.read_csv(root / "description.csv")
    print(f"=== 旗54 実データ：非線形基底＋較正値でメモリの正体を当て直す（{igbp or '全'}）===")
    tally = {}
    n_target = 0
    for _, dd in desc.iterrows():
        ds = str(dd["CSR_DATASET"]); ig = str(dd.get("CSR_IGBP", ""))
        if igbp and igbp.lower() not in ig.lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            daily, meta = load_daily(f, months)
            r = analyze(daily)
        except Exception:
            continue
        if "note" in r:
            continue
        n_target += 1
        v, best = verdict(r)
        key = v.split("が")[0].lstrip("★○―") if best else "どの候補も説明せず"
        tally[key] = tally.get(key, 0) + 1
        print(f"\n  {ds}  (基本 R2={r['r2_0']:.2f} ACF1={r['ac0']:.2f} e-fold={r['ef0']}日)")
        for k, c in r["cands"].items():
            print(f"      {k:<24} ΔACF1={c['dac']:+.2f}  → {c['ac']:+.2f}")
        print(f"    → {v}")
    print(f"\n  === まとめ（適格 {n_target} サイト）===")
    if not tally:
        print("    適格サイトなし（柔軟R²≥0.3・ACF1≥%.2f・e-fold≤%d日）" % (ACF_THR, EFOLD_MAX))
    for k, v in sorted(tally.items(), key=lambda x: -x[1]):
        print(f"    {k:<24} {v}")
    print("\n  読み：『どの候補も説明せず』が多数＝**最も強い交絡制御（非線形基底）を通しても、")
    print("        観測されている水・熱の履歴では説明できない**＝「観測の隙間」は保たれる。")
    print("        候補が説明する側に変われば、旗45 の帰無は線形当てはめの取り残しに埋もれていた")
    print("        だけ＝結論を差し替える。")
    print("  留保：非線形基底は候補ブロックと一部重なるため、候補の寄与は旗45 より小さく出る")
    print("        （＝より保守的）。降水は全COSOREサイトで欠測のままで、Birch はθ増加の代理検定。")


def main():
    p = argparse.ArgumentParser(description="メモリの帰属を非線形基底でやり直す")
    p.add_argument("--cosore-dir"); p.add_argument("--igbp", default="forest")
    p.add_argument("--month", type=int, nargs="+", default=None)
    a = p.parse_args()
    if a.cosore_dir:
        run_real(a.cosore_dir, a.igbp, a.month)
    else:
        run_synth()


if __name__ == "__main__":
    main()
