"""旗55：事前登録した仮説 H1 を検定する——深部地温を測れている所ほど熱慣性が記憶を説明するか。

**予測・統計量・判定規則は `research/PREREGISTRATION_step55.md` に実行前に確定・commit 済み。**
本ファイルはその手続きを実装するだけで、結果を見てから基準を動かさない。

H1：`d_thermal = ΔACF1(熱慣性) − ΔACF1(プラセボ)` は、そのサイトの土壌温度センサの
    **最大深度 max_depth** と正の相関を持つ（Spearman、順列検定 2000 回、片側）。
主判定は**50km クラスタに潰した版**（旗48 の擬似反復の教訓）。

    python research/thermal_depth_step55.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from memory_attribution_step45 import load_daily, build_blocks
from memory_attribution_flex_step54 import flex_basis, _fit, ACF_THR, EFOLD_MAX
from cosore_memory_step40 import _acf_gap, _efold_gap
from out_of_sample_step48 import geo_clusters


def _spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 4:
        return np.nan
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    if rx.std() == 0 or ry.std() == 0:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def _perm_p(x, y, nperm=2000, seed=0):
    """x（深さ）をサイト間で並べ替える片側順列検定。"""
    r = _spearman(x, y)
    if not np.isfinite(r):
        return r, np.nan
    rng = np.random.default_rng(seed)
    null = [_spearman(rng.permutation(x), y) for _ in range(nperm)]
    null = np.array([v for v in null if np.isfinite(v)])
    return r, float((null >= r).mean()) if null.size else np.nan


def site_row(path, months):
    """1サイトの d_thermal と土壌温度センサの深さを返す（適格でなければ None）。"""
    daily, meta = load_daily(path, months)
    if "T_sh" not in daily:
        return None
    y = np.log(daily["Rs"].where(daily["Rs"] > 0)).to_numpy()
    T = daily["T_sh"].to_numpy()
    W = daily["SM_sh"].to_numpy() if "SM_sh" in daily else None
    base = flex_basis(T, W)
    res0, r2_0 = _fit(y, base)
    if res0 is None:
        return None
    ac0, ef0 = _acf_gap(res0, 1), _efold_gap(res0)
    if not (np.isfinite(ac0) and np.isfinite(r2_0)):
        return None
    if r2_0 < 0.3 or ac0 < ACF_THR or ef0 > EFOLD_MAX:      # 旗53/54 と同じ適格条件
        return None

    d = {}
    for name, blk in build_blocks(daily).items():
        res, _ = _fit(y, base + [blk[c].to_numpy() for c in blk.columns])
        if res is None:
            continue
        ac = _acf_gap(res, 1)
        if np.isfinite(ac):
            d[name] = ac0 - ac
    th = next((v for k, v in d.items() if "熱慣性" in k), None)
    pl = next((v for k, v in d.items() if "プラセボ" in k), None)
    if th is None or pl is None:            # プラセボが無いサイトは事前登録通り除外
        return None
    depths = meta.get("T_depths") or []
    if not depths:
        return None
    return {"d_thermal": float(th - pl), "max_depth": float(max(depths)),
            "n_depth": len(depths), "ac0": ac0, "ef0": ef0}


def main():
    p = argparse.ArgumentParser(description="事前登録H1：深部地温と熱慣性の説明力の相関")
    p.add_argument("--cosore-dir", required=True)
    p.add_argument("--igbp", default="forest")
    p.add_argument("--month", type=int, nargs="+", default=None)
    a = p.parse_args()

    root = Path(a.cosore_dir); desc = pd.read_csv(root / "description.csv")
    print("=== 旗55：事前登録 H1 の検定（予測は PREREGISTRATION_step55.md に実行前確定）===")
    print("  H1：熱慣性の説明力 d_thermal = ΔACF1(熱慣性) − ΔACF1(プラセボ) は max_depth と正相関\n")
    print(f"  {'dataset':<32}{'最大深度cm':>10}{'層数':>6}{'d_thermal':>11}{'基本ACF1':>10}{'ef':>5}")
    rows = []
    for _, dd in desc.iterrows():
        ds = str(dd["CSR_DATASET"]); ig = str(dd.get("CSR_IGBP", ""))
        if a.igbp and a.igbp.lower() not in ig.lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            r = site_row(f, a.month)
        except Exception:
            continue
        if r is None:
            continue
        try:
            r["lat"] = float(dd["CSR_LATITUDE"]); r["lon"] = float(dd["CSR_LONGITUDE"])
        except (TypeError, ValueError, KeyError):
            r["lat"] = r["lon"] = np.nan
        r["ds"] = ds; rows.append(r)
        print(f"  {ds:<32}{r['max_depth']:>10.1f}{r['n_depth']:>6}{r['d_thermal']:>+11.3f}"
              f"{r['ac0']:>10.2f}{r['ef0']:>5.0f}", flush=True)

    if len(rows) < 6:
        print(f"\n  適格サイト {len(rows)} 件＝**事前登録の通り検出力不足で判定しない**。")
        return
    df = pd.DataFrame(rows)
    print(f"\n  === 深さの分布（事前登録で『必ず報告する』と決めた項目）===")
    print(f"  最大深度 cm：{sorted(df['max_depth'].round(1).tolist())}")
    print(f"  層数：{sorted(df['n_depth'].tolist())}")

    r1, p1 = _perm_p(df["max_depth"].to_numpy(), df["d_thermal"].to_numpy())
    print(f"\n  --- サイト単位（参考）---  Spearman r={r1:+.2f}  片側 p={p1:.3f}")

    ok = np.isfinite(df["lat"]) & np.isfinite(df["lon"])
    if ok.sum() == len(df):
        g = geo_clusters(df["lat"].to_numpy(), df["lon"].to_numpy(), km=50.0)
        cl = df.assign(cl=g).groupby("cl")[["max_depth", "d_thermal"]].mean()
        r2, p2 = _perm_p(cl["max_depth"].to_numpy(), cl["d_thermal"].to_numpy())
        print(f"  --- **クラスタ単位（主判定）** --- {len(cl)} クラスタ  "
              f"Spearman r={r2:+.2f}  片側 p={p2:.3f}")
    else:
        r2, p2 = r1, p1
        print("  --- 座標欠損のためクラスタ縮約できず。サイト単位を主判定に代用（限界として記録）---")

    print("\n  === 事前登録の判定規則に照らす ===")
    if np.isfinite(r2) and r2 > 0 and np.isfinite(p2) and p2 < 0.05:
        print("  → ★**H1 支持**：深部地温が測れていないことが『説明できなさ』の一因＝")
        print("     観測の隙間の一部は**センサ配置の問題**。")
    elif np.isfinite(r2) and r2 > 0:
        print("  → ○**示唆どまり**：方向は合うが証拠不十分（p≥0.05）。保留として記録する。")
    else:
        print("  → ▲**H1 棄却**：SAVAGE の★は深さでは説明できない＝**1サイトの偶然**として扱い、")
        print("     旗54 の記述から『深部温度が次の一手』という示唆を落とす。")
    print("\n  留保（事前登録で先に認めた通り）：max_depth はセンサ配置の代理であって熱慣性そのものではない。")
    print("  深く測るサイトは他の面でも観測が充実している可能性（交絡）。検出力は低い。")


if __name__ == "__main__":
    main()
