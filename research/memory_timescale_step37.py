"""旗37：呼吸の"記憶"の時間スケールを測る — 違和感駆動で未観測の法則を絞る。

**違和感**：文献（Migliavacca 2011, Stoy 2007, Cranko Page 2022, Cable 2013）は
「RECO のメモリ＝最近の光合成による基質供給、ラグ 1〜5 日」でほぼ収束している。だが旗25 で
**5日累積 Rg を足しても GER 日残差の自己相関が全く落ちなかった**（全生態系・8サイト）。
文献の想定するスケールで説明できていない＝これが違和感。so 窓幅を総当たりして「どの時間スケールなら
記憶が消えるか」を測り、**未観測の法則が"どの時間スケールに棲むか"を特定する**。

やること：観測変数から作れる遅い項（累積 Rg/Ta/θ/P/VPD の窓平均）を窓幅 W=1〜30 日でスイープし、
各 W について GER 日残差の自己相関がどこまで落ちるかを測る（旗25 の一般化＝日次で回帰）。
  ・1〜5 日で落ちる → 基質供給（文献既知, 我々の新規性なし）
  ・10〜30 日で落ちる → 水分/フェノロジーの長い記憶（文献の antecedent moisture: Cable 2013 2〜10週と整合）
  ・**どの窓でも落ちない → アクセスできる観測変数からは作れない＝真に未観測（衛星が要る根拠が確定）**
さらに残差の自己相関関数(ACF)から **e-folding 時間（記憶が 1/e に落ちる日数）** を出し、
「記憶そのものの時間スケール」を文献のラグ（1〜5日 / 2〜10週）と直接比較する。

    python research/memory_timescale_step37.py                      # 合成で検証
    python research/memory_timescale_step37.py --site JP-Tak --qc-max 1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

DRIVERS = ["Rg", "Ta", "VPD", "Ts", "th", "gH", "gLE", "P"]
SLOW_SRC = ["Rg", "Ta", "th", "P", "VPD"]      # 遅い項を作る元（観測変数のみ）
WINDOWS = [1, 2, 3, 5, 7, 10, 14, 21, 30]      # 窓幅[日]（文献のラグ帯を跨ぐ）


def _z(x):
    x = np.asarray(x, float); s = np.nanstd(x)
    return (x - np.nanmean(x)) / s if s > 0 else x - np.nanmean(x)


def _autocorr(x, lag=1):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) <= lag + 2:
        return np.nan
    a, b = x[:-lag], x[lag:]
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def efolding_days(x, maxlag=30):
    """残差 ACF が 1/e(≈0.368) を下回る最初のラグ[日]＝記憶の時間スケール。"""
    acf = [_autocorr(x, l) for l in range(1, maxlag + 1)]
    for i, r in enumerate(acf, start=1):
        if np.isfinite(r) and r < 1 / np.e:
            return i, acf
    return np.nan, acf          # maxlag 内で落ちない＝記憶が非常に長い


def fit_residual(daily, slow_cols=()):
    """日次で GER を回帰し、R² と残差を返す（旗25 と同じ思想＝日次で完結）。"""
    Y = daily["GER"].to_numpy(float)
    base = [v for v in DRIVERS if v in daily]
    zc = {v: _z(daily[v].to_numpy(float)) for v in base}
    cols = [zc[v] for v in base]
    for a, b in [("Ta", "th"), ("Rg", "VPD"), ("Ta", "VPD")]:      # 交互作用
        if a in zc and b in zc:
            cols.append(zc[a] * zc[b])
    for v in ("Ta", "th", "VPD"):                                   # 二次
        if v in zc:
            cols.append(zc[v] ** 2)
    for c in slow_cols:                                             # 遅い項
        if c in daily:
            cols.append(_z(daily[c].to_numpy(float)))
    X = np.column_stack(cols + [np.ones(len(Y))])
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ coef
    ss = np.sum((Y - Y.mean()) ** 2)
    r2 = 1 - np.sum(resid ** 2) / ss if ss > 0 else np.nan
    return r2, resid


def add_slow_window(daily, w):
    """窓幅 w 日の遅い項（観測変数の後方窓平均）を付ける。列名を返す。"""
    d = daily.copy()
    names = []
    for v in SLOW_SRC:
        if v in d:
            nm = f"cum{v}{w}"
            d[nm] = d[v].rolling(w, min_periods=max(1, w // 2)).mean()
            names.append(nm)
    doy = d.index.dayofyear.to_numpy().astype(float)               # フェノロジー
    d["doy_sin"] = np.sin(2 * np.pi * doy / 365.25)
    d["doy_cos"] = np.cos(2 * np.pi * doy / 365.25)
    names += ["doy_sin", "doy_cos"]
    return d, names


def sweep(daily, windows=WINDOWS):
    """窓幅スイープ：各 W で残差自己相関がどこまで落ちるか。"""
    r2_0, res_0 = fit_residual(daily)
    ac_0 = _autocorr(res_0)
    ef_0, _ = efolding_days(res_0)
    rows = [{"w": 0, "r2": r2_0, "ac": ac_0, "ef": ef_0}]
    for w in windows:
        d, names = add_slow_window(daily, w)
        d = d.dropna()
        if len(d) < 40:
            continue
        r2, res = fit_residual(d, names)
        ef, _ = efolding_days(res)
        rows.append({"w": w, "r2": r2, "ac": _autocorr(res), "ef": ef})
    return rows


def make_synth(kind, days=900, seed=0):
    """記憶の時間スケールを仕込んだ合成（窓スイープが正体を当てられるか）。"""
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2001-06-01", periods=days, freq="D")
    Ta = 20 + 5 * np.sin(2 * np.pi * idx.dayofyear / 365.25) + rng.normal(0, 1.5, days)
    Rg = 300 + 80 * np.sin(2 * np.pi * idx.dayofyear / 365.25) + rng.normal(0, 40, days)
    th = np.clip(0.3 + rng.normal(0, 0.05, days), 0.05, 0.6)
    VPD = np.clip(0.5 + 0.06 * (Ta - 15), 0.1, None)
    P = np.clip(rng.normal(0, 1, days), 0, None)
    base = np.exp(0.06 * (Ta - 20))
    if kind == "short":        # 3日スケールの観測性メモリ（基質供給を模す）
        slow = 1.5 * _z(pd.Series(Rg).rolling(3, min_periods=1).mean().to_numpy())
    elif kind == "long":       # 21日スケールの観測性メモリ（水分の長い記憶）
        slow = 1.5 * _z(pd.Series(th).rolling(21, min_periods=1).mean().to_numpy())
    else:                       # 未観測の AR(0.9) メモリ（どの観測窓でも作れない）
        h = np.zeros(days)
        for i in range(1, days):
            h[i] = 0.9 * h[i - 1] + rng.normal(0, 0.4)
        slow = 1.5 * _z(h)
    GER = 2.0 * base + slow + rng.normal(0, 0.15, days)
    return pd.DataFrame({"GER": np.clip(GER, 1e-3, None), "Rg": Rg, "Ta": Ta,
                         "VPD": VPD, "Ts": Ta - 1, "th": th, "P": P,
                         "gH": Rg * 0.2, "gLE": Rg * 0.3}, index=idx)


def _report(rows, tag):
    print(f"\n  === {tag} ===")
    print(f"  {'窓W[日]':>7} {'R²':>6} {'日残差ACF':>9} {'e-folding[日]':>13}")
    for r in rows:
        w = "なし" if r["w"] == 0 else str(r["w"])
        ef = "≥30" if not np.isfinite(r["ef"]) else f"{r['ef']:.0f}"
        print(f"  {w:>7} {r['r2']:>6.3f} {r['ac']:>+9.2f} {ef:>13}")
    base = rows[0]["ac"]
    best = min((r for r in rows[1:]), key=lambda r: r["ac"], default=None)
    if best is None:
        return
    drop = base - best["ac"]
    print(f"  → 最良は W={best['w']}日（ACF {base:+.2f}→{best['ac']:+.2f}, 低下={drop:+.2f}）")
    if best["ac"] < 0.3:
        print(f"     ✅ W={best['w']}日の遅い項で記憶が消える＝**観測から作れる**"
              f"（文献の基質供給1-5日/水分記憶2-10週と照合可）")
    elif drop >= 0.3:
        print(f"     ○ 大きく下がるが残る＝一部は観測の遅い項、一部は未観測")
    else:
        print(f"     ★ **どの窓(1-30日)でも落ちない＝アクセスできる観測変数からは作れない**"
              f"＝真に未観測（衛星プロキシが要る根拠）")


def main():
    p = argparse.ArgumentParser(description="呼吸の記憶の時間スケールを窓スイープで測る")
    p.add_argument("--site")
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--qc-max", type=int, default=None)
    a = p.parse_args()

    if not a.site:
        print("=== 旗37 合成検証：窓スイープで記憶の時間スケールを当てられるか ===")
        for kind, lab in [("short", "3日スケールの観測性メモリ（基質供給を模す）"),
                          ("long", "21日スケールの観測性メモリ（水分の長い記憶）"),
                          ("hidden", "未観測 AR(0.9) メモリ（観測から作れない）")]:
            _report(sweep(make_synth(kind)), lab)
        print("\n  → 3日/21日は該当窓で ACF が落ち、未観測はどの窓でも落ちないのが期待。")
        return

    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    from japanflux_pn.run_robustness import get_site_years
    cfg = AnalysisConfig(qc_max=a.qc_max) if a.qc_max is not None else AnalysisConfig()
    years, mo = get_site_years(a.site)
    ms = sorted(a.month or mo)
    raw = load_raw_all(get_site(a.site), cfg)
    raw = raw[raw.index.month.isin(ms)]
    keep = ["GER"] + [v for v in DRIVERS if v in raw.columns]
    daily = raw[keep].groupby(raw.index.normalize()).mean().dropna()
    if len(daily) < 100:
        print(f"日数不足({len(daily)})"); return
    qtag = f"QC≤{a.qc_max}" if a.qc_max is not None else "gap-fill込み"
    print(f"=== 旗37 実データ {a.site}（呼吸の記憶の時間スケール, {len(years)}年{qtag}, "
          f"{len(daily)}日）===")
    _report(sweep(daily), f"{a.site} GER 日残差")
    print("\n  読み方：文献は『メモリ＝基質供給, ラグ1〜5日』(Migliavacca2011/Stoy2007)、"
          "『antecedent moisture 2〜10週』(Cable2013)。")
    print("    どの窓でも落ちなければ＝これらでは説明できない＝我々のデータが示す違和感の核心。")
    print("  留保：GERは分割派生量(旗32)。日次・夏のみ。e-folding は30日上限。")


if __name__ == "__main__":
    main()
