"""旗24：呼吸の"遅い未観測駆動"を絶対値空間で探す（季節交絡を日平均集約で制御）。

旗23（アノマリ）で呼吸は瞬間駆動と脱結合(R²=0.04)＝残差はノイズと分かった。だが呼吸の意味ある
構造は絶対値側（旗20の応答曲面）に在る。so 絶対値空間で観測駆動（＋交互作用）を当てはめ、
**残差を"日平均"に集約してから自己相関**を見る。狙い：
  ・瞬間の駆動ノイズは日平均で消える → 残る日スケールの残差自己相関＝**遅い未観測駆動**の足跡
    （基質・フェノロジー・深層土壌水は日〜季節スケールで効く）。季節交絡と混ざらず取れる。
  ・日残差 vs 日集約した他の観測変数 → 入れ忘れの"遅い"観測駆動か。
真に未観測（日残差が自己相関）か、観測で説明可（日残差がノイズ）かを分離する。

    python research/respiration_unobserved_abs_step24.py                    # 合成で検証
    python research/respiration_unobserved_abs_step24.py --site JP-Tak --years 1999 ... --month 7 8
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

DRIVERS = ["Rg", "Ta", "VPD", "Ts", "th", "gH", "gLE", "P"]


def _z(x):
    x = np.asarray(x, float); s = np.nanstd(x)
    return (x - np.nanmean(x)) / s if s > 0 else x - np.nanmean(x)


def _fit_resid(Y, drivers: dict):
    """絶対値 GER を観測駆動＋主要交互作用で当てはめ、R² と残差を返す。"""
    names = list(drivers)
    zc = {v: _z(drivers[v]) for v in names}
    cols = [zc[v] for v in names]
    for a, b in [("Ta", "th"), ("Rg", "VPD"), ("Ta", "VPD")]:
        if a in zc and b in zc:
            cols.append(zc[a] * zc[b])
    for v in ("Ta", "th", "VPD"):
        if v in zc:
            cols.append(zc[v] ** 2)
    X = np.column_stack(cols + [np.ones(len(Y))])
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ coef
    ss = np.sum((Y - Y.mean()) ** 2)
    r2 = 1 - np.sum(resid ** 2) / ss if ss > 0 else np.nan
    return r2, resid


def _autocorr1(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)] - np.nanmean(x)
    if len(x) < 3 or x.std() == 0:
        return np.nan
    return float(np.corrcoef(x[:-1], x[1:])[0, 1])


def analyze(df: pd.DataFrame):
    """df: 絶対値、datetime index、GER＋DRIVERS。日平均残差の自己相関などを返す。"""
    Y = df["GER"].to_numpy(float)
    drivers = {v: df[v].to_numpy(float) for v in DRIVERS if v in df}
    r2, resid = _fit_resid(Y, drivers)
    day = pd.Series(resid, index=df.index).groupby(df.index.normalize()).mean()
    ger_day = df["GER"].groupby(df.index.normalize()).mean()
    ac_resid = _autocorr1(day.to_numpy())
    ac_ger = _autocorr1(ger_day.to_numpy())
    # 日残差 vs 日集約観測変数（入れ忘れの遅い駆動候補）
    leftover = []
    for v in drivers:
        vday = df[v].groupby(df.index.normalize()).mean()
        j = pd.concat([day.rename("r"), vday.rename("v")], axis=1).dropna()
        if len(j) > 3:
            leftover.append((v, abs(float(np.corrcoef(j["r"], j["v"])[0, 1]))))
    leftover.sort(key=lambda t: -t[1])
    return {"r2": r2, "n_days": int(day.notna().sum()),
            "ac_resid_day": ac_resid, "ac_ger_day": ac_ger, "leftover": leftover}


def make_synth(hidden=True, days=200, per_day=48, seed=0):
    rng = np.random.default_rng(seed)
    n = days * per_day
    idx = pd.date_range("2001-06-01", periods=n, freq="30min")
    doy = idx.dayofyear.to_numpy()
    hour = idx.hour.to_numpy() + idx.minute.to_numpy() / 60
    Rg = np.clip(np.sin((hour - 6) / 12 * np.pi), 0, None) * 800 + rng.normal(0, 30, n)
    Ta = 20 + 6 * np.sin((hour - 9) / 24 * 2 * np.pi) + 0.03 * (doy - doy.mean()) + rng.normal(0, 1, n)
    th = 0.3 + 0.05 * np.sin(doy / 30) + rng.normal(0, 0.01, n)
    VPD = np.clip(0.5 + 0.06 * (Ta - 15), 0.1, None) + rng.normal(0, 0.1, n)
    # 遅い未観測駆動：日スケールの AR（基質/フェノロジー的）
    hd = np.zeros(days)
    for d in range(1, days):
        hd[d] = 0.85 * hd[d - 1] + rng.normal(0, 0.3)
    hidden_hh = np.repeat(hd, per_day)
    base = np.exp(0.06 * (Ta - 20)) * (th / (0.2 + th))
    GER = 2.0 * base + (0.8 * hidden_hh if hidden else 0) + rng.normal(0, 0.2, n)
    df = pd.DataFrame({"GER": np.clip(GER, 1e-3, None), "Rg": Rg, "Ta": Ta,
                       "VPD": VPD, "Ts": Ta - 1, "th": th, "gH": Rg * 0.2,
                       "gLE": Rg * 0.3, "P": rng.normal(0, 0.5, n)}, index=idx)
    return df


def _report(res, tag):
    print(f"\n  === {tag} ===")
    print(f"  R²(観測駆動＋交互作用, 絶対値)={res['r2']:.3f}  有効日数={res['n_days']}")
    print(f"  日平均GER 自己相関(参考)={res['ac_ger_day']:+.2f}")
    print(f"  ★日平均『残差』の自己相関={res['ac_resid_day']:+.2f}"
          f"（大＝遅い未観測駆動の足跡 / ≈0＝観測＋非加法で説明可）")
    print(f"  日残差 vs 日集約観測（入れ忘れ候補・上位）: " +
          ", ".join(f"{v}:{c:.2f}" for v, c in res['leftover'][:3]))


def main():
    p = argparse.ArgumentParser(description="呼吸の遅い未観測駆動を絶対値・日平均残差で探す")
    p.add_argument("--site")
    p.add_argument("--years", type=int, nargs="+")
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    a = p.parse_args()

    if not a.site:
        print("=== 旗24 合成検証：絶対値・日平均残差で遅い未観測駆動を検出 ===")
        _report(analyze(make_synth(hidden=True)), "隠れ駆動あり（日スケールの遅い未観測）")
        _report(analyze(make_synth(hidden=False)), "隠れ駆動なし（観測＋交互作用で説明可）")
        print("\n  → 期待：ありは日残差自己相関が大、なしは≈0。＝真に未観測(遅い) vs 説明可 を分離。")
        return

    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    cfg = AnalysisConfig()
    raw_all = load_raw_all(get_site(a.site), cfg)
    ms = sorted(a.month)
    frames = []
    used = []
    for y in a.years or []:
        start = pd.Timestamp(year=y, month=ms[0], day=1)
        end = pd.Timestamp(year=y, month=ms[-1], day=1) + pd.offsets.MonthBegin(1)
        r = raw_all[(raw_all.index >= start) & (raw_all.index < end)]
        keep = ["GER"] + [v for v in DRIVERS if v in r.columns]
        r = r[keep].dropna()
        if len(r) > 100:
            frames.append(r); used.append(y)
    if not frames:
        print("有効年なし"); return
    df = pd.concat(frames)
    print(f"=== 旗24 実データ {a.site}（GER 遅い未観測, 絶対値・日平均残差, {len(used)}年, N={len(df)}）===")
    _report(analyze(df), f"{a.site} 呼吸 GER")
    print("\n  読み方：日残差自己相関が大＝観測(非加法込み)でも説明しきれない『遅い未観測駆動』の足跡")
    print("     （基質/フェノロジー/深層土壌水）＝領域知識で仮説化し衛星プロキシで検証の入口(北極星)。")
    print("     ≈0＝呼吸は観測＋非加法でほぼ説明可＝未観測は小さい、と正直に線引き。")


if __name__ == "__main__":
    main()
