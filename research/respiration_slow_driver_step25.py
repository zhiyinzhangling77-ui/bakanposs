"""旗25：呼吸の"遅い成分"は観測から作れる遅い駆動か、真に未観測か（決着）。

旗24 で絶対値の夏 GER は遅い成分に支配され瞬間の観測駆動は8%しか説明せず、日残差自己相関0.90＝
遅い未観測駆動の足跡が出た。だが「遅い」≠「観測不能」：基質(最近の光合成)・フェノロジー・熱履歴は
**観測変数の時間積分**から作れる。so 遅いプロキシを足して日残差自己相関が落ちるかで決着：
  ・累積Rg（過去5日平均＝光/基質履歴）／累積Ta（熱履歴/GDD）／通日 sin,cos（フェノロジー）
  ・落ちる→「観測から作れる遅い項の入れ忘れ」＝真に未観測でない（モデルに遅い項を足せば済む）。
  ・落ちない→真に未観測（衛星プロキシ SIF/SMAP の出番）。
炭素分割(GEP/GER/NEE)は使わない（累積Rg で基質を代理＝分割循環を避ける）。

    python research/respiration_slow_driver_step25.py                    # 合成で検証
    python research/respiration_slow_driver_step25.py --site JP-Tak --years 1999 ... --month 7 8
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

DRIVERS = ["Rg", "Ta", "VPD", "Ts", "th", "gH", "gLE", "P"]
WIN = 5 * 48   # 5日（30分値）


def _z(x):
    x = np.asarray(x, float); s = np.nanstd(x)
    return (x - np.nanmean(x)) / s if s > 0 else x - np.nanmean(x)


def _autocorr1(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)] - np.nanmean(x)
    if len(x) < 3 or x.std() == 0:
        return np.nan
    return float(np.corrcoef(x[:-1], x[1:])[0, 1])


def add_slow(df: pd.DataFrame) -> pd.DataFrame:
    """観測から作れる遅いプロキシ（年内で計算）：累積Rg/Ta(過去5日平均)・通日 sin/cos。"""
    d = df.copy()
    d["cumRg"] = d["Rg"].rolling(WIN, min_periods=WIN // 2).mean()
    d["cumTa"] = d["Ta"].rolling(WIN, min_periods=WIN // 2).mean()
    doy = d.index.dayofyear.to_numpy().astype(float)
    d["doy_sin"] = np.sin(2 * np.pi * doy / 365.25)
    d["doy_cos"] = np.cos(2 * np.pi * doy / 365.25)
    return d


def _fit_ac(daily, use_slow):
    """日次データで GER_day を回帰し、R² と残差の日次自己相関を返す（日次で完結＝日周期の混入なし）。"""
    Y = daily["GER"].to_numpy(float)
    base = [v for v in DRIVERS if v in daily]
    zc = {v: _z(daily[v].to_numpy(float)) for v in base}
    cols = [zc[v] for v in base]
    for a, b in [("Ta", "th"), ("Rg", "VPD"), ("Ta", "VPD")]:
        if a in zc and b in zc:
            cols.append(zc[a] * zc[b])
    for v in ("Ta", "th", "VPD"):
        if v in zc:
            cols.append(zc[v] ** 2)
    if use_slow:
        for v in ("cumRg", "cumTa", "doy_sin", "doy_cos"):
            if v in daily:
                cols.append(_z(daily[v].to_numpy(float)))
    X = np.column_stack(cols + [np.ones(len(Y))])
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ coef
    ss = np.sum((Y - Y.mean()) ** 2)
    r2 = 1 - np.sum(resid ** 2) / ss if ss > 0 else np.nan
    return r2, _autocorr1(resid)


def analyze(df):
    df = df.dropna()
    # 日次に集約してから回帰（日周期を除き、遅い成分＝日々の構造だけを見る）
    daily = df.groupby(df.index.normalize()).mean()
    daily = daily.dropna()
    ac_ger = _autocorr1(daily["GER"].to_numpy())
    r2a, aca = _fit_ac(daily, use_slow=False)
    r2b, acb = _fit_ac(daily, use_slow=True)
    return {"r2_base": r2a, "ac_base": aca, "r2_slow": r2b, "ac_slow": acb,
            "ac_ger": ac_ger, "n_days": len(daily)}


def _report(res, tag):
    print(f"\n  === {tag} ===")
    print(f"  瞬間駆動のみ:  R²={res['r2_base']:.3f}  日残差自己相関={res['ac_base']:+.2f}")
    print(f"  ＋遅いプロキシ: R²={res['r2_slow']:.3f}  日残差自己相関={res['ac_slow']:+.2f}"
          f"（回復={res['r2_slow']-res['r2_base']:+.3f}）")
    drop = res['ac_base'] - res['ac_slow']
    if np.isfinite(drop):
        if res['ac_slow'] < 0.3:
            v = "✅ 遅いプロキシで自己相関がほぼ消えた＝遅い成分は『観測から作れる』（真に未観測でない）"
        elif drop >= 0.3:
            v = "○ 大きく下がったが残る＝一部は観測プロキシで説明、一部は未観測"
        else:
            v = "★ ほぼ落ちない＝観測プロキシで説明できない『真に未観測』の遅い駆動（衛星の出番）"
        print(f"  → 判定: {v}")


def make_synth(kind, days=200, per_day=48, seed=0):
    rng = np.random.default_rng(seed)
    n = days * per_day
    idx = pd.date_range("2001-06-01", periods=n, freq="30min")
    doy = idx.dayofyear.to_numpy(); hour = idx.hour.to_numpy() + idx.minute.to_numpy() / 60
    Rg = np.clip(np.sin((hour - 6) / 12 * np.pi), 0, None) * 800 + rng.normal(0, 30, n)
    Ta = 20 + 6 * np.sin((hour - 9) / 24 * 2 * np.pi) + 0.03 * (doy - doy.mean()) + rng.normal(0, 1, n)
    th = 0.3 + rng.normal(0, 0.02, n)
    VPD = np.clip(0.5 + 0.06 * (Ta - 15), 0.1, None)
    cumRg = pd.Series(Rg).rolling(WIN, min_periods=1).mean().to_numpy()
    base = np.exp(0.06 * (Ta - 20)) * (th / (0.2 + th))
    if kind == "obs_slow":      # 遅い成分＝累積Rg（観測から作れる基質）
        slow = 1.5 * _z(cumRg)
    elif kind == "hidden_slow": # 遅い成分＝観測と無関係な AR（真に未観測）
        hd = np.zeros(days)
        for d in range(1, days):
            hd[d] = 0.85 * hd[d - 1] + rng.normal(0, 0.3)
        slow = 1.5 * np.repeat(hd, per_day)
    GER = 2.0 * base + slow + rng.normal(0, 0.15, n)
    return pd.DataFrame({"GER": np.clip(GER, 1e-3, None), "Rg": Rg, "Ta": Ta, "VPD": VPD,
                         "Ts": Ta - 1, "th": th, "gH": Rg * 0.2, "gLE": Rg * 0.3,
                         "P": rng.normal(0, 0.5, n)}, index=idx)


def main():
    p = argparse.ArgumentParser(description="呼吸の遅い成分は観測から作れるか真に未観測か")
    p.add_argument("--site")
    p.add_argument("--years", type=int, nargs="+")
    p.add_argument("--month", type=int, nargs="+", default=[7, 8])
    p.add_argument("--qc-max", type=int, default=None,
                   help="QCフラグ閾値（例1＝実測寄り。既定None＝gap-fill込み）。"
                        "穴埋め由来かの検証用：--qc-max 1 で遅い未観測が生き残るか。")
    a = p.parse_args()

    if not a.site:
        print("=== 旗25 合成検証：遅い成分は観測プロキシで消えるか ===")
        _report(analyze(add_slow(make_synth("obs_slow"))), "遅い成分＝累積Rg（観測から作れる）")
        _report(analyze(add_slow(make_synth("hidden_slow"))), "遅い成分＝観測無関係AR（真に未観測）")
        print("\n  → 期待：観測性の遅い成分は＋プロキシで自己相関が消え、真に未観測は残る。")
        return

    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    cfg = AnalysisConfig(qc_max=a.qc_max) if a.qc_max is not None else AnalysisConfig()
    raw_all = load_raw_all(get_site(a.site), cfg)
    qtag = f"・QC≤{a.qc_max}(実測寄り)" if a.qc_max is not None else "・gap-fill込み(既定)"
    ms = sorted(a.month)
    frames, used = [], []
    for y in a.years or []:
        start = pd.Timestamp(year=y, month=ms[0], day=1)
        end = pd.Timestamp(year=y, month=ms[-1], day=1) + pd.offsets.MonthBegin(1)
        r = raw_all[(raw_all.index >= start) & (raw_all.index < end)]
        keep = ["GER"] + [v for v in DRIVERS if v in r.columns]
        r = r[keep]
        if len(r.dropna()) > 200:
            frames.append(add_slow(r)); used.append(y)   # 遅いプロキシは年内で計算
    if not frames:
        print("有効年なし"); return
    df = pd.concat(frames)
    print(f"=== 旗25 実データ {a.site}（GER 遅い成分は観測から作れるか, {len(used)}年{qtag}）===")
    _report(analyze(df), f"{a.site} 呼吸 GER")
    print(f"  ※データ点数 N={len(df.dropna())}（QC≤{a.qc_max} で減れば穴埋め依存の目安）"
          if a.qc_max is not None else f"  ※データ点数 N={len(df.dropna())}（gap-fill込み）")
    print("\n  読み方：＋遅いプロキシ(累積Rg=基質/累積Ta=熱履歴/通日=フェノロジー)で日残差自己相関が")
    print("     消える→遅い成分は観測から作れる(モデルに遅い項を足せば済む)＝真に未観測でない。")
    print("     残る→観測プロキシで説明できない真に未観測の遅い駆動＝衛星プロキシ検証の入口(北極星)。")


if __name__ == "__main__":
    main()
