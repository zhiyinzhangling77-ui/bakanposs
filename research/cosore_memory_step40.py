"""旗40：チャンバー呼吸(COSORE)の4日記憶を直接測る＝分割アーティファクト説の最終決着。

段2まで：呼吸の4日記憶は(a)速い基質でない(真SIF棄却)、(b)一部は分割窓アーティファクト(旗39)。だが
「本物の生物記憶」と「分割の平滑さ」を完全に切れなかった——両方がGER(分割派生量)の中に同居するから。

**チャンバーRsは分割を通さない直接測定**(COSORE, Bond-Lamberty 2020)。同じ地点/同型のDBF森林で、
チャンバーRsの日残差に4日記憶があるかを旗37と同じ方法で測る：
  ・記憶あり(ACF高・e-fold~4日) → 4日記憶は**本物の生物物理**(分割アーティファクトでない, 未観測駆動は実在)。
  ・記憶なし(気象で説明) → 4日記憶は**分割の産物だった**。

COSORE data csv を --file で渡す。列: CSR_FLUX_CO2=Rs[µmol/m²/s]、CSR_TIMESTAMP_BEGIN、CSR_T<深さ>=土壌温度、
CSR_SM<深さ>=土壌水分(m³/m³)。土壌温度/水分は自動検出(--tcol/--smcol で上書き可)。

    python research/cosore_memory_step40.py                                  # 合成で検証
    python research/cosore_memory_step40.py --file /mnt/hdd/cosore-0.7.0/datasets/data_d20190504_SAVAGE_hf006-03.csv
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from memory_timescale_step37 import _z, _autocorr, efolding_days


def _pick_soil_temp(cols):
    """CSR_T<深さcm> のうち 5cm に最も近い層。無ければ CSR_TAIR。"""
    cands = []
    for c in cols:
        m = re.fullmatch(r"CSR_T(\d+\.?\d*)", c)
        if m:
            cands.append((abs(float(m.group(1)) - 5), c))
    if cands:
        return sorted(cands)[0][1]
    for c in ("CSR_TAIR", "CSR_TAIR_AMB"):
        if c in cols:
            return c
    return None


def _pick_sm(cols):
    cands = []
    for c in cols:
        m = re.fullmatch(r"CSR_SM(\d+\.?\d*)", c)
        if m:
            cands.append((abs(float(m.group(1)) - 5), c))
    return sorted(cands)[0][1] if cands else None


def load_cosore(path, months=None):
    df = pd.read_csv(path)
    cols = list(df.columns)
    if "CSR_FLUX_CO2" not in cols:
        raise ValueError(f"CSR_FLUX_CO2 が無い。列: {cols[:12]}")
    tcol = "CSR_TIMESTAMP_BEGIN" if "CSR_TIMESTAMP_BEGIN" in cols else "CSR_TIMESTAMP_END"
    ts = pd.to_datetime(df[tcol], errors="coerce")
    out = pd.DataFrame({"Rs": pd.to_numeric(df["CSR_FLUX_CO2"], errors="coerce")}, index=ts)
    st, sm = _pick_soil_temp(cols), _pick_sm(cols)
    if st:
        out["Tsoil"] = pd.to_numeric(df[st], errors="coerce").to_numpy()
    if sm:
        out["SM"] = pd.to_numeric(df[sm], errors="coerce").to_numpy()
    out = out[out.index.notna()]
    if months:
        out = out[out.index.month.isin(months)]
    return out.dropna(subset=["Rs"]), st, sm


def fit_resid(daily, use_drivers=True):
    """日次 Rs を土壌温度・水分(＋2次・交互作用)で回帰し残差。use_drivers=False は平均のみ。"""
    Y = daily["Rs"].to_numpy(float)
    cols = []
    if use_drivers:
        z = {}
        for v in ("Tsoil", "SM"):
            if v in daily:
                z[v] = _z(daily[v].to_numpy(float)); cols.append(z[v])
        if "Tsoil" in z:
            cols.append(z["Tsoil"] ** 2)
        if "Tsoil" in z and "SM" in z:
            cols.append(z["Tsoil"] * z["SM"])
    X = np.column_stack(cols + [np.ones(len(Y))]) if cols else np.ones((len(Y), 1))
    coef, *_ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ coef
    ss = np.sum((Y - Y.mean()) ** 2)
    r2 = 1 - np.sum(resid ** 2) / ss if ss > 0 else np.nan
    return r2, resid


def analyze(df):
    daily = df.groupby(df.index.normalize()).mean().dropna(subset=["Rs"])
    if len(daily) < 60:
        return {"note": f"日数不足({len(daily)})"}
    ac_raw = _autocorr(daily["Rs"].to_numpy())
    r2, res = fit_resid(daily, use_drivers=True)
    ac_res = _autocorr(res); ef, _ = efolding_days(res)
    return {"n_days": len(daily), "ac_raw": ac_raw, "r2": r2,
            "ac_res": ac_res, "ef": ef,
            "has_T": "Tsoil" in daily, "has_SM": "SM" in daily}


def make_synth(kind, days=800, seed=0):
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2005-05-01", periods=days, freq="D")
    doy = idx.dayofyear.to_numpy()
    Tsoil = 12 + 8 * np.sin(2 * np.pi * (doy - 100) / 365.25) + rng.normal(0, 1, days)
    SM = np.clip(0.3 + rng.normal(0, 0.05, days), 0.05, 0.6)
    base = 2.0 * np.exp(0.07 * (Tsoil - 12))              # 温度応答
    if kind == "bio_memory":     # 気象と独立な4-5日AR記憶（本物の生物記憶）
        S = np.zeros(days)
        for i in range(1, days):
            S[i] = 0.8 * S[i - 1] + rng.normal(0, 0.5)
        Rs = base + 1.2 * _z(S) + rng.normal(0, 0.1, days)
    else:                         # 記憶なし（温度・水分＋白色雑音）
        Rs = base * (0.7 + 0.6 * (SM - 0.3) / 0.2) + rng.normal(0, 0.15, days)
    return pd.DataFrame({"Rs": np.clip(Rs, 1e-3, None), "Tsoil": Tsoil, "SM": SM}, index=idx)


def _report(r, tag, st=None, sm=None):
    if "note" in r:
        print(f"  {tag}: {r['note']}"); return
    drv = f"土壌温度={st or '—'} 水分={sm or '—'}"
    print(f"\n  === {tag}（{r['n_days']}日, {drv}）===")
    print(f"  生Rs 日次ACF        = {r['ac_raw']:+.2f}")
    print(f"  ドライバー回帰 R²   = {r['r2']:.3f}")
    print(f"  日残差ACF           = {r['ac_res']:+.2f}   e-fold = {r['ef']}日")
    if r["ac_res"] > 0.5 and r["ef"] >= 3:
        print("  → ★ チャンバーRsの残差に数日スケールの記憶が残る＝**分割を通さない直接測定でも4日記憶は本物**")
        print("     ＝呼吸の未観測の遅い駆動は分割アーティファクトでなく生物物理（弧の最深の裏づけ）")
    elif r["ac_res"] < 0.3:
        print("  → ・ 残差の記憶がほぼ無い＝チャンバーRsは土壌温度・水分でほぼ説明＝")
        print("     フラックスGERの4日記憶は分割の産物だった可能性（直接測定には無い）")
    else:
        print("  → △ 中間（記憶はあるが弱い）")


def main():
    p = argparse.ArgumentParser(description="チャンバー呼吸(COSORE)の4日記憶を直接測る")
    p.add_argument("--file")
    p.add_argument("--month", type=int, nargs="+", default=None, help="対象月(既定=全月)")
    a = p.parse_args()

    if not a.file:
        print("=== 旗40 合成検証：チャンバーRs残差の記憶を検出できるか ===")
        _report(analyze(make_synth("bio_memory")), "生物記憶あり（気象と独立な4-5日AR）")
        _report(analyze(make_synth("no_memory")), "記憶なし（温度・水分＋白色雑音）")
        print("\n  → 生物記憶ありは残差ACF高・e-fold~4、記憶なしは残差ACF≈0 が期待。")
        return

    df, st, sm = load_cosore(a.file, a.month)
    mtag = f"・月={a.month}" if a.month else "・全月"
    print(f"=== 旗40 実データ COSORE チャンバーRs（{Path(a.file).name}{mtag}, N={len(df)}）===")
    _report(analyze(df), Path(a.file).stem, st, sm)
    print("\n  意味：チャンバーRsは渦相関の分割を通さない直接測定。ここに4日記憶があれば＝GERの未観測駆動は")
    print("    分割アーティファクトでなく生物物理で本物。同一/同型DBF森林で我々のJP-Tak GER(旗37)と直接比較。")
    print("  留保：チャンバーは点測定(空間スケール小)・土壌呼吸のみ(地上部呼吸を含まず)・サイト固有。")


if __name__ == "__main__":
    main()
