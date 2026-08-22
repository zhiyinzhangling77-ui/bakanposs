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
    out = pd.DataFrame(
        {"Rs": pd.to_numeric(df["CSR_FLUX_CO2"], errors="coerce").to_numpy()}, index=ts)
    st, sm = _pick_soil_temp(cols), _pick_sm(cols)
    if st:
        out["Tsoil"] = pd.to_numeric(df[st], errors="coerce").to_numpy()
    if sm:
        out["SM"] = pd.to_numeric(df[sm], errors="coerce").to_numpy()
    out = out[out.index.notna()]
    if months:
        out = out[out.index.month.isin(months)]
    return out.dropna(subset=["Rs"]), st, sm


def _acf_gap(x, lag=1):
    """ギャップ(NaN)対応の lag 自己相関＝カレンダー上 lag 日離れた有限ペアだけで相関。"""
    x = np.asarray(x, float); mu = np.nanmean(x)
    a, b = x[:-lag] - mu, x[lag:] - mu
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 10 or a[m].std() == 0 or b[m].std() == 0:
        return np.nan
    return float(np.corrcoef(a[m], b[m])[0, 1])


def _efold_gap(x, maxlag=30):
    for k in range(1, maxlag + 1):
        r = _acf_gap(x, k)
        if np.isfinite(r) and r < 1 / np.e:
            return k
    return maxlag


def fit_resid(daily, use_drivers=True):
    """日次 Rs を土壌温度・水分(＋2次・交互作用)で回帰し残差。欠測は NaN で温存。
    ほぼ全欠測のドライバーは採らない。NaN 行はマスクして回帰。"""
    Y = daily["Rs"].to_numpy(float)
    cols, used = [], []
    if use_drivers:
        z = {}
        for v in ("Tsoil", "SM"):
            if v in daily:
                a = _z(daily[v].to_numpy(float))
                if np.isfinite(a).sum() > 0.5 * len(a):     # 半分以上 finite のみ採用
                    z[v] = a; cols.append(a); used.append(v)
        if "Tsoil" in z:
            cols.append(z["Tsoil"] ** 2)
        if "Tsoil" in z and "SM" in z:
            cols.append(z["Tsoil"] * z["SM"])
    X = np.column_stack(cols + [np.ones(len(Y))]) if cols else np.ones((len(Y), 1))
    mask = np.isfinite(Y) & np.all(np.isfinite(X), axis=1)
    resid = np.full(len(Y), np.nan)
    if mask.sum() < 30:
        return np.nan, resid, used
    coef, *_ = np.linalg.lstsq(X[mask], Y[mask], rcond=None)
    resid[mask] = Y[mask] - X[mask] @ coef
    ss = np.sum((Y[mask] - Y[mask].mean()) ** 2)
    r2 = 1 - np.sum(resid[mask] ** 2) / ss if ss > 0 else np.nan
    return r2, resid, used


def analyze(df):
    daily = df.groupby(df.index.normalize()).mean()
    # 連続日グリッドに張り直す（歯抜けを NaN 明示＝lag が真のカレンダー日になる）
    grid = pd.date_range(daily.index.min(), daily.index.max(), freq="D")
    daily = daily.reindex(grid)
    n_obs = int(daily["Rs"].notna().sum())
    if n_obs < 60:
        return {"note": f"日数不足({n_obs})"}
    ac_raw = _acf_gap(daily["Rs"].to_numpy())
    r2, res, used = fit_resid(daily, use_drivers=True)
    return {"n_days": n_obs, "n_grid": len(grid), "ac_raw": ac_raw, "r2": r2,
            "ac_res": _acf_gap(res), "ef": _efold_gap(res), "used": used,
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
    drv = f"使用ドライバー={r.get('used') or 'なし(平均のみ)'}（温度列={st or '—'} 水分列={sm or '—'}）"
    cov = f"{r['n_days']}日/{r.get('n_grid','?')}日グリッド"
    print(f"\n  === {tag}（{cov}, {drv}）===")
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


def run_batch(cosore_dir, months, igbp_filter="forest", min_days=90):
    """description.csv の全サイトを走査し、残差メモリを一括集計（穴①②を叩く：
    地理的に独立な多数の森林で「メモリは生物物理」が普遍かを直接測定で検証）。"""
    root = Path(cosore_dir)
    desc = pd.read_csv(root / "description.csv")
    rows = []
    print(f"=== 旗40 バッチ：COSORE 全サイトのチャンバーRs残差メモリ（{igbp_filter or '全'}, 月={months or '全'}）===")
    print(f"  {'dataset':<34}{'IGBP':<14}{'国/経度':>8} {'日数':>5} {"R²":>5} {"残差ACF":>7} {'e-fold':>6}  判定")
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); igbp = str(d.get("CSR_IGBP", ""))
        if igbp_filter and igbp_filter.lower() not in igbp.lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            df, st, sm = load_cosore(f, months)
            r = analyze(df)
        except Exception as e:
            print(f"  {ds:<34}{igbp[:12]:<14} SKIP {type(e).__name__}"); continue
        if "note" in r:
            print(f"  {ds:<34}{igbp[:12]:<14} {r['note']}"); continue
        # 厳しい分類（穴①再集計）：ドライバーが季節を除去(R²≥0.3)した上で、短e-fold(≤7日)の残差記憶＝真の~4日メモリ。
        # R²低=生の季節ACFを見ているだけ／e-fold長=季節(我々の発見でない)。
        if r["r2"] < 0.3:
            mem = "駆動弱"                       # 温度で季節を除けず＝残差≈生（判定不能）
        elif r["ac_res"] > 0.4 and r["ef"] <= 7:
            mem = "★短メモリ"                   # 我々の~4日メモリと同型（真）
        elif r["ac_res"] > 0.4:
            mem = "季節(長)"                     # e-fold長＝季節自己相関
        else:
            mem = "·なし"
        r.update({"dataset": ds, "igbp": igbp, "lon": d.get("CSR_LONGITUDE"), "mem": mem})
        rows.append(r)
        print(f"  {ds:<34}{igbp[:12]:<14}{str(d.get('CSR_LONGITUDE'))[:7]:>8} "
              f"{r['n_days']:>5} {r['r2']:>5.2f} {r['ac_res']:>+7.2f} {r['ef']:>5}日  {mem}")
    if not rows:
        print("  該当サイトなし"); return
    from collections import Counter
    c = Counter(r["mem"] for r in rows)
    judged = c["★短メモリ"] + c["季節(長)"] + c["·なし"]     # 駆動弱を除いた判定可能サイト
    ef_short = [r["ef"] for r in rows if r["mem"] == "★短メモリ"]
    print(f"\n  === まとめ（n={len(rows)}, 厳しい基準: R²≥0.3 で季節除去 × 残差ACF>0.4 × e-fold≤7）===")
    print(f"  ★短メモリ(真の~4日)：{c['★短メモリ']}／季節(長e-fold)：{c['季節(長)']}／なし：{c['·なし']}"
          f"／駆動弱(R²<0.3,判定不能)：{c['駆動弱']}")
    if judged:
        print(f"  判定可能({judged})サイト中 ★短メモリ：{c['★短メモリ']}/{judged}"
              f"（中央 e-fold={np.median(ef_short):.0f}日）" if ef_short else "")
    print("  読み方：★短メモリ＝温度で季節を除いた上で残る数日記憶＝我々の~4日メモリと同型（真の生物物理）。")
    print("    これが独立多サイトで多数なら擬似反復でない。季節(長)は季節自己相関＝別物。駆動弱はドライバー欠で判定不能。")


def main():
    p = argparse.ArgumentParser(description="チャンバー呼吸(COSORE)の4日記憶を直接測る")
    p.add_argument("--file")
    p.add_argument("--cosore-dir", help="COSORE ルート（--batch 一括集計）")
    p.add_argument("--igbp", default="forest", help="バッチのIGBPフィルタ（既定 forest, 空で全生態系）")
    p.add_argument("--month", type=int, nargs="+", default=None, help="対象月(既定=全月)")
    a = p.parse_args()

    if a.cosore_dir:
        run_batch(a.cosore_dir, a.month, igbp_filter=(a.igbp or None)); return

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
