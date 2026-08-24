"""旗45：呼吸の~4日メモリの『正体』を仮説駆動で当てにいく（未知→候補の検証へ）。

旗40 で ~4日メモリは本物（チャンバー直接測定）と確定したが、正体は「未知の遅い土壌プロセス」で
止まっている。文献は具体的候補を挙げている：
  (a) **先行降水/湿潤パルス（Birch効果, Cable 2013 antecedent moisture）** — 乾いた土が濡れると
      微生物が一気に活性化し数日尾を引く。
  (b) **先行水分の積算**（同 Cable：草原は過去2週間、灌木林は10週間が効く）。
  (c) **深層土壌水分** — 表層θでは見えない遅い水。
  (d) **熱慣性/位相遅れ** — 深部土壌温度・表層温度のラグ（Stoy 2007 が落葉期にも同ラグ＝物理拡散を示唆）。

方法：旗40と同じ基本モデル（表層Tsoil+表層SM）の日残差にメモリ(ACF1・e-fold)が在るサイトだけを対象に、
候補ブロックを1つずつ追加して**残差メモリがどれだけ潰れるか**を見る。最も潰す候補＝メモリの正体の第一候補。

**過剰適合の防御**：候補と同次元の**プラセボ**（同じ統計・位相だけずらした系列）を必ず併走させる。
プラセボでも同程度メモリが落ちるなら、それは説明でなく自由度の産物＝「説明せず」と判定する。

**検出器の限界（合成で確認済み・正直に）**：水系候補どうし（湿潤パルス／先行水分／深層水分）は互いに
共線なので、「正体=深層水分」の合成でも僅差で「湿潤パルス」が勝った。＝本ツールが信頼できる判別は
**「水の履歴 か／熱慣性 か／どれでもない か」**の3択まで。水の内訳（表層パルス vs 深層 vs 積算）は鋭くない。

    python research/memory_attribution_step45.py                                   # 合成で検証
    python research/memory_attribution_step45.py --cosore-dir /mnt/hdd/cosore-0.7.0 # 実データ(森林)
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import _acf_gap, _efold_gap

WINDOWS = (3, 7, 14, 30)
LAGS = (1, 2, 4, 8)


# ---------- 読み込み：深さ違いの列と降水列まで拾う -------------------------------
def _depth_cols(cols, prefix):
    out = []
    for c in cols:
        m = re.fullmatch(rf"CSR_{prefix}(\d+\.?\d*)", c)
        if m:
            out.append((float(m.group(1)), c))
    return sorted(out)


def _precip_col(cols):
    for c in cols:
        if re.search(r"PRECIP|PPT|RAIN", c, re.I):
            return c
    return None


def load_daily(path, months=None):
    """日次・連続グリッドの Rs と、深さ別の土壌温度/水分・降水を返す。"""
    df = pd.read_csv(path)
    cols = list(df.columns)
    if "CSR_FLUX_CO2" not in cols:
        raise ValueError("CSR_FLUX_CO2 なし")
    tc = "CSR_TIMESTAMP_BEGIN" if "CSR_TIMESTAMP_BEGIN" in cols else "CSR_TIMESTAMP_END"
    ts = pd.to_datetime(df[tc], errors="coerce")
    d = pd.DataFrame({"Rs": pd.to_numeric(df["CSR_FLUX_CO2"], errors="coerce").to_numpy()}, index=ts)

    temps, sms = _depth_cols(cols, "T"), _depth_cols(cols, "SM")
    meta = {"T_depths": [z for z, _ in temps], "SM_depths": [z for z, _ in sms]}
    if temps:                                     # 浅い=5cmに最も近い / 深い=最深
        sh = min(temps, key=lambda x: abs(x[0] - 5))[1]
        d["T_sh"] = pd.to_numeric(df[sh], errors="coerce").to_numpy()
        if temps[-1][1] != sh:
            d["T_dp"] = pd.to_numeric(df[temps[-1][1]], errors="coerce").to_numpy()
    elif "CSR_TAIR" in cols:
        d["T_sh"] = pd.to_numeric(df["CSR_TAIR"], errors="coerce").to_numpy()
    if sms:
        sh = min(sms, key=lambda x: abs(x[0] - 5))[1]
        d["SM_sh"] = pd.to_numeric(df[sh], errors="coerce").to_numpy()
        if sms[-1][1] != sh:
            d["SM_dp"] = pd.to_numeric(df[sms[-1][1]], errors="coerce").to_numpy()
    pc = _precip_col(cols)
    if pc:
        d["P"] = pd.to_numeric(df[pc], errors="coerce").to_numpy()
    meta["precip"] = pc

    d = d[d.index.notna()]
    if months:
        d = d[d.index.month.isin(months)]
    if d.empty:
        raise ValueError("有効行なし")
    agg = {c: ("sum" if c == "P" else "mean") for c in d.columns}
    daily = d.groupby(d.index.normalize()).agg(agg)
    grid = pd.date_range(daily.index.min(), daily.index.max(), freq="D")
    return daily.reindex(grid), meta


# ---------- 候補ブロックの構築 ---------------------------------------------------
def _roll(s, w):
    return s.rolling(w, min_periods=max(2, w // 3)).mean()


def build_blocks(daily):
    """候補ブロック（名前 → 列DataFrame）。データに無い候補は作らない。"""
    b = {}
    sm = daily["SM_sh"] if "SM_sh" in daily else None

    # (a) 先行降水/湿潤パルス（Birch）：降水があれば積算、無ければ θ の"増加分"＝浸潤イベントで代用
    if "P" in daily:
        wet = daily["P"].fillna(0.0)
        lab = "先行降水(Birch)"
    elif sm is not None:
        wet = sm.diff().clip(lower=0).fillna(0.0)      # θの増加＝降雨浸潤の代理
        lab = "湿潤パルス(θ増加=降雨代理)"
    else:
        wet = None
    if wet is not None:
        b[lab] = pd.DataFrame({f"wet{w}": wet.rolling(w, min_periods=1).sum() for w in WINDOWS})

    # (b) 先行水分の積算（Cable 2013）
    if sm is not None:
        b["先行水分(積算θ)"] = pd.DataFrame({f"sm{w}": _roll(sm, w) for w in WINDOWS})

    # (c) 深層土壌水分
    if "SM_dp" in daily:
        b["深層水分"] = pd.DataFrame({"dp": daily["SM_dp"],
                                  **{f"dp{w}": _roll(daily["SM_dp"], w) for w in WINDOWS[:3]}})

    # (d) 熱慣性/位相遅れ：深部温度＋表層温度のラグ
    th = {}
    if "T_dp" in daily:
        th["Tdp"] = daily["T_dp"]
    if "T_sh" in daily:
        for L in LAGS:
            th[f"Tlag{L}"] = daily["T_sh"].shift(L)
    if th:
        b["熱慣性(深部T/ラグT)"] = pd.DataFrame(th)

    # プラセボ：同次元・同統計・位相だけずらした系列（過剰適合のベースライン）
    if sm is not None:
        n = len(sm); sh = max(60, n // 3)
        fake = pd.Series(np.roll(sm.to_numpy(), sh), index=sm.index)
        b["【プラセボ】位相ずらしθ"] = pd.DataFrame({f"pl{w}": _roll(fake, w) for w in WINDOWS})
    return b


# ---------- 当てはめとメモリ測定 -------------------------------------------------
def _fit_resid(y, X):
    """NaN行を落として最小二乗。残差を元のインデックス長に戻す（欠測はNaN）。"""
    M = np.column_stack([X.to_numpy(float), np.ones(len(X))])
    yy = y.to_numpy(float)
    ok = np.isfinite(yy) & np.isfinite(M).all(axis=1)
    if ok.sum() < max(60, 5 * M.shape[1]):
        return None, np.nan
    coef = np.linalg.lstsq(M[ok], yy[ok], rcond=None)[0]
    res = np.full(len(yy), np.nan)
    res[ok] = yy[ok] - M[ok] @ coef
    ss = np.sum((yy[ok] - yy[ok].mean()) ** 2)
    r2 = 1 - np.sum(res[ok] ** 2) / ss if ss > 0 else np.nan
    return res, float(r2)


def _mem(res):
    if res is None:
        return np.nan, np.nan
    return _acf_gap(res, 1), _efold_gap(res)


def analyze(daily, meta):
    if "T_sh" not in daily:
        return {"note": "土壌温度なし"}
    base_cols = [c for c in ("T_sh", "SM_sh") if c in daily]
    X0 = daily[base_cols]
    y = np.log(daily["Rs"].where(daily["Rs"] > 0))     # 呼吸は温度の指数関数＝ln で当てる
    res0, r2_0 = _fit_resid(y, X0)
    if res0 is None:
        return {"note": "点不足"}
    ac0, ef0 = _mem(res0)
    if not (np.isfinite(ac0) and np.isfinite(r2_0)) or r2_0 < 0.3 or ac0 <= 0.4 or ef0 > 7 or ef0 < 2:
        return {"note": f"対象外(基本R2={r2_0:.2f}, ACF1={ac0:.2f}, e-fold={ef0})"}

    out = {"r2_0": r2_0, "ac0": ac0, "ef0": ef0, "cands": {}}
    for name, blk in build_blocks(daily).items():
        X = pd.concat([X0, blk], axis=1)
        res, r2 = _fit_resid(y, X)
        ac, ef = _mem(res)
        out["cands"][name] = {"r2": r2, "ac": ac, "ef": ef,
                              "dac": ac0 - ac if np.isfinite(ac) else np.nan}
    return out


def verdict(res):
    """プラセボを超えて残差メモリを潰した候補があるか。"""
    if "note" in res:
        return "―" + res["note"], None
    cs = res["cands"]
    pl = next((v["dac"] for k, v in cs.items() if "プラセボ" in k), 0.0)
    pl = pl if np.isfinite(pl) else 0.0
    real = {k: v for k, v in cs.items() if "プラセボ" not in k and np.isfinite(v["dac"])}
    if not real:
        return "―候補なし", None
    best = max(real, key=lambda k: real[k]["dac"])
    d = real[best]["dac"]
    if d > pl + 0.15 and real[best]["ac"] < 0.4:
        return f"★{best}が記憶を説明", best
    if d > pl + 0.15:
        return f"○{best}が部分的に説明", best
    return "―どの候補も説明せず(未知のまま)", None


# ---------- 合成検証 --------------------------------------------------------------
def _synth(kind, days=900, seed=0):
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2015-01-01", periods=days, freq="D")
    doy = idx.dayofyear.to_numpy()
    T = 12 + 11 * np.sin(2 * np.pi * (doy - 100) / 365) + rng.normal(0, 1.2, days)
    rain = (rng.random(days) < 0.12) * rng.exponential(12, days)          # 降雨イベント
    sm = np.zeros(days); s = 0.25
    for i in range(days):                                                  # バケツ：降雨で上昇・指数減衰
        s = s + 0.004 * rain[i] - 0.02 * (s - 0.15)
        sm[i] = s
    lnR = 0.06 * T + 2.0 * sm
    if kind == "birch":                       # 湿潤パルスが数日尾を引く（Birch効果）
        pulse = np.zeros(days); p = 0.0
        for i in range(days):
            p = 0.72 * p + 0.05 * rain[i]     # e-fold ≈ 3日
            pulse[i] = p
        lnR += pulse
    elif kind == "deep":                      # 深層水分（遅い）が効く
        dp = pd.Series(sm).rolling(6, min_periods=1).mean().to_numpy()
        lnR += 2.5 * (dp - dp.mean())
    elif kind == "unknown":                   # どの観測候補でもない隠れAR
        h = np.zeros(days); v = 0.0
        for i in range(days):
            v = 0.72 * v + rng.normal(0, 0.08); h[i] = v
        lnR += h
    d = pd.DataFrame({"Rs": np.exp(lnR + rng.normal(0, 0.03, days)), "T_sh": T, "SM_sh": sm}, index=idx)
    if kind == "deep":
        d["SM_dp"] = pd.Series(sm, index=idx).rolling(6, min_periods=1).mean()
    return d, {}


def main():
    p = argparse.ArgumentParser(description="呼吸メモリの正体を候補で当てる")
    p.add_argument("--cosore-dir"); p.add_argument("--igbp", default="forest")
    p.add_argument("--month", type=int, nargs="+", default=None)
    a = p.parse_args()

    if not a.cosore_dir:
        print("=== 旗45 合成検証：仕込んだ正体を当てられるか（プラセボ併走）===")
        for kind, lab in [("birch", "正体=湿潤パルス(Birch)"), ("deep", "正体=深層水分"),
                          ("unknown", "正体=観測外の隠れAR(当ててはいけない)")]:
            d, m = _synth(kind)
            r = analyze(d, m)
            v, _ = verdict(r)
            print(f"  {lab}")
            if "note" in r:
                print(f"    {v}\n"); continue
            print(f"    基本: R2={r['r2_0']:.2f} ACF1={r['ac0']:.2f} e-fold={r['ef0']}日")
            for k, c in r["cands"].items():
                print(f"      {k:<22} ΔR2={c['r2']-r['r2_0']:+.3f}  ACF1→{c['ac']:+.2f} "
                      f"(Δ{c['dac']:+.2f})  e-fold→{c['ef']}日")
            print(f"    → {v}\n")
        print("  期待：Birch/深層は該当候補が★、隠れARは『どの候補も説明せず』。プラセボが勝つなら検出器不良。")
        return

    root = Path(a.cosore_dir); desc = pd.read_csv(root / "description.csv")
    print(f"=== 旗45 実データ：~4日メモリの正体を候補で当てる（{a.igbp}）===")
    tally = {}
    for _, dd in desc.iterrows():
        ds = str(dd["CSR_DATASET"]); igbp = str(dd.get("CSR_IGBP", ""))
        if a.igbp and a.igbp.lower() not in igbp.lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            daily, meta = load_daily(f, a.month)
            r = analyze(daily, meta)
        except Exception as e:
            continue
        if "note" in r:
            continue
        v, best = verdict(r)
        key = v.split("が")[0].lstrip("★○―")[:16] if best else "どの候補も説明せず"
        tally[key] = tally.get(key, 0) + 1
        print(f"\n  {ds}  (基本 R2={r['r2_0']:.2f} ACF1={r['ac0']:.2f} e-fold={r['ef0']}日"
              f" / 深さ T{meta.get('T_depths')} SM{meta.get('SM_depths')} 降水={meta.get('precip')})")
        for k, c in r["cands"].items():
            print(f"      {k:<22} ΔR2={c['r2']-r['r2_0']:+.3f}  ACF1→{c['ac']:+.2f} "
                  f"(Δ{c['dac']:+.2f})  e-fold→{c['ef']}日")
        print(f"    → {v}")
    print("\n  === まとめ（メモリを持つサイトの正体）===")
    if not tally:
        print("    対象サイトなし（基本R2≥0.3・ACF1>0.4・e-fold 2〜7日を満たすサイト）")
    for k, v in sorted(tally.items(), key=lambda x: -x[1]):
        print(f"    {k:<20} {v}")
    print("  読み：★/○＝その候補が残差メモリを潰した＝メモリの正体の第一候補。")
    print("        『どの候補も説明せず』が多数＝正体は観測されている水・熱の履歴では説明できない")
    print("        （＝基質供給や微生物動態など、チャンバー観測の外にある可能性が残る）。")
    print("  注：プラセボ(位相ずらしθ)を超える改善のみ採用＝自由度による見かけの改善を除いてある。")


if __name__ == "__main__":
    main()
