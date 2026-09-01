"""旗100：**~4 日メモリは、我々が駆動に入れていなかった気象か**（事前登録 step100）。

旗99 で「**地点スケールの現象**」と分かった。**その中身の候補**は
**①我々が入れていない気象 ②林分規模の生物過程 ③土壌水文**。
**①は、同じ場所のタワーの気象を駆動に足せば検定できる。**

**旗45 は候補を残差に当てたが、それはすべて「チャンバー自身のセンサから作った量」**
（先行湿潤・積算水分・深層水分・Birch・熱慣性）。
**同じ場所のタワーが測っている気象（Rg・Ta・VPD・P）を駆動に入れたことは一度も無い。**

**事前登録 step100 で固定済み**：
  ・**対象は規則で選ぶ**——**チャンバーが★短メモリ**かつ**タワーと重なる日数 ≥60・暦年 ≥3**
    （**メモリが無い所で「消えた」と言っても意味がない**＝旗95 の対照の教訓）
  ・**統計量はラグ1の自己共分散**（**ACF1 だけでは記憶量を測れない**＝旗61）
  ・**「プラセボ 5 通りの最大値を上回る」かつ「減少率 20% 以上」の両方**（旗75 の `FLOOR`）
  ・**プラセボはタワーの気象を年ごと入れ替えた系列**（旗45/71）

    python research/tower_weather_step100.py                 # 合成で検証（既定）
    python research/tower_weather_step100.py --real --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore
from model_richness_step74 import design, residuals, measure, star
from same_site_arc_step66 import PAIRS, SENSITIVITY_ONLY

MIN_DAYS, MIN_YEARS = 60, 3
LAGS = (0, 1, 2, 3, 4, 5)                 # タワー気象のラグ（事前登録で固定）
WVARS = ("Rg", "Ta", "VPD", "P")          # タワーが測る気象（事前登録で固定）
N_PLACEBO = 5                             # 年ごと入れ替えの通し数
FLOOR = 0.20                              # 絶対的な下限（旗75）


def lag1_cov(x):
    """**ラグ1の自己共分散**＝記憶量（旗61：ACF1 だけでは測れない）。欠測跨ぎは除く。"""
    x = np.asarray(x, float)
    mu = np.nanmean(x)
    a, b = x[:-1] - mu, x[1:] - mu
    m = np.isfinite(a) & np.isfinite(b)
    return float(np.mean(a[m] * b[m])) if m.sum() >= 20 else np.nan


def tower_weather(site):
    """**タワーの気象だけ**を日平均で読む（`gLE`/`gH` の欠測で日を落とさないため）。"""
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    raw = load_raw_all(get_site(site), AnalysisConfig())
    keep = [c for c in WVARS if c in raw.columns]
    if len(keep) < 2:
        return None
    d = raw[keep].groupby(raw.index.normalize()).mean()
    return d.dropna(how="all")


def weather_design(w):
    """タワー気象と**そのラグ**の計画行列（**線形＋二乗**）。"""
    cols = []
    for c in w.columns:
        s = w[c]
        for L in LAGS:
            cols.append(s.shift(L).to_numpy())
        cols.append((s - s.mean()).to_numpy() ** 2)      # 同日の二乗だけ足す
    X = np.column_stack(cols + [np.ones(len(w))])
    return X


def reduction(res0, w):
    """`res0`（基準の残差）を**タワー気象**で説明し直し、**記憶量の減少率**を返す。"""
    X = weather_design(w)
    ok = np.isfinite(X).all(axis=1) & np.isfinite(res0)
    if ok.sum() < MIN_DAYS:
        return np.nan
    r1 = residuals(X[ok], res0[ok], True)
    if r1 is None:
        return np.nan
    c0, c1 = lag1_cov(res0[ok]), lag1_cov(r1)
    if not (np.isfinite(c0) and np.isfinite(c1)) or c0 <= 0:
        return np.nan
    return float(1 - c1 / c0)


def run_pair(site, ds, ch_daily, tw, verbose=True):
    """1 組。**★の確認 → 基準残差 → 気象を足す → プラセボ**。"""
    common = ch_daily.index.intersection(tw.index)
    yrs = pd.Index(common).year.nunique() if len(common) else 0
    if verbose:
        print(f"    重なり {len(common)} 日／{yrs} 暦年")
    if len(common) < MIN_DAYS or yrs < MIN_YEARS:
        return None, "重なり不足"
    ch = ch_daily.reindex(common)
    w = tw.reindex(common)

    # ── チャンバーが★短メモリか（旗74 と同一の物差し）──
    y = ch["Rs"].to_numpy(); T = ch["Tsoil"].to_numpy()
    W = ch["SM"].to_numpy() if "SM" in ch else None
    m = measure(y, T, W, "テンソルビン", True)
    st = star(m)
    if verbose and m:
        print(f"    チャンバー：R²={m['r2']:.2f}／ACF1={m['acf1']:+.2f}／"
              f"e-fold={m['efold']:.0f}日 → {'**★短メモリ**' if st else '★でない'}")
    if not st:
        return None, "チャンバーが★でない"

    # ── 基準の残差（{Ts, SM} のテンソルビン＋外挿）──
    X0 = design("テンソルビン", T, W)
    with np.errstate(divide="ignore", invalid="ignore"):
        ly = np.log(np.where(y > 0, y, np.nan))
    res0 = residuals(X0, ly, True)
    if res0 is None:
        return None, "基準残差を作れない"

    r_obs = reduction(res0, w)
    # ── プラセボ：**タワーの気象を年ごと入れ替える** ──
    pls = []
    for k in range(1, N_PLACEBO + 1):
        ws = tw.copy()
        ws.index = ws.index + pd.Timedelta(days=365 * k)
        wk = ws.reindex(common)
        if wk.notna().any(axis=1).sum() < MIN_DAYS:
            continue
        v = reduction(res0, wk)
        if np.isfinite(v):
            pls.append(v)
    if not np.isfinite(r_obs) or len(pls) < 3:
        return None, "減少率を出せない"
    pl_max = float(np.max(pls))
    won = bool(r_obs > pl_max and r_obs >= FLOOR)
    if verbose:
        print(f"    記憶量の減少率：**実測 {r_obs:+.0%}**／"
              f"プラセボ最大 {pl_max:+.0%}（{len(pls)} 通り・中央 {np.median(pls):+.0%}）")
        print(f"      → {'**説明された**（プラセボ超え かつ 20% 以上）' if won else '**説明されない**'}"
              f"{'' if r_obs >= FLOOR else '（**下限 20% に届かない**）'}")
    return won, None


def load_chambers(cosore_dir):
    root = Path(cosore_dir)
    out = {}
    for site, ds, km in PAIRS:
        if (site, ds) in SENSITIVITY_ONLY:
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            df, st, sm = load_cosore(f, None)
        except Exception:
            continue
        if df is None or "Tsoil" not in df or "Rs" not in df:
            continue
        cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
        d = df[cols].groupby(df.index.normalize()).mean()
        if len(d) >= MIN_DAYS:
            out.setdefault(site, []).append((ds, d))
    return out


def synth(kind, years=6, seed=0):
    """**タワー気象で作った駆動／タワーが測っていない駆動／駆動なし**。"""
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2010-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy()
    Rg = np.clip(150 + 120 * np.sin(2 * np.pi * (doy - 80) / 365)
                 + rng.normal(0, 30, len(idx)), 5, None)
    Ta = 12 + 12 * np.sin(2 * np.pi * (doy - 100) / 365) + rng.normal(0, 2.5, len(idx))
    VPD = np.clip(0.4 + 0.05 * Ta + rng.normal(0, 0.25, len(idx)), 0.02, None)
    P = np.where(rng.random(len(idx)) < 0.2, rng.exponential(4, len(idx)), 0.0)
    tw = pd.DataFrame({"Rg": Rg, "Ta": Ta, "VPD": VPD, "P": P}, index=idx)
    T = 15 + 10 * np.sin(2 * np.pi * (doy - 110) / 365) + rng.normal(0, 1.5, len(idx))
    W = np.clip(0.2 + 0.04 * np.sin(2 * np.pi * (doy - 200) / 365)
                + rng.normal(0, .02, len(idx)), .02, .6)
    if kind == "weather":                    # **タワー気象のラグで作った隠れ駆動**
        hid = (pd.Series(P, index=idx).rolling(4, min_periods=1).sum().to_numpy() / 4
               + 0.004 * pd.Series(Rg, index=idx).shift(2).fillna(Rg.mean()).to_numpy())
        hid = (hid - np.nanmean(hid)) / np.nanstd(hid)
    elif kind == "hidden":                   # **タワーが測っていない駆動**（移動平均）
        hid = np.convolve(rng.normal(0, 1, len(idx)), np.ones(6) / 6, "same")
    else:                                    # hidden_ar：**AR(1) の隠れ駆動**（別の作り方）
        hid = np.zeros(len(idx))
        e = rng.normal(0, 1, len(idx))
        for i in range(1, len(idx)):
            hid[i] = 0.78 * hid[i - 1] + e[i]
    hid = (hid - np.nanmean(hid)) / (np.nanstd(hid) + 1e-12)
    # **否定側の合成も★にする**——**★でなければ「説明されない」枝を通せない**
    # （第1版は振幅が足りず、`hidden`/`noise` が★の関門で落ちていた）。
    amp = 0.30 if kind == "weather" else 0.55
    lrs = -1.0 + 0.06 * T + 2.0 * W + amp * hid + rng.normal(0, .05, len(idx))
    ch = pd.DataFrame({"Rs": np.exp(lrs), "Tsoil": T, "SM": W}, index=idx)
    return ch, tw


def main():
    ap = argparse.ArgumentParser(description="旗100：タワーの気象を駆動に足す")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--cosore-dir", default="/mnt/hdd/cosore-0.7.0")
    a = ap.parse_args()

    print("=== 旗100：~4 日メモリは、我々が駆動に入れていなかった気象か ===")
    print("  **旗45 が当てたのはチャンバー自身のセンサから作った量**——")
    print("  **同じ場所のタワーが測る気象を駆動に入れたことは一度も無い。**")
    print(f"  **統計量はラグ1の自己共分散**（旗61）／"
          f"**プラセボ超え かつ 減少率 {FLOOR:.0%} 以上**で「説明された」（旗75）。")

    if not a.real:
        print("\n  【合成データで道具を検証する】**実データに触れる前に**、")
        print("  **`hidden`・`noise` で「説明された」と出ないか**を必ず見る（過剰適合の試験）。")
        for kind, want in (
                ("weather", "**説明されるべき**"),
                ("hidden", "**説明されてはいけない**（移動平均の隠れ駆動・過剰適合の試験）"),
                ("hidden_ar", "**説明されてはいけない**（AR(1) の隠れ駆動・別の作り方）")):
            print(f"\n  ===== 合成 `{kind}` —— 期待：{want} =====")
            ch, tw = synth(kind)
            v, why = run_pair("合成", kind, ch, tw)
            print(f"  【判定】{v if v is not None else '判定しない（' + str(why) + '）'}")
        print("\n  → **weather→True・hidden→False・hidden_ar→False** なら道具は使える。")
        return

    chambers = load_chambers(a.cosore_dir)
    print(f"\n  読めたチャンバー：{sum(len(v) for v in chambers.values())} 本"
          f"／タワー {len(chambers)} サイト")
    res, why = {}, {}
    for site in sorted(chambers):
        try:
            tw = tower_weather(site)
        except Exception as e:
            print(f"\n  ━━ {site} ━━\n    タワーを読めない {type(e).__name__}: {str(e)[:90]}")
            continue
        if tw is None or tw.empty:
            print(f"\n  ━━ {site} ━━\n    タワーの気象が無い"); continue
        print(f"\n  ━━ {site} ━━（気象 {list(tw.columns)}／{len(tw):,} 日）")
        for ds, ch in chambers[site]:
            print(f"    ── {ds} ──")
            v, w = run_pair(site, ds, ch, tw)
            key = f"{site}／{ds}"
            if v is None:
                why[key] = w
                print(f"      → **判定しない**（{w}）")
            else:
                res[key] = v

    print("\n  === 集計（事前登録の判定規則に当てる）===")
    for k, v in res.items():
        print(f"    {k:<48}{'**説明された**' if v else '説明されない'}")
    for k, w in why.items():
        print(f"    {k:<48}判定しない（{w}）")
    n = len(res)
    print("\n  === 結論 ===")
    if n < 3:
        print(f"  **判定しない**——条件を満たす組が {n} で 3 未満。")
        print("  **どの条件で落ちたか**は上の各行に書いてある。")
    elif sum(res.values()) > n / 2:
        print("  **★~4 日メモリは、我々が入れていなかった気象だった**。")
        print("  ＝**「観測の隙間」は埋まる**。**この研究の題名そのものが変わる。**")
    elif sum(res.values()) <= 1:
        print("  **▲気象では説明されない**——**地点規模の生物過程か土壌水文**に絞られる。")
        print("  ＝**新規観測へ渡す**（`NEW_OBSERVATION_DESIGN.md`）。")
    else:
        print("  **○一部で説明される**——**組ごとに記し、まとめない**。")
    print("\n  留保（事前登録どおり）：")
    print("   ・**タワーとチャンバーは同じ場所ではない**（0.00–0.27 km）。**フェッチと視野が違う。**")
    print("   ・**タワーの気象も穴埋め済み**（旗46）。**駆動側に入るので向きが違う**——そう読む。")
    print("   ・**「気象で説明された」と出ても、どの気象かは言えない**（まとめて足すため）。")
    print("   ・**説明されなくても気象一般が外れたとは言えない**——")
    print("     **タワーが測る 4 変数とその 1–5 日ラグの範囲**である。")


if __name__ == "__main__":
    main()
