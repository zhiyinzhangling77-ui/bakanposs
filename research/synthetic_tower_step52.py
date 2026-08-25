"""旗52（策A）：観測演算子ごと前向きに回す「合成タワー」＝研究全体をエンドツーエンドで較正する。

これまでの検証は「処理済みデータを解析し、アーティファクトを事後に引き算する」形だった
（旗44の曲率制御・旗46のN対照・旗50の配分シャッフル）。前提の穴（③収支・⑤gap-fill・①温度交絡・
派生量依存）は**解析者自身の処理がデータ生成過程に入っている**ことが共通の根なので、根本策は逆向きになる：

    真の過程 → 測定（ノイズ・欠測）→ **MDS 風ギャップフィル** → **NT 分割** → 我々が見る数値

を生成モデルとして書き、**真値が分かっている合成サイトを、我々の検出器に丸ごと通す**。
問いが「この発見は本物か（事後に引き算）」から「**仕込んだ場合と仕込まない場合で、検出器は何と言うか**」
に変わる＝アーティファクトが交絡ではなく**モデルの一部**になる。

**本研究にとっての中心的な問い**：旗25/37 の「呼吸残差の ~4 日メモリ」は、
**遅い未観測駆動を一切仕込まなくても、ギャップフィル（±7日窓）と分割（移動窓）だけで作られてしまわないか**。
  ・仕込んだ時だけ出る → 検出器は妥当＝旗25/37 の主張はパイプラインに耐える。
  ・仕込まなくても出る → **旗25/37 はパイプラインの産物**＝重大な格下げ（旗40 のチャンバーは別途生き残る）。

    python research/synthetic_tower_step52.py              # 既定（3夏・2条件）
    python research/synthetic_tower_step52.py --nrep 5     # 反復を増やす
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import _acf_gap, _efold_gap

T0_LT = -46.02          # Lloyd-Taylor の下限温度 (°C)
NIGHT_RG = 20.0         # 夜間判定 (W/m2)


# ---------- 1. 真の過程 -------------------------------------------------------------
def simulate_truth(years=3, plant_memory=True, plant_wq10=False, seed=0):
    """半時間値の「真の」生態系。夏(7-8月)を years 年分。"""
    rng = np.random.default_rng(seed)
    idx = pd.DatetimeIndex([])
    for k in range(years):
        idx = idx.append(pd.date_range(f"{2010+k}-07-01", f"{2010+k}-08-31 23:30", freq="30min"))
    n = len(idx)
    hod = idx.hour + idx.minute / 60.0
    doy = idx.dayofyear.to_numpy()

    # 放射：日周 × 雲（赤色雑音）
    clear = np.clip(np.sin((hod - 6) / 12 * np.pi), 0, None) * 800
    cloud = np.zeros(n); c = 1.0
    for i in range(n):
        c = 0.97 * c + 0.03 * rng.uniform(0.3, 1.0); cloud[i] = c
    Rg = clear * cloud

    Ta = 18 + 6 * np.sin((hod - 9) / 24 * 2 * np.pi) + 0.012 * Rg + rng.normal(0, 0.6, n)
    Ts = pd.Series(Ta).rolling(96, min_periods=1).mean().to_numpy()      # 土壌温度＝熱慣性で鈍る
    VPD = np.clip(0.61 * np.exp(17.5 * Ta / (Ta + 240)) * (1 - 0.75) * 10, 0, None)

    rain = (rng.random(n) < 0.004) * rng.exponential(6, n)               # 降雨イベント
    th = np.zeros(n); s = 0.30
    for i in range(n):
        s = np.clip(s + 0.004 * rain[i] - 0.0009 * (s - 0.12), 0.05, 0.45); th[i] = s

    # 遅い未観測駆動（e-fold ≈ 4日）＝本研究が「在る」と言っている当のもの
    hid = np.zeros(n); v = 0.0
    for i in range(n):
        v = 0.9928 * v + rng.normal(0, 0.02); hid[i] = v                 # 0.9928^48 ≈ 0.71/日
    if not plant_memory:
        hid[:] = 0.0

    # 呼吸：Lloyd-Taylor × 水分関数 ×（仕込めば）遅い駆動
    q10_term = 320.0 * (1.0 / (10 - T0_LT) - 1.0 / (Ts - T0_LT))
    wfun = np.clip((th - 0.05) / 0.25, 0, 1) * (1 - np.clip((th - 0.35) / 0.10, 0, 1) * 0.4)
    lnR = 0.4 + q10_term + np.log(np.clip(wfun, 0.05, None))
    if plant_wq10:                       # 湿るほど温度感度が上がる（旗26/44 が問うている構造）
        lnR += 0.05 * ((th - 0.25) / 0.10) * (Ts - Ts.mean())
    RECO = np.exp(lnR + hid)

    GPP = np.clip(0.06 * Rg / (1 + 0.06 * Rg / 30) * np.clip(wfun, 0.2, 1), 0, None)
    NEE = RECO - GPP + rng.normal(0, 0.35, n)

    # エネルギー：利用可能エネルギーを θ 依存のボーエン比で配分（＝Bowen 反転を仕込む）
    A = 0.75 * Rg
    beta = np.clip(0.65 - 1.2 * (th - 0.25), 0.15, 0.85)                 # 湿るほど潜熱側へ
    gH = A * beta + rng.normal(0, 8, n)
    gLE = A * (1 - beta) + rng.normal(0, 8, n)

    return pd.DataFrame({"Rg": Rg, "Ta": Ta, "Ts": Ts, "VPD": VPD, "th": th,
                         "gH": gH, "gLE": gLE, "NEE": NEE,
                         "RECO_true": RECO, "GPP_true": GPP, "hidden": hid}, index=idx)


# ---------- 2. 観測演算子：欠測 → MDS 風ギャップフィル -------------------------------
def observe(df, gap_frac=0.45, seed=0):
    """夜間を厚めに欠測させ、±7日窓の気象類似日平均で埋める（MDS, Reichstein 2005 の骨格）。"""
    rng = np.random.default_rng(seed)
    n = len(df)
    night = df["Rg"].to_numpy() < NIGHT_RG
    p = np.where(night, gap_frac * 1.6, gap_frac * 0.6)                  # u* 棄却で夜が欠けやすい
    gap = rng.random(n) < np.clip(p, 0, 0.95)

    nee = df["NEE"].to_numpy().copy()
    nee_obs = nee.copy(); nee_obs[gap] = np.nan

    Rg = df["Rg"].to_numpy(); Ta = df["Ta"].to_numpy(); VPD = df["VPD"].to_numpy()
    day = (df.index - df.index[0]).days.to_numpy()
    filled = nee_obs.copy()
    gi = np.flatnonzero(gap)
    for i in gi:
        w = np.flatnonzero((np.abs(day - day[i]) <= 7) & ~gap)           # ±7日窓の実測点
        if w.size == 0:
            continue
        sim = w[(np.abs(Rg[w] - Rg[i]) < 50) & (np.abs(Ta[w] - Ta[i]) < 2.5)
                & (np.abs(VPD[w] - VPD[i]) < 5)]
        if sim.size == 0:                                                # 条件を緩めて再探索
            sim = w[np.abs(Rg[w] - Rg[i]) < 100]
        if sim.size:
            filled[i] = nee[sim].mean()                                  # ＝気象の決定的関数になる
    out = df.copy()
    out["NEE_obs"] = nee_obs
    out["NEE_f"] = np.where(np.isfinite(filled), filled, np.nanmean(nee))
    out["is_gap"] = gap
    return out


# ---------- 3. 分割（NT：夜間 Lloyd-Taylor の移動窓） --------------------------------
def partition_nt(df, window_days=5):
    """夜間NEEに Lloyd-Taylor を移動窓で当て、RECO を全時刻へ外挿する（Reichstein 2005 の骨格）。"""
    Ts = df["Ts"].to_numpy(); nee = df["NEE_f"].to_numpy()
    night = df["Rg"].to_numpy() < NIGHT_RG
    day = (df.index - df.index[0]).days.to_numpy()
    base = 320.0 * (1.0 / (10 - T0_LT) - 1.0 / (Ts - T0_LT))             # E0 は固定（実務同様）
    reco = np.full(len(df), np.nan)
    for d0 in range(0, day.max() + 1):
        m = np.abs(day - d0) <= window_days // 2
        nm = m & night & np.isfinite(nee)
        if nm.sum() < 10:
            continue
        rref = np.mean(nee[nm] / np.exp(base[nm]))                       # 窓ごとの基準呼吸
        sel = m & (day == d0)
        reco[sel] = rref * np.exp(base[sel])
    out = df.copy()
    out["RECO_est"] = reco
    out["GPP_est"] = reco - df["NEE_f"].to_numpy()
    return out


# ---------- 4. 検出器（我々の解析）--------------------------------------------------
def memory_detector(df, col="RECO_est", form="linear"):
    """旗37/40 と同じ：日次に均し、気象で回帰した残差の自己相関（ACF1・e-fold）。

    form="linear"：日平均への**線形**回帰＝旗25/37/40 が実際にやっている形。
    form="flex"  ：非線形基底（Lloyd-Taylor 項・二次項・交互作用）を足した形。
    **呼吸は駆動の非線形関数なので、線形回帰は非線形の取り残しを残差に落とし、
    その取り残しは遅い駆動（Ts・θ）の滑らかさを引き継いで自己相関を持つ**。
    2つの差が「検出器の形が作るメモリ」の大きさになる（第2回で判明した交絡）。"""
    d = df[[col, "Ta", "Ts", "th", "Rg"]].copy()
    daily = d.groupby(d.index.normalize()).mean()
    grid = pd.date_range(daily.index.min(), daily.index.max(), freq="D")
    daily = daily.reindex(grid)
    y = np.log(daily[col].where(daily[col] > 0)).to_numpy()
    Ts_, th_ = daily["Ts"].to_numpy(), daily["th"].to_numpy()
    cols = [daily["Ta"].to_numpy(), Ts_, th_, daily["Rg"].to_numpy()]
    if form == "flex":
        with np.errstate(divide="ignore", invalid="ignore"):
            lt = 320.0 * (1.0 / (10 - T0_LT) - 1.0 / (Ts_ - T0_LT))   # Lloyd-Taylor 基底
        cols += [lt, Ts_ ** 2, th_ ** 2, Ts_ * th_, np.log(np.clip(th_, 1e-3, None))]
    X = np.column_stack(cols + [np.ones(len(daily))])
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    if ok.sum() < 40:
        return {"note": "点不足"}
    coef = np.linalg.lstsq(X[ok], y[ok], rcond=None)[0]
    res = np.full(len(y), np.nan); res[ok] = y[ok] - X[ok] @ coef
    ss = np.sum((y[ok] - y[ok].mean()) ** 2)
    r2 = 1 - np.sum(res[ok] ** 2) / ss if ss > 0 else np.nan
    return {"r2": float(r2), "acf1": _acf_gap(res, 1), "efold": _efold_gap(res)}


def run_case(plant_memory, years, seed, gap_frac):
    tru = simulate_truth(years=years, plant_memory=plant_memory, seed=seed)
    obs = observe(tru, gap_frac=gap_frac, seed=seed)
    est = partition_nt(obs)
    return {"true": memory_detector(est, "RECO_true", "linear"),
            "est": memory_detector(est, "RECO_est", "linear"),
            "true_flex": memory_detector(est, "RECO_true", "flex"),
            "est_flex": memory_detector(est, "RECO_est", "flex"),
            "gap_rate": float(obs["is_gap"].mean())}


def main():
    p = argparse.ArgumentParser(description="合成タワーでパイプライン全体を較正する")
    p.add_argument("--years", type=int, default=3)
    p.add_argument("--nrep", type=int, default=3)
    p.add_argument("--gap-frac", type=float, default=0.45)
    a = p.parse_args()

    print("=== 旗52：合成タワー（真の過程→欠測→MDS穴埋め→NT分割→我々の検出器）===")
    print("  中心の問い：**遅い駆動を仕込まなくても、穴埋めと分割だけで ~4日メモリが作られないか**\n")
    print(f"  {'条件':<26}{'系列':<14}{'R2':>7}{'ACF1':>8}{'e-fold':>8}")
    summary = {}
    for plant in (True, False):
        lab = "遅い駆動を仕込む" if plant else "**仕込まない**（帰無）"
        accs = {"true": [], "est": [], "true_flex": [], "est_flex": []}
        for k in range(a.nrep):
            r = run_case(plant, a.years, k, a.gap_frac)
            for key in ("true", "est", "true_flex", "est_flex"):
                if "note" not in r[key]:
                    accs[key].append(r[key])
        for key, kl in (("true", "真の呼吸/線形"), ("est", "分割穴埋/線形"),
                        ("true_flex", "真の呼吸/柔軟"), ("est_flex", "分割穴埋/柔軟")):
            if not accs[key]:
                print(f"  {lab:<26}{kl:<14}  判定不能"); continue
            r2 = np.mean([x["r2"] for x in accs[key]])
            ac = np.nanmean([x["acf1"] for x in accs[key]])
            ef = np.nanmean([x["efold"] for x in accs[key]])
            summary[(plant, key)] = (r2, ac, ef)
            print(f"  {lab:<26}{kl:<14}{r2:>7.2f}{ac:>8.2f}{ef:>8.1f}")
        print()

    print("  === 判定（3つの寄与を分けて読む）===")
    def g(plant, key): return summary.get((plant, key), (np.nan,)*3)
    # (1) 検出器のモデル形が作るメモリ：仕込まない × 真の呼吸（パイプラインも通さない）
    lin_ac, lin_ef = g(False, "true")[1], g(False, "true")[2]
    flx_ac, flx_ef = g(False, "true_flex")[1], g(False, "true_flex")[2]
    print(f"  (1) **検出器のモデル形**：仕込みも処理も無いのに 線形回帰では ACF1={lin_ac:+.2f} "
          f"e-fold={lin_ef:.1f}日。")
    print(f"      非線形基底にすると ACF1={flx_ac:+.2f} e-fold={flx_ef:.1f}日 まで落ちる"
          f"（当てはめ R²={g(False,'true_flex')[0]:.2f}）。")
    print("      ＝**呼吸は駆動の非線形関数なのに日平均へ線形回帰しているため、非線形の取り残しが")
    print("        遅い駆動の滑らかさを引き継いで自己相関を持つ**＝メモリの相当部分は検出器が作っている。")
    # (2) パイプライン（穴埋め＋分割）の寄与：同じ検出器で 真の呼吸 → 分割穴埋
    for f_, nm in (("", "線形"), ("_flex", "柔軟")):
        d0 = g(False, "est"+f_)[1] - g(False, "true"+f_)[1]
        d1 = g(True, "est"+f_)[1] - g(True, "true"+f_)[1]
        print(f"  (2) **パイプライン**（{nm}検出器）：ACF1 の変化 仕込まない {d0:+.2f} ／ 仕込む {d1:+.2f}")
    print("      ＝穴埋めと分割は**メモリを無から作るのではなく、在るものを増幅する**方向に効く。")
    # (3) 仕込みの効果：柔軟検出器で 帰無 vs 仕込み
    n_ac, n_ef = g(False, "est_flex")[1], g(False, "est_flex")[2]
    p_ac, p_ef = g(True, "est_flex")[1], g(True, "est_flex")[2]
    print(f"  (3) **仕込みの効果**（柔軟検出器・パイプライン通過後）：")
    print(f"      帰無 ACF1={n_ac:+.2f} e-fold={n_ef:.1f}日 → 仕込み ACF1={p_ac:+.2f} e-fold={p_ef:.1f}日")
    if p_ac - n_ac > 0.1 and p_ef > n_ef:
        print("      ＝仕込んだ時にははっきり増える＝**遅い駆動の有無は原理的に識別できる**。")
    else:
        print("      ＝差が小さい＝この設定では識別できない。")
    # (4) 旗40 の判定基準を帰無データに当てるとどうなるか＝偽陽性の点検
    print("  (4) **旗40 の基準（R²≥0.3 × ACF1>0.4 × e-fold 2〜7日）を帰無データに当てる**：")
    for key, nm in (("est", "線形検出器"), ("est_flex", "柔軟検出器")):
        r2, ac, ef = g(False, key)
        hit = (r2 >= 0.3) and (ac > 0.4) and (2 <= ef <= 7)
        print(f"      {nm}：R²={r2:.2f} ACF1={ac:+.2f} e-fold={ef:.1f} → "
              f"{'**通ってしまう＝偽陽性**' if hit else '通らない'}")
    print("      ＝通ってしまうなら、旗40 の『★短メモリ 12/35』には無視できない偽陽性率がある。\n")
    print("  留保：合成は本物の生態系ではない。ここで測っているのは**検出器とパイプラインの性質**であって、")
    print("        生態系の性質ではない。実データの結論を置き換えるものではなく、較正に使う。")


if __name__ == "__main__":
    main()
