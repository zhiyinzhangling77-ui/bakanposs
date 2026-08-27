"""旗73：**中核主張の検出器の偽陽性率を測る**——閾値 ACF1≥0.64 は何%を許すのか。

本研究で最も強い主張は **A-1「チャンバー呼吸に多日メモリがある（森林 22/45）」**である。
判定は旗53/54 の較正済み検出器：**非線形基底で当てはめ→残差の ACF1≥0.64・e-fold≤7日・R²≥0.3**。

しかし **ACF1≥0.64 という閾値は「帰無が 0.49 を出したから 0.64 にした」**という決め方であり、
**それが何%の偽陽性を許すのかは測っていない**。旗72 で O-information の帰無に見つけたのと
**同じ形の穴**である（「帰無を作ったが、その校正を確かめていない」）。

もし偽陽性率が 15% なら、45 サイト中 **約7件は偶然**であり、**22 は実質 15 前後**になる。
**主張の大きさが変わる**ので、測る価値がある。

## 二つのやり方（**両方やる**）

**(A) 合成（コンテナ内で完結）**：季節＋自己相関のある駆動を作り、**メモリを一切入れず**に
Rs = 非線形関数(T, W) + **白色雑音** を生成する。＝**真のメモリはゼロ**。
検出器をかけ、**閾値ごとの偽陽性率**を測る。同じ設定でメモリを植えた場合の**検出力**も測る。

**(B) 実データに基づく帰無（サイトごと・より忠実）**：各サイトの**本物の Tsoil/SM をそのまま使い**、
実測 Rs に当てた**非線形基底の予測値＋白色雑音**（雑音の大きさは実測残差に合わせる）で
**メモリの無い Rs** を作る。＝**駆動の季節性・自己相関・欠測構造をすべて本物にしたまま**、
メモリだけを消せる。これを何度も作って **★が出る割合＝そのサイトの偽陽性率**を測る。

(B) が本命である。(A) は道具が正しく動くことの確認と、閾値の形を見るためのもの。

**注意**：雑音は**白色**でなければならない。実測残差の自己相関に合わせて雑音を作ると、
**測りたいメモリそのものを帰無に入れてしまう**（旗71 で同型の誤りをした）。

    python research/detector_fpr_step73.py --synth
    python research/detector_fpr_step73.py --cosore-dir /mnt/hdd/cosore-0.7.0 --igbp forest --nrep 40
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore, _acf_gap, _efold_gap
from memory_attribution_flex_step54 import flex_basis, _fit, ACF_THR, EFOLD_MAX

THRESHOLDS = (0.40, 0.50, 0.64, 0.70, 0.80)


def detect(y, T, W):
    """旗53/54 と同じ検出器（非線形基底→残差の ACF1・e-fold・R²）。"""
    with np.errstate(divide="ignore", invalid="ignore"):
        ly = np.log(np.where(y > 0, y, np.nan))
    res, r2 = _fit(ly, flex_basis(T, W))
    if res is None or not np.isfinite(r2):
        return None
    return {"r2": r2, "acf1": _acf_gap(res, 1), "efold": _efold_gap(res)}


def is_star(m, thr=ACF_THR):
    """★＝R²≥0.3 かつ ACF1≥thr かつ e-fold≤7日（旗53/54 の規則）。"""
    if m is None or not np.isfinite(m["r2"]) or m["r2"] < 0.3:
        return None                       # 駆動弱＝判定不能（★でも×でもない）
    if not (np.isfinite(m["acf1"]) and np.isfinite(m["efold"])):
        return None
    return bool(m["acf1"] >= thr and m["efold"] <= EFOLD_MAX)


def null_series(pred, resid_sd, rng):
    """**メモリの無い** Rs：当てはめ予測 ＋ **白色**雑音（対数空間）。"""
    return np.exp(pred + rng.normal(0, resid_sd, len(pred)))


# ---------- (A) 合成 -------------------------------------------------------------
def synth_drivers(n_days, rng):
    """季節＋自己相関のある地温・水分（**メモリは Rs に入れない**）。"""
    t = np.arange(n_days)
    seas = np.sin(2 * np.pi * t / 365.25)
    wT = np.zeros(n_days); wW = np.zeros(n_days)
    for i in range(1, n_days):
        wT[i] = 0.85 * wT[i - 1] + rng.normal(0, 1)     # 気象は自己相関する
        wW[i] = 0.90 * wW[i - 1] + rng.normal(0, 1)
    T = 12 + 10 * seas + 1.2 * wT
    W = np.clip(0.25 + 0.05 * seas + 0.02 * wW, 0.03, 0.6)
    return T, W


def synth_case(n_days, rng, plant_efold=None):
    """Rs を作る。``plant_efold`` が None なら**メモリ無し**（＝帰無）。"""
    T, W = synth_drivers(n_days, rng)
    base = 0.6 + 0.08 * (T - 12) + 2.0 * (W - 0.25) - 3.0 * (W - 0.25) ** 2
    ly = np.log(np.clip(base, 0.05, None))
    if plant_efold:
        phi = np.exp(-1.0 / plant_efold)
        h = np.zeros(n_days)
        for i in range(1, n_days):
            h[i] = phi * h[i - 1] + rng.normal(0, 1)
        ly = ly + 0.25 * h / max(np.std(h), 1e-9)
    ly = ly + rng.normal(0, 0.12, n_days)               # **白色**雑音
    return np.exp(ly), T, W


def run_synth(nrep, n_days, seed=0):
    print("  ── (A) 合成：**メモリを入れない**系列で、閾値ごとの偽陽性率を測る ──")
    print(f"     日数 {n_days}・反復 {nrep}。駆動（地温・水分）は季節＋自己相関つき。")
    for lab, ef in [("**帰無（メモリ無し）**", None), ("参考：メモリを植えた場合（e-fold 4日）", 4.0)]:
        stars = {t: 0 for t in THRESHOLDS}
        acfs, judged = [], 0
        for r in range(nrep):
            rng = np.random.default_rng(seed + r)
            y, T, W = synth_case(n_days, rng, ef)
            m = detect(y, T, W)
            if m is None or not np.isfinite(m["r2"]) or m["r2"] < 0.3:
                continue
            judged += 1
            acfs.append(m["acf1"])
            for t in THRESHOLDS:
                if is_star(m, t):
                    stars[t] += 1
        print(f"     {lab}（判定できた {judged}/{nrep}・ACF1 の中央 "
              f"{np.median(acfs):+.2f}）" if judged else f"     {lab}：判定できず")
        if judged:
            row = "  ".join(f"{t:.2f}→{stars[t]/judged:5.1%}" for t in THRESHOLDS)
            print(f"       閾値ごとの★率： {row}")
    print("     → **帰無の★率が 5% を大きく超える閾値は使ってはいけない**。")
    print("       採用中の 0.64 がそこに入っていないかを見る。\n")


# ---------- (B) 実データに基づく帰無 ------------------------------------------------
def run_real(cosore_dir, igbp, months, nrep, seed=0):
    root = Path(cosore_dir); desc = pd.read_csv(root / "description.csv")
    print(f"  ── (B) 実データに基づく帰無（{igbp or '全'}）──")
    print("     各サイトの**本物の Tsoil/SM をそのまま使い**、実測 Rs に当てた非線形基底の")
    print("     **予測値＋白色雑音**で**メモリの無い Rs** を作る＝駆動の構造は本物のまま。")
    print(f"     反復 {nrep} 回／サイト。★が出た割合＝**そのサイトの偽陽性率**。\n")
    print(f"  {'dataset':<30}{'N':>5}{'実測':>7}{'一律':>5}"
          f"{'帰無中央':>9}{'帰無95%':>9}{'偽陽性率':>9}{'サイト別':>9}")
    tot_fp, tot_site, real_star, cal_star, both = 0.0, 0, 0, 0, 0
    disagree = []
    for _, d in desc.iterrows():
        ds = str(d["CSR_DATASET"]); ig = str(d.get("CSR_IGBP", ""))
        if igbp and igbp.lower() not in ig.lower():
            continue
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            df, st, sm = load_cosore(f, months)
        except Exception:
            continue
        if df is None or "Tsoil" not in df or "Rs" not in df:
            continue
        cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
        daily = df[cols].groupby(df.index.normalize()).mean()
        daily = daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D"))
        y = daily["Rs"].to_numpy(); T = daily["Tsoil"].to_numpy()
        W = daily["SM"].to_numpy() if "SM" in daily else None
        m = detect(y, T, W)
        star = is_star(m)
        if star is None:
            continue
        # 当てはめ予測と残差 sd（**白色**雑音の大きさに使う）
        with np.errstate(divide="ignore", invalid="ignore"):
            ly = np.log(np.where(y > 0, y, np.nan))
        res, _ = _fit(ly, flex_basis(T, W))
        pred = ly - res
        sd = float(np.nanstd(res))
        if not np.isfinite(sd) or sd <= 0:
            continue
        nstar, njud, nacf = 0, 0, []
        for r in range(nrep):
            rng = np.random.default_rng(seed + 7919 * r + len(ds))
            yn = null_series(pred, sd, rng)
            mn = detect(yn, T, W)
            sn = is_star(mn)
            if sn is None:
                continue
            njud += 1
            nacf.append(mn["acf1"])
            nstar += int(sn)
        if njud == 0:
            continue
        fp = nstar / njud
        # **サイトごとの帰無で判定する**（一律の閾値ではなく、そのサイトの帰無の95%点を超えるか）
        q95 = float(np.percentile(nacf, 95))
        star_cal = bool(m["acf1"] >= q95 and m["efold"] <= EFOLD_MAX and m["r2"] >= 0.3)
        tot_fp += fp; tot_site += 1; real_star += int(star); cal_star += int(star_cal)
        both += int(star and star_cal)
        if star != star_cal:
            disagree.append((ds, star, star_cal, m["acf1"], q95))
        print(f"  {ds:<30}{int(np.isfinite(y).sum()):>5}{m['acf1']:>7.2f}"
              f"{'★' if star else '·':>5}{np.median(nacf):>9.2f}{q95:>9.2f}"
              f"{fp:>8.0%}{'★' if star_cal else '·':>9}")
    print(f"\n  === まとめ ===")
    if tot_site == 0:
        print("  判定できたサイトが無い"); return
    mean_fp = tot_fp / tot_site
    print(f"  判定できたサイト {tot_site}")
    print(f"  **一律の閾値（ACF1≥{ACF_THR}）での★：{real_star} 件**")
    print(f"  **偽陽性率の平均 {mean_fp:.1%}** ＝ 偶然の★は期待値 **{mean_fp*tot_site:.1f} 件**")
    print(f"  ＝ 一律の★ {real_star} 件のうち、**実質は約 {real_star - mean_fp*tot_site:.1f} 件**")
    print(f"  **サイトごとの帰無（各サイトの95%点）での★：{cal_star} 件**（両方で★：{both} 件）")
    if disagree:
        print(f"  **判定が食い違ったサイト {len(disagree)} 件**（一律／サイト別）：")
        for ds, a1, a2, ac, q in disagree[:12]:
            print(f"    {ds:<30} {'★' if a1 else '·'}／{'★' if a2 else '·'}"
                  f"  実測 {ac:+.2f} 対 そのサイトの95%点 {q:+.2f}")
        if len(disagree) > 12:
            print(f"    …他 {len(disagree)-12} 件")
    print("\n  読み方：")
    print("   ・偽陽性率が 5% 前後なら、**閾値 0.64 は妥当**であり A-1 の 22/45 はほぼそのまま。")
    print("   ・**15% を超えるなら A-1 の件数は過大**であり、**閾値を上げ直す必要がある**。")
    print("   ・サイトごとに偽陽性率が大きく違うなら、**日数・欠測・駆動の強さで決まっている**")
    print("     ＝**一律の閾値そのものが不適切**という別の問題を示す。")
    print("     合成（A）では駆動の自己相関を強くすると帰無 ACF1 の中央が +0.59 まで上がった")
    print("     （旗53 が報告した帰無 0.49 より高い）＝**帰無の水準は駆動の構造で変わる**。")
    print("     → **サイトごとの帰無の95%点で判定する方が原理的に正しい**。上表の右端がそれ。")
    print("  留保：")
    print("   ・帰無は『**非線形基底で表せる部分＋白色雑音**』である。基底で表せない")
    print("     **系統的な形（未知の非線形）**が実在すれば、それを『メモリ』と誤認する余地は残る。")
    print("     ＝本ツールが測るのは**雑音由来の偽陽性**であって、**モデル誤特定由来ではない**。")
    print("   ・雑音を**白色**にしたのは意図的である。実測残差の自己相関に合わせると")
    print("     **測りたいメモリを帰無に入れてしまう**（旗71 で同型の誤りをした）。")


def main():
    p = argparse.ArgumentParser(description="中核検出器の偽陽性率を測る")
    p.add_argument("--cosore-dir")
    p.add_argument("--igbp", default="forest")
    p.add_argument("--month", type=int, nargs="+", default=None)
    p.add_argument("--nrep", type=int, default=40)
    p.add_argument("--days", type=int, default=700)
    p.add_argument("--synth", action="store_true")
    a = p.parse_args()

    print("=== 旗73：中核主張（A-1 チャンバー多日メモリ）の検出器の偽陽性率 ===")
    print(f"  規則は旗53/54：R²≥0.3・**ACF1≥{ACF_THR}**・e-fold≤{EFOLD_MAX}日。")
    print("  この閾値は『帰無が 0.49 を出したから 0.64』と決めたもので、"
          "**許す偽陽性率は未測定**だった。\n")
    if a.synth or not a.cosore_dir:
        run_synth(a.nrep, a.days)
        if not a.cosore_dir:
            return
    run_real(a.cosore_dir, a.igbp, a.month, a.nrep)


if __name__ == "__main__":
    main()
