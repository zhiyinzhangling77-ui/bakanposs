"""旗144：**骨格の `r` に、年をまたぐ対比はどれだけ混ざっているか**（事前登録 step144）。

**旗143 は `ES-FcO`（3 年）で、報告していた `r̂` の大半が「年の対比」だったことを示した**
（`−0.224 → −0.039`／`−0.283 → −0.199`）。**年の対比はたった 3 点の相関**である。

**同じ汚染が骨格（旗88/89/91/106/107）に乗っていないかを測る。**
これらは **15〜22 年**あるので効かない公算が高いが、**測っていない。**

**測るのは一つだけ**（`PREREGISTRATION_step144.md`）：

  ・`r_raw` ＝ 現行の偏 Spearman（Rg のみ除去）＝**これまで報告してきた値**
  ・`r_yr`  ＝ **Rg に加えて年ダミーも除去**した偏 Spearman（旗143 の `nullwhy` と同じ）
  ・`d = r_yr − r_raw` ＝ **年の対比が `r̂` に足していた分**
  ・**分解**：`r_raw = r_between + r_within`（恒等式・下の `decompose` で厳密に成り立つ）

**新しい検定はしない。CI も p も出さない**——**較正していない CI は出さない**（旗139）。

    python research/year_contrast_step144.py            # 門①（合成 3 本）＝既定
    python research/year_contrast_step144.py --real     # 実データ（/mnt/hdd）
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from moisture_control_atlas_step31 import partial_spearman, _rank
from stratified_bowen_step89 import cell_of, MIN_DAYS, MIN_YEARS, SEASON
from soiltemp_match_step90 import band, SPRING, AUTUMN
from rain_history_probe_step103 import (rain_history, daily_precip,
                                        PRIMARY_THR, REMOTE_MIN)
from evaporation_regime_step36 import daily_energy
from runlog import tee_stdout

SMALL, BIG = 0.05, 0.15                 # 事前登録で固定した 3 段階の境
S88, A88 = (4, 5), (9, 10)              # 旗88 の春・秋（旗89/91/107 は SPRING/AUTUMN）
NA3 = ("US-Wkg", "US-Whs", "US-SRM")
MN3 = ("MN-Hst", "MN-Nkh", "MN-Kbu")
REPS = 200                              # 門①の反復（事前登録で固定）


# ------------------------------------------------------------- 中核の測定
def _dummies(yr):
    """年ダミー（**最後の年を落とす**＝定数項と共線にしない）。"""
    u = np.unique(yr)
    return [(yr == y).astype(float) for y in u[:-1]]


def decompose(y, x, rg, yr):
    """`r_raw = r_between + r_within` を**厳密に**成り立たせて分解する。

    Rg の順位を除去した残差 `xr`・`yr_` を、**年平均成分**と**年内成分**に分ける。
    交差項は年内成分の年内平均が 0 なので**厳密に消える**——だから和は正確に一致する。
    """
    y = np.asarray(y, float); x = np.asarray(x, float); rg = np.asarray(rg, float)
    ok = np.isfinite(y) & np.isfinite(x) & np.isfinite(rg)
    if ok.sum() < 20:
        return None
    ry, rx, rz = _rank(y[ok]), _rank(x[ok]), _rank(rg[ok])
    Z = np.column_stack([rz, np.ones(ok.sum())])
    xr = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    yr_ = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    g = np.asarray(yr)[ok]
    xb = np.zeros_like(xr); yb = np.zeros_like(yr_)
    for u in np.unique(g):
        m = g == u
        xb[m] = xr[m].mean(); yb[m] = yr_[m].mean()
    xw, yw = xr - xb, yr_ - yb
    sx, sy = xr.std(), yr_.std()
    if sx == 0 or sy == 0:
        return None
    return (float(np.mean(xb * yb) / (sx * sy)),      # r_between
            float(np.mean(xw * yw) / (sx * sy)),      # r_within
            float(np.mean(xb ** 2) / (sx ** 2)))      # θ 側の年成分の分散比


def r_set(y, x, rg, yr):
    """`r_raw`・`r_yr`・分解を一度に返す。**判定はしない**（測るだけ）。"""
    y = np.asarray(y, float); x = np.asarray(x, float)
    rg = np.asarray(rg, float); yr = np.asarray(yr)
    r_raw, n = partial_spearman(y, x, [rg])
    r_yr, _ = partial_spearman(y, x, [rg] + _dummies(yr))
    dec = decompose(y, x, rg, yr)
    return {"raw": float(r_raw), "yr": float(r_yr), "n": int(n),
            "d": float(r_yr - r_raw) if np.isfinite(r_raw) and np.isfinite(r_yr) else np.nan,
            "between": dec[0] if dec else np.nan,
            "within": dec[1] if dec else np.nan,
            "vshare": dec[2] if dec else np.nan}


def klass(v):
    """事前登録の 3 段階。**追補 A 以降、当てるのは `r_between`**（`d` ではない）。

    **欠陥 #71**：`d = r_yr − r_raw` には**分母の再正規化**が混ざる
    ——**年の対比が偶然 0 でも、年成分を除けば分散が減って `|r|` は必ず大きくなる。**
    **`r_between` は `r_raw` と同じ分母で測るので、この影響を受けない。**
    """
    if not np.isfinite(v):
        return "?"
    return "無視" if abs(v) < SMALL else ("効く" if abs(v) < BIG else "骨格")


# ------------------------------------------------------------- 印字と集計
class Tally:
    def __init__(self):
        self.rows = []

    def add(self, arena, site, season, half, res, flip):
        self.rows.append(dict(arena=arena, site=site, season=season, half=half,
                              flip=flip, **res))

    def counts(self):
        c = {"無視": 0, "効く": 0, "骨格": 0, "?": 0}
        for r in self.rows:
            k = klass(r["between"])
            if r["flip"] and k != "?":
                k = "骨格"
            c[k] += 1
        return c


def measure_pair(tal, arena, site, season, sub):
    """1 つの部分集合（サイト×季節）について LE と H の両方を測って印字する。"""
    n, ny = len(sub), sub.index.year.nunique()
    if n < MIN_DAYS or ny < MIN_YEARS:
        print(f"      {season:<4} {n:>5} 日／{ny:>2} 年  **下限未満＝判定しない**")
        return None
    out = {}
    print(f"      {season:<4} {n:>5} 日／{ny:>2} 年（年あたり {n / ny:.0f} 日）")
    for k, col, nm in (("le", "gLE", "θ→γLE|Rg"), ("h", "gH", "θ→γH|Rg")):
        res = r_set(sub[col].to_numpy(), sub["th"].to_numpy(),
                    sub["Rg"].to_numpy(), sub.index.year.to_numpy())
        # **追補 A**：符号反転の判定は `r_within` と `r_raw` の間で見る（`r_yr` ではない）。
        flip = (np.isfinite(res["raw"]) and np.isfinite(res["within"])
                and np.sign(res["raw"]) != np.sign(res["within"]))
        tal.add(arena, site, season, k, res, flip)
        out[k] = res
        mark = klass(res["between"])
        if flip and mark != "?":
            mark = "骨格(符号反転)"
        print(f"        {nm:<10} raw {res['raw']:+.3f} → 年除去 {res['yr']:+.3f}"
              f"  **d {res['d']:+.3f}**  [{mark}]"
              f"   分解 between {res['between']:+.3f} + within {res['within']:+.3f}"
              f"（θ 年成分の分散比 {res['vshare']:.2f}）")
    return out


def show_delta(tal, arena, site, sp, au):
    """Δ = r_秋 − r_春 を raw と年除去の両方で出す（旗91/107 の主判定量）。

    **追補 A**：分解は Δ についても成り立つ（**各季節の恒等式の差**なので厳密）——
    `Δ_raw = Δ_between + Δ_within`。**判定は `Δ_between` に、同じ 3 段階の境で当てる。**
    """
    if sp is None or au is None:
        return
    for k, nm in (("le", "Δ θ→γLE"), ("h", "Δ θ→γH")):
        d_raw = au[k]["raw"] - sp[k]["raw"]
        d_yr = au[k]["yr"] - sp[k]["yr"]
        dd = d_yr - d_raw
        d_bet = au[k]["between"] - sp[k]["between"]
        d_win = au[k]["within"] - sp[k]["within"]
        flip = (np.isfinite(d_raw) and np.isfinite(d_win)
                and np.sign(d_raw) != np.sign(d_win))
        tal.add(arena + "/Δ", site, "Δ", k,
                {"raw": d_raw, "yr": d_yr, "n": 0, "d": dd,
                 "between": d_bet, "within": d_win, "vshare": np.nan}, flip)
        mark = klass(d_bet) + ("(符号反転)" if flip else "")
        print(f"        {nm:<10} raw {d_raw:+.3f} → 年除去 {d_yr:+.3f}"
              f"  **d {dd:+.3f}**  [{mark}]"
              f"   分解 between {d_bet:+.3f} + within {d_win:+.3f}")


# ------------------------------------------------------------- 土俵ごとの部分集合
def arena88(tal, site, d):
    print(f"\n    ── 旗88 の土俵（全日・春 4–5 月／秋 9–10 月）──")
    sp = d[[m in S88 for m in d.index.month]]
    au = d[[m in A88 for m in d.index.month]]
    a = measure_pair(tal, "旗88", site, "春", sp)
    b = measure_pair(tal, "旗88", site, "秋", au)
    show_delta(tal, "旗88", site, a, b)


def arena89(tal, site, d):
    print(f"\n    ── 旗89 の土俵（θ高×Rg高 の中を季節で割る）──")
    lab, _, _ = cell_of(d)
    hh = d[lab == "θ高×Rg高"]
    sp = hh[[SEASON[m] == "春" for m in hh.index.month]]
    au = hh[[SEASON[m] == "秋" for m in hh.index.month]]
    a = measure_pair(tal, "旗89", site, "春", sp)
    b = measure_pair(tal, "旗89", site, "秋", au)
    show_delta(tal, "旗89", site, a, b)


def arena91(tal, site, d):
    """旗90/91 の土俵：`θ高×Rg高` の中の **Ts 重なり帯**。"""
    print(f"\n    ── 旗91 の土俵（θ高×Rg高 かつ Ts 重なり帯）──")
    if "Ts" not in d.columns:
        print("      **Ts が無い**＝この土俵は無い"); return
    lab, _, _ = cell_of(d)
    hh = d[(lab == "θ高×Rg高") & d["Ts"].notna()]
    sp = hh[[m in SPRING for m in hh.index.month]]
    au = hh[[m in AUTUMN for m in hh.index.month]]
    if len(sp) == 0 or len(au) == 0:
        print("      春か秋が空＝この土俵は無い"); return
    lo, hi = band(sp["Ts"].to_numpy(), au["Ts"].to_numpy())
    if hi <= lo:
        print("      帯が作れない＝この土俵は無い"); return
    print(f"      帯 Ts ∈ [{lo:.1f}, {hi:.1f}]")
    a = measure_pair(tal, "旗91", site, "春", sp[(sp["Ts"] >= lo) & (sp["Ts"] <= hi)])
    b = measure_pair(tal, "旗91", site, "秋", au[(au["Ts"] >= lo) & (au["Ts"] <= hi)])
    show_delta(tal, "旗91", site, a, b)


def arena107(tal, site, d, P):
    """旗106/107 の土俵：`θ高×Rg高` かつ雨から `遠い`（≥7 日）。"""
    print(f"\n    ── 旗107 の土俵（θ高×Rg高 かつ雨から遠い ≥{REMOTE_MIN} 日）──")
    lab, _, _ = cell_of(d)
    hh = d[lab == "θ高×Rg高"]
    j = hh.join(rain_history(P, PRIMARY_THR), how="left")
    j["usable"] = j["usable"].fillna(False).astype(bool)
    j = j[j["usable"]]
    far = j[j["dry"] >= REMOTE_MIN]
    sp = far[[m in SPRING for m in far.index.month]]
    au = far[[m in AUTUMN for m in far.index.month]]
    a = measure_pair(tal, "旗107", site, "春", sp)
    b = measure_pair(tal, "旗107", site, "秋", au)
    show_delta(tal, "旗107", site, a, b)


# ------------------------------------------------------------- 門①（合成）
def synth(kind, years, seed, ndays=100):
    """**年内の連関**と**年をまたぐ対比**を作り分ける。

      ・`within_only`       —— 年内にだけ真の連関。年平均は θ と γ で**独立**（年分散は小）
      ・`between_only`      —— 年内の連関は 0。**年平均どうしにだけ**真の関係
      ・`no_between_signal` —— 年内に真の連関。**年平均は独立**（＝年の対比は**偶然だけ**・年分散は大）
    """
    rng = np.random.default_rng(seed)
    n = years * ndays
    g = np.repeat(np.arange(years), ndays)
    idx = pd.date_range("2000-01-01", periods=n, freq="D")
    Rg = rng.normal(0, 1, n)
    if kind == "within_only":
        s_y, beta_w, tie = 0.3, -0.8, False
    elif kind == "between_only":
        s_y, beta_w, tie = 1.0, 0.0, True
    else:                                    # no_between_signal
        s_y, beta_w, tie = 1.0, -0.8, False
    th_y = rng.normal(0, s_y, years)
    th_w = rng.normal(0, 1, n)
    th = th_y[g] + th_w
    y_y = -1.0 * th_y if tie else rng.normal(0, s_y, years)
    y = y_y[g] + beta_w * th_w + rng.normal(0, 1, n) + 0.5 * Rg
    return pd.DataFrame({"th": th, "Rg": Rg, "y": y}, index=idx), g


def gate_stats(kind, years, reps=REPS):
    """反復ごとの `r_raw`・`r_yr`・`d`・`r_between`・`r_within` を要約する。

    **追補 A 以降、門①が当てるのは `r_between` と `r_within`**（`d` は参考として残す）。
    """
    ds, raws, yrs_, bets, wins = [], [], [], [], []
    for i in range(reps):
        df, g = synth(kind, years, seed=1000 + i)
        res = r_set(df["y"].to_numpy(), df["th"].to_numpy(), df["Rg"].to_numpy(), g)
        ds.append(res["d"]); raws.append(res["raw"]); yrs_.append(res["yr"])
        bets.append(res["between"]); wins.append(res["within"])
    ds, raws, yrs_, bets, wins = map(
        lambda a: np.asarray(a, float), (ds, raws, yrs_, bets, wins))
    return dict(d_med=float(np.median(ds)), d_sd=float(np.std(ds, ddof=1)),
                d_abs=float(np.median(np.abs(ds))),
                raw_med=float(np.median(raws)), yr_med=float(np.median(yrs_)),
                yr_abs=float(np.median(np.abs(yrs_))),
                bet_med=float(np.median(bets)), bet_sd=float(np.std(bets, ddof=1)),
                bet_abs=float(np.median(np.abs(bets))),
                win_med=float(np.median(wins)), win_abs=float(np.median(np.abs(wins))))


def gates():
    """門①（追補 A の G1'〜G3'）。**合格条件は実データを見る前に固定済み。**

    旧 G1/G2/G3（`d` に当てる版）は**欠陥 #71 で無効**——`d` には分母の再正規化が混ざる。
    **参考として旧条件の数字も併記するが、合否には使わない。**
    """
    print("=== 旗144 門①（合成 3 本・**追補 A の G1'〜G3'**）——**実データより先に走らせる** ===")
    print(f"  反復 {REPS}・年あたり 100 日。**合格条件は `PREREGISTRATION_step144_amendment.md` で固定済み**。")
    print("  **当てる量は `r_between`／`r_within`**（`d` は欠陥 #71 のため合否に使わない）。")
    ok = {}

    g1 = gate_stats("within_only", 20)
    c1 = g1["bet_abs"] < 0.03 and g1["win_med"] < -0.15
    ok["G1'"] = c1
    print(f"\n  G1' `within_only`（20 年）：|r_between| 中央 {g1['bet_abs']:.3f}（< 0.03 が条件）"
          f"／r_within 中央 {g1['win_med']:+.3f}（< −0.15 が条件） → {'合格' if c1 else '**不合格**'}")
    print(f"       （参考・合否外）r_raw 中央 {g1['raw_med']:+.3f}／|d| 中央 {g1['d_abs']:.3f}")

    c2 = True
    print("\n  G2' `between_only`：")
    for yy in (3, 20):
        g2 = gate_stats("between_only", yy)
        c = g2["bet_med"] <= -0.15 and g2["win_abs"] < 0.05
        c2 = c2 and c
        print(f"    {yy:>2} 年：r_between 中央 {g2['bet_med']:+.3f}（≤ −0.15）"
              f"／|r_within| 中央 {g2['win_abs']:.3f}（< 0.05）"
              f"／（参考）r_raw 中央 {g2['raw_med']:+.3f} → {'合格' if c else '**不合格**'}")
    ok["G2'"] = c2

    print("\n  G3' `no_between_signal`（偶然の年対比だけ）：")
    sds, abss = {}, {}
    for yy in (3, 20):
        g3 = gate_stats("no_between_signal", yy)
        sds[yy] = g3["bet_sd"]; abss[yy] = g3["bet_abs"]
        print(f"    {yy:>2} 年：|r_between| 中央 {g3['bet_abs']:.3f}（< 0.05）"
              f"／r_between の SD {g3['bet_sd']:.3f}"
              f"／（参考）r_within 中央 {g3['win_med']:+.3f}・d 中央 {g3['d_med']:+.3f}")
    ratio = sds[3] / sds[20] if sds[20] > 0 else np.inf
    c3 = (abss[3] < 0.05 and abss[20] < 0.05
          and 1.8 <= ratio <= 3.5 and sds[20] < 0.05)
    ok["G3'"] = c3
    print(f"    条件：|r_between| 中央が両方 < 0.05／SD 比 {ratio:.2f}（1.8〜3.5・"
          f"平方根なら sqrt(20/3)=2.58）／SD(20年) {sds[20]:.3f} < 0.05"
          f" → {'合格' if c3 else '**不合格**'}")

    print(f"\n  === 門①のまとめ：{ok} ===")
    if all(ok.values()):
        print("  **3 本とも合格＝実データに進んでよい**（`--real`）。")
    else:
        print("  **落ちた門がある＝実データに進まない**（旗52 の作法）。")
        print("  **G3' が落ちた場合は「年数が多いから安全」という読み自体が成り立たない。**")
    return all(ok.values())


# ------------------------------------------------------------- G3' が落ちた理由の切り分け（合成のみ）
def synth_share(years, s_y, seed, ndays=100):
    """`no_between_signal` の**年分散の取り分だけ**を振る（他は同じ作り）。"""
    rng = np.random.default_rng(seed)
    n = years * ndays
    g = np.repeat(np.arange(years), ndays)
    Rg = rng.normal(0, 1, n)
    th_y = rng.normal(0, s_y, years)
    th_w = rng.normal(0, 1, n)
    th = th_y[g] + th_w
    y_y = rng.normal(0, s_y, years)                 # θ と**独立**＝年の対比は偶然だけ
    y = y_y[g] - 0.8 * th_w + rng.normal(0, 1, n) + 0.5 * Rg
    return th, y, Rg, g


def g3_sweep(reps=REPS):
    """**探索（事前登録外・合成のみ）**：G3' が落ちたのは年数のせいか、年分散の取り分のせいか。

    **`r_between` の大きさは `年分散の取り分 vshare` × `年平均どうしの偶然の相関` である。**
    後者は年数 `k` の平方根で薄まる（G3' の SD 比 2.15 はこれと整合）。**前者は年数と無関係。**
    **＝「年数が多いから安全」が成り立つ条件は、年数ではなく `vshare` が決める**——そのはずである。
    ここでは `vshare` を振って、**20 年で |r_between| 中央 < 0.05 に入る `vshare` の上限**を測る。
    """
    print("\n=== 探索（事前登録外・合成のみ）：G3' が落ちた理由 ===")
    print("  **合否には使わない。次の事前登録（追補 B）の閾値を導くための測定である。**")
    print(f"  反復 {reps}・年あたり 100 日・年の対比は偶然だけ（θ と γ の年平均は独立）。")
    print(f"  {'s_y':>5} {'年':>4} {'vshare':>7} {'|r_bet| 中央':>12} {'r_bet SD':>9}"
          f" {'恒等式の最大誤差':>16}")
    for s_y in (0.15, 0.25, 0.35, 0.5, 0.7, 1.0):
        for years in (3, 20):
            bs, vs, errs = [], [], []
            for i in range(reps):
                th, y, Rg, g = synth_share(years, s_y, seed=7000 + i)
                res = r_set(y, th, Rg, g)
                bs.append(res["between"]); vs.append(res["vshare"])
                errs.append(abs(res["raw"] - (res["between"] + res["within"])))
            bs = np.asarray(bs, float)
            print(f"  {s_y:>5.2f} {years:>4} {np.mean(vs):>7.3f}"
                  f" {np.median(np.abs(bs)):>12.3f} {np.std(bs, ddof=1):>9.3f}"
                  f" {max(errs):>16.2e}")
    print("  **恒等式の最大誤差が 1e-12 台なら `r_raw = r_between + r_within` は数値的に厳密**"
          "（回帰確認）。")


# ------------------------------------------------------------- 実データ
def real(sites, qc_max=None):
    tal = Tally()
    print("=== 旗144：骨格の r に年の対比はどれだけ混ざっているか（実データ）===")
    print("  **測るだけ。判定は差し替えない**（事前登録 step144・穴 #69 の規則）。")
    print(f"  段階（**追補 A：当てるのは r_between**）：|r_between| < {SMALL} 無視／"
          f"{SMALL} 以上 {BIG} 未満 効く／{BIG} 以上 または sign(r_within)≠sign(r_raw) 骨格")
    for s in sites:
        print(f"\n  ━━ {s} ━━")
        try:
            d, _ = daily_energy(s, list(range(1, 13)), qc_max, extra=("Ts",))
        except Exception as e:
            print(f"    読み込み失敗 {type(e).__name__}: {str(e)[:120]}"); continue
        if len(d) < MIN_DAYS or not {"th", "Rg", "gLE", "gH"} <= set(d.columns):
            print(f"    変数不足（{len(d)} 日・列 {list(d.columns)}）"); continue
        print(f"    全日数 {len(d):,}／年数 {d.index.year.nunique()}"
              f"／Ts {'あり' if 'Ts' in d.columns else 'なし'}")
        # **旗88/89 は Ts 無しで読んだ日集合を使う**（extra を足すと dropna が変わる＝旗36 の注意）
        try:
            d0, _ = daily_energy(s, list(range(1, 13)), qc_max)
        except Exception as e:
            print(f"    Ts 無しの読み直しに失敗 {type(e).__name__}"); d0 = d
        arena88(tal, s, d0)
        arena89(tal, s, d0)
        arena91(tal, s, d)
        try:
            P = daily_precip(s, qc_max)
        except Exception as e:
            P = None
            print(f"\n    ── 旗107 の土俵 ──\n      降水が読めない {type(e).__name__}")
        if P is not None and not P.dropna().empty:
            arena107(tal, s, d0, P)
        elif P is not None:
            print(f"\n    ── 旗107 の土俵 ──\n      **降水 P が空**")

    print("\n  === 集計（事前登録の規則に当てる）===")
    c = tal.counts()
    print(f"    判定できたセル {sum(v for k, v in c.items() if k != '?')}"
          f"／無視 {c['無視']}・効く {c['効く']}・**骨格 {c['骨格']}**（判定不能 {c['?']}）")
    bad = [r for r in tal.rows if klass(r["between"]) in ("効く", "骨格") or r["flip"]]
    if bad:
        print(f"    **|r_between| ≥ {SMALL} だったセル（名指しで残す）**：")
        for r in sorted(bad, key=lambda r: -(abs(r["between"])
                                             if np.isfinite(r["between"]) else np.inf)):
            print(f"      {r['arena']:<8}{r['site']:<8}{r['season']:<3}{r['half']:<3}"
                  f" raw {r['raw']:+.3f} = between {r['between']:+.3f}"
                  f" + within {r['within']:+.3f}（年除去 {r['yr']:+.3f}・d {r['d']:+.3f}）"
                  f"{'  **符号反転**' if r['flip'] else ''}")
    # 事前予測の採点
    na88 = [r for r in tal.rows if r["arena"] == "旗88" and r["site"] in NA3
            and r["season"] in ("春", "秋")]
    h1 = bool(na88) and all(abs(r["d"]) < SMALL for r in na88)
    a91 = [r for r in tal.rows if r["arena"] == "旗91" and r["season"] in ("春", "秋")]
    h2 = any(abs(r["d"]) >= SMALL for r in a91)
    cells = [r for r in tal.rows if r["season"] in ("春", "秋") and np.isfinite(r["d"])]
    opp = sum(1 for r in cells if np.sign(r["d"]) != np.sign(r["raw"]))
    h3 = opp >= 8
    print(f"\n    ★事前予測：H1（旗88 の北米 {len(na88)} セルすべて |d| < {SMALL}）"
          f"＝{'当たり' if h1 else '外れ'}")
    print(f"    ★事前予測：H2（旗91 の帯に |d| ≥ {SMALL} が 1 つ以上／対象 {len(a91)} セル）"
          f"＝{'当たり' if h2 else '外れ'}")
    print(f"    ★事前予測：H3（d の符号が r_raw と逆：{opp}/{len(cells)} セル・8 以上で当たり）"
          f"＝{'当たり' if h3 else '外れ'}")

    # 追補 A の H4〜H6（**当てる量は `r_between`**）
    b88 = [r for r in na88 if np.isfinite(r["between"])]
    h4 = bool(b88) and all(abs(r["between"]) < SMALL for r in b88)
    bcells = [r for r in tal.rows if r["season"] in ("春", "秋")
              and np.isfinite(r["between"])]
    top = max(bcells, key=lambda r: abs(r["between"])) if bcells else None
    h5 = top is not None and top["arena"] == "旗91"
    same = sum(1 for r in bcells if np.sign(r["between"]) == np.sign(r["raw"]))
    h6 = same >= 8
    print(f"    ★事前予測：H4（旗88 の北米 {len(b88)} セルすべて |r_between| < {SMALL}）"
          f"＝{'当たり' if h4 else '外れ'}")
    where = ("なし" if top is None else
             "{} {} {} {} {:+.3f}".format(top["arena"], top["site"], top["season"],
                                          top["half"], top["between"]))
    print(f"    ★事前予測：H5（|r_between| 最大が旗91 の帯：{where}）"
          f"＝{'当たり' if h5 else '外れ'}")
    print(f"    ★事前予測：H6（r_between の符号が r_raw と同じ：{same}/{len(bcells)} セル・"
          f"8 以上で当たり）＝{'当たり' if h6 else '外れ'}")

    print("\n  結論（事前登録の集計規則）：")
    if c["骨格"] == 0 and c["効く"] == 0:
        print("    **骨格の r は年の対比に依っていない。旗88/89/91/106/107 の記述は変えない。**")
    elif c["骨格"] == 0:
        print("    **効いているセルがある＝名指しで併記する。判定そのものは変えない。**")
    else:
        print("    **骨格に触るセルがある＝該当の主張に注記を入れ、次の事前登録の対象にする。**")
        print("    **本周では判定を差し替えない**（穴 #69 の規則）。")
    print("\n  留保（事前登録どおり）：")
    print("   ・**年除去の r が「正しい r」ではない。** 生の r は日と年をまとめた連関、")
    print("     年除去の r は**年内の連関**である。**誤りは前者を後者と読むことだけ。**")
    print("   ・**本測定は季節差の実在を再判定しない。**")
    return tal


def main():
    ap = argparse.ArgumentParser(description="旗144：骨格の r に混ざる年の対比を測る")
    ap.add_argument("--real", action="store_true", help="実データ（/mnt/hdd）")
    ap.add_argument("--g3sweep", action="store_true",
                    help="G3' が落ちた理由の切り分け（探索・合成のみ）")
    ap.add_argument("--sites", nargs="+", default=list(NA3) + list(MN3))
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()
    tee_stdout("step144")
    if a.g3sweep:
        g3_sweep()
        return 0
    if not a.real:
        gates()
        return 0
    real(a.sites, a.qc_max)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
