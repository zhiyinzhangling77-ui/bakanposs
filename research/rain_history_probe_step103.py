"""旗103：**手C（降雨イベントからの経過）の前提の下調べ**——**検定はしない**。

`OPEN_QUESTIONS_OPTIONS.md` の**手C**：**同じ θ でも「大きな雨の 3 日後」と
「小さな雨が続いた末の同じ θ」は違う**のではないか。**`P` は測定量**であり、
**旗45/77 は呼吸の Birch 効果に使ったが、Bowen 反転には一度も使っていない。**

**この道具は分布と日数しか出さない**（旗94/98/101 と同じ作法）——
**θ→γLE・θ→γH は一行も計算しない**。**検定の答えを事前登録の前に見ないため**である。

## なぜ「事前登録」ではなく「下調べ」から始めるか（**旗96/97 の授業料**）

**旗96（春は高 VPD のはず → 実際は低かった）**と**旗97（春は深層が乾いているはず →
実際は逆または差が無かった）**は、**どちらも一行で確かめられる分布を見ずに設計に入っていた**。
**旗97 で作法に加えた**——**事前登録の前に、前提にする事実を分布で確かめる**。
**手C の前提**は「**春と秋では雨の入り方が違う**」であり、**それは P の分布で確かめられる**。

## もう一つ、**先に確かめないと事前登録が無駄になる**もの

**旗97 の追補で分かった構造**：**重なり帯は「群間で違う量」を揃えるのに使うと自己矛盾する。**
手C も同じ形を踏む危険がある——**「大きな雨の直後」の日は当然 θ が高い**ので、
**θ を揃えようとすると帯が作れない／帯の中が下限未満になる**。
**そうなるかどうかは、検定を組む前に日数を数えれば分かる。** 本道具はそれを数える。

## 出すもの（サイト × 春/秋）

  1. **P（降水）が読めるか**・有効日数・年数・年降水量・雨日率
  2. **イベントの閾値候補**（5 mm・10 mm・**雨日の 90 パーセンタイル**）ごとのイベント数
  3. **最後のイベントからの日数**（`dryspell`）の分布と、**直前 7/30 日積算 P** の分布
  4. **層別候補の日数**：`直後`（dryspell ≤ 3）／`遠い`（dryspell ≥ 7）が
     **下限（60 日・3 年）に届くか**——**春と秋のそれぞれで**
  5. **θ を揃えられるか**：`直後`と`遠い`の **θ_1 の重なり帯**の幅と、**帯の中の日数**
  6. **dryspell と積算 P が同じ量になっていないか**（Spearman）——
     **同じなら二つ登録する意味が無い**

    python research/rain_history_probe_step103.py            # 合成で道具を検証（既定）
    python research/rain_history_probe_step103.py --check    # dryspell の自己試験だけ
    python research/rain_history_probe_step103.py --real \\
        --sites US-Wkg US-Whs US-SRM
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaporation_regime_step36 import daily_energy
# **既にある道具を呼ぶ**（旗90/96/97 と同じ作法＝しきい値も定義も作り直さない）
from stratified_bowen_step89 import cell_of, MIN_DAYS, MIN_YEARS
from soiltemp_match_step90 import band, SPRING, AUTUMN
from vpd_match_step96 import spearman
from precip_pressure_test_step77 import dryspell, _cum

EVENT_THR = (5.0, 10.0)          # イベントの閾値候補 [mm/日]（**下調べなので併記する**）
PRIMARY_THR = 5.0                # 日数を数えるときの主閾値（**事前登録ではまだ固定しない**）
RECENT_MAX, REMOTE_MIN = 3, 7    # `直後` ≤3 日／`遠い` ≥7 日（**候補**）
WARMUP = 30                      # 記録先頭 30 日は使わない（**下の欠陥 #35**）
CUM_WINS = (7, 30)


def daily_precip(site, qc_max=None):
    """**日降水量 [mm/日]**。**`daily_energy` の `extra` では取れない**。

    理由は二つあり、**どちらも実データで結果を壊す**：
      ・`daily_energy` は全列を **平均**で日集約する。**降水は合計でなければ意味が無い**
        （30 分値の平均は日雨量の 1/48）。
      ・`daily_energy` は最後に `dropna()` する。**列を増やすと落ちる日が変わり、
        旗89 以降のセル定義が動く**（`daily_energy` の docstring にある警告そのもの）。
    ＝**P は別に読み、合計で日集約し、あとからセルの日に貼る。**
    """
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    cfg = AnalysisConfig(qc_max=qc_max) if qc_max is not None else AnalysisConfig()
    raw = load_raw_all(get_site(site), cfg)
    if "P" not in raw.columns:
        return None
    p = raw["P"].groupby(raw.index.normalize()).sum(min_count=1)
    return p.dropna()


def rain_history(P, thr=PRIMARY_THR):
    """**連続した日付の上で**降雨履歴を作る。**戻り値は日付を索引とする DataFrame**。

    **落とし穴（この道具の要点）**：θ高×Rg高 のセルの日は**飛び飛び**である。
    **その部分列の上で「最後の雨からの日数」を数えると、全く別の量になる**——
    間に挟まった日が見えないので、**乾燥期間を過小に数える**。
    ＝**P の全期間を日単位で埋め直してから数え、あとでセルの日に貼る。**

    **道具の欠陥 #35（旗77 の `dryspell` から引き継ぐ性質）**：
    `dryspell` は**記録先頭を 0 から始める**——**「乾燥期間が未知」と「昨日雨が降った」を
    区別しない**。**欠測を跨いだときも値を持ち越すだけで増やさない**（乾燥期間の過小評価）。
    旗77 は残差の説明変数として使ったので影響は小さかったが、
    **手C は `dryspell` そのもので層別する**ので効く。
    ＝ここで **`usable`** を作って落とす：**記録先頭 30 日**（`WARMUP`）・
    **直前 7 日に欠測がある日**・**まだ一度もイベントが起きていない日**は使わない。
    """
    full = pd.date_range(P.index.min(), P.index.max(), freq="D")
    pf = P.reindex(full)
    v = pf.to_numpy()
    ds = pd.Series(dryspell(v, thr), index=full)
    cum = _cum(v, CUM_WINS)                       # **直前 N 日（当日を除く）**＝旗77 と同一
    ev = (pf > thr).fillna(False).to_numpy()
    seen = np.cumsum(ev) > 0                      # **一度でもイベントが起きたか**
    ok7 = (pf.notna().astype(float)          # **bool のまま rolling しない**
           .rolling(7, min_periods=7).sum().eq(7).to_numpy())
    warm = np.arange(len(full)) >= WARMUP
    out = pd.DataFrame({"P": pf.to_numpy(), "dry": ds.to_numpy(),
                        "usable": seen & ok7 & warm}, index=full)
    for w, c in zip(CUM_WINS, cum):
        out[f"cum{w}"] = c
    return out


def _q(x, tag, unit=""):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return f"{tag}：なし"
    return (f"{tag}：n={x.size} 中央 {np.median(x):.1f}{unit} "
            f"[{np.percentile(x, 10):.1f}–{np.percentile(x, 90):.1f}]")


def season_probe(name, sub):
    """**1 季節ぶんの分布**を出す（**検定はしない**）。**層別候補の日数**を返す。"""
    if sub.empty:
        print(f"      {name}：セル内に日が無い")
        return None
    yrs = sub.index.year.nunique()
    print(f"      {name}：セル内 {len(sub)} 日／{yrs} 年"
          f"／うち使える日 {int(sub['usable'].sum())}")
    u = sub[sub["usable"]]
    if u.empty:
        print(f"        **使える日が無い**（欠測か記録先頭）")
        return {"n": 0, "yrs": 0}
    print(f"        {_q(u['dry'], '最後のイベントからの日数', ' 日')}")
    for w in CUM_WINS:
        print(f"        {_q(u[f'cum{w}'], f'直前 {w} 日の積算 P', ' mm')}")
    rec = u[u["dry"] <= RECENT_MAX]
    rem = u[u["dry"] >= REMOTE_MIN]
    for tag, g in (("直後(≤%d 日)" % RECENT_MAX, rec), ("遠い(≥%d 日)" % REMOTE_MIN, rem)):
        ny = g.index.year.nunique() if len(g) else 0
        ok = len(g) >= MIN_DAYS and ny >= MIN_YEARS
        med = f"／θ 中央 {g['th'].median():.3f}" if len(g) else ""
        print(f"        {tag}：{len(g)} 日／{ny} 年"
              f"  {'**下限を満たす**' if ok else '**下限未満**（60 日・3 年）'}{med}")
    return {"n": len(u), "yrs": u.index.year.nunique(),
            "recent": len(rec), "remote": len(rem),
            "recent_ok": len(rec) >= MIN_DAYS and rec.index.year.nunique() >= MIN_YEARS
            if len(rec) else False,
            "remote_ok": len(rem) >= MIN_DAYS and rem.index.year.nunique() >= MIN_YEARS
            if len(rem) else False,
            "rec": rec, "rem": rem}


def matchability(name, st):
    """**θ を揃えられるか**——`直後`と`遠い`の **θ_1 の重なり帯**を測る。

    **旗97 の追補で判った自己矛盾**をここで先に検出する：
    **「大きな雨の直後は θ が高い」のは当たり前**なので、
    **帯が作れない／帯の中が下限未満なら、手C は θ を揃えた形では設計できない。**
    """
    if not st or st.get("recent") in (None, 0) or st.get("remote") in (None, 0):
        print(f"      {name}：**片方が空**＝帯を測れない＝**この形では設計できない**")
        return {"band": None, "n_rec": 0, "n_rem": 0, "ok": False}
    rec, rem = st["rec"], st["rem"]
    lo, hi = band(rec["th"].to_numpy(), rem["th"].to_numpy())
    print(f"      {name}：θ 中央 直後 {rec['th'].median():.3f}／遠い {rem['th'].median():.3f}")
    if hi <= lo:
        print(f"        → **θ の帯が作れない**（[{lo:.3f}, {hi:.3f}]）"
              f"＝**この形では交絡を潰せない**")
        return {"band": (lo, hi), "n_rec": 0, "n_rem": 0, "ok": False}
    a = rec[(rec["th"] >= lo) & (rec["th"] <= hi)]
    b = rem[(rem["th"] >= lo) & (rem["th"] <= hi)]
    ok = (len(a) >= MIN_DAYS and len(b) >= MIN_DAYS
          and a.index.year.nunique() >= MIN_YEARS and b.index.year.nunique() >= MIN_YEARS)
    print(f"        **帯 [{lo:.3f}, {hi:.3f}]**（幅 {hi-lo:.3f}）／帯の中："
          f"直後 {len(a)} 日・{a.index.year.nunique()} 年／"
          f"遠い {len(b)} 日・{b.index.year.nunique()} 年"
          f"  {'**下限を満たす**' if ok else '**下限未満**'}")
    return {"band": (lo, hi), "n_rec": len(a), "n_rem": len(b), "ok": ok}


def run_site(tag, d, P):
    """1 サイト。**セル → 降雨履歴 → 春/秋の分布 → 層別候補の日数 → θ を揃えられるか**。"""
    print(f"\n  ━━ {tag} ━━")
    if P is None or P.dropna().empty:
        print("    **降水 P が無い**＝手C はこのサイトでは設計できない")
        return {"P": False}
    lab, tmed, rmed = cell_of(d)
    hh = d[lab == "θ高×Rg高"].copy()
    print(f"    しきい値（旗89 と同一・作り直さない）：θ={tmed:.3f}／Rg={rmed:.1f}"
          f"／θ高×Rg高 {len(hh)} 日")
    yrs = P.index.year.nunique()
    wet = P[P >= 1.0]
    p90 = float(np.percentile(wet.to_numpy(), 90)) if len(wet) else float("nan")
    print(f"    P：{len(P)} 日／{yrs} 年／年降水 {P.sum()/max(yrs,1):.0f} mm"
          f"／雨日率（≥1 mm）{len(wet)/max(len(P),1):.0%}"
          f"／雨日の 90 パーセンタイル {p90:.1f} mm")
    for t in EVENT_THR:
        n = int((P >= t).sum())
        print(f"      イベント（≥{t:.0f} mm）：{n} 回＝年 {n/max(yrs,1):.1f} 回")
    hist = rain_history(P, PRIMARY_THR)
    j = hh.join(hist, how="left")
    # **P の期間外のセル日は使えない**。join で入った欠損を bool に畳む
    j["usable"] = j["usable"].fillna(False).astype(bool)
    r = spearman(j.loc[j["usable"], "dry"], j.loc[j["usable"], "cum7"])
    print(f"    **dryspell と直前 7 日積算 P の Spearman r = {r:+.2f}**"
          f"{'  ← **ほぼ同じ量**（別々に登録する意味が薄い）' if np.isfinite(r) and abs(r) > 0.90 else ''}")
    print(f"    ── 季節ごと（イベント閾値 {PRIMARY_THR:.0f} mm）──")
    out = {"P": True, "r_dry_cum7": r}
    for nm, mon in (("春", SPRING), ("秋", AUTUMN)):
        st = season_probe(nm, j[[m in mon for m in j.index.month]])
        out[nm] = st
    print(f"    ── θ を揃えられるか（`直後` 対 `遠い`・季節の中だけで）──")
    for nm in ("春", "秋"):
        out[nm + "_match"] = matchability(nm, out.get(nm))
    return out


def synth(kind, years=16, seed=0):
    """**下調べの機械が正しく数えるか**を試す四つ。**検定の合成ではない**。

      ・`separable`  —— 春も秋も雨が十分あり、**θ は雨の直後性だけでは決まらない**
                        → **両季節で下限を満たし、θ の帯も作れる**べき
      ・`confounded` —— **θ がほぼ雨の直後性で決まる**
                        → **θ の帯が作れない／片方が空／帯の中が下限未満**のどれかが出るべき
                        （＝**この形では設計できない**と、事前登録の前に分かる）
                        （＝**旗97 の自己矛盾を、事前登録の前に検出できるか**）
      ・`sparse`     —— **春はほとんど降らない**
                        → **春の `直後` が下限未満**と出るべき
      ・`nodata`     —— **P が無い** → **「P が無い」で落ちる**べき
    """
    rng = np.random.default_rng(seed)
    if kind == "separable":
        # **旗104：この対照だけ 30 年**。16 年では**秋が帯の中で下限に届かない**——
        # 秋の `Rg高` は 9 月にほぼ限られ、セル内の日数がそもそも足りない（下の較正の記録）。
        # **対照の役目は「下限を満たす」枝が到達可能で正しいと示すこと**であって、
        # **実サイトの年数を模すことではない**。**実データ側の年数はいじらない。**
        years = 30
    idx = pd.date_range("2008-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy()
    Rg = np.clip(150 + 120 * np.sin(2 * np.pi * (doy - 80) / 365)
                 + rng.normal(0, 25, len(idx)), 5, None)
    # 降雨：春（前線性）と夏〜秋（モンスーン）に山を持つ確率で発生させる
    spring = np.exp(-0.5 * ((doy - 100) / 30.) ** 2)
    monsoon = np.exp(-0.5 * ((doy - 230) / 40.) ** 2)
    if kind == "sparse":
        # **春をほとんど降らせない**（春の山を消し、春の底も下げる）
        lam = 0.10 * (1 - 0.95 * spring) + 0.30 * monsoon
    elif kind == "separable":
        # **旗104 で較正した**（旗103 は一度も走らせずに書いたので値が合っていなかった）。
        # 旧値は `confounded` と共用の `0.10 + 0.22*spring + 0.30*monsoon` で、
        # **イベントが年 54 回＝平均間隔 6.8 日**。**`遠い`(≥7 日) は構造的に作れなかった**——
        # ＝**陽性の対照が陽性を出せていなかった**。**雨を疎にして間隔を空ける。**
        lam = 0.070 + 0.150 * spring + 0.170 * monsoon
    else:
        lam = 0.10 + 0.22 * spring + 0.30 * monsoon
    wet = rng.random(len(idx)) < np.clip(lam, 0.005, 1)
    P = np.where(wet, rng.gamma(1.3, 7.0, len(idx)), 0.0)
    ds = dryspell(P, PRIMARY_THR)
    recent = np.exp(-np.nan_to_num(ds, nan=30.) / 4.0)               # 直後ほど 1 に近い
    if kind == "confounded":
        th = np.clip(0.10 + 0.30 * recent + rng.normal(0, 0.005, len(idx)), .02, .6)
    elif kind == "separable":
        # **旗104 で較正**：θ を**雨の直後性から切り離す**。
        # 旧式は `0.16 + 0.06*recent + 季節 + 雑音` で、**θ高 のセルが `直後` の日に偏っていた**。
        # **ゆっくり動く独立成分**（深層貯留のような、雨の直後性では決まらない量）を主にする。
        # **これが `separable` の定義そのもの**——「θ は雨の直後性だけでは決まらない」。
        slow = (pd.Series(rng.normal(0, 1, len(idx)))
                .rolling(45, min_periods=1).mean().to_numpy())
        slow = slow / (np.std(slow) + 1e-12)
        th = np.clip(0.16 + 0.015 * recent + 0.065 * slow
                     + 0.040 * np.sin(2 * np.pi * (doy - 200) / 365)
                     + rng.normal(0, 0.020, len(idx)), .02, .6)
    else:
        th = np.clip(0.16 + 0.06 * recent
                     + 0.06 * np.sin(2 * np.pi * (doy - 200) / 365)
                     + rng.normal(0, 0.045, len(idx)), .02, .6)
    d = pd.DataFrame({"th": th, "Rg": Rg}, index=idx)
    if kind == "nodata":
        return d, None
    return d, pd.Series(P, index=idx)


def check():
    """**`dryspell` と `usable` の自己試験**。**答えが手で分かる系列**で確かめる。"""
    print("  【自己試験】**記録先頭・欠測・イベント無しを正しく落とすか**")
    idx = pd.date_range("2020-01-01", periods=60, freq="D")
    p = np.zeros(60)
    p[40] = 20.0                       # **40 日目に初めてのイベント**
    P = pd.Series(p, index=idx)
    h = rain_history(P, PRIMARY_THR)
    bad = h.loc[h.index[:40], "usable"].any()
    print(f"    ① **初イベント前の 40 日**が使えると誤判定：{bad}（False であるべき）")
    print(f"    ② 41 日目の dryspell = {h['dry'].iloc[41]:.0f}（1 であるべき）／"
          f"使える = {bool(h['usable'].iloc[41])}（True であるべき）")
    P2 = P.copy(); P2.iloc[45:48] = np.nan
    h2 = rain_history(P2, PRIMARY_THR)
    print(f"    ③ 欠測直後（46〜54 日目）に使える日があるか："
          f"{bool(h2['usable'].iloc[46:54].any())}（**直前 7 日に欠測がある間は False**）")
    # **セルの日だけで数えると別物になる**ことを見せる
    cell = idx[::3]
    wrong = pd.Series(dryspell(P.reindex(cell).to_numpy(), PRIMARY_THR), index=cell)
    right = h["dry"].reindex(cell)
    dif = int((wrong != right).sum())
    print(f"    ④ **飛び飛びの日の上で数えた dryspell** は "
          f"{dif}/{len(cell)} 日で全期間から作った値と食い違う"
          f"（**0 でないことがこの道具の存在理由**）")
    ok = (not bad) and h["dry"].iloc[41] == 1 and bool(h["usable"].iloc[41]) and dif > 0
    print(f"    → 自己試験 {'**通過**' if ok else '**失敗**'}")
    return ok


def main():
    ap = argparse.ArgumentParser(description="旗103：手C の前提の下調べ（検定はしない）")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--check", action="store_true", help="dryspell の自己試験だけ")
    ap.add_argument("--sites", nargs="+", default=["US-Wkg", "US-Whs", "US-SRM"])
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()

    print("=== 旗103：手C（降雨イベントからの経過）の**前提の下調べ** ===")
    print("  **検定はしない**——θ→γLE・θ→γH は一行も計算しない（旗94/98/101 と同じ作法）。")
    print("  **出すのは分布と日数だけ**：**春と秋で雨の入り方が違うか**と、")
    print("  **`直後` 対 `遠い` の層別が下限に届くか**と、**θ を揃えられるか**。")
    print(f"  下限は旗58 以来と同じ（{MIN_DAYS} 日・{MIN_YEARS} 年）。"
          f"イベントの主閾値 {PRIMARY_THR:.0f} mm/日 は**下調べの値であって、"
          f"事前登録ではまだ固定しない**。")

    if a.check:
        check(); return

    if not a.real:
        print("\n  【合成データで道具を検証する】**実データに触れる前に**、")
        print("  **数え方が正しいか**（自己試験）と、**四つの場合を区別できるか**を見る。")
        check()
        for kind, want in (
                ("separable", "**両季節で下限を満たし、θ の帯も作れる**べき"),
                ("confounded", "**θ の帯が作れない／片方が空／帯の中が下限未満**が出るべき"),
                ("sparse", "**春の `直後` が下限未満**と出るべき"),
                ("nodata", "**「P が無い」で落ちる**べき")):
            print(f"\n  ===== 合成 `{kind}` —— 期待：{want} =====")
            d, P = synth(kind)
            run_site(f"合成-{kind}", d, P)
        print("\n  → **四つを区別できていれば道具は使える**（旗52 以来の作法）。")
        print("     **`confounded` で帯が潰れることを確かめる意味**は、")
        print("     **旗97 の追補（重なり帯は群間で違う量に使うと自己矛盾する）を、")
        print("     今回は事前登録の前に検出するためである。**")
        return

    res = {}
    for s in a.sites:
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max)
            P = daily_precip(s, a.qc_max)
        except Exception as e:
            print(f"\n  ━━ {s} ━━\n    読み込み失敗 {type(e).__name__}: {str(e)[:110]}")
            continue
        res[s] = run_site(s, d, P)

    print("\n  === この下調べが次に決めさせること（**結論ではない**）===")
    print("  ① **前提**：春と秋で `最後のイベントからの日数`・`直前 7/30 日積算 P` の")
    print("     分布が**実際に違うか**。**違わなければ手C の仮説は前提から成り立たない**")
    print("     （**旗96/97 と同じ失敗を、今度は検定を組む前に見つける**）。")
    print("  ② **設計可能性**：`直後`／`遠い` が**下限に届くか**。届かない季節・サイトは")
    print("     **事前登録の対象から外す**（**走らせて「判定しない」を並べない**＝旗97 の作法）。")
    print("  ③ **交絡を潰せるか**：θ の帯が作れるか。**作れないなら、θ を揃える形は捨てて")
    print("     別の形（例：θ を共変量に入れる）で登録するしかない**——**それを先に知る**。")
    print("  ④ **二つの層別変数が同じ量でないか**（dryspell と積算 P の r）。")
    print("\n  留保：")
    print("   ・**`P_F` はギャップフィル済みの可能性がある**（旗46 が `SWC_F_MDS` で見たのと同じ問題）。")
    print("     **手元では埋めた日と観測日を区別できない**——**そう書いて読む**。")
    print("   ・**独立クラスタは 2 つ**（Walnut Gulch・Santa Rita）＝3 サイト≠3 反復。")
    print("   ・**この道具は何も検定していない**。**★も▲も出さない。**")


if __name__ == "__main__":
    main()
