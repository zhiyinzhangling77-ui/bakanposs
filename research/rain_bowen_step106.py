"""旗106：**Bowen 反転は「雨からの日数」で説明できるか**（手C・事前登録 step106）。

## 論理（**秋の中だけで決着する**）

**春が反転しないのは、春の日が「雨から遠い」日ばかりだからか。**
**旗105 で前提は確かめてある**——**秋は春の 2.4〜3.5 倍「雨の直後」の日を持つ。**

- **反転に「最近濡れたこと」が要るなら、秋の `遠い` 層でも反転は消えるはず。**
- **秋の `遠い` 層で反転が残るなら、`遠い` であること自体は反転を妨げない**
  ＝**春が反転しない理由は「雨からの日数」ではない。**

**春は主検定に入れない**——**春の `直後` は US-Wkg 47 日・US-Whs 44 日で下限未満**（旗105）。
**春で層別できるのは US-SRM だけ＝独立クラスタ 1。** **参考として出す。**

## **門（旗104 の失敗を踏まえて先に決めた）**

  ・**門①-a（実データ）**：**そのサイトの秋全体で反転すること。** しなければ**判定しない**。
    **旗96 の秋の値は「重なり帯の中」で測ったもの**であり、**本検定は帯を使わない**
    ——**指定が変わった以上、対照は測り直す。**
  ・**門①-b（合成）**：**期待する枝に実際に到達することを確かめてから実データに触れる。**
    **「期待どおりに出るはず」と書いただけの対照は対照ではない**（旗104 の規則）。

    python research/rain_bowen_step106.py                    # 合成で検証（既定）
    python research/rain_bowen_step106.py --real             # 実データ（/mnt/hdd）
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stratified_bowen_step89 import cell_of, test_cell, reversed_, MIN_DAYS, MIN_YEARS
from soiltemp_match_step90 import SPRING, AUTUMN
from rain_history_probe_step103 import (rain_history, daily_precip,
                                        PRIMARY_THR, RECENT_MAX, REMOTE_MIN)
from evaporation_regime_step36 import daily_energy

SITES = ("US-Wkg", "US-Whs", "US-SRM")     # 旗82/89/96/105 と同一（**独立クラスタは 2**）


def strata(j):
    """`直後`（≤3 日）と `遠い`（≥7 日）に割る。**`usable` でない日は使わない**（欠陥 #35 の手当て）。"""
    u = j[j["usable"]]
    return {"直後": u[u["dry"] <= RECENT_MAX], "遠い": u[u["dry"] >= REMOTE_MIN]}


def show(tag, sub):
    """1 層を測って印字し、**反転したか**を返す。**θ の分布を必ず併記する**（事前登録どおり）。"""
    n, y = len(sub), (sub.index.year.nunique() if len(sub) else 0)
    if n < MIN_DAYS or y < MIN_YEARS:
        print(f"      {tag:<6}{n:>4} 日／{y:>2} 年  **下限未満**")
        return None
    th = sub["th"].to_numpy()
    res = test_cell(sub)
    if res is None:
        print(f"      {tag:<6}{n:>4} 日／{y:>2} 年  測れない")
        return None
    (rl, cl, _), (rh, ch, _) = res["le"], res["h"]
    rv = reversed_(res)
    print(f"      {tag:<6}{n:>4} 日／{y:>2} 年  "
          f"θ 中央 {np.nanmedian(th):.2f} [{np.nanpercentile(th,25):.2f}–"
          f"{np.nanpercentile(th,75):.2f}]  "
          f"θ→γLE {rl:+.2f} [{cl[0]:+.2f},{cl[1]:+.2f}]  "
          f"θ→γH {rh:+.2f} [{ch[0]:+.2f},{ch[1]:+.2f}]  "
          f"→ {'**反転**' if rv else '反転せず'}")
    return bool(rv)


def run_site(tag, d, P, season=AUTUMN, season_name="秋"):
    """1 サイト・1 季節。**門①-a → 層別**。**判定は返り値で返し、印字と分けない**。"""
    print(f"\n  ━━ {tag}（{season_name}）━━")
    if P is None or P.dropna().empty:
        print("    **降水 P が無い**"); return None
    lab, tmed, rmed = cell_of(d)
    hh = d[lab == "θ高×Rg高"]
    hh = hh[[m in season for m in hh.index.month]]
    j = hh.join(rain_history(P, PRIMARY_THR), how="left")
    j["usable"] = j["usable"].fillna(False).astype(bool)

    print(f"    ── 門①-a：{season_name}全体（層別せず）──")
    whole = show("全体", j)
    if not whole:
        print(f"    → **判定しない**（**門①-a を通らない**"
              f"{'：測れない' if whole is None else '：反転しない'}）")
        return None

    print(f"    ── 層別（`直後` ≤{RECENT_MAX} 日／`遠い` ≥{REMOTE_MIN} 日）──")
    st = strata(j)
    a, b = show("直後", st["直後"]), show("遠い", st["遠い"])
    if a is None or b is None:
        print("    → **判定しない**（どちらかの層が下限未満）"); return None
    if a and not b:
        v = "★雨からの日数で説明できる"
    elif a and b:
        v = "▲両層とも反転＝説明しない"
    elif b and not a:
        v = "○予想と逆（`遠い` だけ反転）"
    else:
        v = "○どちらも反転せず"
    print(f"    → **{v}**")
    return v


def synth(kind, years=22, seed=0):
    # **`theta_partial` だけ 30 年**——**この対照の役目は「両層とも下限を超えた状態で
    # 交絡を持たせる」こと**であり、22 年では `遠い` が 61 日と縁すぎる（旗104 と同じ理由・
    # **実データ側の年数はいじらない**）。
    """**四つとも「期待する枝に到達すること」を確かめる**（旗104 の規則）。

      ・`rain_driven` —— **反転が `直後` の日にだけ起きる** → **★を返すべき**
      ・`rain_free`   —— **反転が dryspell に依らず起きる** → **▲を返すべき**
      ・`theta_driven`—— **反転が θ の水準だけで決まり、θ は dryspell と独立**
                         → **▲を返すべき**
      ・`theta_partial` —— **反転が θ の水準だけで決まり、θ が雨の直後性に「部分的に」従う**
                         （＝**実データの形**・旗105）→ **★を返してしまわないかの試験・最重要**
      ・`theta_confounded` —— **同じだが交絡が極端** → **`遠い` が潰れて判定しないはず**
                         （**極端な交絡では道具が自ら判定を拒む**ことの確認）
      ・`no_reversal` —— **そもそも反転しない** → **門①-a で落ちるべき**
    """
    if kind == "theta_partial":
        years = 30
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2003-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy()
    Rg = np.clip(150 + 120 * np.sin(2 * np.pi * (doy - 80) / 365)
                 + rng.normal(0, 25, len(idx)), 5, None)
    # **雨は疎にする**（旗104 の較正：密だと `遠い` が構造的に作れない）
    spring = np.exp(-0.5 * ((doy - 100) / 30.) ** 2)
    monsoon = np.exp(-0.5 * ((doy - 250) / 45.) ** 2)
    lam = 0.055 + 0.110 * spring + 0.150 * monsoon
    wet = rng.random(len(idx)) < np.clip(lam, 0.005, 1)
    P = np.where(wet, rng.gamma(1.3, 7.0, len(idx)), 0.0)
    from precip_pressure_test_step77 import dryspell
    ds = dryspell(P, PRIMARY_THR)
    recent = np.exp(-np.nan_to_num(ds, nan=30.) / 4.0)
    # **θ の作り方は仕込みで変える**。
    # 既定は**雨の直後性から切り離した θ**（旗104 の `separable` と同じ作り）。
    # **`theta_confounded` だけは θ を雨の直後性に強く結びつける**——
    # **実データがその形だから**である（旗105：`直後` の θ 中央は `遠い` より高い。
    # 春 12.32 対 9.78・7.39 対 4.52・6.29 対 5.03）。
    # **第1版の `theta_driven` は θ を切り離したまま作っており、
    # 両層の θ 中央が 0.21 で揃っていた＝交絡の試験になっていなかった**（旗104 と同じ形）。
    if kind == "theta_confounded":
        th = np.clip(0.10 + 0.30 * recent + rng.normal(0, 0.010, len(idx)), .02, .6)
    elif kind == "theta_partial":
        # **交絡を「実データくらい」に弱める**。`recent` の係数を掃引して選んだ
        # （0.015→θ 比 1.02／0.05→1.07／0.12→1.18 だが 0.12 では `遠い` が 51 日で潰れる）。
        # **0.05 は `遠い` が下限をぎりぎり超える最大の交絡**である。
        slow = (pd.Series(rng.normal(0, 1, len(idx))).rolling(45, min_periods=1)
                .mean().to_numpy())
        slow = slow / (np.std(slow) + 1e-12)
        th = np.clip(0.16 + 0.050 * recent + 0.065 * slow
                     + 0.040 * np.sin(2 * np.pi * (doy - 200) / 365)
                     + rng.normal(0, 0.020, len(idx)), .02, .6)
    else:
        slow = (pd.Series(rng.normal(0, 1, len(idx))).rolling(45, min_periods=1)
                .mean().to_numpy())
        slow = slow / (np.std(slow) + 1e-12)
        th = np.clip(0.16 + 0.015 * recent + 0.065 * slow
                     + 0.040 * np.sin(2 * np.pi * (doy - 200) / 365)
                     + rng.normal(0, 0.020, len(idx)), .02, .6)
    thz = (th - th.mean()) / th.std()

    # ── 反転の強さ g（θ が γLE を上げ γH を下げる度合い）を仕込み分けで変える ──
    if kind == "rain_driven":
        g = 1.0 * (ds <= RECENT_MAX)                 # **`直後` の日だけ反転**
    elif kind == "rain_free":
        g = np.full(len(idx), 1.0)                   # **いつでも反転**
    elif kind in ("theta_driven", "theta_partial", "theta_confounded"):
        g = 1.0 * (thz > 0)                          # **θ が高い日だけ反転**（雨の日数は無関係）
    else:                                            # no_reversal
        g = np.zeros(len(idx))
    avail = 0.75 * Rg
    frac = np.clip(0.45 + 0.22 * g * thz + rng.normal(0, 0.05, len(idx)), .05, .95)
    gLE = avail * frac
    gH = avail * (1 - frac)
    d = pd.DataFrame({"th": th, "Rg": Rg, "gLE": gLE, "gH": gH}, index=idx)
    return d, pd.Series(P, index=idx)


def main():
    ap = argparse.ArgumentParser(description="旗106：Bowen 反転は雨からの日数で説明できるか")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--qc-max", type=int, default=None,
                    help="旗103/105 と同じ既定（None＝設定ファイルのまま）")
    a = ap.parse_args()

    print("=== 旗106：Bowen 反転は「雨からの日数」で説明できるか（手C）===")
    print("  **秋の中だけで決着する**——**春は `直後` が下限未満で主検定にできない**（旗105）。")
    print("  **秋の `遠い` 層でも反転が残るなら、春が反転しない理由は雨からの日数ではない。**")
    print(f"  層は **`直後` ≤{RECENT_MAX} 日／`遠い` ≥{REMOTE_MIN} 日**"
          f"（イベント {PRIMARY_THR:.0f} mm・旗103–105 と同一）。")

    if not a.real:
        print("\n  【合成データで検証する】**四つとも「期待する枝に到達するか」を見る**。")
        print("  **`theta_partial` で★が出るなら、交絡と本物を区別できない**——")
        print("  **その場合は★の枝を「説明できた」と読めない**（追補に書く）。")
        want = {"rain_driven": "★", "rain_free": "▲", "theta_driven": "▲",
                "theta_partial": "**★が出たら区別できない（最重要）**",
                "theta_confounded": "`遠い` が潰れて判定しない",
                "no_reversal": "門①-a で落ちる"}
        got = {}
        for k, w in want.items():
            print(f"\n  ===== 合成 `{k}` —— 期待：**{w}** =====")
            d, P = synth(k)
            got[k] = run_site(f"合成/{k}", d, P)
        print("\n  === 合成のまとめ ===")
        for k, w in want.items():
            print(f"    {k:<14}期待 {w:<16}実際 {got[k]}")
        print("\n  **rain_driven→★・rain_free→▲・theta_driven→▲・no_reversal→判定しない**")
        print("  **が揃って初めて、この道具は実データに使える。**")
        return

    verdicts = {}
    for s in SITES:
        try:
            # **道具の欠陥 #37**：第1版は `daily_energy(s, months)` と呼んでいた。
            # **`qc_max` は必須の位置引数**で、**戻り値は `(d, 年数)` の組**である
            # （旗91/93/95/97/103 はすべて `d, _ = daily_energy(s, months, qc_max)` と書いている）。
            # **合成の枝は `synth` が d と P を直接返すので、この行を一度も通らなかった。**
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max)
            P = daily_precip(s, a.qc_max)
        except Exception as e:
            print(f"\n  ━━ {s} ━━\n    読めない {type(e).__name__}: {str(e)[:90]}")
            continue
        verdicts[s] = run_site(s, d, P, AUTUMN, "秋")

    print("\n  === 集計（事前登録の判定規則に当てる）===")
    for s, v in verdicts.items():
        print(f"    {s:<10}{v if v else '判定しない'}")
    ok = [v for v in verdicts.values() if v]
    star = sum(v.startswith("★") for v in ok)
    tri = sum(v.startswith("▲") for v in ok)
    print("\n  === 結論 ===")
    if len(ok) < 2:
        print(f"  **判定しない**——門①-a を通ったサイトが {len(ok)} で 2 未満。")
    elif star >= 2:
        print("  **★手C が説明する**——**「季節」を「雨からの日数」で言い直せた。**")
    elif tri >= 2:
        print("  **▲手C は説明しない**——**`遠い` でも反転する以上、春の不反転の理由ではない。**")
    else:
        print("  **○中間**——**サイトごとに書き、まとめない。**")

    print("\n  ── 参考（**主判定には使わない**）：春の US-SRM ──")
    try:
        d_sp, _ = daily_energy("US-SRM", list(range(1, 13)), a.qc_max)
        run_site("US-SRM", d_sp, daily_precip("US-SRM", a.qc_max), SPRING, "春")
    except Exception as e:
        print(f"    読めない {type(e).__name__}: {str(e)[:80]}")

    print("\n  留保（事前登録どおり）：")
    print("   ・**交絡は解けない**——**`直後` の日は `遠い` の日より湿っている**（旗105）。")
    print("     **言えるのは「`遠い` 層で反転が残るか」までで、「雨からの日数そのものの効果」ではない。**")
    print("   ・**春は主検定に入っていない。** **春の不反転を直接説明したわけではない**——")
    print("     **説明の候補として残るか消えるかだけが決まる。**")
    print("   ・**独立クラスタは 2**（Walnut Gulch＝Wkg/Whs・Santa Rita＝SRM）。")
    print("     **「3 サイトで確かめた」とは書かない。**")
    print(f"   ・**イベントのしきい値 {PRIMARY_THR:.0f} mm は我々が決めた**（旗103）。**恣意性は残る。**")


if __name__ == "__main__":
    main()
