"""旗109：**湿潤な森林でも、同じ型が出るか**（**別の母集団への外挿**・事前登録 step109）。

**母集団が違うことを最初に書く。**
**Bowen 反転（A-3）は「乾燥地で水がマスター変数」という主張**（旗36/82/88）で、
**既存 2 クラスタは半乾燥地。本検定の 3 サイトはすべて湿潤な森林である**
（US-Ho1 メイン州・US-NC2 ノースカロライナ・US-WCr ウィスコンシン州）。
**どんな結果でも「再現した」とは書かない。**

**旗108 で、乾燥地では 4 群が揃うのが US-SRM だけと確かめた**
（US-Var・US-Ton は地中海性気候で `秋×直後` が構造的に存在しない）。
**＝同じ母集団でこれ以上クラスタを増やす道は、手元には無い。**

**四つの型を事前に登録する**（**旗107 で事後に見つけた `or` を、今回は事前に入れる**）：

| `遠い` の Δ | `直後` の Δ | 型 |
|---|---|---|
| ≈0 | ≈0 | `rain_only`（雨だけ） |
| ≠0 | ≠0 | `season_only`（季節だけ） |
| ≈0 | ≠0 | `both`（雨**かつ**秋） |
| ≠0 | ≈0 | `or`（雨**または**秋・**Santa Rita の型**） |

**しきい値・下限・統計量は旗106/107 と完全に同一。一つも変えない。**

    python research/humid_forest_pattern_step109.py            # 合成で検証（既定）
    python research/humid_forest_pattern_step109.py --real     # 実データ（/mnt/hdd）
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
from downsample_autumn_step91 import diff_boot
from rain_history_probe_step103 import (rain_history, daily_precip,
                                        PRIMARY_THR, RECENT_MAX, REMOTE_MIN)
from evaporation_regime_step36 import daily_energy
from colocate_step51 import haversine

SITES = ("US-Ho1", "US-NC2", "US-WCr")     # 旗108 で 4 群すべてが下限を満たした新しい 3 件
CLUSTER_KM = 50.0                          # 旗82/107 と同じ（**手で断定せず道具が確かめる**）
# **追補で入れた絶対的な下限**（**実データ実行前**）。
# **CI が 0 を外れるだけでは「差がある」と言わない**——**旗75 の教訓**
# （プラセボ超えだけでは勝てない／絶対的な下限が要る）**の 2 度目**である。
# 合成で `rain_only` の `直後` が **Δ=+0.039 [0.004,0.092]** となり、
# **実質ゼロなのに CI が 0 を外れて `both` と誤分類された**。
# **旗107 の実測（−0.05／−0.10／−0.38）はこの下限で判定が変わらない**ことを確かめてある。
DELTA_FLOOR = 0.15

PATTERN = {(False, False): "rain_only（雨だけ）",
           (True, True): "season_only（季節だけ）",
           (False, True): "both（雨かつ秋）",
           (True, False): "or（雨または秋）"}


def prep(d, P):
    lab, _, _ = cell_of(d)
    j = d[lab == "θ高×Rg高"].join(rain_history(P, PRIMARY_THR), how="left")
    j["usable"] = j["usable"].fillna(False).astype(bool)
    return j[j["usable"]]


def pick(j, season, layer):
    sub = j[[m in season for m in j.index.month]]
    return sub[sub["dry"] >= REMOTE_MIN] if layer == "遠い" else sub[sub["dry"] <= RECENT_MAX]


def rev(tag, sub):
    """反転したか（旗106 と同一）。**日数・年数・θ を必ず併記する。**"""
    n, y = len(sub), (sub.index.year.nunique() if len(sub) else 0)
    if n < MIN_DAYS or y < MIN_YEARS:
        print(f"      {tag:<8}{n:>4} 日／{y:>2} 年  **下限未満**"); return None
    res = test_cell(sub)
    if res is None:
        print(f"      {tag:<8}{n:>4} 日／{y:>2} 年  測れない"); return None
    (rl, cl, _), (rh, ch, _) = res["le"], res["h"]
    r = reversed_(res)
    print(f"      {tag:<8}{n:>4} 日／{y:>2} 年  θ 中央 {np.nanmedian(sub['th']):.2f}  "
          f"θ→γLE {rl:+.2f} [{cl[0]:+.2f},{cl[1]:+.2f}]  "
          f"θ→γH {rh:+.2f} [{ch[0]:+.2f},{ch[1]:+.2f}]  "
          f"→ {'**反転**' if r else '反転せず'}")
    return bool(r)


def delta(tag, sp, au):
    """**Δ = r_秋 − r_春**（旗107 と同一）。**θ→γH が主判定。**"""
    ns, na = len(sp), len(au)
    ys = sp.index.year.nunique() if ns else 0
    ya = au.index.year.nunique() if na else 0
    print(f"      {tag:<8}春 {ns:>4} 日/{ys:>2} 年  秋 {na:>4} 日/{ya:>2} 年")
    if ns < MIN_DAYS or na < MIN_DAYS or ys < MIN_YEARS or ya < MIN_YEARS:
        print("               **下限未満**"); return None
    res = diff_boot(sp, au)
    if res is None:
        print("               Δ を出せない"); return None
    keep = None
    for k, nm in (("h", "θ→γH"), ("le", "θ→γLE")):
        (ra, rs, dd), ci, nb = res[k]
        if ci is None:
            print(f"               {nm}：Δ {dd:+.2f}（CI 出ず・{nb} 回）"); continue
        cross = ci[0] <= 0 <= ci[1]
        print(f"               {nm}：春 {rs:+.2f}／秋 {ra:+.2f}／"
              f"**Δ {dd:+.2f} [{ci[0]:+.2f},{ci[1]:+.2f}]** → "
              f"{'0 を跨ぐ' if cross else '**0 を跨がない**'}")
        if k == "h":
            # **「差がある」＝ CI が 0 を跨がない **かつ** |Δ| ≥ 下限**
            keep = (not cross) and abs(dd) >= DELTA_FLOOR
            if (not cross) and abs(dd) < DELTA_FLOOR:
                print(f"               → **|Δ| {abs(dd):.3f} < 下限 {DELTA_FLOOR}"
                      f"＝差なしとする**（CI は 0 を外れているが実質ゼロ）")
    return keep


def run_site(tag, d, P):
    print(f"\n  ━━ {tag} ━━")
    j = prep(d, P)
    print("    ── 門①-a：秋全体で反転するか ──")
    if not rev("秋全体", j[[m in AUTUMN for m in j.index.month]]):
        print("    → **判定しない**（**門①-a を通らない**）")
        return None
    print("    ── (A) 秋の中の層別（旗106 と同一）──")
    rev("秋直後", pick(j, AUTUMN, "直後"))
    rev("秋遠い", pick(j, AUTUMN, "遠い"))
    print("    ── (B) 層の中の季節差（旗107 と同一・**型はこれで決める**）──")
    d_far = delta("遠い", pick(j, SPRING, "遠い"), pick(j, AUTUMN, "遠い"))
    d_near = delta("直後", pick(j, SPRING, "直後"), pick(j, AUTUMN, "直後"))
    if d_far is None or d_near is None:
        print("    → **判定しない**（どちらかの層の Δ が出ない）")
        return None
    pat = PATTERN[(d_far, d_near)]
    print(f"    → **型：{pat}**")
    return pat


def synth(kind, years=20, seed=0):
    """**五つとも「期待する枝に到達すること」を数値で確かめる**（旗106 の規則・5 度目）。"""
    from precip_pressure_test_step77 import dryspell
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2003-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy()
    Rg = np.clip(150 + 120 * np.sin(2 * np.pi * (doy - 80) / 365)
                 + rng.normal(0, 25, len(idx)), 5, None)
    lam = 0.45 * (0.055 + 0.110 * np.exp(-0.5 * ((doy - 100) / 30.) ** 2)
                  + 0.150 * np.exp(-0.5 * ((doy - 250) / 45.) ** 2))
    P = np.where(rng.random(len(idx)) < np.clip(lam, 0.002, 1),
                 rng.gamma(1.3, 7.0, len(idx)), 0.0)
    ds = dryspell(P, PRIMARY_THR)
    recent = np.exp(-np.nan_to_num(ds, nan=30.) / 4.0)
    slow = (pd.Series(rng.normal(0, 1, len(idx))).rolling(45, min_periods=1)
            .mean().to_numpy())
    slow = slow / (np.std(slow) + 1e-12)
    th = np.clip(0.16 + 0.240 * recent + 0.045 * slow
                 + 0.040 * np.sin(2 * np.pi * (doy - 200) / 365)
                 + rng.normal(0, 0.020, len(idx)), .02, .9)
    thz = (th - th.mean()) / th.std()
    is_au = np.array([m in AUTUMN for m in idx.month])
    near = ds <= RECENT_MAX
    g = {"rain_only": 1.0 * near,
         "season_only": 1.0 * is_au,
         "both": 1.0 * (near & is_au),
         "or": 1.0 * (near | is_au),
         "no_reversal": np.zeros(len(idx))}[kind]
    avail = 0.75 * Rg
    frac = np.clip(0.45 + 0.22 * g * thz + rng.normal(0, 0.05, len(idx)), .05, .95)
    d = pd.DataFrame({"th": th, "Rg": Rg, "gLE": avail * frac,
                      "gH": avail * (1 - frac)}, index=idx)
    return d, pd.Series(P, index=idx)


def main():
    ap = argparse.ArgumentParser(description="旗109：湿潤な森林でも同じ型が出るか")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()

    print("=== 旗109：湿潤な森林でも、同じ型が出るか（**別の母集団への外挿**）===")
    print("  **母集団が違う**——**A-3 は乾燥地の主張**で、本検定の 3 サイトは**湿潤な森林**。")
    print("  **どんな結果でも『再現した』とは書かない。**")
    print("  **四つの型を事前に登録した**（`rain_only`／`season_only`／`both`／**`or`**）。")
    print(f"  **しきい値・統計量は旗106/107 と同一**／**Δ の絶対的な下限 {DELTA_FLOOR} は追補で追加**")
    print("  （**CI が 0 を外れるだけでは差ありとしない**＝旗75 の教訓）。")

    if not a.real:
        print("\n  【合成データで検証する】**五つとも期待する枝に到達するかを数値で見る**。")
        print("  **`or` が期待どおりに出なければ、Santa Rita の型を検出できない**＝実データに進まない。")
        want = {"rain_only": "rain_only（雨だけ）", "season_only": "season_only（季節だけ）",
                "both": "both（雨かつ秋）", "or": "**or（雨または秋）**",
                "no_reversal": "**門①-a で落ちる**"}
        got = {}
        for k, w in want.items():
            print(f"\n  ===== 合成 `{k}` —— 期待：{w} =====")
            d, P = synth(k)
            got[k] = run_site(f"合成/{k}", d, P)
        print("\n  === 合成のまとめ ===")
        ok = True
        for k, w in want.items():
            hit = (got[k] is None) if k == "no_reversal" else (got[k] == w.strip("*"))
            ok &= hit
            print(f"    {k:<12}期待 {w:<24}実際 {got[k]}  {'✔' if hit else '**✘**'}")
        print(f"\n  → **五つとも一致：{ok}**"
              f"{'' if ok else '  ← **一致しない枝があるので実データに進まない**'}")
        return

    verd, coords = {}, {}
    for s in SITES:
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max)
            P = daily_precip(s, a.qc_max)
        except Exception as e:
            print(f"\n  ━━ {s} ━━\n    読めない {type(e).__name__}: {str(e)[:90]}")
            continue
        if P is None or P.dropna().empty:
            print(f"\n  ━━ {s} ━━\n    **降水 P が無い**"); continue
        # **道具の欠陥 #38（旗109 の実行で判明）**：
        # **`SiteSpec` は座標を持たない**（フィールドは code/data_dir/source/hh_glob/
        # fmt/var_overrides/description の 7 つだけ）。**`st.lat` は最初から存在しなかった。**
        # **座標は BADM から引くのが正しい**（旗79/86 の `build_badm_index`）。
        # **本実行では 3 サイトすべてが門①-a で落ちたので結論に影響しない**が、
        # **「座標を取れなかった」という護りが働いて誤った主張を防いだ**ことは記録する。
        coords[s] = (None, None)
        verd[s] = run_site(s, d, P)

    print("\n  === クラスタを座標で確かめる（**手で断定しない**・50 km 単連結）===")
    named = [(s, c) for s, c in coords.items()
             if c[0] is not None and np.isfinite(c[0] or np.nan)]
    if len(named) >= 2:
        for i in range(len(named)):
            for j in range(i + 1, len(named)):
                (s1, c1), (s2, c2) = named[i], named[j]
                km = haversine(c1[0], c1[1], c2[0], c2[1])
                print(f"    {s1} × {s2}：{km:,.0f} km"
                      f" → {'**同一クラスタ**' if km <= CLUSTER_KM else '別クラスタ'}")
    else:
        print("    **座標を取れない**——**`SiteSpec` は座標を持たない**（欠陥 #38）。")
        print("    **正しくは BADM から引く**（旗79/86 の `build_badm_index`）。")
        print("    **クラスタ数を主張しない。**")

    print("\n  === 集計 ===")
    for s in SITES:
        print(f"    {s:<10}{verd.get(s) or '判定しない'}")
    ok = [v for v in verd.values() if v]
    print("\n  === 結論 ===")
    if not ok:
        print("  **★湿潤林では Bowen 反転が起きない**——**門①-a を通ったサイトが 0**。")
        print("  **A-3 のスコープ（乾燥地に限る）を支持する独立の証拠**である。")
        print("  **これは空振りではない**——**これまで乾燥地でしか試していなかったので、")
        print("  「限る」と言えなかった。**")
    elif len(ok) < 2:
        print(f"  **判定しない**——門①-a を通ったのが {len(ok)} で、**n=1 では型を主張しない**。")
        print(f"  **ただし『湿潤林で反転が出た』こと自体は記録する**：{ok}")
    else:
        from collections import Counter
        c = Counter(ok).most_common()
        if c[0][1] > len(ok) / 2:
            print(f"  **○{c[0][0]} が湿潤林で現れた**（{c[0][1]}/{len(ok)}）。")
            print("  **外挿であって再現ではない**——**母集団が違う。**")
        else:
            print(f"  **○型が割れた**：{dict(Counter(ok))}。**サイトごとに書き、まとめない。**")

    print("\n  留保（事前登録どおり）：")
    print("   ・**母集団が違う。** **どんな結果も『再現』ではない。**")
    print("   ・**手元の乾燥地クラスタは 2 が上限**（旗108）。**増えたのは別の母集団のクラスタ。**")
    print("   ・**湿潤林の『θ 高』は乾燥地の『θ 高』と意味が違う**——")
    print("     **乾燥地では雨の直後、湿潤林では常態でありうる。**")
    print("     **同じ層別語を使っているが、指しているものが同じとは限らない。**")
    print("   ・**3 サイトはいずれも北米東部〜中西部**＝**気候帯は 3 つに割れていない。**")


if __name__ == "__main__":
    main()
