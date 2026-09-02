"""旗107：**層を揃えても、季節差は残るか**（手C はどこまで説明したか・事前登録 step107）。

**旗106 は「秋の中で `遠い` なら反転が消える」ことを 1 クラスタで示した。**
**だが春の不反転が説明されたかは別問題である**——
**`遠い` であることが理由のすべてなら、`遠い` 層の中では春と秋が同じに振る舞うはず。**

**測るのは差そのもの**：**Δ = r_秋 − r_春**（旗91 の `diff_boot`・年ブロック・ブート）。
**二値化しない**（旗59 の教訓：**閾値をまたがなかっただけで値は下がっていた**）。

**事前登録 step107 で固定済み**：
  ・**主軸は `遠い` 層**（春も秋も 6 通りすべてで下限を満たす・旗105）
  ・**主判定は θ→γH の Δ**（旗106 留保②：`遠い` で落ちたのは顕熱側だけだった）
  ・**規則は「クラスタ」で書く**——**旗106 で『サイト』で書いて 1 クラスタを 2 と数えた**
  ・**独立クラスタは 2**：Walnut Gulch（Wkg・Whs）／Santa Rita（SRM）
  ・**クラスタ内で割れたら、そのクラスタは判定しない**

    python research/season_within_stratum_step107.py            # 合成で検証（既定）
    python research/season_within_stratum_step107.py --real     # 実データ（/mnt/hdd）
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stratified_bowen_step89 import cell_of, MIN_DAYS, MIN_YEARS
from soiltemp_match_step90 import SPRING, AUTUMN
from downsample_autumn_step91 import diff_boot
from rain_history_probe_step103 import (rain_history, daily_precip,
                                        PRIMARY_THR, RECENT_MAX, REMOTE_MIN)
from evaporation_regime_step36 import daily_energy

CLUSTERS = {"Walnut Gulch": ("US-Wkg", "US-Whs"), "Santa Rita": ("US-SRM",)}


def prep(d, P):
    """`θ高×Rg高` セルに降雨履歴を貼る（**旗106 と同一**）。"""
    lab, _, _ = cell_of(d)
    hh = d[lab == "θ高×Rg高"]
    j = hh.join(rain_history(P, PRIMARY_THR), how="left")
    j["usable"] = j["usable"].fillna(False).astype(bool)
    return j[j["usable"]]


def pick(j, season, layer):
    m = [mo in season for mo in j.index.month]
    sub = j[m]
    return sub[sub["dry"] >= REMOTE_MIN] if layer == "遠い" else sub[sub["dry"] <= RECENT_MAX]


def show_delta(tag, sp, au):
    """**Δ = r_秋 − r_春** を出す。**日数・年数・θ の分布を必ず併記する。**"""
    ns, na = len(sp), len(au)
    ys = sp.index.year.nunique() if ns else 0
    ya = au.index.year.nunique() if na else 0
    print(f"      {tag}")
    print(f"        春 {ns:>4} 日／{ys:>2} 年  θ 中央 "
          f"{np.nanmedian(sp['th']) if ns else float('nan'):.2f}"
          f"   秋 {na:>4} 日／{ya:>2} 年  θ 中央 "
          f"{np.nanmedian(au['th']) if na else float('nan'):.2f}")
    if ns < MIN_DAYS or na < MIN_DAYS or ys < MIN_YEARS or ya < MIN_YEARS:
        print("        **下限未満**＝この層は判定しない")
        return None
    res = diff_boot(sp, au)
    if res is None:
        print("        Δ を出せない")
        return None
    out = {}
    for k, nm in (("h", "θ→γH"), ("le", "θ→γLE")):
        (ra, rs, dd), ci, nb = res[k]
        if ci is None:
            print(f"        {nm}：春 {rs:+.2f}／秋 {ra:+.2f}／Δ {dd:+.2f}（CI 出ず・{nb} 回）")
            out[k] = None; continue
        cross = ci[0] <= 0 <= ci[1]
        print(f"        {nm}：春 {rs:+.2f}／秋 {ra:+.2f}／"
              f"**Δ {dd:+.2f} [{ci[0]:+.2f},{ci[1]:+.2f}]** → "
              f"{'**0 を跨ぐ＝季節差なし**' if cross else '**0 を跨がない＝季節差が残る**'}")
        out[k] = (not cross, dd)
    return out


def run_site(tag, d, P):
    print(f"\n  ━━ {tag} ━━")
    j = prep(d, P)
    print("    ── 主軸：`遠い` 層（dryspell ≥7 日）の中で春と秋を比べる ──")
    main = show_delta("遠い", pick(j, SPRING, "遠い"), pick(j, AUTUMN, "遠い"))
    print("    ── 参考：`直後` 層（**春が下限を満たすのは US-SRM だけ**・主判定に使わない） ──")
    show_delta("直後", pick(j, SPRING, "直後"), pick(j, AUTUMN, "直後"))
    if main is None or main.get("h") is None:
        print("    → **このサイトは判定しない**")
        return None
    kept, dd = main["h"]
    print(f"    → **{'季節差が残る' if kept else '季節差なし'}**（θ→γH の Δ = {dd:+.2f}）")
    return kept


def synth(kind, years=20, seed=0):
    """**三つとも「期待する枝に到達すること」を数値で確かめる**（旗106 の規則・4 度目）。

      ・`rain_only`   —— **反転は `直後` の日にだけ起き、季節は無関係** → **Δ ≈ 0**
      ・`season_only` —— **反転は秋にだけ起き、雨の日数は無関係** → **Δ が 0 を跨がない**
      ・`both`        —— **雨と季節の両方が要る** → **`遠い` 層では Δ ≈ 0**
        **（追補で訂正）**：**当初「Δ が 0 を跨がない」と書いたのは私の推論の誤り。**
        **`遠い` 層では雨の条件が春も秋も満たされないので、どちらも反転せず差は 0 になる。**
        **＝`遠い` 層の Δ ≈ 0 は `rain_only` と `both` を区別できない。**
    """
    from precip_pressure_test_step77 import dryspell
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2003-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy()
    Rg = np.clip(150 + 120 * np.sin(2 * np.pi * (doy - 80) / 365)
                 + rng.normal(0, 25, len(idx)), 5, None)
    # **雨は旗106 の `theta_real` と同じ疎さ**（実データの層の大きさに合う）
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
    if kind == "rain_only":
        g = 1.0 * near
    elif kind == "season_only":
        g = 1.0 * is_au
    else:                                   # both
        g = 1.0 * (near & is_au)
    avail = 0.75 * Rg
    frac = np.clip(0.45 + 0.22 * g * thz + rng.normal(0, 0.05, len(idx)), .05, .95)
    d = pd.DataFrame({"th": th, "Rg": Rg, "gLE": avail * frac,
                      "gH": avail * (1 - frac)}, index=idx)
    return d, pd.Series(P, index=idx)


def main():
    ap = argparse.ArgumentParser(description="旗107：層を揃えても季節差は残るか")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()

    print("=== 旗107：層を揃えても、季節差は残るか（手C はどこまで説明したか）===")
    print("  **`遠い` 層の中で春と秋を比べる**——**`遠い` が理由のすべてなら、差は消えるはず。**")
    print("  **測るのは差そのもの Δ = r_秋 − r_春**（旗91 の年ブロック・ブート）。")
    print("  **二値化しない**（旗59：閾値をまたがなかっただけで値は下がっていた）。")
    print("  **主判定は θ→γH の Δ**（旗106 留保②）。**規則はクラスタで書く**（旗106 の反省）。")

    if not a.real:
        print("\n  【合成データで検証する】**三つとも期待する枝に到達するかを数値で見る**。")
        want = {"rain_only": "**Δ ≈ 0（跨ぐ）**",
                "season_only": "**Δ が 0 を跨がない**",
                "both": "**Δ ≈ 0（跨ぐ）**——**追補で訂正**"}
        got = {}
        for k, w in want.items():
            print(f"\n  ===== 合成 `{k}` —— 期待：{w} =====")
            d, P = synth(k)
            got[k] = run_site(f"合成/{k}", d, P)
        print("\n  === 合成のまとめ ===")
        for k, w in want.items():
            print(f"    {k:<12}期待 {w:<26}実際 "
                  f"{'季節差が残る' if got[k] else ('季節差なし' if got[k] is False else '判定しない')}")
        print("\n  **rain_only→季節差なし・season_only→残る・both→季節差なし**")
        print("  **が揃って初めて、この道具は実データに使える。**")
        print("\n  **重要（追補）**：**`rain_only` と `both` は `遠い` 層では区別できない**——")
        print("  **どちらも Δ ≈ 0 になる。** **区別するのは `直後` 層である**")
        print("  （`rain_only` は Δ ≈ 0・`both` は Δ ≠ 0）。")
        print("  **実データで `直後` 層の春が下限を満たすのは US-SRM だけ**なので、")
        print("  **Δ ≈ 0 が出ても『手C が説明しきった』とは書けない。**")
        return

    verd = {}
    for s in ("US-Wkg", "US-Whs", "US-SRM"):
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max)
            P = daily_precip(s, a.qc_max)
        except Exception as e:
            print(f"\n  ━━ {s} ━━\n    読めない {type(e).__name__}: {str(e)[:90]}")
            continue
        if P is None or P.dropna().empty:
            print(f"\n  ━━ {s} ━━\n    **降水 P が無い**"); continue
        verd[s] = run_site(s, d, P)

    print("\n  === クラスタでまとめる（**事前登録どおり**）===")
    cl = {}
    for name, members in CLUSTERS.items():
        vals = [verd.get(m) for m in members if verd.get(m) is not None]
        if not vals:
            cl[name] = None; print(f"    {name:<14}判定しない（サイトが 0）"); continue
        if len(set(vals)) > 1:
            cl[name] = None
            print(f"    {name:<14}**判定しない（クラスタ内で割れた）**："
                  f"{ {m: verd.get(m) for m in members} }")
            continue
        cl[name] = vals[0]
        print(f"    {name:<14}{'**季節差が残る**' if vals[0] else '季節差なし'}"
              f"（{len(vals)} サイト一致）")

    ok = [v for v in cl.values() if v is not None]
    print("\n  === 結論 ===")
    if len(ok) < 2:
        print(f"  **判定しない**——判定できたクラスタが {len(ok)} で 2 未満。")
    elif all(not v for v in ok):
        print("  **○雨からの日数と整合する**——**`遠い` 層では季節差が無い。**")
        print("  **ただし『雨かつ季節の両方が要る』世界とも整合し、この層では区別できない**")
        print("  **（追補：合成 `both` も Δ ≈ 0 を返した）。**")
        print("  **『手C が季節差を説明しきった』とは書かない。**")
    elif all(ok):
        print("  **▲手C は季節差を残す**——**説明は部分的。旗106 の★を割り引く。**")
    else:
        print("  **○部分的**——**クラスタごとに書き、まとめない。**")

    print("\n  留保（事前登録どおり）：")
    print("   ・**`遠い` に揃えても、春と秋は他のすべてで違う**（日長・気温・植生・前年の履歴）。")
    print("     **Δ ≠ 0 は『雨以外の何かが残る』までしか言わない。** **それが何かは言えない。**")
    print("   ・**春の `遠い` は秋の 2〜3 倍の日数**（502 対 173 ほか）＝**CI の幅が季節で違う。**")
    print("   ・**独立クラスタは 2 しかない。** **どの結論も弱い。強い言葉を使わない。**")
    print("   ・**`直後` 層の春は US-SRM だけ**＝**参考にしかならない。**")


if __name__ == "__main__":
    main()
