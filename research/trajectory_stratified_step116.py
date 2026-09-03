"""旗116：**手B——θ の変化の向き（軌道）で季節差は説明されるか**（事前登録 step116）。

## 問い（二段）

1. **門①-a**：各クラスタで秋全体（θ高×Rg高）に Bowen 反転が起きるか（旗106/109/111/113 と同一）。
2. **主判定**：反転するクラスタで、**Δθ の符号（上昇/下降）を揃えても季節差（春 vs 秋）が残るか。**

**「軌道で説明される」＝両方の Δθ 層で季節差が消える**（|Δ_H| < DELTA_FLOOR または CI が 0 を跨ぐ）。
**片層だけ消える → ▲部分的。** **両層で残る → ▲軌道でも説明されない。**

**標的は Santa Rita**（旗107 で雨を揃えても残った `or` 型の季節差）。

## 縛り（事前登録 step116）

- **層は 2・クラスタは 3。片層だけ差なしを拾って「軌道で説明される」と書かない（規則どおり ▲部分的）。**
- **門②**：各層の春・秋の θ 中央値を併記し、**Δθ の層別が θ の水準を作り直していないか**を確認する。
- **統計量・しきい値・下限・DELTA_FLOOR・diff_boot を一つも変えない**（旗89/107/109 と同一）。

    python research/trajectory_stratified_step116.py           # 合成で検証（既定）
    python research/trajectory_stratified_step116.py --real    # 実データ（/mnt/hdd）
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
from humid_forest_pattern_step109 import rev, delta, DELTA_FLOOR
from change_rate_probe_step115 import delta_theta, WINDOW, CLUSTER, daily_energy
from evaporation_regime_step36 import daily_energy as _de  # noqa: F401  (実データ経路の明示)
from runlog import tee_stdout

SITES = ("US-Wkg", "US-Whs", "US-SRM", "CN-Du2")   # 旗115 で 4 群が下限を満たした 3 クラスタ
SIGNS = ("上昇", "下降")   # Δθ>0 / Δθ<0


def prep_b(d: pd.DataFrame) -> pd.DataFrame:
    """θ高×Rg高 セルの日に、暦日ベースの Δθ 符号を付ける。"""
    lab, _, _ = cell_of(d)
    j = d[lab == "θ高×Rg高"].copy()
    dth = delta_theta(d)                 # 全 d 上で暦日差、d.index に整列
    j["dth"] = dth.reindex(j.index)
    return j.dropna(subset=["dth"])


def pick_b(j: pd.DataFrame, season, sign: str) -> pd.DataFrame:
    sub = j[[m in season for m in j.index.month]]
    return sub[sub["dth"] > 0] if sign == "上昇" else sub[sub["dth"] < 0]


def _thmed(sub) -> str:
    return f"{np.nanmedian(sub['th']):.3f}" if len(sub) else "—"


def run_site(tag: str, d: pd.DataFrame):
    print(f"\n  ━━ {tag}（{CLUSTER.get(tag, '?')}）━━")
    j = prep_b(d)
    # 門①-a：秋全体で反転するか
    print("    ── 門①-a：秋全体（θ高×Rg高）で反転するか ──")
    au_all = j[[m in AUTUMN for m in j.index.month]]
    if not rev("秋全体", au_all):
        print("    → **判定しない**（**門①-a を通らない**）")
        return None
    # 門②の基準：セル全体の春秋 θ 中央値の差
    sp_all = j[[m in SPRING for m in j.index.month]]
    base_gap = abs(np.nanmedian(au_all["th"]) - np.nanmedian(sp_all["th"])) \
        if len(sp_all) and len(au_all) else float("nan")
    print(f"    （門②基準：セル全体 θ 中央 春 {_thmed(sp_all)}／秋 {_thmed(au_all)}／差 {base_gap:.3f}）")

    # 主判定：各 Δθ 層の中で Δ = r_秋 − r_春
    diffs = {}
    for sign in SIGNS:
        sp = pick_b(j, SPRING, sign)
        au = pick_b(j, AUTUMN, sign)
        print(f"    ── Δθ {sign} 層の季節差（θ→γH が主）──")
        print(f"        θ 中央：春 {_thmed(sp)}／秋 {_thmed(au)}")
        # 門②：層内の春秋 θ 差が基準より大きければ、θ 帯の作り直しを疑う
        if len(sp) and len(au):
            gap = abs(np.nanmedian(au["th"]) - np.nanmedian(sp["th"]))
            if np.isfinite(base_gap) and gap > base_gap + 0.02:
                print(f"        **門②警告**：層内の春秋 θ 差 {gap:.3f} > 基準 {base_gap:.3f}"
                      f"＝Δθ の層別が θ の水準を作り直している疑い")
        diffs[sign] = delta(f"{sign}", sp, au)
    return diffs


def classify(diffs: dict) -> str:
    """判定規則（事前登録どおり）。diffs[sign] は True=差あり／False=差なし／None=測れない。"""
    vals = [diffs.get(s) for s in SIGNS]
    有差 = [s for s in SIGNS if diffs.get(s) is True]
    無差 = [s for s in SIGNS if diffs.get(s) is False]
    測れず = [s for s in SIGNS if diffs.get(s) is None]
    if 測れず:
        return f"判定不能（{'・'.join(測れず)} 層が測れない）"
    if not 有差:
        return "○軌道で説明される（両層とも季節差なし）"
    if len(有差) == 1:
        return f"▲部分的（{有差[0]} 層でのみ季節差が残る）"
    return "▲説明されない（両層で季節差が残る）"


# ───────────────────────── 合成 ─────────────────────────

def synth(kind: str, years: int = 22, seed: int = 0):
    """反転の強さ k をどう置くかで、判定規則の各枝に到達させる。

    - not_explained : k は季節で決まる（秋>春）。Δθ 層に依らない → 両層で季節差。
    - explained     : k は Δθ 符号で決まる（上昇>下降）。季節に依らない → 両層で季節差なし。
    - partial       : k は (秋 かつ 上昇) でだけ高い → 上昇層でのみ季節差。
    - no_reversal   : k=0 → 門①-a を通らない。
    - theta_level   : k は θ 水準で決まり、Δθ 符号が θ 水準と強く相関 → 門②が警告を出す。
    """
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2001-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy()
    Rg = np.clip(150 + 120 * np.sin(2 * np.pi * (doy - 80) / 365)
                 + rng.normal(0, 20, len(idx)), 5, None)
    # θ：ゆっくりした軌道（両符号の Δθ が出る）＋弱い季節（秋にモンスーンの山）＋雨的パルス
    slow = (pd.Series(rng.normal(0, 1, len(idx))).rolling(20, min_periods=1)
            .mean().to_numpy())
    slow = slow / (np.std(slow) + 1e-12)
    monsoon = 0.05 * np.exp(-0.5 * ((doy - 250) / 35.) ** 2)
    pulse = np.where(rng.random(len(idx)) < 0.12, rng.gamma(1.2, 0.03, len(idx)), 0.0)
    th = np.clip(0.28 + 0.06 * slow + monsoon + pulse
                 + rng.normal(0, 0.010, len(idx)), .05, .9)
    is_au = np.array([m in AUTUMN for m in idx.month])
    dth0 = th - pd.Series(th, index=idx).shift(WINDOW).to_numpy()
    up0 = np.nan_to_num(dth0, nan=0.0) > 0
    if kind == "theta_reband":
        # Δθ 上昇層の秋にだけ θ を大きく持ち上げ、層内の春秋 θ 差を基準より広げる
        # ＝Δθ の層別が θ の帯を作り直している状況（門②が警告を出すべき検査用）
        th = np.clip(th + 0.12 * (is_au & up0) + rng.normal(0, 0.005, len(idx)), .05, .9)
    thz = (th - th.mean()) / (th.std() + 1e-12)
    dth = th - pd.Series(th, index=idx).shift(WINDOW).to_numpy()
    up = np.nan_to_num(dth, nan=0.0) > 0

    KHI, KLO = 0.42, 0.06
    if kind == "no_reversal":
        k = np.zeros(len(idx))
    elif kind == "not_explained":
        k = np.where(is_au, KHI, KLO)
    elif kind == "explained":
        k = np.where(up, KHI, KLO)
    elif kind == "partial":
        k = np.where(is_au & up, KHI, KLO)
    elif kind == "theta_reband":
        k = np.where(is_au, KHI, KLO)   # 反転は起きる。門②が θ 帯の作り直しを警告するかを見る
    else:
        raise ValueError(kind)

    avail = 0.75 * Rg
    frac = np.clip(0.45 + k * thz + rng.normal(0, 0.04, len(idx)), .05, .95)
    d = pd.DataFrame({"th": th, "Rg": Rg, "gLE": avail * frac,
                      "gH": avail * (1 - frac), "Ta": 15 + rng.normal(0, 3, len(idx))},
                     index=idx)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", action="store_true", help="実データ（/mnt/hdd）")
    ap.add_argument("--only", default=None)
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()
    tee_stdout("step116")

    print("=== 旗116：手B（θ の変化の向きで季節差は説明されるか）===")
    print(f"  Δθ = θ(t) − θ(t−{WINDOW}日)。符号で層別。**軌道で説明される＝両層で季節差が消える。**")
    print(f"  下限 {MIN_DAYS} 日・{MIN_YEARS} 暦年／DELTA_FLOOR {DELTA_FLOOR}／θ高×Rg高 セル。")

    if not a.real:
        print("\n  ##### 合成検証：各枝が期待どおり分かれるか #####")
        expect = {"not_explained": "▲説明されない（両層で季節差が残る）",
                  "explained": "○軌道で説明される（両層とも季節差なし）",
                  "partial": "▲部分的（上昇 層でのみ季節差が残る）",
                  "no_reversal": "判定しない（門①-a を通らない）",
                  "theta_reband": "（門②が警告を出すか）"}
        ok_all = True
        for kind, exp in expect.items():
            print(f"\n  ===== 合成 `{kind}` =====")
            diffs = run_site(f"合成/{kind}", synth(kind))
            got = "判定しない（門①-a を通らない）" if diffs is None else classify(diffs)
            if kind == "theta_reband":
                print(f"    → 得：{got}（**門②の警告が上に出ていれば合格**）")
                continue
            mark = "✔" if got == exp else "✘"
            ok_all = ok_all and (got == exp)
            print(f"    → 期待 {exp}／実際 {got}  {mark}")
        print(f"\n  → **4 枝（門①落ち含む）すべて一致：{ok_all}**")
        return

    sites = SITES
    if a.only:
        want = {x.strip() for x in a.only.split(",") if x.strip()}
        sites = tuple(s for s in SITES if s in want) or tuple(sorted(want))
    print(f"\n  実データ：{', '.join(sites)}")
    verdict = {}
    for s in sites:
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max)
        except Exception as e:                # noqa: BLE001
            print(f"\n  ━━ {s} ━━  読み込み失敗：{type(e).__name__}: {e}")
            verdict[s] = None
            continue
        diffs = run_site(s, d)
        verdict[s] = None if diffs is None else classify(diffs)

    print("\n  === クラスタごとの結論 ===")
    by_cl: dict[str, list] = {}
    for s, v in verdict.items():
        by_cl.setdefault(CLUSTER.get(s, s), []).append((s, v))
    for cl, items in by_cl.items():
        print(f"    {cl}：")
        for s, v in items:
            print(f"      {s}: {v if v is not None else '門①-a を通らない/測れない'}")


if __name__ == "__main__":
    main()
