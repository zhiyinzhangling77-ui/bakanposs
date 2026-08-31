"""旗96：**春に Bowen 反転が起きないのは、飽差（VPD）が高いからか**（事前登録 step96）。

旗95 で「θ・Rg・Ts・GCC の四つを試して打ち切る」と書いたのは**早すぎた**——
**`VPD` は旗36・89・90・95 のどこにも現れず、一度も層別に使っていなかった**。
**しかも Bowen 比にとって最も直接的な量**（**気孔を閉じさせるのは飽差**）。

**事前登録 step96 で固定済み**：
  ・場所・帯の作り方・統計量は**旗89/90/95 と同一**。**作り直さない。**
  ・**門①（対照）**：**帯の中で秋が反転しなければ、そのサイトは判定しない**
    ——**旗95 で秋を「参考」にしたまま US-Wkg で対照ごと落ちた**（欠陥31）。**今度は門にする。**
  ・**門②（分離可能性）**：**セル内で `r(VPD, Rg) > 0.7` なら判定しない**
    ——**VPD は `es(Ta)(1−RH/100)` で `Ta` は `Rg` と相関する**（旗32）。
    **VPD の帯が Rg の帯を作り直しているだけ**なら意味がない。**閾値 0.7 は固定。**

    python research/vpd_match_step96.py                    # 合成で検証（既定）
    python research/vpd_match_step96.py --real
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaporation_regime_step36 import daily_energy, _fmt
from stratified_bowen_step89 import cell_of, test_cell, reversed_, MIN_DAYS, MIN_YEARS
from soiltemp_match_step90 import band, SPRING, AUTUMN

R_MAX = 0.70          # 門②の閾値（事前登録で固定）


def _rank(x):
    x = np.asarray(x, float)
    return pd.Series(x).rank().to_numpy()


def spearman(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 20:
        return np.nan
    return float(np.corrcoef(_rank(a[ok]), _rank(b[ok]))[0, 1])


def run_site(tag, d, verbose=True):
    """1 サイト。**門②→帯→検定→門①** の順に見る。"""
    if "VPD" not in d.columns or d["VPD"].notna().sum() < MIN_DAYS:
        print(f"    **VPD が無い／足りない**＝判定しない"); return None, "VPD無し"
    lab, tmed, rmed = cell_of(d)
    hh = d[(lab == "θ高×Rg高") & d["VPD"].notna()]
    sp = hh[[m in SPRING for m in hh.index.month]]
    au = hh[[m in AUTUMN for m in hh.index.month]]
    print(f"    θ={tmed:.3f}／Rg={rmed:.1f}（旗89 と同一）"
          f"／θ高×Rg高：春 {len(sp)} 日・秋 {len(au)} 日")
    if len(sp) == 0 or len(au) == 0:
        print("    片方が空＝判定しない"); return None, "片方が空"

    # ── 門②：VPD と Rg を分離できるか ──
    r_vr = spearman(hh["VPD"].to_numpy(), hh["Rg"].to_numpy())
    print(f"    **門②**：セル内 r(VPD, Rg) = {r_vr:+.2f}"
          f"（閾値 {R_MAX:.2f}）")
    if not np.isfinite(r_vr) or abs(r_vr) > R_MAX:
        print(f"      → **VPD と Rg を分離できない＝判定しない**")
        return None, "門②（VPD≈Rg）"

    lo, hi = band(sp["VPD"].to_numpy(), au["VPD"].to_numpy())
    print(f"    VPD 中央値：春 {sp['VPD'].median():.3f}／秋 {au['VPD'].median():.3f}")
    if hi <= lo:
        print(f"    → **VPD の帯が作れない**（[{lo:.3f}, {hi:.3f}]）"
              f"＝**春と秋は VPD で完全に分離**＝判定しない")
        return None, "帯なし"
    spb = sp[(sp["VPD"] >= lo) & (sp["VPD"] <= hi)]
    aub = au[(au["VPD"] >= lo) & (au["VPD"] <= hi)]
    print(f"    **VPD の帯 [{lo:.3f}, {hi:.3f}]**（幅 {hi-lo:.3f}）"
          f"／帯の中：春 {len(spb)} 日・秋 {len(aub)} 日")
    for nm, sub in (("春", spb), ("秋", aub)):
        if len(sub):
            print(f"      {nm}の帯内 中央値：VPD {sub['VPD'].median():.3f}／"
                  f"Rg {sub['Rg'].median():.1f}／"
                  f"Ts {sub['Ts'].median():.1f}／θ {sub['th'].median():.3f}"
                  if "Ts" in sub.columns else
                  f"      {nm}の帯内 中央値：VPD {sub['VPD'].median():.3f}／"
                  f"Rg {sub['Rg'].median():.1f}／θ {sub['th'].median():.3f}")
    if len(spb) and len(aub):
        drg = abs(spb["Rg"].median() - aub["Rg"].median())
        iqr = float(hh["Rg"].quantile(.75) - hh["Rg"].quantile(.25))
        if iqr > 0 and drg > 0.5 * iqr:
            print(f"      ※**帯の中で Rg が揃っていない**（差 {drg:.1f}／IQR {iqr:.1f}）"
                  f"＝**読み方を弱める**（門にはしない・事前登録どおり）")

    res = {}
    for nm, sub in (("春", spb), ("秋", aub)):
        r = test_cell(sub)
        if r is None:
            print(f"      {nm}（帯の中）：日数 {len(sub)}／年 "
                  f"{sub.index.year.nunique() if len(sub) else 0}＝**下限未満**")
            res[nm] = None; continue
        rev = reversed_(r); res[nm] = rev
        print(f"      {nm}（帯の中）：日数 {r['n']:>4}／年 {r['yrs']:>3}  "
              f"{_fmt(r['le'])}  {_fmt(r['h'])}  "
              f"{'**Bowen反転**' if rev else '反転せず'}")
    if res.get("秋") is None or res.get("春") is None:
        print("      → **判定しない**（下限未満）"); return None, "下限未満"
    # ── 門①：対照（秋）が働いているか ──
    if not res["秋"]:
        print("      → **門①：帯の中で秋が反転しない＝対照が働いていない＝判定しない**")
        return None, "門①（対照が落ちた）"
    verdict = "explained" if res["春"] else "not_explained"
    print(f"      → {'**帯の中では春も反転する**＝**VPD で説明された**' if res['春'] else '**帯の中でも春は反転しない**＝**VPD でも説明されない**'}")
    sub_verdict = within_spring(sp)
    if verdict == "explained" and sub_verdict not in ("low_only", None):
        print("      ※**副検定と食い違う**——**explained は VPD の手柄とは限らない。両方を書く。**")
    return verdict, None


def within_spring(sp):
    """**追補の副検定**：**春の中だけ**で VPD を高低に分け、**θ を揃えて**比べる。

    **主検定は「VPD を揃えると春秋の差が消えるか」しか答えない**——
    合成 `wet_overlap`（θ が真因・VPD は無関係）でも **explained** が出た。
    **「VPD が効いているか」は、春の中で VPD を動かして初めて分かる。**
    **「低 VPD 日は湿った日」という交絡は、θ の重なり帯で潰す。**
    """
    print("    ── **副検定：春の中だけで VPD を高低に分ける**（追補）──")
    vm = float(sp["VPD"].median())
    lo_g, hi_g = sp[sp["VPD"] <= vm], sp[sp["VPD"] > vm]
    if len(lo_g) < 10 or len(hi_g) < 10:
        print("      片方が少なすぎる＝判定しない"); return None
    lo, hi = band(lo_g["th"].to_numpy(), hi_g["th"].to_numpy())
    print(f"      VPD 中央値 {vm:.3f}／低 VPD {len(lo_g)} 日・高 VPD {len(hi_g)} 日")
    print(f"      θ 中央値：低 VPD {lo_g['th'].median():.3f}／高 VPD {hi_g['th'].median():.3f}")
    if hi <= lo:
        print(f"      → **θ の帯が作れない**（[{lo:.3f}, {hi:.3f}]）"
              f"＝**低 VPD と湿りを分離できない＝判定しない**")
        return None
    lb = lo_g[(lo_g["th"] >= lo) & (lo_g["th"] <= hi)]
    hb = hi_g[(hi_g["th"] >= lo) & (hi_g["th"] <= hi)]
    print(f"      **θ の帯 [{lo:.3f}, {hi:.3f}]**／帯の中：低 VPD {len(lb)}・高 VPD {len(hb)}")
    r = {}
    for nm, sub in (("低VPDの春", lb), ("高VPDの春", hb)):
        t = test_cell(sub)
        if t is None:
            print(f"      {nm}：日数 {len(sub)}／年 "
                  f"{sub.index.year.nunique() if len(sub) else 0}＝**下限未満**")
            r[nm] = None; continue
        rev = reversed_(t); r[nm] = rev
        print(f"      {nm}：日数 {t['n']:>4}／年 {t['yrs']:>3}  "
              f"{_fmt(t['le'])}  {_fmt(t['h'])}  "
              f"{'**Bowen反転**' if rev else '反転せず'}")
    if r.get("低VPDの春") is None or r.get("高VPDの春") is None:
        print("      → **判定しない**"); return None
    if r["低VPDの春"] and not r["高VPDの春"]:
        print("      → **低 VPD の春だけ反転**＝**VPD が効いている**"); return "low_only"
    if r["低VPDの春"] and r["高VPDの春"]:
        print("      → **どちらも反転**＝**VPD は効いていない**（春でも反転する日がある）")
        return "both"
    if not r["低VPDの春"] and not r["高VPDの春"]:
        print("      → **どちらも反転しない**＝**春では VPD を下げても反転しない**")
        return "neither"
    print("      → **高 VPD の春だけ反転**＝**想定と逆**（解釈しない）"); return "high_only"


def synth(kind, years=16, seed=0):
    """場合を作る。**`collinear` は門②の試験**である。

    **`*_overlap` は作為的に日々のばらつきを大きくした検証用**——
    **現実的な季節振幅（春 1.7 / 秋 0.5）では帯が作れず、肯定側の枝を通せなかった**ため。
    **これは実データへの予測でもある**：**春と秋の VPD が重ならなければ、この検定は動かない。**
    """
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2008-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy()
    Rg = np.clip(150 + 120 * np.sin(2 * np.pi * (doy - 80) / 365)
                 + rng.normal(0, 25, len(idx)), 5, None)
    wet = (np.exp(-0.5 * ((doy - 90) / 35.) ** 2) + np.exp(-0.5 * ((doy - 230) / 35.) ** 2))
    th = np.clip(0.14 + 0.10 * wet + rng.normal(0, 0.035, len(idx)), 0.02, 0.6)
    Ts = 15 + 12 * np.sin(2 * np.pi * (doy - 120) / 365) + rng.normal(0, 3, len(idx))
    if kind == "collinear":
        VPD = 0.5 + 0.010 * Rg + rng.normal(0, 0.05, len(idx))      # **ほぼ Rg の関数**
    else:
        # **春（モンスーン前）に高く、秋（モンスーン後）に低い**——doy 150 前後で最大
        nz = 0.90 if kind.endswith("_overlap") else 0.30
        VPD = np.clip(1.6 + 1.4 * np.sin(2 * np.pi * (doy - 60) / 365)
                      - 1.1 * wet + rng.normal(0, nz, len(idx)), 0.05, None)
    if kind.startswith("season"):
        on = pd.Series(idx.month).isin(AUTUMN).to_numpy()
    elif kind.startswith("wet"):
        on = th >= np.median(th)
    else:                                    # vpd / collinear：**VPD が低い日だけ反転**
        on = VPD <= np.median(VPD)
    beta = np.where(on, 1.6, 0.0)
    gLE = Rg * (0.25 + beta * (th - th.mean())) + rng.normal(0, 8, len(idx))
    gH = Rg * (0.45 - beta * (th - th.mean())) + rng.normal(0, 8, len(idx))
    return pd.DataFrame({"th": th, "Rg": Rg, "Ts": Ts, "VPD": VPD,
                         "gLE": np.clip(gLE, 0, None), "gH": np.clip(gH, 0, None)},
                        index=idx)


def main():
    ap = argparse.ArgumentParser(description="旗96：VPD で春と秋を揃える")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--sites", nargs="+",
                    default=["US-Wkg", "US-Whs", "US-SRM", "MN-Hst", "MN-Nkh", "MN-Kbu"])
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()

    print("=== 旗96：春に Bowen 反転が起きないのは、飽差（VPD）が高いからか ===")
    print("  **旗95 の打ち切りは早かった**——**VPD は一度も層別に使っていなかった。**")
    print(f"  **門①**：帯の中で**秋（対照）が反転しなければ判定しない**（旗95 の欠陥31）。")
    print(f"  **門②**：セル内 **|r(VPD,Rg)| > {R_MAX:.2f} なら判定しない**（VPD が Rg の言い換え）。")

    if not a.real:
        print("\n  【合成データで道具を検証する】**実データに触れる前に**、")
        print("  **`collinear` で門②が発火するか**を必ず見る。")
        for kind, want in (
                ("vpd", "**現実的な季節振幅**＝**帯が作れず判定しない**はず"
                        "（**実データへの予測でもある**）"),
                ("collinear", "**門②が発火すべき**"),
                ("vpd_overlap", "**重なりを作為的に広げた検証用**：春も反転＝explained"),
                ("season_overlap", "**同上**：春は反転しない＝not_explained"),
                ("wet_overlap", "**同上**：**春だけ反転してはいけない**")):
            print(f"\n  ===== 合成 `{kind}` —— 期待：{want} =====")
            v, why = run_site(f"合成-{kind}", synth(kind))
            print(f"  【判定】{v or ('判定しない（' + str(why) + '）')}")
        print("\n  → **collinear→門②・vpd_overlap→explained・season_overlap→not_explained**")
        print("     **なら道具は使える。** `vpd`（現実的な振幅）が「帯なし」で終わるのは正しい挙動であり、")
        print("     **同時に実データへの予測**である——**春と秋の VPD が重ならなければ検定は動かない。**")
        return

    out, why = {}, {}
    for s in a.sites:
        print(f"\n  ━━ {s} ━━")
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max, extra=("Ts", "VPD"))
        except Exception as e:
            print(f"    読み込み失敗 {type(e).__name__}: {str(e)[:110]}")
            why[s] = "読み込み失敗"; continue
        v, w = run_site(s, d)
        out[s] = v
        if w:
            why[s] = w

    print("\n  === 集計（事前登録の判定規則に当てる）===")
    lab = {"explained": "帯の中では春も反転（VPD で説明された）",
           "not_explained": "帯の中でも春は反転しない（VPD でも説明されない）"}
    for s in a.sites:
        v = out.get(s)
        print(f"    {s:<9}{lab.get(v, '判定しない（' + why.get(s, '—') + '）')}")
    vals = [v for v in out.values() if v]
    n = len(vals)
    print("\n  === 結論 ===")
    if n < 2:
        print(f"  **判定しない**——門を通ったサイトが {n} で 2 未満。")
        print("  **どの門で落ちたか**は上の各行に書いてある。")
    elif sum(v == "explained" for v in vals) > n / 2:
        print("  **★VPD で説明された**——**飽差を揃えると、春でも Bowen 反転が起きる**。")
        print("  ＝A-3 を『**θ が高く、かつ飽差が低いとき Bowen 反転が起きる**』と書き換える。")
        print("  **ただし「VPD が原因」とは言わない**——**揃えると差が消える**までである。")
        print("  **VPD は Ta と RH の決定的な関数**なので、**独立な駆動とは書かない。**")
    else:
        print("  **▲VPD でも説明されない**——**測定量（θ・Rg・Ts・GCC・VPD）を使い切った**と確定する。")
        print("  ＝**残るのは新規観測**（`NEW_OBSERVATION_DESIGN.md`）**か別の枠組み**である。")
    print("\n  留保（事前登録どおり）：")
    print("   ・**独立クラスタは 2 つ**（Walnut Gulch・Santa Rita）＝3 サイト≠3 反復。")
    print("   ・**VPD は `es(Ta)(1−RH/100)`**＝旗32 が「独立変数として綺麗でない」と名指しした量。")
    print("     **層別に使っただけであり、独立な駆動とは書かない。**")
    print("   ・**帯に絞れば n は減り相関は減衰する**＝**判定は符号と CI のみ**。")


if __name__ == "__main__":
    main()
