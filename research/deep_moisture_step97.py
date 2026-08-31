"""旗97：**春に Bowen 反転が起きないのは、水が深くまで届いていないからか**（事前登録 step97）。

旗88–96 は **θ・Rg・Ts・GCC・VPD の“水準”**で説明を試みて、すべて外れた。
**だが、試した θ は `SWC_F_MDS_1`（最浅の慣例）だけ**である——
**FLUXNET 版には最初から複数層が入っている**（旗79：US-Wkg 3 層・US-Whs 6 層・US-SRM 8 層）。

**仮説**：**春の前線性降雨は浅く濡らし、夏のモンスーンは深くまで濡らす。**
**「同じ θ」で揃えたつもりが、表層だけ揃えていた**のではないか。

**事前登録 step97 で固定済み**：
  ・**セルは旗89 のまま**（`SWC_F_MDS_1` による θ高×Rg高）＝**比較可能性のため作り直さない**
  ・**門①**：帯の中で**秋が反転しなければ判定しない**（旗95/96 と同じ）
  ・**門②**：`|r(θ_1, θ_deep)| > 0.90` なら判定しない（**深層が浅層の言い換え**）
  ・**重なりが無い場合**は「**春は深層まで濡れない**」という**記述**であって**検定ではない**
  ・**副検定**：春の中だけで θ_deep を高低に分け、**θ_1 の帯で交絡を潰す**（旗96 の教訓）

    python research/deep_moisture_step97.py                # 合成で検証（既定）
    python research/deep_moisture_step97.py --real
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaporation_regime_step36 import daily_energy, _fmt
from stratified_bowen_step89 import cell_of, test_cell, reversed_, MIN_DAYS, MIN_YEARS
from soiltemp_match_step90 import band, SPRING, AUTUMN
from vpd_match_step96 import spearman

R_SAME = 0.90                 # 門②の閾値（事前登録で固定）
MIN_LAYER_DAYS = 1000         # θ_deep に採る層の下限（事前登録で固定）
LAYER_RE = re.compile(r"^SWC_F_MDS_(\d+)$")


def load_layers(site, qc_max=None):
    """**HH ファイルから `SWC_F_MDS_<n>` を全部読み、日平均にする。**

    `daily_energy` の `extra` は **R&K 表記の 11 変数しか通らない**（`load_raw_all` が
    それしか返さない）ので、**層は生ファイルから直接読む**。
    """
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import find_corevars_files
    cfg = AnalysisConfig(qc_max=qc_max) if qc_max is not None else AnalysisConfig()
    files = find_corevars_files(get_site(site))
    parts = []
    for f in files:
        head = pd.read_csv(f, nrows=2)
        cols = [c for c in head.columns if LAYER_RE.match(str(c))]
        if not cols:
            continue
        want = set(["TIMESTAMP_START"] + cols)
        df = pd.read_csv(f, usecols=lambda c: c in want)
        ts = pd.to_datetime(pd.to_numeric(df["TIMESTAMP_START"]).astype("int64").astype(str),
                            format="%Y%m%d%H%M")
        df = df.drop(columns=["TIMESTAMP_START"]).replace(cfg.na_sentinel, np.nan)
        df.index = ts
        parts.append(df)
    if not parts:
        return None
    raw = pd.concat(parts)
    raw = raw[~raw.index.duplicated(keep="first")].sort_index()
    return raw.groupby(raw.index.normalize()).mean()


def pick_deep(layers, th1, verbose=True):
    """**層の棚卸し**を出し、**θ_deep**（有効日数 ≥1000・年数 ≥3 の最大添字）を選ぶ。"""
    idx = sorted((int(LAYER_RE.match(c).group(1)), c) for c in layers.columns)
    if verbose:
        print(f"    {'層':<16}{'有効日':>8}{'年':>5}{'中央値':>9}{'r(θ_1と)':>11}")
    ok = []
    for n, c in idx:
        s = layers[c].dropna()
        r = spearman(th1.reindex(s.index), s) if n != 1 else 1.0
        if verbose:
            print(f"    {c:<16}{len(s):>8}{s.index.year.nunique():>5}"
                  f"{s.median():>9.2f}{r:>11.2f}")
        if len(s) >= MIN_LAYER_DAYS and s.index.year.nunique() >= MIN_YEARS and n != 1:
            ok.append((n, c, r))
    if not ok:
        return None, None
    n, c, r = ok[-1]                      # **添字が最大＝最も深い（慣例）**
    return c, r


def test_pair(a, b, name_a, name_b, key, label):
    """`key` の重なり帯で `a`/`b` を絞り、それぞれ検定する。**帯と結果を返す。**"""
    lo, hi = band(a[key].to_numpy(), b[key].to_numpy())
    print(f"    {label}：{name_a} 中央値 {a[key].median():.2f}／"
          f"{name_b} 中央値 {b[key].median():.2f}")
    if hi <= lo:
        print(f"      → **{key} の帯が作れない**（[{lo:.2f}, {hi:.2f}]）")
        return None, None, (lo, hi)
    ab = a[(a[key] >= lo) & (a[key] <= hi)]
    bb = b[(b[key] >= lo) & (b[key] <= hi)]
    print(f"      **帯 [{lo:.2f}, {hi:.2f}]**（幅 {hi-lo:.2f}）"
          f"／帯の中：{name_a} {len(ab)} 日・{name_b} {len(bb)} 日")
    out = {}
    for nm, sub in ((name_a, ab), (name_b, bb)):
        r = test_cell(sub)
        if r is None:
            print(f"      {nm}：日数 {len(sub)}／年 "
                  f"{sub.index.year.nunique() if len(sub) else 0}＝**下限未満**")
            out[nm] = None; continue
        rev = reversed_(r); out[nm] = rev
        print(f"      {nm}：日数 {r['n']:>4}／年 {r['yrs']:>3}  "
              f"{_fmt(r['le'])}  {_fmt(r['h'])}  "
              f"{'**Bowen反転**' if rev else '反転せず'}")
    return out.get(name_a), out.get(name_b), (lo, hi)


def run_site(tag, d, layers):
    """1 サイト。**棚卸し → 門② → 主検定（門①）→ 副検定**。"""
    if layers is None or layers.empty:
        print("    **層が読めない**＝判定しない"); return None, "層なし", None
    lab, tmed, rmed = cell_of(d)
    hh = d[lab == "θ高×Rg高"].copy()
    print(f"    θ={tmed:.3f}／Rg={rmed:.1f}（旗89 と同一）／θ高×Rg高 {len(hh)} 日")
    print("    ── 層の棚卸し ──")
    deep, r_same = pick_deep(layers, d["th"])
    if deep is None:
        print("    → **下限を満たす深い層が無い**＝判定しない"); return None, "深層なし", None
    print(f"    **θ_deep = {deep}**（θ_1 との r = {r_same:+.2f}）")
    if abs(r_same) > R_SAME:
        print(f"    → **門②：|r| > {R_SAME:.2f}＝深層が浅層の言い換え＝判定しない**")
        return None, "門②（深層≒浅層）", None

    hh["deep"] = layers[deep].reindex(hh.index)
    hh = hh.dropna(subset=["deep"])
    sp = hh[[m in SPRING for m in hh.index.month]]
    au = hh[[m in AUTUMN for m in hh.index.month]]
    print(f"    深層が在る日：春 {len(sp)}・秋 {len(au)}")
    if len(sp) < 10 or len(au) < 10:
        print("    → 片方が少なすぎる＝判定しない"); return None, "日数不足", None

    print("    ── 主検定：**深層 θ を揃えて春と秋を比べる** ──")
    v_sp, v_au, (lo, hi) = test_pair(sp, au, "春", "秋", "deep", "深層 θ")
    if hi <= lo:
        print(f"      → **春と秋の深層 θ が重ならない**＝**春は深層まで濡れない**という"
              f"**記述**である。**仮説の検定ではない**（事前登録どおり）。")
        print(f"      春 10/90 = [{np.percentile(sp['deep'],10):.2f}, "
              f"{np.percentile(sp['deep'],90):.2f}]／"
              f"秋 = [{np.percentile(au['deep'],10):.2f}, "
              f"{np.percentile(au['deep'],90):.2f}]")
        main = None; why = "帯なし（春は深層まで濡れない）"
    elif v_sp is None or v_au is None:
        main = None; why = "下限未満"
    elif not v_au:
        print("      → **門①：帯の中で秋が反転しない＝対照が働いていない＝判定しない**")
        main = None; why = "門①（対照が落ちた）"
    else:
        main = "explained" if v_sp else "not_explained"
        print(f"      → {'**帯の中では春も反転**＝**深層の水で説明された**' if v_sp else '**帯の中でも春は反転しない**＝**深層の水でも説明されない**'}")
        why = None

    # ── 副検定：春の中だけで深層 θ を高低に分ける ──
    print("    ── 副検定：**春の中だけで深層 θ を高低に分ける**（θ_1 の帯で交絡を潰す）──")
    dm = float(sp["deep"].median())
    wet_g, dry_g = sp[sp["deep"] >= dm], sp[sp["deep"] < dm]
    sub = None
    if len(wet_g) < 10 or len(dry_g) < 10:
        print("      片方が少なすぎる＝判定しない")
    else:
        a, b, _ = test_pair(wet_g, dry_g, "深部湿潤の春", "深部乾燥の春", "th", "表層 θ")
        if a is None or b is None:
            print("      → **判定しない**")
        elif a and not b:
            print("      → **深部湿潤の春だけ反転**＝**深い水が効いている**"); sub = "wet_only"
        elif a and b:
            print("      → **どちらも反転**＝**深さは効いていない**"); sub = "both"
        elif not a and not b:
            print("      → **どちらも反転しない**＝**春では深層を濡らしても反転しない**")
            sub = "neither"
        else:
            print("      → **深部乾燥の春だけ反転**＝**想定と逆**（解釈しない）"); sub = "dry_only"
    if main == "explained" and sub not in ("wet_only", None):
        print("      ※**副検定と食い違う**——**explained は深さの手柄とは限らない。両方を書く。**")
    return main, why, sub


def synth(kind, years=16, seed=0):
    """四つの場合。**`identical` は門②の試験**、**`shallow` は副検定の試験**である。"""
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2008-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy()
    Rg = np.clip(150 + 120 * np.sin(2 * np.pi * (doy - 80) / 365)
                 + rng.normal(0, 25, len(idx)), 5, None)
    spring_rain = np.exp(-0.5 * ((doy - 90) / 35.) ** 2)      # 前線性（浅い）
    monsoon = np.exp(-0.5 * ((doy - 230) / 35.) ** 2)         # モンスーン（深い）
    th = np.clip(0.14 + 0.10 * (spring_rain + monsoon) + rng.normal(0, 0.035, len(idx)), .02, .6)
    if kind == "identical":
        deep = th + rng.normal(0, 0.002, len(idx))            # **ほぼ同一**
    else:
        # **深層はモンスーンでしか濡れない**＋独立な雑音
        deep = np.clip(0.10 + 0.12 * monsoon + 0.02 * spring_rain
                       + rng.normal(0, 0.030, len(idx)), .02, .6)
    Ts = 15 + 12 * np.sin(2 * np.pi * (doy - 120) / 365) + rng.normal(0, 3, len(idx))
    if kind == "season":
        on = pd.Series(idx.month).isin(AUTUMN).to_numpy()
    elif kind == "shallow":
        on = th >= np.median(th)
    else:                                     # deep / identical
        on = deep >= np.median(deep)
    beta = np.where(on, 1.6, 0.0)
    gLE = Rg * (0.25 + beta * (th - th.mean())) + rng.normal(0, 8, len(idx))
    gH = Rg * (0.45 - beta * (th - th.mean())) + rng.normal(0, 8, len(idx))
    tow = pd.DataFrame({"th": th, "Rg": Rg, "Ts": Ts,
                        "gLE": np.clip(gLE, 0, None), "gH": np.clip(gH, 0, None)}, index=idx)
    lay = pd.DataFrame({"SWC_F_MDS_1": th, "SWC_F_MDS_3": deep}, index=idx)
    return tow, lay


def main():
    ap = argparse.ArgumentParser(description="旗97：深い層の θ で層別する")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--sites", nargs="+", default=["US-Wkg", "US-Whs", "US-SRM"])
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()

    print("=== 旗97：春に Bowen 反転が起きないのは、水が深くまで届いていないからか ===")
    print("  **試してきた θ は `SWC_F_MDS_1` だけ**——**複数層は最初から手元にあった。**")
    print(f"  **門①**：帯の中で秋が反転しなければ判定しない／"
          f"**門②**：|r(θ_1, θ_deep)| > {R_SAME:.2f} なら判定しない。")
    print("  **重なりが無い場合は「春は深層まで濡れない」という記述であって検定ではない。**")

    if not a.real:
        print("\n  【合成データで道具を検証する】**実データに触れる前に**、")
        print("  **`identical` で門②が発火するか**と、**`shallow` で副検定が誤らないか**を見る。")
        for kind, want in (
                ("deep", "主検定 ★・副検定「深部湿潤群だけ反転」が出るべき"),
                ("shallow", "**副検定で「両群とも反転」**が出るべき（浅層を拾ってはいけない）"),
                ("season", "**▲**が出るべき"),
                ("identical", "**門②が発火すべき**")):
            print(f"\n  ===== 合成 `{kind}` —— 期待：{want} =====")
            tow, lay = synth(kind)
            m, why, sub = run_site(f"合成-{kind}", tow, lay)
            print(f"  【判定】主={m or ('判定しない（'+str(why)+'）')}／副={sub}")
        print("\n  → **deep→explained+wet_only・shallow→副が both・season→▲・identical→門②**")
        print("     **なら道具は使える。**")
        return

    out, why_, sub_ = {}, {}, {}
    for s in a.sites:
        print(f"\n  ━━ {s} ━━")
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max, extra=("Ts",))
            lay = load_layers(s, a.qc_max)
        except Exception as e:
            print(f"    読み込み失敗 {type(e).__name__}: {str(e)[:110]}")
            why_[s] = "読み込み失敗"; continue
        m, w, sub = run_site(s, d, lay)
        out[s] = m; sub_[s] = sub
        if w:
            why_[s] = w

    print("\n  === 集計（事前登録の判定規則に当てる）===")
    lab = {"explained": "帯の中では春も反転（深層の水で説明された）",
           "not_explained": "帯の中でも春は反転しない（深層の水でも説明されない）"}
    for s in a.sites:
        print(f"    {s:<9}主={lab.get(out.get(s), '判定しない（' + why_.get(s, '—') + '）')}"
              f"／副={sub_.get(s)}")
    vals = [v for v in out.values() if v]
    n = len(vals)
    print("\n  === 結論（**追補どおり、判定の主軸は副検定**）===")
    print("  **主検定は、仮説が正しいほど成立しない**——合成で確認済み。")
    print("  **「春と秋で深層 θ が違う」のが仮説そのもの**なので、**重なり帯が作れない**。")
    if n:
        print(f"  参考・主検定の内訳：{ {s: out[s] for s in out if out[s]} }")
    sv = [v for v in sub_.values() if v]
    m = len(sv)
    if m < 2:
        print(f"  **判定しない**——**副検定で判定できたサイトが {m} で 2 未満**。")
        print("  **どの門・どの下限で落ちたか**は上の各行に書いてある。")
    elif sum(v == "wet_only" for v in sv) > m / 2:
        print("  **★深い水が効いている**——**表層を揃えても、深層が濡れた春だけ反転する**。")
        print("  ＝A-3 の季節依存を『**深層まで水が届いた春では反転する**』と書く。")
        print("  **ただし副検定には対照が無い**（春の中だけで秋を使わない）＝**そこは弱い。**")
    elif sum(v == "both" for v in sv) > m / 2:
        print("  **▲深さは効いていない**——**深層が乾いた春でも反転する**（表層で足りる）。")
    elif sum(v == "neither" for v in sv) > m / 2:
        print("  **▲春では深層を濡らしても反転しない**＝**深さでも説明されない**。")
    else:
        print("  **判定が割れた**——各サイトの行をそのまま記録し、**まとめない**。")
    print("\n  留保（事前登録どおり）：")
    print("   ・**深度は列名から分からない**（旗33/80）＝**「_1 より深い層」までしか言えない**。")
    print("   ・**層は本来強く相関する**——**門②を通っても独立とは言えない**。")
    print("   ・**独立クラスタは 2 つ**（Walnut Gulch・Santa Rita）＝3 サイト≠3 反復。")
    print("   ・**深い層ほど欠測が多い可能性**＝**有効日数を棚卸しに出してある**。
   ・**副検定には対照が無い**（春の中だけで秋を使わない）。**代わりに合成 `shallow` で**
     **「浅層を拾わない」ことを確かめてある**が、**主検定の門①ほど強くはない**。
   ・**重なり帯という道具は、“群間で違う量”を揃えるのに使うと自己矛盾する**——
     **旗90(Ts)・旗95(GCC)・旗96(VPD) の「判定しない」の一部も、この構造による可能性がある。**")


if __name__ == "__main__":
    main()
