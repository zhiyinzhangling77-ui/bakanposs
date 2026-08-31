"""旗93：**春と秋の差は、植生の活性（フェノロジー）で説明できるか**（事前登録 step93）。

**派生量には踏み込まない**——`NEE` はこの研究が**独立に測っている 8 つ**に数えている量であり、
**派生量は NEE を分割して作る GPP と GER の方**である（旗32）。**GPP も LAI も使わない。**

**原理的な危険が二つあり、設計で受けている**（事前登録 step93）：
  ① **乾燥地では「緑」と「湿り」が切り離せない**（緑の春＝雨の多かった春）
     → **旗90 と同じ重なり帯で θ を揃えたうえで緑度だけを変える**。
       **帯が作れなければ「水と分離できない」と結論して終える。**
  ② **NEE は θ→気孔→蒸散/光合成 の経路上のコライダーになりうる**
     → **日ごとの NEE では条件付けない**。**年ごとの春の平均 NEE で年をラベル**し、
       **検定は各春の日次データの中で行う**＝**条件付ける階層と検定する階層を分ける**。

**結論：この検定は実行しない**（旗93）。合成検証で、
**緑度と水を区別できる設定（4月のみ×年の上下1/4）は実データで日数が下限を割り**、
**日数の足りる設定（4–5月）は区別に失敗する**と分かった。
＝**区別できる設計には検出力が無く、検出力のある設計は区別できない。**
**下限を緩めて答えを出すことはしない。**

    python research/phenology_match_step93.py                  # 合成で検証（これだけ使う）
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
from soiltemp_match_step90 import band, SPRING, AUTUMN, PLO, PHI

# **追補で固定**：年でラベルするなら、**日付は狭く固定しなければならない**。
# 合成検証で分かった欠陥——**春の中の日付による緑度の変化が、年による変化より大きい**ため、
# 春（3–5月）全体を使うと**枯れた年の 5 月にも「緑の日」が入り**、
# `pheno`（緑度で反転）と `wet`（水で反転）が**どちらも「both」になって区別できなかった**。
# ＝**同じ暦の窓の中で、年だけを変えて比べる**。
WINDOW = (4, 4)          # 4月のみ（主検定）
WINDOW_NAME = "4月"


def greenness_by_year(d):
    """**年ごとの春の平均 NEE**＝緑度の指標（**負ほど活性が高い**）。

    **その年の春の日を全部使う**（θ高×Rg高 に絞らない）——**指標を安定させるため**。
    **日ごとの NEE では条件付けない**（コライダー対策・事前登録どおり）。
    """
    sp = d[[m in SPRING for m in d.index.month]]
    if "NEE" not in sp.columns:
        return None
    g = sp.groupby(sp.index.year)["NEE"].mean().dropna()
    return g if len(g) >= 2 * MIN_YEARS else None


def run_site(site, d, verbose=True):
    if verbose:
        print(f"\n  ━━ {site} ━━")
    need = {"th", "Rg", "gLE", "gH", "NEE"}
    if not need <= set(d.columns):
        print(f"    **変数が足りない**（在るのは {sorted(set(d.columns))}）＝判定しない")
        return None
    g = greenness_by_year(d)
    if g is None:
        print(f"    **緑度を作れない**（春の年数が {2*MIN_YEARS} 未満）＝判定しない")
        return None
    med = float(g.median())
    green_yrs = set(g.index[g <= med])      # NEE が小さい（負が大きい）＝**緑**
    brown_yrs = set(g.index[g > med])
    if verbose:
        print(f"    緑度＝春の平均 NEE（負ほど活性）：中央値 {med:+.2f}／"
              f"緑の年 {len(green_yrs)}・枯れた年 {len(brown_yrs)}"
              f"（範囲 {g.min():+.2f}〜{g.max():+.2f}）")

    lab, tmed, rmed = cell_of(d)
    hh = d[lab == "θ高×Rg高"]
    # **春全体ではなく、暦の窓を固定する**（追補）＝**日付を揃えて年だけを変える**
    sp = hh[[WINDOW[0] <= m <= WINDOW[1] for m in hh.index.month]]
    gr = sp[[y in green_yrs for y in sp.index.year]]
    br = sp[[y in brown_yrs for y in sp.index.year]]
    if verbose:
        print(f"    θ・Rg のしきい値（旗89 と同一）：θ={tmed:.3f}／Rg={rmed:.1f}"
              f"／θ高×Rg高 の**{WINDOW_NAME}**は {len(sp)} 日（緑 {len(gr)}・枯 {len(br)}）")
    if len(gr) == 0 or len(br) == 0:
        print("    片方が空＝判定しない"); return None

    lo, hi = band(gr["th"].to_numpy(), br["th"].to_numpy())
    if hi <= lo:
        print(f"    → **θ の帯幅が 0 以下**（[{lo:.3f}, {hi:.3f}]）"
              f"＝**緑度と水が完全に分離している**＝**水と切り離せない**")
        return None
    grb = gr[(gr["th"] >= lo) & (gr["th"] <= hi)]
    brb = br[(br["th"] >= lo) & (br["th"] <= hi)]
    if verbose:
        print(f"    **θ の重なり帯 [{lo:.3f}, {hi:.3f}]**（幅 {hi-lo:.3f}）／残存："
              f"緑 {len(grb)}日({len(grb)/max(len(gr),1):.0%})・"
              f"枯 {len(brb)}日({len(brb)/max(len(br),1):.0%})")
        print(f"      帯の中の θ 中央値：緑 {grb['th'].median():.3f}／枯 {brb['th'].median():.3f}"
              f"（**近いほど水が揃っている**）")

    out = {}
    for name, sub in ((f"緑の年の{WINDOW_NAME}", grb), (f"枯れた年の{WINDOW_NAME}", brb)):
        res = test_cell(sub)
        if res is None:
            print(f"      {name}：日数 {len(sub)}／年 "
                  f"{sub.index.year.nunique() if len(sub) else 0}＝**判定しない**（下限未満）")
            out[name] = None; continue
        rev = reversed_(res)
        out[name] = rev
        print(f"      {name}：日数 {res['n']:>5}／年 {res['yrs']:>3}  "
              f"{_fmt(res['le'])}  {_fmt(res['h'])}  "
              f"{'**Bowen反転**' if rev else '反転せず'}")
    # 参考：同じ帯の中の秋（**反転する側の対照**）
    au = hh[[m in AUTUMN for m in hh.index.month]]
    aub = au[(au["th"] >= lo) & (au["th"] <= hi)]
    ares = test_cell(aub)
    if ares is not None:
        print(f"      参考・秋（同じ帯）：日数 {ares['n']:>5}／年 {ares['yrs']:>3}  "
              f"{_fmt(ares['le'])}  {_fmt(ares['h'])}  "
              f"{'**Bowen反転**' if reversed_(ares) else '反転せず'}")
    else:
        print(f"      参考・秋（同じ帯）：日数 {len(aub)}＝判定しない")

    g_key, b_key = f"緑の年の{WINDOW_NAME}", f"枯れた年の{WINDOW_NAME}"
    if out.get(g_key) is None or out.get(b_key) is None:
        print("      → **このサイトは判定しない**（帯の中の日数・年数が下限未満）")
        return None
    if out[g_key] and not out[b_key]:
        print(f"      → **緑の年の{WINDOW_NAME}だけ反転**＝**フェノロジーで説明される**")
        return "pheno"
    if not out[g_key] and not out[b_key]:
        print(f"      → **どちらの年の{WINDOW_NAME}でも反転しない**＝**フェノロジーでも説明されない**")
        return "neither"
    if out[g_key] and out[b_key]:
        print(f"      → **どちらの年の{WINDOW_NAME}でも反転する**＝**旗89/90 と食い違う**"
              "＝**帯の作り方を疑うこと**（事前登録どおり）")
        return "both"
    print(f"      → **枯れた年の{WINDOW_NAME}だけ反転**＝**想定外**")
    return "brown_only"


def synth(kind, years=20, seed=0):
    """**水と緑を取り違えていないか**を含めて、道具を三つの場合で試す。

      ・``pheno``  —— **緑度が高い日だけ**反転。**緑度は θ と独立な年変動を持つ**
        ＝**分離できる状況**。→ **緑の春だけ反転**が出るべき。
      ・``season`` —— **秋だけ**反転。→ **どちらの春でも反転しない**が出るべき。
      ・``wet``    —— **θ が高い日だけ**反転し、**緑度は θ の関数**。
        → **帯で θ を揃えれば差が消える**べき。
        **ここで「緑の春だけ反転」が出たら、道具は水と緑を分離できていない＝実データに進まない。**
    """
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2004-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy(); yr = idx.year.to_numpy()
    Rg = np.clip(150 + 120 * np.sin(2 * np.pi * (doy - 80) / 365)
                 + rng.normal(0, 25, len(idx)), 5, None)
    wet_shape = (np.exp(-0.5 * ((doy - 90) / 35.0) ** 2)
                 + np.exp(-0.5 * ((doy - 230) / 35.0) ** 2))
    th = np.clip(0.14 + 0.10 * wet_shape + rng.normal(0, 0.035, len(idx)), 0.02, 0.6)
    mon = pd.Series(idx.month)
    is_sp, is_au = mon.isin(SPRING).to_numpy(), mon.isin(AUTUMN).to_numpy()

    # 緑度：季節の形 × 年ごとの係数
    canopy = np.exp(-0.5 * ((doy - 175) / 70.0) ** 2)        # 生育期の形
    uy = {y: rng.uniform(0.3, 1.7) for y in np.unique(yr)}    # **θ と独立**な年変動
    yfac = np.array([uy[y] for y in yr])
    if kind == "wet":
        # 緑度が **θ の関数**＝水と緑が絡む（分離できない状況）
        green = canopy * (th / th.mean())
    else:
        green = canopy * yfac
    # NEE：緑ほど負（呼吸の底上げも入れる）
    NEE = 2.0 - 9.0 * green + rng.normal(0, 0.8, len(idx))

    if kind == "season":
        on = is_au
    elif kind == "wet":
        on = th >= np.median(th)
    else:                                    # pheno
        on = green >= np.median(green)
    beta = np.where(on, 1.6, 0.0)
    gLE = Rg * (0.25 + beta * (th - th.mean())) + rng.normal(0, 8, len(idx))
    gH = Rg * (0.45 - beta * (th - th.mean())) + rng.normal(0, 8, len(idx))
    return pd.DataFrame({"th": th, "Rg": Rg, "NEE": NEE,
                         "gLE": np.clip(gLE, 0, None), "gH": np.clip(gH, 0, None)},
                        index=idx)


def main():
    p = argparse.ArgumentParser(description="旗93：緑度で春を割る")
    p.add_argument("--real", action="store_true")
    p.add_argument("--sites", nargs="+", default=["US-Wkg", "US-Whs", "US-SRM"])
    p.add_argument("--qc-max", type=int, default=None)
    a = p.parse_args()

    print("=== 旗93：春と秋の差はフェノロジーで説明できるか ===")
    print("  **派生量には踏み込まない**——緑度は **NEE**（独立測定の 8 つの一つ）で代理する。")
    print("  **GPP も LAI も使わない。**")
    print("  危険①**緑と湿りは乾燥地で切り離せない** → **θ の重なり帯で揃える**。")
    print("  危険②**NEE はコライダーになりうる** → **年単位でラベルし、検定は日次**。")

    if a.real:
        print("\n  " + "!" * 68)
        print("  **この検定は実データに使ってはいけない**（旗93 の合成検証で確定）。")
        print("  **緑度で反転する系列（pheno）と、水で反転する系列（wet）を区別できる設定は**")
        print("  **『4月のみ × 年の上下 1/4』だけ**だが、**実データではその窓の日数が下限 60 を割る**")
        print("  （推定：US-Wkg 51・US-Whs 35・US-SRM 37 日）。")
        print("  **窓を 4–5 月に広げると日数は足りるが、乱数の種を変えると 4 回中 2 回で区別に失敗する。**")
        print("  ＝**区別できる設計には検出力が無く、検出力のある設計は区別できない。**")
        print("  **下限を緩めて答えを出すことはしない**（旗58 の教訓）。")
        print("  出力を見ても**結論には使わないこと**。詳細は FLAGS_LOG.md の旗93。")
        print("  " + "!" * 68)
    if not a.real:
        print("\n  【合成データで道具を検証する】**実データに触れる前に**、")
        print("  **`wet`（θ で反転し緑度は θ の関数）で「緑の春だけ反転」が出ないか**を見る。")
        print("  **出たら水と緑を分離できていない＝実データに進まない。**")
        for kind, want in (
                ("pheno", "**緑の春だけ反転**が出るべき（緑度は θ と独立な年変動）"),
                ("season", "**どちらの春でも反転しない**が出るべき"),
                ("wet", "**緑の春だけ反転してはいけない**（水と緑の取り違えの試験）")):
            print(f"\n  ===== 合成 `{kind}` —— 期待：{want} =====")
            print(f"  【判定】{run_site(f'合成-{kind}', synth(kind))}")
        print("\n  → **`pheno` で pheno、`season` で neither、`wet` で pheno 以外**なら道具は使える。")
        return

    out = {}
    for s in a.sites:
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max, extra=("Ts", "NEE"))
        except Exception as e:
            print(f"\n  ━━ {s} ━━\n    読み込み失敗 {type(e).__name__}: {str(e)[:120]}")
            continue
        out[s] = run_site(s, d)

    print("\n  === 集計（事前登録の判定規則に当てる）===")
    lab = {"pheno": f"緑の年の{WINDOW_NAME}だけ反転（フェノロジーで説明される）",
           "neither": f"どちらでも反転しない（説明されない）",
           "both": "どちらでも反転する（旗89/90 と食い違う）",
           "brown_only": "枯れた年だけ反転（想定外）",
           None: "判定できない"}
    for s, v in out.items():
        print(f"    {s:<9}{lab.get(v, str(v))}")
    vals = [v for v in out.values() if v]
    n = len(vals)
    print("\n  === 結論 ===")
    if n < 2:
        print(f"  **判定しない**——判定できたサイトが {n} で 2 未満。")
        print("  **帯が作れなかったのか日数が足りなかったのか**は、上の各サイトの行に書いてある。")
    elif sum(v == "pheno" for v in vals) > n / 2:
        print("  **★フェノロジーで説明された**——**θ を揃えても、緑の春だけ反転する**。")
        print("  ＝A-3 を『**θ が高く、かつ林冠が活性のとき Bowen 反転が起きる**』と書き換える。")
        print("  **ただし「林冠活性が原因」とは言わない**——**緑度を変えると差が消える**までである。")
    elif sum(v == "neither" for v in vals) > n / 2:
        print("  **▲フェノロジーでも説明されない**——**どちらの春でも反転しない**。")
        print("  ＝旗88–91 と合わせ、**観測データの層別による説明は打ち切る**と確定する。")
    elif sum(v == "both" for v in vals) > n / 2:
        print("  **○旗89/90 と食い違う**——**θ を揃えると春でも反転する**ことになる。")
        print("  ＝**本検定ではなく、帯の作り方（旗90 の 10/90 パーセンタイル）を疑う**")
        print("     ——**事前登録でそう決めてある。後から都合よく解釈しない。**")
    else:
        print("  **判定が割れた**——各サイトの行をそのまま記録し、**まとめない**。")
    print("\n  留保（事前登録どおり）：")
    print("   ・**独立クラスタは 2 つ**（Walnut Gulch・Santa Rita）＝3 サイト≠3 反復。")
    print("   ・**`NEE_VUT_REF` は穴埋め済み**（旗46）＝**緑度の指標に穴埋めが入る**。")
    print("   ・**NEE は呼吸も含む**＝緑度の代理として完全ではない（**検出力を下げる向き**）。")
    print("   ・**年を 2 群に割ると群内の年数が半減**＝CI は広がる。**判定は符号と CI のみ。**")


if __name__ == "__main__":
    main()
