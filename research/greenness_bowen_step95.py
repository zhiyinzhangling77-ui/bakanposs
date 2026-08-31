"""旗95：**春に Bowen 反転が起きないのは、林冠が緑でないからか**（事前登録 step95）。

旗88–91 で**季節依存は θ・Rg・Ts では説明できない**と確定し、旗93 は**検出力の壁**で取りやめ、
旗94 で **PhenoCam の日次 GCC** を得た。**日ごとに緑度で層別できる**ようになったので検定する。

**事前登録 step95 で固定済み**（旗94 の下調べが突きつけた三つを含む）：
  ・**主 ROI はタワーの生態系名と一致するもの**（Wkg→GR・Whs→SH・SRM→SH）。**他方は感度確認。**
    **食い違ったら両方報告し、都合の良い方を選ばない。**
  ・**主指標は `gcc_90`（平滑なし・90 パーセンタイル）**。
    **理由は権威ではなく旗94 の発見**——**GCC は濡れと照度に応答する**ので、
    **高パーセンタイルは暗い/濡れたコマの影響を受けにくい**。
  ・**US-Whs はカメラ交換で ROI が 4 期**——**期ごとに中央値で割ってから束ねる**
    （**GCC の絶対値を期をまたいで比べない**）。
  ・場所・帯・統計量は**旗89/90/94 と同一**。**作り直さない。**

    python research/greenness_bowen_step95.py                    # 合成で検証（既定）
    python research/greenness_bowen_step95.py --real \\
        --phenocam-dir /mnt/hdd/PhenoCam
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
from phenocam_probe_step94 import read_summary

# ── 事前登録で固定（旗95）──
COL = "gcc_90"                       # 主指標：平滑なし 90 パーセンタイル
SENS_COLS = ("smooth_gcc_90", "gcc_mean")
PRIMARY = {"US-Wkg": ("kendall", "GR_1000"),      # Kendall Grasslands → 草
           "US-Whs": ("luckyhills", "SH_"),        # Lucky Hills Shrub → 低木（全期）
           "US-SRM": ("srm", "SH_1000")}           # Santa Rita Mesquite → 低木
SENSITIVITY = {"US-Wkg": ("kendall", "SH_1000"),
               "US-SRM": ("srm", "GR_1000")}


def gcc_files(root, site_key, roi_key):
    """`{site}_{roi}…_1day.csv` を集める。**ndvi と simplified は除く**（別製品）。"""
    out = []
    for p in sorted(Path(root).rglob("*_1day.csv")):
        n = p.name.lower()
        if site_key not in n or roi_key.lower() not in n:
            continue
        if "ndvi" in n or "simplified" in n or "transition" in n:
            continue
        out.append(p)
    return out


def load_gcc(paths, col):
    """**期ごとに**（＝ファイルごとに）読み、`(ラベル, 系列)` の一覧で返す。

    **束ねる前に期の中で中央値を取るため、ここでは束ねない**（事前登録どおり）。
    """
    out = []
    for p in paths:
        g, err = read_summary(p)
        if g is None:
            print(f"    ※{p.name}：読めない（{err}）"); continue
        if col not in g.columns:
            print(f"    ※{p.name}：**列 {col} が無い**（在るのは {list(g.columns)}）"); continue
        s = g[col].dropna()
        if len(s):
            out.append((p.name, s))
    return out


def split_pool(sp, eras):
    """**期ごとに GCC 中央値で緑/枯を割ってから束ねる**（事前登録どおり）。"""
    parts = []
    for name, s in eras:
        j = sp.join(s.rename("gcc"), how="inner").dropna(subset=["gcc"])
        if len(j) < 10:
            continue
        j = j.copy()
        j["green"] = j["gcc"] >= j["gcc"].median()   # **期の中で**割る
        j["era"] = name
        parts.append(j)
    if not parts:
        return None
    pool = pd.concat(parts)
    return pool[~pool.index.duplicated(keep="first")].sort_index()


def run_one(tag, sp, au, eras, verbose=True):
    """1 サイト×1 ROI を検定する。**帯は旗90/94 と同一**。"""
    pool = split_pool(sp, eras)
    if pool is None:
        print(f"    {tag}：**結合できる日が無い**"); return None
    g_, b_ = pool[pool["green"]], pool[~pool["green"]]
    if len(g_) < 10 or len(b_) < 10:
        print(f"    {tag}：緑 {len(g_)}・枯 {len(b_)}＝**片方が少なすぎる**"); return None
    lo, hi = band(g_["th"].to_numpy(), b_["th"].to_numpy())
    print(f"    {tag}：期 {len(eras)}／束ねた日 {len(pool)}（緑 {len(g_)}・枯 {len(b_)}）"
          f"／年 {pool.index.year.nunique()}")
    print(f"      GCC 中央値：緑 {g_['gcc'].median():.4f}／枯 {b_['gcc'].median():.4f}")
    print(f"      θ 中央値：緑 {g_['th'].median():.3f}／枯 {b_['th'].median():.3f}"
          f"  ← **緑の日が{'湿っている' if g_['th'].median() > b_['th'].median() else '乾いている'}**"
          f"（危険③がどちらに効いているか）")
    if hi <= lo:
        print(f"      → **θ の帯が作れない**（[{lo:.3f}, {hi:.3f}]）＝**緑と水が分離できない**")
        return None
    gb = g_[(g_["th"] >= lo) & (g_["th"] <= hi)]
    bb = b_[(b_["th"] >= lo) & (b_["th"] <= hi)]
    print(f"      **θ の帯 [{lo:.3f}, {hi:.3f}]**（幅 {hi-lo:.3f}）"
          f"／帯の中：緑 {len(gb)}・枯 {len(bb)}")
    res = {}
    for nm, sub in (("緑の日", gb), ("枯れた日", bb)):
        r = test_cell(sub)
        if r is None:
            print(f"      {nm}：日数 {len(sub)}／年 "
                  f"{sub.index.year.nunique() if len(sub) else 0}＝**判定しない**（下限未満）")
            res[nm] = None; continue
        rev = reversed_(r); res[nm] = rev
        print(f"      {nm}：日数 {r['n']:>4}／年 {r['yrs']:>3}  "
              f"{_fmt(r['le'])}  {_fmt(r['h'])}  "
              f"{'**Bowen反転**' if rev else '反転せず'}")
    # 参考：同じ帯の秋
    aub = au[(au["th"] >= lo) & (au["th"] <= hi)]
    ar = test_cell(aub)
    if ar is not None:
        print(f"      参考・秋（同じ帯）：日数 {ar['n']:>4}／年 {ar['yrs']:>3}  "
              f"{_fmt(ar['le'])}  {_fmt(ar['h'])}  "
              f"{'**Bowen反転**' if reversed_(ar) else '反転せず'}")
    if res.get("緑の日") is None or res.get("枯れた日") is None:
        print("      → **判定しない**"); return None
    if res["緑の日"] and not res["枯れた日"]:
        print("      → **緑の日だけ反転**＝**フェノロジーで説明される**"); return "pheno"
    if not res["緑の日"] and not res["枯れた日"]:
        print("      → **どちらでも反転しない**＝**説明されない**"); return "neither"
    if res["緑の日"] and res["枯れた日"]:
        print("      → **どちらでも反転する**＝**旗89/90 と食い違う＝帯の作り方を疑う**")
        return "both"
    print("      → **枯れた日だけ反転**＝**想定外**"); return "brown_only"


def prep_tower(d):
    lab, tmed, rmed = cell_of(d)
    hh = d[lab == "θ高×Rg高"]
    return (hh[[m in SPRING for m in hh.index.month]],
            hh[[m in AUTUMN for m in hh.index.month]], tmed, rmed)


def synth(kind, years=14, seed=0):
    """**日次の緑度で層別できるか**を三つの場合で試す（旗93 と違い年ラベルではない）。"""
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2010-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy(); yr = idx.year.to_numpy()
    Rg = np.clip(150 + 120 * np.sin(2 * np.pi * (doy - 80) / 365)
                 + rng.normal(0, 25, len(idx)), 5, None)
    wet = (np.exp(-0.5 * ((doy - 90) / 35.) ** 2) + np.exp(-0.5 * ((doy - 230) / 35.) ** 2))
    th = np.clip(0.14 + 0.10 * wet + rng.normal(0, 0.035, len(idx)), 0.02, 0.6)
    canopy = np.exp(-0.5 * ((doy - 200) / 55.) ** 2)
    uy = {y: rng.uniform(0.4, 1.6) for y in np.unique(yr)}     # **θ と独立**な年変動
    if kind == "wet":
        green = 0.33 + 0.12 * (th - th.mean()) + rng.normal(0, 0.001, len(idx))
        on = th >= np.median(th)
    else:
        green = (0.33 + 0.05 * canopy * np.array([uy[y] for y in yr])
                 + rng.normal(0, 0.004, len(idx)))
        on = (pd.Series(idx.month).isin(AUTUMN).to_numpy() if kind == "season"
              else green >= np.median(green))
    beta = np.where(on, 1.6, 0.0)
    gLE = Rg * (0.25 + beta * (th - th.mean())) + rng.normal(0, 8, len(idx))
    gH = Rg * (0.45 - beta * (th - th.mean())) + rng.normal(0, 8, len(idx))
    tow = pd.DataFrame({"th": th, "Rg": Rg,
                        "gLE": np.clip(gLE, 0, None), "gH": np.clip(gH, 0, None)}, index=idx)
    return tow, pd.Series(green, index=idx)


def main():
    ap = argparse.ArgumentParser(description="旗95：緑の日と枯れた日で春を割る")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--phenocam-dir", default="/mnt/hdd/PhenoCam")
    ap.add_argument("--qc-max", type=int, default=None)
    a = ap.parse_args()

    print("=== 旗95：春に Bowen 反転が起きないのは、林冠が緑でないからか ===")
    print(f"  **主指標 `{COL}`（平滑なし 90 パーセンタイル）**——"
          "**GCC は濡れと照度に応答する**ため高パーセンタイルを採る（旗94 の発見）。")
    print("  **主 ROI は生態系名と一致するもの**。**食い違ったら両方報告し、選ばない。**")
    print("  **US-Whs は期ごとに中央値で割ってから束ねる。**")

    if not a.real:
        print("\n  【合成データで道具を検証する】**実データに触れる前に**、")
        print("  **`wet`（θ で反転し緑度は θ の関数）で「緑の日だけ反転」が出ないか**を見る。")
        for kind, want in (("pheno", "**緑の日だけ反転**が出るべき"),
                           ("season", "**どちらでも反転しない**が出るべき"),
                           ("wet", "**緑の日だけ反転してはいけない**")):
            tow, g = synth(kind)
            sp, au, _, _ = prep_tower(tow)
            print(f"\n  ===== 合成 `{kind}` —— 期待：{want} =====")
            print(f"  【判定】{run_one(f'合成-{kind}', sp, au, [('synth', g)])}")
        print("\n  → **`pheno`→pheno・`season`→neither・`wet`→pheno 以外**なら道具は使える。")
        return

    out = {}
    for site in ("US-Wkg", "US-Whs", "US-SRM"):
        print(f"\n  ━━ {site} ━━")
        try:
            d, _ = daily_energy(site, list(range(1, 13)), a.qc_max, extra=("Ts",))
        except Exception as e:
            print(f"    タワーを読めない {type(e).__name__}: {str(e)[:90]}"); continue
        sp, au, tmed, rmed = prep_tower(d)
        print(f"    θ={tmed:.3f}／Rg={rmed:.1f}（旗89 と同一）"
              f"／θ高×Rg高 の春 {len(sp)} 日・秋 {len(au)} 日")
        skey, roi = PRIMARY[site]
        eras = load_gcc(gcc_files(a.phenocam_dir, skey, roi), COL)
        if not eras:
            print(f"    **主 ROI（{skey}/{roi}）の {COL} が読めない**"); continue
        print(f"    **主 ROI：{skey}/{roi}**（{len(eras)} 期）")
        out[site] = run_one("主", sp, au, eras)
        # 感度確認：別 ROI
        if site in SENSITIVITY:
            s2, r2 = SENSITIVITY[site]
            e2 = load_gcc(gcc_files(a.phenocam_dir, s2, r2), COL)
            if e2:
                print(f"    **感度確認：{s2}/{r2}**")
                v2 = run_one("感度", sp, au, e2)
                if v2 != out[site]:
                    print(f"      ※**主 ROI（{out[site]}）と食い違う（{v2}）"
                          f"——両方を記す。選ばない。**")
        # 感度確認：別の列
        for c in SENS_COLS:
            e3 = load_gcc(gcc_files(a.phenocam_dir, skey, roi), c)
            if e3:
                print(f"    **感度確認：列 {c}**")
                v3 = run_one(f"列{c}", sp, au, e3)
                if v3 != out[site]:
                    print(f"      ※**主指標（{out[site]}）と食い違う（{v3}）——両方を記す。**")

    print("\n  === 集計（事前登録の判定規則に当てる）===")
    lab = {"pheno": "緑の日だけ反転（フェノロジーで説明される）",
           "neither": "どちらでも反転しない（説明されない）",
           "both": "どちらでも反転する（旗89/90 と食い違う）",
           "brown_only": "枯れた日だけ反転（想定外）", None: "判定しない"}
    for s, v in out.items():
        print(f"    {s:<8}{lab.get(v, str(v))}")
    vals = [v for v in out.values() if v]
    n = len(vals)
    print("\n  === 結論 ===")
    if n < 2:
        print(f"  **判定しない**——判定できたサイトが {n} で 2 未満。")
    elif sum(v == "pheno" for v in vals) > n / 2:
        print("  **★フェノロジーで説明された**——**θ を揃えても、緑の日だけ反転する**。")
        print("  ＝A-3 を『**θ が高く、かつ林冠が緑のとき Bowen 反転が起きる**』と書き換える。")
        print("  **ただし「緑が原因」とは言わない**——**緑度を揃えると差が消える**までである。")
    elif sum(v == "neither" for v in vals) > n / 2:
        print("  **▲フェノロジーでも説明されない**——**旗88–91 と合わせ、")
        print("     観測データの層別による説明を打ち切ると確定する。**")
    elif sum(v == "both" for v in vals) > n / 2:
        print("  **○旗89/90 と食い違う**——**帯の作り方を疑う**（事前登録どおり）。")
    else:
        print("  **判定が割れた**——各サイトの行をそのまま記録し、**まとめない**。")
    print("\n  留保（事前登録どおり）：")
    print("   ・**独立クラスタは 2 つ**——US-Wkg と US-Whs は 10.45km で同一（Walnut Gulch）、")
    print("     US-SRM が別（Santa Rita）。**3 サイト≠3 反復。**")
    print("   ・**カメラの視野とフラックスのフェッチは同じではない**（旗81 と同型）。")
    print("   ・**GCC は色であって光合成ではない**／**濡れと照度に応答する**。")
    print("   ・**PhenoCam とタワーの対応は名前で付けた・座標未確認**（旗51/79）。")


if __name__ == "__main__":
    main()
