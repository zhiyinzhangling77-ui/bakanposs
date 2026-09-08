#!/usr/bin/env python3
"""旗138 の事前登録のための **前提の事実確認**——**「位相が逆である」を測る**。

旗137 は `ES-FcO`（Finca Oran・IGBP=CRO・冬穀）について、
**「春に作物があり、秋は刈り跡である＝A-3 の 4 クラスタと位相が逆」** と書いた。
**その根拠は `Oran_Site_Info.docx` の逐語（"crop growth cycles … December to June"）だけで、
データで測っていない。** 旗97 の作法（**事前登録の前に、前提にする事実を分布で確かめる**）に従い、
**事前登録を書く前に測る。**

**この道具は A-3 の量（θ→γH・θ→γLE）を一切計算しない。** 測るのは **植生／炭素の季節位相だけ** である。

--------------------------------------------------------------------------
## 測る量（実行前に固定）

サイトごと・季節ごと（**MAM=3–5 月 / SON=9–11 月**、旗88/107 と同じ切り方）に：

- **緑度**：**MOD13Q1 NDVI**（`pixel_reliability == 0` の良質のみ）。**Oran にしか無い。**
- **炭素**：**日次の正味 CO2 フラックス**。**負が吸収。**
  - Oran：30 分 `FC_mass`（`FC_QC ≤ 3`＝Foken の高品質）を日平均。**48 本中 24 本以上**（旗137 と同じ被覆規則）。
  - US-Wkg/US-Whs/US-SRM：FLUXNET `NEE_VUT_REF`（`NEE_VUT_REF_QC ≥ 0.5`）。
  - **単位が違う**（Oran は質量系・FLUXNET は gC m-2 d-1）。
    **だから比較は「そのサイトの中で、どちらの季節の吸収が強いか」という順序だけに使う。**
    **サイト間で大きさを比べない。**

**位相の統計量**：**Δ = (SON の平均) − (MAM の平均)**。**負なら「秋のほうが吸収が強い＝秋緑」。**
**年ごとに出し、年をまたぐ中央値を代表値にする**（年の重みを揃えるため）。

--------------------------------------------------------------------------
## 門①（対照）——**実行前に宣言する。落ちたら「前提は立たなかった」と書く**

- **C0（符号の規約）**：各サイトの**月別平均**を印字し、**成長期に負（吸収）になっている**ことを確かめる。
  **全月が正、または成長期が正なら、符号の規約を取り違えている**＝その先を読まない。
- **C1（Oran・二つの独立な観測が一致するか）**：**MODIS NDVI の位相**と**タワー CO2 の位相**が
  **同じ向き（どちらも春が強い）** を指すこと。**食い違ったら「Oran の位相は確定しない」と書き、
  事前登録の識別検定の枠は仮定のままにする。**
- **C2（クラスタ内の一致）**：**US-Wkg と US-Whs は同じクラスタ（Walnut Gulch）**なので、
  **Δ の符号が一致すること。** 割れたら**この位相検出器はサイト水準で信用できない**＝A-3 側の位相は
  「判定しない」と書く（旗106/107 の作法：**割れたクラスタは判定しない**）。
- **C3（下限）**：**各サイト×季節のセルが 60 日以上・3 暦年以上**（旗105/107 の下限）。
  **満たさないセルは判定に使わない。**

**予測（先に書く・当たっても外れても記録する）**：
**H1＝Oran は Δ>0（春のほうが吸収が強い＝春緑）、US-Wkg/Whs/SRM は Δ<0（秋緑）。**
**確信は中程度。** 北米 3 サイトは夏のモンスーンで緑になるので **SON が MAM より強いだろう**と
考えているが、**SON の 11 月は既に枯れている**ので、**打ち消して Δ≈0 になる可能性がある。**

--------------------------------------------------------------------------
使い方:
    .venv/bin/python research/premise_phase_step138.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ORAN_HH = Path("/mnt/hdd/Eddy data in Spain/Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.csv")
ORAN_DD = Path("/mnt/hdd/Eddy data in Spain/Oran_EddyDaily_MASTER_2018_2020_correct.csv")
MOD13 = Path("/mnt/hdd/Dataset/MOD13Q1_NDVI_EVI/MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv")
AMF = Path("/mnt/hdd/AmeriFlux_FLUXNET")

AMF_SITES = {
    "US-Wkg": "AMF_US-Wkg_FLUXNET_2004-2025_v1.3_r1",
    "US-Whs": "AMF_US-Whs_FLUXNET_2007-2025_v1.3_r1",
    "US-SRM": "AMF_US-SRM_FLUXNET_2004-2025_v1.3_r1",
}

MAM = (3, 4, 5)
SON = (9, 10, 11)
MIN_DAYS = 60
MIN_YEARS = 3


def season_of(month: int) -> str | None:
    if month in MAM:
        return "MAM"
    if month in SON:
        return "SON"
    return None


def phase_table(df: pd.DataFrame, value_col: str, label: str) -> dict:
    """df は列 `year`・`month`・value_col を持つ日次（または合成）データ。

    戻り値に **年ごとの季節平均・日数**と、**Δ の年またぎ中央値**を入れる。
    """
    d = df.dropna(subset=[value_col]).copy()
    d["season"] = d["month"].map(season_of)
    d = d[d["season"].notna()]

    rows = []
    for (yr, sea), g in d.groupby(["year", "season"]):
        rows.append({"year": int(yr), "season": sea, "mean": float(g[value_col].mean()), "n": int(len(g))})
    per = pd.DataFrame(rows)

    print(f"\n--- {label} : 年×季節 ---")
    if per.empty:
        print("  （データ 0 行）")
        return {"label": label, "delta_median": np.nan, "ok_floor": False, "n_years": 0}
    piv_m = per.pivot(index="year", columns="season", values="mean")
    piv_n = per.pivot(index="year", columns="season", values="n")
    for yr in piv_m.index:
        mam = piv_m.at[yr, "MAM"] if "MAM" in piv_m.columns else np.nan
        son = piv_m.at[yr, "SON"] if "SON" in piv_m.columns else np.nan
        nm = piv_n.at[yr, "MAM"] if "MAM" in piv_n.columns else np.nan
        ns = piv_n.at[yr, "SON"] if "SON" in piv_n.columns else np.nan
        dlt = son - mam if (pd.notna(mam) and pd.notna(son)) else np.nan
        print(f"  {int(yr)}  MAM={mam:8.3f} (n={nm if pd.isna(nm) else int(nm)})   "
              f"SON={son:8.3f} (n={ns if pd.isna(ns) else int(ns)})   Δ=SON-MAM={dlt:8.3f}")

    # C3: 下限（各季節の合計日数と、両季節が揃う年数）
    tot = per.groupby("season")["n"].sum().to_dict()
    both = piv_m.dropna().index
    n_years = len(both)
    ok_floor = (tot.get("MAM", 0) >= MIN_DAYS and tot.get("SON", 0) >= MIN_DAYS and n_years >= MIN_YEARS)
    deltas = (piv_m["SON"] - piv_m["MAM"]).dropna() if {"MAM", "SON"} <= set(piv_m.columns) else pd.Series(dtype=float)
    dmed = float(np.median(deltas)) if len(deltas) else np.nan

    print(f"  合計日数 MAM={tot.get('MAM', 0)} / SON={tot.get('SON', 0)}   両季節が揃う年={n_years}")
    print(f"  **Δ の年またぎ中央値 = {dmed:.3f}**   [{', '.join(f'{v:.3f}' for v in deltas)}]")
    print(f"  C3（60 日・3 年）: {'○満たす' if ok_floor else '×満たさない'}")
    return {"label": label, "delta_median": dmed, "ok_floor": ok_floor, "n_years": n_years,
            "deltas": list(map(float, deltas))}


def monthly_print(df: pd.DataFrame, value_col: str, label: str) -> None:
    """C0：符号の規約を目で確かめるための月別平均。"""
    m = df.dropna(subset=[value_col]).groupby("month")[value_col].agg(["mean", "size"])
    txt = "  ".join(f"{int(mo):02d}:{r['mean']:+.2f}(n={int(r['size'])})" for mo, r in m.iterrows())
    print(f"\n[C0] {label} 月別平均: {txt}")
    if len(m):
        print(f"     最小（最も吸収側）= {m['mean'].min():+.3f} @ {int(m['mean'].idxmin()):02d}月 / "
              f"最大 = {m['mean'].max():+.3f} @ {int(m['mean'].idxmax()):02d}月")


# ---------------------------------------------------------------- Oran タワー
def load_oran_flux() -> pd.DataFrame:
    """30 分 `FC_mass`（`FC_QC<=3`）を日平均に落とす。48 本中 24 本以上を要求。"""
    use = ["TIMESTAMP", "FC_mass", "FC_QC"]
    d = pd.read_csv(ORAN_HH, usecols=use, low_memory=False)
    n_raw = len(d)
    d["FC_mass"] = pd.to_numeric(d["FC_mass"], errors="coerce")
    d["FC_QC"] = pd.to_numeric(d["FC_QC"], errors="coerce")
    ts = pd.to_datetime(d["TIMESTAMP"], format="mixed", errors="coerce")
    d["date"] = ts.dt.floor("D")
    d = d[d["date"].notna()]

    # ★ゼロ埋めの検め（旗137 欠陥 #62 の作法）：合計・最大・最小を印字する
    v = d["FC_mass"].dropna()
    print(f"[在庫] Oran 30分 FC_mass: 全 {n_raw} 行 / 有効 {len(v)} 本 / "
          f"平均 {v.mean():.4f} 最小 {v.min():.4f} 最大 {v.max():.4f} 合計 {v.sum():.2f}")
    print(f"[在庫] FC_QC の分布: {d['FC_QC'].value_counts(dropna=False).sort_index().to_dict()}")

    good = d[(d["FC_QC"] <= 3) & d["FC_mass"].notna()]
    print(f"[在庫] FC_QC<=3 かつ値あり: {len(good)} 本（全体の {100*len(good)/n_raw:.1f}%）")

    agg = good.groupby("date")["FC_mass"].agg(["mean", "size"]).reset_index()
    day = agg[agg["size"] >= 24].copy()
    print(f"[在庫] 日平均に落ちた日数: {len(day)} 日（48本中24本以上の規則・落とした日 {len(agg)-len(day)}）")
    day["year"] = day["date"].dt.year
    day["month"] = day["date"].dt.month
    return day.rename(columns={"mean": "fc"})[["date", "year", "month", "fc"]]


def load_oran_ndvi_mod13() -> pd.DataFrame:
    d = pd.read_csv(MOD13, low_memory=False)
    d = d[d["ID"] == "Oran"].copy()
    rel = "MOD13Q1_061__250m_16_days_pixel_reliability"
    ndvi = "MOD13Q1_061__250m_16_days_NDVI"
    d["date"] = pd.to_datetime(d["Date"], errors="coerce")
    n_all = len(d)
    d = d[(d[rel] == 0) & d[ndvi].notna()]
    print(f"[在庫] MOD13Q1 Oran: 全 {n_all} 点 / pixel_reliability==0 は {len(d)} 点 / "
          f"期間 {d['date'].min().date()}..{d['date'].max().date()} / "
          f"NDVI 最小 {d[ndvi].min():.3f} 最大 {d[ndvi].max():.3f}")
    d["year"] = d["date"].dt.year
    d["month"] = d["date"].dt.month
    return d.rename(columns={ndvi: "ndvi"})[["date", "year", "month", "ndvi"]]


def load_oran_daily_ndvi() -> pd.DataFrame:
    """日次マスタの `NDVI_csm`（平滑済み＝**派生量**。参考にしか使わない）。"""
    d = pd.read_csv(ORAN_DD, low_memory=False)
    d["date"] = pd.to_datetime(d["Date"], errors="coerce")
    d = d[d["date"].notna()].copy()
    d["year"] = d["date"].dt.year
    d["month"] = d["date"].dt.month
    v = pd.to_numeric(d.get("NDVI_csm"), errors="coerce")
    d["ndvi_csm"] = v
    print(f"[在庫] 日次マスタ NDVI_csm: {int(v.notna().sum())} 日 / "
          f"最小 {v.min():.3f} 最大 {v.max():.3f}")
    return d[["date", "year", "month", "ndvi_csm"]]


# ---------------------------------------------------------------- AmeriFlux
def load_amf_nee(site: str, folder: str) -> pd.DataFrame:
    p = AMF / folder
    cands = sorted(p.glob("*FLUXMET_DD*.csv"))
    if not cands:
        print(f"[在庫] {site}: FLUXMET_DD が見つからない → 除外")
        return pd.DataFrame(columns=["year", "month", "nee"])
    # ★道具の欠陥 #63（1 走目で踏んだ・**落ちる方**）：
    #   `NEE_VUT_REF` が全サイトにあると決め打ちしていたが、**US-Whs には `NEE_CUT_*` しか無い。**
    #   FLUXNET の VUT/CUT は u* しきい値の決め方が違う（年ごと／記録全体）。
    #   **どちらを使ったかを必ず印字する**（混ぜたことを黙って隠さない）。
    head = pd.read_csv(cands[0], nrows=0)
    for cand in ("NEE_VUT_REF", "NEE_CUT_REF"):
        if cand in head.columns and f"{cand}_QC" in head.columns:
            col = cand
            break
    else:
        print(f"[在庫] {site}: NEE_VUT_REF も NEE_CUT_REF も無い（列 {len(head.columns)}）→ 除外")
        return pd.DataFrame(columns=["year", "month", "nee"])
    print(f"[在庫] {site}: 使う列 = **{col}**"
          f"（VUT の有無: {'有' if 'NEE_VUT_REF' in head.columns else '無'}）")
    d = pd.read_csv(cands[0], usecols=["TIMESTAMP", col, f"{col}_QC"], low_memory=False)
    d = d.rename(columns={col: "NEE_VUT_REF", f"{col}_QC": "NEE_VUT_REF_QC"})
    d = d.replace(-9999, np.nan)
    d["date"] = pd.to_datetime(d["TIMESTAMP"].astype(str), format="%Y%m%d", errors="coerce")
    d = d[d["date"].notna()]
    v = d["NEE_VUT_REF"].dropna()
    print(f"[在庫] {site}: 全 {len(d)} 日 / NEE 有効 {len(v)} / "
          f"平均 {v.mean():.3f} 最小 {v.min():.3f} 最大 {v.max():.3f} / "
          f"期間 {d['date'].min().date()}..{d['date'].max().date()}")
    print(f"       QC の分布（>=0.5 が採用）: 有効 {int(d['NEE_VUT_REF_QC'].notna().sum())} / "
          f"中央値 {d['NEE_VUT_REF_QC'].median():.2f} / >=0.5 は {int((d['NEE_VUT_REF_QC']>=0.5).sum())} 日")
    good = d[(d["NEE_VUT_REF_QC"] >= 0.5) & d["NEE_VUT_REF"].notna()].copy()
    good["year"] = good["date"].dt.year
    good["month"] = good["date"].dt.month
    return good.rename(columns={"NEE_VUT_REF": "nee"})[["date", "year", "month", "nee"]]


def oran_cell_floor() -> None:
    """**追補（位相を測ったあとに足した節・判定には使わない）**——事前登録の**下限の下調べ**。

    旗89 の切り方（**そのサイトの中央値より θ が高く、かつ Rg も高い日＝`θ高×Rg高`**）で、
    **ES-FcO の春と秋に何日残るか**を数える。**60 日・3 暦年（旗105/107）に届かなければ、
    事前登録の主判定はこのセルでは書けない。**
    **ここでは日数と θ の分布しか出さない。θ→γH・θ→γLE は計算しない。**
    """
    use = ["TIMESTAMP", "SWC_1_1_1", "SW_IN", "H", "H_QC", "LE", "LE_QC"]
    d = pd.read_csv(ORAN_HH, usecols=use, low_memory=False)
    for c in use[1:]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    ts = pd.to_datetime(d["TIMESTAMP"], format="mixed", errors="coerce")
    d["date"] = ts.dt.floor("D")
    d = d[d["date"].notna()]

    for c in ("SWC_1_1_1", "SW_IN", "H", "LE"):
        v = d[c].dropna()
        print(f"[在庫] Oran 30分 {c}: 有効 {len(v)} / 平均 {v.mean():.3f} / 最小 {v.min():.3f} / "
              f"最大 {v.max():.3f} / 合計 {v.sum():.1f}")

    day = d.groupby("date").agg(th=("SWC_1_1_1", "mean"), n_th=("SWC_1_1_1", "size"),
                                rg=("SW_IN", "mean"), n_rg=("SW_IN", "size")).reset_index()
    day = day[(day["n_th"] >= 24) & (day["n_rg"] >= 24)].dropna(subset=["th", "rg"])
    day["year"] = day["date"].dt.year
    day["month"] = day["date"].dt.month
    day["season"] = day["month"].map(season_of)

    th_med, rg_med = day["th"].median(), day["rg"].median()
    print(f"\n[下調べ] ES-FcO 日次: {len(day)} 日 / θ 中央値 {th_med:.4f} / Rg 中央値 {rg_med:.2f}")
    day["cell"] = (day["th"] > th_med) & (day["rg"] > rg_med)
    sub = day[day["season"].notna()]
    print("[下調べ] 季節ごとの θ の分布と、`θ高×Rg高` に残る日数:")
    for sea in ("MAM", "SON"):
        g = sub[sub["season"] == sea]
        q = g["th"].quantile([0.25, 0.5, 0.75]).values
        cell = g[g["cell"]]
        per_year = cell.groupby("year").size().to_dict()
        print(f"  {sea}: 全 {len(g)} 日 / θ 四分位 {q[0]:.4f}, {q[1]:.4f}, {q[2]:.4f} / "
              f"`θ高×Rg高` {len(cell)} 日 / 年別 {per_year}")
        ok = len(cell) >= MIN_DAYS and len(per_year) >= MIN_YEARS
        print(f"      下限（60 日・3 年）: {'○満たす' if ok else '×満たさない'}")

    # ★H・LE の外れ値と QC——**LE の生値に物理的にありえない値がある**（最小 −7860 W m-2）。
    #   **事前登録の前に、どれだけ落ちるかを数える。**
    for c, qc in (("H", "H_QC"), ("LE", "LE_QC")):
        v = d[c]
        n_have = int(v.notna().sum())
        n_qc = int(((d[qc] <= 3) & v.notna()).sum())
        n_rng = int(((d[qc] <= 3) & v.between(-200, 800)).sum())
        print(f"[下調べ] {c}: 有効 {n_have} → QC<=3 で {n_qc} → さらに −200..800 W m-2 で {n_rng} "
              f"（QC<=3 のうち範囲外 {n_qc - n_rng} 本）")

    dd = d[(d["H_QC"] <= 3) & d["H"].between(-200, 800)
           & (d["LE_QC"] <= 3) & d["LE"].between(-200, 800)]
    dcov = dd.groupby("date").size()
    keep = dcov[dcov >= 24].index
    kd = pd.DataFrame({"date": keep})
    kd["year"] = kd["date"].dt.year
    kd["month"] = kd["date"].dt.month
    kd["season"] = kd["month"].map(season_of)
    print("[下調べ] **H・LE を同時に QC<=3・範囲内にした日**（48 本中 24 本以上）:")
    for sea in ("MAM", "SON"):
        g = kd[kd["season"] == sea]
        per_year = g.groupby("year").size().to_dict()
        ok = len(g) >= MIN_DAYS and len(per_year) >= MIN_YEARS
        print(f"  {sea}: {len(g)} 日 / 年別 {per_year} / 下限: {'○満たす' if ok else '×満たさない'}")


def main() -> int:
    print("=" * 78)
    print("旗138 前提の事実確認：**位相が逆であることを測る**（A-3 の量は計算しない）")
    print("=" * 78)

    res = {}

    print("\n########## Oran（ES-FcO） ##########")
    oran_flux = load_oran_flux()
    monthly_print(oran_flux, "fc", "Oran タワー FC_mass（負が吸収）")
    res["oran_flux"] = phase_table(oran_flux, "fc", "Oran タワー CO2 フラックス")

    oran_ndvi = load_oran_ndvi_mod13()
    monthly_print(oran_ndvi, "ndvi", "Oran MOD13Q1 NDVI（大きいほど緑）")
    res["oran_ndvi"] = phase_table(oran_ndvi, "ndvi", "Oran MOD13Q1 NDVI")

    oran_csm = load_oran_daily_ndvi()
    res["oran_csm"] = phase_table(oran_csm, "ndvi_csm", "Oran 日次マスタ NDVI_csm（派生量・参考）")

    print("\n########## A-3 の北米クラスタ ##########")
    for site, folder in AMF_SITES.items():
        d = load_amf_nee(site, folder)
        if d.empty:
            continue
        monthly_print(d, "nee", f"{site} NEE_VUT_REF（負が吸収）")
        res[site] = phase_table(d, "nee", f"{site} NEE_VUT_REF")

    # ---------------------------------------------------------------- 門①の判定
    print("\n" + "=" * 78)
    print("門①（対照）の判定")
    print("=" * 78)

    # C1：Oran で、NDVI の位相と CO2 の位相が同じ向きを指すか。
    #     NDVI は「大きいほど緑」・CO2 は「負が吸収」なので、
    #     **春が緑** ⇔ NDVI の Δ(SON-MAM) < 0 かつ CO2 の Δ(SON-MAM) > 0。
    dn = res["oran_ndvi"]["delta_median"]
    df_ = res["oran_flux"]["delta_median"]
    c1_ndvi_spring = dn < 0
    c1_flux_spring = df_ > 0
    c1 = bool(c1_ndvi_spring == c1_flux_spring)
    print(f"[C1] Oran: NDVI Δ={dn:+.4f}（{'春が緑' if c1_ndvi_spring else '秋が緑'}） / "
          f"CO2 Δ={df_:+.4f}（{'春の吸収が強い' if c1_flux_spring else '秋の吸収が強い'}） → "
          f"{'○一致' if c1 else '×食い違い＝Oran の位相は確定しない'}")

    # C2：Walnut Gulch の 2 サイトで Δ の符号が一致するか
    if "US-Wkg" in res and "US-Whs" in res:
        a, b = res["US-Wkg"]["delta_median"], res["US-Whs"]["delta_median"]
        c2 = bool(np.sign(a) == np.sign(b))
        print(f"[C2] Walnut Gulch: US-Wkg Δ={a:+.3f} / US-Whs Δ={b:+.3f} → "
              f"{'○符号一致' if c2 else '×割れた＝このクラスタは判定しない'}")
    else:
        c2 = False
        print("[C2] Walnut Gulch: 片方が読めず判定不能")

    # C3：下限
    print("[C3] 下限（各季節 60 日・両季節が揃う年 3 以上）:")
    for k, v in res.items():
        print(f"     {v['label']}: {'○' if v['ok_floor'] else '×'}（年={v['n_years']}）")

    # ---------------------------------------------------------------- まとめ
    print("\n" + "=" * 78)
    print("前提「位相が逆である」の判定")
    print("=" * 78)
    print("  ※ 単位が違うので、比べるのは **各サイト内での季節の順序（Δ の符号）だけ**。")
    for k in ["oran_flux", "US-Wkg", "US-Whs", "US-SRM"]:
        if k not in res:
            continue
        v = res[k]
        d = v["delta_median"]
        if not np.isfinite(d):
            phase = "判定不能"
        else:
            phase = "秋のほうが吸収が強い（秋緑）" if d < 0 else "春のほうが吸収が強い（春緑）"
        print(f"  {v['label']:38s} Δ={d:+8.3f}  {phase}  下限{'○' if v['ok_floor'] else '×'}")
    print("\n  ★この道具は位相しか測っていない。A-3 の量（θ→γH・θ→γLE）は一切計算していない。")

    print("\n" + "=" * 78)
    print("追補：事前登録の**下限の下調べ**（位相の判定には使わない）")
    print("=" * 78)
    oran_cell_floor()
    return 0


if __name__ == "__main__":
    sys.exit(main())
