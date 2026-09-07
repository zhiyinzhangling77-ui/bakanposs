"""旗136：**旗85 の総当たりを、AmeriFlux の外へ広げる**（下調べ・**検定はしない**）。

## なぜこれを見るのか

**旗86 本体は「★で選ばない同一地点の対」を実際に測り、判定できたのは 8 組（地温だけなら 6 組）だった。**
**旗87 はその n を増やそうとして「COSORE の中では救えない」で閉じた**
（温度列は本当に無い・別名でもない・別ファイルにも無い）。
**＝「n を増やす」手は、COSORE の側では尽きている。**

**だが旗85 は自分の留保として、こう書いていた**：

> **AmeriFlux は南北アメリカだけ**＝**総当たりはアメリカ大陸について完全、他は不完全**。
> 実際 `d20190607_RUEHR`（スイス 47.48/8.37）は**10km 以内に AmeriFlux が無く**、この表に出ていない。

**＝タワー側の登録簿が足りていない。** **これは COSORE の限界ではなく、私の照合先の限界である。**
**手元（`/mnt/hdd`）には AmeriFlux 以外に JapanFlux・ChinaFlux・KoFlux がある。**
**旗85 の総当たりは、それらを一度も見ていない。**

## **この道具は実行可能性だけを出す。★も相関も一致率も計算しない。**

**それは検定の答えそのものであり、事前登録の前に見てはいけない**（旗94/98/101/103/108 と同じ作法）。
**出すのは距離・地温列の有無・日数・年数・タワー側の年範囲だけである。**

## 門①（対照）

**新しい登録簿を足す前に、AmeriFlux だけで旗85 を再現する。**
**旗85 の既知の数（`10km 以内のチャンバー 45 件`）が出なければ、この道具は距離を間違えている。**
**再現してから、はじめて他の登録簿を足す。**

## 既知の落とし穴（この道具が踏まないようにしたもの）

- **ファイル名で当てにいかない**（旗64・旗82・旗86 の欠陥19/21/23 と同型）。
  JapanFlux の BADM は **82 サイト中 81 件が `Site_General_Info`（下線）、1 件だけ `Site-General-Info`
  （ハイフン）** である。ハイフンだけで拾うと **1 件しか取れない**。
- **BADM が 2 系統ある**（旗86 追記②）。JapanFlux には `BADM/` と `BADM_corrected/`（zip・83 件）があり、
  **どちらを採るかをファイルの並び順で決めない**。**両方読んで、食い違う座標を必ず報告する。**
- **読めなかったものを黙って飛ばさない**。**登録簿ごとに「何件読めて何件落ちたか」を出す。**

    .venv/bin/python research/global_colocation_step136.py
"""
from __future__ import annotations

import argparse
import glob
import io
import re
import sys
import warnings
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore
from colocate_step51 import haversine
from same_site_arc_step66 import PAIRS
from runlog import tee_stdout

warnings.filterwarnings("ignore", category=UserWarning)

COSORE = "/mnt/hdd/cosore-0.7.0"
AMF_BADM = "/mnt/hdd/AmeriFlux_FLUXNET"
JPF_BADM = "/mnt/hdd/JAPANFLUX/BADM"
JPF_BADM_CORR = "/mnt/hdd/JAPANFLUX/BADM_corrected"
KRF_BADM = "/mnt/hdd/KoFlux/BADM"
CNF_ROOT = "/mnt/hdd/ChinaFlux"
JPF_DATA = "/mnt/hdd/JAPANFLUX"
KRF_DATA = "/mnt/hdd/KoFlux/EC_data"

PAIRED = {ds for _, ds, _ in PAIRS}          # 旗66/86 で既に対にしたチャンバー

# 旗85 が AmeriFlux 単独で出した既知の数（門①の合格条件）
STEP85_WITHIN_10KM = 45


# ---------------------------------------------------------------- 登録簿を読む

def _wide_badm(path_or_buf, name=""):
    """**テンプレート型の BADM**（`Variable / Description / Units / dataValue1..`）を 1 行に畳む。

    JapanFlux・KoFlux が配るのはこの形で、**AmeriFlux の long 型
    （`SITE_ID / VARIABLE / DATAVALUE`）とは別物**である。
    """
    sh = pd.read_excel(path_or_buf, sheet_name=0)
    if sh.empty or str(sh.columns[0]).strip().lower() != "variable":
        return None
    vcols = [c for c in sh.columns
             if str(c).lower().startswith("datavalue") or str(c).startswith("Unnamed")]
    if not vcols:
        return None
    idx = sh.set_index(sh.columns[0])

    def get(key):
        if key not in idx.index:
            return None
        row = idx.loc[key, vcols]
        if isinstance(row, pd.DataFrame):        # 同名の変数が複数行ある場合
            row = row.iloc[0]
        vals = [v for v in row.tolist() if pd.notna(v)]
        return vals[0] if vals else None

    sid = get("SITE_ID")
    lat, lon = pd.to_numeric(get("LOCATION_LAT"), errors="coerce"), \
        pd.to_numeric(get("LOCATION_LONG"), errors="coerce")
    if sid is None or not (np.isfinite(lat) and np.isfinite(lon)):
        return None
    return {"site": str(sid).strip(), "lat": float(lat), "lon": float(lon),
            "igbp": get("IGBP"), "src": name}


def read_template_dir(root, pattern=r"Site[-_]General[-_]Info", label=""):
    """**ファイル名の綴りに依らず**、テンプレート型 BADM を全部読む。"""
    files = [p for p in glob.glob(f"{root}/**/*.xlsx", recursive=True)
             if "__MACOSX" not in p and re.search(pattern, p, re.I)]
    rows, bad = [], []
    for p in sorted(files):
        try:
            r = _wide_badm(p, Path(p).name)
        except Exception as e:
            bad.append((Path(p).name, f"{type(e).__name__}: {str(e)[:40]}")); continue
        if r is None:
            bad.append((Path(p).name, "SITE_ID か座標が無い"))
        else:
            rows.append(r)
    print(f"  [{label}] xlsx {len(files)} 件 → 座標が取れた {len(rows)} 件・落ちた {len(bad)} 件")
    for n, why in bad[:6]:
        print(f"      落ちた：{n} — {why}")
    return pd.DataFrame(rows)


def read_zip_dir(root, pattern=r"Site[-_]General[-_]Info", label=""):
    """**zip の中の BADM** を読む（JapanFlux の `BADM_corrected`）。"""
    zips = sorted(glob.glob(f"{root}/**/*.zip", recursive=True))
    rows, bad = [], []
    for zp in zips:
        try:
            with zipfile.ZipFile(zp) as z:
                names = [n for n in z.namelist()
                         if n.lower().endswith(".xlsx") and "__MACOSX" not in n
                         and re.search(pattern, n, re.I)]
                for n in names:
                    r = _wide_badm(io.BytesIO(z.read(n)), Path(n).name)
                    if r is not None:
                        rows.append(r)
                    else:
                        bad.append((Path(zp).name, f"{Path(n).name}: 座標が無い"))
        except Exception as e:
            bad.append((Path(zp).name, f"{type(e).__name__}: {str(e)[:40]}"))
    print(f"  [{label}] zip {len(zips)} 件 → 座標が取れた {len(rows)} 件・落ちた {len(bad)} 件")
    for n, why in bad[:6]:
        print(f"      落ちた：{n} — {why}")
    return pd.DataFrame(rows)


def read_long_badm(root, label=""):
    """**long 型 BADM**（AmeriFlux）。旗85 `read_all_sites` と同じ読み方。"""
    cands = sorted([p for p in Path(root).rglob("*")
                    if p.is_file() and p.suffix.lower() in (".xlsx", ".xls", ".csv")
                    and ("BIF" in p.name.upper() or "BADM" in p.name.upper())])
    parts, n_read, bad = [], 0, []
    for f in cands:
        try:
            frames = (list(pd.read_excel(f, sheet_name=None).values())
                      if f.suffix.lower() in (".xlsx", ".xls")
                      else [pd.read_csv(f, low_memory=False, encoding="latin-1",
                                        encoding_errors="replace")])
        except Exception as e:
            bad.append((f.name, f"{type(e).__name__}")); continue
        used = False
        for df in frames:
            cols = {str(c).upper(): c for c in df.columns}
            sid = next((cols[k] for k in ("SITE_ID", "SITEID") if k in cols), None)
            var = next((cols[k] for k in ("VARIABLE", "VARIABLE_NAME") if k in cols), None)
            val = next((cols[k] for k in ("DATAVALUE", "VALUE") if k in cols), None)
            if not (sid and var and val):
                continue
            v = df[var].astype(str).str.upper()
            keep = df[v.isin(["LOCATION_LAT", "LOCATION_LONG", "IGBP"])]
            if keep.empty:
                continue
            parts.append(pd.DataFrame({"site": keep[sid].astype(str),
                                       "var": v[keep.index], "val": keep[val]}))
            used = True
        n_read += bool(used)
    print(f"  [{label}] BADM 候補 {len(cands)} 件 → 座標を含んだ {n_read} 件・開けなかった {len(bad)} 件")
    for n, why in bad[:6]:
        print(f"      開けない：{n} — {why}")
    if not parts:
        return pd.DataFrame()
    allp = pd.concat(parts, ignore_index=True)
    piv = (allp.drop_duplicates(subset=["site", "var"], keep="first")
              .pivot(index="site", columns="var", values="val"))
    out = pd.DataFrame({"site": piv.index,
                        "lat": pd.to_numeric(piv.get("LOCATION_LAT"), errors="coerce"),
                        "lon": pd.to_numeric(piv.get("LOCATION_LONG"), errors="coerce"),
                        "igbp": piv.get("IGBP").values if "IGBP" in piv else None})
    out["src"] = label
    return out.dropna(subset=["lat", "lon"]).reset_index(drop=True)


def compare_versions(a, b, label_a, label_b, tol_km=0.05):
    """**2 系統の BADM の座標を突き合わせる**（旗86 追記② の作法＝黙ってどちらかを採らない）。"""
    if a.empty or b.empty:
        print(f"  （{label_a} と {label_b} の比較は片方が空なので行わない）")
        return
    m = a.set_index("site")[["lat", "lon"]].join(
        b.set_index("site")[["lat", "lon"]], lsuffix="_a", rsuffix="_b", how="inner")
    if m.empty:
        print(f"  （{label_a} と {label_b} に共通サイトが無い）")
        return
    d = np.array([haversine(r.lat_a, r.lon_a, r.lat_b, r.lon_b) for r in m.itertuples()])
    n_diff = int((d > tol_km).sum())
    print(f"  **{label_a} 対 {label_b}**：共通 {len(m)} サイト／"
          f"{tol_km}km を超えて食い違う {n_diff} 件（最大 {d.max():.3f} km）")
    if n_diff:
        for i in np.argsort(-d)[:8]:
            print(f"      {m.index[i]}: {d[i]:.3f} km")
        print(f"  **10km の判定を跨ぐ差があれば『判定できない』として扱う**（旗86 の取り決め）。")


# ------------------------------------------------------------ タワー側の年範囲

def tower_years():
    """**タワー側のデータ年範囲**を、置いてあるファイル名から拾う。**中身は開かない。**

    **推測で埋めない**——年が名前から読めないものは `None` を返し、そう印字する。
    """
    yrs = {}
    for p in glob.glob(f"{JPF_DATA}/FLX_*"):
        m = re.search(r"FLX_([A-Z]{2}-\w+)_JapanFLUX\d+_ALLVARS_(\d{4})_(\d{4})", Path(p).name)
        if m:
            yrs[m.group(1)] = (int(m.group(2)), int(m.group(3)))
    for p in glob.glob(f"{KRF_DATA}/data-*.xlsx"):
        m = re.search(r"data-([A-Z]{2}-\w+)-(\d{4})-(\d{4})", Path(p).name)
        if m:
            yrs[m.group(1)] = (int(m.group(2)), int(m.group(3)))
    for p in sorted(glob.glob(f"{CNF_ROOT}/*/*/")):
        code = Path(p).parts[3]
        ys = [int(d) for d in (Path(x).name for x in glob.glob(f"{p}/*"))
              if re.fullmatch(r"\d{4}", d)]
        if ys:
            yrs[code] = (min(ys), max(ys))
    return yrs


# ------------------------------------------------------------------ 本体

def main():
    ap = argparse.ArgumentParser(description="旗136：総当たりを AmeriFlux の外へ広げる")
    ap.add_argument("--km", type=float, default=10.0)
    ap.add_argument("--cosore-dir", default=COSORE)
    a = ap.parse_args()

    tee_stdout("step136")
    print("=== 旗136：旗85 の総当たりを AmeriFlux の外へ広げる（下調べ・検定はしない）===")
    print("  **★も相関も一致率も計算しない**——事前登録の前に答えを見ないため（旗94/98/108 の作法）。")
    print("  出すのは **距離・地温列の有無・日数・年数・タワー側の年範囲** だけ。\n")

    # --- 1. チャンバー側
    desc = pd.read_csv(Path(a.cosore_dir) / "description.csv")
    ch = desc[["CSR_DATASET", "CSR_LATITUDE", "CSR_LONGITUDE", "CSR_IGBP"]].copy()
    ch["lat"] = pd.to_numeric(ch["CSR_LATITUDE"], errors="coerce")
    ch["lon"] = pd.to_numeric(ch["CSR_LONGITUDE"], errors="coerce")
    ch = ch.dropna(subset=["lat", "lon"]).reset_index(drop=True)
    n_coord = len(ch)
    # **旗85 と同じ母集団に揃える**——旗85 の `main` は
    # `if not (datasets/data_{ds}.csv).exists(): continue` で**時系列の無い登録を落としていた**。
    # **この 1 行を揃えないと門①は通らない**（第1版は揃えず 48 対 45 で外した）。
    root0 = Path(a.cosore_dir)
    ch["has_csv"] = [(root0 / "datasets" / f"data_{x}.csv").exists()
                     for x in ch["CSR_DATASET"]]
    n_nocsv = int((~ch["has_csv"]).sum())
    print(f"  COSORE：座標のあるデータセット **{n_coord} 件**／"
          f"うち時系列 csv があるのは **{n_coord - n_nocsv} 件**"
          f"（**{n_nocsv} 件は description にだけ在って csv が無い**）")
    print(f"      csv の無い登録：{sorted(ch.loc[~ch['has_csv'], 'CSR_DATASET'])}")
    print(f"  **旗85 の母集団は後者（csv のある分）である。門①はそちらで比べる。**\n")
    ch = ch[ch["has_csv"]].reset_index(drop=True)

    # --- 2. 登録簿を読む
    print("--- タワー側の登録簿 ---")
    amf = read_long_badm(AMF_BADM, "AmeriFlux")
    jpf = read_template_dir(JPF_BADM, label="JapanFlux BADM")
    jpc = read_zip_dir(JPF_BADM_CORR, label="JapanFlux BADM_corrected")
    krf = read_template_dir(KRF_BADM, label="KoFlux BADM")
    print()
    print("--- 2 系統の突き合わせ（旗86 追記② の作法）---")
    compare_versions(jpf, jpc, "JapanFlux BADM", "BADM_corrected")

    # **ChinaFlux は座標が手元に無い。黙って飛ばさず、穴として名指しする。**
    cn_dirs = [d for d in sorted(Path(CNF_ROOT).iterdir()) if d.is_dir()]
    print()
    print(f"  **[ChinaFlux] サイト {len(cn_dirs)} 件があるが、座標を持つファイルが手元に無い。**")
    print(f"      各サイトにあるのは `data description(en/zh).docx`（変数辞書）と年別 xlsx だけで、")
    print(f"      **緯度経度は入っていない**（本周に中身を確認した）。")
    print(f"      → **ChinaFlux は総当たりに入れられない。穴として記録し、GATE に積む。**")

    reg = pd.concat([d for d in (amf, jpf, krf) if not d.empty], ignore_index=True)
    reg = reg.drop_duplicates(subset=["site"], keep="first").reset_index(drop=True)
    print()
    print(f"  **総当たりに使う登録簿：{len(reg)} サイト**"
          f"（AmeriFlux {len(amf)} ＋ JapanFlux {len(jpf)} ＋ KoFlux {len(krf)}"
          f"／重複を除いた後）\n")

    # --- 3. 門①（対照）：AmeriFlux だけで旗85 を再現する
    print("--- 門①（対照）：AmeriFlux だけで旗85 を再現する ---")
    d_amf = np.array([[haversine(c.lat, c.lon, t.lat, t.lon) for t in amf.itertuples()]
                      for c in ch.itertuples()])
    near_amf = d_amf.min(axis=1)
    n_within = int((near_amf <= a.km).sum())
    print(f"  {a.km:.0f}km 以内にタワーを持つチャンバー：**{n_within} 件**"
          f"（旗85 は **{STEP85_WITHIN_10KM} 件**）")
    if n_within == STEP85_WITHIN_10KM:
        print("  → **一致。門①合格。距離の計算と BADM の読み方は旗85 と同じである。**\n")
    else:
        print(f"  → **不一致（差 {n_within - STEP85_WITHIN_10KM:+d} 件）。**")
        print("     **一致するまで、以下の新しい数を結論に使ってはいけない。**")
        print("     考えうる原因：BADM の版が違う（旗86 追記②）／COSORE の版が違う／距離式が違う。\n")

    # --- 4. 全登録簿での総当たり
    d_all = np.array([[haversine(c.lat, c.lon, t.lat, t.lon) for t in reg.itertuples()]
                      for c in ch.itertuples()])
    k = d_all.argmin(axis=1)
    ch["near_km"] = d_all.min(axis=1)
    ch["near_tower"] = reg["site"].to_numpy()[k]
    ch["near_src"] = reg["src"].to_numpy()[k]
    ch["amf_km"] = near_amf
    ch["paired"] = ch["CSR_DATASET"].isin(PAIRED)

    within = ch[ch["near_km"] <= a.km]
    new = within[(within["amf_km"] > a.km) & (~within["paired"])]
    print("--- 全登録簿での総当たり ---")
    print(f"  {a.km:.0f}km 以内：**{len(within)} 件**"
          f"（AmeriFlux だけなら {n_within} 件＝**{len(within) - n_within} 件増えた**）")
    print(f"  そのうち **旗66/86 で既に対にしていない、かつ AmeriFlux では届かなかった** もの："
          f"**{len(new)} 件**\n")

    if new.empty:
        print("  **新しく対を作れるチャンバーは無い。**")
        print("  ＝**旗86 の n=8 は、タワー側の登録簿を広げても増えない。**")
    else:
        print("  | チャンバー | 生態系 | タワー | 登録簿 | km | AmeriFlux 最近傍 km |")
        print("  |---|---|---|---|---|---|")
        for r in new.sort_values("near_km").itertuples():
            print(f"  | {r.CSR_DATASET} | {r.CSR_IGBP} | {r.near_tower} | {r.near_src} "
                  f"| {r.near_km:.2f} | {r.amf_km:.0f} |")

    # --- 5. 新しい候補が「検定できる」かの下調べ（★は出さない）
    print()
    print("--- 新しい候補の実行可能性（**★は出さない**）---")
    yrs = tower_years()
    root = Path(a.cosore_dir)
    print("  | チャンバー | 地温列 | 水分列 | 日数 | チャンバー年 | タワー年 | 重なる年 |")
    print("  |---|---|---|---|---|---|---|")
    ok_pairs = 0
    for r in new.sort_values("near_km").itertuples():
        f = root / "datasets" / f"data_{r.CSR_DATASET}.csv"
        if not f.exists():
            print(f"  | {r.CSR_DATASET} | **csv が無い** | | | | | |"); continue
        try:
            df, st, sm = load_cosore(f, None)
        except Exception as e:
            print(f"  | {r.CSR_DATASET} | **読めない**（{type(e).__name__}）| | | | | |"); continue
        if df is None or df.empty:
            print(f"  | {r.CSR_DATASET} | **空** | | | | | |"); continue
        cy = (int(df.index.min().year), int(df.index.max().year))
        ndays = int(df.index.normalize().nunique())
        ty = yrs.get(r.near_tower)
        ov = (max(cy[0], ty[0]), min(cy[1], ty[1])) if ty else None
        ov_s = f"{ov[0]}–{ov[1]}" if ov and ov[0] <= ov[1] else ("**重ならない**" if ty else "**不明**")
        print(f"  | {r.CSR_DATASET} | {st or '**無し**'} | {sm or '無し'} | {ndays} "
              f"| {cy[0]}–{cy[1]} | {f'{ty[0]}–{ty[1]}' if ty else '**不明**'} | {ov_s} |")
        if st and ndays >= 60 and ty and ov[0] <= ov[1]:
            ok_pairs += 1
    # --- 5b. **届かなかった非米州チャンバー**も名指しする（黙って消さない）
    print()
    print("--- 非米州のチャンバーで、どの登録簿にも届かなかったもの ---")
    print("  **「登録簿を広げれば届く」わけではないことを、残りの距離で示す。**")
    nonam = ch[~((ch["lon"] > -170) & (ch["lon"] < -30))]
    far = nonam[nonam["near_km"] > a.km].sort_values("near_km")
    print("  | チャンバー | 生態系 | 最近傍タワー | km |")
    print("  |---|---|---|---|")
    for r in far.itertuples():
        print(f"  | {r.CSR_DATASET} | {r.CSR_IGBP} | {r.near_tower} | {r.near_km:.0f} |")
    print(f"  **{len(far)} 件**。**手元に ICOS（欧州）・OzFlux（豪）・AsiaFlux・ChinaFlux の"
          f"座標が無い**ことが、そのまま残っている。")

    # --- 6. タワー側が読めるか（旗80 と同じ確かめ方・**「列が在る」で止めない**）
    print()
    print("--- タワー側が読めるか（旗86 の教訓：**列が在ることと中身が在ることは別**）---")
    tower_ok = {}
    try:
        from japanflux_pn.sites import get_site, RK_VARS
        from japanflux_pn.config import AnalysisConfig
        from japanflux_pn.preprocess import find_corevars_files, load_raw_all
    except Exception as e:
        print(f"  **旗80 の読み込み器を import できない**（{type(e).__name__}: {e}）")
        get_site = None
    if get_site is not None:
        for code in sorted(set(new["near_tower"])):
            print(f"  ━━ {code} ━━")
            tower_ok[code] = False
            try:
                sp = get_site(code)
                files = find_corevars_files(sp)
            except Exception as e:
                print(f"    **登録が引けない／HH が無い**：{type(e).__name__}: {str(e)[:90]}")
                continue
            head = pd.read_csv(files[0], nrows=2)
            vm = sp.var_map()
            miss = {k: v for k, v in vm.items() if v not in set(head.columns)}
            print(f"    HH {len(files)} 件：{files[0].name}")
            if miss:
                print(f"    **写像した列が無い**：{miss}")
                for k, v in miss.items():
                    alt = [c for c in head.columns
                           if c.split("_")[0] == v.split("_")[0] and not c.endswith("_QC")]
                    print(f"       {k}（{v}）の代案：{alt[:6]}")
            else:
                print("    **写像した列はすべて在る**")
            try:
                raw = load_raw_all(sp, AnalysisConfig())
            except Exception as e:
                print(f"    **読み込み失敗**：{type(e).__name__}: {str(e)[:110]}")
                continue
            n_ok = {v: int(raw[v].notna().sum()) for v in ("GER", "Ts", "th") if v in raw}
            print(f"    読み込み成功：{len(raw):,} 行／"
                  f"{raw.index.min():%Y-%m}〜{raw.index.max():%Y-%m}／有効数 {n_ok}")
            empty = [v for v in ("GER", "Ts", "th") if n_ok.get(v, 0) <= 1000]
            for v in empty:
                # **「無い」と言う前に、名前が違うだけでないかを確かめる**（旗87 と同型）。
                pat = {"Ts": r"TS|TEMP|SOIL", "th": r"SWC|VWC|MOIST|WATER|WTD",
                       "GER": r"RECO|ER_"}[v]
                cand = [c for c in head.columns if re.search(pat, c, re.I)]
                print(f"    ※**{v} は有効値 {n_ok.get(v, 0)}**／"
                      f"ヘッダ中で `{pat}` に当たる列：{cand if cand else '**一つも無い**'}")
            tower_ok[code] = not empty
            if empty:
                print(f"    → **{code} は旗66 の検出器に掛けられない**"
                      f"（`var_overrides` で救えるのは名前違いのときだけ）。")

    print()
    print(f"  **チャンバー側が実行可能な候補：{ok_pairs} 件**"
          "（地温あり・60 日以上・タワーと年が重なる）")
    usable = [r.CSR_DATASET for r in new.itertuples() if tower_ok.get(r.near_tower)]
    print(f"  **そのうちタワー側も同じ検出器に掛けられるもの：{len(usable)} 件** {usable}")
    print(f"  **＝旗86 の『★で選ばない群 n=8（地温だけなら 6）』に足せるのは {len(usable)} 件。**")
    print("  **★率も一致率もここでは出さない。事前登録してから測る。**")


if __name__ == "__main__":
    main()
