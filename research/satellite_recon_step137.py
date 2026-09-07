#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""旗137 の下調べ道具 — /mnt/hdd/Dataset の Oran/TzM 衛星スタックで
A-3（θ→蒸発・春秋の非対称）に手が届くかを **実行可能性だけ** 測る。

**この道具は A-3 の量（θ→γH・θ→γLE の Δ）を一切計算しない。**
それは事前登録を書いてからやる（旗107 の作法）。ここで測るのは次の 5 つだけ:

  §1 身元と座標    — スタックの抽出点が誰なのか。タワーとの距離。
  §2 タワー側の在庫 — 旗36/旗107 が要る量（th=SWC・gH=H・gLE=LE・Rg=SW_IN・P・Ts）が
                     何日あるか。**下限は旗105/107 と同じ「各セル 60 日・3 暦年」。**
  §3 前提の穴      — 作付け・灌漑（CRO は A-3 の 4 クラスタと生態系が違う）。
  §4 衛星側の在庫   — 製品ごとの被覆・刻み・タワー年との重なり。
  §5 門①（対照）   — **下の 3 本は実行前に宣言してある。**

門①の宣言（結果を見る前に書いた）:
  C1（陽性対照・読み取りが正しいか）:
      Oran のタワー日平均 SWC_1_1_1 と SMAP sm_surface 日平均（2018-2020）の Spearman ρ。
      **期待 ρ > +0.5。** 同じ点の近傍土壌水分を二つの方法で測っているので、
      これが崩れるなら日付の対応かサイトの取り違えを疑う。
  C2（整合の鋭い側・季節周期を抜く）:
      同じ対の **1 階差分** の ρ、および衛星側を +45 日ずらしたときの ρ。
      **差分 ρ > 0 かつ ずらし ρ が C1 より明確に低い**なら、季節周期ではなく
      日々の事象で一致している。そうでなければ「季節周期しか合っていない」と書く。
  C3（Fast_*.csv の行順の同定）:
      Fast_ 系 csv は `system:index` が 0/1 だけで **サイト名の列が無い**。
      Fast_OranTzM_MOD16.csv の 2 行を、名前つきの AppEEARS MOD16 結果（2021-2024）と
      日付で突き合わせ、どちらが Oran かを相関で決める。
      **両者の差が小さければ「同定できない」と書く**（推測で割り当てない）。

出力はすべて印字する。数値をこのファイルに書き戻さない。
"""

from __future__ import annotations

import csv
import json
import math
import os
import re
import sys
import zipfile
from collections import defaultdict
from datetime import date, datetime, timedelta

DATASET = "/mnt/hdd/Dataset"
EDDY = "/mnt/hdd/Eddy data in Spain"
SITEINFO = "/mnt/hdd/Site Info"

TOWER_CSV = os.path.join(EDDY, "Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.csv")
MASTER_CSV = os.path.join(EDDY, "Oran_EddyDaily_MASTER_2018_2020_correct.csv")
SMAP_CSV = os.path.join(DATASET, "SMAP_OranTzM.csv")

# 旗105/107 の下限。ここでは変えない。
MIN_DAYS_PER_CELL = 60
MIN_YEARS = 3

NA = {"", "nan", "NaN", "NAN", "-9999", "-9999.0", "NA", "null", "None"}


def f(x):
    """欠測に強い float 変換。空・NAN・-9999 は None。"""
    if x is None:
        return None
    s = str(x).strip()
    if s in NA:
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    if v <= -9990.0:
        return None
    return v


def haversine_km(a, b):
    (la1, lo1), (la2, lo2) = a, b
    r = 6371.0088
    p1, p2 = math.radians(la1), math.radians(la2)
    dp = p2 - p1
    dl = math.radians(lo2 - lo1)
    h = math.sin(dp / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dl / 2) ** 2
    return 2 * r * math.asin(math.sqrt(h))


def spearman(xs, ys):
    """順位相関。tie は平均順位。n<3 なら None。"""
    n = len(xs)
    if n < 3:
        return None, n

    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = rank(xs), rank(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx <= 0 or syy <= 0:
        return None, n
    return sxy / math.sqrt(sxx * syy), n


def season_of(d: date):
    m = d.month
    if m in (3, 4, 5):
        return "spring"
    if m in (9, 10, 11):
        return "autumn"
    if m in (6, 7, 8):
        return "summer"
    return "winter"


# ---------------------------------------------------------------- §1 身元
def sec1_identity():
    print("=" * 78)
    print("§1 身元と座標 — このスタックは誰の上で切られているか")
    print("=" * 78)
    coords = {}
    for root, _dirs, files in os.walk(DATASET):
        for fn in files:
            if fn.endswith("request.json"):
                j = json.load(open(os.path.join(root, fn)))
                for c in j.get("params", {}).get("coordinates", []) or []:
                    coords.setdefault(c["id"], set()).add((c["latitude"], c["longitude"]))
                print(f"  request: {fn}")
                print(f"    layers: {[l.get('product') for l in j['params'].get('layers', [])][:3]}"
                      f" ... ({len(j['params'].get('layers', []))} 層)")
                print(f"    dates : {j['params'].get('dates')}")
    print()
    for k, v in sorted(coords.items()):
        print(f"  抽出点 {k}: {sorted(v)}")
        if len(v) > 1:
            print("    ★ 同じ ID に複数の座標がある（製品ごとに違う）——要注意")

    # 登録簿側の座標（European Fluxes Database Cluster の Site details.html から手写し）
    reg = None
    hp = os.path.join(SITEINFO, "Site details.html")
    if os.path.exists(hp):
        t = re.sub(r"<[^>]+>", "\n", open(hp, errors="replace").read())
        m = re.search(r"Site code:\s*\n\s*(\S+)", t)
        code = m.group(1) if m else "?"
        m = re.search(r"([0-9]+\.[0-9]+)\s*\(lat\)\s*/\s*(-?[0-9]+\.[0-9]+)\s*\(long\)", t)
        if m:
            reg = (float(m.group(1)), float(m.group(2)))
        m = re.search(r"IGBP:\s*\n\s*(\S+)", t)
        igbp = m.group(1) if m else "?"
        m = re.search(r"Mean Annual Precipitation:\s*\n\s*(\S+)", t)
        mapmm = m.group(1) if m else "?"
        print()
        print(f"  登録簿（European Fluxes DB, Site details.html）: code={code} "
              f"coords={reg} IGBP={igbp} MAP={mapmm} mm")
    if reg and "Oran" in coords:
        for c in sorted(coords["Oran"]):
            print(f"  ★ 抽出点 Oran {c} ↔ 登録簿 {reg} = {haversine_km(c, reg)*1000:.0f} m")
    if "Oran" in coords and "TzM" in coords:
        o, z = sorted(coords["Oran"])[0], sorted(coords["TzM"])[0]
        print(f"  ★ Oran ↔ TzM = {haversine_km(o, z):.1f} km（＝別ピクセル・別地点）")
    print("\n  注: TzM は request.json の category が 'ReferenceSite' とだけ書かれている。"
          "\n      TzM にタワーがあるかは、この機械の手元のファイルからは分からない。")
    return coords, reg


def load_master_rain():
    """日次マスタの Rain_mm。**タワー 30 分の P は全期間ゼロで使えない**（欠陥 #62）。"""
    rain = {}
    if not os.path.exists(MASTER_CSV):
        return rain
    with open(MASTER_CSV, newline="", errors="replace") as fh:
        for row in csv.DictReader(fh):
            try:
                dd, mm, yy = row["Date"].split("/")
                d = date(int(yy), int(mm), int(dd))
            except (ValueError, KeyError):
                continue
            v = f(row.get("Rain_mm"))
            if v is not None:
                rain[d] = v
    return rain


# ---------------------------------------------------------------- §2 タワー
def sec2_tower():
    print()
    print("=" * 78)
    print("§2 タワー側の在庫 — 旗36/旗107 が要る量が何日あるか")
    print("=" * 78)
    if not os.path.exists(TOWER_CSV):
        print(f"  !! 見つからない: {TOWER_CSV}")
        return None
    need = {"gH": "H", "gLE": "LE", "Rg": "SW_IN", "th": "SWC_1_1_1",
            "Ts": "TS_1_1_1", "P": "P", "Ta": "TA_1_1_1", "VPD": "VPD",
            "Rn": "NETRAD"}
    halfhour = defaultdict(lambda: defaultdict(list))  # day -> key -> [vals]
    nrow = 0
    bad_ts = 0
    qc_ok = defaultdict(int)
    qc_tot = defaultdict(int)
    with open(TOWER_CSV, newline="", errors="replace") as fh:
        rd = csv.DictReader(fh)
        cols = rd.fieldnames or []
        missing = [k for k, c in need.items() if c not in cols]
        print(f"  列数 {len(cols)} / 要る量のうち欠けている: {missing if missing else 'なし'}")
        for row in rd:
            nrow += 1
            tc = (row.get("TIMECOD") or "").strip()
            if len(tc) != 12 or not tc.isdigit():
                bad_ts += 1
                continue
            d = date(int(tc[:4]), int(tc[4:6]), int(tc[6:8]))
            for k, c in need.items():
                v = f(row.get(c))
                if v is not None:
                    halfhour[d][k].append(v)
            for k, qc in (("gH", "H_QC"), ("gLE", "LE_QC")):
                q = f(row.get(qc))
                if q is not None:
                    qc_tot[k] += 1
                    if q <= 1:
                        qc_ok[k] += 1
    print(f"  30 分行 {nrow} 本 / 時刻が読めない {bad_ts} 本 / 日数 {len(halfhour)}")
    for k in ("gH", "gLE"):
        if qc_tot[k]:
            print(f"  QC: {k} の QC<=1 は {qc_ok[k]}/{qc_tot[k]} = "
                  f"{100*qc_ok[k]/qc_tot[k]:.1f}%")

    # 日集計。被覆規則を先に決める: 24/48 以上の半時間が有効な日だけ採る。
    COV = 24
    daily = {}
    for d, kv in halfhour.items():
        rec = {}
        for k, vals in kv.items():
            if k == "P":
                rec[k] = sum(vals) if len(vals) >= COV else None
            else:
                rec[k] = (sum(vals) / len(vals)) if len(vals) >= COV else None
        daily[d] = rec
    core = ("th", "gH", "gLE", "Rg")
    print(f"\n  日集計の被覆規則: 48 本中 {COV} 本以上が有効な日だけ採る")
    print(f"  {'年':>6} {'季':>7} {'全日':>5} {'th/gH/gLE/Rg 揃う':>18} {'下限60':>7}")
    ok_cells = 0
    tot_cells = 0
    years = sorted({d.year for d in daily})
    for y in years:
        for s in ("spring", "autumn"):
            ds = [d for d in daily if d.year == y and season_of(d) == s]
            full = [d for d in ds if all(daily[d].get(k) is not None for k in core)]
            tot_cells += 1
            hit = len(full) >= MIN_DAYS_PER_CELL
            ok_cells += 1 if hit else 0
            print(f"  {y:>6} {s:>7} {len(ds):>5} {len(full):>18} {'○' if hit else '×':>7}")
    print(f"\n  ★ 春秋のセル {ok_cells}/{tot_cells} が下限 {MIN_DAYS_PER_CELL} 日を満たす"
          f"（暦年は {len(years)} 年・下限 {MIN_YEARS} 年）")
    # ---- 雨の源を検める（★ここで一度、道具が黙って通った）
    print("\n  雨の源の点検（**値を使う前に分布を見る**・旗97 の作法）:")
    pv = [v for kv in halfhour.values() for v in kv.get("P", [])]
    print(f"    30 分 P: 有効 {len(pv)} 本 / 合計 {sum(pv):.1f} mm / 最大 {max(pv) if pv else 0:.2f} mm")
    if pv and sum(pv) == 0.0:
        print("    ★ P は全期間ゼロ。**雨量計の列が死んでいる**（欠測ではなく 0 が入っている）。"
              "\n      ＝この列で旗107 の層（`直後`/`遠い`）は作れない。日次マスタの Rain_mm に移る。")
    rain, rain_src = {}, "30 分 P"
    if not pv or sum(pv) == 0.0:
        rain_src = "日次マスタ Rain_mm"
        rain = load_master_rain()
        ds = sorted(rain)
        print(f"    日次マスタ Rain_mm: {len(rain)} 日 {ds[0]}..{ds[-1]} / "
              f"合計 {sum(rain.values()):.1f} mm")
        for y in years:
            n = sum(1 for d in rain if d.year == y)
            print(f"      {y}: {n} 日 / {sum(v for d, v in rain.items() if d.year == y):.1f} mm")
    else:
        rain = {d: (daily[d].get("P") or 0.0) for d in daily}

    # 層（旗107 の `遠い`／`直後`）が作れるかの粗い見積り: 雨からの経過日数
    print(f"\n  層の見積り（雨の源＝{rain_src}・しきい値 3 mm・**旗133 の感度は今回やらない**）")
    print("    ※ 雨の記録が無い日は、経過日数が決められないので層に入れない")
    last_rain = {}
    prev = None
    for d in sorted(set(daily) | set(rain)):
        if rain.get(d) is None:
            last_rain[d] = None
            continue
        if rain[d] >= 3.0:
            prev = d
        last_rain[d] = (d - prev).days if prev else None
    for lab, lo, hi in (("直後(0-2日)", 0, 2), ("遠い(>=7日)", 7, 10**6)):
        print(f"    層 {lab}:")
        for y in years:
            for s in ("spring", "autumn"):
                n = sum(1 for d in daily
                        if d.year == y and season_of(d) == s
                        and last_rain.get(d) is not None and lo <= last_rain[d] <= hi
                        and all(daily[d].get(k) is not None for k in core))
                print(f"      {y} {s:>7}: {n:>4} 日 "
                      f"{'← 下限60 未満' if n < MIN_DAYS_PER_CELL else ''}")
    return daily


# ---------------------------------------------------------------- §3 前提
def sec3_premise(daily):
    print()
    print("=" * 78)
    print("§3 前提の穴 — この地点は A-3 の 4 クラスタと同じ土俵か")
    print("=" * 78)
    dp = os.path.join(SITEINFO, "Oran_Site_Info.docx")
    if os.path.exists(dp):
        xml = zipfile.ZipFile(dp).read("word/document.xml").decode("utf8", errors="replace")
        txt = " ".join(x.strip() for x in re.sub(r"<[^>]+>", "\n", xml).split("\n") if x.strip())
        print("  Oran_Site_Info.docx の逐語（先頭 700 字）:")
        print("   ", txt[:700].replace(" | ", " "))
        for w in ("irrigat", "Irrigat", "rainfed", "rain-fed"):
            print(f"    語 '{w}' の出現: {txt.count(w)} 回")
    if not daily:
        return
    # 灌漑の指標: 雨ゼロなのに日平均 SWC が跳ねた日を数える（**指標であって証拠ではない**）
    # ★ 雨の源は日次マスタ Rain_mm。タワー 30 分の P は全期間ゼロ（欠陥 #62）。
    rain = load_master_rain()
    ds = sorted(d for d in daily if daily[d].get("th") is not None)
    jumps, judged = [], 0
    for a, b in zip(ds, ds[1:]):
        if (b - a).days != 1:
            continue
        if rain.get(a) is None or rain.get(b) is None:
            continue  # 雨の記録が無い日は判定に使わない
        judged += 1
        dth = daily[b]["th"] - daily[a]["th"]
        if dth >= 1.0 and rain[a] + rain[b] < 0.2:
            jumps.append((b, round(dth, 2)))
    print(f"\n  雨ほぼゼロ(<0.2mm)で日平均 SWC が +1.0%pt 以上跳ねた日: "
          f"{len(jumps)} 日 / 判定できた {judged} 日対")
    bym = defaultdict(int)
    for d, _ in jumps:
        bym[d.month] += 1
    print("    月別:", dict(sorted(bym.items())))
    print("    ★ これは灌漑の**指標**であって証拠ではない（雨量計の欠測・露・"
          "センサ較正でも起きる）。夏(6-8月)に集中していれば灌漑を疑う。")
    print("    ★ 判定できた日対は日次マスタの被覆（2018-01-01..2020-06-25）に縛られる。")


# ---------------------------------------------------------------- §4 衛星
def load_smap():
    per = defaultdict(lambda: defaultdict(list))  # site -> day -> [(surf, root)]
    with open(SMAP_CSV, newline="", errors="replace") as fh:
        for row in csv.DictReader(fh):
            try:
                d = datetime.strptime(row["date_time"].strip(), "%Y-%m-%d %H:%M").date()
            except ValueError:
                continue
            s, r = f(row.get("sm_surface")), f(row.get("sm_rootzone"))
            if s is not None:
                per[row["name"].strip()][d].append((s, r))
    out = {}
    for site, dd in per.items():
        out[site] = {d: (sum(x[0] for x in v) / len(v), len(v)) for d, v in dd.items()}
    return out


def sec4_satellite(smap):
    print()
    print("=" * 78)
    print("§4 衛星側の在庫 — 何が・どの刻みで・タワー年(2018-2020)と重なるか")
    print("=" * 78)
    rows = []

    for site, dd in sorted(smap.items()):
        ds = sorted(dd)
        ov = sum(1 for d in ds if 2018 <= d.year <= 2020)
        rows.append(("SMAP L4 sm_surface/rootzone", site, "3 時間→日",
                     f"{ds[0]}..{ds[-1]}", len(ds), ov))

    # AppEEARS 結果 csv（名前つき）
    app = [("MOD16A2_ET_2018_2024_Oran_TzM", "MOD16A2-2018-2024-Oran-TzM-MOD16A2-061-results.csv",
            "MOD16A2 ET/LE", "8 日", "MOD16A2_061_ET_500m"),
           ("MOD13Q1_NDVI_EVI", "MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv",
            "MOD13Q1 NDVI/EVI", "16 日", "MOD13Q1_061__250m_16_days_NDVI"),
           ("MYD21A1D_LST", "MYD21A1D-LST-MYD21A1D-061-results.csv",
            "MYD21A1D LST", "日", "MYD21A1D_061_LST_1KM")]
    for d, fn, label, cad, vcol in app:
        p = os.path.join(DATASET, d, fn)
        if not os.path.exists(p):
            rows.append((label, "?", cad, "ファイル無し", 0, 0))
            continue
        per = defaultdict(list)
        with open(p, newline="", errors="replace") as fh:
            for row in csv.DictReader(fh):
                try:
                    dt = datetime.strptime(row["Date"].strip(), "%Y-%m-%d").date()
                except ValueError:
                    continue
                v = f(row.get(vcol))
                per[row["ID"].strip()].append((dt, v))
        for site, vs in sorted(per.items()):
            valid = [(dt, v) for dt, v in vs if v is not None and v != 0.0]
            ds = sorted(dt for dt, _ in valid)
            ov = sum(1 for dt in ds if 2018 <= dt.year <= 2020)
            rows.append((label, site, cad,
                         f"{ds[0]}..{ds[-1]}" if ds else "有効値なし", len(ds), ov))

    # 名前つきの GEE 系
    for fn, label, cad, dcol, ncol, vcol in (
            ("OranTzM_S2_NDVI.csv", "Sentinel-2 バンド(B8/B4)", "可変", "date", "name", "B8"),
            ("PML_V2_TimeSeries_Oran_TzM.csv", "PML_V2 ET/Ec/Es/Ei", "8 日", "date", "name", "ET")):
        p = os.path.join(DATASET, fn)
        if not os.path.exists(p):
            continue
        per = defaultdict(list)
        with open(p, newline="", errors="replace") as fh:
            for row in csv.DictReader(fh):
                try:
                    dt = datetime.strptime(row[dcol].strip(), "%Y-%m-%d").date()
                except ValueError:
                    continue
                if f(row.get(vcol)) is not None:
                    per[row[ncol].strip()].append(dt)
        for site, ds in sorted(per.items()):
            ds = sorted(ds)
            ov = sum(1 for dt in ds if 2018 <= dt.year <= 2020)
            rows.append((label, site, cad, f"{ds[0]}..{ds[-1]}", len(ds), ov))

    print(f"  {'製品':<28} {'サイト':<6} {'刻み':<7} {'期間':<24} {'有効':>5} {'18-20':>6}")
    for r in rows:
        print(f"  {r[0]:<28} {r[1]:<6} {r[2]:<7} {r[3]:<24} {r[4]:>5} {r[5]:>6}")

    # 名前の無い Fast_ 系
    print("\n  ★ 名前の無いファイル（`system:index` が 0/1 だけ・サイト列が無い）:")
    for fn in sorted(os.listdir(DATASET)):
        if fn.startswith("Fast_") and fn.endswith(".csv"):
            with open(os.path.join(DATASET, fn), newline="", errors="replace") as fh:
                rr = list(csv.reader(fh))
            hdr = rr[0]
            dates = sorted({h[:10] for h in hdr[1:] if re.match(r"\d{4}-\d\d-\d\d", h)})
            nonempty = [sum(1 for x in row[1:] if f(x) is not None) for row in rr[1:]]
            print(f"    {fn:<32} 行 {len(rr)-1} / 列 {len(hdr)-1} / "
                  f"日付 {dates[0] if dates else '?'}..{dates[-1] if dates else '?'} / "
                  f"行ごとの有効数 {nonempty}")
    print("    → どちらの行が Oran かは、このファイル単体では決まらない（§5 C3 で同定を試みる）")


# ---------------------------------------------------------------- §5 門①
def sec5_gates(daily, smap):
    print()
    print("=" * 78)
    print("§5 門①（対照）— 冒頭で宣言した C1/C2/C3 をそのまま実行する")
    print("=" * 78)
    if daily and "Oran" in smap:
        pairs = []
        for d in sorted(daily):
            th = daily[d].get("th")
            sm = smap["Oran"].get(d)
            if th is not None and sm is not None and 2018 <= d.year <= 2020:
                pairs.append((d, th, sm[0]))
        xs = [p[1] for p in pairs]
        ys = [p[2] for p in pairs]
        rho, n = spearman(xs, ys)
        print(f"  C1 陽性対照: タワー SWC_1_1_1 vs SMAP sm_surface（Oran, 2018-2020）")
        print(f"     n={n} 日  Spearman ρ={rho if rho is None else round(rho, 3)}"
              f"   宣言した期待は ρ>+0.5 → "
              f"{'○ 合格' if (rho or -9) > 0.5 else '× 不合格（読み取りを疑う）'}")

        # C2-a: 1 階差分
        dx, dy = [], []
        for (d1, t1, s1), (d2, t2, s2) in zip(pairs, pairs[1:]):
            if (d2 - d1).days == 1:
                dx.append(t2 - t1)
                dy.append(s2 - s1)
        r2, n2 = spearman(dx, dy)
        print(f"  C2a 差分（季節周期を抜く）: n={n2}  ρ={r2 if r2 is None else round(r2, 3)}")
        # C2-b: 衛星を +45 日ずらす
        thmap = {p[0]: p[1] for p in pairs}
        sx, sy = [], []
        for d, th, _ in pairs:
            sm = smap["Oran"].get(d + timedelta(days=45))
            if sm is not None:
                sx.append(th)
                sy.append(sm[0])
        r3, n3 = spearman(sx, sy)
        print(f"  C2b +45 日ずらし: n={n3}  ρ={r3 if r3 is None else round(r3, 3)}")
        if rho is not None and r3 is not None:
            print(f"     → ずらしで ρ が {rho:.3f} → {r3:.3f}（差 {rho - r3:+.3f}）")
        if r2 is not None:
            print(f"     → 差分でも ρ={r2:.3f}"
                  f"{'（日々の事象で一致している）' if r2 > 0.2 else '（季節周期しか合っていない可能性）'}")
        del thmap
    else:
        print("  C1/C2: タワーか SMAP が読めず実行できない")

    # C3: Fast_ の行順の同定
    print("\n  C3 Fast_*.csv の行順の同定（Fast_MOD16 ↔ 名前つき AppEEARS MOD16）")
    fp = os.path.join(DATASET, "Fast_OranTzM_MOD16.csv")
    ap = os.path.join(DATASET, "MOD16A2_ET_2018_2024_Oran_TzM",
                      "MOD16A2-2018-2024-Oran-TzM-MOD16A2-061-results.csv")
    if os.path.exists(fp) and os.path.exists(ap):
        with open(fp, newline="", errors="replace") as fh:
            rr = list(csv.reader(fh))
        hdr = rr[0]
        fast = {}
        for row in rr[1:]:
            idx = row[0]
            series = {}
            for h, v in zip(hdr[1:], row[1:]):
                m = re.match(r"(\d{4}-\d\d-\d\d)_ET$", h)
                if m and f(v) is not None:
                    series[datetime.strptime(m.group(1), "%Y-%m-%d").date()] = f(v)
            fast[idx] = series
        named = defaultdict(dict)
        with open(ap, newline="", errors="replace") as fh:
            for row in csv.DictReader(fh):
                try:
                    dt = datetime.strptime(row["Date"].strip(), "%Y-%m-%d").date()
                except ValueError:
                    continue
                v = f(row.get("MOD16A2_061_ET_500m"))
                if v is not None:
                    named[row["ID"].strip()][dt] = v
        res = {}
        for fi, fs in sorted(fast.items()):
            for site, ns in sorted(named.items()):
                com = sorted(set(fs) & set(ns))
                if len(com) >= 10:
                    r, n = spearman([fs[d] for d in com], [ns[d] for d in com])
                    res[(fi, site)] = (r, n)
                    print(f"     行 {fi} ↔ {site}: n={n} ρ={r if r is None else round(r, 3)}")
        if len(res) == 4:
            a = res.get(("0", "Oran"), (None,))[0], res.get(("0", "TzM"), (None,))[0]
            b = res.get(("1", "Oran"), (None,))[0], res.get(("1", "TzM"), (None,))[0]
            if None not in a and None not in b:
                m0 = "Oran" if a[0] > a[1] else "TzM"
                m1 = "Oran" if b[0] > b[1] else "TzM"
                gap = min(abs(a[0] - a[1]), abs(b[0] - b[1]))
                if m0 != m1 and gap >= 0.05:
                    print(f"     ★ 同定: 行0={m0} / 行1={m1}（差の小さい方でも {gap:.3f}）")
                else:
                    print(f"     ★ 同定できない（行0={m0}, 行1={m1}, 差 {gap:.3f}）"
                          f"— 推測で割り当てない")
    else:
        print("     ファイルが無く実行できない")


def main():
    print("旗137 下調べ: Oran/TzM 衛星スタック × ES-FcO タワー")
    print("実行:", datetime.now().isoformat(timespec="seconds"))
    print("**この道具は A-3 の量を計算しない。実行可能性だけを測る。**\n")
    sec1_identity()
    daily = sec2_tower()
    sec3_premise(daily)
    smap = load_smap()
    sec4_satellite(smap)
    sec5_gates(daily, smap)
    print("\n" + "=" * 78)
    print("下調べはここまで。**θ→γH・θ→γLE の Δ は計算していない**"
          "（事前登録を書いてから）。")
    print("=" * 78)


if __name__ == "__main__":
    sys.exit(main())
