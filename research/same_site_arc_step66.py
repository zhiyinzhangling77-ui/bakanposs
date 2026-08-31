"""旗66：**同一地点でタワーとチャンバーを同じ物差しで並べる**＝穴②を閉じる。

旗64/65 で分かったこと：**JP-Fhk ↔ `d20200328_UEYAMA_HOKUROKU` が 0.00 km**で、
チャンバー側にデータがあり **★短メモリ**（R²=0.90, ACF1=+0.75, e-fold=7日）。
JP-Fhk は**閉合最良（EBR 0.85）でSIFも取得済み**＝**タワー・衛星・チャンバーが同一地点に揃う唯一の組**。

これまで穴②（同一地点で3観測系を突き合わせていない）は本研究**最大の弱点**として残っていた。
本ツールは、**同じ検出器（旗53/54 の較正済み・非線形基底）・同じ判定基準**で両側を測り、
**できれば同じ期間**（チャンバーの観測期間にタワーを合わせる）で並べる。

  ・両側で「メモリが在る」と出る → **同一地点で分割の有無に依らず再現**＝弧が閉じる。
  ・片側だけ → **分割（タワー側）か点測定（チャンバー側）に固有**の可能性＝そう記す。

対照として **JP-Tef ↔ `UEYAMA_TESHIO`（0.01km）**も並べる（チャンバー側は季節メモリ e-fold=30日）。
**同じ場所で両側とも「短メモリなし」なら、それも整合の証拠**になる。

    python research/same_site_arc_step66.py --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore, _acf_gap, _efold_gap
from memory_attribution_flex_step54 import flex_basis, _fit, ACF_THR, EFOLD_MAX

# 旗64/65/67 で確定した同一地点の対（距離は旗51/67 の値）。
# 旗67 で**登録済み10サイトでなく手元の83サイト全部**と突き合わせた結果、対が 2→7 に増えた。
PAIRS = [("JP-Fhk", "d20200328_UEYAMA_HOKUROKU", 0.00),   # 済（旗66・一致）
         ("JP-Tef", "d20200328_UEYAMA_TESHIO", 0.01),     # 済（旗66・不一致）
         ("JP-Yms", "d20200212_ATAKA", 0.03),             # 済（旗66）
         ("JP-Yms", "d20200328_UEYAMA_YAMASHIRO", 0.69),  # 済（旗66）
         # --- 旗79/80：AmeriFlux 登録で加わった対（距離は BADM 座標から実測）---
         # **これで対が日本だけでなくなる**＝旗43 が名指しした擬似反復の懸念に直接効く。
         ("US-SSH", "d20200212_KAYE_UNE", 0.05),
         ("US-SSH", "d20200212_KAYE_UNW", 0.06),
         ("US-SSH", "d20200212_KAYE_USW", 0.07),
         ("US-SSH", "d20200212_KAYE_USE", 0.07),
         ("US-SSH", "d20200212_KAYE_LNE", 0.10),
         ("US-SSH", "d20200212_KAYE_LNW", 0.12),
         # LSE/LSW は旗79 の表示が 6 件で切れていたが、旗78 の座標（40.66,-77.90）から
         # **同一地点のはず**。**距離は未実測**なので、そう明記して並べる。
         ("US-SSH", "d20200212_KAYE_LSE", float("nan")),
         ("US-SSH", "d20200212_KAYE_LSW", float("nan")),
         ("US-Ha1", "d20190415_VARNER", 0.27),
         ("US-Ha1", "d20190504_SAVAGE_hf006-03", 0.27),
         ("US-Ha1", "d20190504_SAVAGE_hf006-05", 0.27),
         ("US-MMS", "d20190424_ZHANG_maple", 0.00),
         ("US-MMS", "d20190424_ZHANG_oak", 0.00),
         # Turkey Point：**0.83km しか離れていない2本のタワー**なので、
         # **距離だけでは対を決められない**。植林年（1939/1974）で対応させた——
         # そして**距離 0.00km の組み合わせが名前とも一致した**＝独立に同じ答えが出た。
         ("CA-TP4", "d20200417_ARAIN_TP39", 0.00),
         ("CA-TP3", "d20200417_ARAIN_TP74", 0.00),
         ("CA-TPD", "d20200417_ARAIN_TPD", 0.00),
         # 乾燥地（森林ではない＝旗57 からメモリ自体が出ない公算が高い。対照として並べる）
         ("US-Wkg", "d20190617_SCOTT_WKG", 0.01),
         ("US-Whs", "d20190617_SCOTT_WHS", 0.02),
         ("US-SRM", "d20190617_SCOTT_SRM", 3.60),
         # --- 旗86：**★で選ばずに取得した対**（ここが本命）------------------
         # 上の 23 組は**旗78 の★順で取得サイトを選んだ**結果を含む＝**私の選択バイアス入り**。
         # 旗85 は「10km 以内にタワーがある」という**座標だけの条件**で未取得の組を名指しし、
         # そのうち**取得できた 17 本**がここ。**★は一切見ていない**。
         # 距離は全て旗79（欠陥21/22 修正版）の BADM 座標実測。
         # **チャンバーは最近傍のタワーに 1 度だけ割り当てる**（重複計上を避ける）——
         # 例：SZUTU_TWITCHELL は US-Bi2 から 9.79km だが US-Tw3 からは 0.08km なので Tw3 に付ける。
         ("BR-Sa3", "d20190527_GOULDEN", 0.14),            # **熱帯**（アマゾン）
         ("CL-SDF", "d20200419_PEREZ-QUEZADA", 0.00),      # **南半球**（チリ）
         ("US-Uaf", "d20200328_UEYAMA_FAIRBANKS", 0.00),   # **高緯度・永久凍土**
         ("CA-Ca1", "d20200108_JASSAL", 0.55),
         ("US-Me6", "d20190520_RUEHR", 0.00),
         ("US-NC2", "d20200220_GAVAZZI", 0.00),
         ("US-NC4", "d20190918_MIAO", 0.00),
         ("US-WCr", "d20190430_DESAI", 0.00),
         ("US-xSE", "d20190526_PENNINGTON", 1.82),
         ("US-Ton", "d20200214_SZUTU_TONZI", 0.00),
         ("US-Ton", "d20191017_BALDOCCHI", 0.08),
         ("US-Var", "d20200214_SZUTU_VAIRA", 0.01),
         ("US-Bi1", "d20200214_SZUTU_BOULDINA", 0.01),
         ("US-Bi2", "d20200214_SZUTU_BOULDINC", 0.01),
         ("US-Tw3", "d20200214_SZUTU_TWITCHELL", 0.08),
         # **感度確認**：SCOTT_SRM に対して US-SRM(3.60km) より US-SRS(2.20km) が近い。
         # だがチャンバー名 SRM は Santa Rita **Mesquite** ＝ US-SRM を指す。
         # **距離と名前が食い違う**ので両方並べ、**集計では US-SRM 側だけを数える**
         # （同じチャンバーを二度数えないため）。読み方は旗86 の出力で明示する。
         ("US-SRS", "d20190617_SCOTT_SRM", 2.20),
         # **旗79 修正版（欠陥21）で座標が取れた 2 本**。
         # 相手のチャンバーは**私の見込みと違った**——US-Ho1 は DAVIDSON ではなく **SIHI**、
         # US-UMB は **MATHES**。**推測で対を書かず、座標で確かめてから書いた**のが効いた。
         ("US-Ho1", "d20190610_SIHI_H1", 0.01),
         ("US-Ho1", "d20190610_SIHI_H2", 0.01),
         # **同じ座標だが湿地**＝旗51 の「近さは十分条件でない」。**別生態系として読む**。
         ("US-Ho1", "d20190610_SIHI_H2wetland", 0.01),
         # **同一著者の 2 データセットが同じ座標**（d20200221 と d20200224）。
         # **独立ではない可能性が高い**＝旗43 の擬似反復の対象として、
         # クラスタ縮約では **1 つに数える**こと。
         ("US-UMB", "d20200221_MATHES", 0.00),
         ("US-UMB", "d20200224_MATHES", 0.00),
         # **距離は個別 BIF による**。multi-site の LEGACY 版とは食い違いがあり、
         # US-NC4 で 1.11km、US-Me6 で 0.27km、US-WCr で 0.03km ずれる（旗79 が報告する）。
         # **どちらが正しいかは決めていない。10km の判定を跨ぐ差は無い。**
         # **US-HB1 は最近傍チャンバーが 358km ＝対を作れない**（登録もしていない）。
         # インドネシア泥炭3組は**タワー側に土壌温度・水分が無い**（旗68 で本体HH＋補助ファイルを
         # 全部走査して確認）＝同一地点比較に使えないので**外した**。
         # ("ID-PaB", "d20200109_HIRANO_PDB", 0.00),
         # ("ID-PaD", "d20200109_HIRANO_PDF", 0.00),
         # ("ID-Pag", "d20200109_HIRANO_PUF", 0.44),
         ]
# **同じチャンバーを二度数えないための除外集合**。
# SCOTT_SRM は US-SRM（名前が一致・3.60km）と US-SRS（距離が近い・2.20km）の両方に対で
# 並べてあるが、**集計では US-SRM 側だけを数える**。US-SRS 側は
# 「**近さで選ぶか名前で選ぶかで答えが変わるか**」を見る感度確認である。
SENSITIVITY_ONLY = {("US-SRS", "d20190617_SCOTT_SRM")}

# 旗38 の実測（記録済み）：季節制御後の SIF 偏相関と、記憶ACFが落ちたか
SIF_NOTE = {"JP-Fhk": "偏相関 +0.03・記憶ACFは落ちず",
            "JP-Tef": "偏相関 +0.09・記憶ACFは落ちず"}   # 他サイトはSIF未取得


def memory_from_daily(daily, ycol, tcol, wcol):
    """旗53/54 と同じ非線形基底で残差メモリを測る（タワー・チャンバー共通の物差し）。"""
    y = np.log(daily[ycol].where(daily[ycol] > 0)).to_numpy()
    T = daily[tcol].to_numpy()
    W = daily[wcol].to_numpy() if wcol and wcol in daily else None
    res, r2 = _fit(y, flex_basis(T, W))
    if res is None or not np.isfinite(r2):
        return None
    return {"r2": r2, "acf1": _acf_gap(res, 1), "efold": _efold_gap(res),
            "n": int(np.isfinite(res).sum())}


def tower_daily(site, span=None):
    from japanflux_pn.config import AnalysisConfig
    from japanflux_pn.sites import get_site
    from japanflux_pn.preprocess import load_raw_all
    cfg = AnalysisConfig()
    raw = load_raw_all(get_site(site), cfg)
    d = raw[["GER", "Ts", "th"]].copy()
    if span:
        d = d.loc[(d.index >= span[0]) & (d.index <= span[1])]
    if d.empty:
        return None
    daily = d.groupby(d.index.normalize()).mean()
    return daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D"))


def chamber_daily(path, span=None):
    df, st, sm = load_cosore(path, None)
    if "Tsoil" not in df:
        return None, (st, sm)
    cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
    d = df[cols].copy()
    if span:
        d = d.loc[(d.index >= span[0]) & (d.index <= span[1])]
    if d.empty:
        return None, (st, sm)
    daily = d.groupby(d.index.normalize()).mean()
    return daily.reindex(pd.date_range(daily.index.min(), daily.index.max(), freq="D")), (st, sm)


def verdict(m):
    if m is None:
        return "判定不能"
    if not np.isfinite(m["r2"]) or m["r2"] < 0.3:
        return f"駆動弱(R²={m['r2']:.2f})＝判定不能"
    if not np.isfinite(m["acf1"]) or not np.isfinite(m["efold"]):
        return "推定不能"
    if m["acf1"] >= ACF_THR and m["efold"] <= EFOLD_MAX:
        return "★短メモリ"
    if m["acf1"] >= ACF_THR:
        return f"季節メモリ(e-fold={m['efold']:.0f}日)"
    return f"メモリ弱(ACF1={m['acf1']:.2f})"


def main():
    p = argparse.ArgumentParser(description="同一地点でタワーとチャンバーを並べる")
    p.add_argument("--cosore-dir", required=True)
    a = p.parse_args()
    root = Path(a.cosore_dir)

    print("=== 旗66：同一地点でタワー×チャンバーを同じ物差しで並べる（穴②）===")
    print(f"  検出器は旗53/54 の較正済み（非線形基底・R²≥0.3・ACF1≥{ACF_THR}・e-fold≤{EFOLD_MAX}日）。")
    print("  タワー側は GER（分割派生）、チャンバー側は Rs（分割を通さない直接測定）。\n")

    for site, ds, km in PAIRS:
        f = root / "datasets" / f"data_{ds}.csv"
        print(f"  ━━ {site} ↔ {ds}（{km:.2f} km）━━")
        if not f.exists():
            print("    チャンバーのデータファイルが無い\n"); continue
        ch_all, (st, sm) = chamber_daily(f)
        if ch_all is None:
            print("    チャンバー側に土壌温度が無い\n"); continue
        span = (ch_all.index.min(), ch_all.index.max())
        print(f"    チャンバー観測期間：{span[0]:%Y-%m}〜{span[1]:%Y-%m}（列 T={st} SM={sm}）")

        try:
            tw_all = tower_daily(site)
            tw_ov = tower_daily(site, span)
        except Exception as e:
            # **例外の中身も出す**（旗66 第1版は型名だけで、未登録なのか変数欠けなのか分からなかった）
            msg = str(e)
            print(f"    タワー側の読み込み失敗 {type(e).__name__}: {msg[:160]}\n"); continue

        rows = []
        m = memory_from_daily(ch_all, "Rs", "Tsoil", "SM")
        rows.append(("チャンバー Rs（全期間）", m))
        if tw_all is not None:
            rows.append(("タワー GER（全期間）", memory_from_daily(tw_all, "GER", "Ts", "th")))
        if tw_ov is not None and len(tw_ov.dropna(subset=["GER"])) >= 60:
            rows.append(("タワー GER（**同一期間**）", memory_from_daily(tw_ov, "GER", "Ts", "th")))
            rows.append(("チャンバー Rs（同一期間）", m))   # チャンバーは元から同一期間
        else:
            print("    ※タワー側に同一期間のデータが足りず、期間を揃えた比較はできない")

        print(f"    {'系列':<26}{'N':>6}{'R²':>7}{'ACF1':>8}{'e-fold':>8}  判定")
        for lab, mm in rows:
            if mm is None:
                print(f"    {lab:<26}{'—':>6}{'—':>7}{'—':>8}{'—':>8}  判定不能"); continue
            print(f"    {lab:<26}{mm['n']:>6}{mm['r2']:>7.2f}{mm['acf1']:>8.2f}"
                  f"{mm['efold']:>8.0f}  {verdict(mm)}")
        print(f"    衛星SIF（旗38 の記録）：{SIF_NOTE.get(site, '—')}")
        print()

    print("  === 読み方 ===")
    print("  同一地点で**両側とも★短メモリ**＝分割の有無に依らず再現＝**穴②が閉じる**。")
    print("  片側だけ★＝分割（タワー）か点測定（チャンバー）に固有の可能性。")
    print("  両側とも★でない＝その地点にはメモリが無い（それも整合の証拠）。")
    print("  留保：タワーは生態系呼吸（GER, 分割派生）・チャンバーは土壌呼吸（Rs, 点測定）で")
    print("        **測っている対象が違う**。同一地点であっても『同じ量』ではない。")
    print("        チャンバー期間が短い場合、同一期間のタワー側は年数不足で不安定になりうる。")


if __name__ == "__main__":
    main()
