"""旗128：**到達性スクリーン（Unpaywall）を、陽性の対照に当てて較正する。**

## なぜ要るのか

**旗127 は申し送りとして次を書いた**（`NOVELTY_ASSESSMENT.md:626`・`SESSION_STATE.md:1102`）：

  > ⚪ に着手する前に `https://api.unpaywall.org/v2/<DOI>?email=<...>` を叩く。
  > **`is_oa:false` かつ `has_repository_copy:false` なら、機関ミラー探しは無駄である。**

**これは「探索を打ち切ってよい」という規則である。にもかかわらず、この規則は
陽性の対照に一度も当てられていない**——**すなわち「既に本文へ到達できた文献を、
この判定器は通すのか」を確かめていない。**

**旗104（陽性の対照が陽性を出せていなかった＝欠陥 #36）と、
旗126 の作法「合成の対照は、名乗る性質を実際に持っていることを数値で確かめる」の、
文献到達側での同型である。**

## 何をするか

  1. **陽性の対照**＝**本文まで一次到達した実績のある 5 件**に判定器を当てる。
     **判定器が「closed」と言う件が出れば、それは偽陰性であり、規則は成り立たない。**
  2. **本題**＝**残る ⚪ 5 件**に同じ判定器を当てる。
  3. **両方を同じ表に並べる**（**対照を別に走らせて別に読むと、比べ忘れる**）。

## 読み方の注意

- **この判定器は「到達できるか」ではなく「Unpaywall が OA コピーを索引しているか」を測る。**
  **索引していないことは、コピーが存在しないことを意味しない。**
- **`WebFetch` の要約モデルを経由せず、生の JSON を読む**（旗127 までは要約経由だった）。
  **要約モデルが値を取り違える余地を消す。**
"""
from __future__ import annotations

import json
import urllib.error
import urllib.request

from runlog import tee_stdout

# Unpaywall は連絡先の付与を求める。**個人のアドレスは渡さない**——
# 公式ドキュメントが例示に使う窓口アドレスを使う。
EMAIL = "unpaywall@impactstory.org"

# ---------------------------------------------------------------------------
# **陽性の対照**：本文まで一次到達した実績のある文献（旗118〜127）。
#   `reached` は「どの周に、どの経路で本文へ届いたか」。**実績であって推測ではない。**
# ---------------------------------------------------------------------------
CONTROLS = [
    ("10.5194/bg-11-259-2014", "Wang 2014, Biogeosciences 11, 259–268",
     "旗122・Copernicus が直接 OA"),
    ("10.3402/tellusb.v39i1-2.15329", "Dörr & Münnich 1987, Tellus 39B, 114–121",
     "旗127・b.tellusjournals.se → GCS の PDF を Read pages"),
    ("10.1046/j.1365-2486.1998.00128.x", "Davidson, Belk & Boone 1998, GCB 4, 217–227",
     "旗124・harvardforest1.fas.harvard.edu の刊行物索引から PDF 直リンク"),
    ("10.1016/j.agrformet.2018.05.005", "Zhang 2018, AFM 259, 184–195",
     "旗121・Elsevier 本体ではない機関ミラー"),
    ("10.1016/j.agrformet.2006.02.005", "Wen 2006, AFM 137, 166–175",
     "旗123・全文到達"),
]

# ---------------------------------------------------------------------------
# **本題**：残る ⚪（`NOVELTY_ASSESSMENT.md:632` の優先順・旗127 で更新されたもの）。
# ---------------------------------------------------------------------------
OPEN = [
    ("10.1016/0038-0717(77)90058-X", "Bunnell et al. 1977, Soil Biol. Biochem. 9, 33–40",
     "Davidson 1998 の非同定性実演の元モデル"),
    ("10.1029/2000GB001365", "Xu & Qi 2001, GBC 15, 687–696",
     "`Q10 = a − bT + cS_w` の原型"),
    ("10.1016/j.agrformet.2006.08.003", "Gaumont-Guay 2006, AFM 140, 220–235",
     "水分と温度の依存の解釈・UBC Biomet"),
    ("10.1007/s10533-016-0200-1", "Tucker & Reed 2016, Biogeochemistry 128, 155–169",
     "題名が『見かけの負の温度感度は水分の支配』＝A-2 に最も近い"),
    ("10.1111/j.1365-2486.2011.02516.x", "Suseela et al. 2012, GCB 18, 336–348",
     "温度感度への水分効果が季節で変わる"),
]


def query(doi: str) -> dict:
    """Unpaywall を 1 件引く。**落ちても走行を止めず、理由を値として返す。**"""
    url = f"https://api.unpaywall.org/v2/{doi}?email={EMAIL}"
    try:
        with urllib.request.urlopen(url, timeout=30) as fh:
            return json.load(fh)
    except urllib.error.HTTPError as e:
        return {"_error": f"HTTP {e.code}"}
    except Exception as e:                     # 到達不能そのものを結果として残す
        return {"_error": f"{type(e).__name__}: {str(e)[:60]}"}


def verdict(rec: dict) -> str:
    """**判定器の出力を、旗127 の規則の言葉に翻訳する。**"""
    if "_error" in rec:
        return "問い合わせ失敗"
    if rec.get("is_oa"):
        return "通す（OA）"
    if rec.get("has_repository_copy"):
        return "通す（リポジトリ複製あり）"
    return "打ち切れと言う"          # ＝旗127 の規則が「機関ミラー探しは無駄」と結論する枝


def report(title: str, rows: list[tuple[str, str, str]]) -> list[tuple]:
    print(f"\n{'=' * 78}\n{title}\n{'=' * 78}")
    out = []
    for doi, label, note in rows:
        rec = query(doi)
        v = verdict(rec)
        n_loc = len(rec.get("oa_locations") or [])
        print(f"\n  {label}")
        print(f"    doi           : {doi}")
        print(f"    備考          : {note}")
        if "_error" in rec:
            print(f"    **問い合わせ失敗**: {rec['_error']}")
        else:
            print(f"    is_oa         : {rec.get('is_oa')}")
            print(f"    has_repo_copy : {rec.get('has_repository_copy')}")
            print(f"    oa_status     : {rec.get('oa_status')}")
            print(f"    oa_locations  : {n_loc} 件")
            for loc in (rec.get("oa_locations") or []):
                print(f"        - {loc.get('host_type')} / {loc.get('version')} "
                      f"/ {loc.get('url_for_pdf') or loc.get('url')}")
        print(f"    **判定器の出力**: {v}")
        out.append((label, v))
    return out


def main() -> None:
    tee_stdout("step128_oa_screen")

    print("""
旗128：到達性スクリーン（Unpaywall）の較正
  判定器 : is_oa==False かつ has_repository_copy==False なら「機関ミラー探しは無駄」（旗127 の申し送り）
  対照   : **本文まで一次到達した実績のある 5 件**。ここで「打ち切れ」が出たら、それは偽陰性。
  注意   : この判定器が測るのは「Unpaywall が OA コピーを索引しているか」であって
           「到達できるか」ではない。**索引していない ≠ 存在しない。**
""")

    ctrl = report("【門①・陽性の対照】本文へ一次到達した実績のある 5 件", CONTROLS)
    open_ = report("【本題】残る ⚪ 5 件", OPEN)

    fn = [lab for lab, v in ctrl if v == "打ち切れと言う"]
    print(f"\n{'=' * 78}\n総括\n{'=' * 78}")
    print(f"  対照 {len(ctrl)} 件のうち、判定器が **打ち切れと言った件＝偽陰性** : "
          f"**{len(fn)} 件**")
    for lab in fn:
        print(f"      - {lab}")
    if fn:
        print("\n  **＝旗127 の申し送り『closed なら機関ミラー探しは無駄』は成り立たない。**")
        print("  **偽陰性の件は、実際には機関ミラー（Harvard Forest 索引など）で本文に届いている。**")
        print("  **⚪ を Unpaywall の closed だけで閉じてはならない。**")
    else:
        print("\n  対照は全件通った。**この範囲では規則を否定する証拠は出ていない**"
              "（成り立つことの証明ではない）。")
    print(f"\n  本題側で『打ち切れと言う』となった ⚪ : "
          f"**{sum(1 for _, v in open_ if v == '打ち切れと言う')} / {len(open_)} 件**"
          f"（**上の較正により、これは到達不能の根拠にならない**）")


if __name__ == "__main__":
    main()
