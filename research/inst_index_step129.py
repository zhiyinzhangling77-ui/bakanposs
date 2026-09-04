#!/usr/bin/env python3
"""旗129 — 残る ⚪ 5 件に「著者の所属機関の刊行物索引」の手を当てる（欠陥 #52 の解消）。

旗128 が確かめたこと：本文まで一次到達した 5 件のうち 3 件は Unpaywall/OpenAlex で `closed` である。
索引 API に見えない経路（著者の所属機関・研究室・学会の刊行物一覧）が実際に効いた。

## この版が前の草稿（未実行・旗128 の周に書かれた）と違う点

前の草稿は **PDF のファイル名を推測して直に叩いていた**
（`.../boreal_aspen_2006.pdf` など）。**これはこの研究で 3 度出ている失敗型**
（`SESSION_STATE.md`「ファイル名で当てにいく」＝旗64/82/86）で、
作法は **「索引を作ってから引く」** である。**推測した URL が 404 を返しても、
それは「本文が無い」ことを何も意味しない。**

そこで本版は 2 段に分ける：

  段1 **索引を取る** — 所属機関・研究室・リポジトリの *一覧ページ / 検索 API* を叩き、
       ページ内の **全リンクを抽出**する（アンカー文と href の対）。
  段2 **索引の中を照合する** — 著者姓・年・題名の語で絞り、**当たった href だけ**を叩く。

**どの段で落ちたかを分けて記録する**（旗108・欠陥 #40「落ちた理由を一つの籠に入れない」）：
  `INDEX_UNREACHABLE`（索引そのものに届かない）/ `NOT_IN_INDEX`（索引は取れたが載っていない）/
  `LINK_DEAD`（索引に載っていたが本体が落ちる）/ `GOT_PDF`（本体を取れた）。

値は WebFetch の要約モデルを経由させず urllib.request で直接取る（旗128 の申し送り）。

出力: research/logs/step129_<timestamp>.txt（runlog・欠陥 #42）
PDF は research/tmp_pdfs/ に落とす（版管理しない。生成物）。
"""

import hashlib
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from runlog import tee_stdout  # noqa: E402

UA = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/124.0 Safari/537.36"
)

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tmp_pdfs")

# 索引ページ内のリンク照合に使う語。**全部を要求しない**（索引の表記ゆれで落ちるため）。
# `must_any`: どれか 1 つ / `year`: 年（あれば加点）。加点式にして、上位を人が読む。
TARGETS = [
    {
        "tag": "Gaumont-Guay2006",
        "cite": "Gaumont-Guay et al. 2006, Agric. For. Meteorol. 140, 220-235",
        "aff": "UBC Biometeorology / Andy Black lab",
        "must_any": ["gaumont", "gaumont-guay"],
        "year": "2006",
        "indexes": [
            "https://biomet.ubc.ca/publications/",
            "https://ibis.geog.ubc.ca/~achristn/publications.html",
            "https://www.landfood.ubc.ca/andrew-black/",
        ],
    },
    {
        "tag": "Tucker-Reed2016",
        "cite": "Tucker & Reed 2016, Biogeochemistry 128, 155-169",
        "aff": "USGS Southwest Biological Science Center",
        "must_any": ["tucker"],
        "year": "2016",
        "indexes": [
            "https://www.usgs.gov/publications/low-soil-moisture-during-hot-periods-drives-apparent-negative-temperature",
            "https://www.sciencebase.gov/catalog/items?q=Low+soil+moisture+during+hot+periods&format=json&max=20",
            "https://www.usgs.gov/staff-profiles/sasha-reed",
        ],
    },
    {
        "tag": "Xu-Qi2001",
        "cite": "Xu & Qi 2001, Global Biogeochem. Cycles 15, 687-696",
        "aff": "LBNL / UC eScholarship",
        "must_any": ["xu", "ponderosa"],
        "year": "2001",
        "indexes": [
            "https://escholarship.org/search/?q=%22soil-surface%20CO2%20efflux%22%20ponderosa%20Xu%20Qi",
            "https://escholarship.org/api/search?q=ponderosa%20pine%20soil%20CO2%20efflux%20Xu",
            "https://www.osti.gov/api/v1/records?title=Soil-surface%20CO2%20efflux&rows=20",
        ],
    },
    {
        "tag": "Suseela2012",
        "cite": "Suseela et al. 2012, Global Change Biol. 18, 336-348",
        "aff": "MBL / Dukes lab (BACE)",
        "must_any": ["suseela"],
        "year": "2012",
        "indexes": [
            "https://www.mbl.edu/research/publications?search_api_fulltext=Suseela",
            "https://jeffdukes.weebly.com/publications.html",
            "https://docs.lib.purdue.edu/do/search/?q=Suseela&start=0&context=52108",
        ],
    },
    {
        "tag": "Bunnell1977",
        "cite": "Bunnell et al. 1977, Soil Biol. Biochem. 9, 33-40",
        "aff": "UBC (Forestry / IBP Tundra Biome)",
        "must_any": ["bunnell"],
        "year": "1977",
        "indexes": [
            "https://open.library.ubc.ca/search?q=Bunnell%20microbial%20respiration%20soil",
            "https://www.forestry.ubc.ca/faculty/fred-bunnell/",
            "https://scholar.archive.org/search?q=%22Microbial+respiration+and+substrate+weight+loss%22",
        ],
    },
]

A_TAG = re.compile(r"<a\b[^>]*?href=[\"']([^\"'#]+)[\"'][^>]*>(.*?)</a>",
                   re.IGNORECASE | re.DOTALL)
TAG_STRIP = re.compile(r"<[^>]+>")
WS = re.compile(r"\s+")


def fetch(url, timeout=45):
    """1 本叩く。**成否と理由を分けて返す**（例外を握り潰さない）。"""
    req = urllib.request.Request(url, headers={
        "User-Agent": UA,
        "Accept": "text/html,application/xhtml+xml,application/json,application/pdf,*/*",
        "Accept-Language": "en-US,en;q=0.9",
    })
    t0 = time.time()
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return {"ok": True, "status": r.status, "body": r.read(),
                    "ctype": r.headers.get("Content-Type", "?"),
                    "final": r.geturl(), "dt": round(time.time() - t0, 1)}
    except urllib.error.HTTPError as e:
        return {"ok": False, "status": e.code, "why": f"HTTPError {e.reason}",
                "dt": round(time.time() - t0, 1)}
    except Exception as e:
        return {"ok": False, "status": None, "why": f"{type(e).__name__}: {e}",
                "dt": round(time.time() - t0, 1)}


def links_of(body, base):
    """索引ページから (アンカー文, 絶対 URL) を全部拾う。**照合はこの後**。"""
    text = body.decode("utf-8", "replace")
    out = []
    for href, inner in A_TAG.findall(text):
        label = WS.sub(" ", TAG_STRIP.sub(" ", inner)).strip()
        out.append((label, urllib.parse.urljoin(base, href)))
    # JSON API には <a> が無い。URL らしき文字列も拾う（sciencebase / osti 用）
    if not out:
        for m in re.finditer(r'https?://[^\s"\'<>\\]+', text):
            out.append(("", m.group(0)))
    return out


def score(label, url, spec):
    """索引の 1 行が目当てのものかを加点で測る。**閾値で切らず、上位を印字して人が読む。**"""
    hay = (label + " " + urllib.parse.unquote(url)).lower()
    s = 0
    if any(k in hay for k in spec["must_any"]):
        s += 2
    if spec["year"] in hay:
        s += 1
    if hay.endswith(".pdf") or ".pdf" in hay:
        s += 1
    return s


def run_target(spec):
    print(f"\n{'=' * 78}\n### {spec['tag']} — {spec['cite']}\n    所属索引: {spec['aff']}")
    index_ok = False      # 1 本でも HTTP が通ったか
    index_substantive = False  # **中身のある索引が 1 本でもあったか**（欠陥 #53）
    author_in_raw = False      # 著者名が生 HTML に在るのに、リンクの札には無い場合を分ける
    cands = []
    for url in spec["indexes"]:
        r = fetch(url)
        if not r["ok"]:
            print(f"  [索引] ✗ {url}\n         status={r['status']} {r['why']} ({r['dt']}s)")
            continue
        index_ok = True
        ls = links_of(r["body"], r["final"])
        # **欠陥 #53**：HTTP 200 は「中身が取れた」ことを意味しない。
        # JS で描かれるページは 2〜4 KB・リンク 1 本の殻を返す。閾値は本周の実測から。
        # （eScholarship 2018B/1link・scholar.archive.org 4429B/1link・OSTI API 2B
        #  ↔ USGS 職員ページ 483 KB/644link・UBC Open Library 31 KB/24link）
        substantive = len(r["body"]) >= 10000 and len(ls) >= 10
        index_substantive = index_substantive or substantive
        print(f"  [索引] ✓ {url}\n         status={r['status']} ctype={r['ctype']} "
              f"bytes={len(r['body'])} links={len(ls)} "
              f"中身={'あり' if substantive else '**なし（JS の殻）**'} ({r['dt']}s)")
        # 本文語が生ページに出るかも見る（索引が JS で描かれている場合の見分け）
        raw = r["body"].decode("utf-8", "replace").lower()
        hit_words = [k for k in spec["must_any"] if k in raw]
        print(f"         生ページに著者語: {hit_words if hit_words else 'なし'}")
        author_in_raw = author_in_raw or bool(hit_words)
        for label, u in ls:
            sc = score(label, u, spec)
            if sc >= 2:
                cands.append((sc, label[:90], u))
    if not index_ok:
        print(f"  → 判定: INDEX_UNREACHABLE（索引そのものに 1 本も届かない）")
        return "INDEX_UNREACHABLE"
    cands = sorted({(s, l, u) for s, l, u in cands}, reverse=True)[:12]
    if not cands and not index_substantive:
        # **これを NOT_IN_INDEX と書いてはいけない**（欠陥 #53）。情報が無いだけである。
        print("  → 判定: INDEX_EMPTY_OR_JS（HTTP は通ったが中身の無い殻。**載っていないの意ではない**）")
        return "INDEX_EMPTY_OR_JS"
    if not cands and author_in_raw:
        # **著者名は載っている。リンクとして辿れないだけ**——これも「載っていない」ではない。
        print("  → 判定: IN_PAGE_NOT_LINKED（著者名は生 HTML に在るが、辿れるリンクの札に無い）")
        return "IN_PAGE_NOT_LINKED"
    if not cands:
        print(f"  → 判定: NOT_IN_INDEX（中身のある索引を取れたが、著者語がページのどこにも無い）")
        return "NOT_IN_INDEX"
    print(f"  [照合] 候補 {len(cands)} 件（加点順）:")
    for s, l, u in cands:
        print(f"    ({s}) {l}\n        {u}")
    # 段2: PDF らしい候補だけ叩く（上位 4 本まで）
    tried = 0
    for s, l, u in cands:
        if ".pdf" not in u.lower() or tried >= 4:
            continue
        tried += 1
        r = fetch(u)
        if not r["ok"]:
            print(f"  [本体] ✗ {u}\n         status={r['status']} {r['why']}")
            continue
        body = r["body"]
        if body[:5] == b"%PDF-":
            os.makedirs(OUTDIR, exist_ok=True)
            path = os.path.join(OUTDIR, spec["tag"] + ".pdf")
            with open(path, "wb") as f:
                f.write(body)
            print(f"  [本体] ★ PDF を取得: {u}\n         bytes={len(body)} "
                  f"sha256={hashlib.sha256(body).hexdigest()[:16]} → {path}")
            return "GOT_PDF"
        print(f"  [本体] ? PDF ではない: {u} ctype={r['ctype']} bytes={len(body)}")
    if tried == 0:
        print("  → 判定: NOT_IN_INDEX（候補はあるが PDF 直リンクが無い。上の候補を人が見る）")
        return "NOT_IN_INDEX"
    print("  → 判定: LINK_DEAD（索引に載っていたが本体が取れない）")
    return "LINK_DEAD"


def main():
    tee_stdout("step129")
    print("=== 旗129 所属機関索引プローブ（索引を作ってから引く版）===")
    print(time.strftime("%Y-%m-%d %H:%M:%S"))
    print(f"UA: {UA}")
    print("段1=索引を取る / 段2=索引の中を照合する。落ちた段を分けて記録する。")
    verdicts = {}
    for spec in TARGETS:
        try:
            verdicts[spec["tag"]] = run_target(spec)
        except Exception as e:  # 1 件の失敗で残りを落とさない
            print(f"  !! 例外: {type(e).__name__}: {e}")
            verdicts[spec["tag"]] = f"ERROR({type(e).__name__})"
    print(f"\n{'=' * 78}\n### まとめ")
    for tag, v in verdicts.items():
        print(f"  {tag:20s} {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
