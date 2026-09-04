#!/usr/bin/env python3
"""旗132 段2: Tucker & Reed 2016 に CDX 経路を当てる（旗131 で確定した 2 本目の経路）。

経路（旗131）:
  1. CDX で捕獲を全部列挙  2. len 列を先に見る  3. 新しい順に id_ 付きで取る
  4. 完本判定 = %PDF- で始まり、末尾 4KB に %%EOF、宣言長に到達

門①（対照）を段の中に置く:
  陽性 = Davidson 1998（旗124/131 で本文到達済み）— 同じ道具が今も通ることを確かめる
  陰性 = 同ホストの実在しないパス — NO_CAPTURE が出ることを確かめる
★CDX の 503 は NO_CAPTURE と別語（CDX_FAIL）にする（道具の欠陥 #57）。
"""
import hashlib
import json
import os
import sys
import time
import urllib.parse
import urllib.request
import urllib.error

UA = "Mozilla/5.0 (X11; Linux x86_64) research-primary-check/1.0"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tmp_pdfs")
os.makedirs(OUT, exist_ok=True)


def http(url, timeout=90):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, dict(r.headers), r.read()
    except urllib.error.HTTPError as e:
        return e.code, dict(e.headers or {}), (e.read()[:1000] if e.fp else b"")
    except Exception as e:
        return None, {}, repr(e).encode()


def cdx(url, match=None, limit=60):
    q = {
        "url": url,
        "output": "json",
        "filter": "statuscode:200",
        "collapse": "digest",
        "limit": str(limit),
    }
    if match:
        q["matchType"] = match
    full = "http://web.archive.org/cdx/search/cdx?" + urllib.parse.urlencode(q)
    st, _, body = http(full, timeout=90)
    if st != 200:
        return "CDX_FAIL", st, []
    if not body.strip():
        return "NO_CAPTURE", st, []
    try:
        rows = json.loads(body)
    except Exception:
        return "CDX_FAIL", st, []
    if len(rows) < 2:
        return "NO_CAPTURE", st, []
    hdr, data = rows[0], rows[1:]
    recs = [dict(zip(hdr, r)) for r in data]
    return "HAVE", st, recs


def grab(orig, ts, tag):
    """id_ 付きで原本を取り、完本かを判定する。"""
    u = "https://web.archive.org/web/%sid_/%s" % (ts, orig)
    st, hdrs, body = http(u, timeout=180)
    if st != 200:
        return {"verdict": "HTTP_%s" % st, "bytes": len(body)}
    declared = hdrs.get("Content-Length")
    is_pdf = body[:5] == b"%PDF-"
    has_eof = b"%%EOF" in body[-4096:]
    complete = is_pdf and has_eof and (declared is None or len(body) >= int(declared))
    res = {
        "verdict": "GOT_PDF" if complete else ("TRUNCATED" if is_pdf else "NOT_PDF"),
        "bytes": len(body),
        "declared": declared,
        "sha256": hashlib.sha256(body).hexdigest()[:16],
        "head": body[:40].decode("latin-1", "replace"),
    }
    if is_pdf:
        path = os.path.join(OUT, "%s_%s.pdf" % (tag, ts))
        with open(path, "wb") as f:
            f.write(body)
        res["path"] = path
    return res


def run(tag, orig, match=None, take=3):
    print("=" * 74)
    print("[%s] %s" % (tag, orig) + ("  (matchType=%s)" % match if match else ""))
    state, st, recs = cdx(orig, match=match)
    print("  CDX:", state, "http=%s" % st, "n=%d" % len(recs))
    if state != "HAVE":
        return state
    recs.sort(key=lambda r: r.get("timestamp", ""), reverse=True)
    for r in recs[:12]:
        print("   ", r.get("timestamp"), r.get("length"), (r.get("mimetype") or "")[:24],
              (r.get("original") or "")[:110])
    if match:  # prefix 検索は列挙が目的。取得はしない
        return state
    got = 0
    for r in recs:
        if got >= take:
            break
        res = grab(orig, r["timestamp"], tag)
        print("   -> ts=%s %s" % (r["timestamp"], res))
        got += 1
        if res["verdict"] == "GOT_PDF":
            return "GOT_PDF"
        time.sleep(2)
    return "NO_FULLTEXT"


# ---- 段2a: 門①（対照）を先に通す --------------------------------------
print("### 門①（対照）")
DAVIDSON = "http://harvardforest.fas.harvard.edu/sites/harvardforest.fas.harvard.edu/files/publications/pdfs/Davidson_GlobalChangeBiology_1998.pdf"
pos = run("CTRL_POS_davidson1998", DAVIDSON)
neg = run("CTRL_NEG_fake", DAVIDSON.replace("Davidson_GlobalChangeBiology_1998",
                                            "Davidson_GlobalChangeBiology_XXNOTREAL"))
print("\n門①: 陽性=%s（期待 GOT_PDF） / 陰性=%s（期待 NO_CAPTURE）\n" % (pos, neg))

# ---- 段2b: 本命 ---------------------------------------------------------
print("### 本命 Tucker & Reed 2016 (doi 10.1007/s10533-016-0200-1)")
TARGETS = [
    ("springer_pdf", "https://link.springer.com/content/pdf/10.1007/s10533-016-0200-1.pdf", None),
    ("springer_pdf_enc", "http://link.springer.com/content/pdf/10.1007%2Fs10533-016-0200-1.pdf", None),
    ("springer_art", "https://link.springer.com/article/10.1007/s10533-016-0200-1", None),
    ("springer_art_enc", "http://link.springer.com/article/10.1007%2Fs10533-016-0200-1", None),
    # 索引を介さない広い列挙（取得はしない・在庫の有無を見るため）
    ("prefix_springer_thisdoi", "link.springer.com/content/pdf/10.1007/s10533-016-0200-1", "prefix"),
    ("prefix_springer_s10533_2016", "link.springer.com/content/pdf/10.1007/s10533-016-02", "prefix"),
]
results = {}
for tag, u, m in TARGETS:
    try:
        results[tag] = run(tag, u, match=m)
    except Exception as e:
        results[tag] = "ERROR %r" % e
    time.sleep(1)

print("\n### まとめ")
for k, v in results.items():
    print("  %-28s %s" % (k, v))
sys.stdout.flush()
