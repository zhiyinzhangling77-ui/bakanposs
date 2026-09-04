#!/usr/bin/env python3
"""旗132 段1: Tucker & Reed 2016 の PDF 元 URL を当てる。

USGS Publications Warehouse の pubs-services API を生 HTTP で引き、
`links[]` に PDF があるかを見る。**要約モデルを介さない**（旗128 の作法）。
出力は生 JSON の必要部分のみ。判断はしない。
"""
import json
import sys
import urllib.request
import urllib.error

UA = "Mozilla/5.0 (X11; Linux x86_64) research-primary-check/1.0"


def get(url, timeout=45):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()[:2000]
    except Exception as e:
        return None, repr(e).encode()


CANDS = [
    # 記事 ID 直引き（旗128 で抄録に届いた 70169003）
    "https://pubs.usgs.gov/pubs-services/publication/70169003?mimetype=json",
    # 検索（ID が違っていた場合の保険）
    "https://pubs.usgs.gov/pubs-services/publication/?q=Low%20soil%20moisture%20during%20hot%20periods%20apparent%20negative%20temperature%20sensitivity&mimetype=json",
]

for url in CANDS:
    st, body = get(url)
    print("=" * 70)
    print("URL:", url)
    print("status:", st, "bytes:", len(body))
    if st != 200:
        print("body head:", body[:400])
        continue
    try:
        j = json.loads(body)
    except Exception as e:
        print("NOT_JSON:", repr(e), body[:400])
        continue
    recs = j.get("records", [j]) if isinstance(j, dict) else []
    print("n_records:", len(recs))
    for rec in recs[:5]:
        print("-" * 60)
        print("indexId:", rec.get("indexId"), "| year:", rec.get("publicationYear"))
        print("title:", (rec.get("title") or "")[:160])
        print("doi:", rec.get("doi"))
        for ln in rec.get("links", []) or []:
            print("   link:", ln.get("type", {}).get("text") if isinstance(ln.get("type"), dict) else ln.get("type"),
                  "|", ln.get("url"), "|", ln.get("linkFileType"))
sys.stdout.flush()
