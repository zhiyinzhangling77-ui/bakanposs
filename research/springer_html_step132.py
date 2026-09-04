#!/usr/bin/env python3
"""旗132 段3: wayback が持つ Springer 記事 HTML の捕獲に、本文が入っているかを見る。

非 OA の Springer 記事ページは通常「抄録まで」。本文が入っていれば一次到達。
**判定は目視ではなく語で行う**: 節見出し（Methods/Results/Discussion）と
本文にしか出ない語（引用文献リスト・図の説明）を数える。

門①（対照）: 同じ処理を、**本文が入っていると分かっている捕獲**に当てる必要がある。
ここでは対照として、上と同じ Springer ドメインの **OA 論文**の記事ページ捕獲を使う
（s10533-016-0233-5 は wayback に PDF 1.5MB があり OA の公算が高い）。
陽性対照が「本文なし」と出るなら、判定器のほうが壊れている。
"""
import html
import os
import re
import sys
import urllib.request
import urllib.error

UA = "Mozilla/5.0 (X11; Linux x86_64) research-primary-check/1.0"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tmp_pdfs")


def http(url, timeout=120):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, (e.read()[:1000] if e.fp else b"")
    except Exception as e:
        return None, repr(e).encode()


def detag(b):
    s = b.decode("utf-8", "replace")
    s = re.sub(r"(?is)<(script|style|noscript)\b.*?</\1>", " ", s)
    s = re.sub(r"(?s)<!--.*?-->", " ", s)
    s = re.sub(r"(?s)<[^>]+>", " ", s)
    return re.sub(r"\s+", " ", html.unescape(s)).strip()


MARKERS = ["Materials and methods", "Methods", "Results", "Discussion",
           "Acknowledg", "References", "Conclusion", "Study site", "Fig. 1", "Table 1"]
PAYWALL = ["This is a preview of subscription content",
           "Buy this article", "Institutional subscriptions", "log in via an institution",
           "Access this article", "Subscribe and save"]


def score(tag, url, save=None):
    st, body = http(url)
    print("=" * 74)
    print("[%s] %s" % (tag, url))
    print("  http:", st, "bytes:", len(body))
    if st != 200:
        return None
    txt = detag(body)
    print("  text chars:", len(txt))
    hit = [m for m in MARKERS if m.lower() in txt.lower()]
    pw = [p for p in PAYWALL if p.lower() in txt.lower()]
    print("  節見出し hit:", hit)
    print("  paywall 文言 hit:", pw)
    if save:
        p = os.path.join(OUT, save)
        with open(p, "w") as f:
            f.write(txt)
        print("  saved:", p)
    return txt


# 本命
t = score("TUCKER_2025capture",
          "https://web.archive.org/web/20250815010715id_/https://link.springer.com/article/10.1007/s10533-016-0200-1",
          save="tucker2016_springer_20250815.txt")
if t:
    for probe in ["negative temperature sensitivity", "Q10", "soil moisture", "abstract",
                  "Canyonlands", "Moab", "Utah", "volumetric water content"]:
        n = t.lower().count(probe.lower())
        print("   語 '%s': %d 回" % (probe, n))
    print("\n---- 先頭 2500 字 ----")
    print(t[:2500])

# 門①: 陽性対照（同ホスト・同誌・wayback に PDF がある＝OA の公算）
c = score("CTRL_OA_s10533-016-0233-5",
          "https://web.archive.org/web/2018/https://link.springer.com/article/10.1007%2Fs10533-016-0233-5")
sys.stdout.flush()
