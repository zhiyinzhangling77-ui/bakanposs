#!/usr/bin/env python3
"""旗131 — wayback の捕獲を **CDX で全部列挙**し、`id_` 付きで 1 つずつ取る。

## 何のための道具か

旗130 は Gaumont-Guay 2006 の全文リンクを 1 つ特定したが、**その 1 つの捕獲時刻
（20080708193216）が 130 KB で切断**し、Range でも継げなかった（`TRUNCATED`）。
**「この捕獲が壊れている」と「本文への経路が無い」は別である**（欠陥 #40 の作法）。
**捕獲は 1 つとは限らない。CDX API は同じ URL の全捕獲を列挙する。**

## 門①（対照）— **これを省かない**

判定器が「取れる／取れない」を言い分けられることを、**同じ走行の中で**示す：

- `CTRL_POS`（陽性対照）: **本文に一次到達済みの Davidson 1998**（旗124）。
  **同じ CDX 経路で完本が取れなければ、この道具が壊れている**（白の失敗を道具のせいにできない）。
- `CTRL_NEG`（陰性対照）: **実在しないパス**（同じホストの綴り違い）。
  **`NO_CAPTURE` が返らなければ、この道具は何にでも「在る」と言っている。**

## 判定語（旗130 の欠陥 #55・#56 の作法を継ぐ）

  NO_CAPTURE  CDX が 1 件も返さない（＝この URL は捕獲されていない）
  CDX_FAIL    CDX API 自体に届かない（＝情報が無い。「捕獲が無い」ではない）
  NOT_PDF     取れたが PDF ではない（中継 HTML・関門など）
  TRUNCATED   PDF だが末尾 %%EOF が無い／宣言長に届かない（**「取れた」と書かない**）
  GOT_PDF     %PDF- で始まり %%EOF で終わり、宣言長に達した完本

**出力は runlog でファイルに残す。**
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from runlog import tee_stdout  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, "tmp_pdfs")
UA = ("Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/124.0 Safari/537.36")

# (札, 種別, 元 URL, 本文にあるはずの語)
TARGETS = [
    ("Gaumont-Guay2006", "WHITE",
     "http://www.biometeorology.umn.edu/pdf/boreal_aspen_2006.pdf",
     ["soil respiration", "aspen", "temperature sensitivity"]),
    ("CTRL_POS_Davidson1998", "CTRL_POS",
     "http://pdfs.semanticscholar.org/f72a/8e9f706df71f0b7ce7750157f43792a2e66b.pdf",
     ["confounded", "soil respiration"]),
    ("CTRL_NEG_nonexistent", "CTRL_NEG",
     "http://www.biometeorology.umn.edu/pdf/boreal_aspen_2006_XXNOTREAL.pdf",
     []),
]


def http(url, timeout=90, start=None):
    """**部分読みを捨てない**——IncompleteRead は `partial` を持っている。"""
    h = {"User-Agent": UA, "Accept": "*/*"}
    if start:
        h["Range"] = f"bytes={start}-"
    req = urllib.request.Request(url, headers=h)
    with urllib.request.urlopen(req, timeout=timeout) as r:
        declared = r.headers.get("Content-Length")
        try:
            return r.read(), r.status, declared, None
        except Exception as e:
            return getattr(e, "partial", b""), r.status, declared, f"{type(e).__name__}: {e}"


def cdx_list(orig_url):
    """CDX API で捕獲を全部列挙。**新しい順に並べ替えて返す。**"""
    q = urllib.parse.urlencode({
        "url": orig_url,
        "output": "json",
        "fl": "timestamp,original,mimetype,statuscode,length",
        "filter": "statuscode:200",
        "collapse": "digest",
        "limit": "60",
    })
    api = "http://web.archive.org/cdx/search/cdx?" + q
    try:
        body, status, _, why = http(api, timeout=60)
    except Exception as e:
        return None, f"CDX_FAIL({type(e).__name__}: {e})"
    if why:
        return None, f"CDX_FAIL(切断 {why})"
    try:
        rows = json.loads(body.decode("utf-8", "replace") or "[]")
    except Exception as e:
        return None, f"CDX_FAIL(JSON 不正 {e})"
    if not rows:
        return [], "NO_CAPTURE"
    head, data = rows[0], rows[1:]
    idx = {k: i for i, k in enumerate(head)}
    out = []
    for r in data:
        out.append({
            "ts": r[idx["timestamp"]],
            "orig": r[idx["original"]],
            "mime": r[idx["mimetype"]],
            "len": r[idx.get("length", 0)] if "length" in idx else "?",
        })
    out.sort(key=lambda d: d["ts"], reverse=True)
    return out, ("NO_CAPTURE" if not out else "OK")


def fetch_capture(ts, orig):
    """`id_` 付きで原本を取る。**切断したら Range で 4 回まで継ぐ。**"""
    url = f"https://web.archive.org/web/{ts}id_/{orig}"
    buf, declared = b"", None
    for attempt in range(1, 6):
        try:
            chunk, status, dec, why = http(url, start=len(buf) if buf else None)
        except urllib.error.HTTPError as e:
            return url, b"", None, f"HTTPError {e.code}"
        except Exception as e:
            return url, buf, declared, f"{type(e).__name__}: {e}"
        if declared is None and dec and not buf:
            declared = int(dec)
        if not chunk:
            break
        buf += chunk
        if why is None:
            break
        time.sleep(1)
    return url, buf, declared, None


def verdict(buf, declared):
    if not buf:
        return "NO_BYTES"
    if buf[:5] != b"%PDF-":
        return "NOT_PDF"
    if b"%%EOF" not in buf[-4096:]:
        return "TRUNCATED"
    if declared and len(buf) < declared:
        return "TRUNCATED"
    return "GOT_PDF"


def pdf_words(path, words):
    """`pdftotext` は **subprocess 経由でのみ通る**（旗130）。"""
    try:
        txt = subprocess.run(["pdftotext", "-q", path, "-"],
                             capture_output=True, timeout=120).stdout.decode("utf-8", "replace")
    except Exception as e:
        return None, f"pdftotext 失敗: {type(e).__name__}: {e}"
    low = txt.lower()
    return txt, [w for w in words if w.lower() in low]


def main():
    tee_stdout("step131")
    print("=== 旗131 — CDX で全捕獲を列挙し、id_ 付きで 1 つずつ取る ===")
    print(time.strftime("%Y-%m-%d %H:%M:%S"))
    os.makedirs(OUTDIR, exist_ok=True)
    summary = []

    for name, kind, orig, words in TARGETS:
        print("\n" + "=" * 78)
        print(f"### {name} [{kind}]\n  元 URL: {orig}")
        caps, note = cdx_list(orig)
        if caps is None:
            print(f"  {note}")
            summary.append((name, kind, note.split("(")[0], ""))
            continue
        print(f"  CDX: {len(caps)} 件（digest 折り畳み後・statuscode:200）")
        for c in caps[:12]:
            print(f"    {c['ts']}  {c['mime']:<20} len={c['len']}")
        if not caps:
            print("  判定: NO_CAPTURE")
            summary.append((name, kind, "NO_CAPTURE", ""))
            continue
        if kind == "CTRL_NEG":
            # 陰性対照は「捕獲が返ってしまった」こと自体が異常。取りには行かない。
            print("  ⚠ 陰性対照に捕獲が返った——この道具は何にでも『在る』と言っている疑い")
            summary.append((name, kind, "UNEXPECTED_CAPTURE", ""))
            continue

        got = ""
        for i, c in enumerate(caps, 1):
            url, buf, declared, err = fetch_capture(c["ts"], c["orig"])
            v = verdict(buf, declared) if not err else f"ERR({err})"
            print(f"  [{i}/{len(caps)}] {c['ts']} → {len(buf)}B "
                  f"宣言={declared} 判定={v}")
            if v != "GOT_PDF":
                time.sleep(1)
                continue
            path = os.path.join(OUTDIR, f"{name}.pdf")
            with open(path, "wb") as f:
                f.write(buf)
            sha = hashlib.sha256(buf).hexdigest()[:16]
            print(f"    ★ 完本 bytes={len(buf)} sha256={sha} → {path}")
            txt, hit = pdf_words(path, words)
            if txt is None:
                print(f"    {hit}")
            else:
                with open(os.path.join(OUTDIR, f"{name}.txt"), "w") as f:
                    f.write(txt)
                print(f"    本文の語 {hit} / {words}（抽出 {len(txt)} 文字）")
                print("    先頭 400 字: " + " ".join(txt[:400].split()))
            got = c["ts"]
            break
        summary.append((name, kind, "GOT_PDF" if got else "ALL_CAPTURES_FAILED", got))

    print("\n" + "=" * 78)
    print("### まとめ")
    for name, kind, v, ts in summary:
        print(f"  {name:<24} {kind:<9} {v:<20} {ts}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
