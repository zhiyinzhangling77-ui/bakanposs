#!/usr/bin/env python3
"""旗129 段4 — `fatcat`（scholar.archive.org の索引）を **残る ⚪ 5 件全部** に当てる。

## なぜこの索引か

旗128 の欠陥 #52：「手を一度使ったら、残りの ⚪ 全部に当ててから『経路が尽きた』と書く」。
段1〜3 で分かったこと：
  ・所属機関の刊行物一覧は **JS で描かれていて生 HTTP に中身が無い**（5 件中 3 件）。
  ・PDF 直リンク（UMN ミラー・AGU）は **生 HTTP でも 403**＝旗128 の「urllib なら通る」は
    **API に対しては真だが、UA を見るファイルサーバに対しては偽**。

`fatcat` は Unpaywall / OpenAlex とは**別の索引**で、**ウェブアーカイブ内の実ファイル**を持つ。
機械可読（JSON）なので JS の問題が無い。**旗128 が言う「索引 API に見えない経路」の逆側**——
索引 API のうち、まだ当てていないものである。

## 対照（旗128 の作法：判定器は既知の陽性に当ててから使う）

**本文まで一次到達した実績のある 5 件**を同じ関数に通す。
`fatcat` がそのうち何件で実ファイルを返すかを見てから、⚪ 側の結果を読む。

出力: research/logs/step129c_<timestamp>.txt
"""

import json
import os
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from runlog import tee_stdout  # noqa: E402

UA = "bakanposs-research/1.0 (mailto:zhiyinzhangling77@gmail.com)"

# 対照＝本文まで一次到達済み（旗121〜127）。判定器がこれを拾えるかを先に見る。
CONTROLS = [
    ("Zhang2018 AFM 259 (旗121 到達)", "10.1016/j.agrformet.2018.05.005"),
    ("Wang2014 BG 11 (旗122 到達)", "10.5194/bg-11-259-2014"),
    ("Wen2006 AFM 137 (旗123 到達)", "10.1016/j.agrformet.2006.02.005"),
    ("Davidson1998 GCB 4 (旗124 到達)", "10.1046/j.1365-2486.1998.00128.x"),
    ("Doerr&Muennich1987 Tellus 39B (旗127 到達)", "10.3402/tellusb.v39i1-2.15331"),
]

# 残る ⚪。DOI は Crossref/検索で確認できたもののみ。無いものは title 検索に落とす。
OPEN = [
    ("Xu&Qi2001 GBC 15 687-696", "10.1029/2000GB001365"),
    ("Gaumont-Guay2006 AFM 140 220-235", "10.1016/j.agrformet.2006.08.003"),
    ("Tucker&Reed2016 Biogeochem 128 155-169", "10.1007/s10533-016-0200-1"),
    ("Suseela2012 GCB 18 336-348", "10.1111/j.1365-2486.2011.02516.x"),
    ("Bunnell1977 SBB 9 33-40", "10.1016/0038-0717(77)90058-X"),
]


def get_json(url, timeout=45):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return json.loads(r.read().decode("utf-8", "replace")), None
    except urllib.error.HTTPError as e:
        return None, f"HTTP {e.code} {e.reason}"
    except Exception as e:
        return None, f"{type(e).__name__}: {e}"


def probe(label, doi):
    """fatcat に DOI を引いて、**実ファイルの URL があるか**だけを見る。"""
    url = ("https://api.fatcat.wiki/v0/release/lookup"
           f"?doi={urllib.parse.quote(doi)}&expand=files&hide=abstracts,refs")
    data, err = get_json(url)
    if err:
        print(f"  ✗ {label}\n      doi={doi}\n      lookup: {err}")
        return "NO_RECORD" if "404" in err else "API_ERROR"
    files = data.get("files") or []
    urls = []
    for f in files:
        for u in (f.get("urls") or []):
            urls.append((f.get("mimetype"), f.get("size"), u.get("rel"), u.get("url")))
    print(f"  {'★' if urls else '·'} {label}\n      doi={doi}  fatcat_ident={data.get('ident')}  "
          f"files={len(files)}  file_urls={len(urls)}")
    for mt, sz, rel, u in urls[:6]:
        print(f"        [{rel}] {mt} {sz}B  {u}")
    return "HAS_FILE" if urls else "RECORD_NO_FILE"


def main():
    tee_stdout("step129c")
    print("=== 旗129 段4：fatcat（scholar.archive.org 索引）を全 ⚪ に当てる ===")
    print(time.strftime("%Y-%m-%d %H:%M:%S"))
    print("**先に対照（本文到達済み 5 件）に当て、感度を見てから ⚪ を読む。**")

    print("\n--- 対照（本文まで一次到達した実績あり）---")
    c = {}
    for label, doi in CONTROLS:
        c[label] = probe(label, doi)
        time.sleep(0.5)
    hit = sum(1 for v in c.values() if v == "HAS_FILE")
    print(f"\n  感度: {hit}/{len(CONTROLS)}（実ファイルを持つ対照の数）")

    print("\n--- 残る ⚪ 5 件 ---")
    o = {}
    for label, doi in OPEN:
        o[label] = probe(label, doi)
        time.sleep(0.5)

    print("\n=== まとめ ===")
    print("  [対照]")
    for k, v in c.items():
        print(f"    {v:15s} {k}")
    print("  [⚪]")
    for k, v in o.items():
        print(f"    {v:15s} {k}")
    print("\n  **感度が低ければ、⚪ 側の `RECORD_NO_FILE` は到達不能の根拠にならない**"
          "（旗128 と同じ読み方）。")
    return 0


if __name__ == "__main__":
    import urllib.parse  # noqa: E402  (probe 内で使う)
    sys.exit(main())
