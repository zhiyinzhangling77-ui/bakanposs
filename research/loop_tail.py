#!/usr/bin/env python3
"""claude --output-format stream-json の生ログを人が読める形にする。

run_loop.sh から 1 周ごとに呼ばれる。壊れた行・未知の形は黙って飛ばす
（表示の都合でループを落とさない）。
"""
import json
import sys

MAX_TEXT = 1200   # 1 発話あたりの表示上限
MAX_TOOLS = 60    # ツール呼び出しの表示上限


def brief(tool: str, inp: dict) -> str:
    """ツール呼び出しを 1 行に潰す。"""
    if not isinstance(inp, dict):
        return tool
    for key in ("command", "file_path", "pattern", "path", "url"):
        val = inp.get(key)
        if isinstance(val, str):
            val = " ".join(val.split())
            if len(val) > 100:
                val = val[:100] + "…"
            return f"{tool}: {val}"
    return tool


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: loop_tail.py <stream.jsonl>", file=sys.stderr)
        return 2

    texts, tools, result = [], [], None
    try:
        with open(sys.argv[1], encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line or not line.startswith("{"):
                    continue
                try:
                    ev = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if not isinstance(ev, dict):
                    continue

                kind = ev.get("type")
                if kind == "assistant":
                    content = (ev.get("message") or {}).get("content")
                    for block in content if isinstance(content, list) else []:
                        if not isinstance(block, dict):
                            continue
                        if block.get("type") == "text" and block.get("text", "").strip():
                            texts.append(block["text"].strip())
                        elif block.get("type") == "tool_use":
                            tools.append(brief(block.get("name", "?"), block.get("input")))
                elif kind == "result":
                    result = ev
    except OSError as exc:
        print(f"  (ログを読めない: {exc})", file=sys.stderr)
        return 1

    if tools:
        print(f"  ツール呼び出し {len(tools)} 件:")
        for t in tools[:MAX_TOOLS]:
            print(f"    · {t}")
        if len(tools) > MAX_TOOLS:
            print(f"    … 他 {len(tools) - MAX_TOOLS} 件")

    if texts:
        last = texts[-1]
        if len(last) > MAX_TEXT:
            last = last[:MAX_TEXT] + "\n    …(以下略・全文は生ログに)"
        print("  最後の発言:")
        for row in last.splitlines():
            print(f"    {row}")

    if isinstance(result, dict):
        bits = []
        if result.get("subtype"):
            bits.append(str(result["subtype"]))
        if isinstance(result.get("num_turns"), int):
            bits.append(f"{result['num_turns']} turns")
        if isinstance(result.get("duration_ms"), (int, float)):
            bits.append(f"{result['duration_ms'] / 1000:.0f}s")
        if isinstance(result.get("total_cost_usd"), (int, float)):
            bits.append(f"${result['total_cost_usd']:.2f}")
        if bits:
            print(f"  [{' / '.join(bits)}]")
        if result.get("is_error"):
            print("  ⚠ claude はエラーで終了している")

    if not (texts or tools or result):
        print("  (読み取れる出力が無い。生ログを見ること)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
