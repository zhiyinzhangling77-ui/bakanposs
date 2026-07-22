"""Markdown をトップレベル見出し(# 章)ごとに分割する。"""

from __future__ import annotations

import re
from dataclasses import dataclass

H1_RE = re.compile(r"^#\s+(.*\S)\s*$")


@dataclass
class Chapter:
    index: int          # 1 始まりの章番号
    title: str
    body: str           # 見出し行を含む章本文


def _sanitize(title: str, maxlen: int = 60) -> str:
    s = re.sub(r"\s+", "_", title.strip())
    s = re.sub(r'[\\/:*?"<>|#]+', "", s)
    s = s.strip("_.")
    return (s or "section")[:maxlen]


def split_by_h1(markdown: str) -> tuple[str, list[Chapter]]:
    """(冒頭部, 章リスト)を返す。冒頭部は最初の # より前の内容。"""
    lines = markdown.splitlines()
    frontmatter: list[str] = []
    chapters: list[Chapter] = []
    cur_title: str | None = None
    cur_body: list[str] = []
    in_code = False

    def flush():
        if cur_title is not None:
            chapters.append(
                Chapter(
                    index=len(chapters) + 1,
                    title=cur_title,
                    body="\n".join(cur_body).strip() + "\n",
                )
            )

    for line in lines:
        if line.lstrip().startswith("```"):
            in_code = not in_code
        m = None if in_code else H1_RE.match(line)
        if m:
            flush()
            cur_title = m.group(1).strip()
            cur_body = [line]
        elif cur_title is None:
            frontmatter.append(line)
        else:
            cur_body.append(line)
    flush()
    return "\n".join(frontmatter).strip(), chapters


def chapter_filename(ch: Chapter) -> str:
    return f"{ch.index:02d}_{_sanitize(ch.title)}.md"


def parse_selection(spec: str | None, n: int) -> list[int]:
    """--chapters '1,3,5-7' を 1 始まりインデックスのリストに。None なら全章。"""
    if not spec:
        return list(range(1, n + 1))
    picked: set[int] = set()
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            a, _, b = part.partition("-")
            picked.update(range(int(a), int(b) + 1))
        else:
            picked.add(int(part))
    return sorted(i for i in picked if 1 <= i <= n)
