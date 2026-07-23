"""Markdown をトップレベル見出し(# 章)ごとに分割する。"""

from __future__ import annotations

import re
from dataclasses import dataclass

H1_RE = re.compile(r"^#\s+(.*\S)\s*$")
# 見出しに混じる Markdown 強調記号(** * __ `)を除去するための正規表現
_EMPHASIS_RE = re.compile(r"(\*\*|__|\*|`)")


def clean_heading(text: str) -> str:
    """見出しテキストから強調記号を取り除く(YAML title / ファイル名用)。"""
    return _EMPHASIS_RE.sub("", text).strip()


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
            cur_title = clean_heading(m.group(1).strip())
            cur_body = [line]
        elif cur_title is None:
            frontmatter.append(line)
        else:
            cur_body.append(line)
    flush()
    return "\n".join(frontmatter).strip(), chapters


ANY_HEADING_RE = re.compile(r"^(#{1,6})\s+(.*\S)\s*$")

# mode="pattern" で --chapter-pattern 未指定のときに使う既定(章の見出しっぽいもの)
DEFAULT_CHAPTER_PATTERN = r"(?i)^\s*(chapter\s+\d+|第\s*\d+\s*[章節]|\d+[\.\s])"


def split_none(markdown: str, fallback_title: str) -> tuple[str, list[Chapter]]:
    """分割しない: 変換範囲まるごとを1章として返す。

    `--page-range` で1章ずつ変換する運用に最適(marker が節見出しを # にして
    しまい過剰分割する問題を回避できる)。タイトルは最初の見出し、無ければ
    fallback_title(通常は PDF 名)。
    """
    title = fallback_title
    for line in markdown.splitlines():
        m = ANY_HEADING_RE.match(line)
        if m:
            title = clean_heading(m.group(2).strip())
            break
    return "", [Chapter(index=1, title=title, body=markdown.strip() + "\n")]


def split_by_pattern(
    markdown: str, pattern: str | None
) -> tuple[str, list[Chapter]]:
    """見出し(任意レベル)のうち、正規表現に一致するものだけを章境界にする。"""
    rx = re.compile(pattern or DEFAULT_CHAPTER_PATTERN)
    lines = markdown.splitlines()
    frontmatter: list[str] = []
    chapters: list[Chapter] = []
    cur_title: str | None = None
    cur_body: list[str] = []
    in_code = False

    def flush():
        if cur_title is not None:
            chapters.append(
                Chapter(len(chapters) + 1, cur_title, "\n".join(cur_body).strip() + "\n")
            )

    for line in lines:
        if line.lstrip().startswith("```"):
            in_code = not in_code
        m = None if in_code else ANY_HEADING_RE.match(line)
        is_boundary = bool(m and rx.search(m.group(2).strip()))
        if is_boundary:
            flush()
            cur_title = clean_heading(m.group(2).strip())
            cur_body = [line]
        elif cur_title is None:
            frontmatter.append(line)
        else:
            cur_body.append(line)
    flush()
    return "\n".join(frontmatter).strip(), chapters


def split_chapters(
    markdown: str,
    mode: str = "h1",
    pattern: str | None = None,
    fallback_title: str = "section",
) -> tuple[str, list[Chapter]]:
    """分割モードを選んで (冒頭部, 章リスト) を返す。

    mode: "h1"(# ごと・従来) / "none"(分割しない) / "pattern"(章見出し正規表現)
    """
    if mode == "none":
        return split_none(markdown, fallback_title)
    if mode == "pattern":
        fm, chs = split_by_pattern(markdown, pattern)
        if chs:
            return fm, chs
        # パターンに1つも一致しなければ分割せず1章にまとめる(過剰分割を避ける)
        return split_none(markdown, fallback_title)
    return split_by_h1(markdown)


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
