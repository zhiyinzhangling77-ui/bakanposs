"""作成済みの章 Markdown から学習用の派生物を作る(mdを作った後の工程)。

- build_index    : 各章の3行要約+重要語を1枚に集めた _INDEX.md(“地図”)。API不要=0トークン。
- generate_notes : 各章から Claude で凝縮スタディノートを作る(API必要・一度きりの投資)。

トークン効率の考え方: まず _INDEX.md(地図)を常時ロード → 必要な章だけ精読 →
凝縮ノートを一度だけ作って以後はそれで学習、が最小トークンで最大理解。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ChapterInfo:
    path: Path
    title: str
    chapter: int | None
    summary: list[str] = field(default_factory=list)   # 3行要約(各行 "- ...")
    keywords: str = ""                                  # 重要語(1行)
    body: str = ""                                      # 本文(# 見出し以降)


def chapter_files(folder: Path) -> list[Path]:
    """章ファイル(NN_*.md)。_INDEX.md / _conversion_log.md 等(_始まり)は除外。"""
    return sorted(p for p in folder.glob("*.md") if not p.name.startswith("_"))


def parse_chapter(md_text: str) -> tuple[str | None, int | None, list[str], str, str]:
    lines = md_text.splitlines()
    title: str | None = None
    chapter: int | None = None

    # --- YAML frontmatter ---
    i = 0
    if lines and lines[0].strip() == "---":
        i = 1
        while i < len(lines) and lines[i].strip() != "---":
            l = lines[i]
            if l.startswith("title:"):
                title = l.split(":", 1)[1].strip().strip('"')
            elif l.startswith("chapter:"):
                v = l.split(":", 1)[1].strip()
                chapter = int(v) if v.isdigit() else None
            i += 1

    # --- 要約 / 重要語 / 本文 ---
    summary: list[str] = []
    keywords = ""
    body_lines: list[str] = []
    state: str | None = None
    in_body = False
    for l in lines:
        s = l.strip()
        if not in_body and s.startswith("# "):
            in_body = True
        if in_body:
            body_lines.append(l)
            continue
        if s.startswith("### 3行要約"):
            state = "summary"
            continue
        if s.startswith("### 重要語"):
            state = "keywords"
            continue
        if state == "summary" and s.startswith("-"):
            summary.append(s)
        elif state == "keywords" and s and s != "---" and not s.startswith(">"):
            if not keywords:
                keywords = s

    return title, chapter, summary, keywords, "\n".join(body_lines).strip()


def read_chapters(folder: Path) -> list[ChapterInfo]:
    infos: list[ChapterInfo] = []
    for p in chapter_files(folder):
        title, chapter, summary, keywords, body = parse_chapter(
            p.read_text(encoding="utf-8", errors="replace")
        )
        infos.append(
            ChapterInfo(
                path=p, title=title or p.stem, chapter=chapter,
                summary=summary, keywords=keywords, body=body,
            )
        )
    return infos


def build_index(folder: Path, book_title: str | None = None) -> Path:
    """_INDEX.md を作る。Obsidian の [[wikilink]] で各章へ飛べる“地図”。"""
    infos = read_chapters(folder)
    name = book_title or folder.name
    out: list[str] = [f"# 📚 {name} — 章インデックス(要約マップ)\n"]
    out.append(f"章数: {len(infos)}\n")
    for info in infos:
        link = f"[[{info.path.stem}|{info.title}]]"
        out.append(f"## {link}")
        if info.summary:
            out.extend(info.summary)
        else:
            out.append("- (要約なし — `pdf2md run` で要約を付けると充実します)")
        if info.keywords:
            out.append(f"\n**重要語**: {info.keywords}")
        out.append("")  # 空行
    index_path = folder / "_INDEX.md"
    index_path.write_text("\n".join(out).rstrip() + "\n", encoding="utf-8")
    return index_path


# ---------------------------------------------------------------- 凝縮ノート

NOTE_MAX_INPUT_CHARS = 45_000

NOTE_SYSTEM = (
    "あなたは学習コーチです。与えられた章から、学習者が後で見返す"
    "『自分用の凝縮ノート』を日本語で作ります。冗長な引用はせず要点を圧縮します。"
)

NOTE_TEMPLATE = """次の章の本文から、凝縮スタディノートを**日本語**で作ってください。

出力フォーマット(この見出しを使う):
## 重要概念(定義つき)
- 用語 — 1行の定義
## 概念どうしの関係・全体像
- (概念のつながりを2〜5点)
## つまずきやすい点・よくある誤解
- (2〜4点)
## 理解確認クイズ(5問)
1. …
2. …
3. …
4. …
5. …
## クイズの答え
1. …(以下略)

--- 章本文ここから ---
{body}
--- 章本文ここまで ---"""


def generate_note(body: str, model: str = "claude-opus-4-8", max_tokens: int = 2000):
    """1章分の凝縮ノート本文を返す。(text, ok, note)"""
    try:
        import anthropic
    except Exception:
        return "", False, "anthropic SDK 未インストール"

    text = body.strip()
    if len(text) > NOTE_MAX_INPUT_CHARS:
        text = text[:NOTE_MAX_INPUT_CHARS]
    try:
        client = anthropic.Anthropic()
        resp = client.messages.create(
            model=model,
            max_tokens=max_tokens,
            system=NOTE_SYSTEM,
            messages=[{"role": "user", "content": NOTE_TEMPLATE.format(body=text)}],
        )
        out = "".join(b.text for b in resp.content if b.type == "text").strip()
        if not out:
            raise RuntimeError("空の応答")
        return out, True, ""
    except anthropic.AuthenticationError:
        return "", False, "APIキー未設定/無効(ANTHROPIC_API_KEY か ant auth login)"
    except Exception as e:  # noqa: BLE001
        return "", False, f"生成エラー: {e}"


@dataclass
class NotesReport:
    made: int = 0
    skipped: int = 0
    failed: int = 0
    notes: list[str] = field(default_factory=list)


def generate_notes(
    folder: Path,
    model: str = "claude-opus-4-8",
    overwrite: bool = False,
) -> NotesReport:
    """各章 → 凝縮ノートを folder/_notes/NN_..._notes.md に保存。"""
    report = NotesReport()
    notes_dir = folder / "_notes"
    notes_dir.mkdir(exist_ok=True)
    for info in read_chapters(folder):
        dest = notes_dir / f"{info.path.stem}_notes.md"
        if dest.exists() and not overwrite:
            report.skipped += 1
            continue
        text, ok, note = generate_note(info.body or info.path.read_text(
            encoding="utf-8", errors="replace"), model=model)
        if not ok:
            report.failed += 1
            report.notes.append(f"{info.path.name}: {note}")
            continue
        header = (
            f"---\nsource_chapter: \"{info.path.name}\"\ntype: study-note\n---\n\n"
            f"# 凝縮ノート: {info.title}\n\n"
        )
        dest.write_text(header + text + "\n", encoding="utf-8")
        report.made += 1
    return report
