"""各章に付ける「Claude 作の 3行要約 + 重要語10個」を生成する。

Anthropic API(公式 SDK)を使う。ANTHROPIC_API_KEY(または `ant auth login` の
プロファイル)で認証する。キーが無ければ要約はスキップし、章本文はそのまま保存する。
"""

from __future__ import annotations

from dataclasses import dataclass

# 章が長すぎるときに送る先頭文字数(概ね 1万数千トークン相当)。
MAX_INPUT_CHARS = 45_000

SYSTEM_PROMPT = (
    "あなたは学習教材の編集者です。与えられた章のテキストを読み、"
    "日本語で簡潔かつ正確に要約します。誇張や創作はせず、本文にある内容だけを使います。"
)

USER_TEMPLATE = """次の「章」のテキストを読んで、下記フォーマットで**日本語**で出力してください。

- 3行要約: 各行は1文、先頭に "- " を付ける(ちょうど3行)
- 重要語: その章の重要語をちょうど10個、読点区切りで1行

出力フォーマット(この見出しをそのまま使う):
### 3行要約
- (1行目)
- (2行目)
- (3行目)
### 重要語
語1、語2、…、語10

--- 章テキストここから ---
{chapter}
--- 章テキストここまで ---"""


@dataclass
class SummaryResult:
    text: str            # 付与するヘッダ用テキスト(要約 + 重要語)
    ok: bool
    note: str = ""       # スキップ理由など


def _fallback_header(note: str) -> str:
    return (
        "### 3行要約\n"
        "- (要約は未生成: " + note + ")\n"
        "### 重要語\n"
        "(未生成)\n"
    )


def summarize_chapter(
    chapter_text: str,
    model: str = "claude-opus-4-8",
    max_tokens: int = 1024,
) -> SummaryResult:
    try:
        import anthropic
    except Exception:
        note = "anthropic SDK 未インストール"
        return SummaryResult(text=_fallback_header(note), ok=False, note=note)

    body = chapter_text.strip()
    truncated = False
    if len(body) > MAX_INPUT_CHARS:
        body = body[:MAX_INPUT_CHARS]
        truncated = True

    try:
        client = anthropic.Anthropic()  # キー/プロファイルは環境から解決
        resp = client.messages.create(
            model=model,
            max_tokens=max_tokens,
            system=SYSTEM_PROMPT,
            messages=[{"role": "user", "content": USER_TEMPLATE.format(chapter=body)}],
        )
        out = "".join(b.text for b in resp.content if b.type == "text").strip()
        if not out:
            raise RuntimeError("空の応答")
        if truncated:
            out += "\n\n> ※ 章が長いため先頭部分のみから要約しています。"
        return SummaryResult(text=out + "\n", ok=True)
    except anthropic.AuthenticationError:
        note = "APIキー未設定/無効(ANTHROPIC_API_KEY を設定してください)"
        return SummaryResult(text=_fallback_header(note), ok=False, note=note)
    except Exception as e:  # noqa: BLE001 — 要約失敗でも章保存は続ける
        note = f"要約生成エラー: {e}"
        return SummaryResult(text=_fallback_header(note), ok=False, note=note)


def build_chapter_header(
    source_pdf: str,
    chapter_title: str,
    chapter_index: int,
    summary: SummaryResult,
) -> str:
    """章ファイル冒頭に付ける YAML frontmatter + 要約ブロック。"""
    tags = "study/summary" + ("" if summary.ok else " study/summary-missing")
    return (
        "---\n"
        f'source: "{source_pdf}"\n'
        f"chapter: {chapter_index}\n"
        f'title: "{chapter_title}"\n'
        f"tags: [{tags.replace(' ', ', ')}]\n"
        f"summary_generated_by: {'claude' if summary.ok else 'none'}\n"
        "---\n\n"
        "> **この章の要点(Claude 生成)**\n\n"
        f"{summary.text}\n"
        "---\n\n"
    )
