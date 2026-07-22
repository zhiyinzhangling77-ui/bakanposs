# pdf2md — PDF(教科書・論文)→ AIで勉強しやすい Markdown 変換パイプライン

PDF の教科書・論文を、Obsidian で扱いやすく・AI に読ませやすい Markdown に変換する
再利用可能なパイプラインです。

**やること**
- marker-pdf を主バックエンド、MinerU をフォールバック(日本語・中国語 PDF 向け)として変換
- GPU があれば使い、無ければ CPU で動く(重い場合は Colab 手順あり → 下記)
- **まず品質サンプル(先頭5ページ)→ 目視で確認 → 本処理** の順を守る
- 目次(見出し構造)を一覧表示し、`--chapters` で「この章だけ」を指定可能(未指定なら全章)
- トップレベル見出し(`# 章`)ごとに別ファイルへ分割
- 索引・参考文献の羅列・ページ番号などの定型ノイズを除去(**数式 LaTeX と表 Markdown は残す**)
- 各章の冒頭に **Claude 作の「3行要約 + 重要語10個」** を付与
- 失敗ファイルはスキップし、最後に **成功/失敗の一覧**を `_conversion_log.md` に出力
- 既存ファイルを上書きする前には**必ず確認**(自動化したいときだけ `--yes`)

---

## 1. セットアップ(専用 venv)

```bash
cd pdf2md
bash setup_venv.sh                 # CPU / 既定
# GPU(CUDA)がある場合:
PDF2MD_CUDA=1 bash setup_venv.sh   # CUDA 対応 torch を先に入れてから依存を入れる

source .venv/bin/activate
export ANTHROPIC_API_KEY=sk-ant-...   # 要約に必要。未設定なら要約はスキップされ本文は保存されます
```

`marker-pdf` と `mineru` は初回実行時にモデルをダウンロードします(数百 MB〜)。

---

## 2. 使い方(サンプル → 確認 → 本処理)

### ① まず品質サンプル(本処理には進まない)
```bash
python -m pdf2md sample --input "<PDFフォルダ>" --pages 5
# 保存もしたいとき: --save _samples/sample.md
```
表示された Markdown の品質(数式・表・見出し・文字化けの有無)を目で確認します。
**ここで OK になるまで本処理に進みません。**

日本語・中国語で marker の結果が悪ければ `--prefer mineru` を試してください。

### ② 目次(見出し構造)を確認
```bash
python -m pdf2md toc --input "<PDFフォルダ>"
```
各 PDF の章に `[番号]` が付きます。この番号を `--chapters` に使います。

### ③ 本処理(サンプル確認後)
```bash
# 全章を変換
python -m pdf2md run --input "<PDFフォルダ>" --output "<Obsidian保存先>" --reviewed-sample

# 特定の章だけ(例: 1章・3章・5〜7章)
python -m pdf2md run --input "<PDFフォルダ>" --output "<Obsidian保存先>" \
    --reviewed-sample --chapters "1,3,5-7"

# 1冊だけ / 要約なし / モデル変更
python -m pdf2md run --input "<PDFフォルダ>" --output "<保存先>" --reviewed-sample \
    --only "book.pdf" --no-summary --model claude-opus-4-8
```

- `--reviewed-sample` が無いと本処理は始まりません(サンプル確認の徹底のため)。
- 保存先に既存 `.md` があると、上書き確認のプロンプトが出ます(`--yes` で自動承認)。
- 出力は `<保存先>/<書名>/NN_章タイトル.md`、ログは `<保存先>/_conversion_log.md`。

---

## 3. 出力される章ファイルの形

```markdown
---
source: "book.pdf"
chapter: 3
title: "第3章 …"
tags: [study/summary]
summary_generated_by: claude
---

> **この章の要点(Claude 生成)**

### 3行要約
- …
- …
- …
### 重要語
語1、語2、…、語10

---

# 第3章 …
（本文。数式 $…$ と表 | … | はそのまま残ります）
```

---

## 4. Claude 要約について

- 各章のテキストを Anthropic API に送り、`3行要約 + 重要語10個`を生成します(既定モデル `claude-opus-4-8`)。
- 章が長い場合は先頭 45,000 文字までを使い、その旨を注記します。
- コスト重視なら `--model claude-haiku-4-5` などに変更可(モデル ID は環境で有効なもの)。
- `ANTHROPIC_API_KEY` 未設定・API エラー時は、その章の要約だけスキップして
  「未生成」ヘッダを付け、本文は保存します(処理は止まりません)。

---

## 5. GPU が無い/重いときは Colab で

`COLAB.md` を参照してください(コピペで回せる手順を記載)。要点:

1. Colab のランタイムを **GPU** に変更(ランタイム → ランタイムのタイプを変更 → T4 GPU)
2. Google Drive をマウントして PDF フォルダ・保存先を Drive 上に置く
3. このリポジトリを clone し、`PDF2MD_CUDA=1 bash pdf2md/setup_venv.sh` で環境構築
4. `ANTHROPIC_API_KEY` を Colab のシークレットで設定
5. `sample → toc → run` を同じコマンドで実行(保存先を Drive のパスにする)

---

## 6. 仕組み(モジュール構成)

| ファイル | 役割 |
|---|---|
| `pdf2md/convert.py` | marker/MinerU 呼び出し、GPU/CPU 判定、フォールバック |
| `pdf2md/toc.py` | PyMuPDF で目次(埋め込みアウトライン)抽出 |
| `pdf2md/clean.py` | ノイズ除去(索引・参考文献・ページ番号・繰り返しヘッダ) |
| `pdf2md/split.py` | `# 章`ごとの分割、章選択のパース |
| `pdf2md/summarize.py` | Claude 要約(3行要約+重要語10個)とヘッダ生成 |
| `pdf2md/pipeline.py` | sample/toc/run の本体、成功/失敗ログ |
| `pdf2md/cli.py` | コマンドライン(サンプル→確認→本処理のガード、上書き確認) |

## 注意 / 既知の制約
- ノイズ除去はヒューリスティックです。まれに本文の見出し名が「References/文献」等と
  一致すると誤って落とすことがあります。心配なときは `sample` の結果で確認してください。
- `--chapters` は「変換後 Markdown のトップレベル見出し」に対する番号です。埋め込み目次の
  番号と概ね対応しますが、PDF によってはズレることがあります(`toc` と `sample` で確認)。
