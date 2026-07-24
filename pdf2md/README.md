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

`marker-pdf` は初回実行時にモデルをダウンロードします(数百 MB〜)。

**実行時の注意(ディレクトリ)**: `python -m pdf2md ...` は `pdf2md` パッケージが
見えるディレクトリ(=この `pdf2md/` フォルダ)から実行してください。リポジトリの
ルート等から実行すると `No module named pdf2md.__main__` になります。どこからでも
呼びたい場合は一度だけ editable install すると `pdf2md` コマンドが使えます:

```bash
cd pdf2md
./.venv/bin/pip install -e .
# 以後どのディレクトリからでも(venv 有効化中):  pdf2md sample --input ... など
```

**GPU メモリが少ない場合(例: 4GB)**: そのままだと CUDA out of memory になることが
あります。`--device cpu` を付ければ CPU で確実に動きます(自動でも OOM を検知したら
CPU で 1 回再試行します)。

**MinerU(日本語・中国語フォールバック)について**: MinerU は Pillow のバージョンが
marker と衝突するため、**同じ venv には入れません**。使う場合は別 venv を作ります:

```bash
python3 -m venv .venv-mineru
source .venv-mineru/bin/activate
pip install -r requirements-mineru.txt PyMuPDF anthropic
# この venv から実行すると marker が無いので自動的に mineru が使われます
python -m pdf2md sample --input <PDF> --pages 5 --prefer mineru
```

---

## 2. 使い方(サンプル → 確認 → 本処理)

### ① まず品質サンプル(本処理には進まない)
```bash
python -m pdf2md sample --input "<PDFフォルダ>" --pages 5
# 保存もしたいとき: --save _samples/sample.md
```
表示された Markdown の品質(数式・表・見出し・文字化けの有無)を目で確認します。
**ここで OK になるまで本処理に進みません。**

日本語・中国語で marker の結果が悪ければ、上記の別 venv で MinerU を試してください。
低VRAM GPU で落ちるときは `--device cpu`。`--input` はフォルダでも単一 .pdf でも可。

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

**速度の注意(重要)**: `--chapters` は「変換後に書き出す章を絞る」だけで、変換自体は
全ページ実行します。CPU で1冊フルは遅いので、**特定の章だけ速く変換したいときは
`--page-range` を使ってください**。`toc` に出るページ番号(1始まり)で指定します:

```bash
# 例: toc で第2章が p.15〜p.48 と分かったら、その範囲だけ変換(高速)
python -m pdf2md run --input <PDF> --output <保存先> --reviewed-sample \
    --device cpu --page-range "15-48"
```

**分割の注意 & 推奨運用**: 既定の `--split h1` は「`#` 見出しごとに1ファイル」ですが、
marker は**章の“節”見出しも `#`(H1)にする**ことが多く、1つの章が十数ファイルに
過剰分割されがちです。おすすめは次の2択:

- **`--page-range` で1章ずつ変換 + `--split none`**(範囲まるごと1ファイル)。
  → 章の切れ目を自分で決められて、過剰分割も起きない。**いちばん確実**。
  ```bash
  # 第2章(p.15〜48)を1ファイルとして書き出す
  python -m pdf2md run --input <PDF> --output <保存先> --reviewed-sample \
      --device cpu --page-range "15-48" --split none
  ```
- **`--split pattern`**: 「Chapter N」「第N章」「N.」のような**章見出しだけ**を境界にする
  (本が章番号を振っていて marker がそれを保持している場合に有効)。必要なら
  `--chapter-pattern '正規表現'` で境界を自分で指定。

**`--by-toc`(全書を1コマンドで章ごとに)**: PDFに埋め込み目次があれば、**章リストを
手入力せずに**、目次のトップレベル章ごとにそのページ範囲だけを変換し、1章1ファイルで
保存します(タイトルも目次から自動)。全書変換の一番ラクな方法:

```bash
python -m pdf2md run --input "$PDF" --output <保存先> --reviewed-sample \
    --device cpu --no-summary --by-toc
```

注意: 目次に「Part I」「Cover」等が章と同レベルで入っていると、それらも小さな
ファイルになります(不要なら後で削除)。埋め込み目次が無い PDF では使えないので、
その場合は `--page-range` を使ってください。`--by-toc` 使用時は
`--page-range/--split/--title` は無視されます。

- `--reviewed-sample` が無いと本処理は始まりません(サンプル確認の徹底のため)。
- 保存先に既存 `.md` があると、上書き確認のプロンプトが出ます(`--yes` で自動承認)。
- 出力は `<保存先>/<書名>/NN_章タイトル.md`、ログは `<保存先>/_conversion_log.md`。

---

## 2.5 学習用の派生物(mdを作った“後”の工程)

トークンを抑えて深く理解するための2コマンド。**変換で章 .md ができた後**に使う。

```bash
# 地図: 各章の3行要約+重要語を1枚に集めた _INDEX.md(API不要=0トークン)
python -m pdf2md index --dir "<保存先>/<書名>" --title "本のタイトル"

# 凝縮ノート: 各章 → Claude が自分用スタディノート(重要概念/関係/誤解/クイズ)
#   API必要(ANTHROPIC_API_KEY か ant auth login)。一度きりの投資。
python -m pdf2md notes --dir "<保存先>/<書名>" --model claude-opus-4-8
```

- `_INDEX.md` は Obsidian の `[[wikilink]]` で各章へ飛べる“地図”。**これを常時開き、必要な章だけ精読**するとトークンを最小化できる。
- `notes` は `<書名>/_notes/NN_..._notes.md` を作る(既存はスキップ、`--overwrite` で再生成)。
- おすすめの流れ: **要約ON で変換 → `index` で地図 → `notes` で凝縮 → その凝縮ノートで学習/自己テスト**。

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
