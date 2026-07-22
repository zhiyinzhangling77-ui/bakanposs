# Colab で回す手順(GPU が無い / CPU だと重いとき)

Google Colab の無料 GPU(T4)で pdf2md を回すための手順です。各ステップを
Colab のセルにコピペして順に実行してください。

---

## 0. ランタイムを GPU に
Colab メニュー → **ランタイム → ランタイムのタイプを変更 → ハードウェアアクセラレータ = T4 GPU** → 保存。

確認:
```python
!nvidia-smi -L
```

## 1. Google Drive をマウント
PDF フォルダと Obsidian 保存先を Drive 上に置くと、結果が消えません。
```python
from google.colab import drive
drive.mount('/content/drive')
# 例:
#   入力: /content/drive/MyDrive/pdfs
#   出力: /content/drive/MyDrive/ObsidianVault/textbooks
```

## 2. リポジトリを取得して環境構築
```bash
%%bash
git clone <このリポジトリのURL> repo
cd repo/pdf2md
# Colab は既に CUDA 対応 torch が入っていることが多いので、まず素の依存を入れる
pip -q install marker-pdf mineru PyMuPDF anthropic
python - <<'PY'
import shutil, torch
print("marker_single:", bool(shutil.which("marker_single")))
print("mineru:", bool(shutil.which("mineru")))
print("cuda:", torch.cuda.is_available())
PY
```
> torch が CPU 版になってしまう場合のみ:
> `pip install torch --index-url https://download.pytorch.org/whl/cu121` を先に実行。

## 3. Anthropic API キー(要約用)
Colab の左メニュー「🔑 シークレット」に `ANTHROPIC_API_KEY` を登録してから:
```python
import os
from google.colab import userdata
os.environ["ANTHROPIC_API_KEY"] = userdata.get("ANTHROPIC_API_KEY")
```
(要約が不要なら省略可。`run` に `--no-summary` を付ける)

## 4. サンプル → 確認 → 本処理
```bash
%%bash
cd repo
IN=/content/drive/MyDrive/pdfs
OUT=/content/drive/MyDrive/ObsidianVault/textbooks

# ① 品質サンプル(先頭5ページ)
python -m pdf2md sample --input "$IN" --pages 5
```
出力を目視で確認し、良ければ:
```bash
%%bash
cd repo
IN=/content/drive/MyDrive/pdfs
OUT=/content/drive/MyDrive/ObsidianVault/textbooks

# ② 目次
python -m pdf2md toc --input "$IN"

# ③ 本処理(全章)。特定章だけなら --chapters "1,3,5-7"
python -m pdf2md run --input "$IN" --output "$OUT" --reviewed-sample --yes
```
> Colab は対話入力(上書き確認)が使いにくいので、本処理では `--yes` を付けて
> 上書き確認を自動承認しています。保存先を間違えないよう注意してください。

## 5. 結果
- 章ファイル: `<OUT>/<書名>/NN_章.md`
- ログ: `<OUT>/_conversion_log.md`(成功/失敗の一覧)

Drive 上に出力しておけば、そのまま Obsidian(Drive 同期)から開けます。

---

### ヒント
- 日本語・中国語 PDF で marker の結果が悪いときは、`sample`/`run` に `--prefer mineru`。
- 大きな本は時間がかかります。まず `--chapters` で数章だけ試すのがおすすめ。
- Colab のセッションが切れても、Drive 上の出力と `_conversion_log.md` は残ります。
