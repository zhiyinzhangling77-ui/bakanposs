あなたはこのリポジトリの研究（`japanflux_pn` / 旗方式）を 1 周だけ進めます。
**あなたには一切の会話履歴がありません。** 状態はすべてファイルにあります。
このプロセスは 1 周で終了し、次の周は別プロセスとして起動されます。

## 1. 状態を復元する

```bash
cat research/SESSION_STATE.md
cat research/HUMAN_GATES.md
cat research/OPEN_QUESTIONS_OPTIONS.md
cat research/LOOP_PROTOCOL.md
git log --oneline -20
git status --short
```

**`research/FLAGS_LOG.md` は 4000 行超。全文を読まない。**
`grep -n '^## 旗' research/FLAGS_LOG.md | tail -20` で索引を引き、`sed -n 'X,Yp'` で必要な旗だけ。

`research/COLD_START.md` の Step 2（整合性チェック）を実行し、
**SESSION_STATE の最新旗と FLAGS_LOG の最終旗がずれていたら、それを直すことを今周の作業にする。**

## 2. この周でやること — **ちょうど 1 単位**

`LOOP_PROTOCOL.md` の「ループ本体」に従い、**BLOCKED でない最優先の作業を 1 つだけ**やる。
1 単位とは「1 つの旗を閉じる」か「1 つの事前登録を書き上げる」か「1 つの道具を作り切る」。

**このプロセス内で終わらない量に着手しないこと。** 大きすぎるなら分割して、今周は最初の 1 つだけやる。

**あなたはコンテナ内かもしれないし、実データのあるローカルかもしれない。** 最初に確かめること:

```bash
ls /mnt/hdd/ 2>/dev/null && echo "LOCAL: 実データあり" || echo "CONTAINER: 実データなし"
```

- **実データがある**なら、`HUMAN_GATES.md` の **D 分類（GATE-01〜04, 06）はゲートではない**。自分で実行してよい
- **実データが無い**なら、D 分類は従来どおりゲート。合成検証・道具作り・文章に徹する

## 3. 終わり方 — 必ずどちらか

### (a) 作業を 1 単位進めた場合

1. `research/FLAGS_LOG.md` に旗を追記（または今周の記録を該当文書へ）
2. **`research/SESSION_STATE.md` を更新**（これをせずに周を終えない）
3. 失敗した試行・道具の欠陥（番号を継ぐ）・予測の勝敗も記録する
4. コミットして push:
   ```bash
   git add -A
   git commit -m "旗NN <一行>"
   git push -u origin claude/new-branch-creation-heq0i1
   ```
5. **`research/.loop_stop` は作らない。** 何も言わずに終了してよい

### (b) 残りが**全部** BLOCKED だった場合（＝人間の手が要る）

**1 つブロックされただけで (b) にしないこと。** 手元でできる作業が 1 つでもあるなら (a) を選ぶ。

本当に全部ブロックされているときだけ、`research/.loop_stop` を書いて終了する:

```bash
cat > research/.loop_stop <<'EOF'
理由: 残りが全部人間待ち

## お願いすること
### バッチ1: ローカル実行（夜間放置可）
GATE-xx, GATE-yy
<コピペで動くコマンド>
→ 返してほしいもの: <具体的に>

### バッチ2: 判断（PC 不要・会話だけで答えられる）
GATE-zz: <選択肢A> か <選択肢B> か。推奨は <>、理由は <>。

## これが解けると次に動くもの
GATE-xx → 手C・手B
EOF
```

書いたらそれもコミットして push する。ドライバがこれを見てループを止める。

## 禁止（`LOOP_PROTOCOL.md` の再掲・最重要のみ）

- **数値を推測で埋めない。** 実データが要るなら GATE にする。**分からないことを分かったように書かない**
- **事前登録の判定規則を、結果を見てから変えない**（変えるのは次の事前登録から）
- **門①（対照）を省かない**
- **集計としてしか言えないものをサイトごとに述べない**（旗81：一致確率 83%）
- **★選択バイアス（旗85）を忘れて率を語らない**（無作為 25〜29% 対 取得群 42〜50%）
- **一次文献を確認したことにしない**（旗92：この環境では構造的に不可能）
- `MEASURED_ONLY_SPINE.md` / `PREMISE_AUDIT.md` / `FLAGS_LOG.md` の**過去の記録を書き換えない**（追記のみ）
- **main に push しない**
- **PR を作らない**
