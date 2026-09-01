---
name: research-loop
description: この研究リポジトリで自走ループを開始・再開する。ユーザーの手が必要な作業だけを台帳に登録して迂回し、残り全部が人間待ちになった時だけ停止する。コンテキスト消失後の再開にも使う。ユーザーが「ループを回して」「自走して」「続きから再開して」「どこまでやった？」と言ったとき、または /research-loop と打ったときに使う。
---

# research-loop

## まず状態を復元する（会話履歴を信用しない）

```bash
cat docs/LOOP_STATE.md          # 現在地（唯一の再開ポインタ）
cat docs/LOOP_QUEUE.md          # 残作業
cat docs/HUMAN_GATES.md         # 人間待ち一覧
git log --oneline -15
git status --short
```

`SESSION_SUMMARY.md` は 36 KB。**全文を読まないこと。**
`grep -n '^#\{1,3\} ' SESSION_SUMMARY.md` で節の行番号を引き、`sed -n 'X,Yp'` で必要な節だけ読む。
最低限 `§10.3`（禁止 14 項）と `§10.4`（不変条件）は毎回読む。

## 整合性チェック

- `LOOP_STATE.md` の「進行中の作業」が空でない → 前の周が中断されている。まずそれを完結させる
- `LOOP_STATE.md` の「まだファイルに書いていない情報」が空でない → 該当ファイルへ移してコミット
- `git status` が汚れている → `git diff` を読み、完結させるか捨てるか判断
- `LOOP_QUEUE.md` に `DOING` が残っている → それが中断された作業

## 実行

`prompts/02_continuous_loop.md` の `---` 以降に完全な手順がある。**それを読んで従うこと。**

要点だけ再掲:

1. **最初のブロッカーで止まらない。** 人間が要る作業は `HUMAN_GATES.md` に `GATE-xx` として登録し、`LOOP_QUEUE.md` を `BLOCKED(GATE-xx)` にして**次へ進む**
2. **停止するのは残り全部が BLOCKED になった時だけ**
3. **1 周 = 1 コミット + push。** `LOOP_STATE.md` の更新とコミットを済ませずに周を終えない
4. **数値を推測で埋めない。** 分からないなら BLOCKED にする（`§10.3` の禁止 14 項も毎周確認）
5. **コンテキスト残量が尽きる前に引き継ぎ書き出し**を行い、`prompts/00_cold_start.md` で再開できる状態にして止まる

## 停止時

`prompts/02_continuous_loop.md` の「停止時の報告フォーマット」に従う。
GATE はバッチごとにまとめ、ユーザーが席に着く回数を最小化すること。
