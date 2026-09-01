---
name: research-loop
description: この研究（japanflux_pn / 旗方式）の自走ループを開始・再開する。人手が要る作業だけを HUMAN_GATES に登録して迂回し、残り全部が人間待ちになった時だけ停止する。コンテキスト消失後の再開にも使う。ユーザーが「ループを回して」「自走して」「続きから再開して」「どこまでやった？」と言ったとき、または /research-loop と打ったときに使う。
---

# research-loop

## まず状態を復元する（会話履歴を信用しない）

```bash
cat research/SESSION_STATE.md            # 現在地（唯一の再開ポインタ）
cat research/HUMAN_GATES.md              # 人手が要る一覧
cat research/OPEN_QUESTIONS_OPTIONS.md   # 残っている手（A〜F）
git log --oneline -20
git status --short
```

`research/FLAGS_LOG.md` は 4000 行超。**全文を読まない。**
`grep -n '^## 旗' research/FLAGS_LOG.md | tail -20` で索引を引き、`sed -n 'X,Yp'` で必要な旗だけ読む。

## 整合性チェック

- **`SESSION_STATE.md` の最新旗 ≠ `FLAGS_LOG.md` の最終旗** → その差分が失われた作業。
  `git log` と `FLAGS_LOG.md` から復元し、SESSION_STATE を先に更新する
- 事前登録だけあって結果が無い `PREREGISTRATION_stepNN.md` → その旗が中断されている
- `git status` が汚れている → `git diff` を読み、完結させるか捨てるか判断

## 実行

`research/LOOP_PROTOCOL.md` の「## ループ本体」以降に完全な手順がある。**それを読んで従うこと。**

要点だけ再掲:

1. **最初のブロッカーで止まらない。** 人手が要る作業は `HUMAN_GATES.md` に `GATE-xx` として登録し、
   `BLOCKED` にして**手元でできる次へ進む**
2. **停止するのは残り全部が BLOCKED になった時だけ。** 停止時は GATE をバッチにまとめる
3. **1 旗 = 1 コミット + push。`SESSION_STATE.md` を更新せずに旗を閉じない**
4. **数値を推測で埋めない。** 分からないなら GATE にする
5. **事前登録の判定規則を結果を見てから変えない**（変えるのは次の事前登録から）
6. **門①（対照）を省かない。** 集計としてしか言えないものをサイトごとに述べない
7. **コンテキスト残量が尽きる前に引き継ぎ書き出し**をして、`research/COLD_START.md` で再開できる状態にする

台帳を洗い直すときは `research/GATE_AUDIT_PROMPT.md` を使う。
