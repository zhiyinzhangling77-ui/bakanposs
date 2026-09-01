# コールドスタート — 記憶ゼロからの再開

> **使う場面**: コンテキストが上限に達した／セッションが切れた／別のセッションで続きをやる。
> 新しいセッションを開いて、以下の `---` 内をそのまま貼る。**これ以外に何も説明しなくてよい設計**。

---

あなたはこのリポジトリで進行中の研究（`japanflux_pn` / 旗方式）を引き継ぎます。
**あなたには一切の会話履歴がありません。** 状態はすべてファイルに書かれています。

## Step 1 — 状態の復元（この順番に意味がある）

```bash
cat research/SESSION_STATE.md            # 1. 現在地（唯一の再開ポインタ）
cat research/HUMAN_GATES.md              # 2. 何が人間待ちか
cat research/OPEN_QUESTIONS_OPTIONS.md   # 3. 残っている手（A〜F）
cat research/LOOP_PROTOCOL.md            # 4. ループ規約と禁止事項
git log --oneline -20                    # 5. 直近の実作業
```

**`research/FLAGS_LOG.md` は 4000 行超あります。最初から全部読まないこと**（それ自体がコンテキストを食う）。
必要な旗だけを引く:

```bash
grep -n '^## 旗' research/FLAGS_LOG.md | tail -20   # 旗の索引と行番号
sed -n 'X,Yp' research/FLAGS_LOG.md                 # 必要な旗だけ
```

`MEASURED_ONLY_SPINE.md`（骨格の正）と `PREMISE_AUDIT.md`（前提の穴）は、
**主張に触れる作業をするときだけ**読めば足ります。

## Step 2 — 整合性チェック（引き継ぎ事故の検出）

読み終えたら以下を確認し、**食い違いがあれば作業を始める前に報告**してください:

- [ ] **`SESSION_STATE.md` が扱う最新の旗番号と、`FLAGS_LOG.md` の最後の旗番号が一致するか**
      → **ずれていたら、その差分が失われた作業である。** `git log` と `FLAGS_LOG.md` から復元し、
        `SESSION_STATE.md` を先に更新してから次へ進む
      （**実例**: 本書作成時点で SESSION_STATE は旗99・FLAGS_LOG は旗102 でずれていた）
- [ ] `SESSION_STATE.md` に「進行中」と書かれた旗があるか → あればそれを先に閉じる
- [ ] `git status` がクリーンか → 汚れていれば `git diff` を読み、完結させるか捨てるか判断
- [ ] 事前登録だけあって結果が無い `PREREGISTRATION_stepNN.md` はないか
      → あればその旗が中断されている

## Step 3 — 再開

整合性が取れたら `research/LOOP_PROTOCOL.md` の「## ループ本体」以降に従ってください。
自走させるなら `/loop` にその内容を貼る。

## やってはいけないこと

- **私に「何をしていましたか？」と聞かない。** 答えは `SESSION_STATE.md` にあります。
  聞かれること自体が引き継ぎ設計の失敗です
- **記憶がないことを理由に最初からやり直さない。** `FLAGS_LOG.md` の旗は完了済みです
- **中断された作業を「たぶんこうだった」で補完しない。** `git diff` で現物を確認すること
- **記憶が無い状態で `MEASURED_ONLY_SPINE.md` や `PREMISE_AUDIT.md` を書き換えない。**
  骨格と前提は、旗を積んで得られたものです
