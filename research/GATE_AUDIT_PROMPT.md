# 人手が要る箇所の洗い出しプロンプト（ときどき回して台帳を更新する）

> 使い方: 以下の `---` 内をそのまま貼る。`research/HUMAN_GATES.md` が差分更新される。

---

この研究を最後まで完走するために、**私（人間）の手が必要になる箇所を、洗いざらい・網羅的に**列挙してください。

## 前提
- あなたはコンテナ内で動いており、実データ `/mnt/hdd/...` に**到達できません**（`SESSION_STATE.md` 冒頭）
- 一次文献にも届きません（旗92・計 8 ホスト遮断、`LITERATURE_VERIFICATION_TODO.md`）
- 外部サービス（AmeriFlux・PhenoCam・COSORE 等）への認証・DUA 同意も持ちません

## やること

**Step 1 — 棚卸し.** `OPEN_QUESTIONS_OPTIONS.md`（手 A〜F）・`japanflux_pn/RUN_HEAVY.md`・
`NEW_OBSERVATION_DESIGN.md`・`LITERATURE_VERIFICATION_TODO.md`・`CONTACT_DRAFT_TKY.md`・
`FUTURE_PLAN.md`（Phase 1〜5）・直近の旗の「次に効く一手」を起点に、
**卒論の完成まで**に必要な作業をすべて列挙する。

**Step 2 — AUTO / GATE に二分する.**
`AUTO` = コンテナ内で完結（合成検証・道具作り・文章・記録の整理）。`GATE` = 人手が要る。

**Step 3 — GATE を原因で分類する.**
`D` 実行（ローカルの実データ）／`A` 取得（アカウント・DUA）／`P` 外界（人に聞く・現地）／
`J` 判断（研究上の決定でありAIが代行してはならない）／`S` 仕様（締切・書式など外から与えられる）。

**Step 4 — 各 GATE に `LOOP_PROTOCOL.md` の 9 項目を埋める.**

**Step 5 — 依存グラフ.** GATE 間・GATE→AUTO の依存を図示し、
**「これ 1 つ解けば最も多くが動き出す」順に上位 3 つ**を明示する。

**Step 6 — バッチ化.** 私が 1 回 PC に向かえばまとめて片付く GATE をグループ化し、
**私が席に着く回数を最小化**する構成にする。夜間放置できるものは明示する。

## 出力
`research/HUMAN_GATES.md` を**差分更新**（既存 ID は維持、新規は連番、解消済みは `status: resolved`）してコミット。

## 禁止
- `LOOP_PROTOCOL.md` の「絶対にやってはいけないこと」に触れる提案をしない
- **GATE を減らすために結果を推測で埋めない**。分からないものは GATE のまま残す
- 「たぶん手元で動くはず」で AUTO に分類しない。**実データが要るなら必ず D**
