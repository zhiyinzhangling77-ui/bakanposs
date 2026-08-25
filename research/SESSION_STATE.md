# セッション現在地（再開用の状態記録）

コンテキストが切れても再開できるよう、現在地・結論・穴の対処・保留中の作業を1枚に。
詳細は FLAGS_LOG.md（各旗）/ MEASURED_ONLY_SPINE.md（測定量だけの骨格）/ PREMISE_AUDIT.md（前提の監査）/
LIMITATIONS.md（穴と対処）/ PREREGISTRATION_step55.md・step57.md（事前登録）/ LITERATURE_NOTES.md /
CONTACT_DRAFT_TKY.md（外部連絡の草案）/ synthesis.html（発表, artifact 公開済）。
ブランチ：`claude/new-branch-creation-heq0i1`（全作業push済, リモートが正）。実データはユーザーの /mnt/hdd（コンテナ外）。

## 主題
「豊かな因果構造は、貧しく一部は計算で作られた観測にどこまで写り、どこで写らないか」。

## いまの骨格（`MEASURED_ONLY_SPINE.md` が正）

**A. 測定量だけで立つ（骨格の中核）**
1. チャンバー呼吸の多日メモリ：**22/45**、候補は **20/22 で説明せず**（旗40→53/54で較正し直し）
2. チャンバー呼吸の水分依存Q10：**18/36**（旗42→44で温度交絡を制御）。符号の規則は未解明（旗48/55）
3. θ→蒸発の支配と Bowen 反転：派生量ゼロ。ただし統計的重みは**独立地点 n=1**（旗43）

**B. 計算量に依存したまま（明記して述べる）**：タワー側メモリ／SIF棄却／呼吸の相乗／θ→GER
　※旗56 でタワー側を測定量だけ（夜間実測NEE）に移そうとしたが**検定が成立せず**
　（8中7が駆動弱 R²=0.04〜0.22）＝**このデータでは原理的に動かせない**と確定

**C. 恒等式・健全性チェック**：エネルギー系の冗長（旧「背骨」）＝収支 H+LE≈Rn−G の再記述（旗50, 78/87）

## 前提の穴（PREMISE_AUDIT が正）— in-container で叩けるものは全て対処済み
①Q10の温度交絡→旗44で部分的中（2/3→約半数）／②因果十分性→言い直し＋旗47(LPCMCI, **実データ未実行**)／
③背骨は恒等式→**旗50で的中・確定**／④SIFは弱い操作変数→言い直し／⑤gap-fill→旗46（背骨は無傷・**炭素コアは的中**）／
⑥out-of-sample→旗48（自分の後付けが落第）＋**旗55で初の事前登録・示唆を規則通り取り下げ**／
⑦新規性→方法論監査に位置づけ（Q10はDAMMの再発見の可能性）／
**⑨検出器のモデル形（旗52で新発見）**→非線形系への線形当てはめ自体が自己相関残差を作る。旗53/54で制御しても記憶は残存

## ★保留中：旗57（実行待ち）
**中心的主張の事前登録レプリケーション**——メモリ解析は森林でしかやっておらず、**非森林サイトは未使用**。
予測は `PREREGISTRATION_step57.md` に実行前確定・commit 済（1d676ce）：
H1 普及率>0.25（森林0.49）／H2 説明されなさ>0.6（森林0.91）／判定可能<6なら判定しない。
```
python research/replicate_nonforest_step57.py --cosore-dir /mnt/hdd/cosore-0.7.0
```
結果は FLAGS_LOG に追記し、規則通りにスコープを決める（★再現／○部分再現／▲森林に限定）。

## 残る穴（in-container では叩けない）
- **② TKY 同一サイト**：旗51 で「手元のチャンバーとタワーは10km以内に一組も無い」と確認＝外部依存で確定。
  連絡文面の草案は `CONTACT_DRAFT_TKY.md`（**送信はユーザー判断**。送る前に公開済みでないか等の確認手順つき）
- **④残り 呼吸の相乗の GER 依存**：チャンバーはドライバー2変数で3変数以上の相乗が測れない＝**構造的に不可**
- **⑦ 標本の構造的偏り**：夏(7-8月)のみ・短記録・θ深度不統一

## 主要ツール（research/）
旗25/37 memory_timescale・旗31 moisture_control_atlas・旗36 evaporation_regime・旗32 derivation_audit・
旗35 link_provenance・旗38 sif_respiration＋sif_extract_*・旗39 memory_partition・旗40 cosore_memory・
旗41 fdr_correction・旗42 cosore_q10・旗43 fdr_multisite・旗44 q10_confound・旗45 memory_attribution・
旗46 backbone_gapfill・旗47 latent_confounder(LPCMCI)・旗48 out_of_sample・旗50 energy_identity・
旗51 colocate・**旗52 synthetic_tower（合成タワー＝全体較正）**・旗53 chamber_memory_recount・
旗54 memory_attribution_flex・旗55 thermal_depth・旗56 nightly_nee_memory・旗57 replicate_nonforest。
データ：/mnt/hdd に JAPANFLUX・cosore-0.7.0・TROPOMI_grid。

## この研究の作法（維持すること）
- 合成で検出器を検証してから実データに当てる。**帰無条件を必ず作る**。プラセボ（対照）を併走させる。
- 後付けの説明は**事前登録してから**検定する（旗55/57）。外れたら規則通りに取り下げる。
- 撤回・格下げは記録に残す。**自分の道具の欠陥も記録する**（旗52-56 で計6件見つけて直した）。
- 捏造しない。未確認は「未確認」と書く（例：doi.org 遮断で一次確認できなかった文献は引用しない）。
