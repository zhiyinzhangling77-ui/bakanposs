# セッション現在地（再開用の状態記録）

コンテキストが切れても再開できるよう、現在地・結論・穴の対処・保留中の作業を1枚に。
詳細は FLAGS_LOG.md（各旗）/ VALIDITY_SYNTHESIS.md（発見統合）/ LIMITATIONS.md（自己監査）/
LITERATURE_NOTES.md（文献）/ DATA_QUALITY.md（品質）/ synthesis.html（発表, artifact公開済）。
ブランチ：`claude/new-branch-creation-heq0i1`（全作業push済, リモートが正）。実データはユーザーの /mnt/hdd（コンテナ外）。

## 主題
「豊かな因果構造は、貧しく一部は計算で作られた観測にどこまで写り、どこで写らないか」。渦相関11変数の
うち独立測定は8つ、VPD/GER/GEPは派生量（旗32/35）。

## 確定した発見（3層, honest scope）
- **① 背骨**：放射共通駆動の冗長（O-info）。標本21サイトで一貫（日本冷温帯林に偏る＝擬似反復注意）。
  独立測定間リンク Rg→γH/γLE/Ta/Ts が岩盤（旗35）。GEP→NEEは恒等式で因果でない。
- **② 水マスター**：モンゴル半乾燥ステップ（近接3サイト, "乾燥草原一般"でなくこの母集団に限定）で
  θが呼吸(旗31)と蒸発(旗36)を同符号支配＋Bowen反転。θ→γLEは全独立変数で最も堅い。
- **③ 呼吸の未観測駆動**：下記の弧で「本物の~4日土壌プロセス」と確定。

## 未観測駆動の弧（旗25→42）
旗25 未観測の遅い駆動→旗37 約4日(基質と一致だが分割窓とも一致=交絡)→旗38-GOSIF null(弱いテスト)→
旗39 DT/NT判別で一部は分割窓→旗38-真SIF(8日)季節の生物成分は実在4日記憶は未分解→
旗38-段2 near-daily真SIF「4日記憶=速い基質」を棄却→**旗40 チャンバー直接測定で4日記憶=本物の生物物理**。
**結論**：呼吸の~4日記憶は本物の生物物理（分割非依存の直接測定で確認・大陸をまたぐ独立森林の約1/3・中央
e-fold4日=普遍でない）／最近の光合成でない（SIF棄却）／正体は遅い土壌プロセス（深土壌水分・熱慣性・微生物）。

## 自己監査＝穴と対処（LIMITATIONS.md が正）
- ③ SIF棄却=森林n1最悪品質 → **対処済**（良質4forest JP-Fhk等で一貫, 旗38再テスト）。
- ① 呼吸メモリの緩いバッチ基準 → **再集計で正直化**（普遍51/53→独立森林~1/3・~4日, 旗40厳しい基準）。
- ① 背骨/水マスターの擬似反復 → **主張scope限定で対処**（標本母集団に絞る）。
- ⑤ 多重比較 → **FDR対処**（旗41: 旗13/15の6検定BH, q=0.05で0/6・q=0.10で2/6=状態依存は周辺的, ヌルは頑健）。
- ⑥ 共通w' → Bowen反転で頑健。
- **残り**：② TKY同一サイト（外部連絡=Kishimoto-Mo氏）／④ 呼吸相乗・Q10のGER依存／⑤ ブートCI多サイト比較の族補正／⑦ 夏のみ短記録等。

## ★保留中の作業（次にやること）
- **旗42（穴④を叩く）：水分依存Q10をCOSOREチャンバーで直接検証** — `cosore_q10_step42.py` 作成・合成検証済・push済。
  **ユーザーがローカルで実行し出力待ち**：
  `python research/cosore_q10_step42.py --cosore-dir /mnt/hdd/cosore-0.7.0`（森林）／`--igbp Grassland`（草原）。
  読み：多数森林で★水分依存Q10なら本物(旗26の分割由来説を否定)、割れれば分割由来/サイト固有だった。
  出力を貼ってもらったら解釈→FLAGS_LOG/VALIDITY_SYNTHESIS/LIMITATIONS(④)に反映。

## 主要ツール（research/）
旗25/37 memory_timescale, 旗31 moisture_control_atlas, 旗36 evaporation_regime, 旗32 derivation_audit,
旗35 link_provenance, 旗38 sif_respiration + sif_extract_{geotiff,netcdf,ungridded} + sif_coords,
旗39 memory_partition, 旗40 cosore_memory(＋--cosore-dir バッチ), 旗41 fdr_correction, 旗42 cosore_q10。
データ済：ユーザー /mnt/hdd に JAPANFLUX・cosore-0.7.0・TROPOMI_grid(2018-2021 gridded 8day)。site_coords.csv 生成済。

## 実務メモ（発表/論文で使うなら）
- TROPOMI SIF → Caltech(Koehler/Frankenberg)へ連絡（利用ポリシー）。
- COSORE → 各データセット contributor にクレジット＋Bond-Lamberty et al. 2020 GCB 引用。
- 捏造しない・過大主張しない・撤回は正直に、を通した（複数の自己修正が信頼性の実績）。
