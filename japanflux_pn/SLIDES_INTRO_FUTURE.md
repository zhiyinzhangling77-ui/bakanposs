# 発表スライド：冒頭2枚（背景・新規性・課題）＋末尾2枚（未達・計画）

先行研究（`LINEAGE.md` の原典）に実引用を紐づけた、スライド用の文面。各スライド＝
【スライドに載せる文（簡潔）】＋【話す（口頭補足）】＋【図】＋【引用】。捏造なし・正直な限界つき。

---

## スライド 1 —— 背景・重要性・課題
**タイトル案**: 「生態系の炭素・水・エネルギー交換：相関の先の"相互作用の構造"を読む」

【スライドに載せる文】
- **重要性（数値つき）**: 陸は人為 CO₂ の**約 1/4〜1/3 を吸収**（Global Carbon Budget 2024: 陸シンク
  2.3±1.0 vs 総排出 11.1±0.9 GtC/yr）。だが**"将来どれだけ吸うか"はモデルで大きく食い違う**——
  CMIP6 の**陸の炭素–気候フィードバックは −45.1 ± 50.6 PgC/°C**（Arora et al. 2020）で、
  **海（−17.2 ± 5.0）の約 3 倍・ばらつき約 10 倍**（誤差がゼロをまたぐ）。IPCC AR6 第5章も
  「**陸の炭素貯留の変化が炭素循環予測の鍵となる不確実性**」とする。
- **その不確実性の"原因"＝本研究が測る結合**: ①**土壌呼吸の温度感度**（Booth et al. 2012:「光合成代謝の
  温度依存が最大の不確実性」）、②**水ストレス下の光合成**（Humphrey et al. 2021:「**土壌水分が陸シンクの
  経年変動の約 90% を駆動**」）、③**年々変動は温度か水か**（Jung et al. 2017）。
- **正直な留保**: 2100 年 CO₂ 予測**全体**の不確実性は**排出シナリオの方が大きい**（AR6）。ただし
  **炭素循環フィードバックのばらつきは陸が支配的**、という切り分けで話す。
- **課題1（疑似相関・交絡）**: 変数は互いに動くが、**放射など共通原因**が"見かけの連動"を作り、本当の
  因果を覆う（**Reichenbach の共通原因原理**）。主流のフラックス解析（相関・回帰・機械学習）は
  **交絡の明示的制御が不十分**なことが多い。
- **課題2（相互作用の"形"が未検証）**: 温度と水分は**独立に効くのか、組み合わさって効くのか（相乗）**——
  その相互作用の"**形**"を、**関数形を仮定せず観測から検証する枠組みが乏しい**。
  （相互作用モデルは存在するが特定の関数形を仮定し、主に点スケールで検証）。
- **課題3（対象の偏り）**: 情報理論プロセスネットワークの先行は**主に北米のサイト中心**
  （**Ruddell & Kumar 2009 は米国イリノイ州 Bondville のトウモロコシ–大豆畑**、Goodwell & Kumar も米国⚠）。
  **湿潤モンスーンの東アジア、とくに湛水水田は、この枠組みで扱われていない**。

【話す（自分の言葉の台本）】「陸は今、人間が出す CO₂ の約 1/4〜1/3 を吸っています。でも**将来どれだけ
吸うかはモデルで大きく食い違う**。CMIP6 では陸のフィードバックが −45±51 PgC/°C で、**誤差がゼロを
またぐほど大きい**——海（−17±5）の約 3 倍・ばらつき約 10 倍。IPCC も『陸の炭素貯留が炭素循環予測の
鍵となる不確実性』としています。正直に言うと、2100 年 CO₂ 全体の不確実性は排出シナリオの方が大きい。
それでも、**炭素循環フィードバックのばらつきは陸が支配的**。そしてその原因が、①土壌呼吸の温度感度、
②水ストレス下の光合成、③年々変動が温度か水か——**まさに本研究が情報理論で測る"駆動変数とフラックスの
結合"そのもの**なんです。」
【図】**`fig_uncertainty.png`（このスライドの主役）**＝左: 陸 vs 海のフィードバック（陸が大きく・不確実）、
右: 不確実性の3原因＝本研究が測る結合。（導入の概念図 `fig0_concept_network.png` はスライド2冒頭でも可）
【引用】[1] Reichenbach 1956（共通原因）／[2] Runge et al. 2019（交絡・因果の必要性）／
[3] Ruddell & Kumar 2009（プロセスネットワーク）／[4] Goodwell & Kumar 2017（相乗・冗長）。

---

## スライド 2 —— アプローチ・新規性
**タイトル案**: 「情報理論×因果推論で"相互作用の構造"を測り、既存モデルの仮定を吟味する」

【スライドに載せる文】
- **アプローチ**: 相関でなく、**情報理論**（Transfer Entropy／**O-information**：相乗/冗長を符号付きで）
  ＋**因果推論**（**PCMCI**：共通原因・時間遅れを除いた多変量因果）。
- **系譜（本研究はこの"次の一歩"）**: プロセスネットワーク[3] → PID/相乗[4] → 多変量因果[5,6] →
  情報でモデル診断[7]。
- **新規性（正直に3点）**:
  1. 手法・概念は既存。貢献は**新公開 JapanFlux2024＋ChinaFlux＋KoFlux（日中韓）・湿潤モンスーン・
     管理生態系（水田）への系統的拡張**。
  2. **O-information で"湛水が地下の高次相乗を選択的に壊す"機構を特定**（日韓2水田 vs 中国畑の対照）＝
     自然サイト中心の先行では原理的に出せない一点。
  3. **関数形を仮定せず**相互作用構造を**モデルフリーに検出**（DAMM 等は特定の形を仮定）。
- **すでに得た核（1行ずつ）**: 因果骨格は疎で普遍／湛水で地下相乗が崩壊（日韓再現）／光利用の放射脱結合は
  水律速の森林で。

【話す】情報理論の強みは「相互作用の"質"（相乗か冗長か）を符号付きで測れる」こと。O-information の定義は
**Rosas 2019 の式に厳密準拠**。これを東アジア・水田へ広げ、モデルの決めつけを観測で吟味する。
【図】`fig_pipeline.png`（手法の流れ）＋（結果の予告として）`fig6_flooding_mechanism.png` を小さく。
【引用】[3] Ruddell & Kumar 2009／[4] Goodwell & Kumar 2017／[5] Runge et al. 2019（PCMCI）／
[6] Krich et al. 2020（PCMCI×生物圏–大気）／[7] Nearing et al. 2020（情報でモデル診断）／
[8] Rosas et al. 2019（O-information の定義）。

---

## スライド N−1 —— まだできていないこと（正直に）
**タイトル案**: 「今の限界」

【スライドに載せる文（シンプルに4つ）】
- **サイトがまだ少ない**（とくに水田・畑）。
- **「湛水で土壌水分の変動が消える」ことを、まだ直接は確かめていない**（分散を比べていない）。
- **季節は夏だけ**（春・秋・昼夜はこれから）。
- **モデルとはまだ比べていない**（今は観測だけ）。

【話す】機構は見えたが、標本・直接確認・季節・モデル比較が残っている。ここを一つずつ埋めたい。
（＝断定でなく「まだできていない・これから確かめたい」という姿勢で話す）
【図】`fig5b_oinfo_ci.png`（標本が少ないと誤差が大きい、を正直に示す図）。

---

## スライド N —— これからやってみたいこと（学部 → 大学院）
**タイトル案**: 「今後やってみたいこと」

【スライドに載せる文（シンプル・"やってみたい"の姿勢）】
- **まず手元でできる小さな確認から**: 水田と畑で**土壌水分の変動の大きさを比べてみたい**
  （湛水で本当に変動が消えているか）。
- **サイトを少し増やしたい**: 水田・畑をもう数か所（KoFlux ほか）。
- **既存モデルと比べてみたい**: 呼吸・光合成の式（温度×水分・VPD）が、観測の相互作用を再現できるか。
- **気候の幅を広げたい**: 乾燥・大陸のサイトも入れて、湿潤日本の外へ。

【使いたいデータ】いまの日本＋中国＋韓国の 30 分フラックス（＋将来はモデル出力）。
【使いたい手法】いまの情報理論（O-information など）＋因果推論（PCMCI）を、モデルとの比較にも広げる。

【話す】いきなり大きいことでなく、**まず"土壌水分の変動を比べる"という小さな確認**から始めたい。
そこからサイトを増やし、モデルと比べる、という順で進めたい。
（＝まだ着手していないので「やるべき」でなく「やってみたい・こう進めたい」と言う）
【図】（任意）簡単なロードマップ図（学部＝確認と基盤 / 大学院＝モデル比較と気候拡張）。作成可。

---

## 引用番号（スライド脚注用・詳細は LINEAGE.md）
[1] Reichenbach 1956 *The Direction of Time*（共通原因原理）
[2] Runge et al. 2019 *Nat. Commun.*（地球システムの因果推論）
[3] Ruddell & Kumar 2009 *WRR*（エコ水文プロセスネットワーク）
[4] Goodwell & Kumar 2017 *WRR*（時間的情報分割：相乗/固有/冗長・TIPNets）
[5] Runge et al. 2019 *Sci. Adv.*（PCMCI）
[6] Krich et al. 2020 *Biogeosciences*（PCMCI×生物圏–大気）
[7] Nearing et al. 2020 *WRR*（情報理論でモデル診断）
[8] Rosas et al. 2019 *Phys. Rev. E* 100, 032305（O-information）
（背景の呼吸・気孔・モデル評価: DAMM=Davidson 2012／Medlyn 2011／Jung 2020 FLUXCOM は
`MODEL_COMPARISON_LIT.md`・`LINEAGE.md §K` 参照。⚠ p./原文は原著で確定）

**スライド1「重要性・不確実性」の裏づけ（Perplexity 探索で確認, ✅ DOI取得済 / ⚠図表番号・原文は原著で最終確認）:**
- ✅ **Arora, V. K. et al. 2020**. Carbon–concentration and carbon–climate feedbacks in CMIP6 models…
  *Biogeosciences* 17, 4173–4222. DOI:10.5194/bg-17-4173-2020。**引用値**: 陸の炭素–気候フィードバック
  **−45.1 ± 50.6 PgC/°C**、海 **−17.2 ± 5.0**（陸≈3倍・spread≈10倍）。図: **Fig. 5** 系（⚠番号確認）。
- ✅ **Friedlingstein, P. et al. 2006**. Climate–Carbon Cycle Feedback Analysis: C4MIP. *J. Climate* 19, 3337–3353.
  DOI:10.1175/JCLI3800.1。原文: "…there is still a large uncertainty on the magnitude of these sensitivities."
- ✅ **Friedlingstein, P. et al. 2014**. Uncertainties in CMIP5 projections due to carbon cycle feedbacks.
  *J. Climate* 27, 511–526. DOI:10.1175/JCLI-D-13-00177.1。**引用値**: RCP8.5 で 2100 年 CO₂ が **795–1145 ppm**。
- ✅ **Booth, B. B. B. et al. 2012**. High sensitivity of future global warming to land carbon cycle processes.
  *ERL* 7, 024002. DOI:10.1088/1748-9326/7/2/024002。原文: "…temperature dependences of photosynthetic
  metabolism represents the most important uncertainty identified."
- ✅ **Cox, P. M. et al. 2000**. Acceleration of global warming due to carbon-cycle feedbacks. *Nature* 408,
  184–187. DOI:10.1038/35041539。
- ✅ **Humphrey, V. et al. 2021**. Soil moisture–atmosphere feedback dominates land carbon uptake variability.
  *Nature* 592, 65–69. DOI:10.1038/s41586-021-03325-5。原文: "soil moisture drives 90% of the inter-annual
  variability in global land carbon uptake"。
- ✅ **Humphrey, V. et al. 2018**. *Nature* 560, 628–631. DOI:10.1038/s41586-018-0424-4（TWS と CO₂ 成長率）。
- ✅ **Jung, M. et al. 2017**. Compensatory water effects link yearly global land CO₂ sink changes to temperature.
  *Nature* 541, 516–520. DOI:10.1038/nature20780。
- ✅ **Global Carbon Budget 2024** = Friedlingstein, P. et al. 2025. *ESSD* 17, 965–1039. DOI:10.5194/essd-17-965-2025。
  **引用値**: 2023 陸シンク **2.3 ± 1.0**、総排出 **11.1 ± 0.9 GtC/yr**（≈1/4〜1/3）。
- ⚠ **IPCC AR6 WG1 Ch.5**（Canadell et al. 2021）: "Changes in land carbon storage remain the key uncertainty
  in carbon cycle projections."／"Uncertainty in atmospheric CO₂ by 2100 is dominated by emissions scenarios
  rather than … carbon–climate feedbacks."（節・図番号は AR6 PDF で確認）。
- 補足: Baldocchi, Chu & Reichstein 2018 *Agric. For. Meteorol.* 249, 520–533（フラックスの経年変動）／
  Jung et al. 2020 *Biogeosciences* 17, 1343–1365（FLUXCOM: IAV 過小評価）DOI:10.5194/bg-17-1343-2020。
**課題3の裏づけ（先行のサイト）:**
- Ruddell & Kumar 2009 *WRR*: 米国イリノイ州 **Bondville**（トウモロコシ–大豆畑, US Corn Belt, 2003-07）
  ＝本文で確認済み。Goodwell & Kumar 2017/2018 のサイトは ⚠ 原著で確認（米国サイトの見込み）。
