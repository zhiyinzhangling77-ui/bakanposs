# 既存研究が示す大課題3つ（厳密引用・課題を示す図・僕の解決）

論文の導入に使える形。各課題＝**厳密な引用（著者・年・誌・巻・DOI）＋原文（一次資料から抽出）
＋その課題を"示す図"の場所（論文・図番号・何を示すか）**。⚠ は原著で最終確認（未入手PDF）。
図番号は**手持ちPDFで実キャプションを確認済**のものはその旨明記。

---

## 大課題 1 —— 相関では因果がわからない：観測時系列から"共通原因を除いた因果"を推定する必要
**引用**: Runge, J., Bathiany, S., Bollt, E. et al. (2019). *Inferring causation from time series in Earth
system sciences.* **Nature Communications** 10, 2553. DOI:10.1038/s41467-019-10105-3。
**原文（本文より抽出）**: "Although the truism that correlation does not imply causation holds, the key
idea shared by several approaches follows **Reichenbach's common cause principle**: if variables are
dependent then they are either causal to each other … or **driven by a common driver**."
**課題を示す図**:
- **Runge 2019, Fig. 4「Methodological challenges for causal discovery in complex spatio-temporal
  systems」**（自己相関・共通駆動など、相関を誤らせる要因を図解）＝**手持ちPDFでキャプション確認済**。
- 補助: **Fig. 3「Key generic problems in Earth system sciences」**（同、一般的問題の図）。
> 根っこの一次源: **Reichenbach, H. (1956). The Direction of Time.**（共通原因原理）＝Runge の引用元。

## 大課題 2 —— 生態系（生物圏–大気）の過程どうしの"因果・高次の相互作用"が未解明（ペアワイズ相関では届かない）
**引用**:
- Krich, C., Runge, J., Miralles, D.G. et al. (2020). *Estimating causal networks in biosphere–atmosphere
  interaction with the PCMCI approach.* **Biogeosciences** 17, 1033–1061. DOI:10.5194/bg-17-1033-2020。
- Goodwell, A.E., Jiang, P., Ruddell, B.L., Kumar, P. (2020). *Debates—Does Information Theory Provide a
  New Paradigm for Earth Science? Causality, Interaction, and Feedback.* **Water Resources Research** 56,
  e2019WR024940. DOI:10.1029/2019WR024940。
**原文（抽出）**: Krich: "…there are still **substantial unknowns regarding the exact causal dependencies
among the different processes**." Goodwell: 因果的相互作用（相乗・冗長・フィードバック）を情報理論で
特徴づける枠組みを提示（Abstract 参照）。
**課題を示す図**:
- **Goodwell 2020, Fig. 1「Illustration of different types of causal interactions」**（プロセスネット
  ワークのノード＝時系列変数、相互作用の型＝相乗/冗長/固有）＝**手持ちPDFでキャプション確認済**。
- **Krich 2020, Fig. 4/5**（渦相関3サイトの**因果ネットワーク**と、その**月ごとの相互作用構造の変化**）
  ＝手持ちPDFで確認済。さらに Krich には「**Fig. 5 を単純相関で描いた版**」があり、**相関だと偽リンクが
  増える**＝相関の限界を示す比較図（補足）。

## 大課題 3 —— 陸炭素予測の最大の不確実性は過程応答（呼吸・光合成の温度×水）だが、その相互作用の"形"の妥当性が観測から十分検証されていない
**引用**:
- Arora, V.K., Katavouta, A., Williams, R.G. et al. (2020). *Carbon–concentration and carbon–climate
  feedbacks in CMIP6 models…* **Biogeosciences** 17, 4173–4222. DOI:10.5194/bg-17-4173-2020。
- Booth, B.B.B., Jones, C.D., Collins, M. et al. (2012). *High sensitivity of future global warming to
  land carbon cycle processes.* **Environmental Research Letters** 7, 024002.
  DOI:10.1088/1748-9326/7/2/024002。
**原文/数値（Perplexity 探索, ⚠原著で最終確認）**: Booth: "…**temperature dependences of photosynthetic
metabolism represents the most important uncertainty identified**." Arora: 陸の炭素–気候フィードバック
**−45.1 ± 50.6 PgC/°C** vs 海 **−17.2 ± 5.0**（陸≈3倍・spread≈10倍）。
**課題を示す図（⚠原著で番号確認, 未入手PDF）**:
- **Arora 2020, Fig. 5**（land の β_L・γ_L の**モデル間ばらつき**を示す主要図）⚠。
- Booth 2012（陸炭素過程の仮定で**将来昇温レンジが大きく変わる**図）⚠。
- ※ 自作の要約図 `japanflux_pn/slides/fig_uncertainty.png`（陸 vs 海のフィードバック＋原因）で代替可。

---

## そこから僕が解決できそうなこと（大きく3つ・各課題に対応）

### 解決 1（↔ 大課題 1）: 共通原因を段階的に差し引き、"見かけの連動"から本当の因果骨格を取り出す
- 条件付き相互情報量 → **PCMCI** で放射などの共通駆動を除去。**結果**: 疎で普遍的な因果骨格
  （Rg→γLE, GEP→NEE, Ta→Ts）が**6サイト・全バイオームで安定**（自作図 `fig3_causal_skeleton.png` /
  `fig4_robustness.png`）。→ 相関の限界（課題1）を東アジアの実データで具体的に突破。

### 解決 2（↔ 大課題 2）: 系レベルの高次相互作用（相乗/冗長）を測り、土地管理での変化を特定
- **O-information**（Rosas 2019）で、ペアワイズを超えた**系全体の相乗を符号付き**で測定。**結果**:
  自然生態系は地下（θ×温度×呼吸）に相乗、**湛水水田は崩壊、非湛水畑は保持**＝"湛水が相乗を壊す"
  （自作図 `fig6_flooding_mechanism.png` / `fig5_oinfo_synergy.png`）。→「高次・相互作用・管理下」の
  空白（課題2）を突破。

### 解決 3（↔ 大課題 3）: 相互作用の"形"を関数形フリーで検出し、モデルの過程仮定を観測から吟味する材料を出す
- 分離型（掛け算）か相乗型かを、**特定の関数形を仮定せず** O-info で判定（自作図 `fig_q10_schematic.png`
  が問い、`fig6` が答え）。**気候脱結合**も生態系依存として観測から提示（`fig1`/`fig7`）。→ 陸炭素
  不確実性の核（課題3）に、**観測からの相互作用構造の証拠**を与える。

---

## 一言（このパートの背骨）
> **大課題**＝①相関では共通原因を除いた因果が読めない（Runge 2019）、②生態系の過程の因果・高次相互作用が
> 未解明（Krich 2020・Goodwell 2020）、③陸炭素の最大不確実性＝過程応答の相互作用の"形"が観測未検証
> （Arora 2020・Booth 2012）。**僕の解決**＝①共通駆動を除いた**疎で普遍の因果骨格**、②**O-informationで
> 系レベル相乗と湛水崩壊**、③**関数形フリーで相互作用の形を観測から**。厳密引用と図の所在は上表。
