# 先行研究が示す大課題3つ → 本研究の解決3つ（日本語・コピペ可・厳密引用つき）

論文・スライドにそのまま貼れる形。引用は標準形式（著者・年・誌・巻・DOI）。原文（英語）には
自然な和訳を添えた。**図の所在は実PDFで確認済み**（⚠ は無し）。

---

## 大課題 1：相関では因果はわからない — 観測データから「共通原因を除いた因果」を推定する必要がある

近年、地球システム科学では「相関の先へ進み、観測データから**因果**を推定する」動きが立ち上がった。
その出発点は、**「変数が一緒に動くとき、それは互いに因果があるか、あるいは共通の原因に一緒に
動かされているだけか」**という、ライヘンバッハの共通原因原理である。つまり、放射のような共通の
駆動要因があると、見かけの連動（疑似相関）が本当の因果を覆い隠してしまう。

- **引用**: Runge, J. et al. (2019). Inferring causation from time series in Earth system sciences.
  *Nature Communications*, 10, 2553. DOI:10.1038/s41467-019-10105-3
- **原文**: “correlation does not imply causation … the key idea shared by several approaches follows
  Reichenbach’s common cause principle: if variables are dependent then they are either causal to each
  other … or driven by a common driver.”
  （和訳：相関は因果を意味しない…諸手法に共通する鍵はライヘンバッハの共通原因原理である。変数が
  依存するなら、互いに因果があるか、共通の駆動要因に動かされているかのいずれかだ。）
- **課題を示す図**: **Runge et al. 2019, 図4「Methodological challenges for causal discovery（因果探索の
  方法的な難しさ）」**——自己相関や共通駆動など、相関を誤らせる要因を図解（補助として図3「一般的な
  問題」）。※ 根の一次源はライヘンバッハ Reichenbach (1956) *The Direction of Time*。

## 大課題 2：生態系（生物圏–大気）の過程どうしの「因果と高次の相互作用」がまだよくわかっていない

フラックス観測で変数が同時に動くことは分かっても、**どの要因がどの過程を駆動し、それらが
どう組み合わさっているか**（＝因果と相互作用の構造）は未解明である。とくに、2つずつ（ペアワイズ）の
相関では捉えられない、**3つ以上が組んで初めて現れる高次の相互作用（相乗・冗長）やフィードバック**は
ほとんど特徴づけられていない。

- **引用**:
  - Krich, C. et al. (2020). Estimating causal networks in biosphere–atmosphere interaction with the
    PCMCI approach. *Biogeosciences*, 17, 1033–1061. DOI:10.5194/bg-17-1033-2020
  - Goodwell, A. E., Jiang, P., Ruddell, B. L., Kumar, P. (2020). Debates—Does Information Theory Provide
    a New Paradigm for Earth Science? Causality, Interaction, and Feedback. *Water Resources Research*,
    56, e2019WR024940. DOI:10.1029/2019WR024940
- **原文（Krich 2020）**: “there are still substantial unknowns regarding the exact causal dependencies
  among the different processes.”（和訳：異なる過程どうしの正確な因果依存には、まだ大きな未解明が残る。）
- **課題を示す図**:
  - **Goodwell et al. 2020, 図1「Illustration of different types of causal interactions（因果的相互作用の
    型の図解）」**——プロセスネットワークのノード＝時系列変数、相互作用の型＝相乗・冗長・固有。
  - **Krich et al. 2020, 図4・図5**——渦相関3サイトの**因果ネットワーク**と、その**月ごとの相互作用構造の
    変化**。さらに Krich には「図5を単純相関で描いた版」があり、**相関で描くと偽のつながりが増える**
    （＝相関の限界を示す比較）。

## 大課題 3：陸の炭素吸収の予測の最大の不確実性は「過程応答」にあるが、その相互作用の"形"が観測から検証されていない

将来の陸の炭素吸収がモデル間で最も食い違う原因は、**呼吸・光合成の温度・水分への応答**、すなわち
過程どうしの**組み合わさり方（関数形）の表し方**にある。モデルはこの"形"（たとえば温度と水分は
掛け算で別々に効く、という分離型）を決めて使うが、**それが本当に正しいか（組み合わさって効く＝
相乗ではないか）を、関数形を仮定せず観測から確かめた研究は乏しい**。

- **引用**:
  - Arora, V. K. et al. (2020). Carbon–concentration and carbon–climate feedbacks in CMIP6 models and
    their comparison to CMIP5 models. *Biogeosciences*, 17, 4173–4222. DOI:10.5194/bg-17-4173-2020
  - Booth, B. B. B. et al. (2012). High sensitivity of future global warming to land carbon cycle
    processes. *Environmental Research Letters*, 7, 024002. DOI:10.1088/1748-9326/7/2/024002
- **原文・数値（実PDFで確認済み）**:
  - Arora 2020: “the carbon–climate feedback over land (−45.1 ± 50.6 PgC °C⁻¹) is about 3 times larger
    than over ocean (−17.2 ± 5.0 PgC °C⁻¹).”（和訳：陸の炭素–気候フィードバックは −45.1±50.6 PgC/°C で、
    海の −17.2±5.0 の約3倍。ばらつきも桁違いに大きい。）
  - Booth 2012: “The sensitivity of photosynthetic metabolism to temperature emerges as the most
    important uncertainty.”（和訳：光合成代謝の温度感度が、最大の不確実性として浮かび上がる。）
- **課題を示す図**:
  - **Arora et al. 2020, 図5「陸の炭素–濃度／炭素–気候フィードバック係数（CMIP6 各モデル）」**——陸の
    フィードバックのモデル間ばらつきを示す（海は図6）。
  - **Booth et al. 2012, 図1「2100年までの大気CO₂濃度」**——陸の炭素過程の不確実性で**CO₂予測が大きく
    ばらつく**（最大・最小のメンバーを強調）。
- **正直な留保**: 2100年のCO₂予測**全体**の不確実性は排出シナリオの方が大きい（IPCC AR6）。ただし
  **炭素循環フィードバックのばらつきの中では陸が支配的**、という切り分けで述べる。

---

## 本研究が解決できそうな大きな3つ（各課題に対応）

### 解決 1（大課題1へ）：共通原因を段階的に差し引き、"見かけの連動"から本当の因果の骨組みを取り出す
放射などの共通駆動を、条件付き相互情報量→PCMCI で順に除去する。その結果、
**「放射→蒸発」「光合成→正味CO₂」「気温→地温」といった少数の因果の骨組みが、6サイト・全バイオームで
共通して残った**。相関の限界（大課題1）を、東アジアの実データで具体的に突破した。
（自作図: `fig3_causal_skeleton.png`, `fig4_robustness.png`）

### 解決 2（大課題2へ）：系レベルの高次の相互作用（相乗・冗長）を測り、土地管理での変化を突き止める
O-information を用い、ペアワイズを超えた**系全体の相乗を符号つきで**測る。その結果、
**自然生態系は地下（土壌水分×温度×呼吸）に相乗を持ち、湛水した水田ではそれが崩れ、非湛水の畑では
保たれる**——つまり相乗を壊すのは「農業一般」ではなく「湛水」だと特定した。高次・相互作用・管理下という
空白（大課題2）を突破した。（自作図: `fig6_flooding_mechanism.png`, `fig5_oinfo_synergy.png`）

### 解決 3（大課題3へ）：相互作用の"形"を関数形を仮定せず観測から検出し、モデルの過程仮定を吟味する材料を出す
温度×水分が「分離型（掛け算）」か「相乗型」かを、特定の関数形を仮定せず O-information で判定する。
これは、陸炭素予測の最大の不確実性（大課題3）の核心に、**観測からの相互作用構造の証拠**を与える。
あわせて、乾燥時に光合成が放射から脱結合することも生態系依存として示した。
（自作図: `fig_q10_schematic.png` が問い、`fig6`・`fig1`・`fig7` が答え）

---

## まとめの一文（コピペ可）
> 本研究が属する分野の大きな課題は、①相関では共通原因を除いた因果が読めないこと（Runge et al. 2019）、
> ②生態系の過程どうしの因果と高次の相互作用がまだ未解明であること（Krich et al. 2020, Goodwell et al.
> 2020）、③陸の炭素吸収予測の最大の不確実性である過程応答の"相互作用の形"が観測から検証されていない
> こと（Arora et al. 2020, Booth et al. 2012）、の3つである。本研究はこれらに対し、①共通原因を除いた
> 疎で普遍的な因果の骨組みの抽出、②O-information による系レベルの相乗の測定と「湛水が相乗を壊す」機構の
> 特定、③関数形を仮定しない相互作用構造の観測的検出、という3つで応える。
