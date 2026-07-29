# 学習ガイド — この研究を 100% 理解するための知識マップとプロンプト設計

本研究（JapanFlux2024 への Ruddell & Kumar 型プロセスネットワーク解析 + 共通駆動条件付け
+ PCMCI 因果ネットワーク）を根本から理解するために必要な知識を、**依存順**に並べ、各項目を
**あなたのコード/判断のどこに現れるか**に紐づけた。後半に、教科書データベースを持つ AI から
**オフライン保存できる深い解説を効率よく引き出すためのプロンプト設計**を置く。

---

## Part 1. 知識マップ（依存順の学習シラバス）

下の Tier は前提が下から積み上がる。Tier 1→7 の順で学ぶと手戻りが少ない。各項目末尾の
`⟶` は本研究での出現箇所。

### Tier 0 — 数学の土台
- 確率変数、同時分布 p(x,y)、周辺分布、**条件付き分布** p(y|x)、期待値、独立性
- 対数と情報の単位: **nat（自然対数）** と bit（log2）、`H = −Σ p ln p` ⟶ `information_theory.py` は nat
- 総和・ヒストグラム・度数の基本 ⟶ `np.bincount` による同時分布推定

### Tier 1 — 情報理論コア ⟶ `information_theory.py`
- **Shannon エントロピー** H(X)、同時エントロピー H(X,Y)、条件付きエントロピー H(Y|X)
- **相互情報量** I(X;Y) = H(X)+H(Y)−H(X,Y)。対称・非負・「共有情報」の意味 ⟶ `mutual_information_indices`
- **条件付き相互情報量** I(X;Y|Z) = H(X,Z)+H(Y,Z)−H(Z)−H(X,Y,Z) ⟶ `conditional_mutual_information_indices`
- チェインルール、**データ処理不等式**、KL ダイバージェンス
- 正規化: I' = I / log(m)（%表示の根拠）⟶ `config.log_m`

### Tier 2 — 有限標本でのエントロピー推定 ⟶ m=11 ビン, バイアス floor
- **プラグイン（ヒストグラム）推定**と**正バイアス**（なぜ標本が少ないと I/TE が過大に出るか）
- **次元の呪い**: m=11 の 3 次元同時分布 = 11³=1331 セル、n≈1500 で 1 セル 1 点 ⟶ TE floor, 負の drop
- **Miller-Madow 補正** `+(K−1)/(2N)` ⟶ `_entropy_of_indices(correct=True)`
- **KSG（Kraskov-Stögbauer-Grassberger）/ KNN 推定**（ビン無し、小標本に強い）⟶ 次段階 CMIknn の中身
- カーネル推定、適応ビン、記号（順序）エントロピー

### Tier 3 — 情報流と Transfer Entropy ⟶ `information_theory.py`, `network.py`
- **Transfer Entropy**（Schreiber 2000）: T(X→Y) = I(Y_t ; X_past | Y_past)。方向性・非対称性
- **Knuth 2005 形**（R&K 式5）: `T = H(X_{t−τ},Y_{t−1}) + H(Y_t,Y_{t−1}) − H(Y_{t−1}) − H(X_{t−τ},Y_t,Y_{t−1})` ⟶ `transfer_entropy_indices`
- **TE = 条件付き MI**、Gaussian で**Granger 因果と等価**（Barnett 2009）
- **R&K 2009 プロセスネットワーク**: 11 変数、5 日前方窓アノマリ、固定ビン、結合分類
  （type 1 sync / 2 feedback / 3 forcing / 4 uncoupled）、**Tz = T/I** の比 ⟶ `classify_coupling`, `build_network`
- 特徴ラグ τ'、ラグ走査曲線 T(τ) ⟶ `first_significant_peak`, `te_lag_curve`
- **発展**: Goodwell & Kumar 2017 **TIPNets** と **PID（部分情報分解: 冗長/固有/相乗）**（Williams & Beer 2010）
  — 「なぜ同期支配か（冗長性か相乗か）」を分解する、R&K の直系後継

### Tier 4 — 有意性検定と多重比較 ⟶ サロゲート, `peak_min_run`
- 仮説検定、p 値、**ヌル分布**の考え方
- **サロゲートデータ法**: シャッフル置換（時間結合を壊す）、**層内置換**（(X,Z) を保ち X-Y|Z だけ壊す＝条件独立ヌル）、
  block bootstrap、IAAFT ⟶ `surrogate_te_stats`, `surrogate_cmi_stats`
- しきい値 Δ = μ_ss + c·σ_ss、c=2.36 が α=0.01 片側に対応する理由（正規近似）⟶ `config.sig_c`
- **多重比較問題**: 36 ラグ走査で FWER = 1−(1−α)³⁶ ≈ 30%、Bonferroni / FDR、
  「連続有意帯を要求」する経験的補正 ⟶ `peak_min_run`

### Tier 5 — 時系列の因果推論 ⟶ `condition_driver.py`, `causal_network.py`
- **交絡（confounding）・共通駆動・連鎖・合流点（collider）**、DAG、**d 分離** — なぜペアワイズ相関/MI が偽の結合を生むか
- **条件付け**で共通駆動を除く原理 ⟶ I(X;Y|Rg) の drop 診断
- 多変量/条件付き **Granger 因果**（VAR）
- **PCMCI / PCMCI+**（Runge et al. 2019）: PC 条件選択 + **MCI（Momentary Conditional Independence）**、
  自己相関への頑健性、独立性検定の選択（**ParCorr 線形** vs **CMIknn 非線形/KNN** vs GPDC）⟶ `run_pcmci`
- **仮定**: 因果十分性（未観測交絡なし）、定常性、忠実性（faithfulness）、マルコフ性 — 破れると何が起きるか（例: Rg への逆向き偽リンク）
- **代替**: 収束クロスマッピング **CCM**（Sugihara 2012, Takens 埋め込み）— TE/Granger と因果の定義が違う裏取り

### Tier 6 — 時系列解析の基礎 ⟶ アノマリ, 定常性
- 定常性、自己相関関数、スペクトル解析、ウェーブレットコヒーレンス（同期の帯域を見る）
- **日周期の除去**（なぜ 5 日前方窓アノマリか）、季節性 ⟶ `_forward_window_anomaly`
- 埋め込み・ラグ・欠測（listwise deletion, gap 跨ぎ）⟶ `gap_guard`, `step_index`

### Tier 7 — フラックス／生態水文のドメイン知識 ⟶ 変数の意味, 解釈
- **渦相関法（eddy covariance）**の基礎、エネルギー収支（Rn = H + LE + G）
- 11 変数の物理: SW_IN(Rg), TA, VPD, TS, SWC(θ), H(γH), LE(γLE), NEE, GPP(GEP), RECO(GER)
- **炭素フラックス分割**: NEE = RECO − GPP、**昼間法(DT, Lasslop)/夜間法(NT, Reichstein)**、
  **u\* フィルタ**、VUT/CUT、`_vUT` 命名 ⟶ `sites.DEFAULT_VAR_MAP`
- **ギャップフィル（MDS: marginal distribution sampling）**と QC フラグ — なぜ I が膨らむ恐れ ⟶ FINDINGS の caveat
- 陸面-大気結合、蒸発散、**土壌水分-温度-呼吸**の関係、気孔・光合成の生理
- サイト生態の違い: 冷温帯落葉樹林 / 水田（湛水管理）/ 湿原（水位一定）⟶ 3 サイト対比の解釈

---

## Part 2. 効率的に学ぶためのプロンプト設計

教科書 DB を持つ AI から「深い・正確・オフライン保存できる」解説を引き出すコツ。

### 2.1 良いプロンプトの 7 原則
1. **1 プロンプト 1 概念**（詰め込むと浅くなる）。Tier の 1 項目ずつ。
2. **自分の前提を渡す**（現在の理解度・ゴール）→ AI が深さを較正できる。
3. **出力契約を明示**: 構成・深さ・記法・形式・長さ・引用。
4. **順序を指定**: 直感 → 定義 → 数式 → 小さな数値例 → **本研究への接続** → 落とし穴 → 自己テスト。
5. **オフライン形式**: 「1 つの自己完結した Markdown、数式は LaTeX、外部リンク不要」。
6. **出典要求**: 教科書名・章・式番号を明記させ、**標準と本研究の逸脱点**を指摘させる。
7. **反復**: 曖昧に残った箇所だけ深掘りする追撃プロンプトを用意。

### 2.2 マスターテンプレート（コピーして各トピックに使う）

```
あなたは情報理論・因果推論・生態水文学に精通した教師です。手元の教科書データベース
（Cover & Thomas, MacKay, Bossomaier, Runge/Spirtes/Pearl, Aubinet 等）を根拠に解説してください。

# 私の背景
- 専攻/レベル: <例: 学部4年、微積と線形代数は可、測度論は未修>
- 既知: <例: 相互情報量の定義は分かる>
- ゴール: <例: 自分の PCMCI 解析で「共通駆動を条件付ける」原理を説明できるようになる>

# 学びたいトピック
<例: 条件付き相互情報量 I(X;Y|Z) と、それが共通駆動をどう除くか>

# 本研究の文脈（なぜ学ぶか）
30分値フラックス11変数に Ruddell&Kumar 型 TE と PCMCI を適用。放射 Rg が全変数を共通駆動する
ため、ペアワイズ I が偽の同期を生む。I(X;Y|Rg) の減少率で共通駆動を切り分けている。

# 出力への要求（厳守）
1. 直感（比喩1つ）→ 2. 正式な定義（LaTeX、記法表つき）→ 3. 主要性質と証明の要点
→ 4. 手計算できる最小の数値例（3値程度、途中式つき）→ 5. 本研究の該当処理への接続
（上の文脈に具体的に結びつける）→ 6. よくある誤解・落とし穴 3つ → 7. 理解度自己テスト
（設問5・解答つき）→ 8. 出典（教科書名・章・式番号）と、本研究の手法が標準からずれる点。
- 形式: 1つの自己完結した Markdown 文書。外部リンク・画像なしで完結。数式は $...$ / $$...$$。
- 深さ: 学部上級〜院初級。要点を省かず、ただし冗長な前置きはなし。
- 不確かな点は「教科書での扱いが分かれる」と明示。憶測で断定しない。
```

### 2.3 そのまま使える記入済み例（重要トピック 3 本）

**(a) 推定バイアスと Miller-Madow**
```
（マスターテンプレートの背景・形式はそのまま）
# 学びたいトピック
ヒストグラム(プラグイン)エントロピー推定の有限標本バイアスと、Miller-Madow 補正。
なぜ次元が増える(2D MI→3D 条件付きMI)とバイアスが増え、条件付けで見かけ上MIが増える
「負の drop」が起きるのか。KSG(KNN)推定がなぜこれを緩和するか。
# 特に答えてほしい問い
- プラグイン H の期待バイアスは −(K−1)/(2N) 程度、という式の導出の要点
- m=11, 3次元, N≈1500 で具体的にどれくらいのバイアスか概算
- Miller-Madow が主要項をどう打ち消すか、限界は何か
- KSG が次元の呪いをどう回避するか（k近傍の考え方）
```

**(b) Transfer Entropy と Granger・条件付きMI の関係**
```
# 学びたいトピック
Transfer Entropy の定義(Schreiber)と Knuth 形(R&K式5)の同値性。TE が条件付き相互情報
I(Y_t; X_{t−τ} | Y_{t−1}) であること。Gaussian で Granger 因果と等価(Barnett 2009)である理由。
# 特に答えてほしい問い
- Knuth 形の4エントロピー項が I(Y_t;X_{t−τ}|Y_{t−1}) に一致する代数展開
- なぜ Y_{t−1} で条件付けるのか（自己ダイナミクスの除去）
- Tz = T/I が「forcing 対 feedback」を分ける情報理論的意味
- TE の「情報流」を因果と呼べる条件と、呼べない反例
```

**(c) PCMCI+ と共通駆動・多重比較**
```
# 学びたいトピック
PCMCI+ の2段階(PC条件選択 + MCI)。共通駆動と間接経路をどう除くか。自己相関がなぜ
偽リンクを生み、MCI がどう対処するか。ParCorr(線形)と CMIknn(非線形KNN)の違い。
因果十分性・定常性・忠実性の仮定が破れたとき何が起きるか(例:放射Rgへの逆向き偽リンク)。
# 特に答えてほしい問い
- d分離と条件付き独立の対応。共通駆動 Rg を条件集合に入れると偽リンクが消える理屈
- MCI が Bonferroni より検出力を保てる理由
- CMIknn のシャッフル検定が重い理由(計算量)と、knn/sig_samples の役割
- 我々の「逆向き Rg リンク」偽陽性の原因候補（線形近似 / 有限標本 / faithfulness 違反）
```

### 2.4 反復（深掘り）プロンプトの型
```
前回の「<トピック>」の説明のうち、<箇所> がまだ腑に落ちません。
具体的には <疑問>。別の比喩と、<例: 2変数の具体数値> で追加説明してください。
また私の理解「<自分の言葉での要約>」の誤りを指摘してください。
```

### 2.5 オフライン保存の実務
- 各回答を `study/<tier>_<topic>.md` として保存（この Markdown 群がオフライン教材になる）。
- 冒頭に**記法表**と**用語集**を作らせ、全ファイルで記法を統一させる（「以降この記法で」と固定）。
- 最後に **1 枚の要約カード**（定義・式・落とし穴・本研究での意味）を作らせると復習が速い。
- 数式が多い場合は「Pandoc で PDF 化できる純粋 Markdown で」と指定。

---

## Part 3. 学習順序と自己確認

1. **Tier 1→2→4** を先に（情報理論コア・推定バイアス・有意性）。ここが本研究の推定と検定の心臓部。
2. 次に **Tier 3**（TE と R&K）。1,2,4 が入っていれば式5の意味が腑に落ちる。
3. **Tier 5**（因果推論・PCMCI）。ここが「ペアワイズの限界→条件付け→多変量因果」という本研究の物語の核。
4. **Tier 6,7** は並行で可（時系列・ドメイン）。特に Tier 7 は解釈（3 サイト対比）に必須。
5. **自己確認**: 各 Tier 後に、対応するコード（下表）を読んで「式↔実装」を照合する。読めれば理解できている。

| 概念 | 読むコード |
|---|---|
| エントロピー/MI/CMI/TE | `information_theory.py` |
| バイアス補正 | `_entropy_of_indices(correct=True)` |
| サロゲート・多重比較 | `surrogate_*_stats`, `network.first_significant_peak(min_run)` |
| 結合分類 Tz | `network.classify_coupling` |
| 共通駆動条件付け | `condition_driver.py` |
| 多変量因果 | `causal_network.py` |
| アノマリ・欠測 | `preprocess.py` |
| ドメイン(変数/分割) | `sites.py`, `config.RK_VARS` |

---

## Part 4. 標準教科書・原著（AI に根拠として指定する用）

**情報理論**
- Cover & Thomas, *Elements of Information Theory* — エントロピー/MI/チェインルール/DPI の定番
- MacKay, *Information Theory, Inference, and Learning Algorithms* — 直感重視、無料PDF あり

**Transfer Entropy / 情報流**
- Bossomaier et al., *An Introduction to Transfer Entropy* — TE の教科書
- Schreiber 2000 (PRL); Knuth 2005; **Ruddell & Kumar 2009 (WRR)**（本研究の原著）
- **Goodwell & Kumar 2017 (WRR) TIPNets**; Williams & Beer 2010（PID）
- Barnett, Barrett & Seth 2009 (PRL) — TE↔Granger 等価
- Kraskov, Stögbauer & Grassberger 2004 (PRE) — KSG 推定

**因果推論**
- Pearl, *Causality*; Peters, Janzing & Schölkopf, *Elements of Causal Inference*（無料PDF）
- Spirtes, Glymour & Scheines, *Causation, Prediction, and Search* — PC アルゴリズム
- **Runge et al. 2019 (Science Advances / Nature Comms) PCMCI**; Runge 2018（総説）
- Sugihara et al. 2012 (Science) — CCM

**時系列**
- Hamilton, *Time Series Analysis*; Shumway & Stoffer, *Time Series Analysis and Its Applications*

**渦相関・生態水文**
- Aubinet, Vesala & Papale (eds.), *Eddy Covariance: A Practical Guide*
- Reichstein et al. 2005 (NT 分割); Lasslop et al. 2010 (DT 分割); Pastorello et al. 2020 (FLUXNET2015/ONEFlux)
- Bonan, *Ecological Climatology*; Monteith & Unsworth, *Principles of Environmental Physics*; Brutsaert, *Hydrology*

---

> 補足: この学習ガイドは本研究のコードと結論（同期支配 / 共通駆動交絡 / 2 群分離 / PCMCI 有向網）
> に対応づけてある。各 Tier を学んだら FINDINGS.md の該当主張を自分の言葉で再導出できるか試すと、
> 理解の穴が可視化できる。
