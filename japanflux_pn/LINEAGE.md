# 源流文献マップ — 各概念・手法の"元になった文献"を全部たどる

本研究が使う概念・手法の**一次源（origin）**を、アップロードされた 4 本のレビュー/枠組み論文から
抽出して整理した。各項目に **[出所]**＝どの論文がその原典を指しているか、を付す。末尾に
**「この 4 本には無いが我々の手法に必要な原典」**を別掲（正直に区別）。⚠ は DOI/正確な書誌が
要確認。引用前に一次資料で確定すること。

**アップロードされた 4 本（＝"源流をたどる入口"）**:
- **[R] Runge et al. 2019**, *Nature Communications* 10:2553, 因果推論の総説（地球システム）。DOI:10.1038/s41467-019-10105-3
- **[K] Krich et al. 2020**, *Biogeosciences* 17, 1033–1061, PCMCI で生物圏–大気の因果ネットワーク。DOI:10.5194/bg-17-1033-2020
- **[G] Goodwell, Jiang, Ruddell & Kumar 2020**, *WRR* 56, e2019WR024940, 情報理論は新パラダイムか（因果・相互作用）。
- **[N] Nearing, Ruddell, Bennett, Prieto & Gupta 2020**, *WRR*, 同（仮説検定）。DOI:10.1029/2019WR024918

---

## A. 因果の哲学的・理論的基盤（＝我々の「共通原因を差し引く」の根）
| 概念 | 原典 | 出所 |
|---|---|---|
| **共通原因原理（相関≠因果、依存なら因果か共通駆動）** ★我々の交絡除去の思想の根 | **Reichenbach, H. 1956. The Direction of Time.** Univ. California Press | [R]#8 |
| 経路解析（相関と因果の分離の始祖） | Wright, S. 1921. Correlation and causation. *J. Agric. Res.* 20, 557–585 | [R]#20 |
| 構造的因果モデル（SCM）・do計算 | **Pearl, J. 2000/2009. Causality: Models, Reasoning, and Inference.** Cambridge UP | [R]#4, [K], [G]（Pearl 1995 も） |
| 因果探索の古典（PC アルゴリズム） | **Spirtes, Glymour & Scheines 2000/2001. Causation, Prediction, and Search.** MIT Press | [R]#12, [K] |
| 因果推論の現代的定式化 | Peters, Janzing & Schölkopf 2017. Elements of Causal Inference. MIT Press | [R]#13, [K] |
| 多変量時系列因果の原理と問題 | Eichler, M. 2013. Causal inference with multiple time series. *Phil. Trans. R. Soc. A* 371 | [G] |

## B. 情報理論の基盤（我々の全指標の土台）
| 概念 | 原典 | 出所 |
|---|---|---|
| **エントロピー・情報量の創始** | **Shannon, C. E. 1948.** A mathematical theory of communication. *Bell Syst. Tech. J.* | [G], [N] |
| 情報理論の標準教科書 | Cover, T. & Thomas, J. 2006. Elements of Information Theory (2nd). Wiley | [G] |
| 情報量測度の一般化 | Csiszár, I. 1972. A class of measures of informativity… | [N] |

## C. 情報流・方向性（＝我々の Transfer Entropy）
| 概念 | 原典 | 出所 |
|---|---|---|
| **Transfer Entropy（情報流の向き）** ★我々の TE の原典 | **Schreiber, T. 2000. Measuring information transfer. *Phys. Rev. Lett.* 85, 461–464** | [R]#37, [G], [K] |
| **Granger 因果（予測に基づく因果）** | **Granger, C. W. J. 1969. *Econometrica* 37, 424–438** | [R]#9, [G], [K] |
| Granger 因果 = TE（ガウスで等価） | Barnett, Barrett & Seth 2009. *Phys. Rev. Lett.* 103, 238701 | [G] |
| Granger と有向情報理論の関係（レビュー） | Amblard & Michel 2013. *Entropy* 15(1), 113–143 | [G] |
| 「情報流」と「情報転送」の用語区別 | Lizier & Prokopenko 2010 | [G] |

## D. プロセスネットワーク（＝本研究の直接の土台, R&K）
| 概念 | 原典 | 出所 |
|---|---|---|
| **エコ水文プロセスネットワーク（TE でノード・リンク定義）** ★本研究の直接の親 | **Ruddell, B. L. & Kumar, P. 2009. Ecohydrologic process networks 1 & 2. *WRR* 45** (DOI:10.1029/2008WR007279 / …7280) | [K] |
| 情報駆動の自己組織化 | Kumar, P. & Ruddell, B. L. 2010. *Entropy* 12, 2085–2096. DOI:10.3390/e12102085 | [K], [G] |
| フラックスでの TE 応用群 | Ruddell et al. 2015; Gerken et al. 2018; Yu et al. 2019 | [K] |

## E. 部分情報分解・高次相互作用（＝我々の PID / 相乗・冗長）
| 概念 | 原典 | 出所 |
|---|---|---|
| **部分情報分解 PID（冗長/固有/相乗の分解, I_min）** ★我々の PID の原典 | **Williams, P. L. & Beer, R. D. 2010, 2011.** Nonnegative decomposition of multivariate information（arXiv:1004.2515 ⚠） | [G] |
| **時間的情報分割（synergy/uniqueness/redundancy）** ★我々の相乗概念の直接先行 | **Goodwell, A. E. & Kumar, P. 2017a. *WRR* 53, 5920–5942. DOI:10.1002/2016WR020216** | [G] |
| **TIPNets（プロセスネットワーク×PID）** | **Goodwell & Kumar 2017b. *WRR* 53, 5899–5919. DOI:10.1002/2016WR020218** | [G] |
| （気候ストレス応答）動的過程結合 | Goodwell et al. 2018. *PNAS* 115. DOI:10.1073/pnas.1800236115（⚠ [G] 本文/文献で確認） | [G] |

## F. 多変量因果探索（＝我々の PCMCI）
| 概念 | 原典 | 出所 |
|---|---|---|
| **PCMCI（条件選択＋MCI 検定）** ★我々の PCMCI の原典 | **Runge, J. et al. 2019. Detecting and quantifying causal associations in large nonlinear time series datasets. *Sci. Adv.*（[R]#23 は arXiv:1702.07007）** | [R]#23, [K]（Runge et al. 2019a） |
| 因果ネットワーク再構成の理論→実装 | Runge, J. 2018. *Chaos* 28, 075310 | [R]#24 |
| Tigramite（PCMCI 実装ソフト） | Runge, J. TIGRAMITE. github.com/jakobrunge/tigramite | [K] |

## G. 非線形力学系の因果（比較・文脈）
| 概念 | 原典 | 出所 |
|---|---|---|
| **収束クロスマッピング CCM** | **Sugihara, G. et al. 2012. Detecting causality in complex ecosystems. *Science* 338, 496–500** | [R]#11, [G]（Paluš et al. 2018 が改良） |
| 非線形相互依存の検出 | Arnhold et al. 1999. *Physica D* 134, 419–430 | [R]#10 |

## H. 情報理論によるモデル診断・仮説検定（＝我々の「モデルとの対決」の先行）
| 概念 | 原典 | 出所 |
|---|---|---|
| **情報理論でモデル評価・仮説検定** ★我々のモデル対決の直接先行 | **Nearing et al. 2020（本 [N] 論文）** | [N] |
| IT で認識論的/偶然的不確実性を推定 | Gong, Gupta et al. 2013. *WRR* 49, 2253–2273 | [N] |
| モデル構造妥当性の包括評価 | Gupta, Clark et al. 2012. *WRR* 48, W08301 | [N] |
| 複雑さの統計力学・情報理論的視点 | Balasis et al. 2013. *Entropy* 15(11), 4844–4888 | [G] |

## I. フラックス・領域知識（応用の文脈）
| 概念 | 原典 | 出所 |
|---|---|---|
| 生態系炭素フラックスの経年変動レビュー | Baldocchi, Chu & Reichstein 2018. *Agric. For. Meteorol.* 249, 520–533 | [R]#6, [K] |
| データ駆動地球科学・深層学習 | Reichstein et al. 2019. *Nature* 566, 195–204 | [R]#19 |

---

## J. ★重要：この 4 本には「無い」が、我々の手法に必要な原典（別途たどる）
これらは 4 本のレビューに含まれないので、**一次資料を別に当たる**こと（引用の穴になりやすい）。
| 我々が使う手法 | 原典（要一次確認） |
|---|---|
| **O-information（系全体の相乗/冗長, Ω）** ★novelty の核 | **Rosas, F. et al. 2019. *Phys. Rev. E* 100, 032305. DOI:10.1103/PhysRevE.100.032305** |
| **KNN 相互情報量推定（KSG）**（CMIknn の中身） | Kraskov, Stögbauer & Grassberger 2004. *Phys. Rev. E* 69, 066138 ⚠ |
| **サロゲート/順列による有意性検定** | Theiler et al. 1992. *Physica D* 58, 77–94（サロゲートデータ法）⚠ |
| **エントロピー推定のバイアス補正（Miller-Madow）** | Miller, G. 1955; Paninski 2003（バイアス補正の系）⚠ |
| **条件付き独立検定 CMIknn**（PCMCI 非線形版） | Runge, J. 2018. Conditional independence testing based on a nearest-neighbor estimator of CMI. *AISTATS* ⚠ |
| **相互作用情報 II（測度不変の相乗）** | McGill 1954; Bell 2003（多変量相互作用）⚠ |

---

## 使い方（源流のたどり方）
1. **まず A–F の "★" 印**（Reichenbach 共通原因 / Schreiber TE / Ruddell&Kumar プロセスネットワーク /
   Williams&Beer PID / Goodwell&Kumar 相乗 / Runge PCMCI）を読む＝本研究の骨格の原典。
2. 各 [R]/[K]/[G]/[N] の**該当論文の参考文献リスト**から、上表の原典の完全書誌を取る（本ファイルの
   DOI は暫定、原典で確定）。
3. **J の 6 本は別途**（特に Rosas 2019 O-information は novelty の核なので必ず一次確認）。
4. 論文の「関連研究」は **D→E→F→H の流れ**（プロセスネットワーク→PID/相乗→多変量因果→モデル診断）で
   書くと、本研究がこの系譜の"次の一歩"だと示せる。

> 各論文の本文中の該当箇所（式・定義）は、必要なら該当 PDF のページを指定して抜き出す（捏造しない）。
