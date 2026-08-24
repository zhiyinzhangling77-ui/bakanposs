# 文献ノート — 呼吸のメモリ・基質供給・衛星・情報理論×因果（Perplexity 検索・要一次確認）

> **出所と信頼度**：以下は Perplexity 検索で得た整理で、挙がった論文は実在の既知論文だが、
> **数値・主張は一次文献で必ず確認してから引用すること**（本ノートは所在の記録であって引用可能な確定情報ではない）。
> 捏造防止のため「ページ番号・原文引用・具体数値」は一次確認まで断定しない。日付：2026-08。

## 【重大】この分野の先行研究と、我々の位置づけ（誠実に）

我々が「発見」と思ったことの多くに**先行研究がある**。誇張を避けるため明記する：

| 我々の結果 | 先行研究（既存） | 我々の新規性の余地 |
|---|---|---|
| 呼吸に未観測の遅いメモリ（旗25）| **Migliavacca 2011・Cranko Page 2022・Stoy 2007 等＝ほぼコンセンサス**（正体=最近の光合成による基質供給, ラグ1〜5日）| メモリ自体は既知。**独立な検出法（残差自己相関×多生態系）と時間スケール特徴づけ**の切り口 |
| 因果骨格（PCMCI）・干ばつで NEE 脱結合 | **Krich 2020/2021 (Biogeosciences) が PCMCI を FLUXNET に適用済み**。季節で NEE 脱結合・気候帯を越えた機能収束を報告 | 骨格・脱結合は既出。**派生量監査(旗32)・リンク素性仕分け(旗35)・EBR品質重み**は独自寄り |
| O-info 相乗/冗長（旗14/34）| Goodwell & Kumar 2017/2018 (PNAS) が **TIPNet=PID(Williams-Beer)** で相乗冗長を実施 | ~~★O-information の FLUXNET 直接適用は査読研究に見当たらない~~ **【旗49で撤回】Eldhose & Ghosh 2025 ERL が O-Information を NEE に適用済み（ただし月次グリッドのFLUXCOM）。残る差分は半時間タワーデータ＋派生量監査のみ** |
| 乾燥草原=水分律速の呼吸・蒸発（旗31/36）| Budyko 枠組み・Wang 2025 単峰応答・Cable 2013 antecedent moisture | 情報理論での定量化・偏相関での炭素×水の統一署名は独自寄り |

**結論（方針・旗49で更新）**：発見の「新規性」は慎重に。呼吸メモリ=基質供給も、PCMCI×FLUXNET も、
**O-information×NEE も、水分依存Q10(DAMM)も既出**。
**残る独自性は (a)~~O-information の適用~~→**半時間タワー×派生量監査**という適用の仕方、(b) 方法論の監査一式
(旗32/35/39/43/44/46/47/48)、(c) フラックス／衛星／チャンバー3観測系の消去法設計**。既存を必ず引用し、上書きしない。

---

## ① 呼吸(RECO)のメモリ／履歴依存：残差自己相関・時間スケール・機構

**残差に自己相関が残る（温度・水分だけでは足りない）**
- **Migliavacca et al. (2011) Global Change Biology 17(1):390–409** — FLUXNET 104サイト。気候(T,SWC)のみのRECOモデルは
  時間変動の相当部分を説明できず、**日次GPPを共変数に加えると性能向上＝「最近の光合成による基質供給」が主要駆動**。
  最大LAI(LAI_MAX)が基準呼吸の空間変動を説明（フェノロジー）。
- **Cranko Page et al. (2022) Biogeosciences 19(7):1913–1936** — 豪12サイト。環境メモリ(lagged climate)を入れるとNEE予測が
  平均17%向上。SAM解析で気温ラグ≥2日(サイトにより4〜5日)、降水影響は21〜270日と幅広い。TBMがメモリを明示表現せず=炭素
  フィードバック不確実性の主因。
- Ruehr et al. (2010) Biogeochemistry 98:153–170／Savage & Davidson (2003) — 土壌呼吸残差に日〜週スケール自己相関。
- 残差は白色でなく**赤色雑音(遅い低周波成分)**、が一貫した知見。

**時間スケール（自己相関の半減）**
- **Stoy et al. (2007) Plant Cell & Environment 30(6):712–723** — 隣接するマツ林・広葉樹林で光合成→土壌呼吸のラグ τ_PR=1〜5日
  (特に1〜3日)。**落葉期にも同ラグ→全てが最近の光合成でなく、物理的CO₂拡散時間も寄与**（重要な留保）。
- **Cable et al. (2013) Global Change Biology 19(11):3435–3449** — 草原は過去2週間の水分が重要・累積効果6週間、灌木林は10週間。
- Cranko Page 2022／Liu 2019 — AR(1)残差の φ=0.5〜0.8 ＝**自己相関半減は約2〜5日**。
- Ping et al. (2023) Sci. Adv. — CCMで GPP→Ra ラグ 0.70±0.35ヶ月(≈21日)、GPP→Rh 1.21±0.45ヶ月(≈36日)。Rhの方がラグ長い。

**提案されている機構**
1. **基質供給(最近の光合成)**：Migliavacca 2011／Kuzyakov & Gavrichkova (2010) GCB 16:3386–3401（レビュー, ラグ数時間〜数日）
   ／Tang et al. (2005) GCB 11:1295–1304／Stoy 2007。森林で根呼吸が土壌呼吸の30〜50%(Högberg 2001, McDowell 2004)。
2. **フェノロジー**：Migliavacca 2011(LAI_MAX)／Zhang et al. (2018) Agric For Meteorol 259:178–189(GPP-Ts ラグ−25〜25週の
   ヒステリシス)／**Besnard et al. (2019) Nature Communications 10:633**(気候・植生のメモリがNEE季節変動に寄与)。
3. **微生物/根の動態**：Cable 2013／Ryan et al. (2015) GCB 21:4139–4153(antecedent 水分・温度がRECO応答を修飾)。
4. **熱・水の履歴**：Cranko Page 2022(14〜270日加重平均)／Feldman 2021(降雨→植物水分ピーク5日ラグ)／Huxman 2004・Cleverly 2013。

**未解決点**：メモリの生物過程(基質・微生物)vs物理過程(CO₂拡散・熱伝導)の分解／微生物プール回転時間との直接結合は
間接推論止まり／根の炭素供給をフラックスから直接推定する手法は未発達(同位体ラベリングが要る)／フェノロジー×気候メモリの
相互作用は未定量／TBMがメモリを表現しない。

## ② GPP→RECO ラグ（基質供給）と循環回避

- ラグの日数：日周1〜4時間、日〜週1〜5日、月次0.7〜1.2ヶ月(21〜36日, Ping 2023)。Rh>Ra でラグ長。樹高が高いほど長い傾向。
- 生態系差：森林1〜5日、草原/農地/湿地は日周1〜2時間が典型(Han 2014 Soil Biol Biochem, 湿地)。
- モデル誤差説明：Ono et al. (2025) Biogeosciences 22:5833 — Q10(Lloyd-Taylor)残差とGPP/PARが日周・週次で有意コスペクトル
  ＝温度のみで説明できない変動の主要部を基質供給が説明しうる。
- **循環(GPPとRECOが同じNEE由来)の回避**：
  (a) **昼NEE(GPP支配)/夜NEE(RECOのみ)を直接使う**(Ping 2023 CCM)。
  (b) **独立指標SIFを使う**：Shi 2021(SIFでGPP制約しNEE分割, 高VPDで優位)／Wang 2023(SIF≈GPP日周季節)／SMUrF(NASA 2025,
     SIFから4日平均GPP)／NN_SIF 2022(SIF組込NNでNEE分割, 高温・水制約で優位)。
  ※我々の旗25 は累積Rgを基質代理にして分割循環を避けた＝(a)(b)と同じ問題意識。

## ③ 衛星SIF→RECO/呼吸残差（まだ初期段階）

- SIF→GPP はほぼ標準化(NDVI/EVIより高相関)。**SIF→RECO/呼吸残差は明らかに少ない**。
- ケーススタディ：砂漠化鉱山域UAVハイパースペクトル＋SIFで土壌呼吸推定+26.8%向上、SIF>NDVI/NIRv(MDPI Remote Sens 18(10):1475)。
  根に基質輸送→根呼吸の遅い成分をSIFが捉える解釈。Ono 2025で Rh が GPP/PAR に2〜4時間ラグ。
- **未解決**：SIF→呼吸のラグ構造を全球多生態系で評価した研究は稀／SIFが温度残差のどの部分(Q10季節変化/水ストレス/基質)を
  捉えるか未分解／SIFが「基質代理」か「フェノロジー代理」か識別曖昧／衛星SIFの時空間解像度限界／乾燥地の非温度性呼吸を
  SIFがどこまで捉えるか未評価。
- footprint不一致：タワー(数十〜数百m) vs 衛星SIF(数百m〜数km)。対策=ダウンスケール(0.5°→0.1°でGPP相関向上, SIFtotal_01)、
  複数ピクセル平均(サイト限定だとN不足で季節検出できず広めが安定, というトレードオフ)。呼吸特化のfootprint補正理論は発展途上。

## ④ 衛星SMAP L4 根圏土壌水分→呼吸（深度問題は直接未検証）

- SMAP L4_C：表層(0–5cm)＋根圏(0–100cm)土壌水分を日次入力。**根圏→GPP制約(VPDと併用)、表層→Rh(乾燥抑制の非線形)** の役割分担。
- タワー検証：26コアサイトでNEE/GPP/RECOのubRMSE評価。**だが「浅いセンサーで捉えられない深層水分の効果を衛星が補えた」という
  深度依存の分解検証は公開報告に無い**（旗33の深度問題に直接答える実証はまだ）。
- 参考：Endsley et al. (2026) JGR Biogeosciences(SMAP L4C 10年総括)／Jones et al. (2017) Remote Sens Environ／Kimball 2012/2014
  (L4_Cアルゴリズム)／SMAP L4_C Validation Assessment V5(2021)/V8(2026) NASA GMAO／SMAP Handbook。
- 今後：多深度タワーθ×SMAP表層/根圏で呼吸残差を深度別に分解(回帰/PID)。乾燥地・草原(US-ARM, AU-Str, CN-Du2等)で根圏>表層かを検証。

## ⑤ 情報理論×因果探索 on FLUXNET（Ruddell & Kumar 2009 以降）

- **Ruddell & Kumar (2009a,b)** — transfer entropy で eco-hydrologic process network（本研究の出発点）。
- **Goodwell & Kumar (2017)／Goodwell et al. (2018) PNAS** — TIPNet：TEを情報分解し unique/synergistic/redundant に分解。
  2 CZOフラックスタワーで干ばつ初期に LE への相乗情報が増大、後期に脱結合＝**latent driver(深根・土壌水分再充填)の示唆**。
  標高でドライバーが異なる＝未観測(根深・保水力)の示唆。**＝相乗/冗長・latent driver は我々の主題と重なる先行研究**。
- **Krich et al. (2020) Biogeosciences** — **PCMCIを初めてFLUXNETに適用**。地中海サバンナで季節により因果構造変化、
  夏にNEEが気象から脱結合、Rg→NEEは実はRg→T→NEEの間接(共通原因制御)、latent driverの兆候。
- **Krich et al. (2021) Biogeosciences** — 119サイト。**気候帯を越えた機能収束**(温帯ピーク期≈熱帯雨林、干ばつ時は皆
  地中海乾期構造へ収束)＝気象が因果構造を支配。**＝我々の「背骨は生態系普遍」と強く重なる**。
- Yuan et al. (2021) — TEで森林伐採後の EF への情報フロー減少、E3SM評価。Kang et al. (2017) — 情報フローで傾斜地CO₂移流検出。
- 因果アルゴリズム：Runge et al. (2019) Sci Adv(PCMCI)／Runge et al. (2023) J-PCMCI+(共通latent contextを考慮)。
  latent許容：FCI(Spirtes 2001)、生態系適用は限定的。残差分解でlatent代理：Mahecha 2017／Sippel 2017(ICA/NMF)。
- O-information：Rosas et al. (2019)(Ω>0冗長/Ω<0相乗)／Scagliarini 2023,2024(O-info gradient=変数寄与)／Faes 2022(生理ネットに適用)。
  **★O-information を FLUXNET に直接適用した査読研究は Perplexity 検索で確認できず＝我々の O-info アトラス(21サイト)が新規性を持ちうる場所**。
- 手薄領域(＝新規性の候補)：O-info の FLUXNET 適用／latent を明示モデル化する因果発見の生態系適用／非線形・非定常CI／
  空間的 latent driver／**PCMCI + O-information の統合**。

## 主要文献リスト（一次確認用・実在既知論文）
- Ruddell & Kumar (2009a,b) Water Resources Research — process network / transfer entropy
- Goodwell & Kumar (2017) WRR; Goodwell et al. (2018) PNAS — TIPNet 情報分解
- Krich et al. (2020) Biogeosciences 17; Krich et al. (2021) Biogeosciences — PCMCI on FLUXNET
- Runge et al. (2019) Science Advances — PCMCI; Runge et al. (2023) — J-PCMCI+
- Rosas et al. (2019) Phys Rev E — O-information; Scagliarini et al. (2023/2024) — O-info gradients
- Migliavacca et al. (2011) Global Change Biology 17(1):390–409 — RECO 半経験モデル・基質供給
- Cranko Page et al. (2022) Biogeosciences 19(7):1913–1936 — 環境メモリ
- Stoy et al. (2007) Plant Cell & Environment 30(6):712–723 — 光合成-呼吸ラグ 1–5日
- Cable et al. (2013) Global Change Biology 19(11):3435–3449 — antecedent moisture 2–10週
- Kuzyakov & Gavrichkova (2010) GCB 16(12):3386–3401 — 光合成-土壌CO₂ ラグ レビュー
- Tang et al. (2005) GCB 11(8):1295–1304; Zhang et al. (2018) Agric For Meteorol 259:178–189
- Besnard et al. (2019) Nature Communications 10:633 — memory effects on NEE
- Ryan et al. (2015) GCB 21(11):4139–4153; Ruehr et al. (2010) Biogeochemistry 98(1):153–170
- Ping et al. (2023) Science Advances — CCM GPP→Ra/Rh ラグ; Ono et al. (2025) Biogeosciences 22:5833
- Endsley et al. (2026) JGR Biogeosciences; Jones et al. (2017) Remote Sensing of Environment — SMAP L4_C
- Wang et al. (2025) — 呼吸の土壌水分に対する単峰応答(FLUXNET 135サイト)  ※旗31/36 の裏づけ
- Pastorello et al. (2020) Scientific Data; Lasslop et al. (2010) GCB; Reichstein et al. (2005) GCB;
  Papale et al. (2006) Biogeosciences; Stoy et al. (2020) Boundary-Layer Meteorology; Hammerle et al. (2007)
  ※旗29/30/32 の裏づけ（別途 DATA_QUALITY / FLAGS_LOG に記録済み）

---

# 旗49（前提監査⑦）：一次確認による位置づけの訂正 — 2026-08

> **確認の水準（重要）**：出版社サイト（ScienceDirect / IOPscience）と OpenAlex API は本コンテナの
> egress プロキシで遮断されており、**一次PDFは読めていない**。以下は **WebSearch が返した書誌情報と
> 要約に基づく**。書誌（巻・ページ）まで含め、**論文に引用する前に一次で必ず確認すること**。
> 原文引用・独自の数値は一切書かない（捏造防止）。

## ★① 新規性主張(a)は成立しない — O-information はすでに NEE に適用されている
- **Eldhose, E. & Ghosh, S. (2025)「Exploration of synergistic and redundant information sharing from
  hydrometeorological variables to net ecosystem exchange」Environmental Research Letters 20(7),
  doi:10.1088/1748-9326/add8a7（2025-06-06 公開）**。
  **O-Information の枠組みを明示的に用いて**、気温・降水・VPD・陸水貯留量・PAR から NEE への
  相乗/冗長の寄与を分解している。報告された所見：**VPD-PAR は一貫して相乗**、**T-VPD と T-陸水貯留は主に冗長**。
- **使用データ（重要な差分）**：NEE は **FLUXCOM (RS+METEO)**＝機械学習でアップスケールした**月次グリッド積**。
  ほかに Berkeley Earth（気温）・GPCP（降水）・CSR GRACE/GRACE-FO Mascon（陸水貯留）・ERA5（VPD）・
  CERES SYN1deg（PAR）。＝**タワーの半時間渦相関データではない**。
- **我々の主張の訂正**：LITERATURE_NOTES 冒頭の「**O-information(Rosas 2019)を FLUXNET に直接適用した査読研究は
  確認できず＝ここが最も新規性がありうる**」は**誤り＝撤回する**。残る差分は
  「**半時間・サイト単位の渦相関データに適用した点**」と「**派生量監査を伴う点**」のみで、当初主張より遥かに薄い。
  （半時間タワーデータへの O-information 適用例は追加検索でも見つからなかったが、"見つからなかった"は
  "存在しない"ではない。）
- なお彼らの所見（VPD-PAR 相乗／T-VPD 冗長）は、我々の「放射共通駆動による冗長」と**収束的**で、
  対立ではなく補強関係にある。

## ★② 水分依存Q10は「再発見」で確定 — DAMM が機構として先行
- **Davidson, E.A., Samanta, S., Caramori, S.S., Savage, K. (2012)「The Dual Arrhenius and Michaelis-Menten
  kinetics model for decomposition of soil organic matter at hourly to seasonal time scales」
  Global Change Biology 18:371–384**。
  DAMM は Arrhenius（温度）と Michaelis-Menten（基質・酸素）を結合し、**土壌水分が酸素と可溶性炭素の拡散を
  制御することで、呼吸そのものだけでなく"見かけの温度感度"をも変える**という機構を与える。
- ＝**旗26/42/44 の「湿るほど温度感度が上がる」は、この確立した機構の再発見**。
  乾燥時は基質供給が律速するので温度感度が低く、再湿潤で基質制限が解けると温度感度が上がる、という
  DAMM の予測と我々のチャンバー実測（森林 18/36・草原 4/5）は同じ向き。**引用して上書きしないこと**。
- 関連：土壌の乾燥/再湿潤で Q10 が異なるかを直接扱った研究も存在する（Soil Biology & Biochemistry,
  「Is the temperature-dependence of soil respiration Q10 similar during soil drying and rewetting?」）。
  **本文未読**＝要一次確認。我々の結果と最も近い先行になりうる。

## ★③ 旗44 の前提（温度エイリアシング）は文献的に正しい
- **Lloyd, J. & Taylor, J.A. (1994)「On the temperature dependence of soil respiration」
  Functional Ecology 8(3):315–323**。
  指数(Q10)型・通常の Arrhenius 型は呼吸速度の不偏推定を与えず、**呼吸の"実効活性化エネルギーは温度と
  逆方向に変化する"**（＝**低温ほど見かけ Q10 が大きい**）とし、広い温度域で不偏な経験式を提示した。
- ＝旗44 で我々が「Q10 は温度そのものの関数だから、θビンが季節（＝温度帯）にエイリアスすると
  水分効果がゼロでも湿→高Q10 が出る」と述べた前提は、**文献的に支持される**。
  曲率項 e·Tc² を入れて交互作用 d を測った処置は、この文献に沿った正しい方向の対処だった。

## 位置づけへの帰結（誠実版）
1. **「O-information の生態系フラックス適用」は我々の新規性ではない**（撤回）。
2. **水分依存Q10 も新発見ではない**（DAMM の再発見。ただし**分割非依存のチャンバーで多サイト検証し、
   温度交絡を分離した**点は独自の検証価値がある）。
3. 残る独自性は**方法論の監査**——何が測定で何が計算か（旗32/35）、分割窓の分離（旗39）、
   多サイト族補正（旗43）、温度交絡の分離（旗44）、ギャップフィルの影響（旗46）、
   因果十分性の限界の実証（旗47）、後付け検出の out-of-sample（旗48）——と、
   **フラックス／衛星／チャンバーの3観測系を突き合わせた消去法の設計**（旗38/40/45）。
