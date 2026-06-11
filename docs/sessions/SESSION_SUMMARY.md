<<<<<<< HEAD
# Session Summary — Deep Root vs Irrigation in Mediterranean Almond/Cereal EC Sites

Last updated: 2026-05-11
Branch: `claude/quantify-water-divergence-LxPGr`
Repo: `zhiyinzhangling77-ui/bakanposs`

---

## 1. 研究背景・目的

### 当初の研究目的
Mediterranean の 2 つの Eddy Covariance(EC)サイトで観測された **「表層 SWC と蒸散(LE)の解離(decoupling)」**現象を定量化し、その原因を特定する。

- **Oran**(スペイン Albacete 近郊、雨養 winter cereal): 表層 SWC が乾くと LE が顕著に落ちる
- **Tarazona**(同地域、灌漑 almond): 表層 SWC が乾いても LE を維持する

### 当初仮説
「Tarazona アーモンドは **生物学的深根** を持ち、表層センサー(5cm 深)では捉えられない地下水・深層浸透水を利用して蒸散を維持している」

### なぜこの分析を行ったか
- FLUX 時系列の日変化(diurnal cycle)解析で、Tarazona が VPD 高条件で LE を増やす(典型的木本性反応)現象が確認されていた
- 衛星 SWC product(SMAP 等)から ET を推定するモデルが、灌漑農地で系統誤差を出す問題への科学的貢献の可能性
- 「深根」を表層 EC 観測のみで定量化できるかは方法論的に興味深い

### 前提条件
- Oran: 雨養 cereal、生育期 11–6月、SWC 単位は %
- Tarazona: 灌漑 almond、生育期 1–10月、SWC 単位は m³/m³(→ ×100 で% に変換)
- ERA5 で VPD を統一(EC タワー観測の VPD はサイト間で校正不一致)
- エネルギー閉合誤差を slope 補正で対応(Oran 0.74, Tarazona 0.71)

### 仮説の最終的な変遷(重要)
```
深根仮説(初期)
    ↓ v9-v12 で SDS=+0.05 が観測 → 一見支持
    ↓
Tarazona に灌漑あり(ユーザーから情報)
    ↓ 仮説と区別不能化
    ↓
v13: days_since_irrig 階層化 → 8+d で SDS 負の奇妙な値
    ↓ Claude が「8+d は実質"春期"」と気づく
    ↓
v14: 季節分離 + 夏期内 bucket → dose-response 判明
    ↓
最終結論: 灌漑依存(時定数 3-4日)、生物学的深根仮説は棄却
```

---

## 2. 使用データ

### 2.1 EC データ
| サイト | ファイル | 形式 | 期間 | 用途 |
|---|---|---|---|---|
| Oran | `/home/shion-nagamine/Dataset/Eddy data in Spain/Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.xlsx` | xlsx | 2018-01 – 2020-12 | v4 が読込 |
| Oran | `/home/shion-nagamine/Dataset/Eddy data in Spain/Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.csv` | csv | 同上 | v9/v10 が読込 |
| Tarazona | `/home/shion-nagamine/Dataset/Eddy data in Spain/EddyAlmond_Raw5years_withG.xlsx` | xlsx | 2020-06 – 2024-10 | v4 が読込(半時間データ) |
| Tarazona | `/home/shion-nagamine/Dataset/Eddy data in Spain/Daily_Summary_Filtered_forPred_ActEne26.csv` | csv | 同上 | v9/v10/v13/v14 が読込、**`Irrig_mm`, `Rain_mm`, `IrrigRain_mm` を含む** |

### 2.2 衛星・補助データ
| 種別 | パス | 用途 |
|---|---|---|
| Sentinel-2 NDWI Oran | `/mnt/hdd/Dataset/Sentinel2_NDWI/Oran_NDWI_Export.csv` | v9/v10 で NDWI マッチ |
| Sentinel-2 NDWI Tarazona | `/mnt/hdd/Dataset/Sentinel2_NDWI/TzM_NDWI_Export.csv` | 同上 |
| ERA5 2m T/Td | `/mnt/hdd/Dataset/ERA5_2m_Temperature/{year}.nc` (2018-2024) | VPD 統一計算 |
| GRACE-FO | `/mnt/hdd/Dataset/GRACE-FO_TWL/GRCTellus.JPL.200204_202602.GLO.RL06.3M.MSCNv04CRI.nc` | v9/v10 で地下水貯留量参照(空間解像度 ~300km、参考程度) |

### 2.3 サイト座標
- Oran: lat 38.82, lon −1.86
- Tarazona: lat 39.266, lon −1.9397
- GRACE は経度を 0–360° 系に変換(−1.86 → 358.14)

### 2.4 使用カラム
**v4(半時間データ)**
- Oran: `TIMESTAMP`, `SWC_1_1_1`, `SWC_2_1_1`, `SWC_3_1_1`, `H`, `LE`, `G`, `ET`, `VPD`, `FC_mass`(GPP proxy)
- Tarazona: `TIMESTAMP`, `H`, `LE`, `G`, `ET`, `VPD`, `VPD_kPa`, `co2_flux`, `SWC_avg` または `SWC_1_1_1`, `SWC_2_1_1`

**v13/v14(Tarazona 日次)**
- `date`, `Irrig_mm`, `Rain_mm`, `IrrigRain_mm`

### 2.5 前処理(v4 が実施)
1. **物理的 hard-clip**: LE ∈ [-100, 900], H ∈ [-300, 1000], G ∈ [-300, 600] (W/m²)
2. **昼夜分離 MAD ×7**: 昼夜混合では昼間の正常高値が誤除去される → 9-17h と それ以外で分離して MAD 計算、閾値 7
3. **時間帯別 z-score 4σ**: 24時間別の正規分布フィルタ
4. **エネルギー閉合補正**: OLS slope を計算し、`LE_corr = LE/slope`, `H_corr = H/slope`
5. **日次集約**: 30 min → daily(平均または合計)
6. **EF 計算**: `EF = LE_corr / (LE_corr + H_corr)` (closure-corrected)
7. **4象限分類**: SWC<p40 & VPD>p60 で normal/soil_dry/atm_dry/compound に分類
   - 閾値はサイト×is_growing 別に計算

---

## 3. 実施した分析(時系列順)

### 分析 v4(土台、ユーザー作成)

**目的**
半時間データから日次・分類済データセットを生成し、4 象限干ばつ分類で日変化(diurnal cycle)を可視化する。

**実施内容**
- xlsx 読込 → クレンジング → 外れ値除去 → 閉合補正 → 日次集約 → 4 象限分類 → diurnal 図生成
- パッチ追加(Claude 提案):
  ```python
  daily.to_parquet(f"{OUT_DIR}daily_classified_v4.parquet")
  json.dump(closure_slopes, open(f"{OUT_DIR}closure_slopes_v4.json","w"))
  ```

**結果**
- Oran: closure slope=0.7376, R²~0.92
- Tarazona: closure slope=0.7098
- daily 1939 行(Oran 914, Tarazona 689 だが一部期間外で集約後 1939)
- 4 象限分類で diurnal が綺麗に分離

**解釈**
- Oran growing diurnal: soil_dry < normal < compound < atm_dry(典型的浅根反応)
- Tarazona growing diurnal: soil_dry ≈ normal、atm_dry ≈ compound(典型的木本性、ただし灌漑効果含む)

**問題点**
- 直接 hypothesis test の指標を出していなかった → v11 へ

**次の仮説**
SDS のような単一指標でサイト間比較できる定量解析を作る。

---

### 分析 v9: 3軸統合(SWC × NDWI × VPD)

**目的**
「SWC 低 ∧ NDWI 高 ∧ VPD 高」の **同時成立日**を `deep_access` 群として抽出し、深根の直接証拠を探す。

**実施内容**
- NDWI を Sentinel-2 CSV から ±5日でマッチング
- SWC/VPD/NDWI の閾値を自動決定(KDE 谷検出 or unimodal 中央値)
- 6 群分類: deep_access / true_drought / atm_driven / compound_dry / wet_deep / wet_shallow
- GRACE 経度 0–360 変換

**結果**
- Oran deep_access: **0 日**
- Tarazona deep_access: **0 日**(NDWI 閾値=0.109 とした場合)
- Tarazona は v9 では deep_access=86日確認(NDWI 閾値が違うバージョン)

**解釈**
3軸同時成立が稀すぎて、`deep_access` 群がほぼ作れない。

**問題点**
- NDWI 絶対閾値(0.0 や 0.109)はサイト依存・季節依存
- 「SWC 低 ∧ NDWI 高」は自然界では同時に起きにくい(負の相関がある)

**次の仮説**
NDWI を anomaly(月別中央値からの偏差)で評価する v10 へ。

---

### 分析 v10: 交絡対策版

**目的**
v9 の根本問題を修正:
1. LE 絶対閾値 50 W/m²+ NDWI active-canopy filter
2. **Drought Response Ratio**: `LE_dry / LE_wet` の比をサイト間で比較
3. NDWI anomaly 化
4. ラグ相関を 30日移動平均からの偏差で計算
5. 月別 phenology-matched 比較

**結果**
Oran 雨養と Tarazona 灌漑で **絶対 LE 値が桁違い**(Oran ~3 W/m² vs Tarazona ~195 W/m²)を比較してしまう設計問題が露呈。

**解釈**
Cross-site の絶対比較は根本的に意味がない:
- 作物が違う(穀物 vs アーモンド)
- LAI が違う(0-4 vs 通年~4)
- 季節も違う(Oran 冬作 vs Tarazona 通年)
- "深根 vs 浅根" ではなく "種・気候の違い"を見ているだけ

**問題点**
DRR は形は良いが、cross-site でなく within-site で使うべき。

**次の仮説**
各サイト内で SDS を計算し、サイト固有の baseline で評価する。

---

### 分析 v11/v12: within-site SDS(SDS/VAS/DSO/CompoundDrop)

**目的**
v4 の 4 象限分類をそのまま使い、各サイト内で 4 つの干ばつ応答指標を bootstrap 5000回 + 95% CI で計算。

**指標定義**
```
SDS  = (LE_normal − LE_soil_dry) / LE_normal     (≈0 = 深根)
VAS  = (LE_atm_dry − LE_normal) / LE_normal      (高 = VPD駆動)
DSO  = (LE_compound − LE_soil_dry) / LE_soil_dry (高 = 土壌乾燥下VPD駆動)
CompoundDrop = (LE_normal − LE_compound) / LE_normal
```

**実施内容(v11)**
- parquet + json 入力
- 4 指標を LE_corr / EF_corr / Bowen_corr / ET の各変数で計算
- Mann-Whitney pairwise tests
- 8 枚 PNG 出力

**v11 で発覚したバグ**
- `DENOM_FLOOR = 5.0` を全変数一律に適用 → EF(0-1)/Bowen(~1)/ET(<5) で全て NaN
- Fix: 変数別 `DENOM_FLOORS = {"LE_corr": 5.0, "EF_corr": 0.05, "Bowen_corr": 0.1, "ET": 0.3}`
- Bowen は `|Bowen|<20` で極値クリップ追加

**v12 で改善**
- 8 枚の図を自己説明的に再設計
- 休眠期警告(LE_normal < 30 W/m² で注意)
- 自動判定ロジック(CI 重なり判定)

**結果(v12 Growing season, LE_corr)**
| 指標 | Oran | Tarazona |
|---|---|---|
| SDS | +0.312 [+0.21, +0.41] | +0.051 [-0.07, +0.15] |
| VAS | +2.498 | +1.167 |
| DSO | +1.581 | +1.015 |
| CompoundDrop | -0.776 | -0.912 |

→ 一見、深根仮説支持。CI が分離している。

**問題点**
- ユーザー指摘: 「Oran の生育期 11-6月、Tarazona 1-10月で全く違う」→ cross-site での主張は危うい
- さらに後で発覚: 季節間プーリングのアーティファクト(spring の高 SWC が "normal" 群を支配、summer の低 SWC が "soil_dry" を支配)
- v12 の +0.05 SDS は **アーティファクト**だった

**次の仮説**
衛星データで深根仮説を補強できないか。

---

### 分析 v13: 灌漑経過日数で深根 vs 灌漑を切り分け

**契機**
ユーザーから「Tarazona は灌漑あり」の情報。
→ 深根仮説と灌漑仮説が観測上区別不能。Irrig_mm column が CSV にある。

**目的**
「最後の灌漑からの経過日数(days_since_irrig)」でデータを階層化し、灌漑離脱後も LE 維持されるかを検定する。

**実施内容**
- Tarazona daily CSV から `Irrig_mm > 0.5` を灌漑イベント定義
- Days since last irrigation を計算(NaN for never-irrigated)
- Bucket: 0-1d, 2-3d, 4-7d, 8+d
- 各 bucket 内で SDS 計算

**結果**
| Bucket | LE SDS | CI | n_n, n_s |
|---|---|---|---|
| 0-1d | +0.207 | [-0.06, +0.42] | 46, 42 |
| 2-3d | +0.176 | [-0.07, +0.40] | 37, 34 |
| 4-7d | +0.310 | [-0.01, +0.44] | 43, 41 |
| **8+d** | **−0.151** | [-0.47, +0.13] | **249**, 69 |

灌漑統計:
- 138 イベント 2020-06 〜 2024-10
- 月別: 5月=18, 6月=17, 7月=47, 8月=35, 9月=14, 10月=7
- 11月–4月: 灌漑停止(完全に)

**解釈**
8+d バケットの SDS が負・サンプル数が異常に多い(249/69)。
傾き−0.094、相関−0.61。
→ 自動判定が「逆相関、不自然」と判定。

**問題点(Claude が気づく)**
**8+d バケットの 249 サンプルは大半が 1-4月**(灌漑停止期=春期)。
バケット番号 ≈ 季節 になっていた。経過日数の純粋効果を測れていない。

**次の仮説**
季節 × バケットの交絡を分離する v14 へ。

---

### 分析 v14: 季節分離 + 夏期内 bucket(現在の最終版)

**目的**
v13 の交絡を解決:
- Test A: 夏期(7-9月)限定で bucket 解析 → 夏期内で純粋な経過日数効果
- Test B: 季節別(spring 1-4月 / shoulder 5,6,10月 / summer 7-9月)で SDS → 灌漑停止期の純粋テスト

**実施内容**
- 季節定義: SEASONS = {"spring (1-4月)": [1,2,3,4], "shoulder (5,6,10月)": [5,6,10], "summer (7-9月)": [7,8,9]}
- 夏期内 bucket × 変数(LE_corr/EF_corr/ET)で SDS bootstrap
- 季節 × サイト で SDS bootstrap(Oran が Test B の control)
- 月 × bucket の heatmap で交絡可視化
- Verdict 自動判定: |SDS|≤0.15 OR CI が 0 を含む → 深根帯

**結果(v14)**

**Test A: 夏期 bucket(Tarazona LE_corr)**
| Bucket | SDS | CI | n_n, n_s |
|---|---|---|---|
| 0-1d | +0.051 | [-0.12, +0.43] | 28, 25 |
| 2-3d | +0.004 | [-0.29, +0.36] | 24, 21 |
| **4-7d** | **+0.384** | **[+0.11, +0.51]** | 27, 22 |
| **8+d** | **+0.283** | **[+0.14, +0.43]** | 21, 28 |

→ **明確な dose-response**: 灌漑直後 3 日は SDS≈0、4 日以降は SDS が立ち上がり CI が 0 を超える。

**Test B: 季節別 SDS (LE_corr)**
| 季節 | サイト | SDS | CI | n_n, n_s |
|---|---|---|---|---|
| spring | Oran | **+0.689** | [+0.62, +0.72] | 103, 86 |
| spring | Tarazona | +0.616 | [-0.98, +0.83] | 181, **15** (⚠) |
| shoulder | Tarazona | +0.429 | [+0.29, +0.55] | 94, 75 |
| summer | Tarazona | +0.272 | [+0.14, +0.35] | 100, 96 |

**解釈**
1. **Tarazona の夏期は灌漑依存パターン**: 0-3d で SDS≈0、4+d で SDS=0.28-0.38 → 灌漑効果の dose-response
2. **真の深根なら 8+d でも SDS≈0 のはず → 実際は +0.28 → 深根仮説棄却**
3. **Tarazona の shoulder/summer は SDS≈0.27-0.43 で表層感受性あり** → v12 の "+0.05" はアーティファクト
4. **Spring Tarazona は n_s=15 で判定不能**(冬の降雨で SWC>p40 が大半)

**重要な発見: 自動判定ロジックのバグ**
`in_band` 関数が「CI が 0 を含む → 深根」と判定するが、**n_s=15 で CI が広すぎる場合に偽陽性**。
Spring Tarazona の CI [-0.98, +0.83] は単に sample 不足の現れで、深根の証拠ではない。
→ v15 で fix 必要: CI 幅 > 0.5 で判定保留、n_s < 30 で判定保留など

**次の仮説**
- 真の深根仮説は棄却
- **灌漑による SWC-ET decoupling time-constant ~3-4日** が新しい主張
- 衛星 ET 推定モデル(SMAP-based)への応用が contribution として残る

---

## 4. コード変更履歴

### v4 (ユーザー側、既存)
- ファイル: ユーザー側 `Eddy Covariance Flux Analysis v4` script
- 変更箇所: 末尾に 2 行追加(Claude 提案)
  ```python
  import json
  daily.to_parquet(f"{OUT_DIR}daily_classified_v4.parquet")
  with open(f"{OUT_DIR}closure_slopes_v4.json", "w") as f:
      json.dump(closure_slopes, f)
  ```
- 理由: v11 以降の中間ファイル連携のため
- 現在の状態: ユーザーが反映済みの想定。Claude 側 repo にはなし。

### v9 (`analysis_A_v9.py`)
- 場所: `/home/user/bakanposs/analysis_A_v9.py`
- 内容: 3軸統合(SWC×NDWI×VPD)で deep_access 抽出を試行
- 状態: 完成、ただし deep_access=0 で目的未達

### v10 (`analysis_A_v10.py`)
- 場所: `/home/user/bakanposs/analysis_A_v10.py`
- 内容: DRR + NDWI anomaly + lag correlation
- 状態: 完成、ただし cross-site 比較の根本問題が残る

### v11 (`analysis_A_v11.py`)
- 場所: `/home/user/bakanposs/analysis_A_v11.py`
- 内容: 4-class SDS/VAS/DSO/CompoundDrop bootstrap
- 変更履歴:
  - 初版: DENOM_FLOOR=5.0 一律(バグ、EF/Bowen/ET 全 NaN)
  - 修正: `DENOM_FLOORS` を変数別に変更、`BOWEN_CLIP=20` 追加
- 状態: 完成

### v12 (`analysis_A_v12.py`)
- 場所: `/home/user/bakanposs/analysis_A_v12.py`
- 内容: 8 図を自己説明的に再設計、休眠期警告、自動判定
- 状態: 完成。**結果は今思えば季節プーリング artifact**

### v13 (`analysis_A_v13.py`)
- 場所: `/home/user/bakanposs/analysis_A_v13.py`
- 内容: 灌漑経過日数で階層化、Tarazona daily CSV `Irrig_mm` 利用
- 状態: 完成。8+d バケットの交絡が判明、v14 へ。

### v14 (`analysis_A_v14.py`)
- 場所: `/home/user/bakanposs/analysis_A_v14.py`
- 内容: 季節分離 + 夏期内 bucket、Oran spring を control
- 状態: 完成。**現時点の最終解析**
- 既知のバグ: `in_band` ロジックが「CI が 0 を含む」を無条件に深根と判定する偽陽性

### requirements.txt
- 場所: `/home/user/bakanposs/requirements.txt`
- 内容: v4-v14 全部の依存 (numpy, pandas, scipy, matplotlib, xarray, netCDF4, h5netcdf, dask, bottleneck, openpyxl, pyarrow)
- 状態: 完成・push 済

---

## 5. 試して失敗したこと(再試行禁止)

### F1: NDWI 絶対閾値での deep_access 抽出(v9)
- 試行: `NDWI > 0.0` または `> 0.109` を「植被が水を含む」基準にした
- 結果: deep_access 日が 0
- 理由: NDWI 絶対値はサイト・季節依存性が強すぎる。Oran は全体的に負、Tarazona は全体的に正で、固定閾値は意味を成さない
- 教訓: 必ず anomaly(月別中央値からの偏差)を使う

### F2: cross-site 絶対 LE 比較(v10)
- 試行: Oran 雨養日と Tarazona 灌漑日の絶対 LE を Mann-Whitney 検定 → p=1.5e-47
- 結果: 統計的に超有意だが、科学的に無意味
- 理由: Oran 雨養 cereal の dry SWC 日は収穫後/裸地、Tarazona は灌漑活動中で、種・LAI・季節すべてが異なる
- 教訓: cross-site の絶対比較はしない。各サイト内で baseline 比を取る

### F3: 季節間プーリング SDS(v12)
- 試行: Tarazona の生育期 1-10月を全部 1 群として SDS 計算 → SDS=+0.05(深根支持に見えた)
- 結果: v14 で季節分離すると spring/shoulder/summer 各 +0.27〜+0.62 で実は表層感受性あり
- 理由: spring が高 SWC × 低 LE で "normal" 群を支配、summer が低 SWC × 高 LE で "soil_dry" 群を支配 → "soil_dry" の方が見かけ上 LE が高くなる artifact
- 教訓: SWC 閾値や群定義は **季節別に**評価する

### F4: ラグ相関(v9)
- 試行: SWC と NDWI の lag correlation で深根を検出
- 結果: 季節周期に支配されてフラット、解釈不能
- 理由: SWC と NDWI 両方が年周期を持つので Pearson は周期成分で決まる
- 教訓: lag 解析は必ず detrend / deseasonalize してから

### F5: GRACE-FO を point sample で利用(v9)
- 試行: lat=38.82, lon=-1.86 で GRACE TWS を抽出
- 結果: 経度系の不一致で 0 ヶ月分(NaN)
- Fix: 経度 0-360 系に変換(-1.86 → 358.14)
- 教訓: GRACE は ~300km 解像度なのでサイト固有の解析には適さない(参考程度)

### F6: 自動判定の `in_band` 偽陽性(v12, v14 で再発)
- 試行: 「CI が 0 を含む → 深根帯」と単純判定
- 結果: n が小さく CI が広いケース(Tarazona spring n_s=15)で誤った "深根支持" 判定
- 改善案: CI 幅 > 0.5 で保留、n_s < 30 で保留
- 状態: v15 で fix 予定

### F7: DENOM_FLOOR 一律 5.0(v11 初版)
- 試行: 全変数で同じ最小分母閾値
- 結果: EF/Bowen/ET の指標が全 NaN
- Fix: 変数別 DENOM_FLOORS
- 教訓: スケールが違う物理量に同じ定数閾値を使わない

---

## 6. 現在の研究上の論点

### 論点 A: 深根仮説の最終判定
**現在のスタンス**: 棄却寄り
- 夏期 4+d bucket で SDS が有意に正 → 灌漑効果が消えると Tarazona も表層感受性を持つ
- 真の深根なら経過日数によらず SDS≈0 のはず
- ただし春期データ不足(n_s=15)が完全否定を弱める

### 論点 B: 春期の解釈
- Tarazona spring の SDS = +0.62 [-0.98, +0.83] は **判定不能**
- 解釈案 1: 冬の貯留水で SWC>p40 が大半 → 干ばつ条件成立せず → 深根仮説検定不能
- 解釈案 2: 春は花期/葉展開期で transpiration phenology が違う → SDS の意味自体が違う
- どちらにせよ春期では深根の決定的証拠にも反証にもならない

### 論点 C: Irrigation Decoupling Time-constant の意義
- 0-3d で SDS≈0、4+d で SDS が立ち上がる → time-constant ~3-4日
- これは **灌漑農地で衛星 SWC ET 推定モデルの系統誤差を生む現象**
- 論文化の新しい主軸候補

### 論点 D: Phenology vs Hydrology
- Oran と Tarazona の比較は phenology(作物の生育周期)も違う
- "Spring SDS" は両サイトとも生育期だが、Oran は cereal 中期で Tarazona は葉展開期
- 同じ "spring" の SDS=+0.69 vs +0.62 は両者似ているが、phenology 差が交絡しているか

---

## 7. 未解決課題(優先順)

### P1(最優先): 自動判定ロジックのバグ修正
- `in_band()` が CI 幅・サンプル数を考慮しない
- v14 の Tarazona spring 「深根支持」結論は誤り
- 修正案: `(lo <= 0 <= hi) AND (hi - lo < 0.5) AND n_s >= 30`

### P2: Recovery time analysis(灌漑応答時定数の精緻化)
- 現在の 0-1d / 2-3d / 4-7d / 8+d バケットは粗い
- 1日刻みで SDS をプロット → exponential fit → time-constant τ を抽出
- 論文の headline figure 候補

### P3: 春期データ不足の補正
- Tarazona spring で SWC<p40 が n=15 しかない
- Spring 内で別の SWC 閾値(例: p25)を使えば soil_dry サンプル増える可能性
- ただし定義変更による比較性損失とトレードオフ

### P4: Phenology 補正
- Oran vs Tarazona の比較は phenology 交絡あり
- LAI または NDVI で正規化した LE/EF 指標を併用検討
- データに NDVI あれば実装可能

### P5: 衛星補強(別データソース)
- MODIS LST × ERA5 Tair で ΔT(canopy stress)を独立検証
- ECOSTRESS 70m thermal ET で空間検証
- Sentinel-1 SAR で灌漑パルスを Irrig_mm と独立に検出
- どれも追加データ取得が必要

### P6: 多深度 SWC の必要性
- 5cm のみでは drip wetted bulb(20-50cm)を直接観測できない
- 既存タワーに 30cm/50cm SWC センサーあれば再解析価値大
- データ要確認

### P7: 論文骨子の準備
- 仮説変更後の新フレーミング: "Drip irrigation decouples SWC-ET on 3-4 day timescale"
- Introduction, Methods, Results, Discussion の骨組み
- v14 の図を中心に構成

---

## 8. 次セッション開始時にやること(順序)

1. **`SESSION_SUMMARY.md` を必ず最初に読む**(これ)
2. `/home/user/bakanposs/` の状態確認: `git log --oneline -20`, `ls`
3. 中間ファイル状態確認: `daily_classified_v4.parquet`, `closure_slopes_v4.json` が存在するか(v4 実行済みか)
4. v14 の output 確認: `output_analysis_A_v14/fig01-04.png` を視認
5. **ユーザーに優先課題を確認**:
   - P1(verdict bug fix)で v15 を書くか?
   - P2(recovery time analysis)を実装するか?
   - 論文骨子(Path 1)を書き始めるか?
   - 衛星補強(MODIS/ECOSTRESS, P5)を試すか?
6. **絶対に提案しないこと**:
   - cross-site の絶対 LE 比較(F2 で失敗済み)
   - NDWI 絶対閾値(F1 で失敗済み)
   - 季節プーリングなしの SDS(F3 で artifact 確認済み)

---

## 9. Claude への引き継ぎ指示

### 9.1 一行サマリー
"Mediterranean almond(Tarazona, 灌漑)と cereal(Oran, 雨養)の EC データで、当初「深根仮説」を検定したが、季節 + 灌漑経過日数の制御後、Tarazona は **灌漑による 3-4日の SWC-ET decoupling** を示すことが判明。深根仮説は棄却寄り、新しい論文方向性は灌漑による衛星 ET 推定モデルへの含意。"

### 9.2 重要な判断と制約
- **絶対にしてはいけない比較**: Oran と Tarazona の同日絶対 LE 比較(種・LAI・季節すべて違う)
- **必ず seasonal-stratified に**: 全季節プーリングは v12 で artifact を生んだ
- **n_s < 30 のときは結論保留**: spring Tarazona の偽陽性事例あり
- **SWC センサーは 5cm のみ**: drip wetted bulb 直接観測ではない、間接シグナルとして扱う

### 9.3 ユーザーの好み
- 図は **自己説明的**(タイトルに結論、軸ラベル明確、サンプル数表示)を強く好む
- 解析は **段階的に**(各バージョンが前の問題を1つ修正)
- バグ・問題点は**正直に指摘**することを歓迎(誤った "深根支持" の自動判定を指摘した経緯あり)
- コードは **コメントで "なぜそうしたか" を書く**ことを好む
- バージョン管理: `analysis_A_vN.py` と output ディレクトリ `output_analysis_A_vN/` で並走

### 9.4 当面の科学的方向性
- 論文化の方向: **"Drip irrigation decouples surface SWC from canopy ET on a 3-4 day timescale"**
- v14 の Test A 夏期 bucket dose-response が headline figure 候補
- 衛星 ET モデル(SMAP-based)への含意が contribution

### 9.5 注意事項
- v12 の "Tarazona SDS=+0.05 → 深根支持" は **artifact、引用しない**
- v13 の "SDS 逆相関、不自然" の自動判定は **間違い**、8+d が春期だった
- v14 の Tarazona spring "深根支持" 判定も **偽陽性**、n_s=15 が問題

---

## 10. 重要ファイル一覧

### コード(解析A 本流)
| パス | 内容 |
|---|---|
| `/home/user/bakanposs/analysis_A_v9.py` | NDWI 3軸統合(失敗例として保存) |
| `/home/user/bakanposs/analysis_A_v10.py` | DRR + cross-site(失敗例として保存)。**現在 data_loaders import に書き換え済** |
| `/home/user/bakanposs/analysis_A_v11.py` | within-site SDS(DENOM_FLOOR bug fix 済) |
| `/home/user/bakanposs/analysis_A_v12.py` | 自己説明的 figures(8 枚) |
| `/home/user/bakanposs/analysis_A_v13.py` | 灌漑 days_since_irrig 階層化 |
| `/home/user/bakanposs/analysis_A_v14.py` | **解析A 現最終版**: 季節分離 + 夏期 bucket |

### コード(共通ライブラリ・並行解析)
| パス | 内容 |
|---|---|
| `/home/user/bakanposs/data_loaders.py` | **解析A/B/C 共通の EC ローダ**。解析C で発見した v9 ローダの 5 つのバグを修正したクリーン版。`load_oran_ec_clean`, `load_tarazona_ec_clean`, `normalize_swc` を提供。v10 が import 済 |
| `/home/user/bakanposs/analysis_C_v1.py` | **解析C(並行)**: MODIS NDVI を独立軸として導入し、フェノロジー × フラックスを解析。解析A の深根仮説に対する状況証拠を別観点(NDVI vs EF/LE/GPP)から評価 |

### 共通インフラ
| パス | 内容 |
|---|---|
| `/home/user/bakanposs/requirements.txt` | 全環境依存 |
| `/home/user/bakanposs/SESSION_SUMMARY.md` | **このドキュメント** |

### data_loaders.py の主要修正点(解析C 経由で判明)
1. Oran timestamp の自動推定エラー(`%Y/%m/%d` 確定で残り NaT 化)を 3段フォールバック化
2. -9999 等のセンチネル値マスク
3. Oran VPD: hPa → kPa 統一(÷10)、`|VPD|>10 kPa` 外れ値除去
4. Tarazona ET: `ET_avg` ではなく `ET_sum`(日合計 mm)を採用
5. Tarazona VPD: `VPD_kPa` 列を直接採用(`VPD_mean` は Pa で単位不整合)

**最終単位の統一**:
- LE/H/G/Rn: W/m²
- ET: mm/day
- VPD: kPa
- SWC: %(m³/m³ は `normalize_swc` で ×100)

⚠️ **重要**: 解析A v9-v11 は古いローダで実行されていたため、上記バグの影響を受けている可能性あり。次セッションで v15 を書く際は **必ず `data_loaders` を import** すること。

### 解析C(並行)の目的
- 目的1: 生育期/非生育期を NDVI から客観定義(GROWING_MONTHS の妥当性検証)
- 目的2: NDVI vs EF/LE/GPP_proxy を Oran vs Tarazona で比較
- 目的3: サイト間フェノロジー位相の定量化
- 目的4(A 連結): 低 SWC × 高 NDVI 期で Tarazona が EF を保てるか = 深根仮説の状況証拠の再評価
- 目的5(B 連結): EC LE と NDVI×Rn 簡易 GPP プロキシの整合性

→ 次セッションで A/C 結果を統合して論文化の方向を確定する選択肢あり。

### 中間ファイル(ユーザー側で生成、Claude 側 repo にはなし)
| パス | 内容 | 必要なバージョン |
|---|---|---|
| `daily_classified_v4.parquet` | v4 が出力する日次・分類済データ | v11-v14 が読込 |
| `closure_slopes_v4.json` | エネルギー閉合スロープ {Oran:0.7376, Tarazona:0.7098} | v11-v14 が読込 |

### 入力データ(ユーザー側 /home/shion-nagamine/, /mnt/hdd/)
| パス | 内容 |
|---|---|
| `/home/shion-nagamine/Dataset/Eddy data in Spain/Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.xlsx` | Oran 半時間 EC(v4 用) |
| `/home/shion-nagamine/Dataset/Eddy data in Spain/EddyAlmond_Raw5years_withG.xlsx` | Tarazona 半時間 EC(v4 用) |
| `/home/shion-nagamine/Dataset/Eddy data in Spain/Daily_Summary_Filtered_forPred_ActEne26.csv` | Tarazona 日次(`Irrig_mm`, `Rain_mm`, `IrrigRain_mm` 含む。v9/v10/v13/v14 用) |
| `/mnt/hdd/Dataset/Sentinel2_NDWI/Oran_NDWI_Export.csv` | Sentinel-2 NDWI(v9/v10 用) |
| `/mnt/hdd/Dataset/Sentinel2_NDWI/TzM_NDWI_Export.csv` | 同上 Tarazona |
| `/mnt/hdd/Dataset/ERA5_2m_Temperature/{year}.nc` (2018-2024) | ERA5 VPD 計算用 |
| `/mnt/hdd/Dataset/GRACE-FO_TWL/GRCTellus.JPL.200204_202602.GLO.RL06.3M.MSCNv04CRI.nc` | GRACE-FO TWS(参考) |

### 出力(ユーザー側)
| ディレクトリ | 内容 |
|---|---|
| `output_analysis_A_v9/` | v9 figures + CSV |
| `output_analysis_A_v10/` | v10 figures + CSV |
| `output_analysis_A_v11/` | v11 figures + CSV |
| `output_analysis_A_v12/` | v12 figures + CSV |
| `output_analysis_A_v13/` | v13 figures + CSV |
| `output_analysis_A_v14/` | **最終結果**: fig01-04 + CSV 2本 |

### 重要 git コミット
- `e9efc97` v9 added
- `2069be8` v10 + requirements
- `9ca8f73` v11 added
- `f547c3d` requirements consolidated
- `7d90947` v11 fix (DENOM_FLOORS)
- `0c12aea` v12 refactored
- `bcaf1a3` v13 irrigation analysis
- `1dc0ecb` v14 seasonal × bucket
- `2ade97e` (parallel) data_loaders.py 抽出, analysis_C_v1.py 追加, v10 を data_loaders import に書き換え
- `df16969` SESSION_SUMMARY.md 追加

---

## 11. 解析C 関連の補足(並行作業)

### 11.1 解析C の独立性
解析C(`analysis_C_v1.py`)は解析A とは独立した観点で同じ EC データを分析。
深根仮説の代替確認として、NDVI と LE/EF の関係を見る。

### 11.2 解析C で使う追加データ
- **MODIS MOD13Q1 NDVI**: 16-day, 250m
- 取得方法: Google Earth Engine `scripts/gee_extract.js`(別途存在の想定)
- ⚠️ 次セッションで実際の NDVI CSV パスを確認すること

### 11.3 A と C の統合可能性
- A の v14 「灌漑による SWC-ET decoupling」結論
- C の NDVI による phenology 補正
- → 両者を組み合わせれば「LAI/Phenology 差を補正した上での灌漑効果」が示せる可能性

### 11.4 解析B の存在(?)
- `data_loaders.py` のコメントに「解析A/B/C 共通」と記載
- 解析B のスクリプトはまだ未確認(`/home/user/bakanposs/` には `analysis_B_*.py` が存在しないか?要再確認)
- 次セッション開始時に `ls /home/user/bakanposs/` で確認すること
=======
# 統合 SESSION_SUMMARY — 解析A / B / C

**最終更新**: 2026-05-11
**統合元**:
- A: `claude/quantify-water-divergence-LxPGr` / `SESSION_SUMMARY.md` (SHA `499a4c7`)
- B: `claude/compare-ec-satellite-et-ZnENi` / `SESSION_SUMMARY.md` (SHA `e2b7e7a`)
- C: `claude/ndvi-flux-analysis-iaS0b` / `reports/SESSION_SUMMARY.md` (SHA `02ff1d1`)
**統合先ブランチ**: `claude/analysis-presentation-prep-BQRnM`
**統合方針**: 要約ではなく**意思決定履歴を保持したまま** 3 系列を時系列・横断観点で再構成する。各事実には A/B/C いずれの根拠かを明示する。

---

## 1. 研究背景（統合版）

### 1.1 共通の研究対象
スペイン半乾燥地 (Albacete) の Eddy Covariance (EC) 2 サイト:

| サイト | lat / lon | 作物 | 灌漑 | 期間 | 生育期 |
|---|---|---|---|---|---|
| **Oran** | 38.82 / −1.86 | winter cereal (vetch/wheat/pea 輪作) | 雨養 | 2018-01 〜 2020-12 (半時間) | 11–6月 (主成熟 Feb–Jul, 初期生育 Oct–Jan) |
| **Tarazona (TzM)** | 39.266 / −1.9397 | drip-irrigated almond (落葉) | drip | 2020-06 〜 2024-10 (日次) | 1–10月 |

### 1.2 共通の問い
「Tarazona は表層 SWC が乾いても LE を維持する」現象を **何が説明するか**。当初は深根仮説、最終的には drip wet-bulb 仮説に転換。EC × 衛星 × NDVI フェノロジーの 3 軸から独立検証した。

### 1.3 3 解析の分担

| 解析 | ブランチ | 主軸 | 役割 |
|---|---|---|---|
| **A** | `…LxPGr` | 深層水アクセス / SWC×LE 応答曲線 | 仮説生成・within-site SDS の確立。v4→v14 で 7 回反復 |
| **B** | `…ZnENi` | EC vs 衛星 ET (MOD16 / PML / METv3) | 灌漑バイアスの定量化、τ-fit、SMAP root-zone との照合 |
| **C** | `…iaS0b` | NDVI フェノロジー × フラックス | 位相整合、Oran ローダのバグ発見・修正、共通 `data_loaders.py` 提供 |

### 1.4 仮説の変遷（研究の核 — 改竄禁止）

```
深根仮説 (初期, A v9-v12 で +0.05 SDS が観測 → 一見支持)
   │
   ├─ A v13: 灌漑経過日数で階層化 → 8+d で SDS 負の奇妙な値
   │   └─ Claude が「8+d は実質"春期"」と気づく → 季節交絡
   │
   ├─ A v14: 季節分離 + 夏期内 bucket
   │   └─ TzM 夏期で d0-3 SDS≈0、d4+ で SDS が立ち上がる dose-response
   │   → 深根仮説は表層感受性として再解釈
   │
   ├─ B 分析E (H1/H4/H6):
   │   ├─ H1: τ-based 補正で RMSE −49〜−65%、MBE→0
   │   ├─ H4: bias~days_since_irrig が VPD より圧倒的 (ΔAIC > 60)
   │   └─ H6 multi-stratum: TzM d0-3 で r(SWC, SMAP_rz) = **−0.19**
   │       → **depth inversion = drip wet-bulb の直接観測証拠**
   │
   └─ C H1: Tarazona 灌漑ラグ 1-3d vs 8-14d で EF 27% 減衰 (p=2.93e-7)
       → 深根は寄与あっても優先順位は灌漑効果に劣る
```

**最終的な研究主張 (B 由来、A/C で支持)**:
1. PML / METv3 のバイアスは τ ≈ 4-6 d で過渡的に減衰し補正可能。
2. MOD16 は構造的 floor c = −2.3 mm/d を持つが τ-fit の c 項で吸収できる。
3. SMAP root-zone × in-situ SWC の **depth inversion (r=−0.19 at d0-3)** が drip wet-bulb の直接観測証拠。
4. (補助) rainfed Oran spring で SMAP_rz が in-situ SWC を代替できる (SDS 完全一致 +0.43)。

### 1.5 重要な不変条件（絶対変えるな）
- 解析間で同じ EC データを参照するが、**A v9–v11 と B の初期分析は C が発見した Oran ローダのバグ前**に実行された。Oran の絶対値 (Rn, LE, ET) は再評価対象。
- 単位はクリーンローダ出口で統一: W/m² (flux), mm/day (ET), kPa (VPD), % (SWC)。
- `data_loaders.py` の API は固定 (`load_oran_ec_clean`, `load_tarazona_ec_clean`, `normalize_swc`)。
- `analysis_A_v9.py` は「設定ハブ」(PATHS, SITES, GROWING_MONTHS)。壊さない。
- Oran TIMESTAMP は **必ず三段フォールバック**で読む。pandas デフォルト不可。

---

## 2. 現在までの分析フロー（時系列）

時系列ですべて並べる。`[A]` `[B]` `[C]` は所属。

### Phase 0: 土台
- `[A]` **v4** (ユーザー作成): 半時間→日次集計、closure 補正 (Oran slope=0.7376, TzM=0.7098), 4 象限分類 (normal / soil_dry / atm_dry / compound), diurnal 図。`daily_classified_v4.parquet` と `closure_slopes_v4.json` を出力。

### Phase 1: 深根仮説の初期検定 (A 系列)
- `[A]` **v9** — 3 軸 (SWC × NDWI × VPD) の同時成立日として `deep_access` を抽出。結果 0 日 (Oran/TzM とも)。NDWI 絶対閾値・サイト依存性が原因。
- `[A]` **v10** — DRR (Drought Response Ratio) + NDWI anomaly + lag correlation。cross-site 絶対比較が無意味と判明 (Oran ~3 W/m² vs TzM ~195 W/m²)。
- `[A]` **v11** — within-site SDS / VAS / DSO / CompoundDrop を bootstrap 5000 + 95% CI。初版 `DENOM_FLOOR=5.0` 一律バグで EF/Bowen/ET 全 NaN → 変数別 `DENOM_FLOORS` に修正。
- `[A]` **v12** — 8 図を自己説明的に再設計、休眠期警告、自動判定。**TzM SDS=+0.051 が深根支持に見えたが後にアーティファクトと判明**。

### Phase 2: 仮説の崩壊と再構築 (A 系列)
- `[A]` **v13** — ユーザーから「TzM 灌漑あり」情報。`Irrig_mm` で days_since_irrig を計算しバケット化。8+d で SDS=−0.151、n=249 が異常 → Claude が「8+d は春期=灌漑停止期」と気づく。
- `[A]` **v14**（A 現最終版）— 季節分離 (spring / shoulder / summer) + 夏期内 bucket。
  - Test A: TzM 夏期 LE_corr SDS = 0-1d +0.051, 2-3d +0.004, **4-7d +0.384, 8+d +0.283**（dose-response）
  - Test B: spring Oran SDS=+0.689 (n_n=103, n_s=86), TzM=+0.616 (n_s=15 ⚠ 判定不能), shoulder TzM +0.429, summer TzM +0.272
  - **既知バグ**: `in_band()` が CI 幅・サンプル数を考慮せず spring TzM (n_s=15) で偽陽性「深根支持」と判定 → v15 で fix 予定。

### Phase 3: 衛星検証パイプライン (B 系列)
- `[B]` **分析A (B 内)** — `unify_ec_daily.py` / `add_flags.py` / `aggregate_oran_30min.py` で `ec_daily_master_complete.csv` (1,356 行) を生成。`sds_v14_repro.py` で A v14 を再現:
  - Oran spring SDS=+0.43 [+0.36, +0.51], n=202
  - TzM summer SDS=+0.11 [+0.06, +0.17], n=393 (Oran の約 1/4)
  - TzM summer × bucket: d0-3=+0.13, d4-7=+0.01, d8+=0.00 → 灌漑直後だけ感度
- `[B]` **分析B (B 内)** — GEE 7 プロダクトを `unify_satellite.py` で long→daily。TzM all-year: MOD16 MBE=−2.69, RMSE=3.22；PML MBE=−1.45, RMSE=2.10。Oran ではほぼ 0。
- `[B]` **分析C (B 内)** — 灌漑経過日数の指数減衰モデル。初期は bin 中央値で τ 暴走 → raw daily で fit、τ∈(0,60] 制約。TzM summer × NDVI>0.3 で:
  - MOD16: a=−2.31, τ=4.0d, **c=−2.29 (CI 0 含まず → 構造的 floor)**
  - PML: a=−2.81, τ=4.3d, c=−0.57 (CI 0 含む)
  - METv3: a=−4.03, τ=6.0d, c=−0.62 (CI 0 含む)
- `[B]` **分析D (B 内)** — LSA SAF METv3 (~120k NetCDF) を `load_metv3.py` で日積算。SMAP L4 を 6 km バッファで GEE 抽出。`integrate_metv3_smap.py` で `master_full_v2.csv` (53 列 × 1,356 site-days)。
- `[B]` **分析E (B 内, H1/H4/H6)**:
  - **H1**: τ-based 補正で TzM summer×NDVI>0.3 RMSE: MOD16 4.01→1.39 (−65%), PML 2.82→1.44 (−49%), METv3 3.85→1.50 (−61%); MBE 全て ≈0。
  - **H4**: ΔAIC (vs VPD baseline): MOD16=−73, PML=−66, METv3=−153 (Burnham & Anderson decisive evidence)。
  - **H6 multi-stratum**:
    - Oran_spring (n=203): r(SWC,SMAP_rz)=**+0.80**, SDS_in_situ=SDS_smap=**+0.43**
    - TzM_summer_d0-3 (n=281): r=**−0.19** (depth inversion!)
    - TzM_summer_d4-7 (n=75): r=+0.35
    - TzM_summer_d8+ (n=44): r=+0.61

### Phase 4: フェノロジー検証とローダ修正 (C 系列)
- `[C]` **分析1** — A v9 ローダの Rn 中央値が −63.14 W/m² (物理的不可能) と発覚。`pd.to_datetime` の auto-infer が先頭行 `2018/01/01` から `%Y/%m/%d` を確定し、残り 51,690 行 (時刻付き) を NaT 化していた。914 行→ 真夜中のみ残存。
- `[C]` **分析2** — `data_loaders.py` 新設 (三段 TIMESTAMP パース + センチネル値マスク + 単位統一)。Oran 修正後: 有効日 922, **Rn 中央値 +92.82 W/m²**, **ET 中央値 1.526 mm/day**, EF 有効 28→899, 半時間使用率 1.7%→100%。
- `[C]` **分析3 (H1)** — TzM 灌漑ラグ 1-3d (n=120, EF=0.659) vs 8-14d (n=27, EF=0.484), p=2.93e-7。深根は補助、灌漑優位。
- `[C]` **分析4 (H2)** — NDVI 飽和を EVI で検証。Oran slope 0.369→0.603 (高 NDVI で increase=飽和なし)、TzM 0.833→0.528 (やや減少だが軽微)。飽和は主因ではない。
- `[C]` **分析5 (H4)** — Oran で `partial_r(NDVI, H | ALB)` = −0.280 (単純相関 −0.326 の 86% が ALB 制御後も残存)。気孔開閉が主機構、アルベドではない。ALB は percent 形式と判明 (中央値 13.85)。
- `[C]` **分析6 (H7)** — Oran NDVI セカンドピーク (M04=0.52, M12=0.35) を当初「春+冬の 2 作物」と誤解 → ユーザー訂正で「冬小麦 1 サイクルの main_phase / early_phase」へ。
- `[C]` **分析7** — Oran 2 フェーズ解析: main (Feb–Jul, n=217) vs early (Oct–Jan, n=57)。LE 中央値 45.4 vs 11.8 W/m² (4 倍), EF 0.403 vs 0.417 (同等)。`partial_r(LE,Rn|NDVI)` = +0.645 (main) vs **−0.121 (early)** → early phase は低放射レジーム。
- `[C]` **分析8 (Apr–Jun ベンチマーク)** — 季節を揃えて Oran vs TzM の MWU 比較。**実データ実行ログ未収集 (P★★★)**。
- `[C]` **分析9 (年々変動)** — 2018/19/20 別に EF/LE/ET と H1 灌漑ラグ p 値。**実データ実行ログ未収集 (P★★★)**。
- `[C]` **分析10** — `run_analysis_C.py` / `.sh` / `.env.example` / `RUN_ANALYSIS_C.md` を整備。CLI 引数 > 環境変数 > A v9 既定値 > 標準ロケーションの優先順。

---

## 3. 成功した手法

| # | 手法 | 由来 | 効果 |
|---|---|---|---|
| S1 | **within-site SDS bootstrap (5000 reps, 95% CI)** | A v11 | cross-site 絶対比較の罠を回避、各サイト内 baseline で比較 |
| S2 | **季節 × 灌漑バケットの二重層別** | A v13→v14 | 季節 (spring/shoulder/summer) と days_since_irrig を直交させ artifact を消去 |
| S3 | **closure 補正 (OLS slope で LE/H 補正)** | A v4 | Oran 0.7376, TzM 0.7098 で energy balance closure 適用 |
| S4 | **NDWI/NDVI を月別中央値からのアノマリ化** | A v9→v10, C | サイト依存・季節依存を除去 |
| S5 | **NLS exponential fit on raw daily data** | B 分析C | bin 中央値 fit の τ 暴走を回避。τ ∈ (0, 60] 制約必須 |
| S6 | **AIC によるモデル比較 (ΔAIC > 10 で decisive)** | B 分析E H4 | days_since_irrig が VPD を圧倒 (ΔAIC −66〜−153) |
| S7 | **τ-based 補正 (a·exp(−t/τ) + c)** | B 分析E H1 | RMSE −49〜−65%、MBE ≈ 0 で完全に unbiased 化 |
| S8 | **SMAP multi-stratum 検証 (rainfed control × 灌漑バケット)** | B 分析E H6 | Oran spring n=203 で r=+0.80、TzM d0-3 で r=−0.19 (depth inversion) |
| S9 | **TIMESTAMP 三段パース (明示8種→mixed→Julian復元)** | C 分析1 | Oran 半時間データの 84% データ消失を解決 |
| S10 | **共通 `data_loaders.py` モジュール化** | C 分析2 | 単位統一 + バグ修正を A/B/C で共有 |
| S11 | **partial correlation で交絡変数を制御** | C H4 | アルベド寄与を抜いても気孔機構が 86% 残存と確認 |
| S12 | **REPO 相対パス (`REPO = Path(__file__).parent.parent`)** | B コード設計 | スクリプト再現性。生 CSV は git 管理外 |
| S13 | **xarray lazy load + 2 点 selector** | B 分析D | 120k NetCDF × 164 MB を 3.5 時間で処理 |
| S14 | **bootstrap で τ と CI を同時推定 (N=500)** | B 分析C | CI 重なり判定で結論の頑健性を担保 |

---

## 4. 失敗した手法と理由（再試行禁止）

| # | 試行 | 結果 | 原因 | 対策 | 由来 |
|---|---|---|---|---|---|
| F1 | NDWI 絶対閾値 (>0.0 or >0.109) で deep_access 抽出 | deep_access=0 日 | NDWI 絶対値はサイト・季節依存 | anomaly 化必須 | A v9 |
| F2 | cross-site 絶対 LE 比較 (p=1.5e-47) | 統計有意だが科学的無意味 | 種・LAI・季節すべて違う | within-site baseline 比のみ | A v10 |
| F3 | 季節間プーリング SDS (生育期 1-10 月全部 1 群) | TzM SDS=+0.05 (深根支持に見える artifact) | spring の高 SWC が normal 群、summer の低 SWC が soil_dry 群を支配 | 季節別に評価 | A v12 |
| F4 | SWC×NDWI lag correlation (raw) | 周期成分でフラット | 両方が年周期を持つ | detrend / deseasonalize 必須 | A v9 |
| F5 | GRACE-FO を point sample 利用 | 経度系不一致で全 NaN | lon 0–360 vs −180–180 | 0–360 変換 (−1.86 → 358.14)、ただし 300 km 解像度で参考のみ | A v9 |
| F6 | `in_band` を「CI が 0 を含む → 深根帯」で判定 | spring TzM (n_s=15) で偽陽性 | n と CI 幅を考慮していない | `lo≤0≤hi AND (hi-lo)<0.5 AND n_s≥30` | A v12, v14 |
| F7 | `DENOM_FLOOR=5.0` 一律 | EF/Bowen/ET 全 NaN | スケール違う物理量 | 変数別 `DENOM_FLOORS` + `BOWEN_CLIP=20` | A v11 |
| F8 | τ-fit を bin 中央値で実施 | τ 暴走 (>60), CI ±50 | 5-6 点では NLS の自由度不足 | raw daily (n=300-400) で直接 fit | B 分析C |
| F9 | `cell()` で `LE_Wm2` 無条件優先 | Oran SDS 全 NaN | Oran は半時間集約前で n<30 | `LE_Wm2 が <30 件なら ET_mm にフォールバック` | B 分析A |
| F10 | Oran 30 分ファイルを `pd.to_datetime` 直接パース | 99% NaT | 混合フォーマットで format-lock | year + Julian + Time_hours から再構築 | B / C |
| F11 | SMAP を 200/300 m バッファで GEE 抽出 | 全 NaN | SMAP 9 km grid に対し小さすぎ | **6 km バッファ**へ | B 分析D |
| F12 | GEE JS スクリプトに Unicode 矢印 `→` | SyntaxError | GEE JS パーサーが non-ASCII で停止 | ASCII (`->`) のみ | B 分析D |
| F13 | NetCDF を `xr.open_dataset().load()` で全部メモリ | ~20 TB 相当、実行不可 | eager load | lazy load + 2 点 selector | B 分析D |
| F14 | H6 を Oran summer (n=34) のみで検証 | SDS ≈ 0 で inconclusive | post-harvest で植生不活発 | **multi-stratum 化、Oran spring (n=203) を focal に** | B 分析E |
| F15 | ブランチを取り違えて push | 別ブランチに反映されず | 作業前にブランチ未確認 | `git branch` で確認 → `git checkout` | B |
| F16 | pandas auto-format inference for Oran TIMESTAMP | 84% (51,690 行) NaT | 先頭 1 行から format 確定 | 三段フォールバック (S9) | C 分析1 |
| F17 | `df.where(df > -9000)` を DataFrame 全体に | TypeError | 文字列列に数値比較 | 数値カラムだけ選択 | C 分析2 |
| F18 | SW_IN 単位を固定 (kW/m²) で変換 | 一部データで誤判定 | ratio 形式が混在 | 三段判定 (max<5 / <50 / else) | C 分析5 |
| F19 | `ALB` を fraction として扱う | 中央値 13.85 で異常 | percent 形式だった | 自動検出 /100、フォールバックで SW_OUT/SW_IN | C 分析5 |
| F20 | runner で `if __name__` ブロックを exec | SyntaxError | 直接呼び出し手段なし | `exec(open(script).read(), module.__dict__)` で全体実行 | C 分析10 |
| F21 | Oran NDVI 2 ピークを「春小麦 + 冬小麦」と解釈 | フェノロジー誤解 | 地域品種知識不足 | 「冬小麦 1 サイクル: main / early phase」へ | C 分析6 |
| F22 | H1 を「深根アクセスの証拠」と解釈 | 方向逆 | 深根なら lag↑で EF 維持のはず | 「灌漑優位、深根は補助」へ | C 分析3 |
| F23 | A v12 の SDS=+0.05 を深根支持として引用 | artifact | F3 と同じ季節プーリング | v14 で否定済、**引用禁止** | A v12 |
| F24 | A v13 の自動判定「SDS 逆相関、不自然」 | 解釈誤り | 8+d バケットが実質春期だった | 季節 × バケット交絡 | A v13 |
| F25 | A v14 spring TzM「深根支持」判定 | 偽陽性 | n_s=15 で CI [-0.98, +0.83] | F6 の `in_band` バグ | A v14 |

---

## 5. 現在有力な仮説

### H★1: Drip wet-bulb メカニズム (本筋)
ドリップ灌漑が **深さ 10-30 cm に局所湿潤帯 (wet bulb)** を作り、3-4 日持続する。表層 5 cm SWC とも深部 ~1 m SMAP root-zone とも独立に動く。
- **直接証拠**: B 分析E H6 で TzM d0-3 の r(SWC, SMAP_rz) = **−0.19** (depth inversion)
- **間接証拠**: B 分析C で τ=4-6d, MOD16/PML/METv3 共通
- **再現**: C H1 で EF 1-3d (0.659) → 8-14d (0.484) も同方向

### H★2: 灌漑バイアスは指数減衰で補正可能
- PML / METv3: 過渡的減衰のみ (c CI が 0 を含む)。τ-fit でほぼ完全に補正可能 (B H1: MBE→0)
- MOD16: 構造的 floor c = −2.3 mm/d、FLUXNET2015 校正の作物カバー不足由来と推定
- 3 製品共通の τ ≈ 4-6 d → アルゴリズム固有でなく **dry-surface 駆動の普遍的問題**

### H★3: SMAP root-zone は rainfed 条件で SDS 代替可能
- B H6 Oran spring (n=203, rainfed active): SDS_in_situ = SDS_smap = +0.43, r=+0.80
- → EC タワーがない地域でも SMAP-only で SDS を計算でき、広域マッピング (H7) への道筋

### H4 (補助): 深根アクセスは存在しても寄与は補助的
- C H1 で 7 日に EF 27% 減衰 → 主因は灌漑
- A v14 8+d バケットでも SDS=+0.28 → 表層感受性あり
- 完全否定は不可能（n_s 不足、特に spring TzM）が、headline からは外した

### H5 (Oran 側): NDVI–H 負相関は気孔開閉が主、アルベドは補助
- C H4 partial_r で 86% 残存
- C 分析7 main phase で `partial_r(LE,Rn|NDVI)=+0.645`、early phase で −0.121 → 低放射期は別レジーム

---

## 6. 未解決問題

### U1: A v14 偽陽性バグの修正 (F6, F25)
`in_band()` の判定基準を `(lo≤0≤hi) AND (hi−lo<0.5) AND n_s≥30` に書き直し、A v15 として公開する必要。

### U2: A/B 既存値が Oran ローダバグ前 (C 分析1)
v9–v11 の Oran 絶対値 (Rn, LE, ET) は半時間使用率 1.7% で算出 → **再評価対象**。少なくとも以下を再走:
- A v14 の Oran spring SDS
- B 分析A の Oran spring SDS (+0.43 として論文に載せたもの)
- B H6 の Oran spring n=203, r=+0.80

### U3: A v11–v14 の data_loaders 統合
現在 v10 のみ統合 (commit `2ade97e`)。v11/v12/v13/v14 もパターンが同じなら同様に書き換える。

### U4: B ブランチに `analysis_B_*.py` が見当たらない (C 引き継ぎ §7-4)
解析B の実体は `scripts/` 配下と `hypothesis_tests.py` 等に散在。命名規約と所在の整理が必要。

### U5: C 分析8 (Apr–Jun ベンチマーク) 実データ未実行
`run_analysis_C.py` は整備済。実行ログから p 値・中央値・効果サイズを抽出して `reports/analysis_C_report.md §8.6` に埋める。

### U6: C 分析9 (年々変動 2018/19/20) 実データ未実行
同上、`§8.7` に埋める。

### U7: 論文 Abstract 未完 (B P★★★ 1)
Intro/Methods/Results/Discussion/Limitations は揃った。Abstract のみ未完。

### U8: 図品質チェック (B P★★★ 2)
`fig_H1_correction.png`, `fig_H4_aic.png`, `fig_H6_smap_multistratum.png` のラベル・スケール・色を最終提出形式へ。

### U9: 引用整備 (B P★★★ 3)
Mu 2011, Zhang 2019, Trigo 2018, Reichle 2018, Pettorelli 2005, Burnham & Anderson 2002 等の BibTeX 化。

### U10: 仮説 H2 (灌漑タイプ別 τ) 未検証
flood/sprinkler 灌漑サイトを AmeriFlux 等から取得し τ を比較。

### U11: 仮説 H5 (SIGPAC parcel × Sentinel-2 灌漑面積率) 未検証
METv3 の 5 km ピクセル混合 vs アルゴリズム本質の論点 (B 論点3) と紐づく。

### U12: 仮説 H7 (SDS 広域マッピング) 未着手
B H6 で道筋確認済。FLUXNET2015 / ICOS から多サイトの SMAP root-zone を抽出。

### U13: H8 regional ET correction (Júcar 流域水収支検証) 未着手

### U14: 深根アクセス再評価補助 (C P★5)
TzM の非灌漑日 (lag > 14 d) だけで EF 集計。> 0.4 維持なら深根の存在証拠。

### U15: 多深度 SWC センサーの有無確認 (A P6)
5 cm のみでは drip wet-bulb (20-50 cm) を直接観測できない。既存タワーに 30/50 cm センサーあれば再解析価値大。

---

## 7. 次に試すべき実験（優先順位付き）

### 即実行 (★★★)
1. **`run_analysis_C.py` を実データで実行** → 分析8/9 のログ収集 (U5, U6)。
2. **A v15 = `in_band` バグ修正版** (U1) を書く。同時に Oran 関連数値を新ローダで再走 (U2)。
3. **論文 Abstract 執筆** (U7) — B `paper_methods_results.md` の 1.4 Findings preview / 3.7 / 3.8 / 4.1 を引いて。
4. **B P★★★ 図の最終化** (U8) — fig_H1/H4/H6 のラベル整備。

### 中期 (★★)
5. **A v11–v14 の data_loaders 統合** (U3)。grep で旧ローダ呼び出し検出 → v10 と同じパターンで置換。
6. **U2 の影響度評価** — Oran spring SDS が新ローダでも +0.43 を維持するかを最初に確認。
7. **H2 検証** (U10) — AmeriFlux/ICOS から flood/sprinkler サイト 3–5 個を選定 → τ-fit。
8. **論文 BibTeX 整備** (U9)。
9. **C レポート §8.6/8.7 埋め** (U5, U6 のログから抽出)。

### 長期 (★)
10. **H7 SDS 広域マッピング** (U12) — FLUXNET2015 多サイトで SMAP root-zone から SDS。
11. **H5 SIGPAC parcel × S2 灌漑面積率** (U11)。
12. **Sentinel-1 SAR で灌漑検出** — 灌漑記録がない site への運用展開 (H1 の延長)。
13. **H8 Júcar 流域水収支検証** (U13) — 補正前後を MITECO 統計と照合。
14. **深根再評価補助** (U14) — TzM lag>14d だけで EF 集計。
15. **多深度 SWC センサー確認** (U15)。

---

## 8. 関連コード変更履歴

### A ブランチ (`claude/quantify-water-divergence-LxPGr`)
| コミット | 内容 |
|---|---|
| `e9efc97` | v9 added (3 軸 deep_access, 失敗例として保存) |
| `2069be8` | v10 + requirements (DRR + NDWI anomaly + lag) |
| `9ca8f73` | v11 added (SDS bootstrap, 初版バグあり) |
| `f547c3d` | requirements consolidated |
| `7d90947` | v11 fix (`DENOM_FLOORS` 変数別、`BOWEN_CLIP=20`) |
| `0c12aea` | v12 refactored (自己説明 figures, 自動判定) |
| `bcaf1a3` | v13 irrigation analysis (days_since_irrig バケット) |
| `1dc0ecb` | v14 seasonal × bucket (現 A 最終版) |
| `2ade97e` | data_loaders.py 抽出、v10 を import に書き換え、analysis_C_v1.py 追加 |
| `df16969` | SESSION_SUMMARY.md 追加 |

### B ブランチ (`claude/compare-ec-satellite-et-ZnENi`)
| コミット | 内容 |
|---|---|
| `5064656` | aggregate_oran_30min.py: 混合 datetime parser |
| `6694148` | Oran 30min loader fix |
| `69d4587` | unify_satellite.py 初版 |
| `c4950e0` | merge_satellite_ec.py + killer_figures.py |
| `664f588` | figure_C_summer.py (summer × NDVI フィルタ)、paper_outline.md |
| `14e3829` | Oran SDS NaN bug fix (LE_Wm2→ET_mm フォールバック)、tau_fit.py |
| `2bff988` | load_metv3.py, load_smap.py 初版 |
| `66d3298` | Intro 追加、全 path を REPO 相対に統一、integrate_metv3_smap.py |
| `27641e3` | gee_smap_only.js から Unicode 削除 |
| `b14b584` | load_smap.py を新形式 (SMAP_OranTzM.csv 単一) 対応 |
| `48e8ea2` | unify_satellite.py の BASE を /mnt/hdd/Dataset/ に変更 |
| `5f7f73a` | paper draft を 3 製品 (MOD16, PML, METv3) に拡張 |
| `a58ea9f` | sds_vs_bias.py を 3 製品対応 |
| `4cec8a7` | hypothesis_tests.py + analysis_narrative.md |
| `dedca8e` | H1/H4 結果を narrative と paper に反映 |
| `665b2a8` | H6 を multi-stratum に拡張 (Oran spring focal) |
| `d34c572` | paper Methods/Results/Discussion を 3 製品+H6 multi-stratum に拡張 |
| `108b2fe` | analysis_C_v1.py から未使用旧ローダ import 削除 |

### C ブランチ (`claude/ndvi-flux-analysis-iaS0b`)
| コミット | 内容 |
|---|---|
| (該当複数) | `data_loaders.py` 新設 (Oran TIMESTAMP 三段、センチネル値マスク、単位統一) |
| (該当複数) | `analysis_C_v1.py` 1564 行 (全フェーズ統合): irrigation_lag, ndvi_saturation, albedo_feedback, second_peak, crop_split, same_period, interannual |
| `763f8a8` | run_analysis_C.py / .sh / .env.example / RUN_ANALYSIS_C.md (configurable runners) |
| (該当複数) | scripts/adopt_data_loaders.sh, reports/migration_to_data_loaders.md |
| (該当複数) | reports/analysis_C_report.md, reports/SESSION_SUMMARY.md |

### 統合ブランチ (`claude/analysis-presentation-prep-BQRnM`, **本ブランチ**)
| コミット | 内容 |
|---|---|
| (このコミット) | **統合 SESSION_SUMMARY.md** (本ファイル) |

---

## 9. 重要ファイル一覧

### 9.1 共通ライブラリ（C ブランチ起点）
| パス | 役割 |
|---|---|
| `data_loaders.py` | **A/B/C 共通修正ローダ** (root) — `load_oran_ec_clean`, `load_tarazona_ec_clean`, `normalize_swc` |
| `analysis_A_v9.py` | **設定ハブ** — `PATHS`, `SITES`, `GROWING_MONTHS` の参照元。各ブランチで import される |

### 9.2 解析A 本流（A ブランチ）
| パス | 状態 |
|---|---|
| `analysis_A_v9.py` | NDWI 3軸統合 (失敗例として保存) |
| `analysis_A_v10.py` | DRR + NDWI anomaly (data_loaders 統合済 ✓) |
| `analysis_A_v11.py` | within-site SDS (DENOM_FLOOR fix 済、loader 未統合 ⚠) |
| `analysis_A_v12.py` | 自己説明 figures (loader 未統合 ⚠、結果は artifact) |
| `analysis_A_v13.py` | 灌漑階層化 (loader 未統合 ⚠) |
| `analysis_A_v14.py` | 季節 × bucket (**A 現最終**、loader 未統合 ⚠、`in_band` バグあり) |

### 9.3 解析B 本流（B ブランチ）
| パス | 役割 |
|---|---|
| `scripts/unify_ec_daily.py` | EC 日次統合 |
| `scripts/add_flags.py` | 灌漑バケット / 季節 / 干ばつ class |
| `scripts/aggregate_oran_30min.py` | Oran 半時間 → 日次 |
| `scripts/qc_master.py` | 7 項目 QA |
| `scripts/unify_satellite.py` | GEE 7 プロダクト統合 |
| `scripts/load_metv3.py` | METv3 NetCDF → daily |
| `scripts/load_smap.py` | SMAP CSV → daily (両形式対応) |
| `scripts/merge_satellite_ec.py` | EC + 衛星マージ |
| `scripts/integrate_metv3_smap.py` | METv3+SMAP 統合 |
| `scripts/sds_v14_repro.py` | SDS metric (A v14 再現) |
| `scripts/figure_C_summer.py` | 灌漑バケット boxplot (3 製品) |
| `scripts/sds_vs_bias.py` | SDS vs バイアス scatter (3 製品) |
| `scripts/tau_fit.py` | NLS exponential decay (3 製品) |
| `scripts/hypothesis_tests.py` | H1/H4/H6 検証 (380+ 行) |
| `scripts/gee_extract.js` / `scripts/gee_smap_only.js` | GEE 抽出 |
| `paper_outline.md` / `paper_methods_results.md` / `analysis_narrative.md` | 論文関連 |

### 9.4 解析C 本流（C ブランチ）
| パス | 役割 |
|---|---|
| `analysis_C_v1.py` | 解析C 本体 (1564 行) |
| `run_analysis_C.py` / `.sh` / `.env.example` / `RUN_ANALYSIS_C.md` | ランナー |
| `scripts/adopt_data_loaders.sh` | A/B 向け半自動取り込み |
| `reports/migration_to_data_loaders.md` | 移行手順 |
| `reports/analysis_C_report.md` | 解析C 完全レポート |

### 9.5 中間ファイル（git 管理外、再生成可）
| パス | 由来 | 用途 |
|---|---|---|
| `daily_classified_v4.parquet` | A v4 出力 | A v11–v14 入力 |
| `closure_slopes_v4.json` | A v4 出力 (Oran 0.7376, TzM 0.7098) | A v11–v14 入力 |
| `ec_daily_master.csv` / `ec_daily_master_complete.csv` | B 統合 | B 各分析の入力 |
| `master_full.csv` / **`master_full_v2.csv`** | B 統合 (METv3+SMAP 含む 53 col × 1356 行) | B 各分析の入力 |
| `metv3_daily_all.csv` | B (5,114 site-days, 3.5 h 処理) | B 統合用 |
| `smap_daily.csv` | B (5,112 site-days) | B 統合用 |
| `sds_v14_results.csv` / `sds_vs_bias.csv` / `tau_fit_summary.csv` / `fig_C_summer_summary.csv` / `hypothesis_tests_summary.csv` | B 各分析出力 | 論文用 |

### 9.6 元データ（read-only、絶対パスはユーザー環境依存）
| パス | 内容 |
|---|---|
| `/home/shion-nagamine/Dataset/Eddy data in Spain/Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.csv` | Oran 半時間 EC (xlsx 版もあり) |
| `/home/shion-nagamine/Dataset/Eddy data in Spain/EddyAlmond_Raw5years_withG.xlsx` | TzM 半時間 EC (A v4 用) |
| `/home/shion-nagamine/Dataset/Eddy data in Spain/Daily_Summary_Filtered_forPred_ActEne26.csv` | TzM 日次 (Irrig_mm, Rain_mm, IrrigRain_mm 含む) |
| `/mnt/hdd/Dataset/Sentinel2_NDWI/{Oran,TzM}_NDWI_Export.csv` | A v9/v10 用 |
| `/mnt/hdd/Dataset/MOD13Q1_NDVI_EVI/MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv` | C 用 (AppEEARS) |
| `/mnt/hdd/Dataset/ERA5_2m_Temperature/{year}.nc` (2018–2024) | VPD 統一計算 |
| `/mnt/hdd/Dataset/GRACE-FO_TWL/GRCTellus.JPL.200204_202602.GLO.RL06.3M.MSCNv04CRI.nc` | A v9 用 (参考のみ、300km) |
| `/home/shion-nagamine/Dataset/METv3/YYYY/MM/MMDD/YYYYMMDD_HHMM.nc` | LSA SAF METv3 (~120k files) |
| `/mnt/hdd/Dataset/Fast_OranTzM_{MOD16,PML,LAI,LST,ERA5_Daily,CHIRPS}.csv` | GEE wide CSVs |
| `/mnt/hdd/Dataset/OranTzM_S2_NDVI.csv` | GEE Sentinel-2 |
| `/mnt/hdd/Dataset/SMAP_OranTzM.csv` | SMAP 3-hour (新形式単一ファイル) |

### 9.7 図 (B)
| パス | 役割 |
|---|---|
| `figs/fig_C2_bias_by_irrig_seasonal.png` | Fig 5 headline (灌漑バケット boxplot) |
| `figs/fig_F_tau_fit.png` | Fig 7 decay curve |
| `figs/fig_E_sds_vs_bias.png` | Fig 6 SDS × bias |
| `figs/fig_H1_correction.png` | H1 τ-correction 効果 |
| `figs/fig_H4_aic.png` | H4 ΔAIC 比較 |
| `figs/fig_H6_smap_multistratum.png` | H6 multi-stratum |

---

## 10. 次セッション開始用コンテキスト（Claude 向け）

### 10.1 一行サマリー
> スペイン半乾燥地 2 EC サイト (Oran rainfed cereal / Tarazona drip-irrigated almond) で当初「深根仮説」を検定 → 季節 × 灌漑経過日数で層別したら **drip wet-bulb メカニズム (τ ≈ 4-6 d)** が判明。SMAP root-zone × in-situ SWC の depth inversion (r=−0.19) が直接観測証拠。MOD16/PML/METv3 のバイアスが τ-fit で 49–65% RMSE 削減できることまで確立。

### 10.2 最初に必ずやること
```bash
# 1. このファイル (統合 SESSION_SUMMARY) を必ず最初に読む
cat SESSION_SUMMARY.md

# 2. ブランチ確認
git branch
git log --oneline -20

# 3. 各ブランチの SESSION_SUMMARY (元ファイル) も参照可
git show claude/quantify-water-divergence-LxPGr:SESSION_SUMMARY.md      # A 元
git show claude/compare-ec-satellite-et-ZnENi:SESSION_SUMMARY.md         # B 元
git show claude/ndvi-flux-analysis-iaS0b:reports/SESSION_SUMMARY.md      # C 元

# 4. 中間ファイルの存在確認
ls -la *.csv *.parquet *.json figs/ 2>/dev/null
```

### 10.3 絶対にしてはいけないこと
1. **cross-site の絶対 LE 比較**（F2、p=1.5e-47 で意味なし）
2. **NDWI/NDVI 絶対閾値での deep_access 抽出**（F1）
3. **季節プーリングなしの SDS**（F3、artifact 生成）
4. **A v12 の SDS=+0.05 を深根支持として引用**（F23）
5. **A v13 の「8+d 逆相関」を解釈に使う**（F24、実質春期）
6. **A v14 spring TzM「深根支持」判定を引用**（F25、n_s=15 偽陽性）
7. **pandas `pd.to_datetime` デフォルトを Oran TIMESTAMP に使う**（F16、84% NaT）
8. **τ-fit を bin 中央値で**（F8、τ 暴走）
9. **SMAP を 200/300 m バッファで GEE 抽出**（F11、9 km grid と不整合）
10. **GEE JS に non-ASCII**（F12、パーサーエラー）
11. **NetCDF を `.load()` で eager**（F13、20 TB 相当）
12. **生 CSV (master_full*, ec_daily*) を git add**（容量爆発）
13. **deep-root 仮説を肯定する記述を paper に入れる**（反証済み）
14. **main ブランチに直接 push**

### 10.4 不変条件
- `data_loaders.py` API: `load_oran_ec_clean(filepath, verbose=True)`, `load_tarazona_ec_clean(filepath, verbose=True)`, `normalize_swc(df, site_name)`
- 定数: `EF_DENOM_MIN=10.0`, `SENTINEL_THR=-9000.0`, `VPD_MAX_KPA=10.0`
- 単位: W/m² (flux), mm/day (ET), kPa (VPD), % (SWC)
- Oran TIMESTAMP は三段フォールバック (明示8種 → mixed → year+Julian+Time_hours)
- 全スクリプトで `REPO = Path(__file__).parent.parent`、絶対パスは禁止

### 10.5 主張の階層（論文 3+1 段）
1. **PML / METv3 のバイアスは τ ≈ 4-6 d で過渡的に減衰し補正可能**
2. **MOD16 は構造的 floor c = −2.3 mm/d、τ-fit の c 項で吸収できる (RMSE −65%)**
3. **SMAP root-zone × in-situ SWC の depth inversion (r=−0.19 at d0-3) が drip wet-bulb の直接観測証拠**
4. (補助) **rainfed Oran spring で SMAP_rz が in-situ SWC を代替可能** (SDS 完全一致 +0.43)

### 10.6 ユーザーの好み
- 言語: 日本語ベース、技術内容は英語混じり OK
- 図は **自己説明的** (タイトルに結論、軸ラベル明確、n 表示)
- 解析は **段階的** (各バージョンが前バージョンの問題を 1 つ修正)
- バグ・問題点は **正直に指摘** することを歓迎
- コード: **「なぜそうしたか」をコメント**で残す
- バージョン管理: `analysis_A_vN.py` + `output_analysis_A_vN/` で並走
- 実行環境: Claude のクラウド環境ではデータアクセス不可。ユーザーがローカルで実行 → 結果貼り付け

### 10.7 ブランチ運用
| ブランチ | 用途 | 最終 SHA |
|---|---|---|
| `main` | 使わない (直接 push 禁止) | `e322ba4` |
| `claude/quantify-water-divergence-LxPGr` | 解析A | `2ade97e` |
| `claude/compare-ec-satellite-et-ZnENi` | 解析B | `d34c572` (or `108b2fe`) |
| `claude/ndvi-flux-analysis-iaS0b` | 解析C | `763f8a8` |
| `claude/analysis-presentation-prep-BQRnM` | **本ブランチ** — 統合ドラフト・ポスター準備 | (本コミット) |

### 10.8 次セッションでまず確認すべき項目
1. A v14 の `in_band` バグ修正 (U1) — A v15 を書くか?
2. Oran ローダバグの影響範囲再評価 (U2) — A/B の Oran spring 数値が新ローダで維持されるか?
3. A v11–v14 の data_loaders 統合 (U3) — 一括書き換えするか?
4. C 分析8/9 の実データ実行 (U5, U6) — ユーザー環境でランナー実行
5. 論文 Abstract 執筆 (U7) — B の paper_methods_results.md を起点に
6. ポスター原稿の構成 — 背景 → 仮説変遷 → 手法 (A/B/C) → 試行錯誤 (F1–F25) → 結果 (H★1–3) → 課題 → 自分の意見

### 10.9 ポスター発表用の構造ヒント (本統合ブランチの目的)
- **背景知識**: §1 (作物・サイト・問い)
- **分析手法**: §2 (3 系列の時系列) + §3 (成功手法)
- **分析の試行錯誤**: §4 (F1–F25)
- **分析の流れ・なぜこう分析したか**: §1.4 仮説変遷 + §2 Phase 0–4
- **閾値の決め方**:
  - SWC: 各サイト内 p25 (within-site baseline、cross-site 不可、F2)
  - VPD: median (p50)
  - NDWI/NDVI: 月別中央値からの **anomaly** (絶対閾値は F1 で失敗済)
  - SDS: 季節 × bucket 二重層別 (F3 回避)
  - 灌漑バケット: 0-1d / 2-3d / 4-7d / 8+d (B 分析E では d0-3 / d4-7 / d8+ の 3 区分)
  - τ-fit 制約: τ ∈ (0, 60] (F8 回避)
  - 判定: ΔAIC > 10 (B H4)、CI 幅 < 0.5 AND n ≥ 30 (A v15 で fix 予定)
- **意味と自分の意見**: §5 H★1–3 が論文の主張。drip wet-bulb の depth inversion は新規発見で、衛星 ET 補正の運用に直接効く。

---

**EOF — 統合 SESSION_SUMMARY 完成日 2026-05-11**
*次セッションの Claude へ: §10.3 (絶対禁止 14 項) と §10.4 (不変条件) を必ず守ること。本ファイルが最新の意思決定履歴で、A/B/C 各ブランチの元 SESSION_SUMMARY を上書きしない（A/B/C の SHA は §1 先頭参照）。*
>>>>>>> 49a4546eb7e1cf05160fae654111e7da4ca5ce07
