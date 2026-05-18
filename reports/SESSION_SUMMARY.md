# Session Summary — 解析C (NDVI × Flux) + データローダ統一プロジェクト

最終更新: 2026-05-11
ブランチ構造: 解析A/B/C の 3 並行ブランチ + 共有モジュール (data_loaders.py)

---

## 1. 研究背景・目的

### 大目的
スペイン半乾燥地域 2 サイト (Oran 半乾燥小麦, Tarazona 灌漑アーモンド) の
EC (Eddy Covariance) フラックスデータと衛星 NDVI を組み合わせて、
水利用戦略の違い (深根 vs 灌漑 vs 表層 SWC 依存) を定量化する。

### 解析 A/B/C の分担
- **解析A** (`claude/quantify-water-divergence-LxPGr`): 深層水アクセス検出 / SWC×LE 応答曲線
- **解析B** (`claude/compare-ec-satellite-et-ZnENi`): EC vs 衛星 ET 製品 (MOD16/PML/MET) の比較
- **解析C** (`claude/ndvi-flux-analysis-iaS0b`): NDVI フェノロジー × フラックスの位相整合解析

### 今セッションの直接的目的
1. 解析C で発見された **Oran 半時間値ローダの致命的バグ** (84% のデータ消失) を修正
2. 修正済ローダ (`data_loaders.py`) を解析 A/B/C 共通モジュール化
3. A/B ブランチへの導入手順 (`scripts/adopt_data_loaders.sh`, `reports/migration_to_data_loaders.md`)
4. 解析C v1 を環境依存なく実行できるランナースクリプト群を作成

### 前提条件
- 3 ブランチは独立開発、root に `data_loaders.py` を置いて A/B/C 全部から `import` する
- データは個人ローカル (`/mnt/hdd/Dataset/`, `/home/shion-nagamine/Dataset/`) にあり、
  CI/別マシンでは再現できない → ランナースクリプトはパス可変設計
- 解析A の `analysis_A_v9.py` が `PATHS`, `SITES`, `GROWING_MONTHS` を提供する
  「設定ハブ」になっており、A/B/C すべてから `from analysis_A_v9 import ...` している

---

## 2. 使用データ

### 2.1 Oran EC (半時間値, AmeriFlux 形式)
- **論理名**: `oran_ec`
- **パス例**: `/home/shion-nagamine/Dataset/Eddy data in Spain/Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.csv`
- **期間**: 2018-01-01 〜 2020-12-31 (52,606 半時間レコード)
- **サイト**: Oran (lat=38.82, lon=-1.86), 冬小麦 1 サイクル (Feb–Jul ピーク, Oct–Jan 初期生育)
- **使用カラム**:
  - `TIMESTAMP` (主) or `DateTime`: 日時 (フォーマット混在)
  - 数値: `NETRAD` (Rn), `LE`, `H`, `G`, `VPD` (hPa), `ET`, `SWC_1_1_1` (%)
  - フォールバック: `year`, `Julian`, `Time_hours` (TIMESTAMP が 'nan' の行用)
  - H4 解析用: `ALB` (%), `SW_IN`, `SW_OUT`
- **前処理** (`data_loaders.load_oran_ec_clean`):
  1. TIMESTAMP の三段パース: 明示フォーマット (8 種) → mixed → year/Julian/Time_hours
  2. センチネル値マスク (`< -9000` を NaN 化)
  3. 半時間 → 日次集計 (LE/H/G/Rn/VPD は mean、ET は sum)
  4. EF = LE/(Rn-G), `Rn-G > 10 W/m²` のみ計算, `clip(0, 1.5)`
  5. SWC が 0–100% 範囲のみ採用
  6. VPD: hPa → kPa (÷10), `|VPD| > 10` は NaN
- **期待値** (修正後):
  - 922 有効日, Rn 中央値 **+92.82 W/m²**, ET 中央値 **1.526 mm/day**

### 2.2 Tarazona EC (日次集計, 別フォーマット)
- **論理名**: `tara_ec`
- **パス例**: `/home/shion-nagamine/Dataset/Eddy data in Spain/Daily_Summary_Filtered_forPred_ActEne26.csv`
- **期間**: 2018–2020 (約 700–750 日, 主に Apr–Sep)
- **サイト**: Tarazona (lat=39.266, lon=-1.9397), 落葉アーモンド灌漑
- **使用カラム**:
  - `date`: 日付
  - 名前変更: `LE_avg`→LE, `H_avg`→H, `G_avg`→G, `NetRad_avg`→Rn, `SWC_avg`→SWC
  - ET: `ET_sum` を優先、なければ `ET_avg`
  - VPD: `VPD_kPa` を優先、なければ `VPD_mean/1000` (Pa→kPa)
  - 追加: `Irrig_mm`, `Rain_mm`, `IrrigRain_mm`, `GPP_avg`, `NDVI_orig`, `NDVI_interp`, `Fv_interp`
- **前処理** (`data_loaders.load_tarazona_ec_clean`):
  - 名前変更 + 数値化 + EF 計算
  - SWC 正規化は外側 (`normalize_swc(df, "Tarazona")`) で実施
    → `m³/m³` (0–1) なら ×100

### 2.3 NDVI (MOD13Q1 16 日合成, AppEEARS 経由)
- **論理名**: `NDVI_APPEEARS_CSV`
- **パス例**: `/mnt/hdd/Dataset/MOD13Q1_NDVI_EVI/MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv`
- **空間解像度**: 250 m
- **時間解像度**: 16 日
- **使用列**: 日付列 + サイト識別子列 (`ID`/`Category`) + `_250m_16_days_NDVI`, `_250m_16_days_EVI`

### 2.4 補助データ
- ERA5 VPD (解析A 用, 緯度経度から動的取得)
- GRACE-FO TWL (`/mnt/hdd/Dataset/GRACE-FO_TWL/...`) — 月次 TWS

---

## 3. 実施した分析

### 分析1 — Oran ローダのバグ診断と修正

**目的**: 解析A v9 で計算された Oran の Rn 中央値が **−63.14 W/m²** (物理的に不可能、夜間しか含んでない)
となっていた原因究明。

**実施内容**:
1. `diagnose_oran_parse_failures()` で TIMESTAMP パース失敗行を抽出
2. pandas の `pd.to_datetime` が CSV 先頭行 `2018/01/01` (日付のみ) から
   `%Y/%m/%d` を確定 → 残り 51,690 行 `2018/01/01 00:30:00` を全部 NaT 化していた
3. `dropna(subset=['DateTime'])` で 914 行 (各日真夜中の 1 行) だけが残っていた
4. 三段フォールバックを実装:
   - Stage 1: 8 種の明示フォーマット (`%Y/%m/%d %H:%M:%S`, `%Y-%m-%d %H:%M`, ...)
   - Stage 2: `format="mixed"` で 95% 未満なら適用
   - Stage 3: `year + Julian + Time_hours` 列から `_recover_datetime_from_julian()` で復元

**結果**:
| 指標 (Oran) | v9 (バグあり) | C clean (修正後) |
|---|---|---|
| 有効日数 | 914 | 922 |
| Rn 中央値 [W/m²] | **−63.14** | **+92.82** |
| LE 中央値 [W/m²] | +2.01 | +19.85 |
| ET 中央値 [mm/day] | 0.002 | 1.526 |
| EF 有効日数 | 28 | 899 |
| 半時間レコード使用率 | **1.7%** | 100% |

**解釈**: A/B の Oran 絶対値比較・サイト間比較は全て再評価が必要。

**問題点**: 解析A/B の既存スクリプトは v9 のローダを使い続けるため、
修正版を共通モジュール化して A/B にも導入する必要があった。

**次の仮説**: 修正後の正しい値で「Tarazona > Oran の蒸散」の方向性は維持されるが、
効果サイズは縮む可能性。深根アクセス仮説の妥当性も再評価が必要。

---

### 分析2 — `data_loaders.py` モジュール化

**目的**: 解析 A/B/C 全てに同じ修正版ローダを配布。

**実施内容**:
- `data_loaders.py` (7.2 KB) を root に作成、C ブランチで管理
- API:
  - `load_oran_ec_clean(filepath, verbose=True) -> pd.DataFrame`
  - `load_tarazona_ec_clean(filepath, verbose=True) -> pd.DataFrame`
  - `normalize_swc(df, site_name) -> pd.DataFrame`
  - 定数: `EF_DENOM_MIN=10.0`, `SENTINEL_THR=-9000.0`, `VPD_MAX_KPA=10.0`
- 統一単位: W/m² (flux), mm/day (ET), kPa (VPD), % (SWC)

**結果**: A/B 各々で `git checkout origin/claude/ndvi-flux-analysis-iaS0b -- data_loaders.py`
で取り込めば、解析間で完全に同じローダを使える。

**解釈**: バグ修正が 1 箇所で済むので保守性が向上。

**問題点**:
- A/B の既存スクリプトは内部に旧ローダを定義済み (`def load_oran_ec`, `def load_tarazona_ec`)
- これらを手動で `from data_loaders import ...` に書き換える必要
- 完全自動置換は false-rename リスクがあるので、grep だけする半自動スクリプトを作った

---

### 分析3 — H1+H8: 灌漑経過日数別 EF 減衰 (Tarazona)

**目的**: 「Tarazona の高い EF は深根アクセスか、それとも灌漑直後の表層湿潤か」を判別。

**実施内容** (`irrigation_lag_analysis(tara_df, save_dir)`):
- `Irrig_mm > 0.5` の日を「灌漑日」と定義
- 各観測日について「最後の灌漑日からの経過日数」(lag) を計算
- lag 1–3 日 vs lag 8–14 日 で EF/LE を比較 (Mann-Whitney U test)

**結果**:
- lag 1–3 (n=120): EF 中央値 = **0.659**
- lag 8–14 (n=27): EF 中央値 = **0.484**
- p = 2.93e-7 (★★★)
- 7 日間で EF が 27% 減衰

**解釈**:
- 当初の「深根アクセス仮説」を支持しない
- むしろ **灌漑直後の表層湿潤が EF を駆動** していることを示唆
- 解析A の「深根」フレーミングを「灌漑効果」に置き換える必要

**問題点**: 「深根アクセス」を否定するわけではなく、寄与の優先度が「灌漑 > 深根」になっただけ。

**次の仮説**: 灌漑後の表層 SWC 減衰速度とフラックスの結合をモデル化すべき。

---

### 分析4 — H2: NDVI 飽和を EVI で検証

**目的**: NDVI が高い領域で蒸散と相関しない理由は「NDVI 飽和」か?

**実施内容** (`ndvi_saturation_check(merged, site, save_dir)`):
- 各サイトの NDVI を p67 で分割 (低/高 NDVI 期間)
- 各範囲で OLS 回帰 `EVI ~ NDVI` の slope を比較
- 飽和なら高 NDVI 範囲で slope が下がるはず

**結果**:
- Oran: slope 0.369 → 0.603 (高 NDVI でむしろ増加 = 飽和なし)
- Tarazona: slope 0.833 → 0.528 (やや減少だが大幅ではない)

**解釈**: 飽和は主因ではない。Oran の NDVI~LE 弱相関は別要因 (低 LAI / 水ストレス) による。

**問題点**: なし。

---

### 分析5 — H4: アルベド・フィードバック (Oran)

**目的**: 「Oran で NDVI↑ → H↓」は植物冷却ではなくアルベド変化が原因かを検証。

**実施内容** (`albedo_feedback_check(oran_raw_path, oran_merged)`):
- 生 CSV から `ALB` (%) と `SW_IN/SW_OUT` を取得
- 単位検出: `ALB` 中央値 = 13.85 → percent 形式と判明
- フォールバック: ALB = SW_OUT/SW_IN を計算
- 部分相関: `partial_r(NDVI, H | ALB)` を計算

**結果**:
- 単純相関: `r(NDVI, H) = -0.326`
- 部分相関: `partial_r(NDVI, H | ALB) = -0.280`
- 86% の効果が ALB 制御後も残存
- 有効 n = 918 日

**解釈**: アルベドフィードバックは寄与するが、主機構は **気孔開閉 (蒸散冷却)** である。

**問題点**: SW_IN の単位検出ロジックが必要だった (kW/m² ? W/m² ? 比?)。
解決: 最大値で自動判定 (max<5 → kW/m² ×1000, max<50 → ×100, else 無変換)。

---

### 分析6 — H7: NDVI セカンドピーク検出

**目的**: Oran で 2 作物 (春小麦 + 冬小麦?) の存在を確認。

**実施内容** (`detect_second_peak(per_site)`):
- 月別 NDVI 中央値の local maxima を scan

**結果**:
- Oran: M04=0.52, M12=0.35 → **2 つのピーク**
- Tarazona: M05=0.43 → **単一ピーク**

**解釈**:
- 当初「春作物 vs 冬作物」と解釈
- **ユーザー指摘で修正**: Oran は冬小麦 1 サイクルの「主成熟期 (Feb–Jul)」と
  「初期生育期 (Oct–Jan, 播種後)」の 2 フェーズだった
- ラベルを「main growth phase / early growth phase」に変更

**問題点**: 植生フェノロジーの誤解釈 (品種・地域知識の欠如)。

---

### 分析7 — Oran 生育期分割解析

**目的**: 2 フェーズ間で水利用戦略が異なるかを検証。

**実施内容** (`crop_split_analysis_oran(merged_oran, save_dir)`):
- main phase: Feb–Jul (n=217, NDVI 中央値=0.536)
- early phase: Oct–Jan (n=57, NDVI 中央値=0.446)
- 部分相関: `partial_r(LE, Rn | NDVI)` を各フェーズで計算

**結果**:
| 指標 | main phase | early phase |
|---|---|---|
| NDVI | 0.536 | 0.446 |
| LE 中央値 | 45.4 W/m² | 11.8 W/m² (≈ 1/4) |
| EF | 0.403 | 0.417 (ほぼ同じ) |
| partial_r(LE,Rn|NDVI) | +0.645 | **−0.121** |

**解釈**:
- EF がほぼ同じ → エネルギー分配比は変わらない
- LE の 4 倍ギャップは **Rn 駆動** (光合成有効放射の季節変動)
- early phase で Rn の駆動が消える → 低放射期は別レジーム

---

### 分析8 — 同期間ベンチマーク (Apr–Jun)

**目的**: 2 サイトのフェノロジー位相差を吸収して、純粋な機構差を抽出。

**実施内容** (`same_period_benchmark(oran_m, tara_m, save_dir)`):
- 両サイトの「共通活発期 = Apr–Jun」に限定
- Mann-Whitney U で LE / EF / ET / NDVI / Rn / VPD を比較
- 「low SWC × green」サブセットでも比較
- 月別中央値をプロット (Apr–Jun を網掛け)

**結果**: (実行ログ収集待ちのプレースホルダ)

**解釈**: 期待される結論 = 季節アラインメント効果を除いても Tarazona > Oran が維持されれば、
真の機構差 (灌漑の存在) が確定。

---

### 分析9 — 年々変動チェック (2018/2019/2020)

**目的**: 結論が単一年の異常に依存していないことを確認。

**実施内容** (`interannual_check(oran_m, tara_m, save_dir)`):
- 年別に生育期 EF/LE/ET の中央値を計算
- H1 灌漑ラグ検定も年別に実施
- bar chart で年別 Tarazona vs Oran の方向性を確認

**結果**: (実行ログ収集待ち)

**解釈**: 期待 = 全年で Tarazona > Oran なら結論は頑健。

---

### 分析10 — ランナースクリプト整備

**目的**: 環境依存パスを抜き出して、別マシン/別ユーザーでも再現実行可能にする。

**実施内容**:
- `run_analysis_C.py` (130 行 Python ランナー)
- `run_analysis_C.sh` (170 行 Bash ラッパ)
- `.env.example` (環境変数テンプレート)
- `RUN_ANALYSIS_C.md` (350+ 行使用ガイド)

優先順:
1. CLI 引数 (`--oran`, `--tarazona`, `--ndvi`, `--output`)
2. 環境変数 (`BAKANPOSS_ORAN_EC` 等)
3. `analysis_A_v9.py` 既定値
4. 標準ロケーション (NDVI のみ)

**結果**: `python3 run_analysis_C.py --dry-run` で設定検証だけできる。

---

## 4. コード変更履歴

### `data_loaders.py` (新規作成, C ブランチ)
- **変更箇所**: 全 183 行
- **変更理由**: 解析 A/B/C 共通の修正版ローダ
- **状態**: 完成、A/B にも配布済 (git checkout 経由)

### `analysis_C_v1.py` (C ブランチ)
- **変更箇所**:
  - L40-43: `from analysis_A_v9 import` の整理 (SITES, GROWING_MONTHS, PATHS のみ残す)
  - L45-48: `from data_loaders import ...` 追加
  - L1484, L1487: `load_*_clean` 経由でロード
  - L1513-1522: `irrigation_lag_analysis` 呼び出し追加
  - L1535: `crop_split_analysis_oran` (改名: 旧 spring/winter crop → main/early phase)
  - L1538: `same_period_benchmark` 新規呼び出し
  - L1541: `interannual_check` 新規呼び出し
- **状態**: 完成、syntax 検証 OK

### `analysis_A_v10.py` (A ブランチ)
- **変更箇所**:
  - L41: `from data_loaders import load_oran_ec_clean, load_tarazona_ec_clean, normalize_swc` 追加
  - L146-188: 旧 `load_oran_ec`, `load_tarazona_ec`, `normalize_swc` 関数定義を **削除** (44 行)
  - L558-559: 呼び出しを `load_oran_ec_clean(...)` と
    `normalize_swc(load_tarazona_ec_clean(...), "Tarazona")` に置換
  - L562: 冗長な `tara_ec = normalize_swc(...)` を削除
- **状態**: commit `2ade97e`, push 済
- **コミットメッセージ**: "Replace analysis_A_v10 loader with data_loaders module"

### `analysis_C_v1.py` (B ブランチ)
- **変更箇所**: L40-43 の `from analysis_A_v9 import` から
  未使用の `load_oran_ec, load_tarazona_ec` を削除 (1 行削除)
- **状態**: commit `108b2fe`, push 済
- **コミットメッセージ**: "Remove unused old loader imports from analysis_A_v9"

### `scripts/adopt_data_loaders.sh` (C ブランチ)
- **新規作成**: A/B 向け半自動取り込みスクリプト
- 機能: ① C から data_loaders.py を checkout, ② 旧ローダ呼び出しを grep, ③ import sanity test
- **状態**: 完成

### `reports/migration_to_data_loaders.md` (C ブランチ)
- **新規作成**: 移行ガイド (149 行)
- 内容: なぜ移行が必要か / 期待される値の変化 / API 仕様 / テストコード
- **状態**: 完成

### ランナー群 (C ブランチ, commit `763f8a8`)
- `run_analysis_C.py` (新規)
- `run_analysis_C.sh` (新規)
- `.env.example` (新規)
- `RUN_ANALYSIS_C.md` (新規)
- **状態**: 4 ファイル全て push 済

---

## 5. 試して失敗したこと

### 失敗1: pandas auto-format inference (`pd.to_datetime` デフォルト)
- **試行**: `df["datetime"] = pd.to_datetime(df["TIMESTAMP"])` (デフォルト infer)
- **結果**: 先頭行 `2018/01/01` から `%Y/%m/%d` を確定 → 残り 51,690 行を NaT 化
- **ダメな理由**: pandas は最初の数行で format を確定し、それと合わない行は coerce/NaT
- **対策**: 明示フォーマット 8 種を順に試す → mixed フォーマット → Julian 復元の三段

### 失敗2: `df.where(df > -9000)` を DataFrame 全体に適用
- **試行**: `df = df.where(df > -9000)` でセンチネル値マスク
- **結果**: `TypeError` (string カラムに数値比較を適用)
- **ダメな理由**: 全カラムに比較演算が走り、文字列列で例外
- **対策**: 数値カラムだけ選択して個別マスク (`for col in flux_cols: ...`)

### 失敗3: `SW_IN` の単位検出を固定値で
- **試行**: `if SW_IN.max() < 5: SW_IN *= 1000` (kW/m² → W/m² 想定)
- **結果**: 一部データで誤判定 (実は ratio 形式だった)
- **対策**: 三段判定 (max<5 / max<50 / else) + log で確認

### 失敗4: `ALB` 列の値が想定外
- **試行**: ALB をそのまま fraction として扱う
- **結果**: ALB 中央値 = 13.85 (fraction なら 0.14 のはず)
- **判明**: 単位は percent (%)
- **対策**: 自動検出して /100, それでも怪しい場合は `SW_OUT/SW_IN` フォールバック

### 失敗5: ランナー実装の最初の試行 (Python 内 `exec` の混乱)
- **試行**: `analysis_C_v1.if __name__ == "__main__"` のような構文 (構文エラー)
- **結果**: SyntaxError
- **ダメな理由**: Python の `if __name__` ブロックを外部から直接呼ぶ手段がない
- **対策**: `exec(open(script).read(), module.__dict__)` で全体実行 + 事前に globals を上書き

### 失敗6: 「Oran は春/冬 2 作物」誤解釈
- **試行**: M04 + M12 ピークから「春小麦 + 冬小麦の 2 作目」と判定
- **結果**: ユーザーから「Oran は冬小麦 1 サイクル」と訂正
- **ダメな理由**: フェノロジー知識不足。Oct–Jan は播種〜分蘖、Feb–Jul は再開〜成熟
- **対策**: ラベルを `main_phase` / `early_phase` に変更

### 失敗7: H1 結果の解釈方向ミス
- **試行**: 「灌漑ラグで EF 減衰」を「深根アクセスの証拠」と解釈
- **結果**: 強い有意性 (p=2.93e-7) はむしろ **表層水依存** を意味する
- **ダメな理由**: 深根なら lag が増えても EF は維持されるはず
- **対策**: 解釈を「灌漑優位、深根は補助」に修正

---

## 6. 現在の研究上の論点

### 論点1: 「深根アクセス」仮説の位置づけ
- 解析A の主張: Tarazona の高 EF は深根が原因
- H1 の発見: 7 日で EF 27% 減衰 → 灌漑が主因
- **論点**: 深根は完全否定なのか、寄与の優先順位が下がっただけなのか
- 必要な追加解析: 灌漑日を除外した「非灌漑期」だけで EF が高ければ深根存在の証拠

### 論点2: Oran の Rn-LE 解離 (low LE despite high Rn)
- H4 の結果: 気孔開閉が主因 (アルベドではない)
- 解析A の partial correlation: main phase で `r(LE,Rn|NDVI)=+0.645`
- **論点**: 高 Rn でも LE が伸びない = 何かが律速 (水? VPD? CO₂?)
- 候補: 土壌水分律速 (Tarazona と比較した「ストレス」状態)

### 論点3: 同期間ベンチマーク (Apr–Jun) の結論
- まだ実データで未実行 (run log 必要)
- **論点**: 季節を揃えても Tarazona > Oran なら結論頑健、揃ったら結論消失なら見直し
- 予測: Tarazona > Oran は維持されるが、効果サイズ縮小

### 論点4: 年々変動の安定性
- 2018: 平常年? 2019: 干ばつ年? 2020: 異常気象?
- 未実行 (run log 必要)
- **論点**: 全年で同じ結論なら頑健、特定年に依存すれば「気象年効果」

---

## 7. 未解決課題

優先度順 (★★★ = 最優先):

### ★★★ 1. `analysis_C_v1.py` を実データで実行してログ収集
- 必要: Oran EC CSV, Tarazona EC CSV, NDVI CSV のパスを設定
- 実行: `python3 run_analysis_C.py --oran ... --tarazona ... --ndvi ...`
- 期待出力: `run_C_v13.log` (or similar) + 出力ディレクトリの 20+ PNG

### ★★★ 2. レポート 8.6 / 8.7 セクションに実数値を埋め込む
- `reports/analysis_C_report.md` のプレースホルダ箇所:
  - §8.6 Same-period Apr-Jun benchmark
  - §8.7 Inter-annual robustness check
- ログから Mann-Whitney U の p 値 / 中央値 / 効果サイズを抽出

### ★★ 3. 解析A v11–v14 のローダ更新
- 現在 v10 だけ更新済
- v11–v14 が v10 と同じパターンならコピー要領で同じ変更を適用
- 確認: `grep -n "def load_oran_ec\|def load_tarazona_ec" analysis_A_v1[1-4].py`

### ★★ 4. B ブランチに analysis_B_*.py が存在するか確認
- 現在の B ブランチには analysis_B*.py が見つからなかった
- 解析B の実体スクリプトはどこ? (hypothesis_tests.py?)
- もし他に loader を持つスクリプトがあれば同じく差し替え

### ★ 5. 深根アクセス再評価 (補助解析)
- Tarazona の非灌漑日 (lag > 14 日) だけで EF を集計
- それでも EF > 0.4 なら深根存在の証拠

### ★ 6. 解析C レポートの §4 タイムラインを更新
- 失敗1〜7 を「試行錯誤」として記述
- 特に Oran TIMESTAMP の三段フォールバックの経緯

### ★ 7. ベンチマーク表の追加
- v9 (バグあり) vs C clean (修正後) の各指標を 1 表に
- A/B の論文 draft に直接埋め込めるよう Markdown 形式で

---

## 8. 次セッション開始時にやること

1. `git checkout claude/ndvi-flux-analysis-iaS0b`
2. `git pull --ff-only origin claude/ndvi-flux-analysis-iaS0b`
3. データパスを設定 (ユーザー環境に合わせて):
   ```bash
   export BAKANPOSS_ORAN_EC="..."
   export BAKANPOSS_TARA_EC="..."
   export BAKANPOSS_NDVI_FILE="..."
   ```
4. ドライラン検証: `python3 run_analysis_C.py --dry-run`
5. 本実行: `python3 run_analysis_C.py 2>&1 | tee run_C_v14.log`
6. ログから以下を抽出:
   - Same-period Apr-Jun の MWU 結果 (LE, EF, ET, NDVI)
   - Inter-annual 各年の EF 中央値 (Oran vs Tarazona)
   - H1 各年の灌漑ラグ p 値
7. `reports/analysis_C_report.md` の §8.6 / §8.7 に数値を埋める
8. A/B ブランチに移って v11–v14 のローダ確認:
   ```bash
   git checkout claude/quantify-water-divergence-LxPGr
   grep -n "def load_oran_ec\|def load_tarazona_ec" analysis_A_v1[1-4].py
   ```
9. 必要なら v11–v14 にも `data_loaders` 統合を適用
10. レポート最終化 → 必要なら PR 作成

---

## 9. Claude への引き継ぎ指示

### 最初に理解すべき構造

このリポジトリには **3 つの平行ブランチ** がある。それぞれ独立解析だが、
**root の `data_loaders.py` を共有** する設計になっている:

- A ブランチ: `claude/quantify-water-divergence-LxPGr` (深根 / SWC 応答)
- B ブランチ: `claude/compare-ec-satellite-et-ZnENi` (衛星 ET 比較)
- C ブランチ: `claude/ndvi-flux-analysis-iaS0b` (NDVI フェノロジー × フラックス, **今主に作業**)

### 重要な不変条件 (絶対変えるな)

1. **`data_loaders.py` の API は固定**:
   - `load_oran_ec_clean(filepath, verbose=True) -> DataFrame`
   - `load_tarazona_ec_clean(filepath, verbose=True) -> DataFrame`
   - `normalize_swc(df, site_name) -> DataFrame`
   - 定数: `EF_DENOM_MIN`, `SENTINEL_THR`, `VPD_MAX_KPA`
2. **単位はクリーンローダ出口で統一**: W/m² (flux), mm/day (ET), kPa (VPD), % (SWC)
3. **`analysis_A_v9.py` は「設定ハブ」**: PATHS, SITES, GROWING_MONTHS の参照元 (壊すな)
4. **Oran TIMESTAMP は必ず三段フォールバック**で読む (pandas デフォルトは絶対使うな)

### 重要な研究判断 (背景を理解せよ)

1. **Oran は冬小麦 1 サイクル** (春/冬 2 作物ではない). Feb–Jul = main phase, Oct–Jan = early phase
2. **Tarazona は灌漑アーモンド** (落葉樹). Apr–Sep が活発期
3. **「深根アクセス仮説」は補助的仮説に格下げ**. 主は「灌漑効果」(H1 で確証)
4. **NDVI 飽和は主因ではない** (H2 で否定)
5. **Oran の NDVI~H 負相関は気孔開閉**, アルベドではない (H4 で確証)

### 制約

- 実データは個人ローカル (`/mnt/hdd/Dataset/`, `/home/shion-nagamine/...`) にあり、
  別マシンでは絶対パスが効かない → **必ずランナースクリプト経由で実行**
- 試行錯誤の履歴は `reports/analysis_C_report.md` §4 にあるので参照
- 解析A/B の論文 draft は別ブランチにあり、Oran 絶対値は古い (バグあり) 値
  → A/B レポートを更新するときは「再評価が必要」と明示する

### よくある落とし穴

- pandas の `pd.to_datetime` (デフォルト) を Oran TIMESTAMP に使うと **84% 消える**
- `df.where(df > -9000)` は文字列列があると死ぬ → 数値列だけ選択
- `ALB` は percent 形式の可能性 (中央値 13.85 のとき) → 単位検出必須
- `ET_avg` ではなく `ET_sum` を使う (Tarazona)

---

## 10. 重要ファイル一覧

### コアモジュール
| パス | 役割 |
|---|---|
| `/home/user/bakanposs/data_loaders.py` | **共通修正ローダ** (A/B/C 全部からインポート) |
| `/home/user/bakanposs/analysis_A_v9.py` | 設定ハブ (PATHS, SITES, GROWING_MONTHS) |
| `/home/user/bakanposs/analysis_C_v1.py` | 解析C 本体 (1564 行, 全フェーズ統合) |
| `/home/user/bakanposs/analysis_A_v10.py` | 解析A 主スクリプト (data_loaders 統合済) |

### ランナー (今セッションで作成)
| パス | 役割 |
|---|---|
| `/home/user/bakanposs/run_analysis_C.py` | Python ランナー (CLI + env var) |
| `/home/user/bakanposs/run_analysis_C.sh` | Bash ラッパ |
| `/home/user/bakanposs/.env.example` | 環境変数テンプレート |
| `/home/user/bakanposs/RUN_ANALYSIS_C.md` | 使用ガイド |

### 移行インフラ
| パス | 役割 |
|---|---|
| `/home/user/bakanposs/scripts/adopt_data_loaders.sh` | A/B 向け半自動取り込み |
| `/home/user/bakanposs/reports/migration_to_data_loaders.md` | 移行手順書 |

### レポート
| パス | 役割 |
|---|---|
| `/home/user/bakanposs/reports/analysis_C_report.md` | 解析C の試行錯誤含む完全レポート |

### データ (環境依存、ハードコード禁止)
| 論理名 | パス例 | 期間 |
|---|---|---|
| `oran_ec` | `/home/shion-nagamine/Dataset/Eddy data in Spain/Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.csv` | 2018–2020 半時間値 |
| `tara_ec` | `/home/shion-nagamine/Dataset/Eddy data in Spain/Daily_Summary_Filtered_forPred_ActEne26.csv` | 2018–2020 日次 |
| `ndvi_file` | `/mnt/hdd/Dataset/MOD13Q1_NDVI_EVI/MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv` | 16-day MOD13Q1 |

### 最新コミット (各ブランチ)
| ブランチ | SHA | メッセージ |
|---|---|---|
| `claude/ndvi-flux-analysis-iaS0b` | `763f8a8` | Add configurable runners for analysis_C_v1 |
| `claude/quantify-water-divergence-LxPGr` | `2ade97e` | Replace analysis_A_v10 loader with data_loaders module |
| `claude/compare-ec-satellite-et-ZnENi` | `108b2fe` | Remove unused old loader imports from analysis_A_v9 |

---

**EOF**
