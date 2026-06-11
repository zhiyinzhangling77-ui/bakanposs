# Analysis C — NDVI Phenology Reinforcement (PLAN/v1 in progress)

> **Status**: 🔶 v1 in progress (`analysis_C_v1.py` 存在)
> **Role**: 解析A の結論を **NDVI(植被活性)** という別物理量で補強
> **Why**: τ ≈ 3-4d が "actively growing" 状態の特性なら、NDVI で active 状態を客観定義 → 解析A の主張がさらに堅牢化

---

## 1. 目的と研究全体での位置

### 1.1 解析C が答える問い
- **Q1**: 解析A の "active growing" 月定義(Oran Nov-Jun, Tarazona Jun-Sep)は NDVI と整合するか?
- **Q2**: NDVI ピーク時の τ は ≈ 3-4 d か?(生理学的最も active な時期)
- **Q3**: NDVI × Rn(簡易 GPP プロキシ)で EC LE を予測できるか?
- **Q4**: Tarazona の "deep_access" 候補日に NDVI は健康か?(植被水分維持の独立証拠)

### 1.2 解析A への補強関係
```
解析A: EC LE で 'τ ≈ 3-4d, active growing period'
   ↓ "active growing" を月で定義(Nov-Jun, Jun-Sep)
   ↓
解析C: NDVI で active を客観定義
   - 仮定 (active months) と NDVI 観測が整合 → 解析A 強化
   - 不整合 → NDVI-based 再定義で解析A を update
```

### 1.3 論文化への含意
- NDVI 整合性 → 「**vegetation activity ground-truth で fit period を validate**」 と Methods に記述
- NDVI×Rn vs LE 高相関 → 「光合成プロキシで τ メカニズム生理学的に説明可能」

---

## 2. 使用するデータ

### 2.1 衛星 NDVI プロダクト
| プロダクト | 解像度 | 期間 | 推奨度 |
|---|---|---|---|
| **MODIS MOD13Q1** | **250m / 16-day** | 2000-現在 | ★★★ メイン |
| Sentinel-2 NDVI | 10m / 5-day | 2017-現在 | ★★ 高解像度補強 |
| Landsat 8/9 NDVI | 30m / 16-day | 2013-現在 | ★ 補助 |

### 2.2 既に取得済(`analysis_C_v1.py` で言及)
- MODIS MOD13Q1: scripts/gee_extract.js で抽出予定
- パス未確認 → 次セッションで確認必要

---

## 3. 既存実装(`analysis_C_v1.py`)の構造

```
目的1: 生育期/非生育期を NDVI から客観定義
   → 解析A の手動 GROWING_MONTHS フィルタの validate
目的2: NDVI vs EF / LE / GPP_proxy を Oran vs Tarazona で比較
目的3: サイト間フェノロジー位相の定量化
        (ピーク NDVI 日、立ち上がり/枯れ落ちの傾き)
目的4 [A 連結]: 低 SWC × 高 NDVI 期で Tarazona が EF を保てるか
         = 深根アクセスの状況証拠の再評価
目的5 [B 連結]: EC LE と NDVI×Rn(簡易 GPP プロキシ)の整合性チェック
```

→ v1 は **多目的すぎる**。解析A 完成を踏まえ、**目的1+4 に焦点**を絞った v2 が望ましい

---

## 4. v2 への設計提案(解析A 統合版)

### 4.1 簡素化された目的
1. **NDVI で active period 客観定義** → 解析A の Nov-Jun/Jun-Sep を validate
2. **NDVI peak ± Xday window** で τ 再計算 → 解析A 最終結果を強化
3. **NDVI dynamics と EC LE の時間関係**(NDVI lag/lead を見る)

### 4.2 実装テンプレート

```python
"""
analysis_C_v2_ndvi_phenology.py

NDVI で active growing 期間を客観定義し、解析A 結論を強化
"""

# === 1. NDVI 読込 ===
def load_modis_ndvi(filepath, site):
    """MOD13Q1 16-day NDVI を読込み、daily に線形補間"""
    df = pd.read_csv(filepath)
    df["date"] = pd.to_datetime(df["date"])
    # quality filter (VI Quality flags)
    df = df[df["DetailedQA"] < some_threshold]
    return df

# === 2. NDVI 季節曲線抽出 ===
def fit_phenology_curve(ndvi_df, method="logistic"):
    """
    Savitzky-Golay smoothing + logistic fit で
    生育期開始日 (SOS), ピーク日, 終了日 (EOS) を抽出
    """
    ...

# === 3. NDVI-defined active period ===
def define_active_period(pheno, ndvi_threshold=0.4):
    """NDVI > threshold の期間を 'active growing' と定義"""
    ...

# === 4. 解析A 結果と比較 ===
def compare_with_analysis_a():
    """
    Oran: NDVI-active vs Nov-Jun の overlap %
    Tarazona: NDVI-active vs Jun-Sep の overlap %
    """
    ...

# === 5. NDVI peak window で τ 再計算 ===
def tau_at_ndvi_peak(events, ndvi_df, window_days=30):
    """NDVI peak ± 30 日に event を限定して τ fit"""
    # 解析A v27 の fit_with_validation を流用
    ...
```

### 4.3 想定される出力

| ファイル | 内容 |
|---|---|
| `v2_phenology_metrics.csv` | サイト × 年 の SOS/Peak/EOS dates |
| `v2_active_period_validation.csv` | 解析A の active months 定義 validation |
| `v2_ndvi_peak_tau.csv` | NDVI peak window 限定 τ |
| `fig01_ndvi_seasonal_curve.png` | 両サイトの年次 NDVI cycle |
| `fig02_active_period_overlap.png` | 仮定 vs 観測 active period |
| `fig03_ndvi_peak_tau.png` | NDVI peak で τ ≈ 3-4d 確認 |

---

## 5. 期待される結果と論文記述案

### 5.1 Scenario A: NDVI が解析A の active months を validate
> "We verified that the periods used for τ estimation correspond to
> active vegetation states. MODIS MOD13Q1 NDVI showed elevated values
> (>0.4) during the assumed active periods (Nov-Jun for Oran cereal,
> Jun-Sep for Tarazona almond), confirming our phenology-matched
> comparison. Within the NDVI-peak window (±30 days from peak), the
> τ estimates were X.X d for Oran and Y.Y d for Tarazona, consistent
> with the pooled results."

### 5.2 Scenario B: NDVI が解析A 仮定とずれる
> "MODIS NDVI peaked in [actual peak month] for Oran cereal, slightly
> later than the assumed active-period midpoint. Recomputing τ on
> NDVI-defined active periods yielded similar values..."

### 5.3 NDVI × Rn vs LE
> "The product of NDVI and net radiation explained X% of EC LE
> variance at both sites, consistent with the canopy-conductance
> control on ET. The product's seasonal cycle..."

---

## 6. 解析B との並行進行関係

```
解析B (衛星 ET 検証)         解析C (NDVI 補強)
  ↓                            ↓
  τ_satellite                  active period 客観定義
  amplitude ratio              NDVI peak での τ
  ↓                            ↓
       両者統合 → 解析A 主結果の補強
```

### 並行の利点
- 同じ MODIS データ(MOD16 + MOD13)を一度に extract
- GEE script を共通化可能
- 同じ event detection / fit ロジックを使う

---

## 7. データ取得チェックリスト

### Phase 1: MOD13Q1 NDVI 抽出
- [ ] GEE script で Tarazona NDVI 時系列(2018-2024)
- [ ] GEE script で Oran NDVI 時系列(2018-2020)
- [ ] Quality filter(DetailedQA)
- [ ] 16-day → daily 線形補間

### Phase 2: 雲被覆 quality control
- [ ] MOD13Q1 の VI Quality flags 確認
- [ ] 不良データ除去
- [ ] Gap filling

### Phase 3: Phenology metrics
- [ ] Savitzky-Golay smoothing
- [ ] SOS/Peak/EOS extraction (logistic fit or threshold method)

---

## 8. v1 → v2 への改善ポイント

### v1 の問題(あれば、確認後)
- 多目的すぎて出力が散漫
- 解析A の最終結果と integrate されていない可能性

### v2 で focus すべき
1. **解析A の active period の validation** (最重要)
2. **NDVI peak window での τ 再 fit**(解析A 補強)
3. **Phenology metrics の inter-annual variation**(robust性)

省略してよいもの(v1 にあったかも):
- GPP プロキシの詳細(別研究で)
- 単純な NDVI vs EF 散布図(解析A で示せ済)

---

## 9. 解析C Success Criteria

### Minimum(必達)
- [ ] MOD13Q1 NDVI を Tarazona/Oran で時系列化
- [ ] NDVI で active period 客観定義
- [ ] 解析A の active month と一致 / 不一致を定量

### Stretch
- [ ] NDVI peak window 限定 τ で解析A τ を再現
- [ ] NDVI と EC LE の lag correlation
- [ ] Sentinel-2 NDVI(10m)で空間異質性確認

---

## 10. 次セッションで最初にやること

1. 本ファイル `ANALYSIS_C_PLAN.md` を読む
2. `ANALYSIS_A_FINAL.md` で active period 定義を確認
3. `analysis_C_v1.py` を確認(既存実装の機能棚卸し)
4. MOD13Q1 NDVI データの所在を確認
5. v2(整理版)を書くか、v1 を fix するかをユーザに確認

### 質問必要事項
- MOD13Q1 はもう抽出済か? あればパス
- v1 を 100% 書き直すか、増分 patch するか
- 解析B と並行か、順次か

---

## 11. 解析A との連動するパラメータ

解析A から継承する必須項目:
- 座標: Tarazona (39.266, -1.9397), Oran (38.82, -1.86)
- 期間: 2018-2024
- Active months 仮定:
  - Oran cereal: Nov-Jun
  - Tarazona almond: Jun-Sep(peak), Jan-Oct(broad)
- τ ≈ 3-4 d がベンチマーク値
- Amplitude ~4× scaling

→ これらを NDVI で confirm/refute するのが解析C の核心
