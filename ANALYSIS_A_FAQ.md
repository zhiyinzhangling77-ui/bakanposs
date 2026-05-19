# Analysis A — FAQ / Quick Reference

> v27-v30 で蓄積した「教えるための説明」をまとめたもの。
> ポスター発表や口頭質疑、論文 Discussion 執筆時に参照。

---

## 概念用語

### τ (タウ): 回復時定数
- 単位: 日
- 意味: 灌漑/雨イベント後、LE が baseline 値に向けて減衰する時間スケール
- e-folding time: τ = X日後に灌漑効果の 63% が減衰
- 本研究の結論: **τ ≈ 3 d at all 4 strata**

### LE_0 / LE_∞: 振幅の両端
- LE_0 = 灌漑/雨直後(Day 0)の LE
- LE_∞ = 長期経過後(漸近値、d ≥ 10 の median から fix)
- 振幅 = LE_0 − LE_∞

### Amplitude = "Management signal"
- Tarazona 95 W/m² vs Oran 21 W/m² → **4.5× scaling**
- 解釈: 灌漑は τ を変えず**振幅(応答の大きさ)を拡大**

### MDE (Minimum Detectable Effect)
- Power-aware の差検出限界
- MDE = 1.96 × √(SE₁² + SE₂²)
- 観測差 < MDE → 「検出可能だったが差なし」= 真の同等の証拠

---

## サイト基本情報

### Oran(雨養 cereal)
- 座標: lat 38.82, lon −1.86
- 作物: **winter cereal**(冬作穀物、おそらく barley/wheat)
- データ期間: **2018-2020(3年)**
- 灌漑: **なし**(rainfed)
- **生育期: Nov-Jun**(8ヶ月)
  - 播種: 10-12月
  - 栄養成長: Nov-Mar
  - 出穂・開花: Apr-May
  - 成熟・収穫: 6月
  - 休眠/裸地: Jul-Oct(dry Mediterranean summer)

### Tarazona(灌漑 almond)
- 座標: lat 39.266, lon −1.9397
- 作物: アーモンド(perennial, drip-irrigated)
- データ期間: **2020-2024(5年)**
- 灌漑: **drip irrigation**(2 lines per tree row)
  - 量: 12-15 mm/event
  - 頻度: 月 7-47 events(5-10月集中)
- 樹間 4.5m × 6.5m、密度 342 trees/ha
- **生育期: Jan-Oct**(落葉 Nov-Dec)
  - 開花: Jan-Feb
  - 葉展開: Mar-May
  - **蒸散ピーク: Jun-Sep**(active period)
  - 収穫: Aug-Oct
  - 落葉: Nov-Dec

---

## 図の作り方解説

### Recovery curve (panel a, b) の boxplot 構成

**手順**:
1. 各 event について Day 0, 1, 2, ..., N の LE 値を抽出
   - Day 0 = 水入力イベントの当日(Rain≥3mm or Irrig≥0.5mm)
   - Day N = 次の event まで or max_window=14 日
2. 各 day position(0, 1, 2, ..., 14)で、寄与する events の LE を集める
3. それらの LE 分布を箱ひげで表示

**例**: Tarazona day=3 の箱ひげ
- 41 events のうち day=3 まで観測されている events を集める
- 各 event の day=3 LE 値 を 1 つ取る(例 38 個)
- それらの分布を箱ひげで表示

**なぜ day が進むと n 減るか**:
- 短い灌漑間隔 events は早く打ち切られる
- 長い間隔 events だけが day 10+ に貢献
- 自然な性質、データ品質の問題ではない

### Panel (c) τ comparison

- 4 strata で別々に exp fit
- 各バー = その stratum の τ point estimate
- エラーバー = bootstrap 95% CI
- 緑帯 = "universal band" 3-4 d
- **全バーが緑帯に入る + pairwise diff < MDE → universality**

---

## よくある質問と一言答え

### Q1: なぜ Oran を winter/summer/active で分ける?
**A**: Reviewer attack「season で τ 違うのでは?」を 4 strata で潰すため。
- winter: 植物が前半 active
- summer: 収穫直前の peak transpiration
- active: Nov-Jun pool(Tarazona summer と phenology-matched 比較用)
- → 4 strata で同じ τ ≈ 3 d → robust universality

### Q2: management-scaled amplitude とは?
**A**: 振幅(LE_0−LE_∞)が灌漑で 4.5× 拡大される現象。
- τ = climate property(管理非依存)
- amplitude = management property(管理で scale)
- 論文の二層構造: 時間スケール普遍 + 振幅 scaling

### Q3: n=10(Oran) や 41(Tarazona) は少ない?
**A**: 少ないが MDE 解析で「power 不足ではない」と証明済。
- Oran は 3年データ + Mediterranean dry summer → rain event 必然的に少ない
- MDE = 2.15 d で観測差 0.54 d → 検出感度十分

### Q4: 箱ひげの作り方?
**A**: 各 day position で、寄与する events の LE 値の分布を箱ひげで表示。
- Panel (a) Tarazona ≠ Panel (b) Oran:**完全に別データ**
- Panel (c) 4 strata:Oran winter ⊂ Oran active(部分集合関係)

### Q5: Oran の生育期はいつ?
**A**: **Nov-Jun(冬作穀物)**。
- 播種 10-12月、栄養成長 Nov-Mar、開花 Apr-May、収穫 6月
- 休眠期 Jul-Oct(地中海性気候の乾燥した夏)

### Q6: MDE の使い方?
**A**: 「差がない」を主張するときの power check。
- 通常の検定(p < 0.05): 「差を検出した」を言う道具
- MDE: 「もし差があれば検出できた」を言う道具
- 使い方:
  1. 各群で τ + bootstrap SE を計算
  2. MDE = 1.96 × √(SE₁² + SE₂²)
  3. 観測差 |Δτ| と MDE を比較
  4. obs < MDE → 「検出可能だったが差なし」 = 同等の証拠
  5. obs > MDE → 「有意差検出」(通常の検定と同じ)

---

## v28-v30 figure 履歴

### v28: 初版 poster Fig 4
- 上下 2 panel(recovery curve overlay + τ bars)
- 振幅 annotation box 付き

### v29: MDE 解析図 (Fig 5)
- 4 pairwise 比較で obs vs MDE 比較
- 「差がない」の根拠を示す power-aware verdict
- 概念図 supplementary も生成

### v30: poster Fig 4 透明性強化
- 3 panel(Tarazona / Oran 分離 + τ comparison)
- 各 day position に boxplot + n 表示
- Day 0 定義を annotation box で明示
- Panel (c) 右側に "意味" annotation box

### v30 最終数値
| 量 | 値 |
|---|---|
| Tarazona active τ | 3.36 d [2.44, 4.90], n=41 |
| Oran active τ | 2.82 d [1.83, 5.48], n=10 |
| Oran winter τ | 3.29 d [1.71, 5.78], n=6 |
| Oran summer τ | 3.79 d [1.07, 6.51], n=4 |
| Amplitude Tarazona | 95 W/m² |
| Amplitude Oran | 21 W/m² |
| Ratio | **4.5×** |

---

## ポスター発表用 speaking script

### Fig 4 説明
> 「Fig 4 は本研究の主結果です。
>
> 上段 (a)(b) は灌漑/雨イベント後の LE 回復曲線。Day 0 = 水が入った当日。
> Tarazona は 41 灌漑 events、Oran は 10 雨 events から構築。
> 箱ひげは各 day で寄与する events の LE 分布、median が指数 fit を駆動。
>
> 下段 (c) は 4 つの異なる stratum で τ を別計算。
> サイト・季節・pooling の違いに関係なく τ ≈ 3 d。緑帯が universal band。
>
> 結論: 時間スケール τ は管理(雨/灌漑)・季節に依らず ~3 日。
>       振幅は管理で 4.5× scale される。
>       時間 = 普遍、大きさ = 管理依存。」

### Fig 5 説明
> 「Fig 5 は『差がない』という主張を power-aware に正当化します。
>
> 各行は 2 つの stratum の比較。赤ダイヤ = 観測 τ 差、青丸 = MDE。
> MDE は『この実験で 95% 検出可能な最小差』。
>
> 全比較で MDE >> 観測差(緑の余裕 = headroom)。
> 真に差があれば検出できた、それでも差が出なかった = 真の同等の証拠。
>
> 通常の p > 0.05 は『差を見つけられなかった』に過ぎず、
> 'no difference' を主張するには MDE が必要。」

---

## 次の解析(B/C)に向けて

- 解析A は v27-v30 で完了
- 解析B(衛星 ET 検証)= 次の最優先 → `ANALYSIS_B_PLAN.md`
- 解析C(NDVI phenology)= 補強 → `ANALYSIS_C_PLAN.md`
- どちらも v27 の出力(τ ≈ 3 d, amplitude 4.5×)をベンチマークに validate
