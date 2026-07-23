# japanflux_pn — 情報理論プロセスネットワーク解析

Ruddell & Kumar (2009, WRR 45, W03419) の生態水文プロセスネットワーク解析を
JapanFlux2024 (Ueyama et al. 2025, ESSD) の FLUXNET2015 互換 CSV に適用する。
1 サイトに閉じない再利用可能な設計（JP-Tak / JP-Ta2 / JP-BBY / JP-Mse）。

## モジュール

| module | 役割 |
|---|---|
| `config.py` | `AnalysisConfig` + R&K Table 2 の 11 変数正準順序 (`RK_VARS`) |
| `sites.py` | サイトレジストリ（データ位置 + 変数名マッピング）。多サイト拡張の効き所 |
| `preprocess.py` | COREVARS HH 読込 + 5 日前方窓アノマリ + listwise deletion |
| `information_theory.py` | エントロピー / MI / TE(Knuth形) / shuffled surrogate（変数名非依存カーネル） |
| `network.py` | 隣接行列 AI/ATz/Γ 構築、τ' 検出、結合タイプ分類 |
| `viz.py` | ネットワーク図 / タイプ行列 / ラグ診断プロット |
| `run_site.py` | エントリポイント |

## 実行

```bash
python -m japanflux_pn.run_site --site JP-Tak --year 2003 --month 7 --peak-min-run 2
```

出力は `japanflux_pn/outputs/<site>_<year><month>/` に:
隣接行列 CSV 4 種（AI, ATz, Γ, coupling_type）、`meta.csv`、
`network.png` / `coupling_type.png` / `lag_diagnostics.png`。

## 方法論メモ（実装上の判断）

- **5 日アノマリ (R&K §3)**: 各時刻 t から同時刻帯の直後 5 日平均を引く（t を
  day0 に含む前方窓）。窓 5 点が全て揃う時のみ確定。対象月末に窓日数ぶんの
  バッファを読んで末尾損失を防ぐ。1 ヶ月で有効 ~1200 点（推奨 500–1500）。
- **確率推定**: 固定区間ビン m=11。多次元同時分布は混合基数エンコード +
  `np.bincount`。ラグ版の同一変数は同じ bin edges を共有。
- **ラグ三つ組**: 時間順配列上の位置ラグ。実データは `gap_guard` でギャップ跨ぎ
  を除外（`step_index`）、サロゲートは置換で時間意味が消えるため常に off。実測と
  サロゲートで同一 N。
- **TE のバイアス floor**: m=11 の 3 次元ヒストは有限標本で正のバイアス（~0.15
  nats）を持つ。有意性は 0 でなく `μ_ss + c·σ_ss`（c=2.36, α=0.01 片側）に対して判定。

### ⚠ ラグ走査の多重比較（重要）

per-lag のサロゲート検定は α=0.01 で正しく較正されているが、τ を 36 本走査して
「どれか 1 本でも有意なら結合」とすると族全体の偽陽性率が 1−0.99³⁶ ≈ 30% に
膨らむ。独立 iid 合成データでの実測:

| ルール | 偽陽性率 |
|---|---|
| 単一ラグ | 0% |
| 36 ラグ・どれか 1 本 (`peak_min_run=1`, 忠実 R&K) | ~36% |
| 36 ラグ・2 連続有意を要求 (`peak_min_run=2`) | ~0% |

真の情報流ピークは時間的に連続した有意帯として出るため、`peak_min_run=2`
（有意ピークが長さ 2 以上の連続有意帯に属することを要求）で真の結合を保ったまま
単発クロスの偽陽性を除去できる。**既定は 1（忠実な R&K 再現）、実データ解析では
`--peak-min-run 2` を推奨。**

### Tz の解釈上の注意

`Tz = T(τ')/I`。無結合ペアは I≈0（有限標本バイアスのみ）なので、万一 T が偶発的に
有意になると Tz が発散的に大きくなる（偽陽性の指標にもなる）。真の強制結合の Tz は
通常 O(1)–O(数)。`peak_min_run=2` で偽陽性を除けばこの人工的な巨大 Tz も消える。

## テスト

```bash
python -m pytest tests/ -q
```

合成データ（結合ロジスティック写像, R&K §2.2）で TE 実装を検証:
方向非対称性・ピークラグ=1・サロゲート有意性・シャッフルで結合消滅、
前処理のアノマリ/欠測処理、`build_network` の統合（既知結合の検出）。
