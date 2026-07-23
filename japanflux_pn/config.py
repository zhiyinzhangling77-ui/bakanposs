"""解析全体の設定と、R&K Table 2 の 11 変数の正準順序。

すべての隣接行列・図はこの :data:`RK_VARS` の並びを使う。数値パラメータは
:class:`AnalysisConfig` に集約し、サイト間で共有・上書きできるようにする。
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math

# --- R&K Table 2 の正準順序（列/行の並びを全成果物で固定） --------------------
# (R&K 表記, 説明) のタプル。JapanFlux2024 実カラムへの対応は sites.var_map。
RK_VARS: list[str] = [
    "Rg",   # SW_IN_F      下向き短波放射 [W/m^2]
    "Ta",   # TA_F         気温 [degC]
    "VPD",  # VPD_F        飽差 [hPa]
    "Ts",   # TS_F_MDS_1   表層地温 [degC]
    "P",    # P_F          降水 [mm]
    "th",   # SWC_F_MDS_1  表層土壌水分 [%]   (θ)
    "gH",   # H_F_MDS      顕熱 [W/m^2]        (γH)
    "gLE",  # LE_F_MDS     潜熱 [W/m^2]        (γLE)
    "GER",  # RECO_DT_VUT_REF  生態系呼吸 [umol/m^2/s]
    "NEE",  # NEE_VUT_REF      正味生態系交換 [umol/m^2/s]
    "GEP",  # GPP_DT_VUT_REF   総一次生産 [umol/m^2/s]
]

# 人間可読ラベル（図・診断用）
RK_LABELS: dict[str, str] = {
    "Rg": "Rg", "Ta": "Ta", "VPD": "VPD", "Ts": "Ts", "P": "P",
    "th": "θ", "gH": "γH", "gLE": "γLE", "GER": "GER", "NEE": "NEE", "GEP": "GEP",
}


@dataclass(frozen=True)
class AnalysisConfig:
    """1 回の解析を規定するパラメータ一式。

    既定値は Ruddell & Kumar (2009) の設定に合わせてある。
    """

    # --- 前処理 -------------------------------------------------------------
    steps_per_day: int = 48          # 30 分値 (Δt = 30 min)
    anomaly_window_days: int = 5     # 5 日周期アノマリ (R&K §3)
    na_sentinel: float = -9999.0     # FLUXNET2015 欠測センチネル

    # --- 確率推定 -----------------------------------------------------------
    n_bins: int = 11                 # 固定区間ビン数 m (R&K)

    # --- transfer entropy のラグ -------------------------------------------
    lag_min: int = 1                 # 0.5 h  (= 1 step)
    lag_max: int = 36                # 18 h   (= 36 steps)

    # --- サロゲート有意性検定 ----------------------------------------------
    n_surrogates: int = 100          # M
    sig_c: float = 2.36              # Δ = μ_ss + c·σ_ss  (α = 0.01, 片側)

    # 有意ピーク τ' が属すべき「連続有意ラグ帯」の最小長。
    # 1 = 忠実な R&K (per-lag α=0.01, ラグ走査で族全体の偽陽性は膨らむ)。
    # 2 以上 = 時間的に連続した有意帯のみ採用し、単発クロスの偽陽性を抑制。
    # 独立 iid では 36 ラグ走査の偽陽性率が 1→36%, 2→~0%。実データ推奨は 2。
    peak_min_run: int = 1

    # --- ラグ三つ組の構築 ---------------------------------------------------
    # True: 実データでタイムスタンプがちょうど τ/1 step 離れる三つ組のみ採用
    # (ギャップ跨ぎを除外)。サロゲートは置換で時間意味が消えるため常に無視。
    gap_guard: bool = True

    seed: int = 0                    # 再現性のための RNG シード

    @property
    def lags(self) -> list[int]:
        """走査する τ (step 単位) のリスト。"""
        return list(range(self.lag_min, self.lag_max + 1))

    @property
    def log_m(self) -> float:
        """正規化定数 log(m)。T' = T/log(m), I' = I/log(m)。"""
        return math.log(self.n_bins)

    def lag_hours(self, tau_steps: int) -> float:
        """ラグ (step) を時間 [h] に変換。"""
        return tau_steps * 24.0 / self.steps_per_day
