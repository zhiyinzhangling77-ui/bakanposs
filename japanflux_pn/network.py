"""隣接行列の構築と R&K Table 1 の結合タイプ分類。

前処理済み :class:`~japanflux_pn.preprocess.PreprocessResult` から、11 変数間の
- AI  : 相対相互情報量 I' (%)              [対称]
- ATz : 有意な結合の Tz(X→Y, τ')            [有向, 非有意は NaN]
- Γ   : 特徴ラグ τ' [h]                     [有向, 非有意は NaN]
- ctype: 結合タイプ 1/2/3/4                  [有向]
を作る。各有向ペアで TE のラグ曲線とサロゲートしきい値を計算し、「最初の有意な
局所ピークラグ」τ' を検出して Tz = T(τ')/I とタイプを決める。
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .config import AnalysisConfig, RK_VARS
from . import information_theory as it
from .preprocess import PreprocessResult


# ---------------------------------------------------------------------------
# 特徴ラグ検出と分類
# ---------------------------------------------------------------------------
def _run_length_at(sig: np.ndarray) -> np.ndarray:
    """各インデックスが属する連続 True ラン (帯) の長さを返す。"""
    n = len(sig)
    out = np.zeros(n, dtype=int)
    i = 0
    while i < n:
        if sig[i]:
            j = i
            while j < n and sig[j]:
                j += 1
            out[i:j] = j - i
            i = j
        else:
            i += 1
    return out


def first_significant_peak(
    curve: np.ndarray, threshold: np.ndarray, lags: list[int], min_run: int = 1
) -> tuple[int | None, int | None]:
    """最初の有意な局所ピークラグ τ' を返す (R&K "first significant local peak")。

    有意 (curve > threshold) かつ局所最大 (両隣以上) の最初のラグ。厳密な局所
    ピークが無ければ (資格ある帯の) 最初の有意ラグにフォールバック。
    ``min_run`` を 2 以上にすると、長さ min_run 以上の連続有意帯に属するラグのみ
    採用し、単発クロスによる偽陽性を抑える。有意ラグが無ければ (None, None)。
    NaN は非有意扱い。
    """
    c = np.where(np.isfinite(curve), curve, -np.inf)
    sig = c > threshold
    if min_run > 1:
        sig = sig & (_run_length_at(sig) >= min_run)
    n = len(c)
    for i in range(n):
        if not sig[i]:
            continue
        left = c[i] >= c[i - 1] if i > 0 else True
        right = c[i] >= c[i + 1] if i < n - 1 else True
        if left and right:
            return lags[i], i
    if sig.any():
        i = int(np.argmax(sig))
        return lags[i], i
    return None, None


def classify_coupling(t_significant: bool, i_significant: bool, tz: float) -> int:
    """R&K Table 1 に従って結合タイプ 1-4 を返す。

    - Type 1 (同期支配)       : I 有意, T 非有意
    - Type 2 (フィードバック支配): T 有意, Tz < 1
    - Type 3 (強制支配)       : T 有意, Tz > 1
    - Type 4 (非結合)         : 両方非有意
    """
    if t_significant:
        return 3 if tz > 1.0 else 2
    if i_significant:
        return 1
    return 4


# ---------------------------------------------------------------------------
# 結果コンテナ
# ---------------------------------------------------------------------------
@dataclass
class NetworkResult:
    AI: pd.DataFrame           # 相対相互情報量 I' (%)、対称
    ATz: pd.DataFrame          # 有意な Tz(X->Y, τ')、非有意 NaN
    Gamma: pd.DataFrame        # 特徴ラグ τ' [h]、非有意 NaN
    ctype: pd.DataFrame        # 結合タイプ 1/2/3/4 (int)
    mi: pd.DataFrame           # 生 I [nats]、対称
    mi_threshold: pd.DataFrame # Δ(I) [nats]、対称
    te_curves: dict            # (src, dst) -> T(τ) [nats] 配列
    te_threshold: dict         # (src, dst) -> Δ(T, τ) [nats] 配列
    lags: list
    config: AnalysisConfig
    meta: dict = field(default_factory=dict)

    def save(self, outdir) -> None:
        """隣接行列群を CSV で保存。"""
        from pathlib import Path

        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        self.AI.to_csv(outdir / "AI_relative_mutual_information_pct.csv")
        self.ATz.to_csv(outdir / "ATz_transfer_ratio.csv")
        self.Gamma.to_csv(outdir / "Gamma_characteristic_lag_h.csv")
        self.ctype.to_csv(outdir / "coupling_type.csv")
        pd.Series(self.meta).to_csv(outdir / "meta.csv", header=False)


# ---------------------------------------------------------------------------
# メイン: ネットワーク構築
# ---------------------------------------------------------------------------
def build_network(pre: PreprocessResult) -> NetworkResult:
    """前処理結果から 4 種の隣接行列と結合タイプを構築する。"""
    cfg = pre.config
    m = cfg.n_bins
    lags = cfg.lags
    frame = pre.valid_frame
    step_index = pre.step_index
    rng = np.random.default_rng(cfg.seed)

    # 各変数を自身のレンジで bin index 化 (ラグ版は同一 edges を共有)
    bidx = {v: it.digitize_series(frame[v].to_numpy(), m) for v in RK_VARS}

    empty = lambda: pd.DataFrame(
        np.nan, index=RK_VARS, columns=RK_VARS, dtype=float
    )
    AI, ATz, Gamma = empty(), empty(), empty()
    mi_df, mi_thr_df = empty(), empty()
    ctype = pd.DataFrame(4, index=RK_VARS, columns=RK_VARS, dtype=int)
    for v in RK_VARS:
        ctype.loc[v, v] = 0
    te_curves: dict = {}
    te_threshold: dict = {}

    # --- 対称: 相互情報量とその有意性 (無向ペア) ---------------------------
    mi_sig: dict = {}
    for a in range(len(RK_VARS)):
        for b in range(a + 1, len(RK_VARS)):
            x, y = RK_VARS[a], RK_VARS[b]
            I = it.mutual_information_indices(bidx[x], bidx[y], m)
            s = it.surrogate_mi_stats(
                bidx[x], bidx[y], m, cfg.n_surrogates, cfg.sig_c, rng
            )
            for (p, q) in ((x, y), (y, x)):
                mi_df.loc[p, q] = I
                mi_thr_df.loc[p, q] = s["threshold"]
                AI.loc[p, q] = 100.0 * I / cfg.log_m       # I' (%)
                mi_sig[(p, q)] = I > s["threshold"]

    # --- 有向: transfer entropy, τ', Tz, タイプ ---------------------------
    for src in RK_VARS:
        for dst in RK_VARS:
            if src == dst:
                continue
            xi, yi = bidx[src], bidx[dst]
            curve = it.te_lag_curve(
                xi, yi, lags, m, step_index=step_index, gap_guard=cfg.gap_guard
            )
            stats = it.surrogate_te_stats(
                xi, yi, lags, m, cfg.n_surrogates, cfg.sig_c, rng
            )
            thr = stats["threshold"]
            te_curves[(src, dst)] = curve
            te_threshold[(src, dst)] = thr

            tau_prime, j = first_significant_peak(
                curve, thr, lags, min_run=cfg.peak_min_run
            )
            I = mi_df.loc[src, dst]
            i_sig = mi_sig[(src, dst)]
            if tau_prime is not None:
                tz = curve[j] / I if I > 0 else np.inf
                ATz.loc[src, dst] = tz
                Gamma.loc[src, dst] = cfg.lag_hours(tau_prime)
                ctype.loc[src, dst] = classify_coupling(True, i_sig, tz)
            else:
                ctype.loc[src, dst] = classify_coupling(False, i_sig, np.nan)

    meta = dict(pre.meta())
    meta.update({"n_lags": len(lags)})
    return NetworkResult(
        AI=AI, ATz=ATz, Gamma=Gamma, ctype=ctype,
        mi=mi_df, mi_threshold=mi_thr_df,
        te_curves=te_curves, te_threshold=te_threshold,
        lags=lags, config=cfg, meta=meta,
    )
