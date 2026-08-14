"""旗14：O-information をゼロから最小実装し、"相乗 vs 冗長"を答えの分かる合成系で検証。

ペアの相互情報や TE は「2 変数」の話。3 変数以上になると、情報が
  ・**冗長**（共通原因で重複：Z→X1,X2,X3。どれか1本見れば他も分かる）
  ・**相乗**（組にしか宿らない創発：X3=X1⊕X2。どの1本・どの2本の一部を見ても分からず、
    3本揃って初めて意味が出る）
のどちらに支配されるかが問える。これを1数値で表すのが O-information（Rosas 2019）：
  Ω(X) = (n−2)·H(X) + Σ_i [ H(X_i) − H(X_{−i}) ]   （= 全相関 TC − 双対全相関 DTC）
  Ω>0 冗長支配 ／ Ω<0 相乗支配。

答えが分かる2つの合成系で符号を確かめる（北極星：豊かな高次構造を貧しい観測でどう捉えるか）：
  1. 共通原因 Z→X1,X2,X3  → 冗長支配（Ω>0）を期待
  2. XOR   X3 = X1 xor X2 → 相乗支配（Ω<0）を期待
さらに、4変数ヒストは疎で Ω が負にバイアスするので、本体と同じく
  Miller-Madow 補正 ＋ シャッフルサロゲート z で判定する（絶対符号でなく z の符号）。

    python research/o_information_step14.py
"""

from __future__ import annotations

import numpy as np


# ---- 最小の情報量プリミティブ（numpy だけ）--------------------------------
def _digitize(x: np.ndarray, m: int) -> np.ndarray:
    """連続値を等頻度でなく等幅 m ビンへ（本体の digitize_series と同じ思想の最小版）。"""
    x = np.asarray(x, dtype=float)
    lo, hi = x.min(), x.max()
    if hi <= lo:
        return np.zeros(len(x), dtype=np.int64)
    idx = np.floor((x - lo) / (hi - lo) * m).astype(np.int64)
    return np.clip(idx, 0, m - 1)


def _entropy(cols: list[np.ndarray], m: int, correct: bool = False) -> float:
    """列（ビンインデックス）の同時エントロピー H [nats]。correct=Miller-Madow。"""
    code = np.zeros(len(cols[0]), dtype=np.int64)
    for c in cols:
        code = code * m + c
    n = len(code)
    counts = np.bincount(code)
    counts = counts[counts > 0]
    p = counts / n
    H = float(-np.sum(p * np.log(p)))
    if correct:                      # Miller-Madow: +(K-1)/(2N)
        K = len(counts)
        H += (K - 1) / (2.0 * n)
    return H


def total_correlation(cols: list[np.ndarray], m: int, correct: bool = False) -> float:
    """TC = Σ H(X_i) − H(X)。共有された依存の総量（冗長寄り）。"""
    return sum(_entropy([c], m, correct) for c in cols) - _entropy(cols, m, correct)


def o_information(cols: list[np.ndarray], m: int, correct: bool = False) -> float:
    """Ω = (n−2)H(X) + Σ_i[H(X_i) − H(X_{−i})]。Ω>0 冗長支配 / Ω<0 相乗支配。"""
    n = len(cols)
    if n < 3:
        raise ValueError("O-information は 3 変数以上")
    h_all = _entropy(cols, m, correct)
    omega = (n - 2) * h_all
    for i in range(n):
        rest = [cols[j] for j in range(n) if j != i]
        omega += _entropy([cols[i]], m, correct) - _entropy(rest, m, correct)
    return omega


def dual_total_correlation(cols: list[np.ndarray], m: int, correct: bool = False) -> float:
    """DTC = TC − Ω（相乗寄り）。参考表示用。"""
    return total_correlation(cols, m, correct) - o_information(cols, m, correct)


def surrogate_z(cols: list[np.ndarray], m: int, n_surr: int, correct: bool,
                seed: int = 0) -> tuple[float, float, float]:
    """各変数を独立シャッフル（全依存を壊す＝真の Ω=0）したヌルからの z。
    疎性の有限標本バイアスはヌルにも同じだけ乗るので、z の符号で冗長/相乗を正しく判定できる。
    """
    obs = o_information(cols, m, correct)
    rng = np.random.default_rng(seed)
    n = len(cols[0])
    samp = np.empty(n_surr)
    for s in range(n_surr):
        shuf = [c[rng.permutation(n)] for c in cols]
        samp[s] = o_information(shuf, m, correct)
    mu, sigma = float(samp.mean()), float(samp.std())
    z = (obs - mu) / sigma if sigma > 0 else np.nan
    return obs, mu, z


# ---- 答えの分かる合成系 ----------------------------------------------------
def make_redundant(n=6000, m=6, noise=0.3, seed=0):
    """共通原因 Z→X1,X2,X3。3本が同じ Z を共有＝冗長支配（Ω>0）を期待。"""
    rng = np.random.default_rng(seed)
    z = rng.normal(0, 1, n)
    cols = [_digitize(z + rng.normal(0, noise, n), m) for _ in range(3)]
    return cols


def make_synergy(n=6000, m=6, noise=0.05, seed=0):
    """XOR: X3 = sign(X1) xor sign(X2)。どの2本の符号も X3 と独立、3本で初めて決まる
    ＝相乗支配（Ω<0）を期待。連続値＋微ノイズを等幅ビン化しても構造は残る。"""
    rng = np.random.default_rng(seed)
    x1 = rng.normal(0, 1, n)
    x2 = rng.normal(0, 1, n)
    b = ((x1 > 0).astype(int) ^ (x2 > 0).astype(int))     # XOR ビット
    x3 = np.where(b == 1, 1.0, -1.0) + rng.normal(0, noise, n)
    return [_digitize(x1, m), _digitize(x2, m), _digitize(x3, m)]


def make_independent(n=6000, m=6, seed=0):
    """独立 3 本（帰無）。Ω≈0（z 非有意）を期待＝検定の偽陽性チェック。"""
    rng = np.random.default_rng(seed)
    return [_digitize(rng.normal(0, 1, n), m) for _ in range(3)]


def _judge(z: float, c: float = 2.36) -> str:
    if not np.isfinite(z) or abs(z) < c:
        return "有意でない"
    return "冗長支配 (共通駆動)" if z > 0 else "★相乗支配 (創発)"


def main() -> None:
    m = 6
    n_surr = 400
    print("=== 旗14: O-information をゼロから、答えの分かる合成系で検証 ===")
    print(f"  ビン m={m}, Miller-Madow 補正あり, シャッフル {n_surr} 回, |z|≥2.36 で有意\n")
    cases = [
        ("共通原因 Z→X1,X2,X3 (冗長を期待 Ω>0)", make_redundant(m=m), "冗長"),
        ("XOR X3=X1⊕X2       (相乗を期待 Ω<0)", make_synergy(m=m), "相乗"),
        ("独立3本            (帰無 z≈0)",        make_independent(m=m), "無"),
    ]
    print(f"  {'系':<34} {'Ω(MM)':>9} {'TC':>7} {'DTC':>7} {'z':>8}  判定")
    ok = True
    for name, cols, expect in cases:
        omega, mu, z = surrogate_z(cols, m, n_surr, correct=True)
        tc = total_correlation(cols, m, correct=True)
        dtc = dual_total_correlation(cols, m, correct=True)
        verdict = _judge(z)
        print(f"  {name:<34} {omega:9.3f} {tc:7.3f} {dtc:7.3f} {z:8.1f}  {verdict}")
        if expect == "冗長" and not (z > 2.36):
            ok = False
        if expect == "相乗" and not (z < -2.36):
            ok = False
        if expect == "無" and abs(z) >= 2.36:
            ok = False
    print("\n  → 期待どおり: " + ("✅ 共通原因=冗長(z>0)、XOR=相乗(z<0)、独立=非有意"
                                if ok else "⚠ 期待と不一致（ノイズ/ビン/サロゲート数を調整）"))
    print("  意味: 相互情報(2変数)では見えない『情報が組にしか宿る=相乗』を、Ωの符号が捉える。")
    print("        4変数ヒストは疎でΩが負にバイアスするので、絶対値でなく必ず z(サロゲート)で判定。")
    print("        実データは本体 japanflux_pn.oinfo_analysis（同じΩ・MM・zで6サブ系×多年）。")
    print("        → python -m japanflux_pn.oinfo_analysis --site JP-Tak --multiyear")


if __name__ == "__main__":
    main()
