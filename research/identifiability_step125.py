"""旗125：旗44 の `e`（Tc²）と `d`（Tc·θz）は、θ–T が強く共変動しても識別できるか。

## なぜ要るのか（旗124 §4 が登録した自己点検）

**Davidson, Belk & Boone 1998 (GCB 4, 217–227) は、単一サイトでは温度と水分の効果を分けるのは
「統計的に非常に困難で、おそらく不可能」と結論し、Bunnell 型モデルで解が一意でないことを実演した。**
**彼らの非同定性は「主効果どうし」の話で、旗44 が判定に使うのは交互作用 `d` である。**
**だが θz ≈ ρ·Tc/SD_T + √(1−ρ²)·v と書けば `Tc·θz = ρ·Tc²/SD_T + √(1−ρ²)·Tc·v` であり、
第一項は `Tc²` そのもの＝ρ が強いほど `e` と `d` は共線になる。**
**旗44 は「d が拾えること」を θ–T 相関 1 水準で確かめただけで、識別可能性は一度も測っていない。**

## 何をするか

**旗44 の `analyze`・`verdict`・`_boot_d`・`_design` をそのまま import し、しきい値も `nboot` も
`block` も一つも変えずに、合成データだけを差し替えて被覆率と偽陽性率を測る。**

  ・θ–T 相関 ρ を **−0.5 / −0.7 / −0.85 / −0.95** と振る（旗124 の指定）。
  ・**`Var(θ)` を全水準で一定に保つ**——**保たないと `d` の真値（`0.6·SD_θ`）が水準ごとに動き、
    被覆率を水準間で比べられなくなる。**
  ・**帰無腕（真の d=0）**で **被覆率** と **★率（偽陽性）** を測る＝主判定。
  ・**陽性対照腕（真の d=+0.047）**で **検出率** を測る＝**門①**。
    **検出率 80% 未満の水準は、帰無腕の結果を「偽陽性が無い」の証拠に使わない**
    （検出力の無いところで偽陽性が出ないのは当たり前＝旗95 の欠陥 #31 と同型）。
  ・**落ちた反復は理由ごとに数える**（旗108・欠陥 #40）。**一つの籠に入れない。**

## 旗44 の合成との違い（**記録**）

**旗44 の `_synth` は `doy` を一様乱数で引いており、行が時間順に並んでいない
＝`_boot_d` のブロックブートが実質 iid ブートとして走っていた。**
**本ファイルは暦順（30 分刻み）にする。** **旗44 の結果を否定するものではなく、
旗44 の合成が CI のブロック構造を試せていなかったという事実の記録である。**

    .venv/bin/python research/identifiability_step125.py                 # 事前登録どおり
    .venv/bin/python research/identifiability_step125.py --reps 20       # 短縮（動作確認用）

事前登録：`research/PREREGISTRATION_step125.md`（**判定規則は実行前に固定済み**）。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from q10_confound_step44 import T0_LT, _design, analyze  # noqa: E402  ★旗44 の道具をそのまま使う
from runlog import tee_stdout  # noqa: E402

RHOS = (-0.5, -0.7, -0.85, -0.95)
TRUE_D_SCALE = 0.6          # 旗44 の `_synth('true')` と同じ強さ（model 単位では 0.6·SD_θ）
SD_T_NOISE = 2.0
SD_LNR_NOISE = 0.15
STAR_RATE_MAX = 0.05        # 事前登録の線（結果を見て動かさない）
STAR_RATE_WEAK = 0.08
COVERAGE_MIN = 0.90
POWER_MIN = 0.80            # 門①


def _season(n):
    """**暦順**の 30 分刻み（旗44 の `_synth` は一様乱数だった）。"""
    doy = np.arange(n) / 48.0
    return np.sin(2 * np.pi * (doy - 100) / 365.0)


def synth(rho, n, seed, true_d_scale):
    """θ–T 相関を `rho` に合わせ、**`Var(θ)` は水準に依らず一定**に保つ。

    `T = 12 + 12·s + ε_T`, `θ = 0.30 − A·s + ε_θ` のとき
      `Cov(θ,T) = −A·12·Var(s)`, `Var(T) = 144·Var(s) + σ_T²`, `Var(θ) = Var(s)·A² + σ_θ²`
    より `A = −ρ·√(Var(T)·V0) / (12·Var(s))`、`σ_θ² = V0 − Var(s)·A²`。
    `V0` は旗44 既定（A=0.05, σ_θ=0.07）の分散＝**全水準で固定**。
    """
    rng = np.random.default_rng(seed)
    s = _season(n)
    vs = float(np.var(s))
    v0 = vs * 0.05 ** 2 + 0.07 ** 2                    # 旗44 既定の Var(θ)
    var_T = 144.0 * vs + SD_T_NOISE ** 2
    A = -rho * np.sqrt(var_T * v0) / (12.0 * vs)
    var_eps = v0 - vs * A ** 2
    if var_eps <= 0:                                    # 目盛りが物理的に届かない水準
        return None
    T = 12.0 + 12.0 * s + rng.normal(0, SD_T_NOISE, n)
    th = 0.30 - A * s + rng.normal(0, np.sqrt(var_eps), n)
    lnR = 320.0 * (1.0 / (10 - T0_LT) - 1.0 / (T - T0_LT))   # Lloyd-Taylor（旗44 と同一）
    lnR += 0.8 * (th - 0.30) / 0.10                           # 水分は"量"を変える
    if true_d_scale:
        lnR += true_d_scale * ((th - 0.30) / 0.10) * (T - 12.0) / 10.0   # 感度も変える
    R = np.exp(lnR + rng.normal(0, SD_LNR_NOISE, n))
    d_true = true_d_scale * float(np.std(th))            # model 単位の真値
    return T, th, R, d_true, A


def _coefs_and_vif(T, th):
    """`e`・`d` の点推定に使う設計行列と、`Tc·θz` の VIF（他 3 列に対する）。"""
    Tc = T - T.mean()
    thz = (th - th.mean()) / th.std()
    A = _design(Tc, thz, quad=True)          # [Tc, Tc², θz, Tc·θz, 1]
    y = A[:, 3]
    X = np.column_stack([A[:, 0], A[:, 1], A[:, 2], A[:, 4]])
    resid = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    r2 = 1.0 - float(np.var(resid) / np.var(y))
    return A, (1.0 / (1.0 - r2) if r2 < 1 else np.inf)


def _fit_ed(A, R):
    coef = np.linalg.lstsq(A, np.log(R), rcond=None)[0]
    return float(coef[1]), float(coef[3])    # e, d


def run_arm(rho, true_d_scale, reps, n, seed0):
    """1 水準・1 腕。**落ちた反復は理由ごとに数える**（欠陥 #40）。"""
    out = {"n_ok": 0, "cover": 0, "star": 0, "neg": 0, "collapse": 0, "short": 0,
           "ci_none": 0, "unreachable": 0, "ds": [], "es": [], "widths": [], "corrs": [],
           "vifs": [], "d_true": None, "A": None}
    for i in range(reps):
        s = synth(rho, n, seed0 + i, true_d_scale)
        if s is None:            # **この n の季節グリッドでは、Var(θ) 固定のまま ρ に届かない**
            out["unreachable"] += 1
            continue
        T, th, R, d_true, A_amp = s
        out["d_true"], out["A"] = d_true, A_amp
        out["corrs"].append(float(np.corrcoef(th, T)[0, 1]))
        des, vif = _coefs_and_vif(T, th)
        out["vifs"].append(vif)
        e_hat, d_hat = _fit_ed(des, R)
        res = analyze(T, th, R)                        # ★旗44 の手続きそのもの
        if "note" in res:
            if "崩壊" in res["note"]:
                out["collapse"] += 1
            else:
                out["short"] += 1
            continue
        if res["ci"] is None:
            out["ci_none"] += 1
            continue
        lo, hi = res["ci"]
        out["n_ok"] += 1
        out["cover"] += int(lo <= d_true <= hi)
        out["star"] += int(lo > 0)
        out["neg"] += int(hi < 0)
        out["ds"].append(d_hat)
        out["es"].append(e_hat)
        out["widths"].append(hi - lo)
    return out


def _fmt_rate(k, n):
    return f"{k}/{n}={k / n:.0%}" if n else "—(0本)"


def main():
    p = argparse.ArgumentParser(description="旗125：e と d の識別可能性（合成のみ・実データ不要）")
    p.add_argument("--reps", type=int, default=100)
    p.add_argument("--n", type=int, default=20000)
    p.add_argument("--seed", type=int, default=1000)
    a = p.parse_args()
    tee_stdout("step125")

    print("=== 旗125：θ–T 共線のもとで `e`(Tc²) と `d`(Tc·θz) は識別できるか ===")
    print(f"  反復 {a.reps} / n={a.n}（30 分刻み・暦順＝{a.n / 48:.0f} 日）/ 旗44 の analyze をそのまま使用")
    print("  事前登録：research/PREREGISTRATION_step125.md（判定規則は実行前に固定済み）\n")

    null_rows, ctrl_rows = {}, {}
    for rho in RHOS:
        for arm, scale, store in (("帰無(d=0)", 0.0, null_rows),
                                  ("陽性対照(d=+0.047)", TRUE_D_SCALE, ctrl_rows)):
            r = run_arm(rho, scale, a.reps, a.n, a.seed)
            store[rho] = r
            if r["unreachable"] == a.reps:      # **黙って飛ばさない**（旗85 の作法）
                print(f"  ρ目標={rho:+.2f}  {arm:<20} **この n（{a.n / 48:.0f} 日）の季節グリッドでは"
                      f"Var(θ) を固定したまま届かない水準＝{a.reps} 反復すべて生成不能**\n")
                continue
            cm = np.mean(r["corrs"]) if r["corrs"] else float("nan")
            cs = np.std(r["corrs"]) if r["corrs"] else float("nan")
            print(f"  ρ目標={rho:+.2f}  {arm:<20} 実現corr(θ,T)={cm:+.3f}±{cs:.3f}  "
                  f"A={r['A']:.4f}  VIF={np.median(r['vifs']):.1f}  真のd={r['d_true']:+.4f}")
            print(f"     CI が出た反復 {r['n_ok']}/{a.reps}（落ち：当てはめ崩壊 {r['collapse']} / "
                  f"CI不定 {r['ci_none']} / 点不足 {r['short']} / 生成不能 {r['unreachable']}）")
            if r["n_ok"]:
                print(f"     被覆率 {_fmt_rate(r['cover'], r['n_ok'])}   "
                      f"★率 {_fmt_rate(r['star'], r['n_ok'])}   ×率 {_fmt_rate(r['neg'], r['n_ok'])}")
                print(f"     d̂ 中央値={np.median(r['ds']):+.5f}  SD(d̂)={np.std(r['ds']):.5f}  "
                      f"CI幅 中央値={np.median(r['widths']):.5f}  corr(ê,d̂)={np.corrcoef(r['es'], r['ds'])[0, 1]:+.3f}")
            print()

    # ---- 事前登録の判定規則をそのまま当てる（結果を見て線を動かさない） ----
    print("  === 判定（PREREGISTRATION_step125.md の規則をそのまま当てる）===")
    base_sd = np.std(null_rows[RHOS[0]]["ds"]) if null_rows[RHOS[0]]["ds"] else float("nan")
    print(f"  {'ρ':>6} {'被覆率':>8} {'★率':>8} {'対照検出率':>10} {'SD(d̂)倍率':>10}  門①  水準の判定")
    bad = []
    for rho in RHOS:
        nr, cr = null_rows[rho], ctrl_rows[rho]
        if not nr["n_ok"] or not cr["n_ok"]:
            why = "この n では生成不能" if nr["unreachable"] else "CI が出た反復が無い"
            print(f"  {rho:>+6.2f} {'—':>8} {'—':>8} {'—':>10} {'—':>10}  —   判定不能（{why}）")
            bad.append((rho, f"判定不能：{why}"))
            continue
        cov = nr["cover"] / nr["n_ok"]
        star = nr["star"] / nr["n_ok"]
        power = cr["star"] / cr["n_ok"]
        infl = np.std(nr["ds"]) / base_sd if base_sd else float("nan")
        gate = "通過" if power >= POWER_MIN else "**落ち**"
        if power < POWER_MIN:
            v = "門①落ち＝この水準では帰無腕を証拠に使わない"
        elif star > STAR_RATE_MAX or cov < COVERAGE_MIN:
            v = "**破れ**（★率>5% または 被覆率<90%）"
            bad.append((rho, f"★率{star:.0%}/被覆{cov:.0%}"))
        else:
            v = "○名目内"
        print(f"  {rho:>+6.2f} {cov:>8.0%} {star:>8.0%} {power:>10.0%} {infl:>10.1f}  {gate}  {v}")

    print()
    if not bad:
        print("  【結論】○識別できている——Davidson 1998 の非同定性は旗44 の `d` には当たらない。")
        print("          共線は偽陽性ではなく CI 幅（＝検出力の喪失）として現れた。")
    elif len(bad) == 1 and all(("破れ" in b[1] or "判定不能" in b[1]) for b in bad):
        print(f"  【結論】▲弱——破れは {bad[0][0]:+.2f} の 1 水準のみ。断定しない。")
        print("          この水準以上の θ–T 相関を持つ実サイトの ★ は要再点検。")
    else:
        print(f"  【結論】▲識別が崩れる——破れた水準：{[f'{b[0]:+.2f}' for b in bad]}")
        print("          旗44 の 18/36 は共線の強いサイトで水増しされている**可能性がある**。")
    print("\n  ★実データ 36 サイトの corr_thT の分布は**未確認**（この周は /mnt/hdd が無い）。")
    print("    ＝どの水準が実際に起きているかは言えない。分布の測定は実データのある周に回す。")
    print("  ★旗44 の判定規則はこの旗では変えない（旗102 の作法）。")


if __name__ == "__main__":
    main()
