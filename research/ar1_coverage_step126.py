"""旗126：残差が自己相関しているとき、旗44 のブロックブート CI は名目の被覆を保つか。

## なぜ要るのか（旗125 が残した宿題の第一候補）

**旗125 は θ–T 共線を振って旗44 の `d` の被覆率を測ったが、雑音は `N(0, 0.15)` の iid だった。
＝ブロックブート（`block = 48*7`）が本来相手にしている自己相関を、一度も入れていない。**
**旗44 の ★ 18/36 が安全だという旗125 の結論は、iid 雑音のもとでの数である。**

## 何をするか

**旗44 の `analyze`／`_boot_d`／`_design` と、旗125 の `synth`／`_coefs_and_vif`／`_fit_ed` を
そのまま import し、`lnR` に載る観測雑音だけを AR(1) に差し替える。**
**しきい値も `nboot` も `block` も一つも変えない。**

  ・φ = 0 / 0.5 / 0.9 / 0.98 を **ρ=−0.5 に固定して**振る（自己相関を共線から切り離す）。
  ・**周辺 SD を φ に依らず 0.15 に固定する**（固定しないと雑音の大きさ自体が動く）。
  ・**帰無腕（真の d=0）**で被覆率と ★率＝主判定。
  ・**陽性対照腕（門①）の真の d を、そのセル自身の判定境界の 3 倍に較正する**
    ——**旗125 の欠陥 #46（対照が境界の 26 倍で常に 100% 検出）の是正。**
  ・**補助**：`block=1`（iid ブート）の被覆率を並べる＝ブロックブートが仕事をしているかの直接の対照。
  ・**落ちた反復は理由ごとに数える**（旗108・欠陥 #40）。

    .venv/bin/python research/ar1_coverage_step126.py                    # 事前登録どおり
    .venv/bin/python research/ar1_coverage_step126.py --reps 5 --pilot 3 # 短縮（動作確認用）

事前登録：`research/PREREGISTRATION_step126.md`（**判定規則は実行前に固定済み**）。
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
from scipy.signal import lfilter

sys.path.insert(0, str(Path(__file__).resolve().parent))
from identifiability_step125 import SD_LNR_NOISE, _coefs_and_vif, _fit_ed, _season  # noqa: E402
from q10_confound_step44 import T0_LT, _fit_d, analyze  # noqa: E402
from runlog import tee_stdout  # noqa: E402

PHIS = (0.0, 0.5, 0.9, 0.98)
RHO_MAIN = -0.5
PHI_TH_WHITE = 0.0          # 事前登録どおりの軸（θ の雑音は iid＝旗44/125 の生成器のまま）
PHI_TH_RED = 0.99           # ★追補の軸（θ も自己相関する＝実測に近い形。理由は下の追補を見よ）
CROSS = (-0.85, 0.9)        # 補助（共線と自己相関が同時に効くセル・主判定に使わない）
TRUE_D_SCALE_REF = 0.6      # 旗44/125 の対照の強さ（**参考として印字するだけ**）
SD_T_NOISE = 2.0
BLOCK = 48 * 7              # 旗44 と同じ（変えない）
CTRL_MULT = 3.0             # 対照＝判定境界の 3 倍（欠陥 #46 の是正・1〜3 倍の上端）
PILOT_MIN = 5               # 較正パイロットで CI が出た反復がこれ未満なら較正しない
STAR_RATE_MAX = 0.05        # 事前登録の線（結果を見て動かさない）
STAR_RATE_WEAK = 0.08
COVERAGE_MIN = 0.90
POWER_MIN = 0.80            # 門①


def _ar1(rng, n, sd, phi):
    """周辺 SD を `sd` に固定した定常 AR(1)。**φ=0 なら iid（旗125 と同一）。**"""
    if phi == 0.0:
        return rng.normal(0.0, sd, n)
    eta = rng.normal(0.0, sd * np.sqrt(1.0 - phi ** 2), n)
    zi = np.array([phi * rng.normal(0.0, sd)])      # 初期値は定常分布から＝過渡を入れない
    return lfilter([1.0], [1.0, -phi], eta, zi=zi)[0]


def synth_ar1(rho, phi, n, seed, true_d_scale, phi_th=0.0):
    """旗125 の `synth` と同一。**違いは `lnR` の観測雑音が AR(1) であることだけ**（`phi_th=0` のとき）。

    `Var(θ)` は ρ に依らず固定（旗125 の理由）／周辺の雑音 SD は φ に依らず固定（本旗の理由）。

    **`phi_th` は追補**（動作確認の後・判定の実行の前に足した。`PREREGISTRATION_step126.md` の追補を見よ）：
    **θ の雑音を AR(1) にして、回帰子 `Tc·θz` を赤くする。**
    **周辺分散は `phi_th` に依らず保たれるので、`Var(θ)` も `corr(θ,T)` も期待値としては動かない。**
    """
    rng = np.random.default_rng(seed)
    s = _season(n)                                   # 暦順の 30 分刻み（旗125 で直した点）
    vs = float(np.var(s))
    v0 = vs * 0.05 ** 2 + 0.07 ** 2                  # 旗44 既定の Var(θ)
    var_T = 144.0 * vs + SD_T_NOISE ** 2
    A = -rho * np.sqrt(var_T * v0) / (12.0 * vs)
    var_eps = v0 - vs * A ** 2
    if var_eps <= 0:
        return None
    T = 12.0 + 12.0 * s + rng.normal(0, SD_T_NOISE, n)
    th = 0.30 - A * s + _ar1(rng, n, np.sqrt(var_eps), phi_th)
    lnR = 320.0 * (1.0 / (10 - T0_LT) - 1.0 / (T - T0_LT))   # Lloyd-Taylor（旗44 と同一）
    lnR += 0.8 * (th - 0.30) / 0.10
    if true_d_scale:
        lnR += true_d_scale * ((th - 0.30) / 0.10) * (T - 12.0) / 10.0
    eps = _ar1(rng, n, SD_LNR_NOISE, phi)            # ★ここだけが旗125 との違い
    R = np.exp(lnR + eps)
    d_true = true_d_scale * float(np.std(th))
    return T, th, R, d_true, A, eps


def _boot_d_iid(Tc, thz, lnR, nboot=200, seed=0):
    """**補助（判定に使わない）**：ブロックを外した iid ブート。

    **旗44 の `_boot_d` と回数も推定量も同じで、再標本だけを行単位の復元抽出にしたもの。**
    `_boot_d(block=1)` を呼ぶと 1 行ずつの `arange` を 200×n 回作って現実的な時間で終わらないので、
    **同じ意味の再標本をここで直に書く**（**旗44 の道具は改変していない**）。
    """
    rng = np.random.default_rng(seed)
    n = len(Tc)
    ds = []
    for _ in range(nboot):
        idx = rng.integers(0, n, n)
        try:
            ds.append(_fit_d(Tc[idx], thz[idx], lnR[idx], True)[1])
        except np.linalg.LinAlgError:
            continue
    if len(ds) < 30:
        return None
    return (float(np.percentile(ds, 2.5)), float(np.percentile(ds, 97.5)))


def _acf1(x):
    x = np.asarray(x, float)
    x = x - x.mean()
    v = float(np.dot(x, x))
    return float(np.dot(x[1:], x[:-1]) / v) if v > 0 else float("nan")


def _efold_steps(phi):
    return float("nan") if phi <= 0 else -1.0 / np.log(phi)


def run_arm(rho, phi, true_d_scale, reps, n, seed0, with_iid=False, phi_th=0.0):
    """1 セル・1 腕。**落ちた反復は理由ごとに数える**（欠陥 #40）。"""
    out = {"n_ok": 0, "cover": 0, "star": 0, "neg": 0, "collapse": 0, "short": 0,
           "ci_none": 0, "unreachable": 0, "ds": [], "es": [], "widths": [], "corrs": [],
           "vifs": [], "acf_eps": [], "acf_res": [], "acf_x": [], "d_true": None, "A": None,
           "iid_ok": 0, "iid_cover": 0, "iid_star": 0, "iid_widths": []}
    for i in range(reps):
        s = synth_ar1(rho, phi, n, seed0 + i, true_d_scale, phi_th)
        if s is None:
            out["unreachable"] += 1
            continue
        T, th, R, d_true, A_amp, eps = s
        out["d_true"], out["A"] = d_true, A_amp
        out["corrs"].append(float(np.corrcoef(th, T)[0, 1]))
        out["acf_eps"].append(_acf1(eps))
        des, vif = _coefs_and_vif(T, th)
        out["acf_x"].append(_acf1(des[:, 3]))          # ★回帰子 `Tc·θz` 自身の ACF1（追補の要点）
        out["vifs"].append(vif)
        e_hat, d_hat = _fit_ed(des, R)
        out["acf_res"].append(_acf1(np.log(R) - des @ np.linalg.lstsq(des, np.log(R), rcond=None)[0]))
        res = analyze(T, th, R)                      # ★旗44 の手続きそのもの
        if "note" in res:
            out["collapse" if "崩壊" in res["note"] else "short"] += 1
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
        if with_iid:                                  # 補助：ブロックを外した対照（判定に使わない）
            Tc = T - T.mean()
            thz = (th - th.mean()) / th.std()
            ci1 = _boot_d_iid(Tc, thz, np.log(R))
            if ci1 is not None:
                out["iid_ok"] += 1
                out["iid_cover"] += int(ci1[0] <= d_true <= ci1[1])
                out["iid_star"] += int(ci1[0] > 0)
                out["iid_widths"].append(ci1[1] - ci1[0])
    return out


def calibrate(rho, phi, pilot, n, seed0, phi_th=0.0):
    """**そのセル自身の判定境界（CI 半幅の中央値）を測り、対照の真の d をその 3 倍に置く。**

    **較正に使うのは帰無腕の CI 幅だけで、主判定の量（被覆率・★率）は見ない。**
    戻り値 `(半幅, 対照の d_scale, 実現した対照の真値)`。**測れなければ `(nan, None, nan)`。**
    """
    r = run_arm(rho, phi, 0.0, pilot, n, seed0, phi_th=phi_th)
    if r["n_ok"] < PILOT_MIN:
        return float("nan"), None, float("nan"), r
    hw = float(np.median(r["widths"])) / 2.0
    s = synth_ar1(rho, phi, n, seed0, 0.0, phi_th)
    sd_th = float(np.std(s[1])) if s is not None else float("nan")
    scale = CTRL_MULT * hw / sd_th if sd_th > 0 else None
    return hw, scale, (scale * sd_th if scale else float("nan")), r


def _fmt_rate(k, n):
    return f"{k}/{n}={k / n:.0%}" if n else "—(0本)"


def main():
    p = argparse.ArgumentParser(description="旗126：AR(1) 残差でのブロックブート被覆（合成のみ）")
    p.add_argument("--reps", type=int, default=100)
    p.add_argument("--pilot", type=int, default=20)
    p.add_argument("--n", type=int, default=20000)
    p.add_argument("--seed", type=int, default=2000)
    p.add_argument("--no-iid", action="store_true", help="補助の iid ブート対照を省く（時間短縮）")
    a = p.parse_args()
    tee_stdout("step126")
    t0 = time.time()

    print("=== 旗126：残差が AR(1) のとき、旗44 のブロックブート CI は被覆を保つか ===")
    print(f"  反復 {a.reps}（較正パイロット {a.pilot}）/ n={a.n}（30 分刻み・暦順＝{a.n / 48:.0f} 日）")
    print(f"  block={BLOCK} 刻み（{BLOCK / 48:.0f} 日）・nboot=200＝旗44 のまま。analyze もそのまま import。")
    print("  事前登録：research/PREREGISTRATION_step126.md（判定規則は実行前に固定済み）\n")

    print("  --- 前提の事実確認 1/2：自己相関の長さとブロック長の関係 ---")
    for phi in PHIS:
        ef = _efold_steps(phi)
        if np.isnan(ef):
            print(f"    φ={phi:.2f}  e-fold=—（iid）                       ブロック長 {BLOCK} 刻み＝∞ 倍")
        else:
            print(f"    φ={phi:.2f}  e-fold={ef:6.1f} 刻み（{ef / 2:5.1f} 時間）  "
                  f"ブロック長は e-fold の {BLOCK / ef:5.1f} 倍")
    print("    （**ブロック長が相関長より十分長い**ことがブロックブートの効く前提。上で確かめている）\n")

    print("  --- 追補の理由（**動作確認の後・判定の実行の前に足した**）---")
    print("    **動作確認（反復 3）で CI 幅が φ にほとんど反応しなかった**（倍率 1.0〜1.1）。")
    print("    **原因は θ の雑音が iid で、回帰子 `Tc·θz` が白いこと**——白い回帰子は赤い誤差と")
    print("    スペクトルがほとんど重ならないので、自己相関があっても推定量の分散は膨らまない。")
    print(f"    **実測の土壌水分は 30 分刻みで強く自己相関する。** そこで **θ の雑音も AR(1)（φ_θ={PHI_TH_RED}）**")
    print("    にした軸を足す。**事前登録の軸（φ_θ=0）はそのまま残し、両方を報告する。**")
    print("    **判定規則は両軸で同一。** 詳細は PREREGISTRATION_step126.md の追補。\n")

    cells = ([(RHO_MAIN, phi, PHI_TH_WHITE, "主(事前登録)") for phi in PHIS]
             + [(RHO_MAIN, phi, PHI_TH_RED, "主(追補・θも赤)") for phi in (0.0, 0.9, 0.98)]
             + [(CROSS[0], CROSS[1], PHI_TH_WHITE, "補助(交差)"),
                (CROSS[0], CROSS[1], PHI_TH_RED, "補助(交差・θも赤)")])
    null_rows, ctrl_rows, cal_rows = {}, {}, {}
    for rho, phi, phi_th, tag in cells:
        key = (rho, phi, phi_th)
        print(f"  ===== セル ρ={rho:+.2f} φ={phi:.2f} φ_θ={phi_th:.2f}  [{tag}] =====")
        hw, scale, d_ctrl, pil = calibrate(rho, phi, a.pilot, a.n, a.seed + 90000, phi_th=phi_th)
        cal_rows[key] = (hw, scale, d_ctrl)
        if scale is None:
            print(f"    **較正できない**（パイロット {a.pilot} 反復で CI が出たのは {pil['n_ok']} 本＜{PILOT_MIN}）"
                  f"＝このセルは判定しない\n")
            null_rows[key], ctrl_rows[key] = pil, None
            continue
        print(f"    較正：帰無パイロットの CI 半幅 中央値={hw:.6f}（判定境界）"
              f"→ 対照の真の d={d_ctrl:+.6f}（境界の {CTRL_MULT:.0f} 倍・d_scale={scale:.3f}）")
        print(f"          参考：旗44/125 の対照は d_scale={TRUE_D_SCALE_REF}"
              f"＝境界の {TRUE_D_SCALE_REF * d_ctrl / scale / hw:.0f} 倍（欠陥 #46）")
        for arm, sc, store in (("帰無(d=0)", 0.0, null_rows),
                               ("陽性対照", scale, ctrl_rows)):
            r = run_arm(rho, phi, sc, a.reps, a.n, a.seed,
                        with_iid=(sc == 0.0 and not a.no_iid), phi_th=phi_th)
            store[key] = r
            if r["unreachable"] == a.reps:
                print(f"    {arm:<10} **この n では生成不能＝{a.reps} 反復すべて**\n")
                continue
            print(f"    {arm:<10} 実現corr(θ,T)={np.mean(r['corrs']):+.3f}  "
                  f"実現ACF1(雑音)={np.mean(r['acf_eps']):+.3f}  実現ACF1(当てはめ残差)={np.mean(r['acf_res']):+.3f}  "
                  f"**ACF1(回帰子 Tc·θz)={np.mean(r['acf_x']):+.3f}**  "
                  f"VIF={np.median(r['vifs']):.1f}  真のd={r['d_true']:+.6f}")
            print(f"       CI が出た反復 {r['n_ok']}/{a.reps}（落ち：当てはめ崩壊 {r['collapse']} / "
                  f"CI不定 {r['ci_none']} / 点不足 {r['short']} / 生成不能 {r['unreachable']}）")
            if r["n_ok"]:
                print(f"       被覆率 {_fmt_rate(r['cover'], r['n_ok'])}   "
                      f"★率 {_fmt_rate(r['star'], r['n_ok'])}   ×率 {_fmt_rate(r['neg'], r['n_ok'])}")
                sd_d = float(np.std(r["ds"]))
                w = float(np.median(r["widths"]))
                print(f"       d̂ 中央値={np.median(r['ds']):+.6f}  SD(d̂)={sd_d:.6f}  CI幅 中央値={w:.6f}  "
                      f"ブート/真のばらつき={(w / 3.92) / sd_d if sd_d > 0 else float('nan'):.2f}  "
                      f"corr(ê,d̂)={np.corrcoef(r['es'], r['ds'])[0, 1]:+.3f}")
                if r["iid_ok"]:
                    print(f"       [補助・判定に使わない] block=1(iid ブート)  "
                          f"被覆率 {_fmt_rate(r['iid_cover'], r['iid_ok'])}  "
                          f"★率 {_fmt_rate(r['iid_star'], r['iid_ok'])}  "
                          f"CI幅 中央値={np.median(r['iid_widths']):.6f}")
            print()

    # ---- 事前登録の判定規則をそのまま当てる（結果を見て線を動かさない） ----
    print("  === 判定（PREREGISTRATION_step126.md の規則をそのまま当てる）===")
    print("  **主判定は ρ=−0.5 のセルのみ。交差セルは補助で、判定に入れない。**")
    print("  **同じ規則を 2 つの軸に別々に当て、両方を報告する**（追補・上記の理由）。\n")

    verdicts = {}
    for phi_th, axis_name, phis in ((PHI_TH_WHITE, "軸1＝事前登録どおり（θ の雑音は iid・回帰子は白い）", PHIS),
                                    (PHI_TH_RED, f"軸2＝追補（θ の雑音も AR(1) φ_θ={PHI_TH_RED}・回帰子が赤い）",
                                     (0.0, 0.9, 0.98))):
        print(f"  --- {axis_name} ---")
        base = null_rows.get((RHO_MAIN, 0.0, phi_th))
        base_w = float(np.median(base["widths"])) if base and base["widths"] else float("nan")
        print(f"  {'φ':>5} {'被覆率':>8} {'★率':>7} {'×率':>7} {'対照検出率':>10} {'CI幅倍率':>9}  門①  セルの判定")
        bad, undecided, decided = [], [], []
        for phi in phis:
            key = (RHO_MAIN, phi, phi_th)
            nr, cr = null_rows.get(key), ctrl_rows.get(key)
            if not nr or not nr["n_ok"] or not cr or not cr["n_ok"]:
                why = "較正不能" if cr is None else "CI が出た反復が無い"
                print(f"  {phi:>5.2f} {'—':>8} {'—':>7} {'—':>7} {'—':>10} {'—':>9}  —   判定不能（{why}）")
                undecided.append((phi, why))
                continue
            cov, star = nr["cover"] / nr["n_ok"], nr["star"] / nr["n_ok"]
            neg, power = nr["neg"] / nr["n_ok"], cr["star"] / cr["n_ok"]
            infl = np.median(nr["widths"]) / base_w if base_w else float("nan")
            gate = "通過" if power >= POWER_MIN else "**落ち**"
            if power < POWER_MIN:
                v = "門①落ち＝このセルでは帰無腕を証拠に使わない"
                undecided.append((phi, f"門①落ち（対照検出率 {power:.0%}）"))
            else:
                decided.append(phi)
                if star > STAR_RATE_MAX or cov < COVERAGE_MIN:
                    v = "**破れ**（★率>5% または 被覆率<90%）"
                    bad.append((phi, star, cov))
                else:
                    v = "○名目内"
            print(f"  {phi:>5.2f} {cov:>8.0%} {star:>7.0%} {neg:>7.0%} {power:>10.0%} {infl:>9.1f}  {gate}  {v}")

        # ★欠陥 #47 の是正：**「走らなかったセル」を「破れたセル」と同じ籠に入れない**
        #   （旗108 の欠陥 #40 と同型。動作確認で自分が踏んだ）
        print()
        if undecided:
            print(f"    **判定に使えなかったセル {len(undecided)}/{len(phis)}**："
                  + " ／ ".join(f"φ={p:.2f}（{w}）" for p, w in undecided))
            print("      **これは「破れた」ではない。この φ については何も言えない、という意味である。**")
        if not decided:
            v = "判定不能"
            print("    【この軸の判定】判定不能——判定に到達したセルが 1 つも無い。**結論を出さない。**")
        elif not bad:
            v = "○"
            print(f"    【この軸の判定】○ブロックブートは自己相関を吸収できている"
                  f"（判定に到達した {len(decided)}/{len(phis)} セルすべてで名目内）")
        elif len(bad) == 1 and bad[0][1] <= STAR_RATE_WEAK:
            v = "▲弱"
            print(f"    【この軸の判定】▲弱——破れは φ={bad[0][0]:.2f} の 1 水準のみ（★率 {bad[0][1]:.0%}）。断定しない。")
        else:
            v = "▲"
            print(f"    【この軸の判定】▲吸収し切れない——破れたセル：{[f'φ={b[0]:.2f}' for b in bad]}")
        verdicts[phi_th] = (v, bad, undecided, decided)
        print()

    v_w, v_r = verdicts[PHI_TH_WHITE][0], verdicts[PHI_TH_RED][0]
    print(f"  === 総合 ===  軸1（事前登録・白い回帰子）={v_w}   軸2（追補・赤い回帰子）={v_r}")
    if v_w == v_r:
        print("  **両軸が一致した。** 回帰子の自己相関の有無に依らず、同じ結論である。")
    else:
        print("  **両軸が食い違った。** ＝結論は θ（回帰子）の自己相関に依存する。")
        print("  **実測の土壌水分は強く自己相関するので、旗44 の実データに近いのは軸2 の方である。**")
        print("  **ただし軸2 は事前登録の後に足した軸であり、事前登録された判定は軸1 である。両方書く。**")
    if v_r == "▲" or v_r == "▲弱":
        print("  **旗44 の 18/36 は、回帰子が赤いときに水増しされている可能性がある**"
              "——**「可能性がある」までしか書かない**（実サイトの ACF1 が未確認・GATE-22）。")

    print("\n  ★実データ 36 サイトの残差 ACF1 の分布は**未確認**（この周は /mnt/hdd が無い）。")
    print("    ＝どの φ が実際に起きているかは言えない。分布の測定は GATE-22（実データのある周）。")
    print("  ★旗44 の判定規則はこの旗では変えない（旗102 の作法）。")
    print(f"\n  所要 {time.time() - t0:.0f} 秒")


if __name__ == "__main__":
    main()
