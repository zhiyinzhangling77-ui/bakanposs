"""旗91：**旗90 の ▲ は、日数の差を見ていただけではないか**（事前登録 step91）。

**新しい主張を作る道具ではない。自分の結論を自分で試す道具である。**

旗90 は「帯の中でも秋は反転・春は反転しない → Ts でも説明されない」と結論したが、
**帯の中では 3 サイトとも秋の方が日数が多い**（134/185・90/123・63/74）。
**旗89 では逆（反転した秋の方が少ない）だったので検出力で説明できなかったが、
旗90 にはその守りが無い。** ＝**秋を春と同じ日数に間引いて、それでも反転するかを見る。**

**事前登録 step91 で固定済み**：
  ・場所は**旗89/90 と同一**（θ・Rg のしきい値、Ts の重なり帯）＝**作り直さない**
  ・**B = 200 回**・**乱数の種 0**
  ・**p ≥ 0.80 → 旗90 を維持し強める**／**p < 0.50 → 旗90 の ▲ を取り下げる**／
    **0.50 ≤ p < 0.80 → 弱い証拠に格下げして残す**
  ・年数が 3 未満になった抽出は**数えない**（**数えなかった回数も出す**）

    python research/downsample_autumn_step91.py                # 合成で検証（既定）
    python research/downsample_autumn_step91.py --real
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaporation_regime_step36 import daily_energy, _fmt
from stratified_bowen_step89 import cell_of, test_cell, reversed_, MIN_DAYS, MIN_YEARS
from soiltemp_match_step90 import band, SPRING, AUTUMN, PLO, PHI

B, SEED = 200, 0                       # 事前登録で固定（間引き検定・判定には使わない）
HI, LO = 0.80, 0.50                    # 間引き検定の閾値（事前登録で固定）
BDIFF, SEEDD = 2000, 0                 # **差のブートストラップ**（追補で固定＝主検定）


def downsample_p(au, n_target, b=B, seed=SEED):
    """秋を n_target 日に無作為抽出して**旗90 と同一の検定**を b 回。反転した割合を返す。"""
    rng = np.random.default_rng(seed)
    n_rev = n_used = n_skip = 0
    rs = []
    for _ in range(b):
        take = au.iloc[rng.choice(len(au), size=n_target, replace=False)].sort_index()
        res = test_cell(take)
        if res is None:                 # 年数が下限未満になった抽出は**数えない**
            n_skip += 1; continue
        n_used += 1
        n_rev += int(bool(reversed_(res)))
        rs.append((res["le"][0], res["h"][0]))
    p = n_rev / n_used if n_used else float("nan")
    med = (float(np.median([r[0] for r in rs])), float(np.median([r[1] for r in rs]))) if rs else None
    return p, n_used, n_skip, med


def _pr(y, x, z):
    """偏 Spearman（旗31 の実装をそのまま使う）。**`(r, n)` を返すので r だけ取る。**"""
    from moisture_control_atlas_step31 import partial_spearman
    return partial_spearman(y, x, [z])[0]


def diff_boot(spb, aub, b=BDIFF, seed=SEEDD):
    """**Δ = r_秋 − r_春** を**年ブロック・ブートストラップ**する（追補の主検定）。

    **春と秋は同じ年から来ている**ので、**年を単位に一緒に再抽出する**——
    別々に抽出すると、**年ごとの気候変動が差に混ざる**。
    """
    yrs = np.array(sorted(set(spb.index.year) | set(aub.index.year)))
    if len(yrs) < MIN_YEARS:
        return None
    rng = np.random.default_rng(seed)
    obs = {}
    for k, col in (("le", "gLE"), ("h", "gH")):
        rs = _pr(spb[col].to_numpy(), spb["th"].to_numpy(), spb["Rg"].to_numpy())
        ra = _pr(aub[col].to_numpy(), aub["th"].to_numpy(), aub["Rg"].to_numpy())
        obs[k] = (float(ra), float(rs), float(ra - rs))
    draws = {"le": [], "h": []}
    sp_by = {y: g for y, g in spb.groupby(spb.index.year)}
    au_by = {y: g for y, g in aub.groupby(aub.index.year)}
    for _ in range(b):
        pick = rng.choice(yrs, size=len(yrs), replace=True)
        sp_i = pd.concat([sp_by[y] for y in pick if y in sp_by]) if any(y in sp_by for y in pick) else None
        au_i = pd.concat([au_by[y] for y in pick if y in au_by]) if any(y in au_by for y in pick) else None
        if sp_i is None or au_i is None or len(sp_i) < 20 or len(au_i) < 20:
            continue
        for k, col in (("le", "gLE"), ("h", "gH")):
            rs = _pr(sp_i[col].to_numpy(), sp_i["th"].to_numpy(), sp_i["Rg"].to_numpy())
            ra = _pr(au_i[col].to_numpy(), au_i["th"].to_numpy(), au_i["Rg"].to_numpy())
            if np.isfinite(rs) and np.isfinite(ra):
                draws[k].append(ra - rs)
    out = {}
    for k in ("le", "h"):
        if len(draws[k]) < 100:
            out[k] = (obs[k], None, len(draws[k])); continue
        ci = (float(np.percentile(draws[k], 2.5)), float(np.percentile(draws[k], 97.5)))
        out[k] = (obs[k], ci, len(draws[k]))
    return out


def run_site(site, d):
    print(f"\n  ━━ {site} ━━")
    if "Ts" not in d.columns:
        print("    **Ts が無い**＝対象外"); return None
    lab, _, _ = cell_of(d)
    hh = d[(lab == "θ高×Rg高") & d["Ts"].notna()]
    sp = hh[[m in SPRING for m in hh.index.month]]
    au = hh[[m in AUTUMN for m in hh.index.month]]
    if len(sp) == 0 or len(au) == 0:
        print("    春か秋が空＝対象外"); return None
    lo, hi = band(sp["Ts"].to_numpy(), au["Ts"].to_numpy())
    if hi <= lo:
        print("    帯が作れない＝対象外"); return None
    spb = sp[(sp["Ts"] >= lo) & (sp["Ts"] <= hi)]
    aub = au[(au["Ts"] >= lo) & (au["Ts"] <= hi)]
    n_sp, n_au = len(spb), len(aub)
    print(f"    帯 Ts ∈ [{lo:.1f}, {hi:.1f}]／帯の中：春 {n_sp} 日・秋 {n_au} 日")
    if n_au <= n_sp:
        print(f"    → **秋の方が多くない**＝**この弱点は存在しない**＝対象外")
        return None
    full = test_cell(aub)
    if full is None:
        print("    → 秋が下限未満＝対象外"); return None
    print(f"    間引く前の秋（{n_au}日）：{_fmt(full['le'])}  {_fmt(full['h'])}  "
          f"{'**Bowen反転**' if reversed_(full) else '反転せず'}")
    # ── **主検定（追補）**：Δ = r_秋 − r_春 を年ブロック・ブートストラップ ──
    db = diff_boot(spb, aub)
    verdict = None
    if db is None:
        print("    → 年数が下限未満＝**差の検定はできない**")
    else:
        print(f"    **主検定：Δ = r_秋 − r_春**（年ブロック・ブート {BDIFF} 回・種 {SEEDD}）")
        for k, lab in (("le", "θ→γLE|Rg"), ("h", "θ→γH|Rg")):
            (ra, rs, dd), ci, nd = db[k]
            cis = f"[{ci[0]:+.2f},{ci[1]:+.2f}]" if ci else "[CI 不能]"
            sig = "✓" if (ci and (ci[0] > 0 or ci[1] < 0)) else "·"
            print(f"      {lab:<10} 秋 {ra:+.2f}／春 {rs:+.2f}／**Δ {dd:+.2f}** {cis}{sig}"
                  f"（有効 {nd} 回）")
        (_, _, _), ci_le, _ = db["le"]
        if ci_le is None:
            print("      → **Δ_LE の CI が作れない＝判定しない**")
        elif ci_le[0] > 0 or ci_le[1] < 0:
            verdict = "差あり"
            print("      → **Δ_LE の CI が 0 を跨がない＝春と秋は本当に違う**")
        else:
            verdict = "差なし"
            print("      → **Δ_LE の CI が 0 を跨ぐ＝春と秋を区別できない**")
    # ── 参考：元の間引き検定（**判定には使わない**）──
    p, used, skip, med = downsample_p(aub, n_sp)
    if used:
        print(f"    参考・間引き検定（秋を {n_sp} 日に・{B} 回）：p = {p:.0%}"
              f"（有効 {used}／除外 {skip}）"
              f"{'／中央値 γLE %+.2f・γH %+.2f' % med if med else ''}")
        print(f"      ※**合成検証で「ほぼ必ず通る」と分かっている**ので、**判定には使わない**。")
    return verdict


def synth(kind, years=20, seed=0):
    """**弱点を検出できる道具か**を確かめる。

    ``real_diff``（**秋だけ強く反転**）を作り、**p を n の関数として出す**。

    **検証の仕方を 2 度変えた。理由を残す**：
      1. **`power_only`（春秋を同じ強さの強い効果に）** —— 失敗。
         **春も反転してしまい「春が落ちる」前提が再現できない。**
      2. **`marginal`（効果を検出の境界に）** —— 失敗。**行き過ぎて全 n で反転しなくなった。**
         そもそも**合成の帯は春161・秋177 で差が 9%** しかなく、
         **「n=161 では落ちるが n=177 では出る」という窓を作るのはほぼ不可能**だった。
      3. → **効果の強さを調整するのをやめ、`p` を `n` の関数として出す**。
         **n を下げれば p が下がる**なら、**道具は n に反応している**＝弱点を検出できる。

    **ここで設計上の大事なことに気づいた**：**実データの秋は r=+0.63〜+0.78 と大きい**ので、
    **春の日数（134/90/63）でもまず検出できてしまう**＝**この検定は落ちにくい**。
    ＝**「ほぼ必ず通る検定」であり、通っても弱い証拠にしかならない。** **そう明記して使う。**
    """
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2004-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy()
    Rg = np.clip(150 + 120 * np.sin(2 * np.pi * (doy - 80) / 365)
                 + rng.normal(0, 25, len(idx)), 5, None)
    wet = (np.exp(-0.5 * ((doy - 90) / 35.0) ** 2)
           + np.exp(-0.5 * ((doy - 230) / 35.0) ** 2))
    th = np.clip(0.14 + 0.10 * wet + rng.normal(0, 0.035, len(idx)), 0.02, 0.6)
    mon = pd.Series(idx.month)
    is_au = mon.isin(AUTUMN).to_numpy()
    is_sp = mon.isin(SPRING).to_numpy()
    # **地温のばらつきを季節で変える**——第1版は一定にしたため、
    # **帯に残るのは春の方が多くなり、実データと逆の状況しか作れなかった**
    # （実データは US-Wkg で春134・秋185＝**秋の方が多い**）。
    # 実データの Ts 分布は**春が広く（10–90 で 12.8℃）秋が狭い（7.4℃）**ので、
    # **春のばらつきを大きくする**（春は総観規模の変動が大きい）。
    # **この検定が要る状況（秋の方が帯に多く残る）を作るためであり、
    # 実データがこうだと主張するものではない。**
    ts_sd = np.where(is_sp, 9.0, np.where(is_au, 3.0, 6.0))
    Ts = 15 + 12 * np.sin(2 * np.pi * (doy - 120) / 365) + rng.normal(0, 1, len(idx)) * ts_sd
    beta = np.where(is_au, 1.6, 0.0)              # **秋だけ**強く反転する
    nz = 8.0
    gLE = Rg * (0.25 + beta * (th - th.mean())) + rng.normal(0, nz, len(idx))
    gH = Rg * (0.45 - beta * (th - th.mean())) + rng.normal(0, nz, len(idx))
    return pd.DataFrame({"th": th, "Rg": Rg, "Ts": Ts,
                         "gLE": np.clip(gLE, 0, None), "gH": np.clip(gH, 0, None)},
                        index=idx)


def main():
    p_ = argparse.ArgumentParser(description="旗91：秋を春と同じ日数に間引く")
    p_.add_argument("--real", action="store_true")
    p_.add_argument("--sites", nargs="+", default=["US-Wkg", "US-Whs", "US-SRM"])
    p_.add_argument("--qc-max", type=int, default=None)
    a = p_.parse_args()

    print("=== 旗91：旗90 の ▲ は日数の差を見ていただけではないか ===")
    print("  **新しい主張を作る検定ではない。自分の結論を自分で試す検定である。**")
    print(f"  **p ≥ {HI:.0%} → 旗90 を維持し強める**／**p < {LO:.0%} → 旗90 の ▲ を取り下げる**／")
    print(f"  **{LO:.0%} ≤ p < {HI:.0%} → 弱い証拠に格下げして残す**（事前登録で固定）")

    if not a.real:
        print("\n  【合成データで道具を検証する】**実データに触れる前に**、")
        print("  **p を n の関数として出し、道具が n に反応するか**を見る。")
        print("  **n を下げても p が下がらなければ、この道具は日数差を検出できない。**")
        d = synth("real_diff")
        lab, _, _ = cell_of(d)
        hh = d[(lab == "θ高×Rg高") & d["Ts"].notna()]
        sp = hh[[m in SPRING for m in hh.index.month]]
        au = hh[[m in AUTUMN for m in hh.index.month]]
        lo, hi = band(sp["Ts"].to_numpy(), au["Ts"].to_numpy())
        aub = au[(au["Ts"] >= lo) & (au["Ts"] <= hi)]
        full = test_cell(aub)
        print(f"\n  合成 `real_diff`（秋だけ強く反転）：帯の中の秋 {len(aub)} 日／"
              f"春 {len(sp[(sp['Ts']>=lo)&(sp['Ts']<=hi)])} 日")
        print(f"    間引く前：{_fmt(full['le'])}  {_fmt(full['h'])}")
        print(f"    {'n':>6}{'有効':>6}{'除外':>6}{'p':>8}   間引き後の中央値")
        for n_t in (len(aub), 120, 80, 50, 30, 20):
            if n_t > len(aub):
                continue
            p, used, skip, med = downsample_p(aub, n_t, b=60)
            m = f"γLE {med[0]:+.2f}／γH {med[1]:+.2f}" if med else "—"
            print(f"    {n_t:>6}{used:>6}{skip:>6}{p:>8.0%}   {m}")
        print("\n  → **n を下げると p が下がるなら、道具は日数差に反応している。**")
        print("  **だが同時に分かること**：**効果が大きければ、少ない日数でも検出できてしまう**。")
        print("  **実データの秋は r=+0.63〜+0.78 と大きい**ので、"
              "**この検定はほぼ必ず通る**——**通っても弱い証拠にしかならない。**")
        return

    out = {}
    for s in a.sites:
        try:
            d, _ = daily_energy(s, list(range(1, 13)), a.qc_max, extra=("Ts",))
        except Exception as e:
            print(f"\n  ━━ {s} ━━\n    読み込み失敗 {type(e).__name__}: {str(e)[:120]}")
            continue
        out[s] = run_site(s, d)

    print("\n  === 集計（事前登録の判定規則に当てる）===")
    for s, v in out.items():
        print(f"    {s:<9}{v or '対象外／判定しない'}")
    vals = [v for v in out.values() if v]
    n = len(vals)
    print(f"\n  === 結論（追補の規則＝Δ_LE で判定）===")
    if n < 2:
        print(f"  **判定しない**——判定できたサイトが {n} で 2 未満。")
    elif sum(v == "差あり" for v in vals) > n / 2:
        print("  **★春と秋は本当に違う**——**Δ_LE の CI が 0 を跨がない**。")
        print("  ＝**旗90 の ▲ を維持する**。「**日数の差ではない**」と言える。")
    else:
        print("  **▲春と秋を区別できない**——**Δ_LE の CI が 0 を跨ぐ**。")
        print("  ＝**旗90 の『Ts では説明されない』を取り下げ**、")
        print("     **「帯の中では春と秋を区別できない」と述べ直す。**")
    print("\n  留保（事前登録どおり）：")
    print("   ・**これは旗90 の弱点①だけを潰す**。②春の CI がそもそも広い／")
    print("     ③帯に残る春が 12–22% しかない／④モンゴルが判定できず独立クラスタ 2 つ")
    print("     ——**これらは潰れない**。**p が高くても旗90 は旗89 より弱いままである。**")
    print("   ・**無作為抽出は年の構造を壊す**。年数が減れば CI は広がる＝**向きは保守的**。")
    print("   ・**間引き検定 p は判定に使っていない**（合成検証で「ほぼ必ず通る」と判明）。")
    print("   ・**Δ が 0 と違っても「春に効果が無い」とは言えない**——**「秋より小さい」まで**。")
    print("   ・**この検定が示せるのは『秋が示した大きさの効果なら春の日数でも検出できる』まで**。")
    print("     **春の真の効果が秋より小さい（が 0 ではない）場合は、間引きでは分からない**")
    print("     ——**合成検証でこの限界に気づいた。p が高くても『春に効果が無い』とは言えない。**")


if __name__ == "__main__":
    main()
