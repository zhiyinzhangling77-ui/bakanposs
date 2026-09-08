"""旗139：**A-3 の春秋の非対称は、植生の季節進行が作っているのか**（事前登録 step139）。

**`ES-FcO`（Finca Oran・地中海の冬穀）は緑度の位相が逆である**（旗138 で実測：MOD13Q1 も
タワー CO2 も 4 月に極値）。**位相が逆の土地で非対称も逆になるなら、非対称は植生に付いている。
位相を逆にしても非対称が暦のまま（秋側）なら、非対称は植生に付いていない。**

**本ファイルは合成検証の周（旗139a）に書かれた。`ES-FcO` の θ→γH・θ→γLE は一度も計算していない。**

## 事前登録で固定済み（`PREREGISTRATION_step139.md`）
  ・**統計量**：季節ごとの `θ→γH | Rg`・`θ→γLE | Rg`（**旗88/89/91/106/107 と同一の実装**）
  ・**Δ = r_SON − r_MAM**、CI は旗91 の `diff_boot`（年ブロック・ブート）
  ・**主判定は θ→γH**。**セルで切らない**（旗138 の事実③：この地点の秋は 19 日で空になる）
  ・**下限**：各季節 60 日・3 暦年
  ・**二値化しない**（旗59 の教訓）

## ★追補 A（**実データを走らせる前**・`PREREGISTRATION_step139_amendment.md`）
**事前登録は「旗88 と同じ Pearson」と書いたが、これは事実誤りである**——
**旗88 の事前登録の統計量は「Rg を偏らせた偏 Spearman」**（`PREREGISTRATION_step88.md` 53 行目）で、
**旗89/91/106/107 の実装（`moisture_control_atlas_step31.partial_spearman`）もすべて Spearman。**
**＝「旗88 の実装をそのまま呼ぶ」を優先し、偏 Spearman を使う。判定規則・下限・主判定は変えない。**

## 門①（対照）——**実データを走らせる前に宣言済み**（事前登録より）
  ・**G1**：同じ実装を US-SRM の秋に当て、旗88（+0.81 / −0.56）・旗106（`遠い`：+0.61 / −0.39）を再現
  ・**G2**：`H`・`LE` は上向き正の W m-2 か（日中の H が正・夜間が負）
  ・**G3**：位相の前提が本検定の 3 年（2018-2020）で成り立つ（旗138 で確認済み）
  ・**G4**：季節ラベルを年内で無作為に付け替えた偽データで Δ の CI が 0 を跨ぐ
**G1〜G4 は実データの周（旗139b）で走らせる。本周（旗139a）は合成と被覆だけ。**

## 予測 H1（事前登録・**低い確信**）
**▲（植生起因ではない）＝`ES-FcO` の Δ(θ→γH) は負、反転は秋側に出る。**
**通算 16 回中 9 勝の続きとして数える。**

    python research/phase_asymmetry_step139.py                # 合成 3 種（既定）
    python research/phase_asymmetry_step139.py --coverage     # diff_boot の被覆（実データの日数で）
    python research/phase_asymmetry_step139.py --bias         # 被覆が落ちた原因の切り分け（間引き か 3 年か）
    python research/phase_asymmetry_step139.py --gates        # 門① G1〜G4（**実データの判定より先に走らせる**）
    python research/phase_asymmetry_step139.py --real         # 実データ（**旗139a では走らせない**）

## ★門①の合否しきい値（**旗139b・G1〜G4 を一度も走らせる前に書いた**）
事前登録は G1〜G4 が何を見るかを固定したが、**「どこまで合えば合格か」は数で書いていなかった**。
**実行前にここで固定する。結果を見てから動かさない。**
  ・**G1**：**旗88 の秋（9–10 月）の US-SRM** で θ→γLE = +0.81・θ→γH = −0.56 を、**符号一致かつ ±0.02 以内**
    （旗88 の記録が小数第 2 位までなので、丸め幅を許容にする）。**SON（9–11 月）も参考に印字するが合否に使わない**
    ——**旗88 の「秋」は 9–10 月であり、本検定の SON（9–11 月）とは月が違う**（**この食い違いは
    事前登録が見落としていた。合否は旗88 と同じ月で取る**）
  ・**G2**：**日中（10–14 時）の H の中央値が正**、**夜間（0–4 時）の H の中央値が負**、**日中の LE の中央値が正**
  ・**G3**：**2018/2019/2020 の 3 年とも** タワー `FC_mass` の Δ=SON−MAM **> 0**、
    かつ MOD13Q1 NDVI の Δ **< 0**（**旗138 は 7 年/3 年で測った。ここは本検定の 3 年に限って確かめる**）
  ・**G4**：**季節ラベルを年内で無作為に付け替えた偽データ 200 本のうち、
    Δ(θ→γH) の CI が 0 を跨ぐ割合 ≥ 0.90**（**旗139a の帰無で断定が 6% 出たので、10% を許容幅にする**）。
    **seed 0 の 1 本が跨ぐことも併せて印字する**
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from moisture_control_atlas_step31 import partial_spearman, _boot_ci
from downsample_autumn_step91 import diff_boot
from soiltemp_match_step90 import SPRING, AUTUMN
from stratified_bowen_step89 import MIN_DAYS, MIN_YEARS
from runlog import tee_stdout

ORAN_HH = Path("/mnt/hdd/Eddy data in Spain/Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.csv")

# **旗138 の事実④で実測した在庫**（QC<=3 かつ −200..800 W m-2 かつ 48 本中 24 本以上）
REAL_N = {"MAM": (63, 56, 56), "SON": (40, 56, 35)}     # 2018 / 2019 / 2020
FLUX_LO, FLUX_HI = -200.0, 800.0                        # 事前登録の物理範囲
QC_MAX, HALFHOUR_MIN = 3, 24                            # Foken の高品質・48 本中 24 本以上


# ------------------------------------------------------------------ 統計量
def season_r(sub: pd.DataFrame) -> dict:
    """**季節ごとの `θ→γH | Rg`・`θ→γLE | Rg`** と年ブロック CI（**旗88/89 と同一実装**）。"""
    yr = sub.index.year.to_numpy()
    out = {}
    for k, col in (("h", "gH"), ("le", "gLE")):
        out[k] = _boot_ci(sub[col].to_numpy(), sub["th"].to_numpy(),
                          [sub["Rg"].to_numpy()], blocks=yr)
    return out


def _q(a):
    return np.nanpercentile(a, [25, 50, 75])


def compute(sp: pd.DataFrame, au: pd.DataFrame, b: int | None = None,
            seed: int | None = None) -> dict | None:
    """**印字せずに `show_delta` と同じ中身を作る**（被覆の検証が同じ経路を通るように）。"""
    kw = {}
    if b is not None:
        kw["b"] = b
    if seed is not None:
        kw["seed"] = seed
    res = diff_boot(sp, au, **kw)
    if res is None:
        return None
    rs_, ra_ = season_r(sp), season_r(au)
    out = {}
    for k in ("h", "le"):
        (r_sp, ci_sp, _), (r_au, ci_au, _) = rs_[k], ra_[k]
        (ra, rsp, dd), ci_d, nb = res[k]
        out[k] = None if ci_d is None else {
            "r_sp": r_sp, "ci_sp": ci_sp, "r_au": r_au, "ci_au": ci_au,
            "delta": dd, "ci_d": ci_d, "cross": ci_d[0] <= 0 <= ci_d[1], "nb": nb}
    return out


def show_delta(tag: str, sp: pd.DataFrame, au: pd.DataFrame) -> dict | None:
    """**事前登録が「必ず併記する」と書いた 6 つを全部出す**——
    日数・年数・θ の中央値と四分位・Rg の中央値・r と CI・Δ と CI。"""
    ns, na = len(sp), len(au)
    ys = sp.index.year.nunique() if ns else 0
    ya = au.index.year.nunique() if na else 0
    print(f"    ── {tag} ──")
    for nm, sub, n, y in (("MAM", sp, ns, ys), ("SON", au, na, ya)):
        if n:
            tq, rq = _q(sub["th"]), _q(sub["Rg"])
            print(f"      {nm}: {n:>4} 日／{y} 年  θ 四分位 {tq[0]:.2f}/{tq[1]:.2f}/{tq[2]:.2f}"
                  f"   Rg 中央 {rq[1]:.1f} W m-2")
        else:
            print(f"      {nm}: 0 日")
    if ns < MIN_DAYS or na < MIN_DAYS or ys < MIN_YEARS or ya < MIN_YEARS:
        print(f"      **下限未満（各季節 {MIN_DAYS} 日・{MIN_YEARS} 年）＝判定しない**")
        return None

    out = compute(sp, au)
    if out is None:
        print("      Δ を出せない（年数が下限未満）")
        return None
    f = lambda r, c: (f"{r:+.2f} [{c[0]:+.2f},{c[1]:+.2f}]" if c else f"{r:+.2f} [CI 出ず]")
    for k, nm in (("h", "θ→γH"), ("le", "θ→γLE")):
        v = out[k]
        if v is None:
            print(f"      {nm}：Δ の CI が出ず")
            continue
        print(f"      {nm}：MAM {f(v['r_sp'], v['ci_sp'])}   SON {f(v['r_au'], v['ci_au'])}")
        print(f"        **Δ = {v['delta']:+.2f} "
              f"[{v['ci_d'][0]:+.2f},{v['ci_d'][1]:+.2f}]**（{v['nb']} 回）"
              f" → {'**0 を跨ぐ**' if v['cross'] else '**0 を跨がない**'}")
    return out


def _neg(r, ci):
    """**反転がある**＝r<0 かつ CI が 0 を跨がない（θ→γH 側の向き）。"""
    return bool(ci is not None and np.isfinite(r) and ci[1] < 0)


def verdict(res: dict | None) -> tuple[str, str]:
    """**事前登録の判定表をそのまま実装する。結果を見てから変えない。**"""
    if res is None or res.get("h") is None:
        return "判定しない", "下限未満または CI が出なかった"
    h = res["h"]
    rev_sp, rev_au = _neg(h["r_sp"], h["ci_sp"]), _neg(h["r_au"], h["ci_au"])
    if not rev_sp and not rev_au:
        return "判定しない", "どちらの季節でも r(θ→γH) の CI が 0 を跨ぐ＝この地点に反転が無い"
    if h["cross"]:
        which = "MAM" if rev_sp else "SON"
        return "○弱い証拠", f"Δ の CI は 0 を跨ぐ。反転があったのは {which} だけ（Δ は判定しない）"
    if h["delta"] > 0 and rev_sp:
        return "★植生起因と整合", "Δ>0 かつ春に反転＝位相を逆にしたら非対称も逆になった"
    if h["delta"] < 0 and rev_au:
        return "▲植生起因ではない", "Δ<0 かつ秋に反転＝位相を逆にしても非対称は暦のまま"
    return "判定しない", f"Δ の符号（{h['delta']:+.2f}）と反転の季節が表のどの行にも当たらない"


# ------------------------------------------------------------------ 合成
def synth(kind: str, years: int = 3, seed: int = 0, thin: bool = True):
    """**三つとも「期待する枝に実際に到達すること」を数値で確かめる**（旗106 の規則・5 度目）。

      ・`phase_driven`    —— **反転は「緑が濃い季節」に起きる**世界。`ES-FcO` の位相（春が緑）を入れる
                             → **Δ>0・春に反転**（★枝）
      ・`calendar_driven` —— **反転は SON にだけ起きる**世界（緑度と無関係）
                             → **Δ<0・秋に反転**（▲枝）
      ・`none`            —— **反転が無い**世界 → **両季節とも CI が 0 を跨ぐ**（判定しない枝）

    **地中海の圃場に似せる**：Rg は夏に極大、雨は冬春に多く夏は乾く、緑度は 4 月に極大で
    夏秋は刈り跡（旗138 の実測に合わせた）。**`thin=True` なら実データの日数（春 175・秋 131）まで間引く。**
    """
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2018-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy()

    # Rg：夏に極大（地中海）
    Rg = np.clip(210 + 150 * np.sin(2 * np.pi * (doy - 80) / 365)
                 + rng.normal(0, 20, len(idx)), 10, None)
    # 雨：冬春に多く、夏はほぼ降らない
    lam = np.clip(0.055 + 0.130 * np.cos(2 * np.pi * (doy - 20) / 365), 0.004, 1.0)
    P = np.where(rng.random(len(idx)) < lam, rng.gamma(1.3, 6.0, len(idx)), 0.0)
    # θ：雨の記憶（指数減衰）＋ ゆっくりした変動
    recent = np.zeros(len(idx))
    acc = 0.0
    for i, p in enumerate(P):
        acc = acc * 0.90 + p
        recent[i] = acc
    recent = recent / (recent.std() + 1e-12)
    slow = pd.Series(rng.normal(0, 1, len(idx))).rolling(45, min_periods=1).mean().to_numpy()
    slow = slow / (slow.std() + 1e-12)
    th = 10.0 + 4.5 * recent + 1.8 * slow + rng.normal(0, 0.6, len(idx))
    thz = (th - th.mean()) / th.std()

    # 緑度：4 月に極大（旗138 の実測：MOD13Q1 も FC も 4 月）、夏秋は刈り跡
    green = np.exp(-0.5 * ((doy - 105) / 48.0) ** 2)
    is_au = np.isin(idx.month, AUTUMN)

    if kind == "phase_driven":
        gate = green
    elif kind == "calendar_driven":
        gate = 1.0 * is_au
    elif kind == "none":
        gate = np.zeros(len(idx))
    else:
        raise ValueError(kind)

    avail = 0.75 * Rg
    frac = np.clip(0.45 + 0.22 * gate * thz + rng.normal(0, 0.05, len(idx)), 0.05, 0.95)
    d = pd.DataFrame({"th": th, "Rg": Rg, "gLE": avail * frac,
                      "gH": avail * (1 - frac), "green": green, "gate": gate}, index=idx)

    if thin:
        keep = []
        for nm, months, per_year in (("MAM", SPRING, REAL_N["MAM"]),
                                     ("SON", AUTUMN, REAL_N["SON"])):
            sub = d[np.isin(d.index.month, months)]
            for j, y in enumerate(sorted(sub.index.year.unique())):
                gy = sub[sub.index.year == y]
                n = min(per_year[j % len(per_year)], len(gy))
                keep.append(gy.iloc[rng.choice(len(gy), size=n, replace=False)])
        d = pd.concat(keep).sort_index()
    return d


def synth_property(kind: str, d: pd.DataFrame) -> None:
    """**各対照が名乗る性質を持っているかを印字して確かめる**（旗106 の追補の規則）。"""
    sp = d[np.isin(d.index.month, SPRING)]
    au = d[np.isin(d.index.month, AUTUMN)]
    g_sp, g_au = float(np.median(sp["green"])), float(np.median(au["green"]))
    print(f"      [性質] 合成の緑度 中央値：MAM {g_sp:.3f} ／ SON {g_au:.3f}"
          f"（{'MAM > SON＝ES-FcO の位相' if g_sp > g_au else '**位相が入っていない**'}）")
    if kind == "phase_driven":
        ok = g_sp > g_au and np.median(sp["gate"]) > np.median(au["gate"])
        print(f"      [性質] 反転の門は緑度そのもの：gate 中央 MAM {np.median(sp['gate']):.3f}"
              f" ／ SON {np.median(au['gate']):.3f} → {'○' if ok else '×'}")
    elif kind == "calendar_driven":
        r = float(np.corrcoef(d["green"].to_numpy(), d["gate"].to_numpy())[0, 1])
        print(f"      [性質] 反転の門は暦だけ：gate 中央 MAM {np.median(sp['gate']):.3f}"
              f" ／ SON {np.median(au['gate']):.3f}")
        print(f"      [性質] corr(緑度, gate) = {r:+.3f} —— **★追補 B：事前登録は『緑度と反転の有無が"
              f"独立』と書いたが、`ES-FcO` の位相を入れた世界では両者は必ず負の相関になる**")
        print("             （春に緑・秋に反転だから）。**正しい性質は「gate は暦だけで決まる」であり、"
              "それは上の 2 行が示している。**")
    else:
        print(f"      [性質] θ→γH の真値は 0：gate は全日 {float(np.max(np.abs(d['gate']))):.3f}"
              f"（**係数 0.22×gate ＝ 0**）")


# ------------------------------------------------------------------ 被覆
def _truth(kind: str, seed: int = 12345) -> dict:
    """**大標本（60 年・間引きなし）で Δ の真値を出す**（被覆の基準）。"""
    d = synth(kind, years=60, seed=seed, thin=False)
    sp = d[np.isin(d.index.month, SPRING)]
    au = d[np.isin(d.index.month, AUTUMN)]
    out = {}
    for k, col in (("h", "gH"), ("le", "gLE")):
        rsp = partial_spearman(sp[col].to_numpy(), sp["th"].to_numpy(), [sp["Rg"].to_numpy()])[0]
        rau = partial_spearman(au[col].to_numpy(), au["th"].to_numpy(), [au["Rg"].to_numpy()])[0]
        out[k] = (float(rsp), float(rau), float(rau - rsp))
    return out


def coverage(kinds, reps: int, b: int) -> None:
    """**実データの日数（春 175・秋 131・3 年）で `diff_boot` の CI が名目 95% を保つか。**

    **旗135 で、この確認を怠ると被覆が壊れうると分かった。**
    **保たなければ、CI ではなく符号の一致だけを報告する**（事前登録どおり）。
    """
    print("\n  【被覆の検証】**実データの日数（春 175・秋 131・3 年）で `diff_boot` の CI を測る**")
    print(f"  **replicate {reps} 回 × ブート {b} 回**（既定の 2000 回から下げている。"
          "**年が 3 つしかないので年ブロックの取り方は高々 10 通りで、b を増やしても粒度は変わらない**）")
    for kind in kinds:
        tr = _truth(kind)
        print(f"\n  ===== `{kind}` =====")
        print(f"    真値（60 年の大標本）：θ→γH  MAM {tr['h'][0]:+.3f} / SON {tr['h'][1]:+.3f}"
              f" / **Δ {tr['h'][2]:+.3f}**")
        print(f"                          θ→γLE MAM {tr['le'][0]:+.3f} / SON {tr['le'][1]:+.3f}"
              f" / **Δ {tr['le'][2]:+.3f}**")
        hit = {"h": 0, "le": 0}
        nb_ok = {"h": 0, "le": 0}
        sign_ok = {"h": 0, "le": 0}
        width = {"h": [], "le": []}
        dhat = {"h": [], "le": []}
        # **季節ごとの CI の被覆**（**判定表の「反転がある」はこちらで決まる**）
        s_hit = {"h_sp": 0, "h_au": 0, "le_sp": 0, "le_au": 0}
        s_n = {"h_sp": 0, "h_au": 0, "le_sp": 0, "le_au": 0}
        rev = {"sp": 0, "au": 0}          # θ→γH が「反転あり」と判定された回数
        verd = {}
        for i in range(reps):
            d = synth(kind, years=3, seed=1000 + i, thin=True)
            sp = d[np.isin(d.index.month, SPRING)]
            au = d[np.isin(d.index.month, AUTUMN)]
            res = compute(sp, au, b=b, seed=1000 + i)
            if res is None:
                continue
            v, _ = verdict(res)
            verd[v] = verd.get(v, 0) + 1
            if res["h"] is not None:
                rev["sp"] += int(_neg(res["h"]["r_sp"], res["h"]["ci_sp"]))
                rev["au"] += int(_neg(res["h"]["r_au"], res["h"]["ci_au"]))
            for k in ("h", "le"):
                if res[k] is None:
                    continue
                nb_ok[k] += 1
                ci, dd = res[k]["ci_d"], res[k]["delta"]
                width[k].append(ci[1] - ci[0])
                dhat[k].append(dd)
                if ci[0] <= tr[k][2] <= ci[1]:
                    hit[k] += 1
                if np.sign(dd) == np.sign(tr[k][2]) or tr[k][2] == 0:
                    sign_ok[k] += 1
                for j, (season, truth) in enumerate((("sp", tr[k][0]), ("au", tr[k][1]))):
                    ci_s = res[k]["ci_sp" if season == "sp" else "ci_au"]
                    if ci_s is None:
                        continue
                    s_n[f"{k}_{season}"] += 1
                    s_hit[f"{k}_{season}"] += int(ci_s[0] <= truth <= ci_s[1])
        for k, nm in (("h", "θ→γH"), ("le", "θ→γLE")):
            if not nb_ok[k]:
                print(f"    {nm}：CI が一度も出なかった")
                continue
            sd = float(np.std(dhat[k], ddof=1))
            print(f"    {nm} の Δ：**被覆 {hit[k]/nb_ok[k]:.3f}**（{hit[k]}/{nb_ok[k]}・名目 0.95）"
                  f"／CI 幅 中央 {np.median(width[k]):.3f}"
                  f"／**Δ̂ の実際のばらつき 4×SD = {4*sd:.3f}**"
                  f"／点推定の符号一致 {sign_ok[k]/nb_ok[k]:.3f}")
            for season, jp in (("sp", "MAM"), ("au", "SON")):
                key = f"{k}_{season}"
                if s_n[key]:
                    print(f"      {nm} の {jp} 単季 CI：**被覆 {s_hit[key]/s_n[key]:.3f}**"
                          f"（{s_hit[key]}/{s_n[key]}）")
        print(f"    θ→γH で「反転あり」と判定した割合：MAM {rev['sp']/max(reps,1):.3f}"
              f" ／ SON {rev['au']/max(reps,1):.3f}")
        print("    **判定表が返した結論の分布**：" +
              " ／ ".join(f"{k} {n}/{reps}" for k, n in sorted(verd.items(), key=lambda x: -x[1])))


def _delta_hat(d: pd.DataFrame) -> dict:
    """**点推定だけ**（ブートしない）。`Δ = r_SON − r_MAM` を θ→γH・θ→γLE で。"""
    sp = d[np.isin(d.index.month, SPRING)]
    au = d[np.isin(d.index.month, AUTUMN)]
    out = {}
    for k, col in (("h", "gH"), ("le", "gLE")):
        rsp = partial_spearman(sp[col].to_numpy(), sp["th"].to_numpy(), [sp["Rg"].to_numpy()])[0]
        rau = partial_spearman(au[col].to_numpy(), au["th"].to_numpy(), [au["Rg"].to_numpy()])[0]
        out[k] = (float(rsp), float(rau), float(rau - rsp))
    return out


def bias(kinds, reps: int) -> None:
    """**被覆が名目を割った原因を切り分ける**（旗139a・被覆 0.46〜0.755 を見て足した診断）。

    **`calendar_driven` の Δ̂ は −0.55 前後で、60 年の真値 −0.869 から CI の外まで離れていた。**
    **CI の幅が狭すぎる（分散の過小評価）のか、Δ̂ そのものが真値からずれている（偏り）のかを分ける。**
    **同じ seed で `thin=False`（3 年・全日）と `thin=True`（3 年・実データの日数）を対で走らせ、
    60 年の真値と並べる。** **判定規則には一切触れない。診断の印字だけである。**
    """
    print("\n  【被覆が落ちた原因の切り分け】**同じ seed で「3 年・全日」と「3 年・間引き」を対で出し、"
          "60 年の真値と並べる**")
    print(f"  **replicate {reps} 回・ブートしない（点推定だけ）。判定規則には触れない。**")
    for kind in kinds:
        tr = _truth(kind)
        full, thin_ = {"h": [], "le": []}, {"h": [], "le": []}
        for i in range(reps):
            for store, th_ in ((full, False), (thin_, True)):
                d = synth(kind, years=3, seed=1000 + i, thin=th_)
                v = _delta_hat(d)
                for k in ("h", "le"):
                    store[k].append(v[k][2])
        print(f"\n  ===== `{kind}` =====")
        for k, nm in (("h", "θ→γH"), ("le", "θ→γLE")):
            t = tr[k][2]
            mf, mt = float(np.median(full[k])), float(np.median(thin_[k]))
            print(f"    {nm} の Δ：**60 年の真値 {t:+.3f}**"
                  f" ／ 3 年・全日（MAM 276・SON 273）**中央 {mf:+.3f}**（SD {np.std(full[k], ddof=1):.3f}）"
                  f" ／ 3 年・間引き（175・131）**中央 {mt:+.3f}**（SD {np.std(thin_[k], ddof=1):.3f}）")
            print(f"      → 真値からのずれ：全日 {mf - t:+.3f} ／ 間引き {mt - t:+.3f}"
                  f" ／ 間引きが足したぶん {mt - mf:+.3f}")
            # **表示用の 1 走（seed 0）が分布のどこにいるか**——
            # **合成のまとめに印字される数は 1 標本であって、設計の効果量ではない。**
            d0 = synth(kind, years=3, seed=0, thin=True)
            v0 = _delta_hat(d0)[k][2]
            pct = 100.0 * float(np.mean(np.array(thin_[k]) <= v0))
            print(f"      表示用の 1 走（seed 0）の Δ̂ = {v0:+.3f} ＝ 間引き分布の下から {pct:.0f} %"
                  f"（**まとめに出る数は 1 標本であって設計の効果量ではない**）")


def cidiag(kinds, reps: int, b: int) -> None:
    """**被覆が落ちた残りの容疑者＝CI そのものを測る**（旗139a・`--bias` で偏りが小さいと分かった後）。

    **`--bias` は Δ̂ の偏りが小さいことを示した**（`calendar_driven` の θ→γH で +0.073）。
    **それだけでは被覆 0.46 は説明できない。** **残る容疑者は CI の位置と幅である**——
    **年ブロックが 3 つしかないので、ブートの再抽出は高々 10 通りの重複組合せしか作れず、
    percentile CI は Δ̂ を中心に乗らないかもしれない。** 次の 3 つを測る：

      ・**CI の中心 − Δ̂**（0 から離れていれば、CI は点推定の周りに乗っていない）
      ・**CI の半幅 ÷ Δ̂ の実際の SD**（1.96 より小さければ狭すぎる）
      ・**外した向き**（真値が CI より上か下か。片側に偏れば位置、両側なら幅）

    **診断の印字だけである。判定規則には触れない。**
    """
    print("\n  【CI そのものの診断】**偏りだけでは被覆 0.46 を説明できないので、CI の位置と幅を測る**")
    print(f"  **replicate {reps} 回 × ブート {b} 回**")
    for kind in kinds:
        tr = _truth(kind)
        print(f"\n  ===== `{kind}` =====")
        for k, nm in (("h", "θ→γH"), ("le", "θ→γLE")):
            off, half, dh = [], [], []
            lo_miss = hi_miss = n = 0
            for i in range(reps):
                d = synth(kind, years=3, seed=1000 + i, thin=True)
                sp = d[np.isin(d.index.month, SPRING)]
                au = d[np.isin(d.index.month, AUTUMN)]
                res = compute(sp, au, b=b, seed=1000 + i)
                if res is None or res[k] is None:
                    continue
                ci, dd, t = res[k]["ci_d"], res[k]["delta"], tr[k][2]
                n += 1
                off.append(0.5 * (ci[0] + ci[1]) - dd)
                half.append(0.5 * (ci[1] - ci[0]))
                dh.append(dd)
                if t < ci[0]:
                    lo_miss += 1
                elif t > ci[1]:
                    hi_miss += 1
            if not n:
                print(f"    {nm}：CI が一度も出なかった")
                continue
            sd = float(np.std(dh, ddof=1))
            mh = float(np.median(half))
            print(f"    {nm}：**CI の中心 − Δ̂ の中央 {np.median(off):+.4f}**"
                  f"（0 なら CI は点推定に乗っている）")
            print(f"      **CI 半幅 中央 {mh:.3f} ÷ Δ̂ の SD {sd:.3f} = {mh/sd:.2f}**"
                  f"（**95% には 1.96 が要る**）")
            print(f"      外した向き：真値が CI より**下 {lo_miss}/{n}**・**上 {hi_miss}/{n}**"
                  f"（片側なら位置のずれ・両側なら幅の不足）")


# ------------------------------------------------------------------ 実データ
def load_oran_daily() -> pd.DataFrame:
    """`ES-FcO` の 30 分を日次に落とす（**事前登録の手順 1〜2 をそのまま**）。

    **`H_QC<=3` かつ `LE_QC<=3` かつ両者が −200..800 W m-2。日次は 48 本中 24 本以上。**
    **落とした本数と日数を必ず印字する。**  **旗139a では走らせていない。**
    """
    use = ["TIMESTAMP", "SWC_1_1_1", "SW_IN", "H", "H_QC", "LE", "LE_QC"]
    d = pd.read_csv(ORAN_HH, usecols=use, low_memory=False)
    # ★欠陥 #67（旗139b で G4 が捕まえた）：`format` を渡さないと pandas は先頭行の
    # `2018/01/01`（時刻なし）から `%Y/%m/%d` を推定し、**時刻つきの 30 分値を全部 NaT にする**
    # ——52,606 本中 921 本（真夜中だけ）しか残らなかった。**旗138 の実装は `format="mixed"` で正しい。**
    d["ts"] = pd.to_datetime(d["TIMESTAMP"], format="mixed", errors="coerce")
    for c in use[1:]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    n_file = len(d)
    d = d[d["ts"].notna()]
    n_raw = len(d)
    print(f"    [在庫] ファイル {n_file:,} 行 → 時刻が読めた {n_raw:,} 本"
          f"（読めなかった {n_file - n_raw:,} 本は **TIMESTAMP が空**）")
    good = d[(d["H_QC"] <= QC_MAX) & d["H"].between(FLUX_LO, FLUX_HI)
             & (d["LE_QC"] <= QC_MAX) & d["LE"].between(FLUX_LO, FLUX_HI)
             & d["SWC_1_1_1"].notna() & d["SW_IN"].notna()].copy()
    print(f"    [在庫] 30 分 {n_raw:,} 本 → QC<=3 かつ {FLUX_LO:.0f}..{FLUX_HI:.0f} W m-2 "
          f"かつ θ・Rg あり: {len(good):,} 本（{100*len(good)/max(n_raw,1):.1f}%）")
    good["date"] = good["ts"].dt.floor("D")
    day = good.groupby("date").agg(th=("SWC_1_1_1", "mean"), Rg=("SW_IN", "mean"),
                                   gH=("H", "mean"), gLE=("LE", "mean"), n=("H", "size"))
    n_day_all = len(day)
    day = day[day["n"] >= HALFHOUR_MIN]
    print(f"    [在庫] 日 {n_day_all} → 48 本中 {HALFHOUR_MIN} 本以上: {len(day)} 日")
    return day.drop(columns="n")


# ------------------------------------------------------------------ 門①（対照）
G88_SRM_AUTUMN = {"le": +0.81, "h": -0.56}      # 旗88 の記録（US-SRM・秋 9–10 月）
G1_TOL = 0.02                                    # 丸め幅（**実行前に固定**）
G4_CROSS_MIN = 0.90                              # 偽の季節で CI が 0 を跨ぐ割合の下限（**実行前に固定**）


def gate_g1() -> bool:
    """**G1（検出器の健全性）**：同じ実装を US-SRM の秋に当て、旗88 の値を再現する。

    **旗88 の「秋」は 9–10 月である**（`FLAGS_LOG.md` 旗88 の表）。
    **本検定の SON は 9–11 月**なので**月が違う**。**合否は旗88 と同じ 9–10 月で取り、
    SON は参考として印字する**（**この食い違いは事前登録が見落としていた**）。
    """
    from evaporation_regime_step36 import daily_energy
    ok = None
    for months, label, judge in (((9, 10), "旗88 の秋（9–10 月）", True),
                                 (list(AUTUMN), "SON（9–11 月・本検定の季節）", False)):
        try:
            d, nyr = daily_energy("US-SRM", list(months), None)
        except Exception as e:                                   # noqa: BLE001
            print(f"    G1 {label}: 読み込み失敗 {type(e).__name__}: {str(e)[:120]}")
            if judge:
                ok = False
            continue
        yr = d.index.year.to_numpy()
        got = {}
        for k, col in (("le", "gLE"), ("h", "gH")):
            got[k] = _boot_ci(d[col].to_numpy(), d["th"].to_numpy(),
                              [d["Rg"].to_numpy()], blocks=yr)
        print(f"    G1 {label}: {len(d)} 日／{d.index.year.nunique()} 年（健全年 {nyr}）")
        for k, nm, ref in (("le", "θ→γLE", G88_SRM_AUTUMN["le"]),
                           ("h", "θ→γH", G88_SRM_AUTUMN["h"])):
            r, ci, n = got[k]
            cis = f"[{ci[0]:+.2f},{ci[1]:+.2f}]" if ci else "[CI 出ず]"
            mark = ""
            if judge:
                hit = np.isfinite(r) and np.sign(r) == np.sign(ref) and abs(r - ref) <= G1_TOL
                ok = hit if ok is None else (ok and hit)
                mark = f"  旗88 = {ref:+.2f}／差 {r - ref:+.3f} → {'○' if hit else '**×**'}"
            print(f"      {nm} = {r:+.3f} {cis}{mark}")
    print(f"    → G1 は {'○合格' if ok else '**×不合格**'}"
          f"（許容 ±{G1_TOL:.2f}・符号一致・**実行前に固定**）")
    return bool(ok)


def gate_g2() -> bool:
    """**G2（符号と単位）**：`H`・`LE` は上向き正の W m-2 か。**日中の H が正・夜間が負**。"""
    use = ["TIMESTAMP", "H", "H_QC", "LE", "LE_QC"]
    d = pd.read_csv(ORAN_HH, usecols=use, low_memory=False)
    ts = pd.to_datetime(d["TIMESTAMP"], format="mixed", errors="coerce")
    for c in ("H", "H_QC", "LE", "LE_QC"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["hour"] = ts.dt.hour
    g = d[(d["H_QC"] <= QC_MAX) & d["H"].between(FLUX_LO, FLUX_HI)
          & (d["LE_QC"] <= QC_MAX) & d["LE"].between(FLUX_LO, FLUX_HI)]
    day = g[g["hour"].between(10, 13)]
    night = g[g["hour"].between(0, 3)]
    mh_d, mh_n = float(np.median(day["H"])), float(np.median(night["H"]))
    ml_d, ml_n = float(np.median(day["LE"])), float(np.median(night["LE"]))
    print(f"    G2 日中 10–14 時（n={len(day):,}）: H 中央 {mh_d:+.1f} ／ LE 中央 {ml_d:+.1f} W m-2")
    print(f"    G2 夜間 0–4 時 （n={len(night):,}）: H 中央 {mh_n:+.1f} ／ LE 中央 {ml_n:+.1f} W m-2")
    hourly = g.groupby("hour")["H"].median()
    print("      H の時刻別中央値: "
          + " ".join(f"{int(h):02d}:{v:+.0f}" for h, v in hourly.items()))
    ok = mh_d > 0 and mh_n < 0 and ml_d > 0
    print(f"    → G2 は {'○合格' if ok else '**×不合格**'}（日中 H>0・夜間 H<0・日中 LE>0）")
    return bool(ok)


def gate_g3() -> bool:
    """**G3（位相の前提）**：旗138 の事実①が本検定の 3 年（2018-2020）で成り立つか。

    **タワー `FC_mass` は 3 年とも Δ=SON−MAM>0（春の吸収が強い）**、
    **MOD13Q1 NDVI は 3 年とも Δ<0（春が緑）**。**旗138 は 7 年/3 年で測った。ここは 3 年に絞る。**
    """
    from premise_phase_step138 import load_oran_flux, load_oran_ndvi_mod13

    def _per_year(df, col, years=(2018, 2019, 2020)):
        d = df[df["year"].isin(years)].copy()
        d["season"] = np.where(np.isin(d["month"], SPRING), "MAM",
                               np.where(np.isin(d["month"], AUTUMN), "SON", None))
        d = d[d["season"].isin(("MAM", "SON"))]
        out = {}
        for y, g in d.groupby("year"):
            m = g[g["season"] == "MAM"][col]
            s = g[g["season"] == "SON"][col]
            out[int(y)] = (float(s.mean() - m.mean()) if len(m) and len(s) else np.nan,
                           len(m), len(s))
        return out

    fc = _per_year(load_oran_flux(), "fc")
    nd = _per_year(load_oran_ndvi_mod13(), "ndvi")
    ok_fc = bool(fc) and all(np.isfinite(v[0]) and v[0] > 0 for v in fc.values())
    ok_nd = bool(nd) and all(np.isfinite(v[0]) and v[0] < 0 for v in nd.values())
    for nm, tab, want in (("タワー FC_mass（Δ>0 を要求）", fc, ">0"),
                          ("MOD13Q1 NDVI（Δ<0 を要求）", nd, "<0")):
        s = " ／ ".join(f"{y}: Δ={v[0]:+.3f}(n={v[1]}/{v[2]})" for y, v in sorted(tab.items()))
        print(f"    G3 {nm}: {s}")
    ok = ok_fc and ok_nd
    print(f"    → G3 は {'○合格' if ok else '**×不合格**'}"
          f"（FC 3/3 {'○' if ok_fc else '×'}・NDVI 3/3 {'○' if ok_nd else '×'}）")
    return bool(ok)


def gate_g4(d: pd.DataFrame, reps: int, b: int) -> bool:
    """**G4（偽の季節）**：季節ラベルを年内で無作為に付け替えた偽データで Δ の CI が 0 を跨ぐか。

    **同じ年の MAM+SON の日をプールし、季節の枚数を保ったまま札を振り直す**
    ——**日数の非対称（春 175・秋 131）はそのまま残り、季節の中身だけが消える。**
    **跨がなければ、Δ は季節ではなく日数の非対称を拾っている。**
    """
    pool = d[np.isin(d.index.month, SPRING) | np.isin(d.index.month, AUTUMN)]
    n_sp = {y: int(np.isin(g.index.month, SPRING).sum()) for y, g in pool.groupby(pool.index.year)}
    cross = n = 0
    first = None
    for i in range(reps):
        rng = np.random.default_rng(20000 + i)
        sp_parts, au_parts = [], []
        for y, g in pool.groupby(pool.index.year):
            perm = rng.permutation(len(g))
            k = n_sp[y]
            sp_parts.append(g.iloc[perm[:k]])
            au_parts.append(g.iloc[perm[k:]])
        sp, au = pd.concat(sp_parts).sort_index(), pd.concat(au_parts).sort_index()
        res = compute(sp, au, b=b, seed=20000 + i)
        if res is None or res["h"] is None:
            continue
        n += 1
        cross += int(res["h"]["cross"])
        if first is None:
            first = res["h"]
    if not n:
        print("    G4: CI が一度も出なかった → **×不合格**")
        return False
    frac = cross / n
    print(f"    G4 偽の季節 {n} 本（要求どおり年内で付け替え）："
          f"**Δ(θ→γH) の CI が 0 を跨いだ割合 {frac:.3f}**（{cross}/{n}）")
    print(f"      1 本目（seed 20000）：Δ = {first['delta']:+.3f} "
          f"[{first['ci_d'][0]:+.3f},{first['ci_d'][1]:+.3f}] → "
          f"{'**0 を跨ぐ**' if first['cross'] else '**跨がない**'}")
    ok = frac >= G4_CROSS_MIN
    print(f"    → G4 は {'○合格' if ok else '**×不合格**'}"
          f"（下限 {G4_CROSS_MIN:.2f}・**実行前に固定**）")
    return bool(ok)


def run_gates(reps: int, b: int) -> bool:
    """**門① G1〜G4 を全部走らせる。ここで落ちたら実データの判定は読まない。**"""
    print("\n  【門①（対照）】**事前登録が実データの前に宣言した 4 本**"
          "——**合否のしきい値は走らせる前に docstring へ固定した**")
    print("\n  ===== G1（検出器の健全性）=====")
    g1 = gate_g1()
    print("\n  ===== G2（符号と単位）=====")
    g2 = gate_g2()
    print("\n  ===== G3（位相の前提が 2018-2020 で成り立つか）=====")
    g3 = gate_g3()
    print("\n  ===== G4（偽の季節）=====")
    d = load_oran_daily()
    g4 = gate_g4(d, reps, b)
    allok = g1 and g2 and g3 and g4
    print("\n  === 門①のまとめ ===")
    for nm, v in (("G1 検出器", g1), ("G2 符号と単位", g2), ("G3 位相の前提", g3), ("G4 偽の季節", g4)):
        print(f"    {nm:<14}{'○合格' if v else '**×不合格**'}")
    print(f"\n  **4 本とも通ったか：{'○＝実データの判定を読んでよい' if allok else '**×＝読まない**'}**")
    return allok


def run_real() -> None:
    print("\n  【実データ】`ES-FcO`（Finca Oran・2018-2020）")
    d = load_oran_daily()
    sp = d[np.isin(d.index.month, SPRING)]
    au = d[np.isin(d.index.month, AUTUMN)]
    res = show_delta("層なし（GATE-26 待ち）", sp, au)
    v, why = verdict(res)
    print(f"\n    → **{v}**（{why}）")
    print("\n    **これは層で揃えていない春秋の差である。旗107 が示したとおり、この差は")
    print("    雨からの日数で説明されうる。層で揃えた判定は GATE-26（`ES-FcO` の雨）が")
    print("    解けてから行う。**")


# ------------------------------------------------------------------ main
def main() -> int:
    ap = argparse.ArgumentParser(description="旗139：春秋の非対称は植生の季節進行が作っているのか")
    ap.add_argument("--real", action="store_true", help="実データ（旗139b で使う）")
    ap.add_argument("--gates", action="store_true", help="門① G1〜G4（実データの判定より先に走らせる）")
    ap.add_argument("--gate-reps", type=int, default=200, help="G4 の偽の季節の本数")
    ap.add_argument("--coverage", action="store_true", help="diff_boot の被覆を測る")
    ap.add_argument("--bias", action="store_true", help="被覆が落ちた原因の切り分け（間引き か 3 年か）")
    ap.add_argument("--cidiag", action="store_true", help="CI の位置と幅を測る（偏りでは説明できない残り）")
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--boot", type=int, default=600)
    a = ap.parse_args()

    tee_stdout("step139")
    print("=== 旗139：A-3 の春秋の非対称は、植生の季節進行が作っているのか ===")
    print("  **統計量は偏 Spearman（旗88/89/91/106/107 と同一）**——**★追補 A**：")
    print("  事前登録の『旗88 と同じ Pearson』は事実誤り。旗88 の事前登録も実装も Spearman である。")
    print("  **主判定は θ→γH の Δ = r_SON − r_MAM。セルで切らない（旗138 の事実③）。**")

    if a.gates:
        run_gates(a.gate_reps, a.boot)
        return 0
    if a.real:
        run_real()
        return 0

    want = {"phase_driven": "**Δ>0・春に反転 → ★植生起因と整合**",
            "calendar_driven": "**Δ<0・秋に反転 → ▲植生起因ではない**",
            "none": "**両季節とも CI が 0 を跨ぐ → 判定しない**"}
    if a.bias:
        bias(list(want), a.reps)
        return 0
    if a.cidiag:
        cidiag(list(want), a.reps, a.boot)
        return 0
    if a.coverage:
        coverage(list(want), a.reps, a.boot)
        return 0

    print("\n  【合成データで検証する】**三つとも期待する枝に実際に到達するかを数値で見る。**")
    print("  **日数は実データに合わせて間引く**（春 175＝63/56/56・秋 131＝40/56/35・旗138 の事実④）。")
    got = {}
    for kind, w in want.items():
        print(f"\n  ===== 合成 `{kind}` —— 期待：{w} =====")
        d = synth(kind, years=3, seed=0, thin=True)
        synth_property(kind, d)
        sp = d[np.isin(d.index.month, SPRING)]
        au = d[np.isin(d.index.month, AUTUMN)]
        res = show_delta(kind, sp, au)
        v, why = verdict(res)
        got[kind] = v
        print(f"      → **{v}**（{why}）")

    exp = {"phase_driven": "★植生起因と整合", "calendar_driven": "▲植生起因ではない",
           "none": "判定しない"}
    print("\n  === 合成のまとめ ===")
    allok = True
    for kind in want:
        ok = got[kind] == exp[kind]
        allok &= ok
        print(f"    {kind:<16}期待 {exp[kind]:<12}実際 {got[kind]:<12}{'○' if ok else '**×**'}")
    print(f"\n  **三つとも期待どおりに到達したか：{'○＝この道具は実データに使える' if allok else '**×＝使えない。設計を直す**'}**")
    print("  **被覆は `--coverage` で別に測る。被覆が保たれなければ、CI ではなく符号の一致だけを報告する。**")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
