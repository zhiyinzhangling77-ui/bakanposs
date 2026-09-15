"""**旗157**：**20 年規模で `none` の帰無が緩い原因を、`Δ` の年間成分／年内成分に切って測る。**

**事前登録 `research/PREREGISTRATION_step157.md`（35 回目・判定の数を見る前に確定）。**

## 何を測るか

**旗156 が 600 本で確かめた事実**：`none`（θ→γ の真値が 0）・20 年・`thin=False` で
`Δ` の偽陽性は `θ→γH` 0.073・`θ→γLE` 0.097 で、**どちらも名目 0.05 より上**。
**緩さは LE 固有ではない**（McNemar 両側 p=0.1456）。

**本周の問い**：**その偽陽性を運んでいるのは `Δ` の年間成分（`Δ_b`）か。**

  ・`Δ_b` = `r_between`(SON) − `r_between`(MAM)
  ・`Δ_w` = `r_within`(SON)  − `r_within`(MAM)
  ・**取り分 `s` = |Δ_b| / (|Δ_b| + |Δ_w|)** ← **主判定はこれに当てる**
    （`|Δ_b|` そのものは `|Δ̂|` の大きさと交絡するので併記に落とす・事前登録 2 節）

## 門①（事前登録 3 節・**実行前に固定**）

  ・**G-J** 600 本すべてで作り直した `Δ` が旗156 の CSV と差 ≤ 5e-7
  ・**G-K** 全分解で |r_raw − (r_between + r_within)| < 1e-9
  ・**G-L** `calendar_driven` 100 本で `s` の中央値 < 0.5（誤帰属していない対照）
  ・**G-M** 自作の置換ループが旗156 の `p` を `i`=0..4 で完全再現

## 土台の道具は一行も書き換えない

`phase_asymmetry_step139.py`・`year_contrast_step144.py` は **import するだけ**。
帰無・統計量・種・`ALPHA`・`NPERM` には触っていない。

    .venv/bin/python research/year_share_none_step157.py --gates
    .venv/bin/python research/year_share_none_step157.py --all
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import phase_asymmetry_step139 as M                      # noqa: E402  **書き換えない**
from year_contrast_step144 import decompose              # noqa: E402  **書き換えない**
from moisture_control_atlas_step31 import partial_spearman  # noqa: E402
from runlog import tee_stdout                            # noqa: E402

LOGS = Path(__file__).resolve().parent / "logs"
CSV_NONE = LOGS / "step156_pairs_none_20260914_200040.csv"
CSV_CAL = LOGS / "step156_pairs_calendar_driven_20260914_203704.csv"

SEED0 = 7000          # 旗156 の `synth` の種（`M._synth_rate` と同じ）
YEARS = 20            # 旗154/156 の `SYNTH_YEARS`
DELTA_TOL = 5e-7      # G-J：CSV は小数 6 桁で保存されている
IDENT_TOL = 1e-9      # G-K：恒等式の許容
SHARE_MAX_CAL = 0.5   # G-L：`calendar_driven` の `s` 中央値の上限
N_GM = 5              # G-M：置換を突き合わせる本数
N_D5 = 60             # D-5：置換の分解を回す本数（事前登録 4 節）
NPERM_D5 = 500        # D-5：置換の本数（事前登録が「落としてよい・明記する」と書いた側）
RATIO_LO, RATIO_HI = 0.8, 1.25   # D-5 の帯（**恣意的に置いた境・事前登録 4 節に明記**）

_G: dict = {"k_worst": 0.0}


# ------------------------------------------------------------------ 中核
def _pooled(d: pd.DataFrame):
    """**`M._pool_arrays` と同じプール**に、年ラベルと季節マスクを添えて返す。

    `M._pool_arrays` は `cols`・`blocks` しか返さないが、分解には年ラベルが要る。
    **並べ方（`sort_index` と月の選び方）は `M._pool_arrays` の実装をそのまま写した**
    ——ずれると `i_sp`・`i_au` が別の行を指す。**G-M がそのずれを捕まえる。**
    """
    m = np.isin(d.index.month, M.SPRING) | np.isin(d.index.month, M.AUTUMN)
    p = d[m].sort_index()
    return p, p.index.year.to_numpy()


def _split(d: pd.DataFrame):
    return (d[np.isin(d.index.month, M.SPRING)], d[np.isin(d.index.month, M.AUTUMN)])


def _dec_sub(sub: pd.DataFrame, col: str):
    return decompose(sub[col].to_numpy(float), sub["th"].to_numpy(float),
                     sub["Rg"].to_numpy(float), sub.index.year.to_numpy())


def _dec_idx(cols, yrs, idx, col):
    """**行番号で切った部分集合を分解する**（置換の中で使う・D-5）。"""
    return decompose(cols[col][idx], cols["th"][idx], cols["Rg"][idx], yrs[idx])


def parts_for(kind: str, i: int) -> dict | None:
    """**replicate 1 本の `Δ̂`・`Δ_b`・`Δ_w`・`s`・`vshare` を出す。**"""
    d = M.synth(kind, years=YEARS, seed=SEED0 + i, thin=False)
    hat = M._delta_hat(d)
    sp, au = _split(d)
    out = {"i": i, "vshare_sp": np.nan, "vshare_au": np.nan}
    for k, col in (("h", "gH"), ("le", "gLE")):
        ds, da = _dec_sub(sp, col), _dec_sub(au, col)
        if ds is None or da is None:
            return None
        # G-K：恒等式（**季節ごと・系列ごとに全部見る**）
        for sub, dd in ((sp, ds), (au, da)):
            r_ps = float(partial_spearman(sub[col].to_numpy(float),
                                          sub["th"].to_numpy(float),
                                          [sub["Rg"].to_numpy(float)])[0])
            _G["k_worst"] = max(_G["k_worst"], abs((dd[0] + dd[1]) - r_ps))
        db, dw = da[0] - ds[0], da[1] - ds[1]
        out[f"delta_{k}"] = float(hat[k][2])
        out[f"db_{k}"], out[f"dw_{k}"] = float(db), float(dw)
        den = abs(db) + abs(dw)
        out[f"s_{k}"] = float(abs(db) / den) if den > 0 else np.nan
        out["vshare_sp"], out["vshare_au"] = float(ds[2]), float(da[2])
    return out


def table(kind: str, csv: Path, n: int) -> pd.DataFrame:
    """**CSV の `sig`・`p` と、作り直した分解を、同じ `i` の上で突き合わせた表。**"""
    ref = pd.read_csv(csv)
    rows = []
    for i in range(n):
        r = parts_for(kind, i)
        if r is None:
            continue
        m = ref[ref["i"] == i]
        if not len(m):
            continue
        for c in ("sig_h", "sig_le", "delta_h", "delta_le", "p_h", "p_le"):
            r[f"csv_{c}"] = float(m[c].iloc[0])
        rows.append(r)
        if (i + 1) % 100 == 0:
            print(f"      … {i + 1}/{n} 本")
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ 門①
def gate_j(t: pd.DataFrame) -> bool:
    print(f"\n  【G-J】作り直した Δ が旗156 の CSV と一致するか（許容 {DELTA_TOL:.0e}・全 {len(t)} 本）")
    worst = 0.0
    for k in ("h", "le"):
        e = float(np.max(np.abs(t[f"delta_{k}"] - t[f"csv_delta_{k}"])))
        worst = max(worst, e)
        print(f"    θ→γ{k.upper()}：最大差 {e:.3e}")
    ok = worst <= DELTA_TOL
    print(f"    → G-J は {'○合格' if ok else '**×不合格**'}（最大 {worst:.3e}）")
    return ok


def gate_k() -> bool:
    print(f"\n  【G-K】恒等式 r_raw = r_between + r_within（許容 {IDENT_TOL:.0e}）")
    ok = _G["k_worst"] < IDENT_TOL
    print(f"    最大差 {_G['k_worst']:.3e} → {'○合格' if ok else '**×不合格**'}")
    return ok


def gate_l(tc: pd.DataFrame) -> bool:
    print(f"\n  【G-L】`calendar_driven` {len(tc)} 本で `s` の中央値 < {SHARE_MAX_CAL}"
          f"（**年内に入れた信号を年間側に付け替えていないこと**）")
    ok = True
    for k in ("h", "le"):
        med = float(np.nanmedian(tc[f"s_{k}"]))
        hit = med < SHARE_MAX_CAL
        ok &= hit
        print(f"    θ→γ{k.upper()}：s 中央 {med:.3f}"
              f"／|Δ_b| 中央 {np.nanmedian(np.abs(tc[f'db_{k}'])):.3f}"
              f"／|Δ_w| 中央 {np.nanmedian(np.abs(tc[f'dw_{k}'])):.3f}"
              f" → {'○' if hit else '**×**'}")
    print(f"    → G-L は {'○合格' if ok else '**×不合格**'}")
    return ok


def _perm_replay(d: pd.DataFrame, nperm: int, seed: int, want_parts: bool):
    """**`M.perm_delta_p` の中身を写した置換ループ**（`rng` の消費の順番まで同じ）。

    **写し間違いを捕まえるのが G-M である。** `want_parts` のときは
    置換 1 本ごとに `Δ*_b`・`Δ*_w` も分解して返す（D-5）。
    """
    cols, blocks = M._pool_arrays(d)
    p, yrs = _pooled(d)
    obs = M._delta_hat(d)
    rng = np.random.default_rng(seed)
    ge = {"h": 0, "le": 0}
    n_ok = 0
    nb = {"h": [], "le": []}
    nw = {"h": [], "le": []}
    for _ in range(nperm):
        sp_parts, au_parts = [], []
        for idx, k in blocks:
            perm = rng.permutation(idx)
            sp_parts.append(perm[:k])
            au_parts.append(perm[k:])
        i_sp = np.sort(np.concatenate(sp_parts))
        i_au = np.sort(np.concatenate(au_parts))
        star = M._delta_from_idx(cols, i_sp, i_au)
        if not all(np.isfinite(star[k][2]) for k in ("h", "le")):
            continue
        n_ok += 1
        for k in ("h", "le"):
            ge[k] += int(abs(star[k][2]) >= abs(obs[k][2]))
        if want_parts:
            for k, col in (("h", "gH"), ("le", "gLE")):
                ds = _dec_idx(cols, yrs, i_sp, col)
                da = _dec_idx(cols, yrs, i_au, col)
                if ds is None or da is None:
                    continue
                nb[k].append(da[0] - ds[0])
                nw[k].append(da[1] - ds[1])
    out = {k: {"delta": obs[k][2], "p": (1 + ge[k]) / (n_ok + 1), "n": n_ok} for k in ("h", "le")}
    if want_parts:
        for k in ("h", "le"):
            out[k]["sd_b"] = float(np.std(nb[k])) if nb[k] else np.nan
            out[k]["sd_w"] = float(np.std(nw[k])) if nw[k] else np.nan
    return out


def gate_m(t: pd.DataFrame) -> bool:
    print(f"\n  【G-M】自作の置換ループが旗156 の `p` を再現するか"
          f"（`i`=0..{N_GM - 1}・nperm={M.NPERM}・seed=`i`）")
    worst = 0.0
    for i in range(N_GM):
        d = M.synth("none", years=YEARS, seed=SEED0 + i, thin=False)
        r = _perm_replay(d, M.NPERM, i, want_parts=False)
        m = t[t["i"] == i]
        if not len(m):
            continue
        for k in ("h", "le"):
            e = abs(r[k]["p"] - float(m[f"csv_p_{k}"].iloc[0]))
            worst = max(worst, e)
            print(f"    i={i} {k}：写し {r[k]['p']:.6f} ／ 旗156 {float(m[f'csv_p_{k}'].iloc[0]):.6f}"
                  f" ／差 {e:.1e}")
    ok = worst <= DELTA_TOL
    print(f"    → G-M は {'○合格' if ok else '**×不合格＝D-5 は読まない**'}（最大 {worst:.3e}）")
    return ok


# ------------------------------------------------------------------ 判定
def d1(t: pd.DataFrame) -> dict:
    print("\n  ===== D-1（主判定）：`s` は有意群のほうが大きいか"
          "（Mann-Whitney 両側・事前登録 4 節）=====")
    out = {}
    for k in ("h", "le"):
        s = t[f"s_{k}"].to_numpy(float)
        g = t[f"csv_sig_{k}"].to_numpy(float) > 0.5
        a, b = s[g & np.isfinite(s)], s[(~g) & np.isfinite(s)]
        if len(a) < 2 or len(b) < 2:
            print(f"    θ→γ{k.upper()}：群が小さすぎる（{len(a)} 対 {len(b)}）→ 判定しない")
            continue
        u = stats.mannwhitneyu(a, b, alternative="two-sided")
        rrb = 2.0 * float(u.statistic) / (len(a) * len(b)) - 1.0
        ma, mb = float(np.median(a)), float(np.median(b))
        print(f"    θ→γ{k.upper()}：有意群 n={len(a)} の s 中央 {ma:.3f}"
              f" ／ 非有意群 n={len(b)} の s 中央 {mb:.3f}")
        print(f"      **U={u.statistic:.0f}／両側 p = {u.pvalue:.4f}"
              f"／rank-biserial {rrb:+.3f}**"
              f"（Bonferroni 0.025 を {'下回る' if u.pvalue < 0.025 else '**下回らない**'}）")
        if u.pvalue < 0.05 and ma > mb:
            verd = "**偽陽性は Δ の年間成分が不釣り合いに運んでいる**"
        elif u.pvalue >= 0.05:
            verd = "**この標本数では、年間成分の取り分が偽陽性と結びつくとは示せなかった**"
        else:
            verd = "**見立てと逆向き＝運んでいるのは年内成分のほう**"
        print(f"      → {verd}")
        out[k] = {"p": float(u.pvalue), "rrb": rrb, "med_sig": ma, "med_ns": mb,
                  "n_sig": len(a), "n_ns": len(b)}
    return out


def d2(t: pd.DataFrame) -> None:
    print("\n  ===== D-2（併記）：|Δ_b|・|Δ_w|・s の分布 =====")
    for k in ("h", "le"):
        g = t[f"csv_sig_{k}"].to_numpy(float) > 0.5
        for nm, m in (("有意群", g), ("非有意群", ~g)):
            sub = t[m]
            if not len(sub):
                continue
            q = lambda a: (np.nanpercentile(a, 25), np.nanmedian(a), np.nanpercentile(a, 75))
            qb = q(np.abs(sub[f"db_{k}"])); qw = q(np.abs(sub[f"dw_{k}"])); qs = q(sub[f"s_{k}"])
            print(f"    θ→γ{k.upper()} {nm}（n={len(sub)}）："
                  f"|Δ_b| {qb[1]:.3f}（{qb[0]:.3f}–{qb[2]:.3f}）"
                  f"／|Δ_w| {qw[1]:.3f}（{qw[0]:.3f}–{qw[2]:.3f}）"
                  f"／s {qs[1]:.3f}（{qs[0]:.3f}–{qs[2]:.3f}）")


def d3(t: pd.DataFrame) -> dict:
    print("\n  ===== D-3（併記）：Var(Δ̂) = Var(Δ_b) + Var(Δ_w) + 2Cov =====")
    out = {}
    for k in ("h", "le"):
        db = t[f"db_{k}"].to_numpy(float); dw = t[f"dw_{k}"].to_numpy(float)
        dh = t[f"delta_{k}"].to_numpy(float)
        vb, vw = float(np.var(db)), float(np.var(dw))
        cv = float(np.cov(db, dw, bias=True)[0, 1])
        vt = float(np.var(dh))
        print(f"    θ→γ{k.upper()}：Var(Δ̂) {vt:.5f}"
              f"／Var(Δ_b) {vb:.5f}（{vb / vt:5.1%}）"
              f"／Var(Δ_w) {vw:.5f}（{vw / vt:5.1%}）"
              f"／2Cov {2 * cv:+.5f}（{2 * cv / vt:+5.1%}）"
              f"／和の誤差 {abs(vb + vw + 2 * cv - vt):.2e}")
        print(f"      SD(Δ_b) {np.sqrt(vb):.4f}／SD(Δ_w) {np.sqrt(vw):.4f}"
              f"／corr(Δ_b, Δ_w) {cv / np.sqrt(vb * vw):+.3f}")
        out[k] = {"sd_b": float(np.sqrt(vb)), "sd_w": float(np.sqrt(vw))}
    return out


def d5(obs_sd: dict, n: int, nperm: int) -> dict:
    print(f"\n  ===== D-5：置換帰無は年間成分のばらつきを再現するか"
          f"（{n} 本 × {nperm} 置換）=====")
    print(f"    **事前登録は nperm を落としてよいと書いた。落とした：{M.NPERM} → {nperm}。**")
    acc = {"h": {"b": [], "w": []}, "le": {"b": [], "w": []}}
    for i in range(n):
        d = M.synth("none", years=YEARS, seed=SEED0 + i, thin=False)
        r = _perm_replay(d, nperm, i, want_parts=True)
        for k in ("h", "le"):
            acc[k]["b"].append(r[k]["sd_b"])
            acc[k]["w"].append(r[k]["sd_w"])
        if (i + 1) % 10 == 0:
            print(f"      … {i + 1}/{n} 本")
    out = {}
    for k in ("h", "le"):
        nb = float(np.nanmedian(acc[k]["b"])); nw = float(np.nanmedian(acc[k]["w"]))
        ob, ow = obs_sd[k]["sd_b"], obs_sd[k]["sd_w"]
        rb, rw = nb / ob, nw / ow
        print(f"    θ→γ{k.upper()}：")
        print(f"      年間 Δ_b：帰無 SD 中央 {nb:.4f} ／ 600 本の SD {ob:.4f}"
              f" → **比 {rb:.3f}**")
        print(f"      年内 Δ_w：帰無 SD 中央 {nw:.4f} ／ 600 本の SD {ow:.4f}"
              f" → **比 {rw:.3f}**")
        if rb < RATIO_LO and RATIO_LO <= rw <= RATIO_HI:
            v = "**置換帰無は年間成分のばらつきだけを過小に見ている＝緩さの出どころ**"
        elif RATIO_LO <= rb <= RATIO_HI and RATIO_LO <= rw <= RATIO_HI:
            v = "**帰無は両成分を再現している＝緩さの出どころは分解の外**（見立ては外れ）"
        elif rw < RATIO_LO:
            v = "**帰無は全体として狭い＝年間／年内の切り方では説明が付かない**"
        else:
            v = "**判定しない**（事前登録の 4 行のどれにも当たらない）"
        print(f"      → {v}")
        out[k] = {"ratio_b": rb, "ratio_w": rw}
    return out


def score(gates: dict, r1: dict, r5: dict) -> None:
    """**★事前予測 H29・H30・H31**（事前登録 5 節・数を見る前に固定）。"""
    print("\n  ----- ★事前予測の勝敗（**採点の可否も事前登録に書いてある**） -----")
    base = gates["j"] and gates["k"] and gates["l"]
    if not base:
        print("    **G-J・G-K・G-L のどれかが落ちた → H29・H30・H31 はいずれも判定不能**")
        return
    h29 = any(v["p"] < 0.05 and v["med_sig"] > v["med_ns"] for v in r1.values())
    print(f"    **H29（s は有意群で大きい・h/le の少なくとも一方）→ "
          f"{'○当たり' if h29 else '**×外れ**'}**")
    if not gates["m"]:
        print("    **G-M が落ちた → H30・H31 は判定不能**")
        return
    if not r5:
        print("    **D-5 を走らせていない（`--gates` だけ）→ H30・H31 は未採点**")
        return
    h30 = all(v["ratio_b"] < 0.7 for v in r5.values())
    h31 = all(RATIO_LO <= v["ratio_w"] <= RATIO_HI for v in r5.values())
    fmt_b = "／".join("%s %.3f" % (k, v["ratio_b"]) for k, v in r5.items())
    fmt_w = "／".join("%s %.3f" % (k, v["ratio_w"]) for k, v in r5.items())
    print(f"    **H30（Δ_b の比 < 0.70）→ {'○当たり' if h30 else '**×外れ**'}（{fmt_b}）**")
    print(f"    **H31（Δ_w の比が {RATIO_LO}–{RATIO_HI}）→ "
          f"{'○当たり' if h31 else '**×外れ**'}（{fmt_w}）**")
    print(f"    **本周の勝ち星：{int(h29) + int(h30) + int(h31)}/3**")


# ------------------------------------------------------------------ 本体
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gates", action="store_true", help="門①と D-1〜D-3 まで（D-5 は走らせない）")
    ap.add_argument("--all", action="store_true", help="D-5 まで走らせる")
    ap.add_argument("--n-none", type=int, default=600)
    ap.add_argument("--n-cal", type=int, default=100)
    ap.add_argument("--n-d5", type=int, default=N_D5)
    ap.add_argument("--nperm-d5", type=int, default=NPERM_D5)
    a = ap.parse_args()
    if not (a.gates or a.all):
        ap.error("--gates か --all を指定する")
    tee_stdout("step157")
    print("**旗157：`none`(20 年) の帰無の緩さを、Δ の年間成分／年内成分に切って測る**")
    print("**事前登録 `research/PREREGISTRATION_step157.md`（35 回目・判定の数を見る前に確定）**")
    print(f"**季節：MAM {M.SPRING} ／ SON {M.AUTUMN}／年数 {YEARS}／種 {SEED0}+i**")
    print(f"**標本：{CSV_NONE.name}（`none`）・{CSV_CAL.name}（`calendar_driven`）**")

    print(f"\n  ----- 標本を作り直して分解する（`none` {a.n_none} 本）-----")
    t = table("none", CSV_NONE, a.n_none)
    print(f"  ----- 同じく `calendar_driven` {a.n_cal} 本（門① G-L 用）-----")
    tc = table("calendar_driven", CSV_CAL, a.n_cal)
    print(f"\n  **`vshare` 中央（`none` {len(t)} 本）：MAM {np.nanmedian(t['vshare_sp']):.3f}"
          f" ／ SON {np.nanmedian(t['vshare_au']):.3f}**"
          "（**θ 側の量なので h と le で同じ値になる**・前提の確認 P-4）")

    print("\n  ===== 門①（対照）——**落ちた門に依存する数は読まない**（事前登録 3 節）=====")
    g = {"j": gate_j(t), "k": gate_k(), "l": gate_l(tc), "m": gate_m(t)}
    print(f"\n  **門①：G-J {'○' if g['j'] else '**×**'}／G-K {'○' if g['k'] else '**×**'}"
          f"／G-L {'○' if g['l'] else '**×**'}／G-M {'○' if g['m'] else '**×**'}**")

    if not (g["j"] and g["k"] and g["l"]):
        print("\n  **G-J・G-K・G-L のどれかが落ちた＝分解そのものを疑う。D-1〜D-5 を読まない。**")
        score(g, {}, {})
        return 0

    r1 = d1(t)
    d2(t)
    sd = d3(t)
    r5 = {}
    if a.all:
        if g["m"]:
            r5 = d5(sd, a.n_d5, a.nperm_d5)
        else:
            print("\n  **G-M が落ちた＝置換の写しが違う。D-5 は走らせない。**")
    score(g, r1, r5)

    out = LOGS / "step157_parts_none.csv"
    try:
        t.to_csv(out, index=False)
        print(f"\n  【記録】replicate ごとの分解 → {out}")
    except Exception as e:
        print(f"\n  （**残せない**：{type(e).__name__}）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
