"""旗99：**~4 日メモリは、どの空間スケールの現象か**（事前登録 step99）。

**機構（基質の転流・微生物・根）は手元で区別できない。だが空間スケールは測れる**——
**同じ日の残差が、チャンバー間でどれだけ揃うか**を見る。

| 揃い方 | 意味 |
|---|---|
| **同一地点で揃い、地点間では揃わない** | **地点スケール**（気象・土壌水文）＝**局所生物ではない** |
| **同一地点でも揃わない** | **チャンバー・スケール**（微生物・根の微小生息場所） |
| **地点間でも揃う** | **広域の気象**——**我々が測っていない気象量** |

**事前登録 step99 で固定済み**：
  ・**残差は旗74 のテンソルビン＋外挿**（**形も交互作用も仮定しない最も豊かな駆動**）
    ——**同一地点のチャンバーは似た Ts・SM を持つ**ので、**誤特定が共通なら残差は自明に揃う**。
    **それを最大限つぶす。**
  ・**プラセボは ±30/60/90 日の多重シフト**（旗71：単一シフトは帰無になっていなかった）。
    **「6 通りの最大値を上回る」**を「揃う」の条件にする。
  ・**KAYE 80% 以上 かつ 地点間 30% 以下 → 地点スケール**／**KAYE 30% 以下 → チャンバー・スケール**／
    **両方 80% 以上 → 広域の気象**／**それ以外 → 中間と記す**。
  ・**同一地点の腕は n=1 地点（KAYE）**。**「28 例で確かめた」とは書かない。**

    python research/chamber_sync_step99.py                    # 合成で検証（既定）
    python research/chamber_sync_step99.py --real --cosore-dir /mnt/hdd/cosore-0.7.0
"""
from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosore_memory_step40 import load_cosore
from colocate_step51 import haversine
from model_richness_step74 import design, residuals
from vpd_match_step96 import spearman
from runlog import tee_stdout      # **出力を最初からファイルに残す**（旗110 の反省）

MIN_DAYS, MIN_YEARS = 60, 3
SHIFTS = (-90, -60, -30, 30, 60, 90)      # プラセボ（事前登録で固定）
HI, LO = 0.80, 0.30                       # 判定の閾値（事前登録で固定）
SAME_KM = 1.0


def chamber_residual(daily):
    """**旗74 のテンソルビン＋外挿残差**を、このチャンバー単体に当てる。"""
    if daily is None or "Tsoil" not in daily:
        return None
    y = daily["Rs"].to_numpy()
    T = daily["Tsoil"].to_numpy()
    W = daily["SM"].to_numpy() if "SM" in daily else None
    X = design("テンソルビン", T, W)
    if X is None:
        return None
    with np.errstate(divide="ignore", invalid="ignore"):
        ly = np.log(np.where(y > 0, y, np.nan))
    r = residuals(X, ly, True)
    if r is None or np.isfinite(r).sum() < MIN_DAYS:
        return None
    return pd.Series(r, index=daily.index).dropna()


def pair_sync(a, b):
    """**同じ日の相関 r_obs** と、**±30/60/90 日ずらしたプラセボの最大値**を返す。"""
    common = a.index.intersection(b.index)
    if len(common) < MIN_DAYS or pd.Index(common).year.nunique() < MIN_YEARS:
        return None
    r_obs = spearman(a.reindex(common).to_numpy(), b.reindex(common).to_numpy())
    pl = []
    for s in SHIFTS:
        bs = b.copy()
        bs.index = bs.index + pd.Timedelta(days=s)
        c2 = a.index.intersection(bs.index)
        if len(c2) < MIN_DAYS:
            continue
        v = spearman(a.reindex(c2).to_numpy(), bs.reindex(c2).to_numpy())
        if np.isfinite(v):
            pl.append(v)
    if not np.isfinite(r_obs) or len(pl) < 3:
        return None
    return {"n": len(common), "yrs": int(pd.Index(common).year.nunique()),
            "r": float(r_obs), "pl_max": float(np.max(pl)),
            "pl_med": float(np.median(pl)), "sync": bool(r_obs > np.max(pl))}


def summarize(pairs, label):
    ok = [p for p in pairs if p]
    if not ok:
        print(f"    {label}：**判定できる対が無い**"); return None
    frac = np.mean([p["sync"] for p in ok])
    print(f"    {label}：**{len(ok)} 組／揃う {sum(p['sync'] for p in ok)} 組"
          f"＝{frac:.0%}**")
    print(f"      r_obs 中央値 {np.median([p['r'] for p in ok]):+.3f}／"
          f"プラセボ中央値 {np.median([p['pl_med'] for p in ok]):+.3f}／"
          f"共通日 中央値 {int(np.median([p['n'] for p in ok]))}")
    return frac


def run(groups, tag=""):
    """`groups`：地点ラベル → {名前: 残差 Series}。**同一地点と地点間を分けて集計する。**"""
    same, cross = {}, []
    for gid, members in groups.items():
        ps = []
        for x, y in itertools.combinations(sorted(members), 2):
            p = pair_sync(members[x], members[y])
            if p:
                ps.append(p)
        if ps:
            same[gid] = ps
    keys = sorted(groups)
    for i, j in itertools.combinations(range(len(keys)), 2):
        for x in groups[keys[i]].values():
            for y in groups[keys[j]].values():
                p = pair_sync(x, y)
                if p:
                    cross.append(p)
    print(f"\n  ── 集計{tag} ──")
    fr_same = {}
    for gid, ps in sorted(same.items(), key=lambda kv: -len(kv[1])):
        f = summarize(ps, f"同一地点 [{gid}]")
        if f is not None:
            fr_same[gid] = (f, len(ps))
    fr_cross = summarize(cross, "地点間（全部）")
    return fr_same, fr_cross, cross


def load_real(cosore_dir):
    root = Path(cosore_dir)
    desc = pd.read_csv(root / "description.csv")
    recs = []
    for _, r in desc.iterrows():
        ds = str(r["CSR_DATASET"])
        f = root / "datasets" / f"data_{ds}.csv"
        if not f.exists():
            continue
        try:
            la, lo = float(r["CSR_LATITUDE"]), float(r["CSR_LONGITUDE"])
        except (TypeError, ValueError, KeyError):
            continue
        if not (np.isfinite(la) and np.isfinite(lo)):
            continue
        try:
            df, st, sm = load_cosore(f, None)
        except Exception:
            continue
        if df is None or "Tsoil" not in df or "Rs" not in df:
            continue
        cols = ["Rs", "Tsoil"] + (["SM"] if "SM" in df else [])
        daily = df[cols].groupby(df.index.normalize()).mean()
        if len(daily) < MIN_DAYS:
            continue
        res = chamber_residual(daily)
        if res is None:
            continue
        recs.append({"ds": ds, "lat": la, "lon": lo, "res": res})
    # 座標で地点をまとめる（単連結）
    n = len(recs); parent = list(range(n))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    for i, j in itertools.combinations(range(n), 2):
        if haversine(recs[i]["lat"], recs[i]["lon"], recs[j]["lat"], recs[j]["lon"]) <= SAME_KM:
            parent[find(i)] = find(j)
    groups, coords = {}, {}
    for i in range(n):
        g = find(i)
        groups.setdefault(g, {})[recs[i]["ds"]] = recs[i]["res"]
        coords[g] = (recs[i]["lat"], recs[i]["lon"])
    # 地点ラベルを読みやすくする
    named = {}
    for g, m in groups.items():
        lab = sorted(m)[0].split("_")[1] if "_" in sorted(m)[0] else str(g)
        named[f"{lab}({len(m)}本)"] = m
        coords[f"{lab}({len(m)}本)"] = coords[g]
    return named, coords


def dist_bins(cross_pairs, groups, coords):
    """**地点間を距離で 3 段に分ける**（事前登録どおり）。"""
    print("    ── 地点間を距離で分ける ──")
    keys = sorted(groups)
    bins = {"< 100 km": [], "100–1000 km": [], "> 1000 km": []}
    for i, j in itertools.combinations(range(len(keys)), 2):
        a, b = keys[i], keys[j]
        d = haversine(*coords[a], *coords[b])
        key = "< 100 km" if d < 100 else ("100–1000 km" if d < 1000 else "> 1000 km")
        for x in groups[a].values():
            for y in groups[b].values():
                p = pair_sync(x, y)
                if p:
                    bins[key].append(p)
    for k, v in bins.items():
        summarize(v, f"      {k}")


def synth(kind, n_site=3, n_ch=4, years=6, seed=0):
    """**地点共通／チャンバー独立／全地点共通**の未観測駆動を仕込む。"""
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2012-01-01", periods=365 * years, freq="D")
    doy = idx.dayofyear.to_numpy()
    glob = pd.Series(np.convolve(rng.normal(0, 1, len(idx)), np.ones(4) / 4, "same"), index=idx)
    groups, coords = {}, {}
    for s in range(n_site):
        site = pd.Series(np.convolve(rng.normal(0, 1, len(idx)), np.ones(4) / 4, "same"), index=idx)
        m = {}
        for c in range(n_ch):
            T = 15 + 10 * np.sin(2 * np.pi * (doy - 100) / 365) + rng.normal(0, 1.5, len(idx))
            W = np.clip(0.2 + 0.05 * np.sin(2 * np.pi * (doy - 200) / 365)
                        + rng.normal(0, .02, len(idx)), .02, .6)
            own = pd.Series(np.convolve(rng.normal(0, 1, len(idx)), np.ones(4) / 4, "same"), index=idx)
            hid = {"site": site, "global": glob}.get(kind, own)
            lrs = -1.0 + 0.06 * T + 2.0 * W + 0.25 * hid.to_numpy() + rng.normal(0, .05, len(idx))
            daily = pd.DataFrame({"Rs": np.exp(lrs), "Tsoil": T, "SM": W}, index=idx)
            r = chamber_residual(daily)
            if r is not None:
                m[f"S{s}C{c}"] = r
        groups[f"地点{s}({len(m)}本)"] = m
        coords[f"地点{s}({len(m)}本)"] = (40.0 + 10 * s, -100.0)
    return groups, coords


def verdict(fr_same, fr_cross):
    print("\n  === 結論（事前登録の判定規則に当てる）===")
    if not fr_same:
        print("  **判定しない**——同一地点で判定できる対が無い。"); return
    gid, (f, npair) = max(fr_same.items(), key=lambda kv: kv[1][1])
    print(f"  主たる地点：**{gid}**（{npair} 組・揃う割合 {f:.0%}）"
          f"／地点間 {fr_cross:.0%}" if fr_cross is not None else "")
    if npair < 10:
        print(f"  **判定しない**——同一地点の対が {npair} 組で 10 未満。"); return
    c = fr_cross if fr_cross is not None else 0.0
    if f >= HI and c <= LO:
        print("  **★地点スケール**——**未観測の駆動はチャンバーより広く、地点より狭い**")
        print("  （気象・土壌水文）。**局所の生物ではない。**")
    elif f <= LO:
        print("  **★チャンバー・スケール**——**微小生息場所の現象**。")
        print("  **「観測の隙間」は空間的にも小さい。**")
    elif f >= HI and c >= HI:
        print("  **★広域の気象**——**我々が測っていない気象量**。**最も大きな示唆。**")
    else:
        print("  **○中間**——**割合をそのまま記し、まとめない。**")
    print("\n  留保（事前登録どおり）：")
    print("   ・**同一地点の腕は n=1 地点**。**「28 例で確かめた」とは書かない。**")
    print("   ・**同一地点のチャンバーは同じ気象を受ける**——**テンソルビンで最大限つぶしたが、")
    print("     引き残しが共通しうる**ことは消えない。")
    print("   ・**プラセボを上回っても機構は分からない**。**空間スケールが絞れるだけ。**")


def main():
    ap = argparse.ArgumentParser(description="旗99：チャンバー間の残差の同期")
    ap.add_argument("--real", action="store_true")
    ap.add_argument("--cosore-dir", default="/mnt/hdd/cosore-0.7.0")
    a = ap.parse_args()

    tee_stdout("step99")
    print("=== 旗99：~4 日メモリは、どの空間スケールの現象か ===")
    print("  **残差は旗74 のテンソルビン＋外挿**（誤特定の共通残りを最大限つぶす）。")
    print(f"  **プラセボは ±30/60/90 日の多重シフト**——**6 通りの最大値を上回る**を「揃う」とする。")
    print(f"  **KAYE {HI:.0%} 以上 かつ 地点間 {LO:.0%} 以下 → 地点スケール**／"
          f"**KAYE {LO:.0%} 以下 → チャンバー・スケール**／**両方 {HI:.0%} 以上 → 広域の気象**。")

    if not a.real:
        print("\n  【合成データで道具を検証する】**実データに触れる前に**、")
        print("  **`chamber`（独立な駆動）で「揃う」と出ないか**を必ず見る。")
        print("  **出たら駆動の誤特定が残差に共通して残っている＝実データに進まない。**")
        for kind, want in (("site", "**同一地点で揃い、地点間では揃わない**べき"),
                           ("chamber", "**同一地点でも揃わない**べき"),
                           ("global", "**地点間でも揃う**べき")):
            print(f"\n  ===== 合成 `{kind}` —— 期待：{want} =====")
            g, c = synth(kind)
            fs, fc, _ = run(g, f"（{kind}）")
            print(f"  【判定】同一地点 {[f'{k}:{v[0]:.0%}' for k, v in fs.items()]}"
                  f"／地点間 {fc:.0%}" if fc is not None else "")
        print("\n  → **site→同一高・地点間低／chamber→同一も低／global→両方高** なら道具は使える。")
        return

    groups, coords = load_real(a.cosore_dir)
    print(f"\n  残差を作れたチャンバー：{sum(len(m) for m in groups.values())} 本"
          f"／地点 {len(groups)}")
    fs, fc, cross = run(groups)
    dist_bins(cross, groups, coords)
    verdict(fs, fc)


if __name__ == "__main__":
    main()
