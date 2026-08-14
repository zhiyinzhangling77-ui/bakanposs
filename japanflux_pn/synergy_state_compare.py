"""旗15 の確認：高次構造の状態依存を"多生態系"で。JP-Tak で探索的に見えた
「乾いた年ほど炭素/光合成の冗長が弱まる（r(z,VPD)<0）」を、独立な他生態系でも
同じ向きに出るかで確認する（p-hack でなくクロスサイト）。

本体 `synergy_state.build`（年ごと各サブ系 z ＋ 生 VPD/土壌水分/気温）を各サイトで作り、
サブ系ごとに r(z, VPD) をサイト横断で並べる。判定は旗13 と同じ検出力ロジック：
r は各サイトの"年数"に効くので、n≥n_min 年のサイトだけで符号一致を採り、少年数は保留。
さらに層別併合検定（各サイト内で VPD を並べ替え）で共通の向きの p 値を 1 本出す。

事前指定の確認対象（JP-Tak で見えた順）：炭素コア・光合成の冗長 vs VPD（期待 r<0）。
呼吸は JP-Tak でヌルだったので「ヌルの再現」を確認。

    python -m japanflux_pn.synergy_state_compare --heatmap synergy_state_cross.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig
from .oinfo_analysis import SUBSYSTEMS, OINFO_BINS_DEFAULT
from .synergy_state import build
from .climate_response import _spearman, _spearman_p

DEFAULT_SITES = ["JP-Tak", "JP-Ta2", "JP-BBY", "CN-HaM", "MN-Kbu"]
STATE = "VPD_mean"
# 事前指定の確認対象（サブ系プレフィックス, 期待符号）
PRESPEC = [("炭素コア", "-"), ("光合成", "-"), ("呼吸", "0")]


def collect(sites: list[str], obins: int, config: AnalysisConfig
            ) -> dict[str, pd.DataFrame]:
    per = {}
    for s in sites:
        print(f"\n----- {s} 収集中 -----", flush=True)
        try:
            df = build(s, obins, config)
        except Exception as e:
            print(f"  {s}: SKIP {type(e).__name__}: {e}", flush=True)
            continue
        if len(df) < 4 or STATE not in df.columns:
            print(f"  {s}: 有効年 {len(df)} < 4 → 除外", flush=True)
            continue
        per[s] = df
        print(f"  {s}: 有効年 {len(df)}", flush=True)
    return per


def _sub_col(df: pd.DataFrame, prefix: str) -> str | None:
    for c in df.columns:
        if str(c).startswith(prefix):
            return c
    return None


def pooled_perm(per: dict[str, pd.DataFrame], prefix: str, n_min: int,
                n_perm: int = 20000, seed: int = 0) -> dict:
    """検出力サイト(n≥n_min年)を層別に束ね、r(z,VPD)の共通の向きの p を出す。"""
    xs, ys, rs, used = [], [], [], []
    for s, df in per.items():
        col = _sub_col(df, prefix)
        if col is None or len(df) < n_min:
            continue
        z = df[col].to_numpy(dtype=float); v = df[STATE].to_numpy(dtype=float)
        r = _spearman(z, v)
        if np.isfinite(r):
            xs.append(z); ys.append(v); rs.append(r); used.append(s)
    k = len(used)
    if k < 2:
        return {"k": k, "used": used, "per_r": dict(zip(used, rs)), "mean_r": np.nan,
                "p_two": np.nan}
    mean_r = float(np.mean(rs))
    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(n_perm):
        acc = sum(_spearman(x, rng.permutation(y)) for x, y in zip(xs, ys))
        if abs(acc / k) >= abs(mean_r) - 1e-12:
            ge += 1
    return {"k": k, "used": used, "per_r": dict(zip(used, rs)),
            "mean_r": mean_r, "p_two": (ge + 1) / (n_perm + 1)}


def draw_heatmap(rmat: pd.DataFrame, path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    jp_path = "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"
    jp = fm.FontProperties(fname=jp_path) if Path(jp_path).exists() else None
    m = rmat
    fig, ax = plt.subplots(figsize=(1.8 + 1.0 * m.shape[1], 0.7 * m.shape[0] + 1.6))
    im = ax.imshow(m.to_numpy(dtype=float), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=20, ha="right")
    ax.set_yticks(range(m.shape[0])); ax.set_yticklabels(m.index, fontproperties=jp)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = m.iat[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=9,
                        color="white" if abs(v) >= 0.55 else "#222")
    ax.set_title("高次構造の状態依存の生態系間比較（r( z , VPD )）", fontproperties=jp)
    cb = fig.colorbar(im, ax=ax, fraction=0.03); cb.set_label("Spearman r", fontproperties=jp)
    fig.savefig(path, bbox_inches="tight", dpi=130); plt.close(fig)


def report(sites: list[str], obins: int = OINFO_BINS_DEFAULT,
           config: AnalysisConfig | None = None, heatmap=None, out_csv=None,
           n_min: int = 10, n_perm: int = 20000) -> pd.DataFrame:
    config = config or AnalysisConfig()
    per = collect(sites, obins, config)
    if len(per) < 2:
        print(f"\n有効サイト {len(per)} < 2。比較できません。")
        return pd.DataFrame()
    used = list(per.keys())
    nyr = {s: len(per[s]) for s in used}
    subs = list(SUBSYSTEMS.keys())

    rmat = pd.DataFrame(index=subs, columns=used, dtype=float)
    for s in used:
        df = per[s]
        for sub in subs:
            col = _sub_col(df, sub.split()[0])
            if col is not None:
                rmat.loc[sub, s] = _spearman(df[col].to_numpy(dtype=float),
                                             df[STATE].to_numpy(dtype=float))

    print(f"\n===== 高次構造の状態依存 r(z, VPD) の生態系間比較（{len(used)} サイト）=====")
    print("  年数: " + "  ".join(f"{s}:{nyr[s]}" for s in used))
    print("  z<0=相乗/z>0=冗長。冗長サブ系で r<0＝乾いた年ほど冗長が弱まる（旗13の高次版）\n")
    print(f"  {'サブシステム':<22}" + "".join(f"{s:>10}" for s in used))
    for sub in subs:
        cells = "".join(f"{rmat.loc[sub, s]:>+10.2f}" if np.isfinite(rmat.loc[sub, s])
                        else f"{'--':>10}" for s in used)
        print(f"  {sub:<22}{cells}")

    print(f"\n=== 事前指定の確認（検出力サイト n≥{n_min}年 の符号一致＋層別併合 p, {n_perm}回）===")
    print(f"  参考: p<0.05 に要る |r| ≈ n=5→0.88, n=10→0.63, n=11→0.60, n=21→0.43")
    for prefix, sign in PRESPEC:
        pt = pooled_perm(per, prefix, n_min, n_perm)
        name = next((s for s in subs if s.startswith(prefix)), prefix)
        if pt["k"] < 2:
            print(f"\n  ■ {name}: 検出力サイト<2 → 併合不能"); continue
        per_r = "  ".join(f"{s}:{pt['per_r'][s]:+.2f}" for s in pt["used"])
        # 各検出力サイトの符号一致（期待符号）
        if sign == "0":
            hits = sum(1 for s in pt["used"] if abs(pt["per_r"][s]) < 0.4)
            agree_txt = f"ヌル一致(|r|<0.4) {hits}/{pt['k']}"
        else:
            want_neg = (sign == "-")
            hits = sum(1 for s in pt["used"]
                       if (pt["per_r"][s] < 0) == want_neg)
            agree_txt = f"符号一致 {hits}/{pt['k']}"
        print(f"\n  ■ {name}（期待 {sign}, 検出力{pt['k']}生態系: {per_r}）")
        print(f"    {agree_txt}   併合 平均r={pt['mean_r']:+.3f}  両側p={pt['p_two']:.4f}")
        if sign == "0":
            verd = ("✅ ヌルを再現（状態非依存＝構造的）" if pt["p_two"] >= 0.05
                    else "△ 併合では有意（単一サイトのヌルと不一致・要精査）")
        else:
            sig = pt["p_two"] < 0.05 and ((pt["mean_r"] < 0) == (sign == "-"))
            verd = ("✅ 併合で有意＝乾燥で冗長が弱まるが生態系をまたいで確か（旗13の高次版を確認）"
                    if sig else "△ 併合で非有意（探索的トレンドは多生態系では確証に届かず）")
        print(f"    → {verd}")

    print("\n  ※これは旗15の探索的トレンドの『事前指定』確認。単一サイトの事後選択でなく、")
    print("    独立サイトで同符号かつ併合有意なら初めて主張できる。")

    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        rmat.to_csv(out_csv); print(f"\n[saved] {out_csv}")
    if heatmap:
        draw_heatmap(rmat.astype(float), heatmap); print(f"[saved] {heatmap}")
    return rmat


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="高次構造の状態依存の生態系間比較")
    p.add_argument("--sites", nargs="+", default=DEFAULT_SITES)
    p.add_argument("--obins", type=int, default=OINFO_BINS_DEFAULT)
    p.add_argument("--n-min", type=int, default=10)
    p.add_argument("--n-perm", type=int, default=20000)
    p.add_argument("--heatmap", default=None)
    p.add_argument("--out-csv", default=None)
    a = p.parse_args(argv)
    report(a.sites, a.obins, heatmap=a.heatmap, out_csv=a.out_csv,
           n_min=a.n_min, n_perm=a.n_perm)


if __name__ == "__main__":
    main()
