"""多サイト横断：高次情報構造（相乗 vs 冗長）が生態系をまたいで普遍かを見る。

本体 `oinfo_analysis.scan_years` を再利用し、各サブシステムの O-information z を
複数生態系で並べる。JP-Tak で見つけた「呼吸系だけ相乗支配・他は冗長支配」が、
生態系を越える**普遍の高次署名**か、その生態系固有かを検定する（旗13 の高次版）。

各サブ系 × サイトで、健全全年の z を集計し
  ・mean z（＜0 相乗寄り／＞0 冗長寄り）
  ・相乗年/冗長年（|z|≥2.36 を有意として）
を出す。全サイトで符号が一致すれば**普遍の高次構造**、飛べば生態系固有。

注（旗13 との検出力の違い）：O-info の z は各年 ~3000 点のサロゲート比なので**各年が
well-powered**。だから年数が少ないサイトでも「その年の相乗/冗長」は信頼できる。年数は
「符号の安定性をどれだけ確信できるか」に効く（少ない＝安定性の確証が弱い）だけ。

    python -m japanflux_pn.oinfo_compare --heatmap oinfo_cross.png --out-csv oinfo_cross.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from .config import AnalysisConfig
from .oinfo_analysis import scan_years, SUBSYSTEMS, OINFO_BINS_DEFAULT


DEFAULT_SITES = ["JP-Tak", "JP-Ta2", "JP-BBY", "CN-HaM", "MN-Kbu"]
SIG = 2.36


def site_summary(site: str, obins: int, config: AnalysisConfig) -> pd.DataFrame:
    """1 サイトの各サブ系: mean z・相乗年・冗長年・年数。"""
    long = scan_years(site, obins, config)
    rows = []
    if long.empty:
        return pd.DataFrame(rows)
    for name, g in long.groupby("subsystem", sort=False):
        z = g["z"].to_numpy(dtype=float)
        rows.append({"subsystem": name, "mean_z": float(np.nanmean(z)),
                     "n_syn": int((z <= -SIG).sum()), "n_red": int((z >= SIG).sum()),
                     "n_years": int(len(g))})
    return pd.DataFrame(rows)


def classify(mean_z: float, n_syn: int, n_red: int, n_years: int) -> str:
    """そのサイト・サブ系の支配傾向（60% 規則、本体 report_years と同じ）。"""
    if n_years == 0:
        return "—"
    if n_syn >= 0.6 * n_years:
        return "相乗"
    if n_red >= 0.6 * n_years:
        return "冗長"
    return "混在"


def collect(sites: list[str], obins: int, config: AnalysisConfig
            ) -> dict[str, pd.DataFrame]:
    per = {}
    for s in sites:
        print(f"\n----- {s} O-info 収集中 -----", flush=True)
        try:
            df = site_summary(s, obins, config)
        except Exception as e:
            print(f"  {s}: SKIP {type(e).__name__}: {e}", flush=True)
            continue
        if df.empty:
            print(f"  {s}: 有効年なし → 除外", flush=True)
            continue
        per[s] = df
    return per


def draw_heatmap(zmat: pd.DataFrame, path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
    jp_path = "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"
    jp = fm.FontProperties(fname=jp_path) if Path(jp_path).exists() else None
    m = zmat
    vmax = float(np.nanmax(np.abs(m.to_numpy(dtype=float)))) or 1.0
    fig, ax = plt.subplots(figsize=(1.8 + 1.0 * m.shape[1], 0.7 * m.shape[0] + 1.6))
    im = ax.imshow(m.to_numpy(dtype=float), cmap="RdBu", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=20, ha="right")
    ax.set_yticks(range(m.shape[0])); ax.set_yticklabels(m.index, fontproperties=jp)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            v = m.iat[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:+.0f}", ha="center", va="center", fontsize=9,
                        color="white" if abs(v) >= 0.55 * vmax else "#222")
    ax.set_title("高次情報構造の生態系間比較（O-info 平均z　青=冗長 / 赤=相乗）",
                 fontproperties=jp)
    cb = fig.colorbar(im, ax=ax, fraction=0.03)
    cb.set_label("平均 z（＞0 冗長 / ＜0 相乗）", fontproperties=jp)
    fig.savefig(path, bbox_inches="tight", dpi=130); plt.close(fig)


def report(sites: list[str], obins: int = OINFO_BINS_DEFAULT,
           config: AnalysisConfig | None = None, heatmap=None, out_csv=None,
           n_min: int = 10) -> pd.DataFrame:
    config = config or AnalysisConfig()
    per = collect(sites, obins, config)
    if len(per) < 2:
        print(f"\n有効サイト {len(per)} < 2。比較できません。")
        return pd.DataFrame()
    used = list(per.keys())
    subs = list(SUBSYSTEMS.keys())

    # 平均 z 表（行=サブ系, 列=サイト）と分類表
    zmat = pd.DataFrame(index=subs, columns=used, dtype=float)
    cls = pd.DataFrame(index=subs, columns=used, dtype=object)
    nyr = {}
    for s in used:
        d = per[s].set_index("subsystem")
        nyr[s] = int(d["n_years"].max()) if len(d) else 0
        for sub in subs:
            if sub in d.index:
                r = d.loc[sub]
                zmat.loc[sub, s] = r["mean_z"]
                cls.loc[sub, s] = classify(r["mean_z"], r["n_syn"], r["n_red"], r["n_years"])

    print(f"\n===== 高次情報構造の生態系間比較（{len(used)} サイト）=====")
    print("  年数: " + "  ".join(f"{s}:{nyr[s]}年" for s in used))
    print("\n  各サブ系の平均 z（＜0 相乗 / ＞0 冗長）と支配傾向\n")
    print(f"  {'サブシステム':<22}" + "".join(f"{s:>12}" for s in used))
    for sub in subs:
        cells = ""
        for s in used:
            v = zmat.loc[sub, s]; c = cls.loc[sub, s]
            cells += f"{v:>+8.1f}{('('+str(c)+')'):>5}" if np.isfinite(v) else f"{'--':>12}"
        print(f"  {sub:<22}{cells}")

    # 普遍性判定：符号（相乗/冗長）が全サイトで一致するか
    print("\n=== 高次署名は生態系をまたいで普遍か（相乗/冗長の符号一致）===")
    for sub in subs:
        signs = [cls.loc[sub, s] for s in used if cls.loc[sub, s] in ("相乗", "冗長")]
        mixed = [s for s in used if cls.loc[sub, s] == "混在"]
        if not signs:
            print(f"  {sub:<22} 判定不能（全サイト混在/欠）"); continue
        syn = signs.count("相乗"); red = signs.count("冗長")
        if syn == len(signs) and not mixed:
            verd = "✅ 全生態系で相乗＝普遍の創発署名"
        elif red == len(signs) and not mixed:
            verd = "✅ 全生態系で冗長＝普遍の共通駆動署名"
        elif syn == 0 or red == 0:
            verd = f"○ 符号一致（相乗{syn}/冗長{red}）だが混在サイトあり: {','.join(mixed)}"
        else:
            verd = f"× 符号が生態系で飛ぶ＝固有（相乗{syn}/冗長{red}, 混在{len(mixed)}）"
        print(f"  {sub:<22} {verd}")

    print("\n  読み方: 全サイトで同符号＝その高次署名（呼吸=相乗 等）は生態系を越える普遍構造。")
    print("         符号が飛ぶ＝その生態系固有。年数が少ないサイトは符号は信頼できるが安定性の確証は弱い。")

    # 検出力版：年数の多いサイト(n_years≥n_min)の平均zの符号で普遍性を採る。
    # 上の分類(60%規則)は年数が少ないと"混在"に落ちやすく、また n=5 サイトの反転に
    # 引っ張られる。旗13 の教訓と同じく、確信を持てるサイトだけで符号を採り直す。
    print(f"\n=== 検出力版の普遍性（年数 n≥{n_min} のサイトの平均z符号で採る）===")
    for sub in subs:
        pw = [(s, zmat.loc[sub, s]) for s in used
              if nyr[s] >= n_min and np.isfinite(zmat.loc[sub, s])]
        if len(pw) < 2:
            print(f"  {sub:<22} 検出力サイト<2 → 保留"); continue
        neg = [s for s, v in pw if v < 0]   # 相乗寄り
        pos = [s for s, v in pw if v > 0]   # 冗長寄り
        detail = " ".join(f"{s}:{v:+.1f}" for s, v in pw)
        if not pos:
            verd = f"✅ 検出力{len(pw)}/{len(pw)} 相乗寄り＝普遍の創発署名"
        elif not neg:
            verd = f"✅ 検出力{len(pw)}/{len(pw)} 冗長寄り＝普遍の共通駆動署名"
        elif len(neg) >= len(pos):
            verd = f"○ 相乗寄り優勢（相乗{len(neg)}/冗長{len(pos)}）: 例外 {','.join(pos)}"
        else:
            verd = f"○ 冗長寄り優勢（相乗{len(neg)}/冗長{len(pos)}）: 例外 {','.join(neg)}"
        print(f"  {sub:<22} {verd}   [{detail}]")
    print(f"\n  ※O-info の z は各年~3000点のサロゲート比＝各年 well-powered。年数は『符号の"
          f"安定性』の確証に効くだけ（旗13 の年数=検出力とは別の意味）。")

    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        zmat.to_csv(out_csv)
        print(f"\n[saved] {out_csv}")
    if heatmap:
        draw_heatmap(zmat.astype(float), heatmap)
        print(f"[saved] {heatmap}")
    return zmat


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="高次情報構造（相乗/冗長）の生態系間比較")
    p.add_argument("--sites", nargs="+", default=DEFAULT_SITES)
    p.add_argument("--obins", type=int, default=OINFO_BINS_DEFAULT)
    p.add_argument("--heatmap", default=None)
    p.add_argument("--out-csv", default=None)
    p.add_argument("--n-min", type=int, default=10,
                   help="検出力版の普遍性で符号を採る最小年数（既定10）")
    a = p.parse_args(argv)
    report(a.sites, a.obins, heatmap=a.heatmap, out_csv=a.out_csv, n_min=a.n_min)


if __name__ == "__main__":
    main()
