"""発表用の図を、これまでの解析結果の数値から生成する（実データ不要・数値埋め込み）。

英語ラベル（学会標準）。出力は japanflux_pn/slides/*.png。説明は SLIDES ガイド参照。
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import networkx as nx

# 日本語フォント（IPAGothic）。ラテン文字・数字・矢印もこのフォントで統一。
_JP_PATH = "/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf"
fm.fontManager.addfont(_JP_PATH)
JP = fm.FontProperties(fname=_JP_PATH)

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 130, "font.size": 13,
    "axes.titlesize": 15, "axes.labelsize": 13, "axes.spines.top": False,
    "axes.spines.right": False, "figure.autolayout": True,
    "font.family": JP.get_name(),
})
OUT = Path("japanflux_pn/slides")
OUT.mkdir(parents=True, exist_ok=True)

BLUE, RED, GREEN, ORANGE, GREY = "#1f6fb2", "#c0392b", "#2e8b57", "#e08a1e", "#888"

# 3図で共通の意味づけ配色（統一パレット）
C_RADIATION = "#f0b73e"   # 放射（共通原因そのもの）
C_APPARENT  = BLUE        # 見かけのつながり（共通原因の影／条件付けで崩れる）
C_APP_LIGHT = "#bcd4e8"   # 同・条件付け後（薄色）
C_REAL      = GREEN       # 直接・頑健に残る因果
C_REAL_LIGHT = "#b8e0b8"  # 同・薄色
C_ARTIFACT  = RED         # 物理的にありえない（アーティファクト）

# 変数の日本語表示名
JVAR = {
    "Rg": "日射", "Ta": "気温", "VPD": "飽差", "Ts": "地表面温度", "P": "降水",
    "th": "土壌水分", "θ": "土壌水分", "gH": "顕熱", "gLE": "潜熱",
    "GER": "呼吸", "NEE": "正味CO2", "GEP": "光合成",
}


# ---------------------------------------------------------------------------
# Fig 1 (star): light-use coupling I(Rg;GPP) vs atmospheric dryness (JP-Tak, 21 yr)
# ---------------------------------------------------------------------------
def fig_climate():
    yrs = [1999,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,
           2012,2013,2014,2015,2016,2017,2018,2019,2020,2021]
    vpd = [3.03,4.34,3.63,3.61,1.97,4.72,2.83,2.55,2.92,3.24,1.56,
           2.52,3.91,2.53,3.22,2.63,2.31,4.26,2.53,2.79,2.34]
    igp = [19.7,16.3,17.5,21.1,27.9,20.8,25.2,27.3,24.1,21.6,21.9,
           22.9,20.6,26.6,19.3,21.6,24.8,14.6,24.3,24.0,20.3]
    vpd, igp = np.array(vpd), np.array(igp)
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.scatter(vpd, igp, s=55, color=BLUE, zorder=3, edgecolor="white")
    m, b = np.polyfit(vpd, igp, 1)
    xs = np.linspace(vpd.min(), vpd.max(), 50)
    ax.plot(xs, m*xs+b, color=RED, lw=2, zorder=2)
    for yr, xx, yy, dx, dy, note in [
        (2003, 1.97, 27.9, 0.15, 0.6, "2003\ncool summer"),
        (2018, 4.26, 14.6, -0.05, 1.2, "2018\nheatwave")]:
        ax.scatter([xx],[yy], s=120, facecolor="none", edgecolor=RED, lw=2, zorder=4)
        ax.annotate(note, (xx, yy), xytext=(xx+dx, yy+dy), fontsize=11,
                    color=RED, ha="center", fontweight="bold")
    ax.set_xlabel("Atmospheric dryness  VPD  (hPa, summer mean)")
    ax.set_ylabel("Light-use coupling  I(Rg ; GPP)  (%)")
    ax.set_title("Drier summers → photosynthesis decouples from light\n"
                 "JP-Tak forest, 21 years   (Spearman r = -0.68, p = 0.001)")
    fig.savefig(OUT/"fig1_climate_response.png"); plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 2: common driver — PID redundancy fraction (flux = shared Rg info, thermal = direct)
# ---------------------------------------------------------------------------
def fig_pid():
    red = [("gH<-GEP",100),("NEE<-gLE",100),("gH<-NEE",99),
           ("gLE<-GEP",97),("VPD<-gLE",86)]
    uniq = [("Ts<-Ta",9),("GER<-Ta",10),("GER<-Ts",13),("th<-GER",12)]
    labels = [l for l,_ in red] + [l for l,_ in uniq]
    vals = [v for _,v in red] + [v for _,v in uniq]
    cols = [BLUE]*len(red) + [ORANGE]*len(uniq)
    y = np.arange(len(labels))[::-1]
    fig, ax = plt.subplots(figsize=(7.6, 5.0))
    ax.barh(y, vals, color=cols)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.axvline(60, color=GREY, ls="--", lw=1)
    ax.set_xlabel("Redundancy fraction  R / I(Y;X)  (%)  — share explained by radiation Rg")
    ax.set_title("Flux/carbon coupling is shared radiation info (redundant);\n"
                 "temperature–respiration coupling is direct (unique).  JP-Tak")
    ax.set_xlim(0, 105)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=BLUE, label="Flux/carbon: redundant with Rg"),
                       Patch(color=ORANGE, label="Thermal: unique (direct)")],
              loc="lower right", fontsize=11, frameon=False)
    fig.savefig(OUT/"fig2_pid_redundancy.png"); plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 3: PCMCI causal skeleton (JP-Tak forest)
# ---------------------------------------------------------------------------
def fig_skeleton():
    # 円は小さめ（0.6倍）、文字は大きく。配置は少し詰める。
    pos = {"Rg":(0,2.1), "VPD":(-2.0,1.0), "Ta":(0.35,1.0), "gLE":(-2.2,-0.65),
           "gH":(-0.5,-0.8), "Ts":(1.55,0.05), "GEP":(2.35,1.4), "NEE":(3.0,-0.15)}
    edges = [("Rg","VPD"),("Rg","Ta"),("Rg","gLE"),("Rg","gH"),("Rg","Ts"),
             ("Ta","Ts"),("VPD","gLE"),("GEP","NEE")]
    # 日本語（正式名称）＋アルファベット略号の2段表示
    labels = {"Rg":"日射\nRg", "VPD":"飽差\nVPD", "Ta":"気温\nTa",
              "gLE":"潜熱\nγLE", "gH":"顕熱\nγH", "Ts":"地表面温度\nTs",
              "GEP":"光合成\nGEP", "NEE":"正味CO2\nNEE"}
    G = nx.DiGraph(); G.add_nodes_from(pos); G.add_edges_from(edges)
    fig, ax = plt.subplots(figsize=(14.0, 9.5))
    # 統一配色（青＋緑のみ）：日射＝共通原因(濃い青)、光合成・正味CO2＝直接残る炭素(緑)、その他気象＝薄青
    ncol = ["#5b9bd5" if n=="Rg" else (C_REAL_LIGHT if n in ("GEP","NEE")
            else C_APP_LIGHT) for n in G.nodes]
    NS = 15600  # 円は 0.6 倍（26000→15600）
    nx.draw_networkx_nodes(G, pos, node_size=NS, node_color=ncol,
                           edgecolors="#555", linewidths=2, ax=ax)
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=26,
                            font_family=JP.get_name(), ax=ax)
    nx.draw_networkx_edges(G, pos, arrowstyle="-|>", arrowsize=26,
                           width=2.5, edge_color=BLUE,
                           node_size=NS, ax=ax)
    ax.set_title("共通原因（日射）を差し引いて残る因果の骨組み（PCMCI）\n"
                 "日射が気象・エネルギーを動かし、光合成が正味炭素を動かす",
                 fontsize=18)
    ax.set_xlim(-3.2, 3.9); ax.set_ylim(-1.6, 2.9)
    ax.axis("off")
    fig.savefig(OUT/"fig3_causal_skeleton.png", bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 4: multi-year robustness — link consistency (core vs artifact)
# ---------------------------------------------------------------------------
def fig_robustness():
    def jp(link):  # "Rg->gH" -> "日射→顕熱"
        a, b = link.split("->")
        return f"{JVAR[a]}→{JVAR[b]}"
    core = [("GEP->NEE",100),("Rg->gH",95),("Rg->gLE",95),("Rg->VPD",95),
            ("Ta->Ts",90),("Rg->Ta",86),("Rg->Ts",76),("Ta->GER",71)]
    art = [("gH->Rg",57),("gLE->Rg",38),("GEP->Rg",14),("NEE->Rg",14)]
    labels = [jp(l) for l,_ in core] + [jp(l) for l,_ in art]
    vals = [v for _,v in core] + [v for _,v in art]
    # 統一配色（青＋緑のみ）：緑＝毎年安定して残るリンク、青＝物理的にありえない（→日射）
    cols = [C_REAL]*len(core) + [C_APPARENT]*len(art)
    y = np.arange(len(labels))[::-1]
    n = len(labels)
    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    ax.barh(y, vals, color=cols)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.set_xlabel("そのリンクが出現した年の割合（％）  ― JP-Tak・21年")
    ax.set_title("同じ因果リンクは年をまたいで安定して現れる。\n"
                 "→日射のリンクは頻度によらず物理的にありえない（＝誤検出）",
                 fontsize=13.5)
    ax.set_xlim(0, 105); ax.set_ylim(-0.7, n-0.3)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=C_REAL, label="毎年安定して残るリンク（コア）"),
                       Patch(color=C_APPARENT, label="→日射（物理的にありえない）＝誤検出")],
              loc="lower right", fontsize=11, frameon=False)
    fig.savefig(OUT/"fig4_robustness.png", bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 5 (novelty): higher-order synergy of soil–respiration system, by ecosystem
# ---------------------------------------------------------------------------
def fig_oinfo():
    eco = ["Forest\n(JP-Tak)", "Rice paddy\n(JP-Mse)", "Bog\n(JP-BBY)"]
    meanz = [-7.1, 2.7, 2.3]
    frac = ["19/21 yr", "0/8 yr", "1/5 yr"]
    cols = [GREEN, RED, GREY]
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    bars = ax.bar(eco, meanz, color=cols, width=0.6)
    ax.axhline(0, color="k", lw=0.8)
    ax.axhspan(-2.36, 2.36, color="#eee", zorder=0)
    ax.text(2.35, 0, "not significant", va="center", ha="right",
            fontsize=9, color=GREY)
    for bar, f in zip(bars, frac):
        yy = bar.get_height()
        ax.text(bar.get_x()+bar.get_width()/2, yy + (0.6 if yy>0 else -0.9),
                f, ha="center", fontsize=11, fontweight="bold")
    ax.set_ylabel("Higher-order coupling  (mean z of O-information)\n"
                  "<0 = synergy (emergent),  >0 = redundancy")
    ax.set_title("Soil–respiration system {Rg,Ta,θ,GER}:\n"
                 "synergy in natural forest, but collapses in managed paddy")
    ax.set_ylim(-11, 6)
    fig.savefig(OUT/"fig5_oinfo_synergy.png"); plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 2b (clearer): coupling BEFORE vs AFTER removing the common driver (Rg)
# ---------------------------------------------------------------------------
def fig_conditioning():
    # (pair, I(X;Y) before, I(X;Y|Rg) after)  — JP-Tak
    flux = [("顕熱–光合成",10.9,1.4),("潜熱–光合成",11.2,2.3),
            ("顕熱–潜熱",8.7,1.1),("正味CO2–光合成",21.9,6.9)]
    therm = [("気温–地表面温度",19.1,21.8),("気温–呼吸",17.5,20.7),
             ("地表面温度–呼吸",13.6,18.7),("土壌水分–呼吸",7.8,12.6)]
    labels = [l for l,_,_ in flux] + [l for l,_,_ in therm]
    before = [b for _,b,_ in flux] + [b for _,b,_ in therm]
    after  = [a for _,_,a in flux] + [a for _,_,a in therm]
    x = np.array([0,1,2,3, 4.8,5.8,6.8,7.8])
    w = 0.38
    FS = 26  # 目盛り・軸見出しの基準サイズ（従来比およそ2倍）
    fig, ax = plt.subplots(figsize=(15.5, 10.5))
    # 統一配色：左＝見かけ(青)、右＝直接残る(緑)。濃色＝条件付け前、薄色＝日射で条件付け後。
    cols_before = [C_APPARENT]*4 + [C_REAL]*4
    cols_after  = [C_APP_LIGHT]*4 + [C_REAL_LIGHT]*4
    ax.bar(x-w/2, before, w, color=cols_before)
    ax.bar(x+w/2, after,  w, color=cols_after, edgecolor="#7f7f7f")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=FS)
    ax.tick_params(axis="y", labelsize=FS)
    ax.set_ylabel("つながりの強さ  I ＝ 相互情報量 / log(11)  （％）", fontsize=FS)
    ax.set_ylim(0, 34)
    ax.axvspan(-0.6, 3.6, color="#eef4fa", zorder=0)
    ax.axvspan(4.2, 8.4, color="#eef7ef", zorder=0)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=GREY, label="日射で条件付け前  I(X;Y)"),
                       Patch(facecolor="#dddddd", edgecolor="#7f7f7f",
                             label="日射で条件付け後  I(X;Y | 日射)")],
              loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2,
              fontsize=FS-3, frameon=False)
    ax.text(1.5, 31.3, "エネルギー・炭素フラックス", ha="center",
            fontsize=30, color=C_APPARENT, fontweight="bold")
    ax.text(6.3, 31.3, "気温・土壌・呼吸", ha="center",
            fontsize=30, color=C_REAL, fontweight="bold")
    for xi, b, a in zip(x[:4], before[:4], after[:4]):
        ax.annotate("", xy=(xi+w/2, a+0.6), xytext=(xi-w/2, b-0.6),
                    arrowprops=dict(arrowstyle="->", color="#555", lw=1.8))
    ax.text(1.3, 15.5, "日射で説明できる分を\n除くとほぼ消える", color=C_APPARENT,
            ha="center", fontsize=27, fontweight="bold")
    fig.savefig(OUT/"fig2b_conditioning.png", bbox_inches="tight"); plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 5b (honest): fraction of years synergistic + Wilson 95% CI (shows sample size)
# ---------------------------------------------------------------------------
def _wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0, 1.0
    p = k / n
    c = (k + z*z/2) / (n + z*z)
    h = z/(n+z*z) * np.sqrt(k*(n-k)/n + z*z/4)
    return p, max(0, c-h), min(1, c+h)


def fig_oinfo_ci():
    data = [("Forest\n(JP-Tak)", 19, 21, GREEN),
            ("Rice paddy\n(JP-Mse)", 0, 8, RED),
            ("Bog\n(JP-BBY)", 1, 5, GREY)]
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for i, (name, k, n, col) in enumerate(data):
        p, lo, hi = _wilson(k, n)
        ax.bar(i, p*100, color=col, width=0.6, zorder=2)
        ax.errorbar(i, p*100, yerr=[[ (p-lo)*100 ],[ (hi-p)*100 ]], fmt="none",
                    ecolor="#333", capsize=8, lw=1.8, zorder=3)
        ax.text(i, hi*100+3, f"{k}/{n} yr", ha="center", fontsize=12,
                fontweight="bold")
    ax.set_xticks(range(3)); ax.set_xticklabels([d[0] for d in data])
    ax.set_ylabel("Years showing higher-order synergy  (%)\n"
                  "in the soil–respiration system")
    ax.set_ylim(0, 112)
    ax.set_title("Synergy: consistent in forest, absent in paddy\n"
                 "(whiskers = 95% CI; wide when few years)", fontsize=14)
    ax.axhline(50, color=GREY, ls=":", lw=1)
    fig.savefig(OUT/"fig5b_oinfo_ci.png"); plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 5c (honest, 2 forests): belowground synergy in two subsystems x four sites
# ---------------------------------------------------------------------------
def fig_oinfo_twosub():
    sites = ["Deciduous\nforest\n(JP-Tak)", "Evergreen\nforest\n(JP-Ta2)",
             "Rice paddy\n(JP-Mse)", "Bog\n(JP-BBY)"]
    # fraction of years synergistic (%)
    resp = [19/21*100, 6/11*100, 0/8*100, 1/5*100]      # {Rg,Ta,θ,GER}
    soil = [11/21*100, 10/11*100, 0/8*100, 1/5*100]     # {Ta,Ts,θ,GER}
    x = np.arange(4); w = 0.38
    fig, ax = plt.subplots(figsize=(8.4, 5.4))
    ax.bar(x-w/2, resp, w, label="Respiration control  {Rg,Ta,θ,GER}", color=BLUE)
    ax.bar(x+w/2, soil, w, label="Soil–thermal  {Ta,Ts,θ,GER}", color=ORANGE)
    ax.axhline(50, color=GREY, ls=":", lw=1)
    ax.axvspan(-0.5, 1.5, color="#eaf5ea", zorder=0)   # forests
    ax.axvspan(1.5, 3.5, color="#f7ecec", zorder=0)    # managed / bog
    ax.text(0.5, 104, "Natural forests", ha="center", color=GREEN,
            fontsize=12, fontweight="bold")
    ax.text(2.5, 104, "Managed paddy / bog", ha="center", color=RED,
            fontsize=12, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(sites, fontsize=10)
    ax.set_ylabel("Years showing higher-order synergy  (%)")
    ax.set_ylim(0, 116)
    ax.set_title("Belowground synergy: present in BOTH natural forests\n"
                 "(different subsystem each), absent in managed paddy", fontsize=14)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2,
              fontsize=10.5, frameon=False)
    fig.savefig(OUT/"fig5c_two_forests.png"); plt.close(fig)


PURPLE, TAN = "#7b4fa3", "#b8934a"


# ---------------------------------------------------------------------------
# Fig 6 (NEW centerpiece): flooding — not cultivation — collapses belowground synergy
#   pooled synergy-year fraction of {Rg,Ta,θ,GER} by ecosystem/management group
# ---------------------------------------------------------------------------
def fig_flooding():
    # (label, synergy-years, total-years, color, member note)
    groups = [
        ("Natural\nforest", 19+17+10, 21+20+11, GREEN, "JP-Tak,Tef,Spp"),
        ("Natural\ngrassland", 5+5+4+6, 5+5+4+10, GREEN, "MN-Kbu,Nkh,Skt,CN-HaM"),
        ("Cropland\n(non-flooded)", 2+7, 2+13, ORANGE, "CN-CnR,HbC"),
        ("Rice paddy\n(FLOODED)", 0+0, 8+7, RED, "JP-Mse,KR-CRK"),
    ]
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    for i, (name, k, n, col, note) in enumerate(groups):
        p = 100.0 * k / n
        ax.bar(i, p, color=col, width=0.62, zorder=2)
        ax.text(i, p + 2.5, f"{k}/{n} yr", ha="center", fontsize=12, fontweight="bold")
        ax.text(i, -7, note, ha="center", fontsize=8.5, color="#666")
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g[0] for g in groups], fontsize=11)
    ax.set_ylabel("Years showing belowground higher-order synergy  (%)\n"
                  "soil–respiration system {Rg, Ta, θ, GER}")
    ax.set_ylim(-12, 100)
    ax.axhline(0, color="k", lw=0.8)
    ax.axvspan(-0.5, 1.5, color="#eaf5ea", zorder=0)
    ax.axvspan(1.5, 2.5, color="#fdf2e6", zorder=0)
    ax.axvspan(2.5, 3.5, color="#f7e6e6", zorder=0)
    ax.set_title("It is FLOODING, not cultivation, that collapses the synergy\n"
                 "managed dryland cropland keeps it; only the flooded paddy loses it",
                 fontsize=13.5)
    ax.annotate("flooding pins soil\nmoisture θ → θ carries\nno information →\nsynergy collapses",
                xy=(3, 3), xytext=(2.55, 55), fontsize=9.5, color=RED,
                ha="left", arrowprops=dict(arrowstyle="->", color=RED, lw=1.3))
    fig.savefig(OUT/"fig6_flooding_mechanism.png"); plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 7 (NEW): cross-site light-use decoupling — forests decouple, non-water-limited don't
# ---------------------------------------------------------------------------
def fig_climate_crosssite():
    # (site, r(VPD), significant, group-color)
    rows = [
        ("JP-Yms (F,6y)",  -1.00, True,  GREEN),
        ("JP-Mse (Paddy,8y)", -0.79, True, RED),
        ("JP-Tak (F,21y)", -0.68, True,  GREEN),
        ("RU-Ege (Bor,8y)",-0.62, False, GREY),
        ("JP-Ynf (F,8y)",  -0.60, False, GREEN),
        ("JP-Spp (F,11y)", -0.57, False, GREEN),
        ("JP-Fhk (F,18y)", -0.54, True,  GREEN),
        ("RU-SkP (Bor,7y)",-0.54, False, GREY),
        ("JP-Tef (F,21y)", -0.44, True,  GREEN),
        ("JP-MBF (F,8y)",  -0.45, False, GREEN),
        ("CN-HaM (Grass,10y)", -0.35, False, TAN),
        ("JP-Tmd (F,19y)", -0.12, False, GREEN),
        ("TH-Mae (Trop,8y)", 0.43, False, PURPLE),
        ("MN-Kbu (Grass,5y)", 0.60, False, TAN),
        ("JP-BBY (Bog,5y)", 0.90, False, BLUE),
    ]
    rows = sorted(rows, key=lambda r: r[1])
    labels = [r[0] for r in rows]
    vals = [r[1] for r in rows]
    cols = [r[3] for r in rows]
    y = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    ax.barh(y, vals, color=cols, zorder=2)
    for yi, (lab, v, sig, c) in zip(y, rows):
        if sig:
            ax.text(v - 0.03, yi, "★", va="center", ha="right", color="k", fontsize=12)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9.5)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("Spearman  r :  I(Rg;GPP) vs summer dryness (VPD)\n"
                  "← more negative = stronger decoupling under drought       ★ = significant (p<0.05)")
    ax.set_xlim(-1.1, 1.05)
    ax.set_title("Light-use decoupling is a water-limited-ecosystem trait\n"
                 "forests & paddy decouple; grassland, bog, wet tropics do not",
                 fontsize=13.5)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=GREEN, label="Forest"), Patch(color=RED, label="Rice paddy"),
                       Patch(color=TAN, label="Grassland"), Patch(color=GREY, label="Boreal"),
                       Patch(color=PURPLE, label="Tropical"), Patch(color=BLUE, label="Bog")],
              loc="upper left", fontsize=9, frameon=False, ncol=2)
    fig.savefig(OUT/"fig7_climate_crosssite.png"); plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 5 (updated): higher-order synergy by ecosystem — now cross-biome + 2nd paddy
# ---------------------------------------------------------------------------
def fig_oinfo_crossbiome():
    data = [  # (label, syn-years, total, color)
        ("Deciduous\nforest\nJP-Tak", 19, 21, GREEN),
        ("Evergreen\nforest\nJP-Ta2", 6, 11, GREEN),
        ("Boreal\nlarch\nRU-SkP", 6, 7, GREEN),
        ("Alpine\ngrassland\nCN-HaM", 6, 10, "#66a366"),
        ("Steppe\nMN-Kbu", 5, 5, "#66a366"),
        ("Cropland\nCN-HbC", 7, 13, ORANGE),
        ("Paddy\nJP-Mse", 0, 8, RED),
        ("Paddy\nKR-CRK", 0, 7, RED),
    ]
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    for i, (name, k, n, col) in enumerate(data):
        p = 100.0 * k / n
        ax.bar(i, p, color=col, width=0.68, zorder=2)
        ax.text(i, p + 2.5, f"{k}/{n}", ha="center", fontsize=10.5, fontweight="bold")
    ax.set_xticks(range(len(data)))
    ax.set_xticklabels([d[0] for d in data], fontsize=9)
    ax.set_ylabel("Years showing higher-order synergy  (%)\n{Rg, Ta, θ, GER}")
    ax.set_ylim(0, 108)
    ax.axhline(50, color=GREY, ls=":", lw=1)
    ax.axvspan(-0.5, 4.5, color="#eaf5ea", zorder=0)
    ax.axvspan(4.5, 5.5, color="#fdf2e6", zorder=0)
    ax.axvspan(5.5, 7.5, color="#f7e6e6", zorder=0)
    ax.text(2, 102, "Natural (forest, grassland, boreal)", ha="center", color=GREEN,
            fontsize=10.5, fontweight="bold")
    ax.text(6.5, 102, "Flooded paddy", ha="center", color=RED, fontsize=10.5,
            fontweight="bold")
    ax.set_title("Belowground synergy generalizes across natural biomes,\n"
                 "survives dryland cropping, and collapses only under flooding",
                 fontsize=13.5)
    fig.savefig(OUT/"fig5_oinfo_synergy.png"); plt.close(fig)


# ---------------------------------------------------------------------------
# 概念図 (schematic) — 発表の 動機/背景/方法 用（実データでなく模式）
# ---------------------------------------------------------------------------
def fig_concept_network():
    """§0 動機: ドライバ↔フラックスの相互作用ネットワーク概念図。"""
    fig, ax = plt.subplots(figsize=(8.4, 5.2))
    drivers = {"Rg": (0, 3), "Ta": (0, 2), "VPD": (0, 1), "θ": (0, 0)}
    fluxes = {"GPP": (4, 2.4), "Reco": (4, 1.2), "LE": (4, 0)}
    for n, (x, y) in {**drivers, **fluxes}.items():
        c = "#ffd24d" if n in drivers else "#b8e0b8"
        ax.add_patch(plt.Circle((x, y), 0.34, color=c, ec="#555", zorder=3))
        ax.text(x, y, n, ha="center", va="center", fontsize=12, zorder=4)
    rng = np.random.default_rng(1)
    for _, (x0, y0) in drivers.items():
        for _, (x1, y1) in fluxes.items():
            ax.annotate("", xy=(x1-0.34, y1), xytext=(x0+0.34, y0),
                        arrowprops=dict(arrowstyle="-|>", color="#1f6fb2",
                                        alpha=0.5, lw=1.2))
    # 相互作用（高次）を示す波線ハイライト
    ax.annotate("interactions?\n(synergy / non-additive)", xy=(2, 1.5),
                xytext=(2, 3.5), ha="center", fontsize=12, color=RED,
                fontweight="bold",
                arrowprops=dict(arrowstyle="-[,widthB=6", color=RED, lw=1.5))
    ax.text(0, 3.7, "Drivers", ha="center", fontsize=12, color="#b8860b", fontweight="bold")
    ax.text(4, 3.2, "Fluxes", ha="center", fontsize=12, color=GREEN, fontweight="bold")
    ax.set_xlim(-1, 5.2); ax.set_ylim(-0.8, 4.2); ax.axis("off")
    ax.set_title("Which drivers drive which fluxes — and do they interact?\n"
                 "we measure the interaction structure itself (information flow)",
                 fontsize=13.5)
    fig.savefig(OUT/"fig0_concept_network.png"); plt.close(fig)


def fig_q10_schematic():
    """§1 背景: 呼吸の 分離型(乗法) vs 相乗型(非加法) を、式つき・違いを強調して図示。"""
    T = np.linspace(5, 30, 100)
    DRY, WET = "#c48a3a", "#1f6fb2"
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.0, 5.4), sharey=True)

    # 左: 分離型 R=R0 e^{kT} g(θ) — 同じ傾き k、θ は倍率だけ（高θ = 低θ ×1.8, 相似）
    lowL = np.exp(0.085*(T-5))
    a1.plot(T, 1.0*lowL, color=DRY, lw=2.8, label="dry soil (low θ)")
    a1.plot(T, 1.8*lowL, color=WET, lw=2.8, label="wet soil (high θ)")
    a1.set_title("Separable — the common operational default", fontsize=12.5)
    a1.text(0.97, 0.90, r"$R = R_0\,e^{kT}\,g(\theta)$", transform=a1.transAxes,
            ha="right", fontsize=15)
    a1.text(0.97, 0.80, "k (temp-sensitivity)\nsame for dry & wet",
            transform=a1.transAxes, ha="right", fontsize=10, color="#333")

    # 右: 相乗型 R=R0 e^{k(θ)T} — 傾き k(θ) が θ で変わる。乾は ほぼ平ら、湿は急上昇
    a2.plot(T, np.exp(0.015*(T-5)), color=DRY, lw=2.8, label="dry soil (low θ)")
    a2.plot(T, np.exp(0.11*(T-5)),  color=WET, lw=2.8, label="wet soil (high θ)")
    a2.set_title("Interaction / synergy (we detect it model-free)", fontsize=12.5, color=RED)
    a2.text(0.5, 0.90, r"$R = R_0\,e^{\,k(\theta)\,T}$", transform=a2.transAxes,
            ha="center", fontsize=15, color=RED)
    a2.text(0.5, 0.80, "k(θ) changes with θ → curves fan out",
            transform=a2.transAxes, ha="center", fontsize=10, color=RED)
    a2.annotate("dry: warming\nbarely matters", xy=(24, np.exp(0.015*19)),
                xytext=(14, 2.3), fontsize=10.5, color=DRY,
                arrowprops=dict(arrowstyle="->", color=DRY, lw=1.4))
    a2.annotate("wet: warming\nmatters a lot", xy=(23.5, np.exp(0.11*18.5)),
                xytext=(15.5, 5.8), fontsize=10.5, color=WET,
                arrowprops=dict(arrowstyle="->", color=WET, lw=1.4))

    for a in (a1, a2):
        a.set_xlabel("Soil temperature  Ts  (°C)")
        a.legend(frameon=False, fontsize=10, loc="lower right")
        a.set_yticks([]); a.set_ylim(0, 9)
    a1.set_ylabel("Respiration  Reco")
    fig.suptitle("Do temperature and soil moisture combine separably, or synergistically?\n"
                 r"criterion: $\partial^2\ln R/\partial T\,\partial\theta = 0$ (separable)  vs  "
                 r"$>0$ (synergy)   — we test this with O-information",
                 fontsize=13.5)
    fig.text(0.5, 0.015,
             "Interaction models exist (e.g. DAMM, microbial models) but assume a fixed form; "
             "we detect the interaction structure model-free, in the coupled network, and under management",
             ha="center", fontsize=8.5, color="#666")
    fig.tight_layout(rect=(0, 0.04, 1, 0.92))
    fig.savefig(OUT/"fig_q10_schematic.png"); plt.close(fig)


def fig_pipeline():
    """§3 方法: 解析パイプライン図。"""
    fig, ax = plt.subplots(figsize=(9.8, 3.2))
    steps = [("Flux + met\n(30-min, 11 vars)", "#cfe3f3"),
             ("Preprocess\n5-day anomaly\nlistwise", "#e8e8e8"),
             ("Info theory + causal\nTE · O-information · PCMCI", "#ffe6b3"),
             ("Network · synergy\ncausal skeleton", "#b8e0b8")]
    x = 0
    for i, (txt, c) in enumerate(steps):
        ax.add_patch(plt.Rectangle((x, 0), 2.2, 1.4, color=c, ec="#555"))
        ax.text(x+1.1, 0.7, txt, ha="center", va="center", fontsize=10.5)
        if i < len(steps)-1:
            ax.annotate("", xy=(x+2.55, 0.7), xytext=(x+2.2, 0.7),
                        arrowprops=dict(arrowstyle="-|>", color="#333", lw=1.8))
        x += 2.55
    ax.set_xlim(-0.2, x); ax.set_ylim(-0.3, 1.7); ax.axis("off")
    ax.set_title("Analysis pipeline", fontsize=13.5)
    fig.savefig(OUT/"fig_pipeline.png"); plt.close(fig)


def fig_positioning3():
    """3円ベン図（日本語）: 手法×データ×問い の共通部分＝本研究 ＋ 役立つこと。"""
    import matplotlib.font_manager as fm
    from matplotlib.patches import Circle, FancyBboxPatch
    jp = fm.FontProperties(fname="/usr/share/fonts/opentype/ipafont-gothic/ipag.ttf")
    fig, ax = plt.subplots(figsize=(10.2, 8.4))
    cen = (0, 0.5); R = 2.35
    c1 = (0, 1.9); c2 = (-1.75, -1.0); c3 = (1.75, -1.0)
    cols = ["#1f6fb2", "#2e8b57", "#d98a1e"]
    for (cx, cy), col in zip((c1, c2, c3), cols):
        ax.add_patch(Circle((cx, cy), R, facecolor=col, edgecolor=col, lw=2, alpha=0.22))
    def T(x, y, s, size, color="#222", w="normal", ha="center"):
        ax.text(x, y, s, fontproperties=jp, fontsize=size, color=color,
                ha=ha, va="center", fontweight=w)
    # 円のラベル（外側）
    T(0, 4.55, "① 情報理論・因果推論", 15, cols[0], "bold")
    T(0, 4.0, "相互情報量・PCMCI・O-information", 10.5, cols[0])
    T(-4.35, -2.15, "② 東アジアのフラックス観測", 14.5, cols[1], "bold", ha="left")
    T(-4.35, -2.7, "日中韓の渦相関・水田など管理生態系", 10, cols[1], ha="left")
    T(4.35, -2.15, "③ 生態系の過程の相互作用", 14.5, cols[2], "bold", ha="right")
    T(4.35, -2.7, "炭素・水・エネルギー／相乗・因果／気候・管理", 9.5, cols[2], ha="right")
    # 中心＝本研究
    T(0, 0.95, "本研究", 17, RED, "bold")
    T(0, 0.15, "東アジアの管理生態系のフラックスから、\n過程の相互作用構造を情報理論で読む",
      10.3, "#7a1f14")
    # 下：役立つこと
    ax.add_patch(FancyBboxPatch((-4.7, -4.65), 9.4, 1.35, boxstyle="round,pad=0.03",
                 facecolor="#fbeeea", edgecolor=RED, lw=1.4))
    T(0, -3.62, "これが読めると", 12.5, RED, "bold")
    T(0, -4.12, "① 生態系の状態を診断（気候ストレス・管理・レジームの“指紋”）　"
      "② モデルの過程表現を検証　③ 陸炭素予測の不確実性に観測から迫る", 10.3, "#333")
    ax.set_xlim(-5.2, 5.2); ax.set_ylim(-5.0, 5.0); ax.axis("off")
    ax.set_aspect("equal")
    ax.set_title("手法 × データ × 問い の交わりに本研究がある",
                 fontproperties=jp, fontsize=15)
    fig.savefig(OUT/"fig_positioning3_jp.png", bbox_inches="tight"); plt.close(fig)


def fig_positioning():
    """スライドA: 2つの研究の流れの交点＝本研究、＋突破点。"""
    from matplotlib.patches import Ellipse, FancyBboxPatch
    fig, ax = plt.subplots(figsize=(10.6, 5.8))
    ax.add_patch(Ellipse((3.4, 3.5), 5.4, 3.6, facecolor="#dce8f5",
                         edgecolor="#1f6fb2", lw=2, alpha=0.8))
    ax.add_patch(Ellipse((6.6, 3.5), 5.4, 3.6, facecolor="#e6f0e0",
                         edgecolor="#2e8b57", lw=2, alpha=0.8))
    ax.text(2.0, 4.7, "① Earth-system causal\ninference & information theory",
            ha="center", fontsize=11, color="#1f6fb2", fontweight="bold")
    ax.text(2.0, 3.55, "Runge 2019 · Krich 2020\nGoodwell & Kumar 2020\nRuddell & Kumar 2009 · Rosas 2019",
            ha="center", fontsize=8.5, color="#1f6fb2")
    ax.text(8.0, 4.7, "② Land-carbon uncertainty\n& process models",
            ha="center", fontsize=11, color="#2e8b57", fontweight="bold")
    ax.text(8.0, 3.55, "Arora 2020 · Booth 2012\nDAMM (Davidson 2012)\nFLUXCOM (Jung 2020)",
            ha="center", fontsize=8.5, color="#2e8b57")
    # intersection
    ax.text(5.0, 4.15, "THIS STUDY", ha="center", fontsize=12.5, color=RED, fontweight="bold")
    ax.text(5.0, 3.35,
            "read the interaction\nstructure from data,\ntest model process\nassumptions",
            ha="center", fontsize=9.2, color="#7a1f14")
    # breakthroughs bar
    ax.add_patch(FancyBboxPatch((0.4, 0.35), 9.6, 1.15, boxstyle="round,pad=0.02",
                 facecolor="#fbeeea", edgecolor=RED, lw=1.3))
    ax.text(5.2, 1.28, "Breakthroughs", ha="center", fontsize=10.5, color=RED, fontweight="bold")
    ax.text(5.2, 0.78,
            "East Asia & humid monsoon   ·   flooded rice paddy (management)   ·   "
            "system-level synergy (O-information)   ·   function-form-free   ·   6 sites / 4 biomes",
            ha="center", fontsize=9.3, color="#333")
    ax.set_xlim(0, 10); ax.set_ylim(0, 5.6); ax.axis("off")
    ax.set_title("Where this study sits: the tools of ① applied to the problem of ②",
                 fontsize=13)
    fig.savefig(OUT/"fig_positioning.png"); plt.close(fig)


def fig_uncertainty():
    """スライド1: 陸の炭素-気候フィードバックは海より大きく・不確実（CMIP6, Arora 2020）
    ＋その不確実性の原因＝本研究が測る結合、を1枚で。"""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.4, 5.0),
                                 gridspec_kw={"width_ratios": [1.15, 1]})
    # 左: land vs ocean feedback γ with error bars
    names = ["Ocean", "Land"]
    val = [-17.2, -45.1]; err = [5.0, 50.6]; col = [BLUE, GREEN]
    y = [0, 1]
    a1.barh(y, val, color=col, height=0.5, zorder=2)
    a1.errorbar(val, y, xerr=err, fmt="none", ecolor="#222", capsize=7, lw=1.8, zorder=3)
    a1.axvline(0, color="k", lw=0.8)
    a1.set_yticks(y); a1.set_yticklabels(names, fontsize=12)
    a1.text(-17.2, 0.32, "−17.2 ± 5.0", ha="center", fontsize=10.5, color=BLUE)
    a1.text(-45.1, 1.32, "−45.1 ± 50.6", ha="center", fontsize=10.5, color=GREEN, fontweight="bold")
    a1.set_xlabel("Carbon–climate feedback  γ  (PgC per °C)")
    a1.set_title("Land feedback: ~3× larger, spread ~10× wider\n(CMIP6, Arora et al. 2020)",
                 fontsize=12)
    a1.set_ylim(-0.6, 1.7)
    a1.annotate("huge model spread\n(crosses zero)", xy=(-45.1-50.6, 1), xytext=(-88, 0.55),
                fontsize=9.5, color=GREEN,
                arrowprops=dict(arrowstyle="->", color=GREEN, lw=1.2))
    # 右: why uncertain -> our 3 couplings
    a2.axis("off")
    a2.text(0.5, 0.96, "Why is the land sink so uncertain?", ha="center",
            fontsize=12.5, fontweight="bold", transform=a2.transAxes)
    boxes = [
        ("① Soil-respiration temperature sensitivity\n(how much CO₂ soils release when warm)",
         "Booth et al. 2012"),
        ("② Photosynthesis under water stress\n(soil moisture drives ~90% of sink IAV)",
         "Humphrey et al. 2021"),
        ("③ Is the year-to-year sink driven by\ntemperature or by water?",
         "Jung et al. 2017"),
    ]
    yb = 0.80
    for txt, cite in boxes:
        a2.add_patch(plt.Rectangle((0.03, yb-0.13), 0.94, 0.15, transform=a2.transAxes,
                     facecolor="#f0f4ea", edgecolor="#9bbf7f"))
        a2.text(0.06, yb-0.03, txt, fontsize=9.7, va="top", transform=a2.transAxes)
        a2.text(0.95, yb-0.115, cite, fontsize=8, ha="right", color="#666",
                style="italic", transform=a2.transAxes)
        yb -= 0.205
    a2.annotate("", xy=(0.5, 0.16), xytext=(0.5, 0.20), transform=a2.transAxes,
                arrowprops=dict(arrowstyle="-|>", color=RED, lw=2))
    a2.text(0.5, 0.13, "= the driver–flux couplings\nTHIS study measures (TE / O-information)",
            ha="center", fontsize=10.5, color=RED, fontweight="bold", va="top",
            transform=a2.transAxes)
    fig.suptitle("Land carbon uptake is a top uncertainty in climate projection — "
                 "and its causes are the couplings we measure", fontsize=12.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT/"fig_uncertainty.png"); plt.close(fig)


if __name__ == "__main__":
    fig_climate(); fig_pid(); fig_skeleton(); fig_robustness()
    fig_conditioning(); fig_oinfo_ci(); fig_oinfo_twosub()
    fig_flooding(); fig_climate_crosssite(); fig_oinfo_crossbiome()
    fig_concept_network(); fig_q10_schematic(); fig_pipeline(); fig_uncertainty()
    fig_positioning(); fig_positioning3()
    for p in sorted(OUT.glob("*.png")):
        print(f"[fig] {p}  ({p.stat().st_size//1024} KB)")
