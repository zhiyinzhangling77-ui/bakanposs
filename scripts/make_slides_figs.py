"""発表用の図を、これまでの解析結果の数値から生成する（実データ不要・数値埋め込み）。

英語ラベル（学会標準）。出力は japanflux_pn/slides/*.png。説明は SLIDES ガイド参照。
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 130, "font.size": 13,
    "axes.titlesize": 15, "axes.labelsize": 13, "axes.spines.top": False,
    "axes.spines.right": False, "figure.autolayout": True,
})
OUT = Path("japanflux_pn/slides")
OUT.mkdir(parents=True, exist_ok=True)

BLUE, RED, GREEN, ORANGE, GREY = "#1f6fb2", "#c0392b", "#2e8b57", "#e08a1e", "#888"


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
    pos = {"Rg":(0,2), "VPD":(-1.6,1), "Ta":(0.2,1), "gLE":(-1.9,-0.2),
           "gH":(-0.6,-0.3), "Ts":(1.1,0.2), "GEP":(1.7,1.1), "NEE":(2.2,0)}
    edges = [("Rg","VPD"),("Rg","Ta"),("Rg","gLE"),("Rg","gH"),("Rg","Ts"),
             ("Ta","Ts"),("VPD","gLE"),("GEP","NEE")]
    G = nx.DiGraph(); G.add_nodes_from(pos); G.add_edges_from(edges)
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    ncol = ["#ffd24d" if n=="Rg" else ("#b8e0b8" if n in ("GEP","NEE")
            else "#cfe3f3") for n in G.nodes]
    nx.draw_networkx_nodes(G, pos, node_size=1500, node_color=ncol,
                           edgecolors="#555", ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=12, ax=ax)
    nx.draw_networkx_edges(G, pos, arrowstyle="-|>", arrowsize=20,
                           width=2, edge_color="#1f6fb2",
                           node_size=1500, ax=ax)
    ax.set_title("Causal skeleton after removing common driving (PCMCI)\n"
                 "Radiation drives climate/energy; photosynthesis drives net carbon")
    ax.axis("off")
    fig.savefig(OUT/"fig3_causal_skeleton.png"); plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 4: multi-year robustness — link consistency (core vs artifact)
# ---------------------------------------------------------------------------
def fig_robustness():
    core = [("GEP->NEE",100),("Rg->gH",95),("Rg->gLE",95),("Rg->VPD",95),
            ("Ta->Ts",90),("Rg->Ta",86),("Rg->Ts",76),("Ta->GER",71)]
    art = [("gH->Rg",57),("gLE->Rg",38),("GEP->Rg",14),("NEE->Rg",14)]
    labels = [l for l,_ in core] + [l for l,_ in art]
    vals = [v for _,v in core] + [v for _,v in art]
    cols = [GREEN]*len(core) + [RED]*len(art)
    y = np.arange(len(labels))[::-1]
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    ax.barh(y, vals, color=cols)
    ax.set_yticks(y); ax.set_yticklabels(labels)
    ax.axvline(70, color=GREY, ls="--", lw=1)
    ax.set_xlabel("Fraction of years the causal link appears (%)  — JP-Tak, 21 yr")
    ax.set_title("Real causal links are stable across years;\n"
                 "physically impossible links (-> Rg) are sporadic = artifacts")
    ax.set_xlim(0, 105)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=GREEN, label="Physical core links"),
                       Patch(color=RED, label="Into radiation Rg (impossible) = artifact")],
              loc="lower right", fontsize=11, frameon=False)
    fig.savefig(OUT/"fig4_robustness.png"); plt.close(fig)


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
    flux = [("gH-GEP",10.9,1.4),("gLE-GEP",11.2,2.3),
            ("gH-gLE",8.7,1.1),("NEE-GEP",21.9,6.9)]
    therm = [("Ta-Ts",19.1,21.8),("Ta-GER",17.5,20.7),
             ("Ts-GER",13.6,18.7),("th-GER",7.8,12.6)]
    labels = [l for l,_,_ in flux] + [l for l,_,_ in therm]
    before = [b for _,b,_ in flux] + [b for _,b,_ in therm]
    after  = [a for _,_,a in flux] + [a for _,_,a in therm]
    x = np.array([0,1,2,3, 4.8,5.8,6.8,7.8])
    w = 0.38
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    ax.bar(x-w/2, before, w, label="Coupling I(X;Y)  (before)", color=BLUE)
    ax.bar(x+w/2, after,  w, label="After removing radiation  I(X;Y | Rg)",
           color="#bcd4e8", edgecolor="#7fa8cc")
    ax.set_xticks(x); ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Coupling strength  I  (%)")
    ax.set_ylim(0, 30)
    ax.set_title("Removing the common driver (radiation Rg):\n"
                 "flux couplings collapse — temperature–respiration survive",
                 fontsize=14)
    ax.axvspan(-0.6, 3.6, color="#eef4fa", zorder=0)
    ax.axvspan(4.2, 8.4, color="#fdf2e6", zorder=0)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2,
              fontsize=10.5, frameon=False)
    ax.text(1.5, 27.5, "Energy / carbon fluxes\n(the sun's shadow)", ha="center",
            fontsize=11, color=BLUE, fontweight="bold")
    ax.text(6.3, 27.5, "Temperature – respiration\n(real, direct)", ha="center",
            fontsize=11, color=ORANGE, fontweight="bold")
    for xi, b, a in zip(x[:4], before[:4], after[:4]):
        ax.annotate("", xy=(xi+w/2, a+0.6), xytext=(xi-w/2, b-0.6),
                    arrowprops=dict(arrowstyle="->", color=RED, lw=1.4))
    ax.text(1.5, 4.5, "collapse", color=RED, ha="center", fontsize=11,
            style="italic")
    fig.savefig(OUT/"fig2b_conditioning.png"); plt.close(fig)


if __name__ == "__main__":
    fig_climate(); fig_pid(); fig_skeleton(); fig_robustness(); fig_oinfo()
    fig_conditioning()
    for p in sorted(OUT.glob("*.png")):
        print(f"[fig] {p}  ({p.stat().st_size//1024} KB)")
