"""ネットワーク図と診断プロット (R&K Figure 6 相当)。

- :func:`draw_network`   : 有意な結合のみ、矢印付き、結合タイプで線種を変える有向グラフ
- :func:`plot_type_matrix`: 結合タイプ 1/2/3/4 のヒートマップ
- :func:`plot_lag_diagnostics`: 主要ペアの T' と I' のラグ依存性
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm
import networkx as nx
import numpy as np

from .config import RK_LABELS, RK_VARS
from .network import NetworkResult

# 結合タイプ → 線種・色 (図の凡例)
TYPE_STYLE = {
    1: {"style": "dotted", "color": "#7f7f7f", "label": "Type 1 (sync)"},
    2: {"style": "solid", "color": "#1f77b4", "label": "Type 2 (feedback)"},
    3: {"style": "dashed", "color": "#d62728", "label": "Type 3 (forcing)"},
}


def _circular_layout(nodes: list[str]) -> dict:
    n = len(nodes)
    ang = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, n, endpoint=False)
    return {v: (np.cos(a), np.sin(a)) for v, a in zip(nodes, ang)}


def draw_network(net: NetworkResult, path: str | Path, title: str | None = None):
    """有意な結合 (Type 2/3) を矢印で描く。線種はタイプ、線幅は Tz、注記は τ'。"""
    pos = _circular_layout(RK_VARS)
    G = nx.DiGraph()
    G.add_nodes_from(RK_VARS)

    fig, ax = plt.subplots(figsize=(9, 9))
    nx.draw_networkx_nodes(G, pos, node_color="#eeeeee", edgecolors="#333333",
                           node_size=1500, ax=ax)
    nx.draw_networkx_labels(G, pos, labels=RK_LABELS, font_size=11, ax=ax)

    tz_vals = net.ATz.to_numpy(dtype=float)
    tz_max = np.nanmax(tz_vals) if np.isfinite(tz_vals).any() else 1.0
    for src in RK_VARS:
        for dst in RK_VARS:
            if src == dst:
                continue
            ctype = int(net.ctype.loc[src, dst])
            if ctype not in (2, 3):
                continue
            tz = net.ATz.loc[src, dst]
            lag = net.Gamma.loc[src, dst]
            style = TYPE_STYLE[ctype]
            width = 1.0 + 3.0 * (tz / tz_max if tz_max > 0 else 0.0)
            ax.annotate(
                "", xy=pos[dst], xytext=pos[src],
                arrowprops=dict(
                    arrowstyle="-|>", color=style["color"],
                    linestyle=style["style"], linewidth=width,
                    shrinkA=22, shrinkB=22,
                    connectionstyle="arc3,rad=0.12",
                ),
            )
            mx, my = (pos[src][0] + pos[dst][0]) / 2, (pos[src][1] + pos[dst][1]) / 2
            ax.text(mx, my, f"{lag:.1f}h", fontsize=7, color=style["color"],
                    ha="center", va="center")

    handles = [
        Line2D([0], [0], color=s["color"], linestyle=s["style"], lw=2, label=s["label"])
        for s in TYPE_STYLE.values()
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=9, framealpha=0.9)
    ax.set_title(title or _default_title(net), fontsize=12)
    ax.set_axis_off()
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_type_matrix(net: NetworkResult, path: str | Path):
    """結合タイプ行列 (1/2/3/4) のヒートマップ。行=source, 列=target。"""
    cmap = ListedColormap(["#ffffff", "#bbbbbb", "#1f77b4", "#d62728", "#f5f5f5"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5, 4.5], cmap.N)
    mat = net.ctype.to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    ax.imshow(mat, cmap=cmap, norm=norm, aspect="equal")
    labels = [RK_LABELS[v] for v in RK_VARS]
    ax.set_xticks(range(len(RK_VARS)), labels, rotation=45, ha="right")
    ax.set_yticks(range(len(RK_VARS)), labels)
    ax.set_xlabel("target (Y)")
    ax.set_ylabel("source (X)")
    for i in range(len(RK_VARS)):
        for j in range(len(RK_VARS)):
            v = int(mat[i, j])
            if v in (2, 3):
                ax.text(j, i, str(v), ha="center", va="center", color="white", fontsize=8)
    handles = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor=c, markersize=12, label=lab)
        for c, lab in [("#bbbbbb", "1 sync"), ("#1f77b4", "2 feedback"),
                       ("#d62728", "3 forcing"), ("#f5f5f5", "4 uncoupled")]
    ]
    ax.legend(handles=handles, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    ax.set_title(f"Coupling types — {net.meta.get('site')} "
                 f"{net.meta.get('year')}-{net.meta.get('month'):02d}", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def plot_lag_diagnostics(net: NetworkResult, pairs: list[tuple[str, str]],
                         path: str | Path):
    """主要ペアの T'(τ) と有意しきい値 Δ(T')、および I' を描く (R&K Fig.6 相当)。"""
    cfg = net.config
    lags_h = [cfg.lag_hours(t) for t in net.lags]
    n = len(pairs)
    ncol = 2
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(6 * ncol, 3.2 * nrow), squeeze=False)

    for k, (src, dst) in enumerate(pairs):
        ax = axes[k // ncol][k % ncol]
        curve = net.te_curves.get((src, dst))
        thr = net.te_threshold.get((src, dst))
        if curve is None:
            ax.set_visible(False)
            continue
        tp = curve / cfg.log_m
        th = thr / cfg.log_m
        ax.plot(lags_h, tp, "-o", ms=3, color="#d62728", label="T'(X→Y)")
        ax.plot(lags_h, th, "--", color="#999999", label="Δ(T')")
        Ip = net.AI.loc[src, dst] / 100.0  # I' (0-1)
        ax.axhline(Ip, color="#1f77b4", ls=":", label="I'(X,Y)")
        gamma = net.Gamma.loc[src, dst]
        if np.isfinite(gamma):
            ax.axvline(gamma, color="#333333", lw=0.8)
            ax.text(gamma, ax.get_ylim()[1] * 0.92, f"τ'={gamma:.1f}h",
                    fontsize=7, ha="left")
        ax.set_title(f"{RK_LABELS[src]} → {RK_LABELS[dst]}", fontsize=10)
        ax.set_xlabel("lag τ [h]")
        ax.set_ylabel("normalized info")
        ax.legend(fontsize=7)

    for k in range(n, nrow * ncol):
        axes[k // ncol][k % ncol].set_visible(False)
    fig.suptitle(f"Lag diagnostics — {net.meta.get('site')} "
                 f"{net.meta.get('year')}-{net.meta.get('month'):02d}", fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _default_title(net: NetworkResult) -> str:
    m = net.meta
    return (f"Process network — {m.get('site')} {m.get('year')}-{m.get('month'):02d}  "
            f"(n={m.get('n_points')}, m={m.get('n_bins')})")
