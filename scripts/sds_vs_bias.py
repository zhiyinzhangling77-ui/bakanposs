"""Unify SDS (in-situ drought sensitivity) with satellite ET bias.

The paper's logic:
  SDS measures how strongly EC ET tracks surface SWC.
  Satellite ET retrievals (MOD16, PML) effectively assume that ET is
  closely tied to surface moisture and LAI/canopy state.
  Therefore strata where SDS is small (decoupled) should be the strata
  where the satellite product underestimates EC the most.

For each (site, season, irrig_bucket) cell with sufficient data we
compute:
  SDS         = 1 - mean(LE | dry SWC) / mean(LE | normal SWC)
  bias_MOD16  = mean(MOD16_ET - EC_ET)
  bias_PML    = mean(PML_ET   - EC_ET)
and plot bias vs SDS — the expected pattern is a positive slope (low
SDS → strongly negative bias).

Reads: master_full.csv
Writes:
  sds_vs_bias.csv
  figs/fig_E_sds_vs_bias.png
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

REPO = Path(__file__).parent.parent
CSV = REPO / "master_full.csv"
OUT_CSV = REPO / "sds_vs_bias.csv"
FIGS = REPO / "figs"
FIGS.mkdir(exist_ok=True)

NDVI_GATE = 0.3
MIN_N = 20
RNG = np.random.default_rng(42)
N_BOOT = 1000


def sds_from_arrays(le: np.ndarray, swc: np.ndarray,
                    p25: float, p75: float) -> float:
    dry = swc < p25
    norm = (swc >= p25) & (swc <= p75)
    if dry.sum() < 5 or norm.sum() < 5:
        return np.nan
    m_norm = np.nanmean(le[norm])
    if not np.isfinite(m_norm) or m_norm <= 0:
        return np.nan
    return 1.0 - np.nanmean(le[dry]) / m_norm


def bootstrap_sds(le: np.ndarray, swc: np.ndarray,
                  p25: float, p75: float) -> tuple[float, float, float]:
    idx = np.arange(len(le))
    samples = []
    for _ in range(N_BOOT):
        b = RNG.choice(idx, size=len(idx), replace=True)
        samples.append(sds_from_arrays(le[b], swc[b], p25, p75))
    s = np.array(samples)
    s = s[np.isfinite(s)]
    if len(s) < 50:
        return np.nan, np.nan, np.nan
    return (float(np.median(s)),
            float(np.quantile(s, 0.025)),
            float(np.quantile(s, 0.975)))


def cell(df: pd.DataFrame, label: str) -> dict | None:
    if len(df) < MIN_N:
        return None
    # Prefer LE if mostly populated; otherwise fall back to ET (Oran daily lacks LE)
    le_series = df["LE_Wm2"] if "LE_Wm2" in df else df["ET_mm"]
    if le_series.notna().sum() < 30:
        le_series = df["ET_mm"]
    le = le_series.to_numpy()
    swc = df["SWC"].to_numpy()
    ok = np.isfinite(le) & np.isfinite(swc)
    if ok.sum() < MIN_N:
        return None
    le = le[ok]; swc = swc[ok]
    p25 = float(np.nanquantile(swc, 0.25))
    p75 = float(np.nanquantile(swc, 0.75))
    sds_med, sds_lo, sds_hi = bootstrap_sds(le, swc, p25, p75)
    out = {
        "stratum": label,
        "n": int(ok.sum()),
        "SDS": sds_med, "SDS_lo": sds_lo, "SDS_hi": sds_hi,
    }
    for col, label_p in [("MOD16_ET", "MOD16"), ("PML_ET", "PML")]:
        if col not in df:
            continue
        bias = (df[col] - df["ET_mm"]).dropna()
        if len(bias) >= 5:
            out[f"bias_{label_p}_mean"] = float(bias.mean())
            out[f"bias_{label_p}_n"] = int(len(bias))
            out[f"bias_{label_p}_se"] = float(bias.std(ddof=1)
                                              / np.sqrt(len(bias)))
    return out


def main() -> None:
    df = pd.read_csv(CSV, parse_dates=["date"])
    df = df[df["NDVI"].fillna(0) > NDVI_GATE].copy()

    rows = []

    for season in ["spring", "summer", "fall"]:
        sub = df[(df.site == "Oran") & (df.season == season)]
        c = cell(sub, f"Oran_{season}")
        if c:
            c["site"] = "Oran"
            c["bucket"] = season
            rows.append(c)

    tzm = df[df.site == "TzM"]
    for season in ["spring", "summer", "fall"]:
        sub = tzm[tzm.season == season]
        c = cell(sub, f"TzM_{season}")
        if c:
            c["site"] = "TzM"
            c["bucket"] = season
            rows.append(c)

    tzm_summer = tzm[tzm.season == "summer"]
    for bk in ["irrig_d0-3", "irrig_d4-7", "irrig_d8plus"]:
        sub = tzm_summer[tzm_summer.irrig_bucket == bk]
        c = cell(sub, f"TzM_summer_{bk}")
        if c:
            c["site"] = "TzM"
            c["bucket"] = bk
            rows.append(c)

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    pd.set_option("display.float_format", "{:+.3f}".format)
    print(out.to_string(index=False))
    print(f"\nwrote {OUT_CSV}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    for ax, (col, label, c) in zip(
            axes, [("bias_MOD16_mean", "MOD16", "tab:purple"),
                   ("bias_PML_mean", "PML", "tab:green")]):
        for site, marker in [("Oran", "o"), ("TzM", "s")]:
            sub = out[out.site == site].dropna(subset=["SDS", col])
            if not len(sub):
                continue
            ax.errorbar(sub["SDS"], sub[col],
                        xerr=[sub["SDS"] - sub["SDS_lo"],
                              sub["SDS_hi"] - sub["SDS"]],
                        yerr=sub.get(col.replace("_mean", "_se"), 0),
                        fmt=marker, color=c if site == "TzM" else "gray",
                        markersize=10, alpha=0.85, capsize=3, label=site)
            for _, r in sub.iterrows():
                ax.annotate(r["bucket"],
                            xy=(r["SDS"], r[col]),
                            xytext=(5, 5), textcoords="offset points",
                            fontsize=8)
        ax.axhline(0, color="k", lw=0.7, ls="--")
        ax.axvline(0, color="k", lw=0.7, ls="--")
        ax.set_xlabel("SDS  (1 − mean LE|dry / mean LE|normal)")
        ax.set_ylabel(f"mean ({label} − EC) ET  [mm/d]")
        ax.set_title(f"{label}: bias vs in-situ drought sensitivity")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
    fig.suptitle("Strata with low SDS (decoupling) show the largest "
                 "satellite underestimation",
                 fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "fig_E_sds_vs_bias.png", dpi=150)
    plt.close(fig)
    print(f"wrote {FIGS / 'fig_E_sds_vs_bias.png'}")


if __name__ == "__main__":
    main()
