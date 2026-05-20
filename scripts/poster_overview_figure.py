"""Poster overview figure: NDVI, GPP, ET dynamics at two flux towers.

Each indicator is loaded from its INDEPENDENT source so that satellite
lines remain plotted on EC-missing days. master_full_v2.csv (EC-joined)
was discarding satellite values on EC-gap days; we bypass that here.

Sources:
  EC tower fluxes/NDVI      : data/master_full_v2.csv      (EC days only)
  Sentinel-2 NDVI           : data/OranTzM_S2_NDVI.csv     (B4/B8 -> NDVI;
                              full S2 revisit, independent of EC)
  Meteosat MGPP (10-day GPP): data/mgpp_decadal_all.csv    (independent)
  Meteosat ETv3 (daily ET)  : data/metv3_daily_all.csv     (independent)

Layout (3 rows x 2 cols):
  Row 1: NDVI (EC tower vs Sentinel-2)
  Row 2: GPP  (EC tower vs Meteosat MGPP)
  Row 3: ET   (EC tower vs Meteosat ETv3)

Spatial / temporal scales are annotated in each legend entry.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

REPO = Path(__file__).parent.parent
EC_CSV    = REPO / "data" / "master_full_v2.csv"
METV3_CSV = REPO / "data" / "metv3_daily_all.csv"
MGPP_CSV  = REPO / "data" / "mgpp_decadal_all.csv"
S2_CSV    = REPO / "data" / "OranTzM_S2_NDVI.csv"
OUT_DIR = REPO / "figures" / "poster"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SITES = ["Oran", "TzM"]
SITE_LABEL = {
    "Oran": "Oran  ·  rainfed vetch  ·  2018–2020",
    "TzM":  "Tarazona  ·  drip-irrig. almond  ·  2020–2024",
}

C_EC      = "#1a1a1a"
C_EC_RAW  = "#c8c8c8"
C_METV3   = "#d1493c"
C_MGPP    = "#2a8c63"
C_S2      = "#1e6fb3"
SHADE     = "#f6e7b8"
GRID      = "#dddddd"


def smooth(s: pd.Series, w: int = 7) -> pd.Series:
    return s.rolling(w, center=True, min_periods=max(3, w // 2)).mean()


def smooth_break(s: pd.Series, gap_days: int = 7, w: int = 7) -> pd.Series:
    """Centered 7-day mean, but NaN inside any run of >gap_days missing days."""
    out = s.rolling(w, center=True, min_periods=max(3, w // 2)).mean()
    miss = s.isna()
    streak = miss.groupby((~miss).cumsum()).cumsum()
    return out.where(~(streak > gap_days))


def break_gaps(dates, values, max_gap_days: int = 15):
    """Return (dates, values) with NaN inserted where consecutive valid
    observations are more than max_gap_days apart, so matplotlib breaks
    the line rather than bridging gaps."""
    d = pd.to_datetime(pd.Series(dates))
    v = pd.Series(values, index=d.index)
    keep = v.notna() & d.notna()
    d = d[keep].reset_index(drop=True)
    v = v[keep].reset_index(drop=True)
    order = d.argsort()
    d = d.iloc[order].reset_index(drop=True)
    v = v.iloc[order].reset_index(drop=True)
    if len(d) == 0:
        return d.to_numpy(), v.to_numpy()
    gaps = d.diff().dt.days
    out_d, out_v = [], []
    for i in range(len(d)):
        if i > 0 and gaps.iloc[i] > max_gap_days:
            out_d.append(d.iloc[i - 1] + pd.Timedelta(days=1))
            out_v.append(np.nan)
        out_d.append(d.iloc[i])
        out_v.append(v.iloc[i])
    return (pd.to_datetime(out_d).to_numpy(),
            np.asarray(out_v, dtype=float))


def reindex_daily(df: pd.DataFrame, d0, d1) -> pd.DataFrame:
    full = pd.date_range(d0, d1, freq="D")
    return df.set_index("date").reindex(full).rename_axis("date").reset_index()


def shade_springs(ax, date_min, date_max):
    for y in range(date_min.year, date_max.year + 1):
        a = pd.Timestamp(f"{y}-03-01")
        b = pd.Timestamp(f"{y}-06-30")
        if b < date_min or a > date_max:
            continue
        ax.axvspan(max(a, date_min), min(b, date_max),
                   color=SHADE, alpha=0.55, zorder=0)


def style_time_ax(ax, d0, d1):
    ax.set_xlim(d0, d1)
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.grid(which="major", color=GRID, lw=0.8, zorder=0)


def load_s2_ndvi() -> pd.DataFrame:
    """Compute NDVI from B4/B8 and aggregate to one value per date+site."""
    raw = pd.read_csv(S2_CSV, parse_dates=["date"])
    raw = raw.rename(columns={"name": "site"})
    raw["NDVI_s2"] = (raw["B8"] - raw["B4"]) / (raw["B8"] + raw["B4"])
    raw = raw.dropna(subset=["NDVI_s2"])
    return (raw.groupby(["site", "date"], as_index=False)["NDVI_s2"]
               .mean())


def _to_daily(series_dt: pd.Series, values: pd.Series,
              d0, d1, interp_limit_days: int) -> pd.Series:
    """Sparse (date, value) -> daily-indexed Series with limited linear
    interpolation across short gaps."""
    s = pd.Series(values.values, index=pd.to_datetime(series_dt.values))
    s = s[~s.index.duplicated(keep="first")].sort_index()
    full = pd.date_range(d0, d1, freq="D")
    s = s.reindex(full)
    s = s.interpolate(method="linear", limit=interp_limit_days,
                      limit_direction="both")
    return s


def phase_stats(ec_daily: pd.Series, sat_daily: pd.Series,
                max_lag_days: int = 120,
                smooth_w: int = 30):
    """Compute lag and peak-DOY statistics between EC and Satellite.

    Both inputs are daily-indexed pd.Series. They are smoothed with a
    ``smooth_w``-day centered rolling mean before peak detection so that
    Sentinel-2's noisy revisits don't dominate.

    Returns dict with:
      dpeak_doy : mean(DOY_max(Sat) - DOY_max(EC)) over years where both
                  have a real annual peak (positive => Sat peaks later)
      n_years   : number of years contributing to dpeak_doy
      lag       : argmax tau in days for corr(EC[t], Sat[t+tau]) flipped
                  so positive => Sat is later than EC (consistent with
                  dpeak sign)
      r_max     : Pearson correlation at the best lag
    """
    ec_s  = ec_daily.rolling(smooth_w, center=True,
                              min_periods=max(5, smooth_w // 3)).mean()
    sat_s = sat_daily.rolling(smooth_w, center=True,
                              min_periods=max(5, smooth_w // 3)).mean()
    ec_s, sat_s = ec_s.align(sat_s, join="outer")

    # ---- cross-correlation lag scan ----
    best_r, best_tau = -np.inf, 0
    for tau in range(-max_lag_days, max_lag_days + 1):
        shifted = sat_s.shift(tau)
        df = pd.concat([ec_s, shifted], axis=1).dropna()
        if len(df) < 90:
            continue
        r = df.iloc[:, 0].corr(df.iloc[:, 1])
        if not np.isnan(r) and r > best_r:
            best_r, best_tau = r, tau
    # shift(+tau): sat value at index i moves to i+tau, i.e., sat appears
    # tau days LATER. If best_tau>0 makes them match, Sat was originally
    # EARLIER than EC. Flip the sign so positive => Sat later than EC.
    lag = -best_tau
    # Flag lags that land on the boundary -> not a real peak inside window
    lag_at_boundary = abs(best_tau) >= max_lag_days - 1

    # ---- annual peak DOY ----
    years_ec  = set(ec_s.dropna().index.year)
    years_sat = set(sat_s.dropna().index.year)
    diffs = []
    for y in sorted(years_ec & years_sat):
        ec_y  = ec_s[ec_s.index.year == y].dropna()
        sat_y = sat_s[sat_s.index.year == y].dropna()
        if len(ec_y) < 60 or len(sat_y) < 60:
            continue
        diffs.append(sat_y.idxmax().dayofyear - ec_y.idxmax().dayofyear)
    dpeak = float(np.mean(diffs)) if diffs else float("nan")
    dpeak_per_year = diffs

    return {
        "dpeak_doy": dpeak, "n_years": len(diffs),
        "lag": lag, "r_max": best_r,
        "lag_at_boundary": lag_at_boundary,
        "dpeak_per_year": dpeak_per_year,
    }


def add_phase_annotation(ax, st: dict, loc: str = "upper right") -> None:
    """Small text box: Δpeak / lag (positive => Sat later) and r_max."""
    if loc == "upper right":
        x, y, ha, va = 0.985, 0.965, "right", "top"
    else:
        x, y, ha, va = 0.015, 0.965, "left", "top"
    dp = st["dpeak_doy"]
    dp_txt = f"{dp:+.0f} d" if not np.isnan(dp) else "n/a"
    boundary_mark = "*" if st.get("lag_at_boundary") else ""
    txt = (f"$\\Delta$peak = {dp_txt}\n"
           f"lag = {st['lag']:+d} d{boundary_mark}\n"
           f"r = {st['r_max']:.2f}")
    ax.text(x, y, txt, transform=ax.transAxes, ha=ha, va=va,
            fontsize=8.4, family="monospace",
            bbox=dict(boxstyle="round,pad=0.32", fc="white",
                      ec="#bbbbbb", lw=0.6, alpha=0.92))


def main() -> None:
    ec_full    = pd.read_csv(EC_CSV,    parse_dates=["date"])
    metv3_full = pd.read_csv(METV3_CSV, parse_dates=["date"])
    mgpp_full  = pd.read_csv(MGPP_CSV,  parse_dates=["date"])
    mgpp_full  = mgpp_full[(mgpp_full["GPP_mgpp"] > 0) &
                           (mgpp_full["GPP_mgpp"] < 30)].copy()
    s2_full    = load_s2_ndvi()

    fig, axes = plt.subplots(
        3, 2, figsize=(13.5, 9.0),
        gridspec_kw={"hspace": 0.32, "wspace": 0.18},
    )
    panel_labels = [["a", "b"], ["c", "d"], ["e", "f"]]
    stat_rows = []  # collects phase stats for CSV export

    for j, site in enumerate(SITES):
        ec_site    = ec_full   [ec_full   ["site"] == site].sort_values("date")
        metv3_site = metv3_full[metv3_full["site"] == site].sort_values("date")
        mgpp_site  = mgpp_full [mgpp_full ["site"] == site].sort_values("date")
        s2_site    = s2_full   [s2_full   ["site"] == site].sort_values("date")

        # x-range = EC observation window per site
        d0, d1 = ec_site["date"].min(), ec_site["date"].max()

        # reindex each source to its OWN daily continuous frame
        ec_d    = reindex_daily(ec_site,    d0, d1)
        metv3_d = reindex_daily(metv3_site[(metv3_site["date"] >= d0) &
                                            (metv3_site["date"] <= d1)],
                                 d0, d1)
        mgpp_d  = mgpp_site[(mgpp_site["date"] >= d0) &
                            (mgpp_site["date"] <= d1)].copy()
        s2_d    = s2_site[(s2_site["date"] >= d0) &
                          (s2_site["date"] <= d1)].copy()

        # Build daily-indexed series for phase stats (interp small gaps)
        ec_daily_index = pd.DatetimeIndex(ec_d["date"])
        ec_ndvi_dly = pd.Series(ec_d["NDVI"].values, index=ec_daily_index)
        ec_gpp_dly  = pd.Series(ec_d["GPP_gC_m2_d"].values, index=ec_daily_index)
        ec_et_dly   = pd.Series(ec_d["ET_mm"].values, index=ec_daily_index)
        s2_dly      = _to_daily(s2_d["date"], s2_d["NDVI_s2"], d0, d1,
                                 interp_limit_days=20)
        mgpp_dly    = _to_daily(mgpp_d["date"], mgpp_d["GPP_mgpp"], d0, d1,
                                 interp_limit_days=15)
        metv3_dly   = pd.Series(metv3_d["ET_mm"].values,
                                 index=pd.DatetimeIndex(metv3_d["date"]))

        st_ndvi = phase_stats(ec_ndvi_dly, s2_dly)
        st_gpp  = phase_stats(ec_gpp_dly,  mgpp_dly)
        st_et   = phase_stats(ec_et_dly,   metv3_dly)
        for var, st in [("NDVI", st_ndvi), ("GPP", st_gpp), ("ET", st_et)]:
            stat_rows.append({"site": site, "indicator": var, **st})

        # ---- Row 1: NDVI ----
        ax = axes[0, j]
        shade_springs(ax, d0, d1)
        s2_x, s2_y = break_gaps(s2_d["date"], s2_d["NDVI_s2"],
                                 max_gap_days=20)
        ax.plot(s2_x, s2_y,
                "-o", color=C_S2, lw=1.0, ms=3.2, mec="white", mew=0.4,
                alpha=0.9, zorder=2,
                label="Sentinel-2 NDVI  ·  10 m  ·  ~5-day")
        ax.plot(ec_d["date"], smooth_break(ec_d["NDVI"]),
                color=C_EC, lw=1.7, zorder=3,
                label="EC tower NDVI  ·  ~100 m  ·  daily (7-d)")
        ax.set_ylabel("NDVI (–)", fontsize=10.5)
        ax.set_ylim(-0.05, 0.85)
        style_time_ax(ax, d0, d1)
        ax.set_title(SITE_LABEL[site], fontsize=12.0, pad=8)
        if j == 0:
            ax.legend(loc="upper left", frameon=False, fontsize=9.5)
        add_phase_annotation(ax, st_ndvi)
        ax.text(-0.085, 1.07, panel_labels[0][j], transform=ax.transAxes,
                fontsize=14, fontweight="bold")

        # ---- Row 2: GPP ----
        ax = axes[1, j]
        shade_springs(ax, d0, d1)
        ax.plot(ec_d["date"], ec_d["GPP_gC_m2_d"], color=C_EC_RAW,
                lw=0.6, alpha=0.9, zorder=1)
        m_d, m_v = break_gaps(mgpp_d["date"], mgpp_d["GPP_mgpp"],
                              max_gap_days=15)
        ax.plot(m_d, m_v,
                "-o", color=C_MGPP, lw=1.2, ms=3.6, mec="white", mew=0.5,
                alpha=0.95, zorder=2,
                label="Meteosat MGPP  ·  0.05° (~5 km)  ·  10-day")
        ax.plot(ec_d["date"], smooth_break(ec_d["GPP_gC_m2_d"]),
                color=C_EC, lw=1.7, zorder=3,
                label="EC tower GPP  ·  ~100 m  ·  daily (7-d)")
        ax.set_ylabel(r"GPP (gC m$^{-2}$ d$^{-1}$)", fontsize=10.5)
        ymax = max(16, np.nanpercentile(ec_d["GPP_gC_m2_d"], 99) + 2)
        ax.set_ylim(-1, ymax)
        style_time_ax(ax, d0, d1)
        if j == 0:
            ax.legend(loc="upper left", frameon=False, fontsize=9.5)
        add_phase_annotation(ax, st_gpp)
        ax.text(-0.085, 1.07, panel_labels[1][j], transform=ax.transAxes,
                fontsize=14, fontweight="bold")

        # ---- Row 3: ET (uses independent Meteosat source!) ----
        ax = axes[2, j]
        shade_springs(ax, d0, d1)
        ax.plot(metv3_d["date"], smooth_break(metv3_d["ET_mm"]),
                color=C_METV3, lw=1.4, alpha=0.95, zorder=2,
                label="Meteosat ETv3  ·  0.05° (~5 km)  ·  daily (7-d)")
        ax.plot(ec_d["date"], smooth_break(ec_d["ET_mm"]),
                color=C_EC, lw=1.7, zorder=3,
                label="EC tower ET  ·  ~100 m  ·  daily (7-d)")
        ax.set_ylabel("ET (mm day$^{-1}$)", fontsize=10.5)
        ax.set_ylim(-0.3, 8.5)
        style_time_ax(ax, d0, d1)
        if j == 0:
            ax.legend(loc="upper left", frameon=False, fontsize=9.5)
        add_phase_annotation(ax, st_et)
        ax.text(-0.085, 1.07, panel_labels[2][j], transform=ax.transAxes,
                fontsize=14, fontweight="bold")

    fig.text(0.5, 0.005,
             "Shaded band: Mar–Jun (Mediterranean growing season). "
             "$\\Delta$peak = mean(DOY$_{max}$(Sat) – DOY$_{max}$(EC)) over "
             "years (positive => Sat later). lag = argmax$_\\tau$ "
             "Pearson(EC[t], Sat[t+$\\tau$]), sign-flipped (positive => Sat "
             "later). r = correlation at the best lag. * = lag at the ±120 d "
             "scan boundary (no real peak inside window — interpret with "
             "care).",
             ha="center", fontsize=8.0, style="italic", color="#555555")

    out_png = OUT_DIR / "overview.png"
    out_pdf = OUT_DIR / "overview.pdf"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")

    stats_df = pd.DataFrame([
        {k: v for k, v in r.items() if k != "dpeak_per_year"}
        for r in stat_rows
    ])
    stats_csv = REPO / "data" / "phase_stats.csv"
    stats_df.to_csv(stats_csv, index=False, float_format="%.3f")
    print(f"wrote {stats_csv}")
    print(stats_df.to_string(index=False))

    # Supplementary figure: per-year peak-DOY differences
    _write_peak_doy_figure(stat_rows)


def _write_peak_doy_figure(stat_rows):
    """6 sub-panels (site x indicator) of per-year Δpeak DOY."""
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 5.0),
                              sharey=True,
                              gridspec_kw={"hspace": 0.35, "wspace": 0.12})
    layout = {
        ("Oran", "NDVI"): (0, 0), ("Oran", "GPP"): (0, 1),
        ("Oran", "ET"):   (0, 2),
        ("TzM",  "NDVI"): (1, 0), ("TzM",  "GPP"): (1, 1),
        ("TzM",  "ET"):   (1, 2),
    }
    ind_color = {"NDVI": C_S2, "GPP": C_MGPP, "ET": C_METV3}

    for row in stat_rows:
        i, j = layout[(row["site"], row["indicator"])]
        ax = axes[i, j]
        diffs = row["dpeak_per_year"]
        if diffs:
            xs = np.arange(1, len(diffs) + 1)
            ax.bar(xs, diffs, color=ind_color[row["indicator"]],
                   edgecolor="black", lw=0.6, alpha=0.85)
            mean = float(np.mean(diffs))
            ax.axhline(mean, color="black", lw=1.0, ls="--")
            ax.set_xticks(xs)
        ax.axhline(0, color="#555555", lw=0.7)
        ax.set_title(f"{row['site']} · {row['indicator']}", fontsize=10.5)
        if j == 0:
            ax.set_ylabel("Δpeak DOY (Sat – EC, d)", fontsize=9.5)
        if i == 1:
            ax.set_xlabel("year index", fontsize=9.5)
        ax.grid(axis="y", color="#dddddd", lw=0.6)

    ymax = max(60, max((abs(d) for r in stat_rows for d in r["dpeak_per_year"]),
                       default=10) * 1.15)
    for ax in axes.flat:
        ax.set_ylim(-ymax, ymax)

    fig.suptitle("Per-year peak DOY differences (Sat – EC). "
                  "Dashed line = multi-year mean.",
                  fontsize=11, y=1.02)
    out_png = OUT_DIR / "phase_peak_doy.png"
    out_pdf = OUT_DIR / "phase_peak_doy.pdf"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()
