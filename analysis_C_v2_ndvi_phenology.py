"""
analysis_C_v2_ndvi_phenology.py
=================================
解析A (v27) の結論を NDVI フェノロジーで補強・検証する整理版。

v1 (analysis_C_v1.py) は多目的に拡張されすぎたため、ANALYSIS_C_PLAN.md §4.2 に
従って 3 つの焦点に絞り直したのが本 v2:

  [Q1] NDVI で active period を客観定義し、解析A の月仮定
       (Oran=Nov–Jun cereal, Tarazona=Jun–Sep almond peak) を validate
  [Q2] NDVI peak ±30 d window 内の event だけで τ を再 fit し、
       解析A の τ ≈ 3.0–3.8 d を生理学的に追認
  [Q3] NDVI と EC LE の lag correlation で位相関係を確認

解析A v27 仕様 (ANALYSIS_A_FINAL.md §3.2) を完全踏襲:
  - LE(d) = LE_∞ + (LE_0 - LE_∞) · exp(-d/τ), 2-parameter fit
  - LE_∞ は d ≥ 10 d の median (event を超えた背景)
  - LE_0 bound は adaptive (Oran rainfed の低 LE で fit 不能になる F6 を回避)
  - τ ∈ [0.5, 60] d
  - bootstrap B=2000, validation filter (τ at boundary / R²<0.3 を除外)
  - MDE = 1.96·√(SE1² + SE2²) で pairwise power 明示

event 定義 (ANALYSIS_A_FINAL.md §3.2):
  - Oran: Rain > 3 mm/day
  - Tarazona: Irrig > 0.5 mm/day
  - window: d=0–14 d post-event

入力データ:
  EC daily : data_loaders 経由 (修正版ローダ)
  NDVI     : MOD13Q1 250m / 16-day AppEEARS CSV (両サイト同梱)
"""

import json
import sys
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

warnings.filterwarnings("ignore")

# v1 からは I/O ユーティリティのみ再利用 (load_ndvi_csv, smooth_ndvi, match_ndvi)
# 解析判定ロジックは v2 で全部書き直し
sys.path.insert(0, str(Path(__file__).parent))
from analysis_C_v1 import load_ndvi_csv, smooth_ndvi, match_ndvi, SITE_ID_IN_CSV
from analysis_A_v9 import SITES, GROWING_MONTHS, PATHS as A_PATHS
from data_loaders import (
    load_oran_ec_clean, load_tarazona_ec_clean, normalize_swc,
)


# ===========================================================================
# 0. 設定
# ===========================================================================
BASE_NC = Path("/mnt/hdd/Dataset")
NDVI_APPEEARS_CSV = BASE_NC / "MOD13Q1_NDVI_EVI/MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv"
SAVE_DIR = Path("./output_analysis_C_v2")
SAVE_DIR.mkdir(exist_ok=True)

# 解析A v27 final から (ANALYSIS_A_FINAL.md §3 / §2.1)
A_RESULTS = {
    "Oran_active": {"tau": 2.82, "ci": (1.82, 5.34), "se": 0.90, "n": 10,
                     "months": [11, 12, 1, 2, 3, 4, 5, 6]},
    "Tarazona_active": {"tau": 3.36, "ci": (2.46, 4.90), "se": 0.62, "n": 41,
                         "months": [6, 7, 8, 9]},
}

NDVI_THR_ACTIVE = 0.30           # plan §4.2 推奨値
EVENT_WINDOW_DAYS = 14           # 解析A v27 仕様
PEAK_WINDOW_DAYS = 30            # plan §4.2 (NDVI peak ±X)
TAU_BOUND = (0.5, 60.0)          # 解析A v27 仕様
BOOTSTRAP_N = 2000
RNG = np.random.default_rng(20260518)


# ===========================================================================
# 1. NDVI フェノロジー指標抽出 (SOS / Peak / EOS)
# ===========================================================================
def extract_phenology_metrics(ndvi_sg, ndvi_threshold=NDVI_THR_ACTIVE):
    """Savitzky-Golay 平滑化済 NDVI から年別に SOS/Peak/EOS を抽出.

    SOS (Start of Season): 春の上昇局面で NDVI が threshold を最初に越えた日
    Peak: その年の NDVI 最大日 (smoothed)
    EOS (End of Season):    SOS 以降に NDVI が最後に threshold を越えていた日
    """
    df = ndvi_sg.copy()
    df["year"] = pd.to_datetime(df["date"]).dt.year
    df = df.dropna(subset=["NDVI_sg"]).sort_values("date").reset_index(drop=True)

    out = []
    for year, g in df.groupby("year"):
        g = g.sort_values("date").reset_index(drop=True)
        if len(g) < 5:
            continue
        above = g["NDVI_sg"] >= ndvi_threshold
        if not above.any():
            continue
        first_above = above.idxmax()
        last_above  = above[::-1].idxmax()
        peak_i = g["NDVI_sg"].idxmax()
        out.append({
            "year": int(year),
            "sos_date":  g.loc[first_above, "date"],
            "peak_date": g.loc[peak_i, "date"],
            "eos_date":  g.loc[last_above, "date"],
            "ndvi_peak": float(g.loc[peak_i, "NDVI_sg"]),
            "n_above":   int(above.sum()),
            "season_days": int((g.loc[last_above, "date"]
                                 - g.loc[first_above, "date"]).days),
        })
    return pd.DataFrame(out)


def validate_active_period(pheno, assumed_months, site):
    """NDVI で active と判定された期間と、解析A 仮定の月集合の overlap.

    Returns dict with overlap_pct, n_days_ndvi, n_days_assumed, etc.
    """
    rows = []
    for _, r in pheno.iterrows():
        year = r["year"]
        sos, eos = pd.Timestamp(r["sos_date"]), pd.Timestamp(r["eos_date"])
        ndvi_days = pd.date_range(sos, eos, freq="D")
        # 解析A 仮定: 月集合に対応する日 (その年内)
        yr_days = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")
        assumed_days = yr_days[yr_days.month.isin(assumed_months)]
        inter = ndvi_days.intersection(assumed_days)
        union = ndvi_days.union(assumed_days)
        rows.append({
            "year": year,
            "ndvi_active_days": len(ndvi_days),
            "assumed_active_days": len(assumed_days),
            "intersection_days": len(inter),
            "iou": len(inter) / max(1, len(union)),
            "ndvi_in_assumed_pct": 100.0 * len(inter) / max(1, len(ndvi_days)),
            "assumed_in_ndvi_pct": 100.0 * len(inter) / max(1, len(assumed_days)),
        })
    return pd.DataFrame(rows)


# ===========================================================================
# 2. Event detection (解析A v27 仕様)
# ===========================================================================
def _oran_rain_from_raw(filepath):
    """Oran AmeriFlux 生 CSV から日次積算降水量を抽出.

    `load_oran_ec_clean` は降水を捨てているので、event 検出のためだけに
    raw を読み直す. AmeriFlux 命名規約で P, P_RAIN, PRECIP_TOT を試す.
    """
    df = pd.read_csv(filepath, low_memory=False)
    src_col = "TIMESTAMP" if "TIMESTAMP" in df.columns else "DateTime"
    dt = pd.to_datetime(df[src_col].astype(str).str.strip(),
                          format="mixed", errors="coerce")
    df = df.assign(datetime=dt).dropna(subset=["datetime"])
    rain_candidates = [c for c in df.columns
                        if c.upper() in ("P", "P_RAIN", "P_RAIN_1_1_1",
                                          "PRECIP", "PRECIP_TOT", "RAIN")]
    if not rain_candidates:
        return None
    col = rain_candidates[0]
    s = pd.to_numeric(df[col], errors="coerce")
    s = s.where(s > -9000, np.nan)  # sentinel
    df["_p"] = s
    df["date"] = df["datetime"].dt.normalize()
    daily = df.groupby("date", as_index=False)["_p"].sum(min_count=1)
    daily = daily.rename(columns={"_p": "Rain_mm"})
    return daily


def detect_events(ec_df, site, oran_raw_path=None):
    """rain (Oran) または irrigation (Tarazona) event の発生日を返す.

    解析A v27 仕様:
      - Oran     : Rain_mm > 3.0
      - Tarazona : Irrig_mm > 0.5
    Oran で日次 EC に Rain が無ければ、生 CSV から抽出してマージする.
    """
    df = ec_df.copy()
    if site == "Oran":
        col_candidates = ["Rain_mm", "Rain", "P", "Precip_mm"]
        thr = 3.0
    else:
        col_candidates = ["Irrig_mm", "Irrigation_mm", "Irrig"]
        thr = 0.5
    col = next((c for c in col_candidates if c in df.columns), None)

    if col is None and site == "Oran" and oran_raw_path is not None:
        try:
            rain_daily = _oran_rain_from_raw(oran_raw_path)
            if rain_daily is not None and not rain_daily.empty:
                df = df.merge(rain_daily, on="date", how="left")
                col = "Rain_mm"
                print(f"  [events Oran] rain merged from raw CSV: "
                        f"n_days_with_P={int(df[col].gt(0).sum())}")
        except Exception as e:
            print(f"  [events Oran] raw rain extract failed: {e}")

    if col is None:
        print(f"  [events {site}] event 列なし — Rain/Irrig 不在")
        return pd.DatetimeIndex([])
    df[col] = pd.to_numeric(df[col], errors="coerce")
    events = df.loc[df[col] > thr, "date"]
    return pd.DatetimeIndex(pd.to_datetime(events).sort_values().unique())


def build_event_pool(events, ec_df, le_col="LE", window=EVENT_WINDOW_DAYS):
    """各 event について d=0..window の LE 系列を pool.

    重複 (次 event が早く来る) は最初の event に帰属させる
    (per-event 短窓を保証するため次 event 直前で打ち切り).
    """
    df = ec_df[["date", le_col]].copy()
    df["date"] = pd.to_datetime(df["date"])
    df = df.dropna(subset=[le_col]).sort_values("date").reset_index(drop=True)
    pool = []
    for i, ev in enumerate(events):
        ev = pd.Timestamp(ev)
        next_ev = events[i + 1] if i + 1 < len(events) else None
        end = ev + pd.Timedelta(days=window)
        if next_ev is not None:
            end = min(end, pd.Timestamp(next_ev) - pd.Timedelta(days=1))
        seg = df[(df["date"] >= ev) & (df["date"] <= end)].copy()
        if len(seg) < 3:
            continue
        seg["d"] = (seg["date"] - ev).dt.days.astype(float)
        seg["event_id"] = i
        pool.append(seg)
    if not pool:
        return pd.DataFrame(columns=["date", le_col, "d", "event_id"])
    return pd.concat(pool, ignore_index=True)


# ===========================================================================
# 3. Recovery curve fit (解析A v27 仕様)
# ===========================================================================
def _recovery_model(d, le0, tau, le_inf):
    return le_inf + (le0 - le_inf) * np.exp(-d / tau)


def fit_recovery(d, le, le_inf=None,
                  le0_bound=None, tau_bound=TAU_BOUND):
    """2-parameter exp recovery fit. Returns dict (None if failed).

    le_inf が None なら d ≥ 10 d の median を使う (解析A v27 仕様).
    le0_bound が None なら observed range で adaptive (CLAUDE.md F6).
    """
    d = np.asarray(d, dtype=float)
    le = np.asarray(le, dtype=float)
    msk = np.isfinite(d) & np.isfinite(le)
    d, le = d[msk], le[msk]
    if len(d) < 3:
        return None

    if le_inf is None:
        tail = le[d >= 10]
        le_inf = float(np.median(tail)) if len(tail) >= 2 else float(np.percentile(le, 25))

    if le0_bound is None:
        # adaptive: observed minimum to 2× observed maximum
        lo = float(min(le_inf - 5.0, np.nanmin(le)))
        hi = float(max(le_inf + 5.0, 2.0 * np.nanmax(le)))
        le0_bound = (lo, hi)

    def f(d, le0, tau):
        return _recovery_model(d, le0, tau, le_inf)

    try:
        p0 = (float(np.nanmax(le)), 3.0)
        popt, pcov = curve_fit(
            f, d, le, p0=p0,
            bounds=([le0_bound[0], tau_bound[0]],
                    [le0_bound[1], tau_bound[1]]),
            maxfev=5000,
        )
    except Exception as e:
        return {"success": False, "err": str(e)}

    le0, tau = float(popt[0]), float(popt[1])
    pred = f(d, *popt)
    ss_res = float(np.sum((le - pred) ** 2))
    ss_tot = float(np.sum((le - le.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

    # validation filter (CLAUDE.md 規約 #5)
    at_boundary = (
        abs(tau - tau_bound[0]) < 1e-3 or abs(tau - tau_bound[1]) < 1e-3
        or abs(le0 - le0_bound[0]) < 1e-3 or abs(le0 - le0_bound[1]) < 1e-3
    )
    return {
        "success": True, "le0": le0, "tau": tau, "le_inf": le_inf,
        "r2": r2, "n": len(d), "at_boundary": bool(at_boundary),
    }


def bootstrap_tau_ci(pool, le_col="LE", B=BOOTSTRAP_N,
                      le_inf=None):
    """Event 単位 (block) bootstrap で τ の 95% CI.

    重要: 日次データを直接 resample すると同 event 内の系列相関を無視するので、
    event_id 単位で resample → 各 event 内は全日採用 (= block bootstrap).
    Validation filter で fit 失敗・τ at boundary・R²<0.3 を捨てる.
    """
    if pool.empty:
        return None
    if le_inf is None:
        tail = pool.loc[pool["d"] >= 10, le_col]
        le_inf = float(tail.median()) if len(tail) >= 2 else float(pool[le_col].quantile(0.25))

    event_ids = pool["event_id"].unique()
    taus, le0s, r2s = [], [], []
    n_excluded = 0
    for _ in range(B):
        sample_ids = RNG.choice(event_ids, size=len(event_ids), replace=True)
        sub = pd.concat([pool[pool["event_id"] == eid] for eid in sample_ids],
                         ignore_index=True)
        fit = fit_recovery(sub["d"].values, sub[le_col].values, le_inf=le_inf)
        if (fit is None or not fit.get("success")
                or fit.get("at_boundary") or fit.get("r2", 0) < 0.30):
            n_excluded += 1
            continue
        taus.append(fit["tau"])
        le0s.append(fit["le0"])
        r2s.append(fit["r2"])
    if len(taus) < 20:
        return {"valid": False, "n_valid": len(taus), "n_excluded": n_excluded}
    taus = np.array(taus)
    le0s = np.array(le0s)
    return {
        "valid": True,
        "tau_med":  float(np.median(taus)),
        "tau_ci":   (float(np.percentile(taus, 2.5)),
                     float(np.percentile(taus, 97.5))),
        "tau_se":   float(taus.std(ddof=1)),
        "le0_med":  float(np.median(le0s)),
        "le0_ci":   (float(np.percentile(le0s, 2.5)),
                     float(np.percentile(le0s, 97.5))),
        "r2_med":   float(np.median(r2s)),
        "n_valid":  int(len(taus)),
        "n_excluded": int(n_excluded),
    }


# ===========================================================================
# 4. NDVI peak window で τ 再計算
# ===========================================================================
def filter_events_by_ndvi_peak(events, pheno, window=PEAK_WINDOW_DAYS):
    """各 event を、その年の NDVI peak から ±window 日に入っていれば残す."""
    if pheno.empty:
        return pd.DatetimeIndex([])
    peak_by_year = {int(r["year"]): pd.Timestamp(r["peak_date"])
                     for _, r in pheno.iterrows()}
    keep = []
    for ev in events:
        ev = pd.Timestamp(ev)
        pk = peak_by_year.get(int(ev.year))
        if pk is None:
            continue
        if abs((ev - pk).days) <= window:
            keep.append(ev)
    return pd.DatetimeIndex(sorted(set(keep)))


# ===========================================================================
# 5. NDVI–LE lag correlation
# ===========================================================================
def ndvi_le_lag_correlation(merged, max_lag=20):
    """NDVI を ±max_lag 日 shift して LE との相関 r(lag).

    lag > 0: NDVI が LE に先行 (NDVI leads)
    lag < 0: NDVI が LE に遅行 (LE leads)
    """
    df = merged.dropna(subset=["NDVI", "LE"]).copy()
    df = df.sort_values("date").reset_index(drop=True)
    if len(df) < 30:
        return pd.DataFrame(columns=["lag", "r", "n"])
    rows = []
    for lag in range(-max_lag, max_lag + 1):
        if lag == 0:
            x, y = df["NDVI"], df["LE"]
        elif lag > 0:
            x = df["NDVI"].iloc[:-lag].reset_index(drop=True)
            y = df["LE"].iloc[lag:].reset_index(drop=True)
        else:
            x = df["NDVI"].iloc[-lag:].reset_index(drop=True)
            y = df["LE"].iloc[:lag].reset_index(drop=True)
        n = min(len(x), len(y))
        if n < 20:
            continue
        r, _ = stats.pearsonr(x[:n], y[:n])
        rows.append({"lag": lag, "r": float(r), "n": int(n)})
    return pd.DataFrame(rows)


# ===========================================================================
# 6. プロット
# ===========================================================================
def plot_phenology_overlap(per_site_pheno, validation, save_dir):
    """両サイト × 各年の NDVI active 期間 と 解析A 仮定 months の重ね描き."""
    fig, axes = plt.subplots(2, 1, figsize=(12, 5), sharex=False)
    for ax, site in zip(axes, ["Oran", "Tarazona"]):
        pheno = per_site_pheno[site]
        if pheno.empty:
            ax.text(0.5, 0.5, f"no phenology data for {site}",
                     ha="center", va="center", transform=ax.transAxes)
            continue
        assumed = GROWING_MONTHS[site]
        for _, r in pheno.iterrows():
            y = int(r["year"])
            sos = pd.Timestamp(r["sos_date"])
            eos = pd.Timestamp(r["eos_date"])
            pk  = pd.Timestamp(r["peak_date"])
            ax.barh(y, (eos - sos).days, left=sos.toordinal(),
                     color="green", alpha=0.55, height=0.6,
                     label="NDVI active" if _ == 0 else None)
            ax.plot(pk.toordinal(), y, "v", color="darkgreen",
                     markersize=8, label="NDVI peak" if _ == 0 else None)
        # 解析A 仮定 months: その年の仮定区間
        for y in sorted(pheno["year"].unique()):
            for m in assumed:
                start = pd.Timestamp(f"{y}-{m:02d}-01")
                if m == 12:
                    end = pd.Timestamp(f"{y + 1}-01-01") - pd.Timedelta(days=1)
                else:
                    end = pd.Timestamp(f"{y}-{m + 1:02d}-01") - pd.Timedelta(days=1)
                ax.barh(y + 0.35, (end - start).days, left=start.toordinal(),
                         color="grey", alpha=0.30, height=0.30)
        v = validation[site]
        if not v.empty:
            iou_mean = v["iou"].mean()
            ttl_suffix = f"(NDVI∩assumed/∪ = {iou_mean:.2f})"
        else:
            ttl_suffix = ""
        ax.set_title(f"{site}: NDVI-defined active (green) vs Analysis-A assumed months (grey) {ttl_suffix}",
                       fontsize=10)
        ax.set_ylabel("year")
        # X 軸を日付化
        years_present = sorted(pheno["year"].unique())
        xmin = pd.Timestamp(f"{min(years_present)}-01-01").toordinal()
        xmax = pd.Timestamp(f"{max(years_present)}-12-31").toordinal()
        ax.set_xlim(xmin, xmax)
        ticks = pd.date_range(f"{min(years_present)}-01-01",
                                f"{max(years_present)}-12-31",
                                freq="3MS")
        ax.set_xticks([t.toordinal() for t in ticks])
        ax.set_xticklabels([t.strftime("%Y-%m") for t in ticks],
                            rotation=45, ha="right", fontsize=8)
    fig.suptitle("Fig 02: Active-period validation — NDVI vs Analysis-A assumed months",
                  fontsize=11, y=0.995)
    fig.tight_layout()
    fp = save_dir / "fig02_active_period_overlap.png"
    fig.savefig(fp, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {fp}")


def plot_ndvi_seasonal(per_site_merged, per_site_pheno, save_dir):
    """両サイトの年次 NDVI 季節曲線 + SOS/Peak/EOS marker."""
    fig, axes = plt.subplots(2, 1, figsize=(12, 5), sharex=False)
    for ax, site in zip(axes, ["Oran", "Tarazona"]):
        m = per_site_merged[site]
        if m is None or m.empty or "NDVI" not in m.columns:
            ax.text(0.5, 0.5, "no NDVI", ha="center", va="center",
                     transform=ax.transAxes)
            continue
        m = m.dropna(subset=["NDVI"]).sort_values("date")
        ax.plot(m["date"], m["NDVI"], "o", ms=2, alpha=0.4, color="grey",
                 label="raw")
        if "NDVI_sg" in m.columns:
            sg = m.dropna(subset=["NDVI_sg"])
            ax.plot(sg["date"], sg["NDVI_sg"], "-", lw=1.0, color="green",
                     label="Savitzky-Golay")
        ax.axhline(NDVI_THR_ACTIVE, ls="--", color="orange", lw=0.8,
                    label=f"thr={NDVI_THR_ACTIVE}")
        pheno = per_site_pheno[site]
        for _, r in pheno.iterrows():
            ax.axvline(pd.Timestamp(r["peak_date"]), color="darkgreen",
                         ls=":", lw=0.7, alpha=0.8)
        ax.set_ylabel("NDVI")
        ax.set_title(f"{site}: NDVI seasonal cycle (peaks marked)",
                       fontsize=10)
        ax.legend(fontsize=8, loc="upper right")
        ax.set_ylim(-0.05, 0.9)
    axes[-1].set_xlabel("date")
    fig.suptitle("Fig 01: NDVI seasonal cycle (MOD13Q1)", fontsize=11, y=0.995)
    fig.tight_layout()
    fp = save_dir / "fig01_ndvi_seasonal_curve.png"
    fig.savefig(fp, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {fp}")


def plot_tau_at_peak(tau_results, save_dir):
    """解析A v27 結果と NDVI peak window τ を並べた forest plot."""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    rows = []
    # 解析A v27 (validated, from ANALYSIS_A_FINAL.md)
    for label, key in [("Oran (A v27, Nov–Jun pool)", "Oran_active"),
                          ("Tarazona (A v27, Jun–Sep pool)", "Tarazona_active")]:
        d = A_RESULTS[key]
        rows.append({"label": label, "tau": d["tau"],
                      "lo": d["ci"][0], "hi": d["ci"][1], "n": d["n"],
                      "source": "A_v27"})
    # NDVI peak window τ
    for site, res in tau_results.items():
        if res is None or not res.get("valid"):
            rows.append({"label": f"{site} (C v2, NDVI peak ±{PEAK_WINDOW_DAYS}d)",
                          "tau": np.nan, "lo": np.nan, "hi": np.nan,
                          "n": res.get("n_events", 0) if res else 0,
                          "source": "C_v2"})
            continue
        rows.append({"label": f"{site} (C v2, NDVI peak ±{PEAK_WINDOW_DAYS}d)",
                      "tau": res["tau_med"],
                      "lo": res["tau_ci"][0], "hi": res["tau_ci"][1],
                      "n": res.get("n_events", res["n_valid"]),
                      "source": "C_v2"})
    df = pd.DataFrame(rows)
    yy = np.arange(len(df))[::-1]
    colors = ["steelblue" if s == "A_v27" else "darkorange"
                for s in df["source"]]
    for y, (_, r), c in zip(yy, df.iterrows(), colors):
        if not np.isfinite(r["tau"]):
            ax.text(0.5, y, "fit failed (n insufficient or τ at boundary)",
                     va="center", fontsize=8, color="red")
            continue
        ax.plot([r["lo"], r["hi"]], [y, y], color=c, lw=2.5)
        ax.plot(r["tau"], y, "o", color=c, ms=8)
        ax.text(r["hi"] + 0.3, y, f"τ={r['tau']:.2f}d, n={r['n']}",
                  va="center", fontsize=8)
    ax.axvspan(3.0, 3.8, color="green", alpha=0.10,
                 label="A v27 universal band 3.0–3.8 d")
    ax.set_yticks(yy)
    ax.set_yticklabels(df["label"], fontsize=9)
    ax.set_xlabel("τ [days]  (95% bootstrap CI)")
    ax.set_xlim(0, max(df["hi"].max(skipna=True) + 2, 8))
    ax.set_title("Fig 03: τ at NDVI-peak window vs Analysis-A v27 (universal band shaded)",
                   fontsize=10)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fp = save_dir / "fig03_ndvi_peak_tau.png"
    fig.savefig(fp, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {fp}")


def plot_ndvi_le_lag(per_site_lag, save_dir):
    """NDVI–LE lag correlation の cross-correlogram."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, site in zip(axes, ["Oran", "Tarazona"]):
        d = per_site_lag[site]
        if d.empty:
            ax.text(0.5, 0.5, "insufficient data", ha="center", va="center",
                     transform=ax.transAxes)
            continue
        ax.plot(d["lag"], d["r"], "-o", ms=4, color="darkgreen", lw=1.2)
        # 最大 r とその lag
        i = d["r"].abs().idxmax()
        best_lag = int(d.loc[i, "lag"])
        best_r   = float(d.loc[i, "r"])
        ax.axvline(best_lag, color="red", ls="--", lw=0.8,
                    label=f"max |r| @ lag={best_lag}d (r={best_r:.2f})")
        ax.axvline(0, color="grey", lw=0.5)
        ax.axhline(0, color="grey", lw=0.5)
        ax.set_xlabel("lag (NDVI - LE) [days]\n(+: NDVI leads)")
        ax.set_title(f"{site}", fontsize=10)
        ax.legend(loc="best", fontsize=8)
    axes[0].set_ylabel("Pearson r(NDVI, LE)")
    fig.suptitle("Fig 04: NDVI ↔ LE lag correlation", fontsize=11)
    fig.tight_layout()
    fp = save_dir / "fig04_ndvi_le_lag.png"
    fig.savefig(fp, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {fp}")


# ===========================================================================
# 7. サイト個別ランナー
# ===========================================================================
def run_site_C2(site, ec_df, ndvi_path):
    """1 サイト分の処理: NDVI 読込 → smooth → match → phenology metrics."""
    print(f"\n{'='*60}\n  解析C v2 — {site}\n{'='*60}")

    ndvi_raw = load_ndvi_csv(ndvi_path, site_name=site,
                                site_id_in_csv=SITE_ID_IN_CSV[site],
                                reliability_max=1)
    ndvi_sg = smooth_ndvi(ndvi_raw, window=5, poly=2)
    # smooth_ndvi は date + NDVI_sg を返す
    # merge: ndvi_raw に NDVI_sg を結合してから match
    ndvi_full = ndvi_raw.merge(ndvi_sg, on="date", how="left")
    merged = match_ndvi(ec_df, ndvi_full, max_days=10)

    pheno = extract_phenology_metrics(ndvi_sg, ndvi_threshold=NDVI_THR_ACTIVE)
    print(f"  [{site}] phenology metrics extracted for "
            f"{len(pheno)} year(s)")
    if not pheno.empty:
        for _, r in pheno.iterrows():
            print(f"    year {int(r['year'])}: SOS={pd.Timestamp(r['sos_date']).date()}  "
                    f"Peak={pd.Timestamp(r['peak_date']).date()} (NDVI={r['ndvi_peak']:.2f})  "
                    f"EOS={pd.Timestamp(r['eos_date']).date()}  season={r['season_days']}d")
    return merged, ndvi_full, pheno


def run_tau_at_peak(site, ec_df, pheno, oran_raw_path=None):
    """site の EC データから event 検出 → NDVI peak window で filter → τ fit."""
    events_all = detect_events(ec_df, site, oran_raw_path=oran_raw_path)
    print(f"\n  [{site}] all events detected: n={len(events_all)}")
    if len(events_all) == 0:
        return None

    events_at_peak = filter_events_by_ndvi_peak(events_all, pheno,
                                                  window=PEAK_WINDOW_DAYS)
    print(f"  [{site}] events within NDVI peak ±{PEAK_WINDOW_DAYS}d: "
            f"n={len(events_at_peak)} / {len(events_all)}")

    pool = build_event_pool(events_at_peak, ec_df, le_col="LE",
                              window=EVENT_WINDOW_DAYS)
    print(f"  [{site}] pooled (event, day) obs: n={len(pool)}")
    if pool.empty or pool["event_id"].nunique() < 3:
        print(f"  [{site}] ⚠️ events too few for tau fit")
        return {"valid": False, "n_events": pool["event_id"].nunique()
                                                if not pool.empty else 0}

    # point estimate fit
    fit_point = fit_recovery(pool["d"].values, pool["LE"].values)
    if fit_point and fit_point.get("success"):
        print(f"  [{site}] point fit: τ={fit_point['tau']:.2f}d  "
                f"LE0={fit_point['le0']:.1f}  LE∞={fit_point['le_inf']:.1f}  "
                f"R²={fit_point['r2']:.2f}  "
                f"{'at-boundary' if fit_point['at_boundary'] else 'OK'}")

    # bootstrap CI
    boot = bootstrap_tau_ci(pool, le_col="LE", B=BOOTSTRAP_N,
                              le_inf=fit_point["le_inf"] if fit_point else None)
    if boot and boot.get("valid"):
        boot["n_events"] = int(pool["event_id"].nunique())
        print(f"  [{site}] bootstrap τ = {boot['tau_med']:.2f} "
                f"[{boot['tau_ci'][0]:.2f}, {boot['tau_ci'][1]:.2f}] d  "
                f"(n_valid={boot['n_valid']}/{BOOTSTRAP_N}, "
                f"R²_med={boot['r2_med']:.2f})")
    else:
        print(f"  [{site}] ⚠️ bootstrap failed: {boot}")
    return boot


# ===========================================================================
# 8. Main
# ===========================================================================
def main():
    print("=" * 60)
    print("解析C v2 — NDVI phenology reinforcement of Analysis A")
    print("=" * 60)

    # --- EC データロード (修正版ローダ) ---
    print("\n[STEP 1] EC データ読込")
    oran_ec = normalize_swc(load_oran_ec_clean(A_PATHS["oran_ec"]), "Oran")
    tara_ec = normalize_swc(load_tarazona_ec_clean(A_PATHS["tara_ec"]),
                              "Tarazona")
    # Oran は Rain_mm が無いことが多いので Tarazona 用 CSV から merge する手は無い
    # (別サイトなので) — Oran は event 列が無ければ event 検出スキップする設計

    # --- NDVI + phenology metrics ---
    print("\n[STEP 2] NDVI 読込 → smooth → phenology metrics")
    per_site_merged = {}
    per_site_pheno = {}
    for site, ec in [("Oran", oran_ec), ("Tarazona", tara_ec)]:
        try:
            merged, ndvi_full, pheno = run_site_C2(site, ec, NDVI_APPEEARS_CSV)
            per_site_merged[site] = merged
            per_site_pheno[site] = pheno
        except Exception as e:
            print(f"  ⚠️ {site}: {e}")
            per_site_merged[site] = pd.DataFrame()
            per_site_pheno[site] = pd.DataFrame()

    # --- 出力1: phenology metrics CSV ---
    all_pheno = []
    for site, p in per_site_pheno.items():
        if p.empty:
            continue
        pp = p.copy()
        pp["site"] = site
        all_pheno.append(pp)
    if all_pheno:
        out = pd.concat(all_pheno, ignore_index=True)
        for c in ["sos_date", "peak_date", "eos_date"]:
            out[c] = pd.to_datetime(out[c]).dt.strftime("%Y-%m-%d")
        fp = SAVE_DIR / "v2_phenology_metrics.csv"
        out.to_csv(fp, index=False)
        print(f"\n  saved: {fp}")

    # --- 出力2: active period validation CSV ---
    print("\n[STEP 3] NDVI active vs Analysis-A assumed months")
    validation = {}
    for site, pheno in per_site_pheno.items():
        if pheno.empty:
            validation[site] = pd.DataFrame()
            continue
        v = validate_active_period(pheno, GROWING_MONTHS[site], site)
        v["site"] = site
        validation[site] = v
        if not v.empty:
            print(f"  [{site}] IoU per year:")
            for _, r in v.iterrows():
                print(f"    {int(r['year'])}: IoU={r['iou']:.2f}  "
                        f"ndvi∩assumed/ndvi={r['ndvi_in_assumed_pct']:.0f}%  "
                        f"ndvi∩assumed/assumed={r['assumed_in_ndvi_pct']:.0f}%")
    val_df = pd.concat([v for v in validation.values() if not v.empty],
                         ignore_index=True) if validation else pd.DataFrame()
    if not val_df.empty:
        fp = SAVE_DIR / "v2_active_period_validation.csv"
        val_df.to_csv(fp, index=False)
        print(f"  saved: {fp}")

    # --- 出力3: NDVI peak window τ ---
    print("\n[STEP 4] τ at NDVI peak window (±{}d)".format(PEAK_WINDOW_DAYS))
    tau_results = {}
    for site, ec in [("Oran", oran_ec), ("Tarazona", tara_ec)]:
        pheno = per_site_pheno.get(site, pd.DataFrame())
        if pheno.empty:
            tau_results[site] = None
            continue
        raw = A_PATHS["oran_ec"] if site == "Oran" else None
        tau_results[site] = run_tau_at_peak(site, ec, pheno, oran_raw_path=raw)

    # CSV 出力
    rows = []
    for site, r in tau_results.items():
        a = A_RESULTS.get(f"{site}_active", {})
        if r and r.get("valid"):
            row = {
                "site": site,
                "n_events_peak_window": r.get("n_events"),
                "n_bootstrap_valid": r["n_valid"],
                "tau_C_v2": r["tau_med"],
                "tau_C_v2_lo": r["tau_ci"][0],
                "tau_C_v2_hi": r["tau_ci"][1],
                "tau_C_v2_se": r["tau_se"],
                "le0_C_v2": r["le0_med"],
                "r2_med": r["r2_med"],
                "tau_A_v27": a.get("tau"),
                "tau_A_v27_lo": a.get("ci", (None, None))[0],
                "tau_A_v27_hi": a.get("ci", (None, None))[1],
                "se_A_v27": a.get("se"),
            }
            # MDE for C_v2 vs A_v27
            if a.get("se") is not None:
                mde = 1.96 * np.sqrt(r["tau_se"] ** 2 + a["se"] ** 2)
                row["MDE"] = float(mde)
                row["obs_diff"] = float(abs(r["tau_med"] - a["tau"]))
                row["significant"] = row["obs_diff"] > mde
            rows.append(row)
        else:
            rows.append({"site": site, "n_events_peak_window":
                            r.get("n_events", 0) if r else 0,
                          "tau_C_v2": None, "fit_status": "failed"})
    if rows:
        fp = SAVE_DIR / "v2_ndvi_peak_tau.csv"
        pd.DataFrame(rows).to_csv(fp, index=False)
        print(f"\n  saved: {fp}")

    # --- 出力4: lag correlation ---
    print("\n[STEP 5] NDVI–LE lag correlation")
    per_site_lag = {}
    for site, m in per_site_merged.items():
        if m is None or m.empty or "NDVI" not in m.columns:
            per_site_lag[site] = pd.DataFrame()
            continue
        lag = ndvi_le_lag_correlation(m, max_lag=20)
        per_site_lag[site] = lag
        if not lag.empty:
            best = lag.loc[lag["r"].abs().idxmax()]
            print(f"  [{site}] max |r| = {best['r']:.2f} at lag={int(best['lag'])}d "
                    f"(NDVI leads by {int(best['lag'])}d if positive)")
    lag_concat = []
    for site, d in per_site_lag.items():
        if d.empty:
            continue
        dd = d.copy()
        dd["site"] = site
        lag_concat.append(dd)
    if lag_concat:
        fp = SAVE_DIR / "v2_ndvi_le_lag.csv"
        pd.concat(lag_concat, ignore_index=True).to_csv(fp, index=False)
        print(f"  saved: {fp}")

    # --- merged CSV (再現用) ---
    for site, m in per_site_merged.items():
        if m is None or m.empty:
            continue
        fp = SAVE_DIR / f"v2_{site}_merged.csv"
        m.to_csv(fp, index=False)
        print(f"  saved: {fp}")

    # --- 図 ---
    print("\n[STEP 6] 図を生成")
    plot_ndvi_seasonal(per_site_merged, per_site_pheno, SAVE_DIR)
    plot_phenology_overlap(per_site_pheno, validation, SAVE_DIR)
    plot_tau_at_peak(tau_results, SAVE_DIR)
    plot_ndvi_le_lag(per_site_lag, SAVE_DIR)

    # --- 最終サマリ ---
    print("\n" + "=" * 60)
    print("★ 解析C v2 サマリ")
    print("=" * 60)
    print("[Q1] NDVI active vs Analysis-A assumed months:")
    for site, v in validation.items():
        if v.empty:
            continue
        print(f"  {site}: mean IoU = {v['iou'].mean():.2f} "
                f"(n_years={len(v)})")
    print("\n[Q2] NDVI peak window τ (compared to Analysis-A v27):")
    for site, r in tau_results.items():
        a = A_RESULTS.get(f"{site}_active", {})
        if r and r.get("valid"):
            mde_str = ""
            if a.get("se") is not None:
                mde = 1.96 * np.sqrt(r["tau_se"] ** 2 + a["se"] ** 2)
                diff = abs(r["tau_med"] - a["tau"])
                sig = "SIG" if diff > mde else "NS (consistent)"
                mde_str = f"  |Δτ|={diff:.2f}d vs MDE={mde:.2f}d → {sig}"
            print(f"  {site}: τ_C_v2={r['tau_med']:.2f}d "
                    f"[{r['tau_ci'][0]:.2f}, {r['tau_ci'][1]:.2f}]  "
                    f"vs τ_A_v27={a.get('tau','?')}d{mde_str}")
        else:
            print(f"  {site}: fit failed")
    print("\n[Q3] NDVI–LE lag (max |r|):")
    for site, d in per_site_lag.items():
        if d.empty:
            continue
        best = d.loc[d["r"].abs().idxmax()]
        print(f"  {site}: r={best['r']:.2f} at lag={int(best['lag'])}d")

    print(f"\n[完了] 出力ディレクトリ: {SAVE_DIR}/")


if __name__ == "__main__":
    main()
