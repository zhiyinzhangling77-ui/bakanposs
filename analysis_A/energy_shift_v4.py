"""
Eddy Covariance Flux Analysis - v4
====================================
v3からの変更点と根拠：

【外れ値判定の改善】
  - MAD法を昼間・夜間に分けて適用（昼夜混合では昼間の正常な高値が誤除去される）
  - FLUXNETの調査（Falge et al., 2001）で除去率9〜65%、平均35%は許容範囲と報告
    → 30〜40%は異常ではないが、昼夜分離でより精度の高い判定が可能
  - MAD閾値を5.0→7.0に緩和（夜間の低値分布が昼間の閾値を不当に下げる問題への対処）
  - 時間帯別zスコアのσを3.5→4.0に緩和（同理由）
  - 物理的絶対限界（第一段階）を残した上で統計フィルター（第二段階）を適用

【生育期定義の修正（文献根拠あり）】
  Oran（穀物・冬作）:
    根拠: Gomez-Gimenez et al. (2024, ADS), Russell (1990) via ScienceDirect
    播種: 10〜12月 → 生育開始: 11〜1月 → 収穫: 6月
    → GROWING_MONTHS = [11,12,1,2,3,4,5,6]  (従来の5月終了を6月まで延長)

  Tarazona（アーモンド・木本）:
    根拠: Agro Invest Spain (2022), Copernicus (2024)
    開花: 1〜2月 → 葉面積最大・蒸散旺盛: 4〜9月 → 収穫: 8〜10月 → 落葉: 11〜12月
    → GROWING_MONTHS = [1,2,3,4,5,6,7,8,9,10]  (葉展開〜収穫後まで)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy import stats

# ========= CONFIG =========
ORAN_PATH  = "/home/shion-nagamine/Dataset/Eddy data in Spain/Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.xlsx"
TARA_PATH  = "/home/shion-nagamine/Dataset/Eddy data in Spain/EddyAlmond_Raw5years_withG.xlsx"
OUT_DIR    = "./"

DAYTIME_START  = 9
DAYTIME_END    = 17
MIN_LE_DAYTIME = 10  # W/m²

# ---- 生育期定義（サイト別・文献根拠あり） ----
# Oran: スペイン冬作穀物 播種11-12月、収穫6月
#   出典: Gomez-Gimenez et al. 2024 (Sentinel-2 phenometrics, Andalusia)
#         Russell 1990 via Cantero-Martinez et al. 1995 (NE Spain barley)
ORAN_GROWING  = [11, 12, 1, 2, 3, 4, 5, 6]   # 旧: [11,12,1,2,3,4,5]

# Tarazona: アーモンド 開花1-2月、収穫8-10月、落葉11-12月
#   出典: Agro Invest Spain 2022 (harvest Aug-Oct, 6-8 month maturation)
#         Copernicus 2024 (winter flowering, frost risk)
TARA_GROWING  = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]   # 旧: [11,12,1,2,3,4,5]

# ---- 外れ値判定パラメータ ----
# 物理的絶対限界（第一段階）- 機器故障・伝送エラー由来の明らかな異常値のみ除去
# これを超える値はどの気候条件でも起こりえない
HARD_LE_MIN, HARD_LE_MAX = -100, 900
HARD_H_MIN,  HARD_H_MAX  = -300, 1000
HARD_G_MIN,  HARD_G_MAX  = -300, 600

# 統計フィルター（第二段階）
# 閾値を緩和: 昼夜混合MADでは昼間の正常高値が誤除去される問題を昼夜分離で解消
# 根拠: FLUXNET平均除去率35% (Falge et al. 2001) → 30-40%は許容範囲
MAD_THRESHOLD  = 7.0   # 旧: 5.0 → 緩和
ZSCORE_HOUR    = 4.0   # 旧: 3.5 → 緩和

DTYPE_COLORS = {
    'normal':   '#2196F3',
    'soil dry': '#FF9800',
    'atm dry':  '#4CAF50',
    'compound': '#F44336',
}


# =========================================================
# 1. LOAD & CLEAN
# =========================================================

def clean_numeric(df, cols):
    df = df.replace(['#N/A', 'NAN', 'nan', 'NaN', ' ', ''], np.nan)
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df


def hard_clip(df, site_name):
    """第一段階: 物理的絶対限界による除去（機器故障・伝送エラーのみ対象）"""
    n_before = len(df)
    mask = (
        df['LE'].between(HARD_LE_MIN, HARD_LE_MAX) &
        df['H'].between(HARD_H_MIN,  HARD_H_MAX)  &
        df['G'].between(HARD_G_MIN,  HARD_G_MAX)
    )
    n_hard = (~mask).sum()
    df = df[mask].copy()
    print(f"  [{site_name}] 物理限界除去: {n_hard}行 ({100*n_hard/n_before:.1f}%)")
    return df


def mad_filter_diurnal(series, threshold=MAD_THRESHOLD):
    """昼夜分離MAD法: Trueが正常値。"""
    ok = pd.Series(True, index=series.index)
    for is_day in [True, False]:
        idx = series.index[is_day_mask.loc[series.index] == is_day]
        s = series.loc[idx].dropna()
        if len(s) < 20:
            continue
        med = s.median()
        mad = (s - med).abs().median()
        if mad == 0:
            continue
        z = (s - med).abs() / (mad * 1.4826)
        ok.loc[idx[z >= threshold]] = False
    return ok


def hourly_zscore_filter(df, col, n_sigma=ZSCORE_HOUR):
    """時間帯別±n_sigmaフィルター: Trueが正常値。"""
    mask = pd.Series(True, index=df.index)
    for h in range(24):
        idx = df.index[df['hour'] == h]
        if len(idx) < 10:
            continue
        s = df.loc[idx, col]
        m, sd = s.mean(), s.std()
        if sd == 0 or np.isnan(sd):
            continue
        bad = (s < m - n_sigma * sd) | (s > m + n_sigma * sd)
        mask.loc[idx[bad]] = False
    return mask


def clip_outliers_v4(df, site_name):
    """
    第一段階: 物理的絶対限界（機器故障のみ）
    第二段階: 昼夜分離MAD + 時間帯別zスコア
    根拠: FLUXNET平均除去率35% (Falge et al., 2001, ScienceDirect Topics)
    """
    global is_day_mask

    # 第一段階
    df = hard_clip(df, site_name)

    # 昼夜マスクを作成（第二段階で使用）
    is_day_mask = (df['hour'] >= DAYTIME_START) & (df['hour'] <= DAYTIME_END)

    n_before = len(df)
    for col in ['LE', 'H', 'G']:
        if col not in df.columns:
            continue
        ok_mad  = mad_filter_diurnal(df[col])
        ok_hour = hourly_zscore_filter(df, col)
        df = df[ok_mad & ok_hour].copy()
        # is_day_maskをdfに合わせて更新
        is_day_mask = (df['hour'] >= DAYTIME_START) & (df['hour'] <= DAYTIME_END)

    n_after = len(df)
    pct = 100 * (n_before - n_after) / n_before
    print(f"  [{site_name}] 統計フィルター除去: {n_before}→{n_after}行 ({pct:.1f}%) ← FLUXNET平均35%と比較")
    return df


# --- Oran ---
oran = pd.read_excel(ORAN_PATH)
oran = clean_numeric(oran, ['SWC_1_1_1','SWC_2_1_1','SWC_3_1_1',
                             'H','LE','G','ET','VPD','FC_mass'])
oran['TIMESTAMP'] = pd.to_datetime(oran['TIMESTAMP'])
oran = oran.assign(
    date       = oran['TIMESTAMP'].dt.date,
    hour       = oran['TIMESTAMP'].dt.hour,
    month      = oran['TIMESTAMP'].dt.month,
    time_str   = oran['TIMESTAMP'].dt.strftime('%H:%M'),
    SWC_mean   = oran[['SWC_1_1_1','SWC_2_1_1','SWC_3_1_1']].mean(axis=1),
    GPP_proxy  = -oran['FC_mass'],
    # ← 修正: 6月を生育期に追加
    is_growing = oran['TIMESTAMP'].dt.month.isin(ORAN_GROWING),
)
oran = oran[['TIMESTAMP','date','hour','month','time_str','is_growing',
             'H','LE','G','ET','VPD','SWC_mean','GPP_proxy']]
oran['site'] = 'Oran'
print(f"\nOran 生育期月: {sorted(ORAN_GROWING)}")
print(f"  根拠: スペイン穀物播種11-12月、収穫6月 (Gomez-Gimenez et al. 2024; Russell 1990)")
print(f"Oran 読込行数: {len(oran)}")
oran = clip_outliers_v4(oran, 'Oran')

# --- Tarazona ---
tara = pd.read_excel(TARA_PATH)
tara = clean_numeric(tara, ['H','LE','G','ET','VPD','VPD_kPa',
                             'co2_flux','SWC_avg','SWC_1_1_1','SWC_2_1_1'])
tara['TIMESTAMP'] = pd.to_datetime(tara['TIMESTAMP'])
swc = (tara['SWC_avg'] if 'SWC_avg' in tara.columns
       else tara[['SWC_1_1_1','SWC_2_1_1']].mean(axis=1))
vpd = tara['VPD'] if 'VPD' in tara.columns else tara['VPD_kPa']
tara = tara.assign(
    date       = tara['TIMESTAMP'].dt.date,
    hour       = tara['TIMESTAMP'].dt.hour,
    month      = tara['TIMESTAMP'].dt.month,
    time_str   = tara['TIMESTAMP'].dt.strftime('%H:%M'),
    SWC_mean   = swc,
    VPD        = vpd,
    GPP_proxy  = -tara['co2_flux'],
    # ← 修正: 1-10月を生育期に
    is_growing = tara['TIMESTAMP'].dt.month.isin(TARA_GROWING),
)
tara = tara[['TIMESTAMP','date','hour','month','time_str','is_growing',
             'H','LE','G','ET','VPD','SWC_mean','GPP_proxy']]
tara['site'] = 'Tarazona'
print(f"\nTarazona 生育期月: {sorted(TARA_GROWING)}")
print(f"  根拠: アーモンド開花1-2月、収穫8-10月 (Agro Invest Spain 2022; Copernicus 2024)")
print(f"Tarazona 読込行数: {len(tara)}")
tara = clip_outliers_v4(tara, 'Tarazona')

df = pd.concat([oran, tara], ignore_index=True)


# =========================================================
# 2. ENERGY BALANCE CLOSURE → サイト別補正係数
# =========================================================

closure_slopes = {}

for site in df['site'].unique():
    sub = df[(df['site'] == site) &
             (df['hour'] >= DAYTIME_START) &
             (df['hour'] <= DAYTIME_END)].dropna(subset=['H','LE','G'])
    rn_proxy = sub['H'] + sub['LE'] + sub['G']
    hle      = sub['H'] + sub['LE']
    mask     = rn_proxy.notna() & hle.notna()
    slope, intercept, r, p, se = stats.linregress(rn_proxy[mask], hle[mask])
    closure_slopes[site] = slope
    print(f"\n[{site}] エネルギー収支閉合 OLS slope = {slope:.3f}  R²={r**2:.3f}")

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(rn_proxy[mask], hle[mask], alpha=0.1, s=4, color='steelblue')
    lim = [rn_proxy[mask].quantile(0.01), rn_proxy[mask].quantile(0.99)]
    ax.plot(lim, lim, 'r--', lw=1.5, label='1:1 line')
    xs = np.linspace(lim[0], lim[1], 200)
    ax.plot(xs, slope * xs + intercept, 'k-', lw=1.2,
            label=f'OLS slope={slope:.2f}  R²={r**2:.2f}')
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_title(f'{site} — Energy Balance Closure (daytime)', fontsize=12)
    ax.set_xlabel('Rn − G  [W m⁻²]'); ax.set_ylabel('H + LE  [W m⁻²]')
    ax.legend(); ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(f"{OUT_DIR}{site}_energy_closure_v4.png", dpi=300)
    plt.close()
    print(f"  Saved: {site}_energy_closure_v4.png")


def apply_closure_correction(df, slopes):
    df = df.copy()
    for site, slope in slopes.items():
        idx = df['site'] == site
        df.loc[idx, 'LE_corr'] = df.loc[idx, 'LE'] / slope
        df.loc[idx, 'H_corr']  = df.loc[idx, 'H']  / slope
    return df

df = apply_closure_correction(df, closure_slopes)


# =========================================================
# 3. DAILY AGGREGATION
# =========================================================

daily = df.groupby(['site','date','is_growing']).agg(
    H         = ('H',         'mean'),
    LE        = ('LE',        'mean'),
    LE_corr   = ('LE_corr',   'mean'),
    H_corr    = ('H_corr',    'mean'),
    G         = ('G',         'mean'),
    ET        = ('ET',        'sum'),
    VPD       = ('VPD',       'mean'),
    SWC_mean  = ('SWC_mean',  'mean'),
    GPP_proxy = ('GPP_proxy', 'mean'),
).reset_index()

daily['EF_corr']    = daily['LE_corr'] / (daily['H_corr'] + daily['LE_corr']).replace(0, np.nan)
daily['Bowen_corr'] = daily['H_corr']  / daily['LE_corr'].replace(0, np.nan)


# =========================================================
# 4. 干ばつ分類（サイト×季節別閾値）
# =========================================================

thresholds = daily.groupby(['site','is_growing']).apply(
    lambda g: pd.Series({
        'swc_th': g['SWC_mean'].quantile(0.4),
        'vpd_th': g['VPD'].quantile(0.6),
    })
).reset_index()

print("\n干ばつ閾値（サイト×季節別）:")
print(thresholds.to_string(index=False))

daily = daily.merge(thresholds, on=['site','is_growing'])

def classify(row):
    if row['SWC_mean'] < row['swc_th'] and row['VPD'] > row['vpd_th']:
        return 'compound'
    elif row['SWC_mean'] < row['swc_th']:
        return 'soil dry'
    elif row['VPD'] > row['vpd_th']:
        return 'atm dry'
    else:
        return 'normal'

daily['drought_type'] = daily.apply(classify, axis=1)

print("\n日数集計（サイト×季節×干ばつタイプ）:")
ct = daily.groupby(['site','is_growing','drought_type']).size().unstack(fill_value=0)
print(ct.to_string())

# n<20の場合に警告
for (site, is_growing, dtype), n in daily.groupby(['site','is_growing','drought_type']).size().items():
    if n < 20:
        label = '生育期' if is_growing else '非生育期'
        print(f"  ⚠ {site} / {label} / {dtype}: n={n}日 → 解釈注意")


# =========================================================
# 5. HELPER: 半時間別日変化
# =========================================================

SEASON_LABELS = {True: 'Growing Season', False: 'Non-Growing Season'}

def get_diurnal_sub(site, dtype, is_growing, daytime_only=False,
                    min_le=None, exclude_neg_H=False):
    cond = ((daily['site'] == site) &
            (daily['drought_type'] == dtype) &
            (daily['is_growing'] == is_growing))
    days = daily[cond]['date']
    sub  = df[(df['site'] == site) & (df['date'].isin(days))].copy()
    if daytime_only:
        sub = sub[(sub['hour'] >= DAYTIME_START) & (sub['hour'] <= DAYTIME_END)]
    if min_le is not None:
        sub = sub[sub['LE'] > min_le]
    if exclude_neg_H:
        sub = sub[sub['H_corr'] >= 0]
    return sub, len(days)


def make_diurnal(sub, col, min_n=50):
    if len(sub) < min_n:
        return None
    d = sub.groupby('time_str')[[col]].mean()
    d.index = pd.to_datetime(d.index, format='%H:%M')
    d = d.sort_index()
    d[col] = d[col].rolling(3, center=True).mean()
    return d[col]


# =========================================================
# 6〜8. プロット共通ヘルパー
# =========================================================
# 一覧PNG用にデータをキャッシュする辞書
# key: (plot_type, site, tag)  value: list of (x, y, color, label)
_plot_cache = {}

def _draw_diurnal_ax(ax, lines_data, title, ylabel, hour_interval=2,
                     hline=None, ylim=None):
    """axes に折れ線を描画する共通処理。"""
    for x, y, color, label in lines_data:
        ax.plot(x, y, color=color, label=label, linewidth=1.8)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=hour_interval))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')
    if hline is not None:
        ax.axhline(hline['y'], color='gray', linestyle='--',
                   lw=0.8, label=hline['label'])
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel('Time of Day', fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    ax.tick_params(labelsize=7)


# =========================================================
# 6. LE日変化
# =========================================================

for site in df['site'].unique():
    for is_growing, season_label in SEASON_LABELS.items():
        tag = 'growing' if is_growing else 'nongrowing'
        season_str = (f'Growing (Oran:Nov–Jun / Tara:Jan–Oct)'
                      if is_growing else
                      f'Non-Growing (Oran:Jul–Oct / Tara:Nov–Dec)')
        lines_data = []
        for dtype in ['normal','soil dry','atm dry','compound']:
            sub, n_days = get_diurnal_sub(site, dtype, is_growing)
            s = make_diurnal(sub, 'LE')
            if s is None:
                continue
            lines_data.append((s.index, s.values,
                                DTYPE_COLORS[dtype], f'{dtype} (n={n_days}d)'))
        if not lines_data:
            continue

        # 個別PNG
        fig, ax = plt.subplots(figsize=(11, 5))
        _draw_diurnal_ax(ax, lines_data,
                         title=f'{site} — LE Diurnal  [{season_str}]',
                         ylabel='LE (W m⁻²)', hour_interval=2)
        fig.tight_layout()
        fig.savefig(f"{OUT_DIR}{site}_LE_diurnal_{tag}_v4.png", dpi=300)
        plt.close()
        print(f"Saved: {site}_LE_diurnal_{tag}_v4.png")

        # キャッシュ保存
        _plot_cache[('LE', site, tag)] = dict(
            lines=lines_data,
            title=f'{site}\nLE [{season_label}]',
            ylabel='LE (W m⁻²)', hour_interval=2)


# =========================================================
# 7. EF日変化（閉合補正後、昼間のみ）
# =========================================================

for site in df['site'].unique():
    slope = closure_slopes[site]
    for is_growing, season_label in SEASON_LABELS.items():
        tag = 'growing' if is_growing else 'nongrowing'
        lines_data = []
        for dtype in ['normal','soil dry','atm dry','compound']:
            sub, n_days = get_diurnal_sub(site, dtype, is_growing,
                                          daytime_only=True, min_le=MIN_LE_DAYTIME,
                                          exclude_neg_H=(site=='Tarazona'))
            if len(sub) < 50:
                continue
            denom = (sub['LE_corr'] + sub['H_corr']).replace(0, np.nan)
            sub = sub.copy()
            sub['EF_corr'] = sub['LE_corr'] / denom
            sub = sub[(sub['EF_corr'] >= 0) & (sub['EF_corr'] <= 1.3)]
            s = make_diurnal(sub, 'EF_corr')
            if s is None:
                continue
            lines_data.append((s.index, s.values,
                                DTYPE_COLORS[dtype], f'{dtype} (n={n_days}d)'))
        if not lines_data:
            continue

        # 個別PNG
        fig, ax = plt.subplots(figsize=(11, 5))
        _draw_diurnal_ax(ax, lines_data,
                         title=f'{site} — EF corr. (slope={slope:.2f})  [{season_label}]',
                         ylabel='EF = LE_corr/(LE_corr+H_corr)',
                         hour_interval=1,
                         hline=dict(y=0.5, label='EF=0.5'),
                         ylim=(0, 1.2))
        fig.tight_layout()
        fig.savefig(f"{OUT_DIR}{site}_EF_corr_{tag}_v4.png", dpi=300)
        plt.close()
        print(f"Saved: {site}_EF_corr_{tag}_v4.png")

        _plot_cache[('EF', site, tag)] = dict(
            lines=lines_data,
            title=f'{site}\nEF corr. [{season_label}]',
            ylabel='EF',
            hour_interval=1,
            hline=dict(y=0.5, label='EF=0.5'),
            ylim=(0, 1.2))


# =========================================================
# 8. Bowen比日変化（昼間 + H>=0）
# =========================================================

for site in df['site'].unique():
    for is_growing, season_label in SEASON_LABELS.items():
        tag = 'growing' if is_growing else 'nongrowing'
        lines_data = []
        for dtype in ['normal','soil dry','atm dry','compound']:
            sub, n_days = get_diurnal_sub(site, dtype, is_growing,
                                          daytime_only=True, min_le=MIN_LE_DAYTIME,
                                          exclude_neg_H=True)
            if len(sub) < 50:
                continue
            sub = sub.copy()
            sub['Bowen_corr'] = sub['H_corr'] / sub['LE_corr'].replace(0, np.nan)
            sub = sub[(sub['Bowen_corr'] > -1) & (sub['Bowen_corr'] < 10)]
            s = make_diurnal(sub, 'Bowen_corr')
            if s is None:
                continue
            lines_data.append((s.index, s.values,
                                DTYPE_COLORS[dtype], f'{dtype} (n={n_days}d)'))
        if not lines_data:
            continue

        # 個別PNG
        fig, ax = plt.subplots(figsize=(11, 5))
        _draw_diurnal_ax(ax, lines_data,
                         title=f'{site} — Bowen Ratio corr. (H≥0)  [{season_label}]',
                         ylabel='Bowen = H_corr/LE_corr',
                         hour_interval=1,
                         hline=dict(y=1.0, label='Bowen=1 (H=LE)'))
        fig.tight_layout()
        fig.savefig(f"{OUT_DIR}{site}_Bowen_corr_{tag}_v4.png", dpi=300)
        plt.close()
        print(f"Saved: {site}_Bowen_corr_{tag}_v4.png")

        _plot_cache[('Bowen', site, tag)] = dict(
            lines=lines_data,
            title=f'{site}\nBowen [{season_label}]',
            ylabel='Bowen',
            hour_interval=1,
            hline=dict(y=1.0, label='Bowen=1'))


# =========================================================
# 8b. エネルギー収支閉合図もキャッシュ（散布図）
# =========================================================
# closure散布図データを再計算してキャッシュ
_closure_cache = {}
for site in df['site'].unique():
    sub = df[(df['site'] == site) &
             (df['hour'] >= DAYTIME_START) &
             (df['hour'] <= DAYTIME_END)].dropna(subset=['H','LE','G'])
    rn  = (sub['H'] + sub['LE'] + sub['G']).values
    hle = (sub['H'] + sub['LE']).values
    slope = closure_slopes[site]
    _closure_cache[site] = dict(rn=rn, hle=hle, slope=slope)


# =========================================================
# 9. 全図一覧PNG
# =========================================================
# レイアウト:
#   行: plot_type × season  (LE_growing, LE_nongrowing,
#                             EF_growing, EF_nongrowing,
#                             Bowen_growing, Bowen_nongrowing,
#                             Closure)
#   列: Oran / Tarazona
#
# 7行×2列 = 14パネル

sites    = ['Oran', 'Tarazona']
ptypes   = ['LE', 'EF', 'Bowen']
seasons  = ['growing', 'nongrowing']

# 行リスト: (plot_type_or_closure, season_tag)
row_defs = ([(pt, s) for pt in ptypes for s in seasons]
            + [('Closure', None)])   # 最終行

n_rows = len(row_defs)   # 7
n_cols = len(sites)      # 2

fig_all, axes_all = plt.subplots(
    n_rows, n_cols,
    figsize=(22, n_rows * 3.8),
    constrained_layout=True
)

for r, (ptype, season_tag) in enumerate(row_defs):
    for c, site in enumerate(sites):
        ax = axes_all[r, c]

        if ptype == 'Closure':
            # 散布図
            cd = _closure_cache[site]
            rn, hle, slope = cd['rn'], cd['hle'], cd['slope']
            p1, p99 = np.percentile(rn, 1), np.percentile(rn, 99)
            ax.scatter(rn, hle, alpha=0.05, s=2, color='steelblue', rasterized=True)
            ax.plot([p1, p99], [p1, p99], 'r--', lw=1.2, label='1:1')
            xs = np.linspace(p1, p99, 200)
            # intercept再計算
            _, intercept, r2, _, _ = stats.linregress(rn, hle)
            ax.plot(xs, slope * xs + intercept, 'k-', lw=1.0,
                    label=f'slope={slope:.2f}')
            ax.set_xlim(p1, p99); ax.set_ylim(p1, p99)
            ax.set_title(f'{site}\nEnergy Closure', fontsize=9)
            ax.set_xlabel('Rn−G [W m⁻²]', fontsize=8)
            ax.set_ylabel('H+LE [W m⁻²]', fontsize=8)
            ax.legend(fontsize=7); ax.grid(alpha=0.3)
            ax.tick_params(labelsize=7)
        else:
            key = (ptype, site, season_tag)
            if key not in _plot_cache:
                ax.set_visible(False)
                continue
            cd = _plot_cache[key]
            _draw_diurnal_ax(
                ax,
                cd['lines'],
                title=cd['title'],
                ylabel=cd['ylabel'],
                hour_interval=cd.get('hour_interval', 2),
                hline=cd.get('hline'),
                ylim=cd.get('ylim'),
            )

# 列タイトル（上部）
for c, site in enumerate(sites):
    axes_all[0, c].set_title(
        f'【{site}】\n' + axes_all[0, c].get_title(),
        fontsize=10, fontweight='bold'
    )

fig_all.suptitle(
    'Eddy Covariance Flux Overview  —  Oran (Cereal) vs Tarazona (Almond)\n'
    'Closure-corrected | Outlier: hard-clip + diurnal MAD(×7) + hourly z-score(4σ)\n'
    'Growing season: Oran Nov–Jun (winter cereal) | Tarazona Jan–Oct (almond)',
    fontsize=11, y=1.01
)

overview_path = f"{OUT_DIR}overview_all_v4.png"
fig_all.savefig(overview_path, dpi=180, bbox_inches='tight')
plt.close(fig_all)
print(f"\nSaved: overview_all_v4.png  ← 全図一覧")


# =========================================================
# 9. サマリー統計
# =========================================================

summary = daily.groupby(['site','is_growing','drought_type']).agg(
    n_days       = ('date',        'count'),
    LE_corr_mean = ('LE_corr',     'mean'),
    H_corr_mean  = ('H_corr',      'mean'),
    EF_corr      = ('EF_corr',     'mean'),
    Bowen_corr   = ('Bowen_corr',  'mean'),
    VPD_mean     = ('VPD',         'mean'),
    SWC_mean     = ('SWC_mean',    'mean'),
    GPP_mean     = ('GPP_proxy',   'mean'),
).round(3).reset_index()

summary['season'] = summary['is_growing'].map(
    {True: 'growing', False: 'non-growing'})
summary = summary.drop(columns='is_growing')

summary.to_csv(f"{OUT_DIR}drought_summary_v4.csv", index=False)
print("\n=== Summary Table ===")
print(summary.to_string(index=False))
print(f"\nSaved: drought_summary_v4.csv")

# 除去率の根拠メモを出力
print("\n" + "="*60)
print("【外れ値除去率の根拠】")
print("  Falge et al. (2001) via ScienceDirect Topics:")
print("  FLUXNETサイト調査 → 除去・欠損率 9〜65%, 平均35%")
print("  → 本解析の除去率が30-40%であれば許容範囲内")
print("  昼夜分離MAD (threshold=7.0) + 時間帯別z-score (σ=4.0) を適用")
print("="*60)
print("✅ All done.")

# === v11連携用: 中間ファイル出力 ===
import json
daily.to_parquet(f"{OUT_DIR}daily_classified_v4.parquet")
with open(f"{OUT_DIR}closure_slopes_v4.json", "w") as f:
    json.dump(closure_slopes, f)
print(f"\n[v11連携] daily_classified_v4.parquet と closure_slopes_v4.json を保存")