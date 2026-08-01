"""サイトレジストリ: サイトコード → データ位置 + 変数名マッピング。

多サイト拡張 (JP-Ta2 / JP-BBY / JP-Mse) の効き所。数値カーネルは変数名を
知らないので、サイトごとの実カラム名・ディレクトリ名の違いはここで吸収する。
FLUXNET2015 の基本名は共通だが ``TS_F_MDS_1`` / ``SWC_F_MDS_1`` の層インデックス
や、湿原・水田で欠ける変数はサイト依存になりうる。
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

from .config import RK_VARS


# R&K 表記 → JapanFlux2024 の既定カラム名。
# JapanFlux2024 は FLUXNET2015 とほぼ同名だが 3 点だけ規約が違う:
#   - 土壌系は層インデックス無し (TS_F_MDS / SWC_F_MDS)
#   - 炭素フラックスは可変 u* 閾値の基準値を `_vUT` と表記 (FLUXNET2015 の
#     `_VUT_REF` 相当、v は小文字・_REF 無し)。USTAR05/50/95 は分位変種なので不使用
#   - GER/GEP は daytime 分割法 (DT) で統一
# サイト固有に上書きしたい場合は SiteSpec.var_overrides で差し替える。
DEFAULT_VAR_MAP: dict[str, str] = {
    "Rg":  "SW_IN_F",
    "Ta":  "TA_F",
    "VPD": "VPD_F",
    "Ts":  "TS_F_MDS",
    "P":   "P_F",
    "th":  "SWC_F_MDS",
    "gH":  "H_F_MDS",
    "gLE": "LE_F_MDS",
    "GER": "RECO_DT_vUT",
    "NEE": "NEE_vUT",
    "GEP": "GPP_DT_vUT",
}


# --- A3 拡張: AmeriFlux/ICOS/KoFlux "BASE" 形式 (韓国・中国の 30 分値) -----------
# 位置修飾子 _H_V_R (水平_鉛直_反復, 例 `_1_1_1`) が付く FLUXNET 系の別規約。
# JapanFlux2024 との実務上の差は 3 点:
#   1. VPD 列が無い → 気温 TA と相対湿度 RH から導出する (下の @VPD マーカ)。
#   2. gap-fill 済みフラックスは `_PI_` 表記 (H_PI/LE_PI/NEE_PI/GPP_PI/RECO_PI)。
#      気象 (SW_IN/TA/TS/SWC/P) は実測のみで _F/_MDS 相当が無い。
#   3. QC フラグは `_SSITC_TEST_` (乱流検定 0/1/2)。値列名+"_QC" は無い。
# 欠測センチネル・TIMESTAMP_START 形式は JapanFlux と同じ (-9999, YYYYMMDDHHMM)。
VPD_FROM_TA_RH = "@VPD_from_TA_RH"   # var_map 上で「導出変数」を示すマーカ

DEFAULT_VAR_MAP_BASE: dict[str, str] = {
    "Rg":  "SW_IN_1_1_1",
    "Ta":  "TA_1_1_1",
    "VPD": VPD_FROM_TA_RH,          # TA + RH から導出
    "Ts":  "TS_1_1_1",
    "P":   "P_1_1_1",
    "th":  "SWC_1_1_1",
    "gH":  "H_PI_1_1_1",
    "gLE": "LE_PI_1_1_1",
    "GER": "RECO_PI_1_1_1",
    "NEE": "NEE_PI_1_1_1",
    "GEP": "GPP_PI_1_1_1",
}

# BASE 形式の相対湿度列 (VPD 導出の入力) と、乱流 QC フラグの割り当て。
BASE_RH_COL = "RH_1_1_1"
# 各 R&K 変数 → SSITC 検定フラグ列。気象・土壌は QC 無し (None=常に採用)。
BASE_QC_MAP: dict[str, str | None] = {
    "Rg": None, "Ta": None, "VPD": None, "Ts": None, "P": None, "th": None,
    "gH":  "H_SSITC_TEST_1_1_1",
    "gLE": "LE_SSITC_TEST_1_1_1",
    "GER": "FC_SSITC_TEST_1_1_1",   # 派生炭素は乱流CO2フラックスの検定を代理に
    "NEE": "FC_SSITC_TEST_1_1_1",
    "GEP": "FC_SSITC_TEST_1_1_1",
}

# --- A3 拡張: ChinaFlux 形式 (65 サイト・30 分 HH xlsx) -------------------------
# 列名は独自 (DR/Ta_above/NEP_sc/ER_DT…)、時刻は Year/Month/Day/hour/min の 5 列、
# ヘッダ行の次に「単位行」がある (自動スキップ)。gap-fill 済みの storage-corrected 値
# (_sc) を採る。NEP=−NEE の符号逆・炭素の単位 (mg CO2/m^2/s) は、情報量が単調変換・
# 符号反転・単位に不変なので解析結果に影響しない。ER_DT/GPP_DT は大文字/小文字が
# サイトで揺れる (ER_dt 等) ため、reader 側で大小無視して解決する。
CHINAFLUX_TIME_COLS = ["Year", "Month", "Day", "hour", "min"]

DEFAULT_VAR_MAP_CHINAFLUX: dict[str, str] = {
    "Rg":  "DR",         # 下向き短波放射
    "Ta":  "Ta_above",   # 林冠上の気温
    "VPD": "VPD",
    "Ts":  "Ts",         # 地温
    "P":   "Rain",       # 降水 ("P" は気圧なので使わない)
    "th":  "SW1",        # 浅層土壌水分
    "gH":  "H_sc",       # 顕熱 (貯熱補正)
    "gLE": "LE_sc",      # 潜熱 (貯熱補正)
    "GER": "ER_DT",      # 生態系呼吸 (昼分割; ER_dt も許容)
    "NEE": "NEP_sc",     # NEP = −NEE (符号逆; 情報量に不変)
    "GEP": "GPP_DT",     # 総一次生産 (昼分割; GPP_dt も許容)
}

# 形式名 → 既定 var_map。SiteSpec.fmt で選択する。
_FORMAT_MAPS: dict[str, dict[str, str]] = {
    "japanflux": DEFAULT_VAR_MAP,
    "base": DEFAULT_VAR_MAP_BASE,
    "chinaflux": DEFAULT_VAR_MAP_CHINAFLUX,
}

# データソース名 → HH ファイル glob。SiteSpec.source で選択する。
_SOURCE_GLOBS: dict[str, str] = {
    "ALLVARS": "**/*ALLVARS_HH_*.csv",
    "COREVARS": "**/*COREVARS_HH_*.csv",
}


@dataclass(frozen=True)
class SiteSpec:
    """1 サイトのデータ位置とメタ情報。"""

    code: str                        # 例 "JP-Tak"
    data_dir: str                    # サイトのルート (配布 ID サブフォルダを内包)
    # データソース: "ALLVARS" (既定, 全変数・上位互換) / "COREVARS" (代表列のみ)。
    # JapanFlux2024 は R&K の 11 変数について両者で同じ参照列名・同じ結果になるが、
    # ALLVARS は NT/DT 分割・土壌の深さ別・u* 閾値の選択肢を残すため既定に採る。
    source: str = "ALLVARS"
    # 30 分値ファイルの glob を明示したい場合に指定 (None なら source から導出)。
    # JapanFlux2024 は <root>/<配布ID>/DATA/<SOURCE>/FLX_<code>_JapanFLUX2024_
    # <SOURCE>_HH_*.csv のように配布 ID フォルダが挟まるため ** で再帰探索する。
    # BASE 形式 (韓国・中国) は AMF_<code>_BASE_HH_*.csv 等を明示する。
    hh_glob: str | None = None
    # データ規約: "japanflux" (既定, FLUXNET2015 名) / "base" (AmeriFlux/ICOS/KoFlux)。
    fmt: str = "japanflux"
    var_overrides: dict[str, str] = field(default_factory=dict)
    description: str = ""

    @property
    def glob(self) -> str:
        """HH ファイル探索 glob。明示指定を優先し、無ければ source から導出。"""
        if self.hh_glob:
            return self.hh_glob
        return _SOURCE_GLOBS.get(self.source, _SOURCE_GLOBS["ALLVARS"])

    def var_map(self) -> dict[str, str]:
        """R&K 表記 → 実カラム名。形式ごとの既定に上書きを適用したもの。"""
        m = dict(_FORMAT_MAPS.get(self.fmt, DEFAULT_VAR_MAP))
        m.update(self.var_overrides)
        return m

    def columns(self) -> list[str]:
        """正準順序での実カラム名リスト (11 本)。"""
        m = self.var_map()
        return [m[v] for v in RK_VARS]


# 実データはユーザーのローカルディスク (/mnt/hdd) 上。JP_Tak のみアンダースコア。
SITES: dict[str, SiteSpec] = {
    "JP-Tak": SiteSpec(
        code="JP-Tak",
        data_dir="/mnt/hdd/JAPANFLUX/JP_Tak",
        description="Takayama deciduous broadleaf forest (24-yr record)",
    ),
    "JP-Ta2": SiteSpec(
        code="JP-Ta2",
        data_dir="/mnt/hdd/JAPANFLUX/JP-Ta2",
        description="Takayama evergreen coniferous forest",
    ),
    "JP-BBY": SiteSpec(
        code="JP-BBY",
        data_dir="/mnt/hdd/JAPANFLUX/JP-BBY",
        description="Bibai bog (wetland)",
    ),
    "JP-Mse": SiteSpec(
        code="JP-Mse",
        data_dir="/mnt/hdd/JAPANFLUX/JP-Mse",
        description="Mase paddy (rice)",
    ),
    # --- A3 (日中韓) 拡張テンプレート ---------------------------------------
    # AmeriFlux/ICOS/KoFlux BASE 形式 (30 分値, `_1_1_1` 位置修飾子)。
    # 使い方: data_dir を実データのルートに、hh_glob を実ファイル名に
    # 合わせて変更 (または新エントリを複製) してから inspect_site を実行する。
    # VPD は TA+RH から自動導出、QC は _SSITC_TEST_ を使う。
    "A3-base": SiteSpec(
        code="A3-base",
        data_dir="/mnt/hdd/A3/SITE",     # ← 実パスに変更
        hh_glob="**/*BASE_HH_*.csv",     # ← 実ファイル名に合わせる
        fmt="base",
        description="A3 (Korea/China) BASE-format 30-min site — edit data_dir/glob",
    ),
}


# ---------------------------------------------------------------------------
# ローカルディスクからの自動発見 (/mnt/hdd/JAPANFLUX に大量サイトがある場合)
# ---------------------------------------------------------------------------
JAPANFLUX_ROOT = "/mnt/hdd/JAPANFLUX"
# JapanFlux2024 ファイル名 FLX_<code>_JapanFLUX2024_COREVARS_HH_*.csv から
# サイトコード (例 "JP-Tak") を取り出す。国コード2文字 + "-" + サイト名。
_CODE_RE = re.compile(r"FLX_([A-Za-z]{2}-[A-Za-z0-9]+)_", re.IGNORECASE)


def _code_from_filename(name: str) -> str | None:
    m = _CODE_RE.search(name)
    return m.group(1) if m else None


@lru_cache(maxsize=8)
def discover_japanflux_sites(
    root: str = JAPANFLUX_ROOT,
    glob: str = "**/*ALLVARS_HH_*.csv",
) -> dict[str, SiteSpec]:
    """``root`` 以下を走査し、見つかった ALLVARS HH から SiteSpec を組む。

    ファイル名からサイトコードを抽出し、``root`` 直下のフォルダを data_dir に採る。
    手登録の :data:`SITES` に無いサイトも ``get_site`` から使えるようにする。
    root が存在しない環境 (このリポジトリのCI等) では空 dict を返す。
    """
    root_p = Path(root)
    found: dict[str, SiteSpec] = {}
    if not root_p.exists():
        return found
    for f in sorted(root_p.glob(glob)):
        code = _code_from_filename(f.name)
        if code is None or code in found:
            continue
        try:
            site_root = root_p / f.relative_to(root_p).parts[0]
        except (ValueError, IndexError):
            site_root = f.parent
        found[code] = SiteSpec(
            code=code, data_dir=str(site_root),
            description="auto-discovered JapanFlux2024 site",
        )
    return found


CHINAFLUX_ROOT = "/mnt/hdd/ChinaFlux"


@lru_cache(maxsize=8)
def discover_chinaflux_sites(root: str = CHINAFLUX_ROOT) -> dict[str, SiteSpec]:
    """ChinaFlux ルート直下の各サイトフォルダを ChinaFlux 形式サイトとして登録する。

    <root>/<CN-XXX>/<生態系名>/<年>/XXX_YYYY_HH_ALL.xlsx。data_dir はサイトフォルダ、
    30 分 HH は ``**/*_HH_ALL.xlsx`` で再帰探索。生態系名フォルダから粗いタイプも拾える。
    """
    root_p = Path(root)
    found: dict[str, SiteSpec] = {}
    if not root_p.exists():
        return found
    for d in sorted(p for p in root_p.iterdir() if p.is_dir()):
        code = d.name                         # 例 "CN-AaG"
        eco = ""
        subs = [s for s in d.iterdir() if s.is_dir()]
        if subs:
            eco = subs[0].name                # 例 "Ansai_grassland"
        found[code] = SiteSpec(
            code=code, data_dir=str(d), fmt="chinaflux",
            hh_glob="**/*_HH_ALL.xlsx",
            description=f"ChinaFlux site ({eco})" if eco else "ChinaFlux site",
        )
    return found


def get_site(code: str) -> SiteSpec:
    """サイトを取得。手登録 :data:`SITES` を優先し、無ければローカル自動発見を試す。"""
    if code in SITES:
        return SITES[code]
    disc = discover_japanflux_sites()
    if code in disc:
        return disc[code]
    cn = discover_chinaflux_sites()
    if code in cn:
        return cn[code]
    known = sorted(set(SITES) | set(disc) | set(cn))
    raise KeyError(f"unknown site {code!r}; known: {known}")


def resolve_qc_columns(header: list[str], site: SiteSpec) -> dict[str, str | None]:
    """各 R&K 変数の QC 列を解決する。

    値列 + ``_QC`` が実ヘッダにあればそれを、無ければ None（＝QC 無しで常に実測扱い）。
    派生炭素 (RECO/GPP) は自前 QC を持たないため NEE の QC を代理に使う。
    BASE 形式では `_SSITC_TEST_` 乱流検定フラグ (:data:`BASE_QC_MAP`) を割り当てる。
    """
    hset = set(header)
    if site.fmt == "base":
        return {v: (c if c in hset else None) for v, c in BASE_QC_MAP.items()}
    vmap = site.var_map()
    nee_qc = vmap["NEE"] + "_QC"
    out: dict[str, str | None] = {}
    for v in RK_VARS:
        own = vmap[v] + "_QC"
        if own in hset:
            out[v] = own
        elif v in ("GER", "GEP") and nee_qc in hset:
            out[v] = nee_qc          # 派生炭素は NEE の品質を継承
        else:
            out[v] = None
    return out
