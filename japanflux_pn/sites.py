"""サイトレジストリ: サイトコード → データ位置 + 変数名マッピング。

多サイト拡張 (JP-Ta2 / JP-BBY / JP-Mse) の効き所。数値カーネルは変数名を
知らないので、サイトごとの実カラム名・ディレクトリ名の違いはここで吸収する。
FLUXNET2015 の基本名は共通だが ``TS_F_MDS_1`` / ``SWC_F_MDS_1`` の層インデックス
や、湿原・水田で欠ける変数はサイト依存になりうる。
"""

from __future__ import annotations

from dataclasses import dataclass, field

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


@dataclass(frozen=True)
class SiteSpec:
    """1 サイトのデータ位置とメタ情報。"""

    code: str                        # 例 "JP-Tak"
    data_dir: str                    # サイトのルート (配布 ID サブフォルダを内包)
    # COREVARS 30 分値ファイルの glob。JapanFlux2024 は
    #   <root>/<配布ID>/DATA/COREVARS/FLX_<code>_JapanFLUX2024_COREVARS_HH_*.csv
    # のように配布 ID フォルダが挟まるため、ルートから再帰的に探す (** を使用)。
    corevars_hh_glob: str = "**/*COREVARS_HH_*.csv"
    var_overrides: dict[str, str] = field(default_factory=dict)
    description: str = ""

    def var_map(self) -> dict[str, str]:
        """R&K 表記 → 実カラム名。既定に上書きを適用したもの。"""
        m = dict(DEFAULT_VAR_MAP)
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
}


def get_site(code: str) -> SiteSpec:
    if code not in SITES:
        raise KeyError(f"unknown site {code!r}; known: {sorted(SITES)}")
    return SITES[code]


def resolve_qc_columns(header: list[str], site: SiteSpec) -> dict[str, str | None]:
    """各 R&K 変数の QC 列を解決する。

    値列 + ``_QC`` が実ヘッダにあればそれを、無ければ None（＝QC 無しで常に実測扱い）。
    派生炭素 (RECO/GPP) は自前 QC を持たないため NEE の QC を代理に使う。
    """
    hset = set(header)
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
