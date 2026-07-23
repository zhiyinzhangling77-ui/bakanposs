"""JapanFlux2024 transfer-entropy process networks (Ruddell & Kumar 2009).

情報理論に基づく生態水文プロセスネットワーク解析を JapanFlux2024 の
FLUXNET2015 互換 CSV に適用するためのパッケージ。数値カーネル
(:mod:`information_theory`) は変数名やサイトを一切知らず、前処理
(:mod:`preprocess`) と :mod:`sites` レジストリがサイト依存性を吸収する。
"""

from .config import AnalysisConfig, RK_VARS

__all__ = ["AnalysisConfig", "RK_VARS"]
