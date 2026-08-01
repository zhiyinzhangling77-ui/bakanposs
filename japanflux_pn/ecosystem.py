"""各サイトの BADM から IGBP 土地被覆コードを抽出し、生態系タイプで分類する。

サイトコードからの推測でなく、配布メタデータ (BADM) の IGBP を一次情報として使い、
「森林 / 草原 / 水田 (CRO) / 湿原 …」でグループ分けする。batch_oinfo などの結果を
生態系タイプで層別し、反復数を客観的に数えるための土台。

    python -m japanflux_pn.ecosystem                 # 全サイトの IGBP 一覧
    python -m japanflux_pn.ecosystem --sites JP-Mse JP-Tak CN-HaM

BADM の書式はデータセットで揺れるため、寛容に走査する: サイトフォルダ内の BADM 風
csv を全部読み、(1)「IGBP」というラベル行の値、(2) セル中に現れる既知 IGBP コード、
の順で拾う。見つからなければ "?" を返す（その場合 --dump で中身を確認できる）。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .sites import JAPANFLUX_ROOT

# IGBP コード → (日本語ラベル, 粗い生態系グループ)
IGBP_INFO: dict[str, tuple[str, str]] = {
    "ENF": ("常緑針葉樹林", "森林"),
    "EBF": ("常緑広葉樹林", "森林"),
    "DNF": ("落葉針葉樹林", "森林"),
    "DBF": ("落葉広葉樹林", "森林"),
    "MF":  ("混交林", "森林"),
    "CSH": ("閉低木林", "低木"),
    "OSH": ("開低木林", "低木"),
    "WSA": ("疎林サバンナ", "サバンナ"),
    "SAV": ("サバンナ", "サバンナ"),
    "GRA": ("草原", "草原"),
    "WET": ("湿原", "湿地"),
    "CRO": ("農地(水田/畑)", "農地"),
    "CVM": ("農地/自然モザイク", "農地"),
    "URB": ("市街", "その他"),
    "SNO": ("雪氷", "その他"),
    "BSV": ("裸地/疎植生", "その他"),
    "WAT": ("水域", "その他"),
}
IGBP_CODES = set(IGBP_INFO)

# BADM がローカルに無い場合のフォールバック分類。
# ★ = 以前アップされた BADM で確定 (SITE_METADATA.md)。それ以外は AsiaFlux/JapanFlux の
# 公表サイト情報に基づく暫定（IGBP 細分類は要確認だが「森林/草原/水田/湿地」の粗グループは
# 確度が高い）。BADM をダウンロードできれば igbp_for_site が自動でこれを上書きする。
KNOWN_IGBP: dict[str, str] = {
    # --- BADM 確定 (★) ---
    "JP-Tak": "DBF", "JP-Ta2": "ENF", "JP-Mse": "CRO", "JP-BBY": "WET",
    # --- 草原 (確度高: モンゴル・内モンゴル・青海) ---
    "CN-HaM": "GRA", "MN-Kbu": "GRA", "MN-Hst": "GRA", "MN-Nkh": "GRA",
    "MN-Skt": "GRA", "MN-Udg": "GRA",
    "CN-In1": "GRA", "CN-In2": "GRA", "CN-In3": "GRA", "CN-In4": "GRA",
    "CN-In5": "GRA", "CN-In6": "GRA", "CN-In7": "GRA", "CN-In8": "GRA",
    # --- シベリア落葉針葉樹 (カラマツ, 確度高) ---
    "RU-Ege": "DNF", "RU-SkP": "DNF", "RU-NeF": "DNF",
}
# 確定 (BADM) と暫定 (文献) の区別
BADM_CONFIRMED = {"JP-Tak", "JP-Ta2", "JP-Mse", "JP-BBY"}

# BADM 風ファイルを探す glob (JapanFlux2024: 名称は多様)
_BADM_GLOBS = ("**/*BADM*.csv", "**/*Site_General*.csv", "**/*BIF*.csv")


def _iter_badm_files(data_dir: str, code: str | None = None) -> list[Path]:
    """サイトの BADM 風 csv を探す。data_dir 下だけでなく、その親 (JAPANFLUX root)
    直下の「サイトコードを含む BADM 兄弟フォルダ」も探す。

    自動発見サイトの data_dir は ALLVARS 専用フォルダ
    (FLX_<code>_JapanFLUX2024_ALLVARS_...) で、BADM は root 直下の別配布フォルダ
    (FLX_<code>_JapanFLUX2024_BADM_...) にあるため、両方を走査する。
    """
    seen: list[Path] = []
    roots = [Path(data_dir)]
    if code:                                   # root 直下の <code> を含む兄弟フォルダ
        parent = Path(data_dir).parent
        roots += [p for p in sorted(parent.glob(f"*{code}*"))
                  if p.is_dir() and p != Path(data_dir)]
    for root in roots:
        for g in _BADM_GLOBS:
            for f in sorted(root.glob(g)):
                if f not in seen:
                    seen.append(f)
    return seen


def _igbp_from_frame(df: pd.DataFrame) -> str | None:
    """DataFrame から IGBP を拾う。ラベル行優先、無ければ既知コードの直接出現。"""
    # (1) 「IGBP」ラベルを含むセルの隣/同行の値
    vals = df.astype(str)
    for r in range(len(vals)):
        row = vals.iloc[r].tolist()
        for c, cell in enumerate(row):
            if cell.strip().upper() in ("IGBP", "IGBP_CLASS", "LAND_COVER"):
                # 同行の後続セルから IGBP コードを探す
                for other in row[c + 1:] + row[:c]:
                    tok = other.strip().upper()
                    if tok in IGBP_CODES:
                        return tok
    # (2) 表中どこかに現れる既知 IGBP コード (単独セル)
    for r in range(len(vals)):
        for cell in vals.iloc[r].tolist():
            tok = cell.strip().upper()
            if tok in IGBP_CODES:
                return tok
    return None


def igbp_for_site(data_dir: str, code: str | None = None) -> str:
    """サイトの BADM 群から IGBP コードを抽出 (見つからねば '?')。"""
    for f in _iter_badm_files(data_dir, code):
        try:
            df = pd.read_csv(f, header=None, dtype=str, on_bad_lines="skip")
        except Exception:  # noqa: BLE001
            continue
        igbp = _igbp_from_frame(df)
        if igbp:
            return igbp
    return "?"


def dump_badm(data_dir: str, code: str | None = None, max_files: int = 6,
              max_rows: int = 25) -> None:
    """診断用: 見つかった BADM ファイルと中身の先頭を表示する。"""
    files = _iter_badm_files(data_dir, code)
    print(f"### BADM 探索 (data_dir={data_dir}, code={code})")
    print(f"  親: {Path(data_dir).parent}")
    if not files:
        print("  BADM 風ファイルが見つかりません。root 直下の配布フォルダ名を確認:")
        for p in sorted(Path(data_dir).parent.glob(f"*{code or ''}*"))[:20]:
            print(f"    {p.name}")
        return
    for f in files[:max_files]:
        print(f"\n--- {f} ---")
        try:
            df = pd.read_csv(f, header=None, dtype=str, on_bad_lines="skip")
            print(df.head(max_rows).to_string(max_cols=8))
        except Exception as e:  # noqa: BLE001
            print(f"  (読めません: {e})")


def classify(root: str = JAPANFLUX_ROOT,
             sites: list[str] | None = None) -> pd.DataFrame:
    from .rank_sites import _resolve_sites
    resolved = _resolve_sites(root)
    codes = sites or sorted(resolved)
    rows = []
    for code in codes:
        spec = resolved.get(code)
        igbp = igbp_for_site(spec.data_dir, code) if spec is not None else "?"
        src = "BADM"
        if igbp == "?" and code in KNOWN_IGBP:      # ローカル BADM 無し → フォールバック
            igbp = KNOWN_IGBP[code]
            src = "確定" if code in BADM_CONFIRMED else "暫定(文献)"
        elif igbp == "?" and spec is None:
            src = "(未登録)"
        lab, grp = IGBP_INFO.get(igbp, ("(不明)", "?"))
        rows.append({"site": code, "igbp": igbp, "type": grp,
                     "label": lab, "source": src})
    return pd.DataFrame(rows)


def report(root: str = JAPANFLUX_ROOT, sites: list[str] | None = None,
           csv: str | None = None) -> pd.DataFrame:
    df = classify(root, sites)
    print(f"### サイト → IGBP → 生態系グループ ({len(df)} サイト)\n")
    print(f"  {'site':<10} {'IGBP':>5}  {'グループ':<8} {'出所':<10} 詳細")
    for _, r in df.sort_values(["type", "site"]).iterrows():
        print(f"  {r['site']:<10} {r['igbp']:>5}  {r['type']:<8} "
              f"{r.get('source',''):<10} {r['label']}")
    print("\n=== 生態系グループ別サイト数 ===")
    for grp, g in df.groupby("type"):
        print(f"  {grp:<8} {len(g):>3}  ({', '.join(g['site'])})")
    if csv:
        Path(csv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(csv, index=False)
        print(f"\n[output] {csv}")
    return df


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="classify sites by BADM IGBP")
    p.add_argument("--root", default=JAPANFLUX_ROOT)
    p.add_argument("--sites", nargs="+", default=None)
    p.add_argument("--csv", default=None)
    p.add_argument("--dump", default=None,
                   help="診断: 指定サイトの BADM ファイルと中身を表示")
    args = p.parse_args(argv)
    if args.dump:
        from .rank_sites import _resolve_sites
        spec = _resolve_sites(args.root).get(args.dump)
        if spec is None:
            print(f"unknown site {args.dump!r}")
        else:
            dump_badm(spec.data_dir, args.dump)
        return
    report(args.root, args.sites, args.csv)


if __name__ == "__main__":
    main()
