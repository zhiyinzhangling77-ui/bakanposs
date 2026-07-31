"""ローカルの JapanFlux ルートを走査し、利用可能なサイトを一覧する。

ローカル (/mnt/hdd が見える環境) で:

    python -m japanflux_pn.discover
    python -m japanflux_pn.discover --root /mnt/hdd/JAPANFLUX

各サイトのコード・観測期間・COREVARS HH ファイル数を表示する。手登録済み
(sites.SITES) か自動発見かも示す。ここに出たコードはそのまま
``inspect_site --site <code>`` / ``run_site --site <code>`` に渡せる。
"""

from __future__ import annotations

import argparse

import pandas as pd

from .sites import SITES, JAPANFLUX_ROOT, discover_japanflux_sites
from .preprocess import find_corevars_files


def _span(path) -> str:
    try:
        df = pd.read_csv(path, usecols=["TIMESTAMP_START"])
        ts = pd.to_datetime(df["TIMESTAMP_START"].astype("int64").astype(str),
                            format="%Y%m%d%H%M")
        return f"{ts.min():%Y-%m} .. {ts.max():%Y-%m}"
    except Exception:  # noqa: BLE001
        return "(span unreadable)"


def report(root: str = JAPANFLUX_ROOT) -> None:
    disc = discover_japanflux_sites(root)
    # 指定 root で発見したものを土台に、手登録サイト (curated data_dir) で上書き。
    resolved = dict(disc)
    resolved.update(SITES)
    codes = sorted(resolved)
    print(f"### JapanFlux root: {root}")
    print(f"### {len(codes)} site(s)  (手登録 {len(SITES)} / 自動発見 {len(disc)})\n")
    if not codes:
        print("  (COREVARS HH が見つかりません。--root を確認してください)")
        return
    for code in codes:
        tag = "registered" if code in SITES else "auto"
        try:
            site = resolved[code]
            files = find_corevars_files(site)
            lo = _span(files[0])
            hi = _span(files[-1]) if len(files) > 1 else lo
            span = lo if len(files) == 1 else f"{lo.split('..')[0].strip()} .. {hi.split('..')[-1].strip()}"
            print(f"  {code:<10} [{tag:^10}] files={len(files):<3} span={span}")
            print(f"             dir: {site.data_dir}")
        except Exception as e:  # noqa: BLE001
            print(f"  {code:<10} [{tag:^10}] (読み取り不可: {e})")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="discover local JapanFlux sites")
    p.add_argument("--root", default=JAPANFLUX_ROOT)
    args = p.parse_args(argv)
    report(args.root)


if __name__ == "__main__":
    main()
