#!/usr/bin/env bash
# 夜通しバッチ：TROPOMI ungridded L2 SIF を月ごとに DL→untar→タワー集約→削除→結合。
# near-daily で「呼吸の4日記憶が速い基質か」を JP-Tak で最終検証する（段2）。
#
# 使い方：
#   1. まず 2019-07 を1ヶ月手で落として `--list` で変数名/構造を確認し、下の 4 変数を埋める。
#   2. bash research/sif_ungridded_overnight.sh   （nohup 推奨: nohup bash ... > sif.log 2>&1 &）
#
# ディスク節約：各月 DL→展開→抽出→即削除。ピークは1ヶ月分（~10GB）に収まる。

set -u
cd "$(dirname "$0")/.."          # リポジトリルートへ

# ===== --list で確認して埋める（推測NG）=====================================
VAR="sif"                        # SIF 値の変数名（例 sif / SIF_743）
LATV="lat"                       # 緯度変数
LONV="lon"                       # 経度変数
TIMEV="TIME"                     # 時刻変数（無ければ空 "" にして日付はファイル名依存）
NCGLOB="*.nc"                    # 展開後の NetCDF の glob（find で確認: *.nc / *.nc4 等）
RADIUS=20                        # タワー半径[km]
SITES="JP-Tak JP-BBY"            # 2018+ と重なるサイトのみ
# 生育期 6-9月 × 2018-2021（7-8月だけにするなら 06,09 を削る）
MONTHS="2018-06 2018-07 2018-08 2018-09 2019-06 2019-07 2019-08 2019-09 \
        2020-06 2020-07 2020-08 2020-09 2021-06 2021-07 2021-08 2021-09"
# ===========================================================================

BASE="ftp://fluo.gps.caltech.edu/data/tropomi/ungridded/SIF740nm"
mkdir -p L2work sif_month

for M in $MONTHS; do
  echo "===== $M  $(date +%H:%M) ====="
  TAR="L2work/TROPO_SIF_${M}.tar.gz"
  if ! wget -q -O "$TAR" "$BASE/TROPO_SIF_${M}.tar.gz"; then
    echo "  DL失敗（存在しない月かも）: $M"; rm -f "$TAR"; continue
  fi
  tar xzf "$TAR" -C L2work 2>/dev/null || { echo "  untar失敗: $M"; rm -rf L2work/*; continue; }
  NCDIR=$(dirname "$(find L2work -name "$NCGLOB" 2>/dev/null | head -1)")
  [ -z "$NCDIR" ] && { echo "  .nc 見つからず: $M"; rm -rf L2work/*; continue; }
  TARG=""; [ -n "$TIMEV" ] && TARG="--time-var $TIMEV"
  python research/sif_extract_ungridded.py --coords site_coords.csv --ncdir "$NCDIR" \
      --glob "$NCGLOB" --var "$VAR" --lat-var "$LATV" --lon-var "$LONV" $TARG --radius "$RADIUS"
  for s in $SITES; do [ -f "${s}_sif.csv" ] && mv "${s}_sif.csv" "sif_month/${s}_${M}.csv"; done
  rm -rf L2work/*                # 容量解放
done

echo "===== 月次CSVを結合 ====="
python - "$SITES" <<'PY'
import sys, glob, pandas as pd
for s in sys.argv[1].split():
    fs = sorted(glob.glob(f"sif_month/{s}_*.csv"))
    if not fs:
        print(f"  {s}: 月次CSVなし"); continue
    df = pd.concat([pd.read_csv(f) for f in fs])
    df = df.groupby("date")["sif"].mean().reset_index().sort_values("date")
    df.to_csv(f"{s}_sif.csv", index=False)
    print(f"  {s}: {len(df)} 日 → {s}_sif.csv")
PY

echo "===== 検証（near-daily で4日記憶が分解できる）====="
for s in $SITES; do
  [ -f "${s}_sif.csv" ] && python research/sif_respiration_step38.py --site "$s" --sif "${s}_sif.csv" --qc-max 1 --deseason
done
echo "===== 完了 $(date) ====="
