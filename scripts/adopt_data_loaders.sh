#!/usr/bin/env bash
# adopt_data_loaders.sh — 解析 A/B ブランチに data_loaders.py を採用するヘルパー
#
# 想定使用法 (A ブランチの作業者):
#   git checkout claude/quantify-water-divergence-LxPGr
#   bash scripts/adopt_data_loaders.sh
#
# このスクリプトは:
#   1. C ブランチから data_loaders.py を取り込む
#   2. analysis_A_v*.py / analysis_B*.py 内のローダ呼び出し箇所を grep で
#      洗い出し、置換候補を表示 (実際の編集は手動)
#   3. data_loaders がインポート可能か健全性テストを走らせる
#
# 安全のため、コードの自動書き換えはしない (見落としや false rename を
# 避けるため)。差分は出すので手動で適用する。

set -euo pipefail

C_BRANCH="claude/ndvi-flux-analysis-iaS0b"
TARGET_FILE="data_loaders.py"

cd "$(git rev-parse --show-toplevel)"

current_branch=$(git rev-parse --abbrev-ref HEAD)
echo "現在のブランチ: $current_branch"
if [ "$current_branch" = "$C_BRANCH" ]; then
    echo "ERROR: あなたは既に $C_BRANCH にいます。A または B のブランチで実行してください。"
    exit 1
fi

# 0. C ブランチ存在チェック
if ! git rev-parse --verify "origin/$C_BRANCH" >/dev/null 2>&1; then
    echo "ERROR: origin/$C_BRANCH が見つかりません。git fetch を実行してください。"
    exit 1
fi

# 1. data_loaders.py を取り込み
echo
echo "STEP 1: $C_BRANCH から $TARGET_FILE を取り込み"
if [ -f "$TARGET_FILE" ]; then
    echo "  既存の $TARGET_FILE があります。バックアップして上書きします。"
    cp "$TARGET_FILE" "${TARGET_FILE}.bak.$(date +%Y%m%d_%H%M%S)"
fi
git checkout "origin/$C_BRANCH" -- "$TARGET_FILE"
echo "  $TARGET_FILE を取り込みました"

# 2. 既存ローダの使用箇所を検出
echo
echo "STEP 2: 既存ローダ呼び出しの検出 (手動置換が必要な箇所)"
echo "----"
echo "  load_oran_ec / load_tarazona_ec / normalize_swc を呼んでいる行:"
grep -rn -E "load_oran_ec|load_tarazona_ec|normalize_swc" \
    --include="*.py" \
    --exclude="$TARGET_FILE" \
    . 2>/dev/null \
    | grep -v -E "(^|/)data_loaders\.py:" \
    | grep -v "scripts/adopt_data_loaders.sh" \
    || echo "  (該当なし — 既に切り替え済みか?)"
echo "----"

cat <<'EOF'

  置換ガイド (各ファイルで手動編集):

    # Before:
    from analysis_A_v9 import load_oran_ec, load_tarazona_ec, normalize_swc
    oran_ec = load_oran_ec(PATHS["oran_ec"])
    tara_ec = load_tarazona_ec(PATHS["tara_ec"])

    # After:
    from data_loaders import (load_oran_ec_clean, load_tarazona_ec_clean,
                                normalize_swc)
    oran_ec = load_oran_ec_clean(PATHS["oran_ec"])
    tara_ec = normalize_swc(load_tarazona_ec_clean(PATHS["tara_ec"]), "Tarazona")

EOF

# 3. data_loaders.py の健全性テスト (ファイル パスは環境依存なので
#    インポートチェックのみ実施)
echo
echo "STEP 3: data_loaders.py のインポートチェック"
python3 - <<'PYEOF'
import importlib.util, sys
spec = importlib.util.spec_from_file_location("data_loaders", "data_loaders.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
required = ["load_oran_ec_clean", "load_tarazona_ec_clean", "normalize_swc",
            "EF_DENOM_MIN", "SENTINEL_THR"]
missing = [n for n in required if not hasattr(m, n)]
if missing:
    print(f"  ERROR: 不足: {missing}")
    sys.exit(1)
print("  OK: 必須シンボル全部あり:", required)
PYEOF

# 4. 期待される動作の説明
echo
cat <<'EOF'
STEP 4: 動作確認 (任意)
  実データで動かすにはサイトの CSV パスが必要なため、本スクリプトは
  実行しない。手動で以下を試してください:

    python3 -c "
    from data_loaders import load_oran_ec_clean, load_tarazona_ec_clean
    oran = load_oran_ec_clean('/your/path/to/Oran_Ameriflux_Cereal_ASV_CLEAN_2018_2020.csv')
    tara = load_tarazona_ec_clean('/your/path/to/Daily_Summary_Filtered_forPred_ActEne26.csv')
    assert 80 < oran['Rn'].median() < 110, f'Oran Rn med 異常: {oran[\"Rn\"].median()}'
    assert tara['ET'].median() > 2, 'Tarazona ET が小さすぎる'
    print('OK')
    "

  詳細は reports/migration_to_data_loaders.md を参照。

完了: data_loaders.py の取り込みが終わりました。
       上の手動置換ガイドに従って各 analysis スクリプトを更新してください。
EOF
