#!/usr/bin/env bash
# Run analysis_C_v1.py with configurable data paths
#
# Usage:
#   ./run_analysis_C.sh                       # Use env vars or defaults
#   ./run_analysis_C.sh --help                # Show options
#   ./run_analysis_C.sh --oran /path --tarazona /path
#
# Environment variables:
#   BAKANPOSS_ORAN_EC    - Path to Oran EC CSV
#   BAKANPOSS_TARA_EC    - Path to Tarazona EC CSV
#   BAKANPOSS_NDVI_FILE  - Path to NDVI CSV
#   BAKANPOSS_OUTPUT_DIR - Output directory

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Parse command line arguments
ORAN_EC=""
TARA_EC=""
NDVI_FILE=""
OUTPUT_DIR=""
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case $1 in
        --oran)
            ORAN_EC="$2"
            shift 2
            ;;
        --tarazona)
            TARA_EC="$2"
            shift 2
            ;;
        --ndvi)
            NDVI_FILE="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --help)
            cat << 'EOF'
Usage: ./run_analysis_C.sh [OPTIONS]

OPTIONS:
    --oran PATH         Path to Oran EC CSV file
    --tarazona PATH     Path to Tarazona EC CSV file
    --ndvi PATH         Path to NDVI CSV file
    --output DIR        Output directory
    --dry-run           Show configuration without running
    --help              Show this message

ENVIRONMENT VARIABLES:
    BAKANPOSS_ORAN_EC    - Path to Oran EC CSV (default: from analysis_A_v9)
    BAKANPOSS_TARA_EC    - Path to Tarazona EC CSV (default: from analysis_A_v9)
    BAKANPOSS_NDVI_FILE  - Path to NDVI CSV (default: ./MOD13Q1-*.csv or /mnt/hdd/...)
    BAKANPOSS_OUTPUT_DIR - Output directory (default: ./output_analysis_C_v1)

EXAMPLE:
    # Using environment variables
    export BAKANPOSS_ORAN_EC="/data/Oran_EC.csv"
    export BAKANPOSS_TARA_EC="/data/Tarazona_EC.csv"
    ./run_analysis_C.sh

    # Using command line arguments
    ./run_analysis_C.sh \
        --oran /data/Oran_EC.csv \
        --tarazona /data/Tarazona_EC.csv \
        --ndvi /data/NDVI.csv \
        --output ./results
EOF
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Resolve paths from environment or defaults
if [[ -z "$ORAN_EC" ]]; then
    ORAN_EC="${BAKANPOSS_ORAN_EC:-}"
fi

if [[ -z "$TARA_EC" ]]; then
    TARA_EC="${BAKANPOSS_TARA_EC:-}"
fi

if [[ -z "$NDVI_FILE" ]]; then
    NDVI_FILE="${BAKANPOSS_NDVI_FILE:-}"
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="${BAKANPOSS_OUTPUT_DIR:-./output_analysis_C_v1}"
fi

# Try to detect from analysis_A_v9 if not set
if [[ -z "$ORAN_EC" ]] || [[ -z "$TARA_EC" ]]; then
    python3 << 'PYEOF'
import sys
sys.path.insert(0, '.')
try:
    from analysis_A_v9 import PATHS
    print(f"ORAN_EC={PATHS['oran_ec']}")
    print(f"TARA_EC={PATHS['tara_ec']}")
except Exception as e:
    print(f"# Could not load from analysis_A_v9: {e}", file=sys.stderr)
PYEOF
fi

# Display configuration
echo "=============================================================================="
echo "Analysis C Configuration"
echo "=============================================================================="
echo "Oran EC:     $ORAN_EC"
echo "Tarazona EC: $TARA_EC"
echo "NDVI file:   $NDVI_FILE"
echo "Output dir:  $OUTPUT_DIR"
echo "=============================================================================="

# Validate paths
ERRORS=0

if [[ -z "$ORAN_EC" ]]; then
    echo "❌ Oran EC path not configured"
    ERRORS=$((ERRORS + 1))
elif [[ ! -f "$ORAN_EC" ]]; then
    echo "❌ Oran EC file not found: $ORAN_EC"
    ERRORS=$((ERRORS + 1))
fi

if [[ -z "$TARA_EC" ]]; then
    echo "❌ Tarazona EC path not configured"
    ERRORS=$((ERRORS + 1))
elif [[ ! -f "$TARA_EC" ]]; then
    echo "❌ Tarazona EC file not found: $TARA_EC"
    ERRORS=$((ERRORS + 1))
fi

if [[ -z "$NDVI_FILE" ]]; then
    # Try to find NDVI file automatically
    if [[ -f "MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv" ]]; then
        NDVI_FILE="./MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv"
    elif [[ -f "/mnt/hdd/Dataset/MOD13Q1_NDVI_EVI/MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv" ]]; then
        NDVI_FILE="/mnt/hdd/Dataset/MOD13Q1_NDVI_EVI/MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv"
    else
        echo "❌ NDVI file not configured or found"
        ERRORS=$((ERRORS + 1))
    fi
elif [[ ! -f "$NDVI_FILE" ]]; then
    echo "❌ NDVI file not found: $NDVI_FILE"
    ERRORS=$((ERRORS + 1))
fi

if [[ $ERRORS -gt 0 ]]; then
    echo ""
    echo "Please configure paths via:"
    echo "  - Command line arguments (--oran, --tarazona, --ndvi, --output)"
    echo "  - Environment variables (BAKANPOSS_*)"
    echo "  - analysis_A_v9.py defaults"
    exit 1
fi

echo "✓ All paths verified"
echo ""

if [[ $DRY_RUN -eq 1 ]]; then
    echo "(Dry-run mode: configuration validated but not executing)"
    exit 0
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run analysis_C_v1.py
echo "=============================================================================="
echo "Starting Analysis C v1"
echo "=============================================================================="
echo ""

python3 << PYEOF
import sys
from pathlib import Path

# Setup paths
SCRIPT_DIR = Path('.')
sys.path.insert(0, str(SCRIPT_DIR))

# Override paths in analysis_A_v9
import analysis_A_v9
analysis_A_v9.PATHS["oran_ec"] = Path("$ORAN_EC")
analysis_A_v9.PATHS["tara_ec"] = Path("$TARA_EC")

# Override analysis_C_v1 globals
import analysis_C_v1
analysis_C_v1.NDVI_APPEEARS_CSV = Path("$NDVI_FILE")
analysis_C_v1.SAVE_DIR = Path("$OUTPUT_DIR")

# Execute the main analysis block from analysis_C_v1
exec(open("analysis_C_v1.py").read(), analysis_C_v1.__dict__)
PYEOF

echo ""
echo "✓ Analysis complete: $OUTPUT_DIR"
