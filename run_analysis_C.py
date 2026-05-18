#!/usr/bin/env python3
"""
Direct runner for analysis_C_v1.py with environment variable configuration.

Usage:
    python3 run_analysis_C.py                          # Use defaults or env vars
    python3 run_analysis_C.py --help                   # Show options
    python3 run_analysis_C.py --oran /path/to/oran.csv --tarazona /path/to/tara.csv

Environment variables:
    BAKANPOSS_ORAN_EC    - Path to Oran EC CSV (default: from analysis_A_v9)
    BAKANPOSS_TARA_EC    - Path to Tarazona EC CSV (default: from analysis_A_v9)
    BAKANPOSS_NDVI_FILE  - Path to NDVI CSV (default: auto-detect)
    BAKANPOSS_OUTPUT_DIR - Output directory (default: ./output_analysis_C_v1)
"""

import argparse
import sys
from pathlib import Path
import os

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))


def get_path_config(args):
    """Resolve paths from args, environment, or defaults."""

    # Try to load defaults from analysis_A_v9
    oran_default = None
    tara_default = None
    try:
        from analysis_A_v9 import PATHS as A_PATHS
        oran_default = A_PATHS.get("oran_ec")
        tara_default = A_PATHS.get("tara_ec")
    except (ImportError, KeyError):
        pass

    paths = {
        "oran_ec": Path(args.oran) if args.oran else Path(
            os.environ.get("BAKANPOSS_ORAN_EC") or oran_default or ""
        ),
        "tara_ec": Path(args.tarazona) if args.tarazona else Path(
            os.environ.get("BAKANPOSS_TARA_EC") or tara_default or ""
        ),
        "ndvi_file": None,
        "output_dir": Path(args.output or os.environ.get(
            "BAKANPOSS_OUTPUT_DIR") or "./output_analysis_C_v1"
        ),
    }

    # Find NDVI file
    if args.ndvi:
        paths["ndvi_file"] = Path(args.ndvi)
    elif os.environ.get("BAKANPOSS_NDVI_FILE"):
        paths["ndvi_file"] = Path(os.environ["BAKANPOSS_NDVI_FILE"])
    else:
        # Try to auto-detect
        for candidate in [
            SCRIPT_DIR / "MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv",
            Path("/mnt/hdd/Dataset/MOD13Q1_NDVI_EVI/MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv"),
        ]:
            if candidate.exists():
                paths["ndvi_file"] = candidate
                break

    return paths


def validate_paths(paths):
    """Check that all required files exist."""
    errors = []

    if not paths["oran_ec"] or not paths["oran_ec"].exists():
        errors.append(f"Oran EC: {paths['oran_ec']} (not found or not configured)")

    if not paths["tara_ec"] or not paths["tara_ec"].exists():
        errors.append(f"Tarazona EC: {paths['tara_ec']} (not found or not configured)")

    if not paths["ndvi_file"] or not paths["ndvi_file"].exists():
        errors.append(f"NDVI file: {paths['ndvi_file']} (not found or not configured)")

    return errors


def main():
    parser = argparse.ArgumentParser(
        description="Run analysis_C_v1 with configurable data paths",
        epilog=__doc__,
    )
    parser.add_argument("--oran", help="Path to Oran EC CSV")
    parser.add_argument("--tarazona", help="Path to Tarazona EC CSV")
    parser.add_argument("--ndvi", help="Path to NDVI CSV")
    parser.add_argument("--output", help="Output directory")
    parser.add_argument("--dry-run", action="store_true", help="Validate config only")

    args = parser.parse_args()
    paths = get_path_config(args)

    # Display configuration
    print("=" * 70)
    print("Analysis C Configuration")
    print("=" * 70)
    print(f"Oran EC:     {paths['oran_ec']}")
    print(f"Tarazona EC: {paths['tara_ec']}")
    print(f"NDVI file:   {paths['ndvi_file']}")
    print(f"Output dir:  {paths['output_dir']}")
    print("=" * 70)

    # Validate
    errors = validate_paths(paths)
    if errors:
        print("\n❌ Configuration errors:")
        for error in errors:
            print(f"   {error}")
        sys.exit(1)

    print("\n✓ All paths verified")

    if args.dry_run:
        print("\n(Dry-run: config validated)")
        return

    # Create output directory
    paths["output_dir"].mkdir(parents=True, exist_ok=True)

    # Configure and run
    print(f"\n{'='*70}")
    print("Starting Analysis C v1")
    print(f"{'='*70}\n")

    import analysis_A_v9
    analysis_A_v9.PATHS["oran_ec"] = paths["oran_ec"]
    analysis_A_v9.PATHS["tara_ec"] = paths["tara_ec"]

    import analysis_C_v1
    analysis_C_v1.NDVI_APPEEARS_CSV = paths["ndvi_file"]
    analysis_C_v1.SAVE_DIR = paths["output_dir"]

    # Run the analysis by executing the main block
    exec(open(str(SCRIPT_DIR / "analysis_C_v1.py")).read())


if __name__ == "__main__":
    main()
