# Running Analysis C v1

This document explains how to run the comprehensive phenology × flux analysis (analysis_C_v1.py) with your data.

## Quick Start

### Using Python Runner (Recommended)

```bash
# 1. Configure your data paths as environment variables
export BAKANPOSS_ORAN_EC="/path/to/Oran_EC.csv"
export BAKANPOSS_TARA_EC="/path/to/Tarazona_EC.csv"
export BAKANPOSS_NDVI_FILE="/path/to/NDVI.csv"
export BAKANPOSS_OUTPUT_DIR="./output_analysis_C_v1"

# 2. Run the analysis
python3 run_analysis_C.py

# Or use command-line arguments directly
python3 run_analysis_C.py \
    --oran /path/to/Oran_EC.csv \
    --tarazona /path/to/Tarazona_EC.csv \
    --ndvi /path/to/NDVI.csv \
    --output ./output_analysis_C_v1
```

### Using Bash Wrapper

```bash
./run_analysis_C.sh \
    --oran /path/to/Oran_EC.csv \
    --tarazona /path/to/Tarazona_EC.csv \
    --ndvi /path/to/NDVI.csv \
    --output ./output_analysis_C_v1
```

### Using .env Configuration File

```bash
# 1. Copy the example configuration
cp .env.example .env

# 2. Edit .env with your paths
# vim .env

# 3. Source the configuration and run
source .env
python3 run_analysis_C.py
```

## Input Data Requirements

### Oran EC Data (Eddy Covariance)

**File Format**: CSV with half-hourly measurements

**Required Columns**:
- `TIMESTAMP`: Date/time in format `YYYY/MM/DD HH:MM:SS` or `YYYY/MM/DD`
  - Fallback formats supported: `YYYY-MM-DD`, ISO 8601, Julian dates
- `NETRAD`: Net radiation [W/m²]
- `LE`: Latent energy flux [W/m²]
- `H`: Sensible heat flux [W/m²]
- `G`: Ground heat flux [W/m²]
- `VPD`: Vapor pressure deficit [hPa]
- `ET`: Evapotranspiration [mm or mm/day]
- `SWC_1_1_1`: Soil water content [%]

**Expected**: ~52,606 half-hourly records (2018-2020), parsed into ~922 daily averages

### Tarazona EC Data (Daily Aggregated)

**File Format**: CSV with daily aggregated measurements

**Required Columns**:
- `date`: Date in format `YYYY-MM-DD`
- `LE_avg`: Latent energy (daily average) [W/m²]
- `H_avg`: Sensible heat (daily average) [W/m²]
- `G_avg`: Ground heat (daily average) [W/m²]
- `NetRad_avg`: Net radiation (daily average) [W/m²]
- `SWC_avg`: Soil water content (daily average) [%]
- `VPD_mean`: Vapor pressure deficit [Pa or kPa]
- `ET_sum` or `ET_avg`: Evapotranspiration [mm or other unit]

**Optional Columns** (preserved in output):
- `Irrig_mm`: Irrigation amount [mm/day]
- `Rain_mm`: Rainfall [mm/day]
- `GPP_avg`: Gross primary productivity [µmol/m²/s]
- `NDVI_orig`, `NDVI_interp`: NDVI values for sanity checks

**Expected**: ~700-750 daily records (2018-2020 growing seasons)

### NDVI Data (Satellite Vegetation)

**File Format**: CSV from Google Earth Engine (MOD13Q1 or similar)

**Source**: MODIS MOD13Q1 16-day composite (250 m resolution)

**Required Columns** (site-specific):
- At minimum: columns for date and NDVI/EVI values per site
  - e.g., `date`, `Oran_NDVI`, `Tarazona_NDVI`
  - or extracted via AppEEARS with site IDs in column names

**Expected**: ~70-90 16-day observations per site

## Configuration Options

### Command-Line Arguments

```
--oran PATH         Path to Oran EC CSV file
--tarazona PATH     Path to Tarazona EC CSV file
--ndvi PATH         Path to NDVI CSV file
--output DIR        Output directory for results
--dry-run           Validate configuration without running
--help              Show help message
```

### Environment Variables

```
BAKANPOSS_ORAN_EC    - Oran EC CSV path
BAKANPOSS_TARA_EC    - Tarazona EC CSV path
BAKANPOSS_NDVI_FILE  - NDVI CSV path
BAKANPOSS_OUTPUT_DIR - Output directory
```

### Auto-Detection

If paths are not specified, the runner attempts to:
1. Load defaults from `analysis_A_v9.py` PATHS dict
2. Look for NDVI file in:
   - Current directory: `MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv`
   - `/mnt/hdd/Dataset/MOD13Q1_NDVI_EVI/MOD13Q1-NDVI-EVI-MOD13Q1-061-results.csv`

## Output Files

Analysis results are saved to the output directory (default: `./output_analysis_C_v1/`) and include:

### Data Exports

- `C_Oran_merged.csv` - Oran EC + NDVI merged daily data
- `C_Tarazona_merged.csv` - Tarazona EC + NDVI merged daily data

### Diagnostic Plots

- `NDVI_timeseries_Oran.png` - Oran NDVI seasonal cycle with Savitzky-Golay smooth
- `NDVI_timeseries_Tarazona.png` - Tarazona NDVI seasonal cycle
- `NDVI_vs_flux_Oran.png` - NDVI scatter plots vs LE, EF, ET, ET
- `NDVI_vs_flux_Tarazona.png` - Same for Tarazona
- `growing_phenology_Oran.png` - NDVI threshold definition and phase classification
- `growing_phenology_Tarazona.png` - Same for Tarazona

### Hypothesis Test Results

- `H1_H8_irrigation_lag_analysis.png` - EF decay with days-since-irrigation (Tarazona)
- `H2_ndvi_saturation_Oran.png` - EVI-based saturation check
- `H2_ndvi_saturation_Tarazona.png` - Same for Tarazona
- `H4_albedo_feedback_Oran.png` - Partial correlation: NDVI, Albedo → H
- `H7_second_peak_detection.png` - Monthly NDVI peaks (crop/phase detection)
- `crop_split_Oran_main_vs_early.png` - Main (Feb-Jul) vs Early (Oct-Jan) growth phase
- `same_period_benchmark_Apr_Jun.png` - Site comparison restricted to Apr-Jun
- `interannual_check.png` - Year-by-year EF/LE/ET robustness (2018/2019/2020)

### Summary Figures

- `final_summary_LE.png` - Combined LE comparison with statistical annotations
- `monthly_ndvi_seasonal.png` - Month-by-month NDVI boxplots
- `partial_correlation_summary.png` - NDVI and Rn independent contributions to LE/H

### Console Output

The script prints to stdout a summary of:
- Data loading diagnostics (TIMESTAMP parse success, unit confirmations)
- NDVI phenology thresholds and peak DOY for each site
- Cross-site comparisons (dry-canopy LE/EF/ET with p-values)
- H1 irrigation lag results
- H2 saturation slope comparisons
- H4 partial correlations
- H7 second-peak detection

## Common Issues

### Missing Data Files

```
❌ Oran EC: /path/to/oran.csv (not found or not configured)
```

**Solution**:
- Verify the file path exists and is readable
- Use `--oran /correct/path.csv` or export `BAKANPOSS_ORAN_EC=/correct/path.csv`

### NDVI File Not Found

The script looks in standard locations. If your NDVI file is elsewhere:

```bash
python3 run_analysis_C.py --ndvi /your/custom/path/NDVI.csv
```

### Timestamp Parsing Issues

The script diagnoses and handles multiple timestamp formats:
- Explicit formats: `%Y/%m/%d %H:%M:%S`, `%Y-%m-%d`, etc.
- Mixed format fallback for heterogeneous columns
- Julian date recovery: reconstructs from `year`, `Julian`, `Time_hours` columns

If warnings appear but processing continues, this is expected and handled.

### Unit Mismatches

The script auto-detects and reports:
- VPD in hPa vs kPa (converts to kPa)
- ET in mm/day vs mm/30min (uses ET_sum when available)
- SWC in % vs m³/m³ (converts to % if needed)

Check console output for unit confirmations like:
```
[Tarazona loader] 741 日分 (ET←ET_sum, VPD←VPD_kPa)
```

## Validation Without Running

To check configuration without executing analysis:

```bash
python3 run_analysis_C.py --dry-run
```

This validates all paths and reports errors without computing results.

## Direct Execution (Advanced)

If you prefer to run `analysis_C_v1.py` directly:

```python
import sys
from pathlib import Path
sys.path.insert(0, '.')

# Configure paths
import analysis_A_v9
analysis_A_v9.PATHS["oran_ec"] = Path("/your/oran.csv")
analysis_A_v9.PATHS["tara_ec"] = Path("/your/tarazona.csv")

import analysis_C_v1
analysis_C_v1.NDVI_APPEEARS_CSV = Path("/your/ndvi.csv")
analysis_C_v1.SAVE_DIR = Path("./output_analysis_C_v1")

# Execute the main analysis block
exec(open("analysis_C_v1.py").read())
```

## Performance Notes

- Full analysis runtime: ~10-30 minutes depending on hardware
- Memory usage: ~2-4 GB (keeps all intermediate DataFrames in memory)
- Disk space: ~500 MB for output figures and CSV exports

## References

- **Main Script**: `analysis_C_v1.py` (1500+ lines, comprehensive phenology-flux linkage analysis)
- **Data Loaders**: `data_loaders.py` (corrected Oran/Tarazona EC loaders with bug fixes)
- **Configuration**: `analysis_A_v9.py` (site coordinates, data paths, constants)
- **Report**: `reports/analysis_C_report.md` (detailed results and interpretation)
