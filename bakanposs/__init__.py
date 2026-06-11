"""bakanposs — Tarazona/Oran EC and satellite ET analysis package.

Subpackages and modules:
  - bakanposs.loaders             : shared EC data loaders (Oran/Tarazona)
  - bakanposs.analysis_a          : Analysis A (LE recovery τ, deep-access test)
  - bakanposs.analysis_c.v1_legacy: Analysis C v1 (multipurpose NDVI explorer)
  - bakanposs.analysis_c.v2_phenology: Analysis C v2 (NDVI phenology validation)

Entry points (also runnable as `python -m bakanposs.<module>`):
  - python -m bakanposs.analysis_a
  - python -m bakanposs.analysis_c.v2_phenology
"""

__version__ = "0.2.0"
