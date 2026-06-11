#!/usr/bin/env python3
"""Unified runner for Analysis C — dispatch to v1 or v2 by --version.

Equivalent to:
  v1: `python -m bakanposs.analysis_c.v1_legacy`
  v2: `python -m bakanposs.analysis_c.v2_phenology`
"""
import argparse
import runpy
import sys


def main():
    ap = argparse.ArgumentParser(description="Run Analysis C")
    ap.add_argument("--version", "-v", choices=["v1", "v2"], default="v2",
                    help="Which analysis variant to run (default: v2)")
    args = ap.parse_args()

    target = {
        "v1": "bakanposs.analysis_c.v1_legacy",
        "v2": "bakanposs.analysis_c.v2_phenology",
    }[args.version]
    print(f"[run_analysis_C] dispatching to {target}", file=sys.stderr)
    runpy.run_module(target, run_name="__main__")


if __name__ == "__main__":
    main()
