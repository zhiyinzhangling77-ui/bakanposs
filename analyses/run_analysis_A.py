#!/usr/bin/env python3
"""Thin wrapper for Analysis A (v9) — invoke the module so its __main__ runs.

Equivalent to: `python -m bakanposs.analysis_a`
"""
import runpy

if __name__ == "__main__":
    runpy.run_module("bakanposs.analysis_a", run_name="__main__")
