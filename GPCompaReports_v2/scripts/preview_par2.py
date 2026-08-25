#!/usr/bin/env python3
"""Build the PAR2 report only (fast preview).

Loads the full data store (needed for cross-GPCR analyses like CFR), then
filters store.gpcr_ids to ['par2_human'] before generating individual
report pages so only PAR2 gets rendered.

Outputs to GPCompaReports_v2/output_par2_preview/ by default.
"""
from __future__ import annotations

import argparse
import shutil
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
V2_ROOT = HERE.parent
PROJECT_ROOT = V2_ROOT.parent
sys.path.insert(0, str(PROJECT_ROOT))

from jinja2 import Environment, FileSystemLoader

from GPCompaReports_v2.analysis.data_loader import GPCRDataStore
from GPCompaReports_v2.analysis.cross_gpcr_analysis import run_cross_gpcr_analysis
from GPCompaReports_v2.analysis.tm_domain_analysis import run_tm_domain_analysis
from GPCompaReports_v2.analysis.cfr_analysis import run_cfr_analysis
from GPCompaReports_v2.analysis.variant_correlation import run_variant_analysis
from GPCompaReports_v2.website.page_generators.gpcr_report_page import generate_all_reports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path,
                        default=V2_ROOT / 'output_par2_preview')
    parser.add_argument('--gpcr', default='par2_human',
                        help='Receptor id to render (default par2_human).')
    args = parser.parse_args()

    start = time.time()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for sd in ['reports', 'data', 'static/css', 'static/js']:
        (args.output_dir / sd).mkdir(parents=True, exist_ok=True)

    static_dir = V2_ROOT / 'static'
    for src in static_dir.rglob('*'):
        if src.is_file():
            rel = src.relative_to(static_dir)
            dst = args.output_dir / 'static' / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)

    print('Loading data store…')
    store = GPCRDataStore()
    store.load_all()

    if args.gpcr not in store.gpcr_ids:
        sys.exit(f"Receptor {args.gpcr!r} not found in batch.")

    print('Running cross-GPCR analysis (full set)…')
    analysis = {}
    analysis['cross_gpcr'] = run_cross_gpcr_analysis(store)
    analysis['tm_domain'] = run_tm_domain_analysis(store)
    analysis['cfr'] = run_cfr_analysis(store)
    cfr_table = analysis['cfr'].get('cfr_table')
    analysis['variant'] = run_variant_analysis(store, cfr_table)

    # Also seed the conservation cache for this receptor if a fetch script
    # exists on disk (we just check the cache file presence — preview build
    # will fall back to the variant-only conservation column otherwise).
    for cache_name in (f'conservation_{args.gpcr}.json', f'alphamissense_{args.gpcr}.json'):
        cache = V2_ROOT / 'output' / 'data' / cache_name
        if cache.exists():
            dst = args.output_dir / 'data' / cache.name
            shutil.copy2(cache, dst)
            print(f"Reused cache: {dst.name}")

    print(f"\nGenerating report for {args.gpcr}…")
    store.gpcr_ids = [args.gpcr]
    env = Environment(loader=FileSystemLoader(str(V2_ROOT / 'templates')))
    generate_all_reports(env, store, args.output_dir, analysis_results=analysis)

    print(f"\nDone in {time.time() - start:.1f}s")
    print(f"Open: {args.output_dir / 'reports' / (args.gpcr + '.html')}")


if __name__ == '__main__':
    main()
