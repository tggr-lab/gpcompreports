#!/usr/bin/env python3
"""GPCompReports v2 — Master script to generate the complete static website.

Usage:
    python3 generate_site.py
    python3 generate_site.py --output /path/to/output
    python3 generate_site.py --limit 3   # test-build only 3 report pages

The v1 legacy site lives in ../GPCompReports/ and is independent of this
build.
"""

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from GPCompReports_v2.website.site_generator import SiteGenerator


def main():
    parser = argparse.ArgumentParser(description='Generate GPCompReports v2 static website')
    parser.add_argument('--output', type=str, default=None,
                        help='Output directory (default: GPCompReports_v2/output)')
    parser.add_argument('--batch-dir', type=str, default=None,
                        help='Batch analysis directory')
    parser.add_argument('--metadata', type=str, default=None,
                        help='Path to class_A_all.csv')
    parser.add_argument('--limit', type=int, default=None,
                        help='Only generate N individual reports (for testing)')
    args = parser.parse_args()

    start = time.time()

    generator = SiteGenerator(
        batch_dir=args.batch_dir,
        metadata_csv=args.metadata,
        output_dir=args.output,
        limit=args.limit,
    )
    generator.run()

    elapsed = time.time() - start
    print(f"\nTotal time: {elapsed:.1f}s")


if __name__ == '__main__':
    main()
