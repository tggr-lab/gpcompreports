#!/usr/bin/env python3
"""GPCompReports — Master script to generate the complete static website.

Usage:
    python generate_site.py
    python generate_site.py --output /path/to/output
"""

import argparse
import sys
import time
from pathlib import Path

# Add parent to path for package imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from GPCompReports.website.site_generator import SiteGenerator


def main():
    parser = argparse.ArgumentParser(description='Generate GPCompReports static website')
    parser.add_argument('--output', type=str, default=None,
                        help='Output directory (default: GPCompReports/output)')
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
