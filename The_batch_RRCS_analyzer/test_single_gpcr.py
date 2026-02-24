#!/usr/bin/env python3
"""
Test script to process a single GPCR
"""

import os
import sys
import logging
from pathlib import Path
from datetime import datetime

# Import the necessary modules
from single_gpcr_processor import SingleGPCRProcessor

def setup_test_logging():
    """Setup logging for test run."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s - %(message)s'
    )
    return logging.getLogger('TestSingleGPCR')

def create_output_dirs(base_dir):
    """Create output directory structure for test."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    dirs = {
        'root': base_dir / f"test_single_gpcr_{timestamp}",
        'reports': base_dir / f"test_single_gpcr_{timestamp}" / "reports",
        'csv_data': base_dir / f"test_single_gpcr_{timestamp}" / "csv_data",
        'plots': base_dir / f"test_single_gpcr_{timestamp}" / "plots",
        'logs': base_dir / f"test_single_gpcr_{timestamp}" / "logs",
        'summary': base_dir / f"test_single_gpcr_{timestamp}" / "summary",
        'rrcs_matrices': base_dir / f"test_single_gpcr_{timestamp}" / "rrcs_matrices",
        'variants': base_dir / f"test_single_gpcr_{timestamp}" / "variants"
    }

    for dir_path in dirs.values():
        dir_path.mkdir(parents=True, exist_ok=True)

    return dirs

def main():
    """Run a single GPCR test."""
    logger = setup_test_logging()

    # Define paths
    data_dir = Path("/media/yamir/TGGR_HDD/The Ultimate RRCS database project/alphafold_multistate_gpcr/gpcrdb_alphafold_models")
    output_dir = Path("./test_output")

    # Test GPCR
    gpcr_name = "5ht1a_human"

    logger.info(f"Testing with GPCR: {gpcr_name}")
    logger.info(f"Data directory: {data_dir}")

    # Check if PDB files exist
    active_pdb = data_dir / "active" / f"{gpcr_name}_active.pdb"
    inactive_pdb = data_dir / "inactive" / f"{gpcr_name}_inactive.pdb"

    if not active_pdb.exists():
        logger.error(f"Active PDB not found: {active_pdb}")
        return 1

    if not inactive_pdb.exists():
        logger.error(f"Inactive PDB not found: {inactive_pdb}")
        return 1

    logger.info(f"✓ Found active PDB: {active_pdb}")
    logger.info(f"✓ Found inactive PDB: {inactive_pdb}")

    # Create output directories
    output_dirs = create_output_dirs(output_dir)
    logger.info(f"Output directory: {output_dirs['root']}")

    # Prepare GPCR info
    gpcr_info = {
        'name': gpcr_name,
        'receptor': '5HT1A',
        'species': 'Human',
        'active_pdb': str(active_pdb),
        'inactive_pdb': str(inactive_pdb)
    }

    # Create processor and run analysis
    logger.info("\n" + "="*60)
    logger.info("Starting analysis...")
    logger.info("="*60 + "\n")

    try:
        processor = SingleGPCRProcessor(
            gpcr_info=gpcr_info,
            output_dirs=output_dirs,
            logger=logger
        )

        results = processor.run_complete_analysis()

        logger.info("\n" + "="*60)
        logger.info("ANALYSIS COMPLETE!")
        logger.info("="*60)
        logger.info(f"Summary:")
        logger.info(f"  - Total contacts: {results['summary']['total_contacts']}")
        logger.info(f"  - Significant changes: {results['summary']['significant_changes']}")
        logger.info(f"  - Variants found: {results['summary']['variants_found']}")
        logger.info(f"  - UniProt ID: {results['summary']['uniprot_id']}")
        logger.info(f"\nFiles generated: {len(results['files_generated'])}")
        for file in results['files_generated']:
            logger.info(f"  - {file}")

        logger.info(f"\n✓ Test completed successfully!")
        logger.info(f"✓ Output directory: {output_dirs['root']}")

        return 0

    except Exception as e:
        logger.error(f"\n✗ Test failed: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())
