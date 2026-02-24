#!/usr/bin/env python3
"""
GPCRdb AlphaFold-Multistate Model Downloader
Downloads active and inactive state models for human GPCRs
"""

import requests
from bs4 import BeautifulSoup
import time
import os
from pathlib import Path
import json
from datetime import datetime
import zipfile
import io

class GPCRdbDownloader:
    def __init__(self, output_dir="gpcrdb_alphafold_models"):
        self.output_dir = Path(output_dir)
        self.base_url = "https://gpcrdb.org"
        self.download_stats = {
            "successful": [],
            "failed": [],
            "total": 0
        }

        # Create organized directory structure
        self.active_dir = self.output_dir / "active"
        self.inactive_dir = self.output_dir / "inactive"
        self.active_dir.mkdir(parents=True, exist_ok=True)
        self.inactive_dir.mkdir(parents=True, exist_ok=True)

    def load_receptor_list(self):
        """
        Load receptor list from available files
        """
        # Try different filenames
        possible_files = [
            "all_receptors.txt",
            "class_a_receptors.txt", 
            "class_a_receptors_manual.txt"
        ]

        for filename in possible_files:
            filepath = Path(filename)
            if filepath.exists():
                with open(filepath, 'r') as f:
                    receptors = [line.strip() for line in f if line.strip()]
                print(f"✓ Loaded {len(receptors)} receptors from {filename}")
                return receptors

        return None

    def download_model(self, receptor, state):
        """
        Download a single model for a receptor in a specific state
        GPCRdb serves files as ZIP archives containing PDB files
        """
        url = f"{self.base_url}/structure/homology_models/{receptor}_{state}_full/download_pdb"

        # Determine output directory and filename
        if state == "active":
            output_file = self.active_dir / f"{receptor}_{state}.pdb"
        else:
            output_file = self.inactive_dir / f"{receptor}_{state}.pdb"

        # Skip if already exists
        if output_file.exists():
            file_size = output_file.stat().st_size / 1024
            if file_size > 1:  # More than 1 KB
                print(f"⊙ {receptor} ({state}): Already downloaded")
                return "exists"

        try:
            response = requests.get(url, timeout=30)

            if response.status_code == 200 and len(response.content) > 1000:
                # GPCRdb returns a ZIP file containing the PDB
                try:
                    # Read ZIP from memory
                    zip_buffer = io.BytesIO(response.content)
                    with zipfile.ZipFile(zip_buffer, 'r') as zip_ref:
                        # Get the PDB file (should be only one file in the archive)
                        pdb_files = [f for f in zip_ref.namelist() if f.endswith('.pdb')]

                        if pdb_files:
                            # Extract the PDB content
                            pdb_content = zip_ref.read(pdb_files[0])

                            # Save with our naming convention
                            with open(output_file, 'wb') as f:
                                f.write(pdb_content)

                            file_size = len(pdb_content) / 1024  # KB
                            print(f"✓ {receptor} ({state}): {file_size:.1f} KB")
                            return True
                        else:
                            print(f"✗ {receptor} ({state}): No PDB file in ZIP")
                            return False

                except zipfile.BadZipFile:
                    # If it's not a ZIP, try saving as-is (fallback for plain PDB)
                    with open(output_file, 'wb') as f:
                        f.write(response.content)

                    file_size = len(response.content) / 1024
                    print(f"✓ {receptor} ({state}): {file_size:.1f} KB (plain PDB)")
                    return True
            else:
                print(f"✗ {receptor} ({state}): No model available (status {response.status_code})")
                return False

        except Exception as e:
            print(f"✗ {receptor} ({state}): Error - {e}")
            return False

    def download_all(self, receptor_list, states=["active", "inactive"]):
        """
        Download all models for all receptors
        """
        total = len(receptor_list) * len(states)
        self.download_stats["total"] = total
        current = 0

        print(f"\nStarting download of {total} models...")
        print(f"Output directory: {self.output_dir.absolute()}")
        print(f"{'='*70}\n")

        start_time = time.time()

        for receptor in receptor_list:
            for state in states:
                current += 1
                print(f"[{current}/{total}] ", end="")

                success = self.download_model(receptor, state)

                if success == True:
                    self.download_stats["successful"].append(f"{receptor}_{state}")
                elif success == "exists":
                    self.download_stats["successful"].append(f"{receptor}_{state}")
                else:
                    self.download_stats["failed"].append(f"{receptor}_{state}")

                # Be nice to the server
                if success == True:  # Only delay for actual downloads
                    time.sleep(0.5)

        elapsed = time.time() - start_time
        self.save_report()
        self.print_summary(elapsed)

    def save_report(self):
        """
        Save download report as JSON
        """
        report = {
            "timestamp": datetime.now().isoformat(),
            "statistics": {
                "total_attempted": self.download_stats["total"],
                "successful": len(self.download_stats["successful"]),
                "failed": len(self.download_stats["failed"])
            },
            "successful_downloads": self.download_stats["successful"],
            "failed_downloads": self.download_stats["failed"]
        }

        report_file = self.output_dir / "download_report.json"
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)

        print(f"\nReport saved to: {report_file}")

    def print_summary(self, elapsed_time):
        """
        Print download summary
        """
        successful = len(self.download_stats["successful"])
        failed = len(self.download_stats["failed"])
        total = self.download_stats["total"]

        print("\n" + "="*70)
        print("DOWNLOAD SUMMARY")
        print("="*70)
        print(f"Total attempted:  {total}")
        print(f"Successful:       {successful} ({successful/total*100:.1f}%)")
        print(f"Failed:           {failed} ({failed/total*100:.1f}%)")
        print(f"Time elapsed:     {elapsed_time/60:.1f} minutes")
        print(f"\nFiles saved in:   {self.output_dir.absolute()}")
        print(f"  - Active:       {self.active_dir} ({len(list(self.active_dir.glob('*.pdb')))} files)")
        print(f"  - Inactive:     {self.inactive_dir} ({len(list(self.inactive_dir.glob('*.pdb')))} files)")
        print("="*70)


def main():
    print("="*70)
    print("GPCRdb AlphaFold-Multistate Model Downloader")
    print("="*70)
    print()

    # Initialize downloader
    downloader = GPCRdbDownloader(output_dir="gpcrdb_alphafold_models")

    # Load receptor list
    receptor_list = downloader.load_receptor_list()

    if not receptor_list:
        print("ERROR: No receptor list found!")
        print("\nPlease run: python3 get_class_a_list.py")
        print("This will create a receptor list file.\n")
        return

    print(f"Found {len(receptor_list)} receptors to download\n")

    # Confirm before starting
    print("This will download:")
    print(f"  • {len(receptor_list)} receptors")
    print(f"  • 2 states each (active + inactive)")
    print(f"  • Total: {len(receptor_list) * 2} PDB files")
    print(f"  • Estimated size: ~50-100 MB total")
    print(f"  • Estimated time: ~{len(receptor_list) * 2 * 0.5 / 60:.0f} minutes")

    response = input("\nContinue? [y/N]: ")
    if response.lower() != 'y':
        print("Cancelled.")
        return

    # Download all models
    downloader.download_all(receptor_list, states=["active", "inactive"])


if __name__ == "__main__":
    main()
