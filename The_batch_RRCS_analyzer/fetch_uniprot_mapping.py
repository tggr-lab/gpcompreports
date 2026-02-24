#!/usr/bin/env python3
"""
Fetch UniProt IDs for all GPCRs from GPCRdb API
Creates uniprot_mapping.txt file
"""

import requests
import time
from pathlib import Path
import sys

def get_gpcr_list(data_dir):
    """Get list of all GPCR names from active directory."""
    active_dir = Path(data_dir) / "active"
    gpcr_names = []

    for pdb_file in active_dir.glob("*_active.pdb"):
        gpcr_name = pdb_file.stem.replace('_active', '')
        gpcr_names.append(gpcr_name)

    return sorted(gpcr_names)

def fetch_uniprot_from_gpcrdb(gpcr_name):
    """
    Fetch UniProt ID from GPCRdb API.

    Args:
        gpcr_name: GPCR name like "5ht1a_human"

    Returns:
        UniProt ID or None
    """
    # GPCRdb API endpoint for protein details
    url = f"https://gpcrdb.org/services/protein/{gpcr_name}"

    try:
        response = requests.get(url, timeout=10)

        if response.status_code == 200:
            data = response.json()

            # The UniProt ID is in the 'accession' field
            uniprot_id = data.get('accession')

            if uniprot_id:
                return uniprot_id
            else:
                print(f"  ⚠ No UniProt ID found in response for {gpcr_name}")
                return None
        else:
            print(f"  ✗ API error for {gpcr_name}: HTTP {response.status_code}")
            return None

    except requests.exceptions.Timeout:
        print(f"  ✗ Timeout for {gpcr_name}")
        return None
    except requests.exceptions.RequestException as e:
        print(f"  ✗ Request failed for {gpcr_name}: {str(e)}")
        return None
    except Exception as e:
        print(f"  ✗ Error for {gpcr_name}: {str(e)}")
        return None

def main():
    """Main function to fetch all UniProt mappings."""
    data_dir = "/media/yamir/TGGR_HDD/The Ultimate RRCS database project/alphafold_multistate_gpcr/gpcrdb_alphafold_models"
    output_file = "uniprot_mapping.txt"

    print("="*60)
    print("UniProt Mapping Fetcher for GPCRdb")
    print("="*60)

    # Get all GPCR names
    print("\nScanning data directory...")
    gpcr_list = get_gpcr_list(data_dir)
    print(f"Found {len(gpcr_list)} GPCRs")

    # Fetch UniProt IDs
    print("\nFetching UniProt IDs from GPCRdb API...")
    print("(This may take a few minutes with rate limiting)\n")

    mappings = {}
    success_count = 0
    failed_count = 0

    for idx, gpcr_name in enumerate(gpcr_list, 1):
        print(f"[{idx}/{len(gpcr_list)}] {gpcr_name}...", end=" ")
        sys.stdout.flush()

        uniprot_id = fetch_uniprot_from_gpcrdb(gpcr_name)

        if uniprot_id:
            mappings[uniprot_id] = gpcr_name
            print(f"✓ {uniprot_id}")
            success_count += 1
        else:
            print(f"✗ Failed")
            failed_count += 1

        # Rate limiting: wait 0.5 seconds between requests
        if idx < len(gpcr_list):
            time.sleep(0.5)

    # Save to file
    print(f"\n{'='*60}")
    print("Saving mappings to file...")

    with open(output_file, 'w') as f:
        f.write("# UniProt ID to GPCRdb ID mapping\n")
        f.write("# Format: UNIPROT_ID  GPCRDB_ID\n")
        f.write(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Total mappings: {len(mappings)}\n")
        f.write("#\n")

        for uniprot_id, gpcr_name in sorted(mappings.items()):
            f.write(f"{uniprot_id}\t{gpcr_name}\n")

    # Print summary
    print(f"{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Total GPCRs:        {len(gpcr_list)}")
    print(f"Successful:         {success_count}")
    print(f"Failed:             {failed_count}")
    print(f"Success rate:       {success_count/len(gpcr_list)*100:.1f}%")
    print(f"\nMapping saved to:   {output_file}")
    print(f"{'='*60}")

    return 0 if failed_count == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
