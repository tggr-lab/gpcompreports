#!/usr/bin/env python3
"""
Fetch Class A GPCR list by scraping GPCRdb structure models page
"""

import requests
from bs4 import BeautifulSoup
import re
import json

def scrape_receptor_list():
    """
    Scrape the GPCRdb structure models page for all receptors
    """
    print("Scraping GPCRdb structure models page...")

    url = "https://gpcrdb.org/structure/homology_models"

    try:
        headers = {
            'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36'
        }
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()

        # Find all receptor names in the HTML
        receptors = set()

        # Method 1: Look for patterns in the entire HTML
        pattern = r'([a-z0-9]+_human)_(?:active|inactive)_full'
        matches = re.findall(pattern, response.text)
        receptors.update(matches)

        # Method 2: Parse with BeautifulSoup
        soup = BeautifulSoup(response.text, 'html.parser')

        # Look for links, buttons, or data attributes
        for element in soup.find_all(['a', 'button', 'div', 'option']):
            text = element.get_text() + str(element.get('href', '')) + str(element.get('value', ''))
            matches = re.findall(pattern, text)
            receptors.update(matches)

        # Method 3: Look in script tags for JSON data
        for script in soup.find_all('script'):
            if script.string:
                matches = re.findall(pattern, script.string)
                receptors.update(matches)

        if receptors:
            receptor_list = sorted(list(receptors))

            # Filter for Class A if we can identify them
            # Class A typically includes: adrenergic, serotonin, dopamine, opioid, etc.
            print(f"\nFound {len(receptor_list)} human receptors")

            # Save all receptors
            with open("all_receptors.txt", "w") as f:
                for receptor in receptor_list:
                    f.write(f"{receptor}\n")
            print(f"✓ Saved all receptors to: all_receptors.txt")

            # Save as JSON too
            with open("all_receptors.json", "w") as f:
                json.dump(receptor_list, f, indent=2)
            print(f"✓ Saved to: all_receptors.json")

            # Show sample
            print(f"\nFirst 20 receptors found:")
            for receptor in receptor_list[:20]:
                print(f"  - {receptor}")
            print(f"  ... and {len(receptor_list) - 20} more")

            return receptor_list
        else:
            print("\nNo receptors found. The page structure may have changed.")
            return None

    except Exception as e:
        print(f"Error: {e}")
        return None

def create_manual_list():
    """
    Create a comprehensive manual list of Class A receptors
    Based on GPCRdb classification
    """
    print("\nCreating comprehensive Class A receptor list...")

    class_a_receptors = [
        # Adenosine receptors
        "aa1r_human", "aa2ar_human", "aa2br_human", "aa3r_human",

        # Adrenergic receptors
        "adrb1_human", "adrb2_human", "adrb3_human",
        "ada1a_human", "ada1b_human", "ada1d_human",
        "ada2a_human", "ada2b_human", "ada2c_human",

        # Serotonin receptors
        "5ht1a_human", "5ht1b_human", "5ht1d_human", "5ht1e_human", "5ht1f_human",
        "5ht2a_human", "5ht2b_human", "5ht2c_human",
        "5ht4r_human", "5ht5a_human", "5ht6r_human", "5ht7r_human",

        # Dopamine receptors
        "drd1_human", "drd2_human", "drd3_human", "drd4_human", "drd5_human",

        # Histamine receptors
        "hrh1_human", "hrh2_human", "hrh3_human", "hrh4_human",

        # Acetylcholine receptors (muscarinic)
        "acm1_human", "acm2_human", "acm3_human", "acm4_human", "acm5_human",

        # Opioid receptors
        "oprm_human", "oprd_human", "oprk_human", "oprx_human",

        # Cannabinoid receptors
        "cnr1_human", "cnr2_human",

        # Chemokine receptors
        "ccr1_human", "ccr2_human", "ccr3_human", "ccr4_human", "ccr5_human",
        "ccr6_human", "ccr7_human", "ccr8_human", "ccr9_human", "ccr10_human",
        "cxcr1_human", "cxcr2_human", "cxcr3_human", "cxcr4_human", "cxcr5_human", "cxcr6_human",
        "cx3c1_human", "ccrl2_human", "ackr3_human",

        # Endothelin receptors
        "ednra_human", "ednrb_human",

        # Protease-activated receptors
        "par1_human", "par2_human", "par3_human", "par4_human",

        # Angiotensin receptors
        "agtr1_human", "agtr2_human",

        # Bradykinin receptors
        "bkrb1_human", "bkrb2_human",

        # Melanocortin receptors
        "mc1r_human", "mc2r_human", "mc3r_human", "mc4r_human", "mc5r_human",

        # Melatonin receptors
        "mtnr1a_human", "mtnr1b_human",

        # Neuropeptide Y receptors
        "npy1r_human", "npy2r_human", "npy4r_human", "npy5r_human",

        # Oxytocin and vasopressin receptors
        "oxyr_human", "v1ar_human", "v1br_human", "v2r_human",

        # Prostaglandin receptors
        "pe2r1_human", "pe2r2_human", "pe2r3_human", "pe2r4_human",
        "ptafr_human", "ptgdr_human", "ptgdr2_human", "ptgfr_human", "ptgir_human", "ptger_human",

        # Orphan receptors
        "gpr3_human", "gpr6_human", "gpr12_human", "gpr18_human", "gpr19_human",
        "gpr20_human", "gpr22_human", "gpr25_human", "gpr27_human", "gpr35_human",
        "gpr37_human", "gpr55_human", "gpr84_human", "gpr132_human",
    ]

    class_a_receptors.sort()

    with open("class_a_receptors_manual.txt", "w") as f:
        for receptor in class_a_receptors:
            f.write(f"{receptor}\n")

    print(f"✓ Created manual list with {len(class_a_receptors)} Class A receptors")
    print(f"✓ Saved to: class_a_receptors_manual.txt")

    return class_a_receptors

if __name__ == "__main__":
    print("="*70)
    print("GPCRdb Receptor List Fetcher")
    print("="*70)
    print()

    # Try scraping first
    receptor_list = scrape_receptor_list()

    # Create manual list as backup
    print()
    manual_list = create_manual_list()

    if receptor_list:
        print(f"\n{'='*70}")
        print("SUCCESS: Use 'all_receptors.txt' with the download script")
        print(f"{'='*70}")
    else:
        print(f"\n{'='*70}")
        print("Use 'class_a_receptors_manual.txt' with the download script")
        print("This contains a curated list of major Class A receptors")
        print(f"{'='*70}")
