# GPCRdb AlphaFold-Multistate Model Downloader

Systematic download of AlphaFold-Multistate predicted structure models for Class A human GPCRs from GPCRdb.

## Overview

This project downloads AI-predicted GPCR structures in two activation states (ACTIVE and INACTIVE) from GPCRdb. These models were generated using a modified AlphaFold protocol by Heo & Feig (2022) and represent state-of-the-art predictions for GPCR structures.

### What You'll Get
- **~113 Class A human GPCRs**
- **2 states per receptor** (active + inactive)
- **~226 total PDB files**
- **Organized by activation state**
- **Total size**: ~50-100 MB
- **Download time**: ~5-10 minutes

## Quick Start

### Prerequisites

1. **Python 3.12+** with virtual environment
2. **Required packages**: `requests`, `beautifulsoup4`

Your environment is already set up at: `~/venv_gpcrdb`

### Step-by-Step Usage

```bash
# 1. Navigate to the project directory
cd "/media/yamir/TGGR_HDD/The Ultimate RRCS database project/alphafold_multistate_gpcr"

# 2. Activate virtual environment
source ~/venv_gpcrdb/bin/activate

# 3. (Optional) Get updated receptor list from GPCRdb
python3 get_class_a_list.py
# This attempts to scrape GPCRdb and creates class_a_receptors_manual.txt as backup

# 4. Download all models
python3 download_gpcrdb_models.py
# Will prompt for confirmation before starting
# Shows real-time progress with ✓ (success), ✗ (failed), ⊙ (exists)

# 5. Deactivate when done
deactivate
```

## Output Structure

```
gpcrdb_alphafold_models/
├── active/
│   ├── par2_human_active.pdb
│   ├── adrb2_human_active.pdb
│   ├── 5ht2a_human_active.pdb
│   └── ... (all active state models)
├── inactive/
│   ├── par2_human_inactive.pdb
│   ├── adrb2_human_inactive.pdb
│   ├── 5ht2a_human_inactive.pdb
│   └── ... (all inactive state models)
└── download_report.json
```

## Features

### Robust Download System
- ✓ **Resume capability**: Skips already downloaded files
- ✓ **Rate limiting**: 0.5s delay between downloads (server-friendly)
- ✓ **Error handling**: Continues downloading even if some fail
- ✓ **Progress tracking**: Real-time counter and status symbols
- ✓ **Comprehensive logging**: JSON report with all successes/failures

### File Organization
- Models organized by activation state (active/inactive directories)
- Clear naming convention: `{receptor}_{state}.pdb`
- Easy to find specific receptors or compare states

### Download Report (JSON)
```json
{
  "timestamp": "2026-02-02T11:30:00",
  "statistics": {
    "total_attempted": 226,
    "successful": 210,
    "failed": 16
  },
  "successful_downloads": ["par2_human_active", "par2_human_inactive", ...],
  "failed_downloads": ["gpr999_human_active", ...]
}
```

## Receptor List

The manual list (`class_a_receptors_manual.txt`) includes 113 Class A GPCRs organized by family:

- **Aminergic**: Adrenergic (α1A-D, α2A-C, β1-3), Serotonin (5-HT1A-F, 2A-C, 4, 5A, 6, 7), Dopamine (D1-5), Histamine (H1-4), Muscarinic (M1-5)
- **Opioid**: μ, δ, κ, nociceptin receptors
- **Chemokine**: CCR1-10, CXCR1-6, CX3CR1, CCRL2, ACKR3
- **Lipid**: Cannabinoid (CB1-2), Prostaglandin (EP1-4, FP, DP, DP2, IP, TP), PAF
- **Peptide**: Angiotensin (AT1-2), Bradykinin (B1-2), Endothelin (ETA-B), Melanocortin (MC1-5), Melatonin (MT1-2), NPY (Y1-2,Y4-5), Vasopressin/Oxytocin (OT, V1A-B, V2), PAR (PAR1-4)
- **Nucleotide**: Adenosine (A1, A2A, A2B, A3)
- **Orphan**: GPR3, GPR6, GPR12, GPR18-20, GPR22, GPR25, GPR27, GPR35, GPR37, GPR55, GPR84, GPR132

## Troubleshooting

### "No receptor list found" Error
**Solution**: Run `python3 get_class_a_list.py` first, or ensure `class_a_receptors_manual.txt` exists.

### Some downloads fail (404 errors)
**Expected**: Not all receptors have AlphaFold-Multistate models available on GPCRdb. Typical success rate: 85-95%.

### Download interrupted
**Solution**: Just run `python3 download_gpcrdb_models.py` again. It will skip existing files and resume from where it stopped.

### Network timeout
**Solution**: The script has 30-second timeout per file. If your connection is slow, you may get some timeouts. Just re-run the script to retry failed downloads.

### Virtual environment issues on external drive
If you need to recreate the venv on the external drive (which doesn't support symlinks):
```bash
python3 -m venv --copies venv
source venv/bin/activate
pip install requests beautifulsoup4
```

## Technical Details

### URL Pattern
```
https://gpcrdb.org/structure/homology_models/{receptor_name}_{state}_full/download_pdb
```

Example:
- `https://gpcrdb.org/structure/homology_models/par2_human_active_full/download_pdb`
- `https://gpcrdb.org/structure/homology_models/adrb2_human_inactive_full/download_pdb`

### Model Quality
- **Method**: AlphaFold-Multistate (Heo & Feig, 2022)
- **Training data**: Experimental structures up to July 2022
- **Accuracy**: Median RMSD 1.12Å (active), 1.41Å (inactive)
- **Version**: 2022 (newer than Zenodo's 2021 version)

### File Validation
- Minimum file size: 1 KB (rejects error pages)
- HTTP status: 200 (success) required
- Content type: PDB format text files
- Typical file size: 10-100 KB per model

## References

1. **AlphaFold-Multistate**: Heo, L. & Feig, M. (2022) "Multi-state modeling of G-protein coupled receptors at experimental accuracy" *Proteins* 90(11):1873-1885
2. **GPCRdb**: https://gpcrdb.org
3. **GPCRdb Paper**: Kooistra, A.J. et al. (2023) "GPCRdb in 2023" *Nucleic Acids Research*
4. **GitHub**: https://github.com/huhlim/alphafold-multistate

## Common Use Cases

### Download only specific receptors
Edit the receptor list file before running the download script:
```bash
# Create a custom list
echo "par2_human" > my_receptors.txt
echo "adrb2_human" >> my_receptors.txt

# Modify download script to read "my_receptors.txt"
```

### Compare active vs inactive states
Both states are in separate directories for easy comparison:
```bash
# Get file sizes to compare
ls -lh gpcrdb_alphafold_models/active/par2_human_active.pdb
ls -lh gpcrdb_alphafold_models/inactive/par2_human_inactive.pdb

# Use molecular viewers (PyMOL, Chimera, VMD) for structural comparison
```

### Find models for specific receptor families
```bash
# Find all serotonin receptors
ls gpcrdb_alphafold_models/active/5ht*

# Find all chemokine receptors
ls gpcrdb_alphafold_models/active/cc*
ls gpcrdb_alphafold_models/active/cxcr*
```

## Support

- **GPCRdb Help**: https://gpcrdb.org/documentation
- **Issues**: Check download_report.json for failed downloads
- **Questions**: Refer to GPCRdb API documentation: https://gpcrdb.org/services/

## License

The downloaded PDB files are provided by GPCRdb under their terms of use. These are AI-predicted models for research purposes.

---

**Last Updated**: February 2026
**Author**: Yamir (Petah Tikva, Israel)
**Project**: The Ultimate RRCS Database Project
