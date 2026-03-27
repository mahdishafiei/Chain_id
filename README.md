# Chain ID — Antibody/Antigen Chain Identifier

Identify antibody and antigen chains from PDB/CIF structure files. For antibody chains, annotate CDR regions using IMGT, Kabat, or Chothia numbering schemes.

## Features

- Parse PDB and mmCIF structure files
- Classify each chain as **antibody** (heavy/light) or **antigen**
- Report chain length (number of residues)
- Identify CDR1, CDR2, CDR3 regions with residue-level mapping for antibody chains
- Display full chain sequences on demand
- Supports IMGT (default), Kabat, and Chothia numbering schemes

## Prerequisites

- Python 3.8+
- HMMER (required by ANARCI for antibody numbering)

```bash
# macOS
brew install hmmer

# Ubuntu/Debian
sudo apt install hmmer

# Conda
conda install -c bioconda hmmer
```

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Basic usage — summarize all chains
python chain_id.py structure.pdb

# Use mmCIF file
python chain_id.py structure.cif

# Choose numbering scheme (imgt, kabat, chothia)
python chain_id.py structure.pdb --scheme kabat

# Show full sequence for specific chains
python chain_id.py structure.pdb --show-seq A B

# Show full sequence for all chains
python chain_id.py structure.pdb --show-seq all

# Show CDR residue-level details
python chain_id.py structure.pdb --cdr-details
```

## Output Example

```
=== Chain Summary ===
Chain  Type              Length
─────────────────────────────────
A      Antibody (Heavy)  220
B      Antibody (Light)  214
C      Antigen           305

=== CDR Regions (IMGT) ===
Chain A (Heavy):
  CDR-H1: GFTFSSY (positions 27-38)
  CDR-H2: ISYDGSN (positions 56-65)
  CDR-H3: ARDLGYYDS (positions 105-117)

Chain B (Light):
  CDR-L1: QSISSY (positions 27-38)
  CDR-L2: AAS (positions 56-65)
  CDR-L3: QQSYSTPLT (positions 105-117)
```
