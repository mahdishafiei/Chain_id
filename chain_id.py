#!/usr/bin/env python3
"""
Chain ID — Identify antibody/antigen chains from PDB/CIF files.

Classifies each chain as antibody (heavy/light) or antigen,
reports chain lengths, CDR regions, and optionally displays sequences.
"""

import argparse
import sys
from pathlib import Path

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import PPBuilder, is_aa
from Bio.Data.IUPACData import protein_letters_3to1

# Map 3-letter uppercase codes to 1-letter
_3to1 = {k.upper(): v for k, v in protein_letters_3to1.items()}


def three_to_one(resname):
    """Convert 3-letter amino acid code to 1-letter."""
    return _3to1[resname.strip().upper()]

try:
    from abnumber import Chain as AbChain, ChainParseError
except ImportError:
    print("Error: 'abnumber' package is required. Install with: pip install abnumber")
    sys.exit(1)


def extract_chains(filepath):
    """Parse a PDB or CIF file and return chain IDs with their sequences and residue info."""
    filepath = Path(filepath)
    suffix = filepath.suffix.lower()

    if suffix in ('.cif', '.mmcif'):
        parser = MMCIFParser(QUIET=True)
    elif suffix in ('.pdb', '.ent'):
        parser = PDBParser(QUIET=True)
    else:
        print(f"Error: Unsupported file format '{suffix}'. Use .pdb or .cif")
        sys.exit(1)

    structure = parser.get_structure('structure', str(filepath))
    model = structure[0]  # Use first model

    chains = {}
    for chain in model:
        chain_id = chain.get_id()
        residues = []
        sequence = []

        for residue in chain:
            if is_aa(residue, standard=True):
                het, resseq, icode = residue.get_id()
                if het.strip():  # Skip HETATM amino acids (modified residues)
                    continue
                try:
                    aa = three_to_one(residue.get_resname())
                except KeyError:
                    continue
                icode_str = icode.strip()
                res_label = f"{resseq}{icode_str}" if icode_str else str(resseq)
                residues.append({
                    'resseq': resseq,
                    'icode': icode_str,
                    'label': res_label,
                    'aa': aa,
                })
                sequence.append(aa)

        if sequence:
            chains[chain_id] = {
                'sequence': ''.join(sequence),
                'residues': residues,
                'length': len(sequence),
            }

    return chains


def classify_chain(sequence, scheme='imgt'):
    """Try to classify a sequence as antibody. Returns AbChain object or None."""
    try:
        ab = AbChain(sequence, scheme=scheme)
        return ab
    except ChainParseError:
        return None
    except Exception:
        return None


def get_chain_type_label(ab_chain):
    """Human-readable label for antibody chain type."""
    ctype = ab_chain.chain_type
    labels = {'H': 'Heavy', 'K': 'Light (Kappa)', 'L': 'Light (Lambda)'}
    return labels.get(ctype, ctype)


def get_cdr_residue_mapping(ab_chain, residues):
    """Map CDR regions to original PDB residue numbers.

    Returns a dict with CDR names as keys and lists of (pdb_label, scheme_pos, aa) as values.
    """
    regions = ab_chain.regions  # OrderedDict{region_name: OrderedDict{Position: aa}}

    # The variable domain sequence from abnumber
    vd_seq = str(ab_chain.seq)
    full_seq = ''.join(r['aa'] for r in residues)

    # Find where the variable domain starts in the full chain
    vd_start = full_seq.find(vd_seq)
    if vd_start == -1:
        vd_start = 0

    vd_residues = residues[vd_start:vd_start + len(vd_seq)]

    # Iterate through all regions, tracking position in variable domain sequence
    vd_idx = 0
    cdrs = {}
    for region_name, region_dict in regions.items():
        if 'CDR' not in region_name:
            # Still need to advance vd_idx for framework regions
            for pos, aa in region_dict.items():
                vd_idx += 1
            continue

        cdr_residues = []
        for pos, aa in region_dict.items():
            pdb_label = vd_residues[vd_idx]['label'] if vd_idx < len(vd_residues) else '?'
            cdr_residues.append((pdb_label, str(pos), aa))
            vd_idx += 1

        if cdr_residues:
            cdrs[region_name] = cdr_residues

    return cdrs


def print_summary(chain_data, scheme):
    """Print the chain summary table."""
    print(f"\n{'=' * 50}")
    print(f"  Chain Summary (numbering scheme: {scheme.upper()})")
    print(f"{'=' * 50}")
    print(f"{'Chain':<8} {'Type':<25} {'Length':>6}")
    print(f"{'─' * 42}")

    for cid, info in chain_data.items():
        ctype = info['type_label']
        length = info['length']
        print(f"{cid:<8} {ctype:<25} {length:>6}")

    print()


def print_cdr_summary(chain_data, scheme):
    """Print CDR regions for antibody chains."""
    ab_chains = {cid: info for cid, info in chain_data.items() if info['is_antibody']}
    if not ab_chains:
        return

    print(f"{'=' * 50}")
    print(f"  CDR Regions ({scheme.upper()})")
    print(f"{'=' * 50}")

    for cid, info in ab_chains.items():
        ab = info['ab_chain']
        ctype = get_chain_type_label(ab)
        short_type = 'H' if ab.chain_type == 'H' else 'L'

        print(f"\nChain {cid} ({ctype}):")
        print(f"  CDR-{short_type}1: {ab.cdr1_seq}")
        print(f"  CDR-{short_type}2: {ab.cdr2_seq}")
        print(f"  CDR-{short_type}3: {ab.cdr3_seq}")

    print()


def print_cdr_details(chain_data, scheme):
    """Print CDR regions with residue-level details."""
    ab_chains = {cid: info for cid, info in chain_data.items() if info['is_antibody']}
    if not ab_chains:
        return

    print(f"{'=' * 50}")
    print(f"  CDR Residue Details ({scheme.upper()})")
    print(f"{'=' * 50}")

    for cid, info in ab_chains.items():
        ab = info['ab_chain']
        cdrs = info['cdr_mapping']
        ctype = get_chain_type_label(ab)
        short_type = 'H' if ab.chain_type == 'H' else 'L'

        print(f"\nChain {cid} ({ctype}):")
        for region_name, res_list in cdrs.items():
            cdr_num = region_name.replace('CDR', '')
            display_name = f"CDR-{short_type}{cdr_num}"
            seq = ''.join(aa for _, _, aa in res_list)
            pdb_labels = [pdb_label for pdb_label, _, _ in res_list]
            pos_range = f"{pdb_labels[0]}-{pdb_labels[-1]}" if pdb_labels else "N/A"
            print(f"  {display_name}: {seq}  (PDB residues {pos_range})")
            for pdb_label, scheme_pos, aa in res_list:
                print(f"    PDB {pdb_label:>6} | {scheme_pos}: {aa}")

    print()


def print_sequence(chain_data, show_chains):
    """Print full sequences for requested chains."""
    if show_chains == ['all']:
        targets = list(chain_data.keys())
    else:
        targets = show_chains

    print(f"{'=' * 50}")
    print(f"  Chain Sequences")
    print(f"{'=' * 50}")

    for cid in targets:
        if cid not in chain_data:
            print(f"\nChain {cid}: not found in structure")
            continue
        info = chain_data[cid]
        seq = info['sequence']
        print(f"\nChain {cid} ({info['type_label']}, {info['length']} residues):")

        # Print sequence in lines of 60 with position markers
        for i in range(0, len(seq), 60):
            chunk = seq[i:i + 60]
            # Add spacing every 10 residues
            spaced = ' '.join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
            print(f"  {i + 1:>5}  {spaced}")

    print()


def main():
    parser = argparse.ArgumentParser(
        description='Identify antibody/antigen chains from PDB/CIF files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python chain_id.py structure.pdb
  python chain_id.py structure.cif --scheme kabat
  python chain_id.py structure.pdb --show-seq A B
  python chain_id.py structure.pdb --show-seq all --cdr-details
        """
    )
    parser.add_argument('input', help='PDB or CIF file path')
    parser.add_argument('--scheme', choices=['imgt', 'kabat', 'chothia'], default='imgt',
                        help='Numbering scheme for CDR annotation (default: imgt)')
    parser.add_argument('--show-seq', nargs='+', metavar='CHAIN',
                        help='Show full sequence for specified chains (use "all" for all)')
    parser.add_argument('--cdr-details', action='store_true',
                        help='Show residue-level CDR details')

    args = parser.parse_args()

    # Check file exists
    if not Path(args.input).exists():
        print(f"Error: File not found: {args.input}")
        sys.exit(1)

    # Parse structure
    print(f"Parsing {args.input}...")
    chains = extract_chains(args.input)

    if not chains:
        print("No protein chains found in the structure.")
        sys.exit(1)

    print(f"Found {len(chains)} chain(s). Classifying...")

    # Classify each chain
    chain_data = {}
    for cid, info in chains.items():
        ab = classify_chain(info['sequence'], scheme=args.scheme)
        is_ab = ab is not None

        if is_ab:
            type_label = f"Antibody ({get_chain_type_label(ab)})"
            cdr_mapping = get_cdr_residue_mapping(ab, info['residues'])
        else:
            type_label = "Antigen"
            cdr_mapping = {}

        chain_data[cid] = {
            'sequence': info['sequence'],
            'residues': info['residues'],
            'length': info['length'],
            'is_antibody': is_ab,
            'type_label': type_label,
            'ab_chain': ab,
            'cdr_mapping': cdr_mapping,
        }

    # Output
    print_summary(chain_data, args.scheme)
    print_cdr_summary(chain_data, args.scheme)

    if args.cdr_details:
        print_cdr_details(chain_data, args.scheme)

    if args.show_seq:
        print_sequence(chain_data, args.show_seq)


if __name__ == '__main__':
    main()
