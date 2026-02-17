# -*- coding: utf-8 -*-
"""
@created: Feb 13 2026
@modified: Feb 16 2026
@author: Yoann Pradat

    Institut Gustave Roussy
    U981 & Prism Center
    114 rue Edouard Vaillant, Villejuif, 94800 France

Anonymize a BAM file by setting the sequenced reads to the reference at the positions of known SNPs.
"""

import argparse
import sys
import pysam

def load_sites(sites_file):
    """
    Parses a tab-delimited file containing SNP positions.
    Expected format (bcftools query output): CHROM  POS  REF  [ALT]
    """
    sites_dict = {}
    print(f"Loading sites from {sites_file}...")
    n_sites = 0

    with open(sites_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("CHROM"):
                continue

            parts = line.split('\t')
            if len(parts) < 3:
                continue

            chrom = parts[0]
            try:
                # VCF is 1-based, Pysam is 0-based.
                pos_1based = int(parts[1])
                pos_0based = pos_1based - 1
            except ValueError:
                continue

            ref_base = parts[2]
            alt_base = parts[3] if len(parts) > 3 else "."

            # Define keys for both chr-prefixed and non-prefixed versions
            if chrom.startswith("chr"):
                keys_to_update = [(chrom, pos_0based), (chrom[3:], pos_0based)]
            else:
                keys_to_update = [(chrom, pos_0based), (f"chr{chrom}", pos_0based)]
            n_sites += 1

            for key in keys_to_update:
                if key in sites_dict:
                    existing_ref, existing_alt = sites_dict[key]
                    if existing_ref != ref_base:
                        pass # Ref mismatch context

                    current_alts = existing_alt.split(',')
                    if alt_base not in current_alts:
                        new_alt_str = existing_alt + "," + alt_base
                        sites_dict[key] = (ref_base, new_alt_str)
                else:
                    sites_dict[key] = (ref_base, alt_base)

    print(f"Loaded {n_sites} unique sites to sanitize.")
    return sites_dict

def sanitize_bam(input_bam, output_bam, sites_dict, preview_limit=0,
                 shift_coord=False, chrom=None, start_fasta=None):
    """
    Reads input BAM, reverts SNPs to reference, writes output.

    If shift_coord is True:
      - Assumes Input BAM coordinates are already shifted relative to start_fasta.
      - Adds (start_fasta - 1) to the read position to look up the SNP in sites_dict.
    """

    try:
        bam_in = pysam.AlignmentFile(input_bam, "rb")
    except ValueError:
        print(f"Error: Could not open {input_bam}. Is it a BAM file?")
        sys.exit(1)

    # 1. Coordinate Logic
    shift_offset = 0
    if shift_coord:
        if not (chrom and start_fasta):
            print("Error: --shift_coord requires --chrom and --start_fasta (to map back to genomic coords).")
            sys.exit(1)

        # If BAM is shifted, Read Pos 0 corresponds to Genomic Pos (START_FASTA - 1)
        shift_offset = start_fasta - 1
        print(f"Shift Coord Mode: Input BAM is shifted. Adding {shift_offset}bp to look up SNPs.")

    # Open output BAM (Keep header exactly as is from input)
    bam_out = pysam.AlignmentFile(output_bam, "wb", template=bam_in)

    print(f"Processing {input_bam} -> {output_bam}...")

    count_total = 0
    count_modified = 0
    preview_count = 0

    for read in bam_in:
        count_total += 1
        if read.is_unmapped:
            bam_out.write(read)
            continue

        aligned_pairs = read.get_aligned_pairs(matches_only=True)
        read_modified = False
        seq_list = list(read.query_sequence)
        original_seq_preview = None
        modifications_log = []

        if preview_count < preview_limit:
            original_seq_preview = read.query_sequence

        is_reverse_strand = read.is_reverse
        read_strand = "-" if is_reverse_strand else "+"

        for query_pos, ref_pos in aligned_pairs:
            if ref_pos is None or query_pos is None:
                continue

            # Calculate Genomic Coordinate
            if shift_coord:
                genomic_ref_pos = ref_pos + shift_offset
            else:
                genomic_ref_pos = ref_pos

            # Lookup SNP using genomic coordinate
            site_key = (read.reference_name, genomic_ref_pos)

            if site_key in sites_dict:
                ref_base_genomic, alt_base_genomic = sites_dict[site_key]
                target_base_read = ref_base_genomic

                current_base = seq_list[query_pos]
                if current_base != target_base_read:
                    seq_list[query_pos] = target_base_read
                    read_modified = True

                    if preview_count < preview_limit:
                        modifications_log.append({
                            'idx': query_pos,
                            'chrom': read.reference_name,
                            'pos': genomic_ref_pos,
                            'bam_pos': ref_pos,
                            'old': current_base,
                            'new': target_base_read,
                            'ref_genomic': ref_base_genomic,
                            'alt_genomic': alt_base_genomic,
                            'strand': read_strand
                        })

        if read_modified:
            count_modified += 1
            read.query_sequence = "".join(seq_list)

            try:
                read.set_tag('MD', None)
                read.set_tag('NM', None)
            except KeyError:
                pass

        # Print Preview
        if read_modified and preview_count < preview_limit:
            print("-" * 80)
            print(f"MODIFIED READ: {read.query_name}")
            disp_start = read.reference_start + 1
            print(f"Coords (BAM):  {read.reference_name}:{disp_start} (Strand: {read_strand})")

            if shift_coord:
                print(f"Coords (Gen):  {read.reference_name}:{disp_start + shift_offset}")

            print("Modifications:")
            for mod in modifications_log:
                print(f"  - Read Idx {mod['idx']}: Changed '{mod['old']}' to '{mod['new']}'")
                print(f"    [SNP Genomic: {mod['chrom']}:{mod['pos']+1} ({mod['ref_genomic']}/{mod['alt_genomic']})]")

            preview_count += 1

        bam_out.write(read)

    bam_in.close()
    bam_out.close()
    print(f"Done. Processed {count_total} reads. Modified {count_modified} reads.")


def main(args):
    # 1. Load Sites
    sites = load_sites(args.sites_file)

    # 2. Run Sanitization
    sanitize_bam(
        args.input_bam,
        args.output_bam,
        sites,
        preview_limit=args.preview,
        shift_coord=args.shift_coord,
        chrom=args.chrom,
        start_fasta=args.start_fasta # Passed correctly
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--input_bam", required=True)
    parser.add_argument("--output_bam", required=True)
    parser.add_argument("--sites_file", required=True)
    parser.add_argument("--preview", type=int, default=0)

    # Arguments for coordinate mapping
    parser.add_argument("--shift_coord", action="store_true")
    parser.add_argument("--chrom", type=str)
    # Changed from --start to --start_fasta for clarity
    parser.add_argument("--start_fasta", type=int, help="The genomic coordinate corresponding to position 0 in the BAM (used for offset).")

    args = parser.parse_args()

    main(args)
