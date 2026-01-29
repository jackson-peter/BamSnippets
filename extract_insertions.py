#!/usr/bin/env python3

import pysam
import sys

bam_path = sys.argv[1]

bam = pysam.AlignmentFile(bam_path, "rb")

# Print header
print("\t".join([
    "chrom",
    "pos",
    "read_name",
    "insertion_length",
    "insertion_sequence",
    "strand",
    "mapq"
]))

for read in bam.fetch(until_eof=True):
    if read.is_unmapped:
        continue
    if read.cigartuples is None:
        continue

    ref_pos = read.reference_start  # 0-based
    read_pos = 0

    for cigar_op, length in read.cigartuples:
        # CIGAR operations:
        # 0 = M, 1 = I, 2 = D, 3 = N, 4 = S, 5 = H, 6 = P, 7 = =, 8 = X

        if cigar_op == 0 or cigar_op == 7 or cigar_op == 8:
            # Match or mismatch
            ref_pos += length
            read_pos += length

        elif cigar_op == 1:
            # Insertion relative to reference
            insertion_seq = read.query_sequence[read_pos:read_pos + length]
            strand = "-" if read.is_reverse else "+"

            print("\t".join(map(str, [
                bam.get_reference_name(read.reference_id),
                ref_pos + 1,  # convert to 1-based
                read.query_name,
                length,
                insertion_seq,
                strand,
                read.mapping_quality
            ])))

            read_pos += length

        elif cigar_op == 2 or cigar_op == 3:
            # Deletion or skipped region
            ref_pos += length

        elif cigar_op == 4:
            # Soft clip
            read_pos += length

        elif cigar_op == 5:
            # Hard clip
            pass

bam.close()