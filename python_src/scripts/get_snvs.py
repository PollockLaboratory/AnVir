import sys
from pprint import pprint
from itertools import product
from ref_sequence import RefSeqWindows, RefSeq

# the tab delimited file containing variants/variant ID
variants = sys.argv[1]

# load the reference and split into windows
ref_seq_windows = RefSeqWindows(reference=sys.argv[2], window_length=14)

# a bit redundant, but I want to load a data structure to query the
# reference with 1 based coords.
ref_seq = RefSeq(reference=sys.argv[2])

contig = ref_seq_windows.contig

with open(variants) as f:
    # the sequences start at col 5 (0-based)
    # the first seq is the "prev" and
    # the final seq is the "next"
    # the rest are the "deviants"
    for i in range(2): next(f) # skip header lines
    for line in f:

        variant_id = line.rstrip().split()[0]

        # cols 5 (0-based) through N contain prev/next/deviants
        seqs = line.rstrip().split()[5:]

        # not snv (prev+next+#deviants) - so for 14-mers 16 signifies snv
        if len(seqs) != 16: continue 

        # for an snv the variant will be the last nucleotide
        # at the first deviant and shift left as the window slides
        prev_seq = seqs[0]
        variant = seqs[1][-1]
        next_seq = seqs[-1]

        prev_coord = ref_seq_windows.ref_coords[prev_seq]
        next_coord = ref_seq_windows.ref_coords[next_seq]

        # get all possible pairs of (prev x next) seq coordinates
        # in the case of snvs the end of prev and the start of next
        # should have a difference of 2 (eg prev.end=1, next.start=3)
        prev_next_product = list(product(prev_coord, next_coord))
        for pairs in prev_next_product:
            # pairs is a tuple of 2 (start,end) coordinates
            # pairs[0] := (start, end) of prev
            # pairs[1] := (start, end) of next
            # there should be only one pair that satisfies this condition
            if pairs[1][0] - pairs[0][1] == 2:
                region = str(pairs[0][1] + 1)
                ref_nucleotide = ref_seq.query(int(region),int(region))
                print('\t'.join([contig, region, region,
                                 variant_id, ref_nucleotide, variant]))
