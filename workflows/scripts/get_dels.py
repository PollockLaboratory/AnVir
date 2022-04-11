import sys
from itertools import product
from ref_sequence import RefSeqWindows, RefSeq

variants = sys.argv[1]

ref_seq_windows = RefSeqWindows(reference=sys.argv[2], window_length=14)
ref_seq = RefSeq(reference=sys.argv[2])
contig = ref_seq_windows.contig

with open(variants) as f:
    for i in range(2): next(f) # skip header
    for line in f:
        A = line.rstrip().split()
        variant_id = A[0]
        count = A[2]

        # cols 5 (0-based) through N contain prev/next/deviants
        seqs = A[5:]

        # 13 deviants + prev/next
        if len(seqs) == 15:
            prev_seq = seqs[0]
            next_seq = seqs[-1]

            prev_coord = ref_seq_windows.ref_coords[prev_seq]
            next_coord = ref_seq_windows.ref_coords[next_seq]

            prev_next_product = list(product(prev_coord, next_coord))
            for pairs in prev_next_product:
                # pairs is a tuple of 2 (start,end) coordinates
                # pairs[0] := (start, end) of prev
                # pairs[1] := (start, end) of next
                # the region of the variant will begin 1 past the
                # end of prev and 1 before the beginning of next
                start = pairs[0][1] + 1
                end = pairs[1][0] - 1
                if start < end:
                    ref_sequence = ref_seq.query(start, end)
                    print('\t'.join(map(str, [contig, start, end,
                                            variant_id, ref_sequence,
                                            "DEL", count])))
