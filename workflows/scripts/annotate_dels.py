import sys
from collections import namedtuple
from itertools import zip_longest
from copy import deepcopy
from Bio.Seq import translate
from ref_sequence import RefSeq

class Change:
    """
    Simplify printing amino acid or
    nucleotide changes with this 
    structure and repr method
    """
    def __init__(self, pos: int, ref: str, changed: str):
        self.pos = pos
        self.ref = ref
        self.changed = changed
    def __repr__(self):
        return f'{self.pos}{self.ref}>{self.changed}'


Variant = namedtuple('Variant', ['chrom', 'start', 'end', 'id', 'ref',
                                 'var', 'count', 'prot_start',
                                 'prot_end', 'gene', 'prot_name'])

if __name__ == '__main__':
    ref = sys.argv[1]

    ## load the intersected dels/annotations
    variants = []
    for line in sys.stdin:
        A = line.rstrip().split('\t')
        variants.append(
            Variant(chrom=A[0], start=int(A[1]), end=int(A[2]),
                    id=A[3], ref=A[4], var=A[5], count=A[6],
                    prot_start=int(A[8]), prot_end=int(A[9]),
                    gene=A[10], prot_name=A[11]))

    ## annotate each snp
    ref_seq = RefSeq(ref)
    for v in variants:
        ref_amino_seq = translate(ref_seq.query(v.prot_start, v.prot_end))
        alt_amino_seq = translate(
            ref_seq.query(start=v.prot_start, end=v.prot_end,
                          change=((v.start, v.end), 'del')))

        changes = []
        for i, (ref, alt) in enumerate(zip_longest(ref_amino_seq, alt_amino_seq, fillvalue='_')):
            if ref != alt:
                changes.append(Change(i+1, ref, alt))
        print('\t'.join(
            map(str, [v.chrom, v.start, v.end,
                      v.id, v.ref, v.var, v.count,
                      v.gene, v.prot_name, *changes])))
        
