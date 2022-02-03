import sys
from collections import namedtuple
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

if __name__ == "__main__":
    ref = sys.argv[1]
    ## load the intersected snps/proteins
    snps = []
    for line in sys.stdin:
        A = line.rstrip().split('\t')
        snps.append(
            Variant(chrom=A[0], start=int(A[1]), end=int(A[2]),
                    id=A[3], ref=A[4], var=A[5], count=A[6],
                    prot_start=int(A[8]), prot_end=int(A[9]),
                    gene=A[10], prot_name=A[11]))


    ## annotate each snp
    ref_seq = RefSeq(ref)
    for S in snps:

        # snp position relative to protien start
        rel_pos = S.start - S.prot_start
        aa_index = rel_pos//3

        # TODO we don't really need to query the whole protein sequence
        # using aa index we could compute just the query of the changed codon
        ref_amino_seq = translate(ref_seq.query(S.prot_start, S.prot_end))
        changed_amino_seq = translate(ref_seq.query(
            S.prot_start, S.prot_end, change=((S.start, S.end), S.var)))

        if (ref:=ref_amino_seq[aa_index]) != (changed:=changed_amino_seq[aa_index]):
            aa_change = Change(aa_index+1, ref, changed)
            print('\t'.join(
                map(str, [S.chrom, S.start, S.end,
                          S.id, S.ref, S.var, S.count,
                          S.gene, S.prot_name, aa_change])))

