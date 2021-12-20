import sys
from copy import deepcopy
from collections import namedtuple
from Bio.Seq import translate
from ref_sequence import RefSeq

"""
* applies snps to reference sequence TODO intersect snps with proteins first
* gets amino acid sequence of ref/snp'd protiens
* get for each protein, note position (1-based), orig amino acid,
  and changed amino acid (eg 1C>Q, 25N>B, ...)
"""

class Change:
    """
    Simplify printing amino acid changes with this structure and repr method
    """
    def __init__(self, pos: int, ref: str, changed: str):
        self.pos = pos
        self.ref = ref
        self.changed = changed
    def __repr__(self):
        return f'{self.pos}{self.ref}>{self.changed}'

if __name__ == '__main__':
    proteins = sys.argv[1]
    snps = sys.argv[2]
    ref = sys.argv[3]

    ### get snp positions/bases
    with open(snps) as f:
        changes = []
        for line in f:
            A = line.rstrip().split()
            changes.append((int(A[1]), A[3]))

    ### load ref sequence and apply snps
    ref_seq = RefSeq(ref)
    changed_seq = deepcopy(ref_seq)
    changed_seq.apply_changes(changes)

    ### get the amino acid changes and annotate at each protein region
    # for each prot region:
    #    translate both seqs to amino acids
    #    for each amino acid:
    #        if different add pos/changed seq to list
    with open(proteins) as f:
        for line in f:
            A = line.rstrip().split()
            start, end = int(A[1]), int(A[2])
            ref_amino_seq = translate(ref_seq.query(start, end))
            changed_amino_seq = translate(changed_seq.query(start, end))
            changes = []
            for i, (r, c) in enumerate(zip(ref_amino_seq, changed_amino_seq)):
                if r != c: changes.append(Change(i+1, r, c))
            if changes:
                print('\t'.join([line.rstrip(), *map(str, changes)]))
