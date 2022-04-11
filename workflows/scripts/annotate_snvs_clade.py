import sys
from pprint import pprint
from collections import defaultdict

snvs = sys.argv[1]
clades = sys.argv[2]

# its bad, but I love it
clades_table = defaultdict(lambda: defaultdict(list))

# map site -> alt -> list of clades
with open(clades) as f:
    next(f) # skip header
    for line in f:
        clade, _, site, alt = line.rstrip().split('\t')
        clades_table[site][alt].append(clade)

# pprint(clades_table)

# for a given snv get the list of clades that it fits into
with open(snvs) as f:
    for line in f:
        A = line.rstrip().split('\t')
        site = A[1]
        alt = A[5]
        clades = ':'.join(clades_table[site][alt])
        print('\t'.join([*A, clades]))
