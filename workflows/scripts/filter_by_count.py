
import sys
from collections import defaultdict

class Value:
    def __init__(self, id: str, ref: str, var: str, count: int):
        self.id = id
        self.ref = ref
        self.var = var
        self.count = count
    def __lt__(self, other):
        """perform lt comp with the count fields of operands"""
        return self.count < other.count


bed = sys.argv[1]

# load the region as keys and the rest as values
regions = defaultdict(list)
with open(bed) as f:
    for line in f:
        A = line.rstrip().split()
        R = '\t'.join(A[:3])
        regions[R].append(Value(id=A[3], ref=A[4], var=A[5], count=int(A[6])))

# for each region, get the key with the max count
for r, v in regions.items():
    max_v = max(v)
    print('\t'.join([r, max_v.id, max_v.ref,
                     max_v.var, str(max_v.count)]))
