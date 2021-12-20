from ref_sequence import RefSeq
from itertools import zip_longest
import sys

def grouper(iterable, n, fillvalue=None):
    """
    Collect data into non-overlapping fixed-length chunks or blocks.
    Taken from itertools docs.
    """
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)



ref_seq = RefSeq(reference=sys.argv[1])

codons = grouper(ref_seq.query(29558,29674), 3, '-')


for codon in list(codons):
    print(''.join(codon))
# start codon
# ATG
# stop codon
# TAA, TAG, TGA

 # 266,805
 # 806,2719
 # 2720,8554
 # 8555,10054
 # 10055,10972
 # 10973,11842
 # 11843,12091
 # 12092,12685
 # 12686,13024
 # 13025,13441
 # 13442,13468
 # 13468,16236
 # 16237,18039
 # 18040,19620
 # 19621,20658
 # 20659,21552

 # 21563,25384
 # 25393,26220
 # 26245,26472
 # 26523,27191
 # 27202,27387
 # 27394,27759
 # 27756,27887
 # 27894,28259
 # 28274,29533
 # 29558,29674
