import sys
from itertools import islice
from collections import defaultdict, deque
from pprint import pprint
from Bio.Seq import translate

## Functions
###############################################################################
def sliding_window(iterable, n):
    """
    sliding_window('ABCDEFG', 4) -> ABCD BCDE CDEF DEFG
    Taken from python itertools example docs
    """
    it = iter(iterable)
    window = deque(islice(it, n), maxlen=n)
    if len(window) == n:
        yield tuple(window)
    for x in it:
        window.append(x)
        yield tuple(window)

def grouper(iterable, n, fillvalue=None):
    """
    Collect data into non-overlapping fixed-length chunks or blocks.
    Taken from itertools example docs.
    """
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


## Classes
###############################################################################
class RefSeqWindows:
    """
    Used to make kmer sequence queries with a sliding window acrross the ref.
    A window size of 14 (default) will yeild few non unique sequences.
    """
    def __init__(self, reference: str, window_length=14,):
        # NOTE: i'm assuming that the ref genome is just one contiguous sequence
        # get sequences in a sliding window along with the coordinates
        fasta = [f.rstrip() for f in open(reference).readlines()[1:]]
        fasta = ''.join(fasta)
        coordinates = list(range(1, len(fasta)+1)) # 1-based indexing
        fasta = sliding_window(fasta, window_length)
        coordinates = sliding_window(coordinates, window_length)

        # create a seq -> list(coordinates[start, end]) mapping (can be non unique)
        self.ref_coords = defaultdict(list)
        for f, c in zip(fasta, coordinates):
            seq = ''.join([i for i in f if i])
            coords = ((x:=[i for i in c if i])[0], x[-1])
            self.ref_coords[seq].append(coords)
        self.contig = open(reference).readline().split()[0][1:]

    def print_non_unique(self) -> None:
        """
        look at the non unique kmers and where they are in the ref
        """
        for key in self.ref_coords:
            if len(self.ref_coords[key])>1:
                print(key, self.ref_coords[key])

class RefSeq:
    """
    Used to make interval (1-based) queries to retrieve sequence from the reference
    """
    def __init__(self, reference: str):
        # NOTE: i'm assuming that the ref genome is just one contiguous sequence
        # get sequences in a sliding window along with the coordinates
        self.fasta = [f.rstrip() for f in open(reference).readlines()[1:]]
        self.fasta = ''.join(self.fasta)

    def query(self, start: int, end: int,
              change: tuple[tuple[int, int], str]=None) -> str:
        """
        make a 1-based closed interval query of the fasta.
        optionally apply change (closed interval) to returned
        query sequence (non-mutating).
        """

        if change:
            f = list(self.fasta[start-1:end])
            cs = change[0][0] - start
            ce = change[0][1] - start
            if len(change[1]) != (ce - cs + 1):
                raise ValueError("change interval != change string")
            f[cs:ce+1] = change[1]
            return ''.join(f)
        return self.fasta[start-1:end]
            

    def apply_changes(self, changes: list[tuple[int, str]]) -> None:
        """
        **This method mutates the reference data structure**
        given a list of genomic (1-based) positions
        and their changed bases, apply changes accross the ref
        """
        f = list(self.fasta)
        for pos, base in changes:
            f[pos-1] = base
        self.fasta = ''.join(f)
            

if __name__ == '__main__':
    ref_seq = RefSeqWindows(reference=sys.argv[1], window_length=14)
    ref_seq.print_non_unique()
