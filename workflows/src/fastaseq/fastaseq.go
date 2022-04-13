package fastaseq

import (
	"bufio"
	"os"
	"strings"
	"annotation/utils"
)
/*
   =============================================================================
   This package provides a datastructure for querying the viral
   genome by k-mer sequence windows, or contiunously (essentially simple
   wrapper over a string).

   /// Sketch to load sliding windows into interval map
   * get fasta path, kmer size
   * declare kmer2coords mapping string to pair of ints
   * parse fasta into single string -> ref string
   * for i = 0 to len(ref) - k
       window = ref[i:i+k] (0-based half-open interval)
       kmer2coords[window] = interval{i+1, i+k} (1-based closed interval)
   =============================================================================
*/

type interval struct {
	start int
	end int
}

// Data structure to store a single fasta record (viral genome) and enable
// queries keyed by kmer sliding widnows.
type FastaSeqWindowed struct {
	Contig string
	K int
	Kmer2coords map[string]interval
}
func (this *FastaSeqWindowed)Query(seq string) (int, int) {
	i := this.Kmer2coords[seq]
	return i.start, i.end
}

// Load reference genome (single fasta record) for kmer window queries.
func LoadRefWindowed(fasta_path string, k int) *FastaSeqWindowed {
	f, err := os.Open(fasta_path)
	utils.Check(err)
	defer f.Close()

	Ref := new(FastaSeqWindowed)
	Ref.Kmer2coords = make(map[string]interval)
	scanner := bufio.NewScanner(f)
	var sb strings.Builder

	scanner.Scan()
	Ref.Contig = scanner.Text()[1:] // discard the '>'
	Ref.K = k

	for scanner.Scan() {
		sb.WriteString(scanner.Text())
	}
	seq := sb.String()
	
	for i := 0; i < len(seq) - k; i++ {
		// remember: sliced sequence is 0-based, half-open
		// genomic interval is 1-based, closed
		Ref.Kmer2coords[seq[i:i+k]] = interval{i+1, i+k}
	}
	return Ref
}
