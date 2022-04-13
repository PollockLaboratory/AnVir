package fastaseq

import (
	// "fmt"
	"bufio"
	"os"
	"strings"
	"annotation/utils"
)
// =============================================================================
// This package provides a datastructure for querying the viral
// genome by k-mer sequence windows, or contiunously (essentially simple
// wrapper over a string).
// =============================================================================


// ============================================================================
/// Helper structs
// ============================================================================

type Interval struct {
	Start int
	End   int
}

// ============================================================================
/// Windowed Reference
// ============================================================================

// Data structure to store a single fasta record (viral genome) and enable
// queries keyed by kmer sliding widnows.
type WindowedReference struct {
	Contig      string
	K           int
	Kmer2coords map[string][]Interval
}

// Given a kmer sequence and the genomic positions, add to
// the data struture.  If location doesn't already exist,
// instantiate a slice keyed by the seq.
func (fs *WindowedReference)addseq(seq string, loc Interval) {
	if l, exists := fs.Kmer2coords[seq]; exists {
		fs.Kmer2coords[seq] = append(l, loc)
	} else {
		fs.Kmer2coords[seq] = []Interval{loc}
	}
}

// Query the data structure with a kmer seq.  I'm not enforcing that
// k-mer is the right size, however the query won't be successful if
// k isn't correct.  Returns start/end genomic coords (1-based, closed)
func (fs *WindowedReference)Query(seq string) []Interval {
	return fs.Kmer2coords[seq]
}

// Load reference genome (single fasta record) for kmer window queries.
// Serves as the constructor for WindowedReference objects.
func LoadWindowedReference(fasta_path string, k int) *WindowedReference {
	f, err := os.Open(fasta_path)
	utils.Check(err)
	defer f.Close()

	Ref := new(WindowedReference)
	Ref.Kmer2coords = make(map[string][]Interval)
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
		// Ref.Kmer2coords[seq[i:i+k]] = Interval{i+1, i+k}
		Ref.addseq(seq[i:i+k], Interval{i+1, i+k})
	}
	return Ref
}
