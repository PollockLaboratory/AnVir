package utils

import (
	"bufio"
	"fmt"
	"os"

	"github.com/biogo/biogo/align"
	"github.com/biogo/biogo/align/matrix"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/seq/linear"
)

// ============================================================================
/// Structs/Related functions
// ============================================================================

// used to store pairs of objects
type Pair[T, U any] struct {
	Fst T
	Snd U
}

// used to store genomic intervals.
// yeah I know that its essentially a Pair[int, int]...
type Interval struct {
	Start int
	End   int
}

// Given two lists of intervals get all possible combinations between
// A and B, discarding combinations if B occurs before A.
func IntervalCartesionProduct(A []Interval, B []Interval) []Pair[Interval, Interval]{
	// we expect len of A/B to be small
	intervals := make([]Pair[Interval, Interval], 0, len(A)*len(B))

	for _, a := range A {
		for _, b := range B {
			if b.Start > a.End {
				intervals = append(intervals, Pair[Interval, Interval]{a, b})
			}
		}
	}
	return intervals
}


// ============================================================================
/// Functions
// ============================================================================

// quick error checking pattern
func Check(e error) {
	if e != nil {
		panic(e)
	}
}


// I wonder what this function does?
func LineCount(filename string) int {
	f, err := os.Open(filename)
	Check(err)
	defer f.Close()
	scanner := bufio.NewScanner(f)
	var i int
	for i = 0; scanner.Scan(); i++ {}
	return i
}


// min/max of two ints because of course its not in the standard lib
func Max(a int, b int) int {
	if a > b { return a }
	return b
}
func Min(a int, b int) int {
	if a < b { return a }
	return b
}

// needlman wunch alignment for DNA or amino acid sequence
// (specified by dna flag). We expect sequences to be short.
func AlignSequences(ref string, alt string, dna bool) (string, string, error){

	// prepare sequences
	aseq := &linear.Seq{Seq: alphabet.BytesToLetters([]byte(ref))}
	bseq := &linear.Seq{Seq: alphabet.BytesToLetters([]byte(alt))}
	if dna {
		aseq.Alpha = alphabet.DNAredundant
		bseq.Alpha = alphabet.DNAredundant
	} else {
		aseq.Alpha = alphabet.Protein
		bseq.Alpha = alphabet.Protein
	}

	// init the similarity matrix and linear
	// component of the affine gap penalty
	var mat [][]int
	if dna {
		mat = matrix.NUC_4_4
	} else {
		mat = matrix.BLOSUM65
	}
	const lin_gap = -1
	for i := range mat {
		mat[i][0] = lin_gap
	}
	for j := range mat[0] {
		mat[0][j] = lin_gap
	}

	needle := align.NWAffine{Matrix: mat, GapOpen: -10}
	aln, err := needle.Align(aseq, bseq)

	if err != nil {
		return "", "", err
	}

	fa := align.Format(aseq, bseq, aln, '-')
	ref_alignment := fmt.Sprintf("%s", fa[0])
	alt_alignment := fmt.Sprintf("%s", fa[1])

	return ref_alignment, alt_alignment, nil
}

func SuffixPrefixOverlap(A string, B string) string {
	// assuming A, B are the same length
	// Find the max overlapping suffix of A and prefix of B.
	n := len(A)
	for i := range(A) {
		if A[i:] == B[:n-i] {
			return A[i:]
		}
	}
	return ""
}
