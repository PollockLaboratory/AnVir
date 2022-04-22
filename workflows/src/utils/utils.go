package utils

import (
	"bufio"
	"fmt"
	"os"
	"os/exec" // execute bash commands
	"strings"
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

// max of two ints because of course its not in the standard lib
func max(a int, b int) int {
	if a > b { return a }
	return b
}

func AlignSequences(ref string, alt string) (string, string){
	width := max(len(ref), len(alt))

	// use emboss needleman-wunsch implementation
	out, err := exec.Command("bash", "-c", fmt.Sprintf(
		`needle -asequence asis:%s -bsequence asis:%s -auto -awidth3 %d -stdout |
			grep -v '#' | grep -v '^$'`,
		ref, alt, width)).Output()
	Check(err)

	A := strings.Split(string(out), "\n")
	ref_alignment := strings.Fields(A[0])[2]
	alt_alignment := strings.Fields(A[2])[2]

	return ref_alignment, alt_alignment
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
