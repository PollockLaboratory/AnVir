package fastaseq_test
import (
	// "fmt"
	// "bufio"
	// "os"
	"testing"
	"annotation/fastaseq"
)

/*
   Test the LoadRefWindowed function with the following input
   test_ref_windows.fa
   >CONTIG_NAME
   ATCGAATT
   TGAATGTA
 */
func TestLoadRefWindows(t *testing.T)() {
	Ref := fastaseq.LoadRefWindowed("test_ref_windows.fa", 3)
	
	// Test contig field
	a := "CONTIG_NAME"
	b := Ref.Contig
	if a != b {
		t.Errorf("contig check: exptected %s but got %s", a, b)
	}

	// Test ATC: unique kmer
	a1, a2 := 1, 3
	b1, b2 := Ref.Query()
	// TODO Need to make the query return a list of possible regions? yes


	// Test GAA: non-unique kmer
}
