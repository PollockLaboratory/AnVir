package fastaseq_test
import (
	"path/filepath"
	"reflect"
	"testing"
	"annotation/fastaseq"
	. "annotation/utils"
)

// ============================================================================
/// Windowed Ref
// ============================================================================

// Test the WindowedReference struct/methods with the following fasta:
// test_ref_windows.fa (in same working directory)
func TestWindowedReference(t *testing.T)() {
	Ref := fastaseq.LoadWindowedReference("test_ref_windows.fa", 3)

	t.Run("Contig Parsed", func(t *testing.T)() {
		correct := "CONTIG_NAME"
		result := Ref.Contig
		if correct != result {
			t.Errorf("contig check: exptected %s but got %s", correct, result)
		}
	})
	t.Run("Unique Kmer", func(t *testing.T)() {
		result := Ref.Query("ATC")
		correct := Interval{1, 3}
		if len(result) != 1 {
			t.Errorf(
				"Query must have only contain 1 interval.  %d found",
				len(result))
		} else if result[0] != correct {
			t.Errorf(
				"correct = %+v; result[0] = %+v", correct, result[0])
		}
	})
	t.Run("Non-Unique Kmer", func(t *testing.T)() {
		result := Ref.Query("GAA")
		correct := []Interval{{4, 6}, {10, 12}}
		if !reflect.DeepEqual(result, correct) {
			t.Errorf("\nQuery: 'GAA'\nresult=%v\ncorrect=%v,", result, correct)
		}
		result = Ref.Query("AAT")
		correct = []Interval{{5, 7}, {11, 13}}
	})
}

func BenchmarkWindowedReference(b *testing.B) {
	path, err := filepath.Abs("../../NC_045512.2.fasta")
	Check(err)
	// for i := 0; i < 100; i++ { // give me challenge
		fastaseq.LoadWindowedReference(path, 14)
	// }
}

// ============================================================================
/// Contiguous Ref
// ============================================================================
func TestContinuousReference(t *testing.T)() {
	Ref := fastaseq.LoadContiguousReference("test_ref_windows.fa")

	t.Run("Contig Parsed", func(t *testing.T)() {
		correct := "CONTIG_NAME"
		result := Ref.Contig
		if correct != result {
			t.Errorf("\ncontig check: exptected %s but got %s", correct, result)
		}
	})
	t.Run("Correct sequence", func(t *testing.T)() {
		correct := "ATCGAATTTGAATGTA"
		result := Ref.Seq
		if correct != result {
			t.Errorf("\ncorrect = %s\nresult = %s\n", correct, result)
		}
	})
	t.Run("Interval query", func(t *testing.T)() {
		correct := "TTTGAATGTA"
		result := Ref.Query(7, 16)
		if correct != result {
			t.Errorf("\ncorrect = %s\nresult = %s\n", correct, result)
		}
	})
	t.Run("Single nucleotide query", func(t *testing.T)() {
		correct := "G"
		result := Ref.Query(10, 10)
		if correct != result {
			t.Errorf("\ncorrect = %s\nresult = %s\n", correct, result)
		}
	})
	t.Run("Single nucleotide query at end of ref", func(t *testing.T)() {
		correct := "A"
		result := Ref.Query(16, 16)
		if correct != result {
			t.Errorf("\ncorrect = %s\nresult = %s\n", correct, result)
		}
	})
}
