package fastaseq_test
import (
	"reflect"
	"testing"
	"annotation/fastaseq"
)

// Test the WindowedReference struct/methods with the following fasta:
// test_ref_windows.fa (in same working directory)
func TestWindowedReference(t *testing.T)() {
	Ref := fastaseq.LoadWindowedReference("test_ref_windows.fa", 3)

	t.Run("Contig Parsed", func(t *testing.T)() {
		a := "CONTIG_NAME"
		b := Ref.Contig
		if a != b {
			t.Errorf("contig check: exptected %s but got %s", a, b)
		}
	})
	t.Run("Unique Kmer", func(t *testing.T)() {
		result := Ref.Query("ATC")
		correct := fastaseq.Interval{1, 3}
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
		correct := []fastaseq.Interval{{4, 6}, {10, 12}}
		if !reflect.DeepEqual(result, correct) {
			t.Errorf("\nQuery: 'GAA'\nresult=%v\ncorrect=%v,", result, correct)
		}
		result = Ref.Query("AAT")
		correct = []fastaseq.Interval{{5, 7}, {11, 13}}
	})
}
