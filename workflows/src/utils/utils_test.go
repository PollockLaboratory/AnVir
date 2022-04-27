package utils_test

import (
	"fmt"
	. "annotation/utils"
	"reflect"
	"testing"
)

func TestInervalCartesianProduct(t *testing.T) {
	t.Run("2x2 product - no filter", func(t *testing.T) {
		A := []Interval{{Start:1, End:2}, {Start:3, End:4}}
		B := []Interval{{Start:5, End:6}, {Start:7, End:8}}

		result := IntervalCartesionProduct(A, B)
		correct := []Pair[Interval, Interval]{
			{Fst: A[0], Snd: B[0]},
			{Fst: A[0], Snd: B[1]},
			{Fst: A[1], Snd: B[0]},
			{Fst: A[1], Snd: B[1]},
		}
		if !reflect.DeepEqual(result, correct) {
			t.Errorf("\nCorrect:\n%+v\nResult:\n%+v\n", correct, result)
		}
	})
	t.Run("1x2 product - no filter", func(t *testing.T) {
		A := []Interval{{Start:1, End:2}}
		B := []Interval{{Start:5, End:6}, {Start:7, End:8}}

		result := IntervalCartesionProduct(A, B)
		correct := []Pair[Interval, Interval]{
			{Fst: A[0], Snd: B[0]},
			{Fst: A[0], Snd: B[1]},
		}
		if !reflect.DeepEqual(result, correct) {
			t.Errorf("\nCorrect:\n%+v\nResult:\n%+v\n", correct, result)
		}
	})
	t.Run("2x1 product - no filter", func(t *testing.T) {
		A := []Interval{{Start:1, End:2}, {Start:3, End:4}}
		B := []Interval{{Start:5, End:6}}

		result := IntervalCartesionProduct(A, B)
		correct := []Pair[Interval, Interval]{
			{Fst: A[0], Snd: B[0]},
			{Fst: A[1], Snd: B[0]},
		}
		if !reflect.DeepEqual(result, correct) {
			t.Logf("A: %+v\n", A)
			t.Logf("B: %+v\n", B)
			t.Errorf("\nCorrect:\n%+v\nResult:\n%+v\n", correct, result)
		}
	})
	// TODO bounds check test
	t.Run("2x1 product - yes filter", func(t *testing.T) {
		A := []Interval{{Start:1, End:2}, {Start:5, End:6}}
		B := []Interval{{Start:3, End:4}}

		result := IntervalCartesionProduct(A, B)
		correct := []Pair[Interval, Interval]{
			{Fst: A[0], Snd: B[0]},
		}
		if !reflect.DeepEqual(result, correct) {
			t.Logf("A: %+v\n", A)
			t.Logf("B: %+v\n", B)
			t.Errorf("\nCorrect:\n%+v\nResult:\n%+v\n", correct, result)
		}
	})
}

func TestSuffixPrefixOverlap(t *testing.T) {
	A := "QRSTUAAT"
	B := "AATUVXYZ"
	correct := "AAT"
	result := SuffixPrefixOverlap(A, B)
	if result != correct {
		t.Errorf("\ncorrect: %s\nresult: %s\n", correct, result)
	}
}

func BenchmarkAlignSequences(b *testing.B) {
	ref := "AAATGGGATCGATC"
	alt := "AAACCCGGGATCGATCG"
	fmt.Println(AlignSequences(ref, alt, true))

	ref = "AAARRREQGHILKM*"
	alt = "AAALLKEQGHIL*"
	fmt.Println(AlignSequences(ref, alt, false))
}

