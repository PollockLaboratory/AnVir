// TODO put test data for all packages into a single directory
package amino_test

import (
	"fmt"
	"path/filepath"
	"reflect"
	"sort"
	"testing"

	. "annotation/amino"
	. "annotation/fastaseq"
	. "annotation/utils"
)

func compare[T any](result T, correct T, t *testing.T) {
	if !reflect.DeepEqual(result, correct) {
		t.Errorf("\ncorrect: %+v\nresult: %+v\n", correct, result)
	}
}

func TestRelativeCoords(t *testing.T) {
	// yes this is probably a stupid test :)
	// 0123456789
	// ABCDEFGHIJ
	// ^    ^
	// 10   15
	// relative
	// 15-10 = 5
	result := RelativeCoords(10, 15)
	correct := 5
	compare(result, correct, t)
}

func TestCodonPos(t *testing.T) {
	// abs:  25       33
	//        v       v  
	// rel:   012345678
	// gene:  |||___|||
	// codon:  0  1  2
	t.Run("start:25,pos:32", func(t *testing.T){
		correct := 2
		result := CodonPos(25, 32)
		compare(result, correct, t)
	})
	t.Run("start:25,pos:28", func(t *testing.T){
		correct := 1
		result := CodonPos(25, 28)
		compare(result, correct, t)
	})
	t.Run("start:25,pos:27", func(t *testing.T){
		correct := 0
		result := CodonPos(25, 27)
		compare(result, correct, t)
	})
}

func TestCodonRange(t *testing.T) {
	// 1  4  7  10 13 16 19
    // v  v  v  v  v  v  v
	// |||___|||___|||___|||
	//  0  1  2  3  4  5  6
	r1, r2 := CodonRange(1, 5, 13)
	c1, c2 := 1, 4
	compare(r1, c1, t)
	compare(r2, c2, t)
}

func TestCodonIndex2Genomic(t *testing.T) {
	// gstart
	// v
	// |||___|||
	// 123456789 <- genomic coordinates
	//  0  1  2  <- cindex
	// eg gstart = 1, cindex = 1
	// then start = 1 + 1*3 = 4
	correct := 4
	result := CodonIndex2Genomic(1, 1)
	compare(result, correct, t)
}

func TestGetCodonSeq(t *testing.T) {
	// example:
	//
	//    4 <-gstart (genomic start pos of gene)
	//    v
	// AAABBBCCCDDDEEE
	// FFFGGGHHHIIIJJJ
	//     ^     ^
	//     5     7 <- codon indices
	ref := LoadContiguousReference("test_data/test.fa")
	correct := "GGGHHHIII"
	result := CodonSeq(ref, 4, 5, 7)
	compare(correct, result, t)
}

func TestGetCodonTable(t *testing.T) {
	path, err := filepath.Abs("../../data/dna_codon_table.tsv")
	Check(err)
	table := GetCodonTable(path)

	// test by inspection only
	// for ease of comparison store keys in sorted order
	keys := make([]string, 0, len(table))
	for k := range table {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for _, k := range keys {
		fmt.Println(k, string(table[k]))
	}
}

func TestApplySNP(t *testing.T) {
	codon_seq := "ABC"
	
	t.Run("pos=0", func(t *testing.T) {
		correct := "XBC"
		pos := 0
		alt_seq := "X"
		result := ApplySNP(codon_seq, pos, alt_seq[0])
		compare(result, correct, t)
	})
	t.Run("pos=1", func(t *testing.T) {
		correct := "AXC"
		pos := 1
		alt_seq := "X"
		result := ApplySNP(codon_seq, pos, alt_seq[0])
		compare(result, correct, t)
	})
	t.Run("pos=2", func(t *testing.T) {
		correct := "ABX"
		pos := 2
		alt_seq := "X"
		result := ApplySNP(codon_seq, pos, alt_seq[0])
		compare(result, correct, t)
	})
}

func TestApplyDEL(t *testing.T) {
	codon_seq := "ABCDEFGHIJKLMNOPQRSTUVWX"
	t.Run("pos=0;length=3", func(t *testing.T) {
		correct := "DEFGHIJKLMNOPQRSTUVWX"
		pos := 0
		length := 3
		result := ApplyDEL(codon_seq, pos, length)
		compare(result, correct, t)
	})
	t.Run("pos=5;length=4", func(t *testing.T) {
		correct := "ABCDEJKLMNOPQRSTUVWX"
		pos := 5
		length := 4
		result := ApplyDEL(codon_seq, pos, length)
		compare(result, correct, t)
	})
	t.Run("pos=5;length=4", func(t *testing.T) {
		correct := "ABCDE"
		pos := 5
		length := 19
		result := ApplyDEL(codon_seq, pos, length)
		compare(result, correct, t)
	})
}

func TestApplyINS(t *testing.T) {
	codon_seq := "ABCDEFGHIJKLMNOPQRSTUVWX"
	t.Run("pos=0;hello", func(t *testing.T) {
		correct := "AhelloBCDEFGHIJKLMNOPQRSTUVWX"
		pos := 0
		result := ApplyINS(codon_seq, pos, "hello")
		compare(result, correct, t)
	})
	t.Run("pos=5;hello", func(t *testing.T) {
		correct := "ABCDEFhelloGHIJKLMNOPQRSTUVWX"
		pos := 5
		result := ApplyINS(codon_seq, pos, "hello")
		compare(result, correct, t)
	})
	t.Run("pos=23;hello", func(t *testing.T) {
		correct := "ABCDEFGHIJKLMNOPQRSTUVWXhello"
		pos := 23
		result := ApplyINS(codon_seq, pos, "hello")
		compare(result, correct, t)
	})
}

func TestApplyCompound(t *testing.T) {
	codon_seq := "ABCDEFGHIJKLMNOPQRSTUVWX"
	t.Run("pos=7;ref=HIJK;alt=#IJ#", func(t *testing.T) {
		correct := "ABCDEFG#IJ#LMNOPQRSTUVWX"
		pos := 6
		ref := "HIJK"
		alt := "#IJ#"
		result := ApplyCompound(codon_seq, pos, ref, alt)
		compare(result, correct, t)
	})
	t.Run("pos=7;ref=H--IJK;alt=H##IJ#", func(t *testing.T) {
		correct := "ABCDEFGH##IJ#LMNOPQRSTUVWX"
		pos := 6
		ref := "H--IJK"
		alt := "H##IJ#"
		result := ApplyCompound(codon_seq, pos, ref, alt)
		compare(result, correct, t)
	})
	t.Run("pos=7;ref=HIJKLMN;alt=HIJ--#-", func(t *testing.T) {
		correct := "ABCDEFGHIJ#OPQRSTUVWX"
		pos := 6
		ref := "HIJKLMN"
		alt := "HIJ--#-"
		result := ApplyCompound(codon_seq, pos, ref, alt)
		compare(result, correct, t)
	})
}

func TestDNA2AminoAcid(t *testing.T) {
	// small subset of codon table
	codon_table := make(map[string]byte)
	codon_table["TTT"] = 'F'
	codon_table["TTG"] = 'L'
	codon_table["TTC"] = 'F'
	codon_table["TTA"] = 'L'
	codon_table["TGT"] = 'C'
	codon_table["TGG"] = 'W'

	t.Run("divisible by 3", func(t *testing.T) {
		dna := "TTTTTGTTCTTATGTTGG"
		correct := "FLFLCW"
		result := DNA2AminoAcid(dna, codon_table)
		compare(result, correct, t)
	})
	t.Run("not divisible by 3", func(t *testing.T) {
		dna := "TTTTTGTTCTTATGTTG"
		correct := "FLFLC"
		result := DNA2AminoAcid(dna, codon_table)
		compare(result, correct, t)
	})
}

func TestAnnotateChanges(t *testing.T) {
	// TODO create small synthetic data
	ref_fasta := "../../data/NC_045512.2.fasta"
	vcf_file := "../classify_variants/test_data/benchmark.vcf"
	outfile := "test_data/out.vcf"
	AnnotateChanges(ref_fasta, vcf_file, outfile)
}
