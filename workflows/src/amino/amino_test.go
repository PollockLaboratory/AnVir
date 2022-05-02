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
	// Example for tests
	// 5  8  11 14 17 20 23
    // v  v  v  v  v  v  v
	// |||___|||___|||___|||
	//  0  1  2  3  4  5  6
	t.Run("case:within_gene;9-17", func(t *testing.T) {
		r1, r2 := CodonRange(5, 25, 9, 17)
		c1, c2 := 1, 4
		compare(r1, c1, t)
		compare(r2, c2, t)
	})
	t.Run("case:v_start_before_gene;4-17", func(t *testing.T) {
		r1, r2 := CodonRange(5, 25, 4, 17)
		c1, c2 := 0, 4
		compare(r1, c1, t)
		compare(r2, c2, t)
	})
	t.Run("case:v_end_after_gene;14-27", func(t *testing.T) {
		r1, r2 := CodonRange(5, 25, 14, 27)
		c1, c2 := 3, 6
		compare(r1, c1, t)
		compare(r2, c2, t)
	})
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
		result, _ := DNA2AminoAcid(dna, codon_table)
		compare(result, correct, t)
	})
	t.Run("not divisible by 3", func(t *testing.T) {
		dna := "TTTTTGTTCTTATGTTG"
		correct := "FLFLC"
		result, _ := DNA2AminoAcid(dna, codon_table)
		compare(result, correct, t)
	})
}

/// TODO
func TestGetGeneIntervals(t *testing.T) {

}

func TestGetChanges(t *testing.T) {
	t.Run("case_single_AA_changed", func(t *testing.T) {
		ref_aa := "F"
		alt_aa := "L"
		cstart := 5
		result := GetChanges(ref_aa, alt_aa, cstart)
		correct := []Change{{From: "F", To: "L", At: 5}}
		compare(result, correct, t)
	})
	t.Run("case_single_AA_no_change", func(t *testing.T) {
		ref_aa := "F"
		alt_aa := "F"
		cstart := 5
		result := GetChanges(ref_aa, alt_aa, cstart)
		correct := []Change{}
		compare(result, correct, t)
	})
	t.Run("case_del_all", func(t *testing.T) {
		ref_aa := "YES"
		alt_aa := ""
		cstart := 5
		result := GetChanges(ref_aa, alt_aa, cstart)
		correct := []Change{
			{From: "Y", To: "del", At: 5},
			{From: "E", To: "del", At: 6},
			{From: "S", To: "del", At: 7},
		}
		compare(result, correct, t)
	})
	t.Run("case_ins_all", func(t *testing.T) {
		
		ref_aa := ""
		alt_aa := "YES"
		cstart := 5
		result := GetChanges(ref_aa, alt_aa, cstart)
		correct := []Change{
			{From: "ins", To: "Y", At: 5},
			{From: "ins", To: "E", At: 5},
			{From: "ins", To: "S", At: 5},
		}
		compare(result, correct, t)
	})
	t.Run("compound1", func(t *testing.T) {
		// 56789
		// ARNDCQE
		// AB--RTE
		ref_aa := "ARNDCQE"
		alt_aa := "AB--RTE"
		cstart := 5
		result := GetChanges(ref_aa, alt_aa, cstart)
		correct := []Change{
			{From: "R", To: "B", At: 6},
			{From: "N", To: "del", At: 7},
			{From: "D", To: "del", At: 8},
			{From: "C", To: "R", At: 9},
			{From: "Q", To: "T", At: 10},
		}
		compare(result, correct, t)
	})
	t.Run("compound2", func(t *testing.T) {
		//         56--789
		ref_aa := "AB--RTE"
		alt_aa := "ARNDCQE"
		cstart := 5
		result := GetChanges(ref_aa, alt_aa, cstart)
		correct := []Change{
			{From: "B", To: "R", At: 6},
			{From: "ins", To: "N", At: 6},
			{From: "ins", To: "D", At: 6},
			{From: "R", To: "C", At: 7},
			{From: "T", To: "Q", At: 8},
		}
		compare(result, correct, t)
	})

}

func BenchmarkAnnotateChanges(t *testing.B) {
	ref_fasta := "../../data/NC_045512.2.fasta"
	vcf_file := "../../output/variants_genes.vcf"
	outfile := "test_data/out.vcf"
	genes_bed := "../../data/genes.bed.gz"
	codons_file := "../../data/dna_codon_table.tsv"
	AnnotateChanges(ref_fasta, vcf_file,
		genes_bed, codons_file, outfile)
}
