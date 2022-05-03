// TODO put test data for all packages into a single directory
package amino_test

import (
	"fmt"
	"path/filepath"
	"reflect"
	"sort"
	"testing"

	. "annotation/amino"
	"annotation/fastaseq"
	. "annotation/fastaseq"
	. "annotation/utils"
	"annotation/vcf"
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
		correct := []Change{{From: "F", To: "L", At: 6}}
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
			{From: "Y", To: "del", At: 6},
			{From: "E", To: "del", At: 7},
			{From: "S", To: "del", At: 8},
		}
		compare(result, correct, t)
	})
	t.Run("case_ins_all", func(t *testing.T) {
		ref_aa := ""
		alt_aa := "YES"
		cstart := 5
		result := GetChanges(ref_aa, alt_aa, cstart)
		correct := []Change{
			{From: "ins", To: "Y", At: 6},
			{From: "ins", To: "E", At: 6},
			{From: "ins", To: "S", At: 6},
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
			{From: "R", To: "B", At: 7},
			{From: "N", To: "del", At: 8},
			{From: "D", To: "del", At: 9},
			{From: "C", To: "R", At: 10},
			{From: "Q", To: "T", At: 11},
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
			{From: "B", To: "R", At: 7},
			{From: "ins", To: "N", At: 7},
			{From: "ins", To: "D", At: 7},
			{From: "R", To: "C", At: 8},
			{From: "T", To: "Q", At: 9},
		}
		compare(result, correct, t)
	})
}

func TestAminoAcidChanges(t *testing.T) {
	
	ref_fasta := "../../data/NC_045512.2.fasta"
	genes_bed := "../../data/genes.bed.gz"
	codons_file := "../../data/dna_codon_table.tsv"

	ref := fastaseq.LoadContiguousReference(ref_fasta)
	gene_intervals := GetGeneIntervals(genes_bed)
	codon_table := GetCodonTable(codons_file)

	// test case vcf records
	SNP1 := "NC_045512.2	24499	3357	T	C	.	.	VARTYPE=SNP;END=24499;COUNT=3585;KMERS=AGTGTTTTAAATGA,GTGTTTTAAATGAC,TGTTTTAAATGACA,GTTTTAAATGACAT,TTTTAAATGACATC,TTTAAATGACATCC,TTAAATGACATCCT,TAAATGACATCCTT,AAATGACATCCTTT,AATGACATCCTTTC,ATGACATCCTTTCA,TGACATCCTTTCAC,GACATCCTTTCACG,ACATCCTTTCACGT,CATCCTTTCACGTC,ATCCTTTCACGTCT;GENE=S"
	SNP2 := "NC_045512.2	22317	3361	G	T	.	.	VARTYPE=SNP;END=22317;COUNT=1281;KMERS=TTATTTGACTCCTG,TATTTGACTCCTGT,ATTTGACTCCTGTT,TTTGACTCCTGTTG,TTGACTCCTGTTGA,TGACTCCTGTTGAT,GACTCCTGTTGATT,ACTCCTGTTGATTC,CTCCTGTTGATTCT,TCCTGTTGATTCTT,CCTGTTGATTCTTC,CTGTTGATTCTTCT,TGTTGATTCTTCTT,GTTGATTCTTCTTC,TTGATTCTTCTTCA,TGATTCTTCTTCAG;GENE=S"
	DEL := "NC_045512.2	26158	3465	GTTA	DEL	.	.	VARTYPE=DEL;END=26161;COUNT=21688;KMERS=GTTCATCCGGAGTT,TTCATCCGGAGTTA,TCATCCGGAGTTAT,CATCCGGAGTTATC,ATCCGGAGTTATCC,TCCGGAGTTATCCA,CCGGAGTTATCCAG,CGGAGTTATCCAGT,GGAGTTATCCAGTA,GAGTTATCCAGTAA,AGTTATCCAGTAAT,GTTATCCAGTAATG,TTATCCAGTAATGG,TATCCAGTAATGGA,ATCCAGTAATGGAA;GENE=ORF3a"
	INS := "NC_045512.2	22205	3505	INS	CGGCAGGCT	.	.	VARTYPE=INS;END=22206;COUNT=2149;KMERS=TAATTTAGTGCGTG,AATTTAGTGCGTGC,ATTTAGTGCGTGCG,TTTAGTGCGTGCGG,TTAGTGCGTGCGGC,TAGTGCGTGCGGCA,AGTGCGTGCGGCAG,GTGCGTGCGGCAGG,TGCGTGCGGCAGGC,GCGTGCGGCAGGCT,CGTGCGGCAGGCTA,GTGCGGCAGGCTAT,TGCGGCAGGCTATC,GCGGCAGGCTATCT,CGGCAGGCTATCTC,GGCAGGCTATCTCC,GCAGGCTATCTCCC,CAGGCTATCTCCCT,AGGCTATCTCCCTC,GGCTATCTCCCTCA,GCTATCTCCCTCAG,CTATCTCCCTCAGG,TATCTCCCTCAGGG,ATCTCCCTCAGGGT;GENE=S"
	COMPOUND_DEL := "NC_045512.2	28273	108	ATGTCTGAT	-TGTCTCTA	.	.	VARTYPE=COMPOUND;END=28283;COUNT=1060102;KMERS=ACGAACAAACTAAA,CGAACAAACTAAAT,GAACAAACTAAATG,AACAAACTAAATGT,ACAAACTAAATGTC,CAAACTAAATGTCT,AAACTAAATGTCTC,AACTAAATGTCTCT,ACTAAATGTCTCTA,CTAAATGTCTCTAA,TAAATGTCTCTAAA,AAATGTCTCTAAAT,AATGTCTCTAAATG,ATGTCTCTAAATGG,TGTCTCTAAATGGA,GTCTCTAAATGGAC,TCTCTAAATGGACC,CTCTAAATGGACCC,TCTAAATGGACCCC,CTAAATGGACCCCA,TAAATGGACCCCAA,AAATGGACCCCAAA,AATGGACCCCAAAA;GENE=N"
	COMPOUND_INS := "NC_045512.2	11082	5381	G---	TTTT	.	.	VARTYPE=COMPOUND;END=11084;COUNT=445;KMERS=TTGTTCTTTTTTTT,TGTTCTTTTTTTTT,GTTCTTTTTTTTTT,TTCTTTTTTTTTTT,TCTTTTTTTTTTTT,CTTTTTTTTTTTTT,TTTTTTTTTTTTTA,TTTTTTTTTTTTAT,TTTTTTTTTTTATG,TTTTTTTTTTATGA,TTTTTTTTTATGAA,TTTTTTTTATGAAA,TTTTTTTATGAAAA,TTTTTTATGAAAAT,TTTTTATGAAAATG,TTTTATGAAAATGC,TTTATGAAAATGCC,TTATGAAAATGCCT,TATGAAAATGCCTT;GENE=ORF1a"
	COMPOUND_SNP := "NC_045512.2	28874	3413	GCAGTAGGGGAAC	TCAGTAGGGGAAT	.	.	VARTYPE=COMPOUND;END=28888;COUNT=2809;KMERS=TTCAACTCCAGGCA,TCAACTCCAGGCAT,CAACTCCAGGCATC,AACTCCAGGCATCA,ACTCCAGGCATCAG,CTCCAGGCATCAGT,TCCAGGCATCAGTA,CCAGGCATCAGTAG,CAGGCATCAGTAGG,AGGCATCAGTAGGG,GGCATCAGTAGGGG,GCATCAGTAGGGGA,CATCAGTAGGGGAA,ATCAGTAGGGGAAT,TCAGTAGGGGAATT,CAGTAGGGGAATTT,AGTAGGGGAATTTC,GTAGGGGAATTTCT,TAGGGGAATTTCTC,AGGGGAATTTCTCC,GGGGAATTTCTCCT,GGGAATTTCTCCTG,GGAATTTCTCCTGC,GAATTTCTCCTGCT,AATTTCTCCTGCTA,ATTTCTCCTGCTAG,TTTCTCCTGCTAGA,TTCTCCTGCTAGAA;GENE=N"

	t.Run("SNP1", func(t *testing.T) {
		// SNP occurs at pos 24499 in gene S.
		// S starts at 21563 (1-based).
		// Codon offset is (24499-21563)/3 = 978.  
		//
		// Relative offset pos in the codon is 
		// (24499-21563)%3 = 2 ie last pos in the codon
		// 
		// Codon ref seq is ref(start:24497, end:24499) = GAT
		//
		// Ref = T, Alt = C -- Alt Codon seq = GAC
		// ref AA is D
		// alt AA is D ==> no change

		rec, err := vcf.ParseVCFRecord(SNP1)
		Check(err)

		res_changes, res_frameshift :=
			AminoAcidChanges(ref, rec, gene_intervals, codon_table)
		corr_changes, corr_frameshift := ".", "false"
		compare(res_changes, corr_changes, t)
		compare(res_frameshift, corr_frameshift, t)
	})
	t.Run("SNP2", func(t *testing.T) {
		// SNP occurs at pos 22317 in gene S.
		// S starts at 21563 (1-based).
		// Codon offset is (22317-21563)/3 = 251.  
		//
		// Relative offset pos in the codon is 
		// (22317-21563)%3 = 1 ie middle pos in codon
		// 
		// Codon ref seq is ref(start:22316, end:22318) = GGT
		//
		// Ref = G, Alt = T -- Alt Codon seq = GTT
		// ref AA is G
		// alt AA is V ==> 1-based codon notation 252G>V

		rec, err := vcf.ParseVCFRecord(SNP2)
		Check(err)

		res_changes, res_frameshift :=
			AminoAcidChanges(ref, rec, gene_intervals, codon_table)
		corr_changes, corr_frameshift := "252G>V", "false"
		compare(res_changes, corr_changes, t)
		compare(res_frameshift, corr_frameshift, t)
	})
	t.Run("DEL", func(t *testing.T) {
		// DEL occurs from 26158 to 26161 in gene ORF3a.
		// ORF3a starts at 25393 (1-based).
		// Codon start is (26158-25393)/3 = 255.  (0-based)
		//   - offset pos within codon (26158-25393) % 3 = 0
		// Codon end is (26161-25393)/3 = 256.  
		//   - offset pos within codon (26161-25393)%3 = 0
		//
		// Codon ref seq is ref(start:26158, end:26163) = GTTAAT
		//
		// Ref = GTTA, Alt = DEL -- Alt Codon seq = ----AT
		// ref AA is VN
		// alt AA is -- ==> 1-based codon notation 256V>del,257N>del
		// Frameshift = true

		rec, err := vcf.ParseVCFRecord(DEL)
		Check(err)

		res_changes, res_frameshift :=
			AminoAcidChanges(ref, rec, gene_intervals, codon_table)
		corr_changes, corr_frameshift := "256V>del,257N>del", "true"
		compare(res_changes, corr_changes, t)
		compare(res_frameshift, corr_frameshift, t)
	})
	t.Run("INS", func(t *testing.T) {
		// INS occurs AFTER 22205 in gene S
		// S starts at 21563 (1-based).
		// Codon start is (22205-21563)/3 = 214.  (0-based)
		//   - offset pos within codon (22205-21563)%3 = 0
		// 
		// Codon ref seq is ref(start:22205, end:22207) = GAT
		//
		// Ref = ins, Alt = CGGCAGGCT; Alt Codon seq = GCGGCAGGCTAT
		// ref AA is D
		// alt AA is AAGY ==> 1-based codon notation 215D>A,215ins>A,215ins>G,215ins>Y
		// Frameshift = False
		rec, err := vcf.ParseVCFRecord(INS)
		Check(err)

		res_changes, res_frameshift :=
			AminoAcidChanges(ref, rec, gene_intervals, codon_table)
		corr_changes, corr_frameshift := "215D>A,215ins>A,215ins>G,215ins>Y", "false"
		compare(res_changes, corr_changes, t)
		compare(res_frameshift, corr_frameshift, t)
	})
	t.Run("COMPOUND_DEL", func(t *testing.T) {
		// Compound variant occurs AFTER 28273 and BEFORE 28283	in gene N
		// N starts at 28274 (1-based) (TSS so start codon shouldn't be affected?)
		// Codon start is (28274-28274)/3 = 0.  (0-based)
		//   - offset pos within codon (28274-28274)%3 = 0
		// Codon end is (28282-28274)/3 = 2.  (0-based)
		//   - offset pos within codon (28283-28274)%3 = 2
		// 
		// Codon ref seq is ref(start:28274, end:28282) = ATGTCTGATAAT
		//
		// Ref = ATGTCTGAT, Alt = -TGTCTCTA; Alt Codon seq = TGTCTCTAAAT
		// ref AA is MSDN
		// alt AA is CL* ==> 1-based codon notation 1M>C,2S>L,3D>*,4N>del
		// Frameshift = True
		rec, err := vcf.ParseVCFRecord(COMPOUND_DEL)
		Check(err)

		res_changes, res_frameshift :=
			AminoAcidChanges(ref, rec, gene_intervals, codon_table)
		corr_changes, corr_frameshift := "1M>C,2S>L,3D>*,4N>del", "true"
		compare(res_changes, corr_changes, t)
		compare(res_frameshift, corr_frameshift, t)
	})
	t.Run("COMPOUND_INS", func(t *testing.T) {
		// Compound variant occurs AFTER 11082 and BEFORE 11084	in gene ORF1a
		// ORF1a starts at 266 (1-based)
		// Codon start is (11082-266)/3 = 3605.  (0-based)
		//   - offset pos within codon (11082-266)%3 = 1
		// Codon end is (11084-266)/3 = 3606.  (0-based)
		//   - offset pos within codon (11084-266)%3 = 0
		// 
		// Codon ref seq is ref(start:11081, end:11086) = TTGTAT
		//                                           siii
		// Ref = G---, Alt = TTTT; Alt Codon seq = TTTTTTTAT
		// ref AA is LY
		// alt AA is FFY ==> 1-based codon notation 3606L>F,3606ins>F
		// Frameshift = false
		rec, err := vcf.ParseVCFRecord(COMPOUND_INS)
		Check(err)

		res_changes, res_frameshift :=
			AminoAcidChanges(ref, rec, gene_intervals, codon_table)
		corr_changes, corr_frameshift := "3606L>F,3606ins>F", "false"
		compare(res_changes, corr_changes, t)
		compare(res_frameshift, corr_frameshift, t)
	})
	t.Run("COMPOUND_SNP", func(t *testing.T) {
		// Compound variant occurs AFTER 28874 and BEFORE 28888	in gene N
		// N starts at 28274 (1-based)
		// Codon start is (28874-28274)/3 = 200.  (0-based)
		//   - offset pos within codon (28874-28274)%3 = 0
		// Codon end is (28888-28274)/3 = 204.  (0-based)
		//   - offset pos within codon (28888-28274)%3 = 2
		// 
		// Codon ref codon seq is ref(start:28874, end:28888) = AGCAGTAGGGGAACT
		// Ref = GCAGTAGGGGAAC, Alt = TCAGTAGGGGAAT; Alt Codon seq = ATCAGTAGGGGAATT
		//         201
		//           v
		// ref AA is SSRGT
		// alt AA is ISRGI ==> 1-based codon notation 201S>I,205T>I
		// Frameshift = false
		rec, err := vcf.ParseVCFRecord(COMPOUND_SNP)
		Check(err)

		res_changes, res_frameshift :=
			AminoAcidChanges(ref, rec, gene_intervals, codon_table)
		corr_changes, corr_frameshift := "201S>I,205T>I", "false"
		compare(res_changes, corr_changes, t)
		compare(res_frameshift, corr_frameshift, t)
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
