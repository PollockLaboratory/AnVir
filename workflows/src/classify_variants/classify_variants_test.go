package classify_variants_test

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"reflect"
	"runtime"
	"strings"
	"testing"

	"annotation/classify_variants"
	. "annotation/utils"
)

func compare_strings(correct string, result string, t *testing.T) {
	if result != correct {
		t.Errorf("\ncorrect = %s\nresult = %s\n", correct, result)
	}
}
// ============================================================================
/// Unit testing
// ============================================================================

func TestMergeDeviants(t *testing.T) {
	t.Run("k=3,len(variant_seq)=9", func(t *testing.T) {
		// test sequence: AAABCDEFGGG
		// where AAA and GGG are the "anchor sequences" 
		k := 3
		variant_seq := []string{
			"AAA", "AAB", "ABC",
			"BCD", "CDE", "DEF",
			"EFG", "FGG", "GGG",
		}
		correct := "AABCDEFGG"
		result := classify_variants.MergeDeviants(variant_seq, k)
		compare_strings(correct, result, t)
	})
	t.Run("k=3,len(variant_seq)=5", func(t *testing.T) {
		// test sequence: AAABCCC
		k := 3
		variant_seq := []string{"AAA", "AAB", "ABC", "BCC", "CCC"}
		correct := "AABCC"
		result := classify_variants.MergeDeviants(variant_seq, k)
		compare_strings(correct, result, t)
	})
	t.Run("k=4,len(variant_seq)=4", func(t *testing.T) {
		// test sequence: AAA_CCC (deletion)
		k := 3
		variant_seq := []string{"AAA", "AAC", "ACC", "CCC"}
		correct := "AACC"
		result := classify_variants.MergeDeviants(variant_seq, k)
		compare_strings(correct, result, t)
	})
}

// ============================================================================
/// Full Variant Classification test
// ============================================================================
// TODO what if I have inserted seq that is same as prev/next? it shouldnt occur,
//      because it wouldnt be a deviant sequence (it maps to the ref)
// TODO does the variant caller have the ability to detect tandem duplications?
// TODO do some brainstorming and add a case in the classify variant package if needed


// Test variant classification functionality on a toy dataset.
// See ./test_data/notes.org for information/reason on test cases.
func TestClassifyVariant(t *testing.T) {
	test_fasta, _ := filepath.Abs("test_data/test_ref.fa")
	test_variants, _ := filepath.Abs("test_data/test_variants.tsv")
	path, _ := filepath.Abs("test_data/out.vcf")
	out, err := os.Create(path)
	Check(err)
	classify_variants.GetVariants(test_variants, test_fasta, 5, out)

	/// Simple variants -----------------------------------------
	t.Run("SNP@6-6", func(t *testing.T) {
		out, _ := exec.Command(
			"bcftools", "view", "-i", "ID=\"1\"", "-H", path).CombinedOutput()
		Check(err)
		correct := []string{
			"contig", "6", "1", "T", "t", ".", ".",
			"TYPE=SNP;END=6;COUNT=1;KMERS=ATCGA,TCGAt,CGAtA,GAtAT,AtATG,tATGG,ATGGC",
		}
		result := strings.Fields(string(out))
		if !reflect.DeepEqual(result, correct) {
			t.Errorf("\nCORRECT:\n%s\nRESULT\n%s", correct, result)
		}
	})
	t.Run("DEL:len=1@18-18", func(t *testing.T) {
		out, _ := exec.Command(
			"bcftools", "view", "-i", "ID=\"2\"", "-H", path).CombinedOutput()
		Check(err)
		correct := []string{
			"contig", "18", "2", "T", "DEL", ".", ".",
			"TYPE=DEL;END=18;COUNT=33;KMERS=CGCAT,GCATT,CATTA,ATTAG,TTAGA,TAGAT",
		}
		result := strings.Fields(string(out))
		if !reflect.DeepEqual(result, correct) {
			t.Errorf("\nCORRECT:\n%s\nRESULT\n%s", correct, result)
		}
	})
	t.Run("DEL:len=6@22-27", func(t *testing.T) {
		out, _ := exec.Command(
			"bcftools", "view", "-i", "ID=\"3\"", "-H", path).CombinedOutput()
		Check(err)
		correct := []string{
			"contig", "22", "3", "ATTCGA", "DEL", ".", ".",
			"TYPE=DEL;END=27;COUNT=42;KMERS=TTTAG,TTAGT,TAGTC,AGTCG,GTCGG,TCGGG",
		}
		result := strings.Fields(string(out))
		if !reflect.DeepEqual(result, correct) {
			t.Errorf("\nCORRECT:\n%s\nRESULT\n%s", correct, result)
		}
	})
	t.Run("INS@len=1:12-13", func(t *testing.T) {
		out, _ := exec.Command(
			"bcftools", "view", "-i", "ID=\"4\"", "-H", path).CombinedOutput()
		Check(err)
		correct := []string{
			"contig", "12", "4", "INS", "a", ".", ".",
			"TYPE=INS;END=13;COUNT=11;KMERS=TGGCG,GGCGa,GCGaC,CGaCG,GaCGC,aCGCA,CGCAT",
		}
		result := strings.Fields(string(out))
		if !reflect.DeepEqual(result, correct) {
			t.Errorf("\nCORRECT:\n%s\nRESULT\n%s", correct, result)
		}
	})
	t.Run("INS@len=3:17-18", func(t *testing.T) {
		out, _ := exec.Command(
			"bcftools", "view", "-i", "ID=\"5\"", "-H", path).CombinedOutput()
		Check(err)
		correct := []string{
			"contig", "17", "5", "INS", "abc", ".", ".",
			"TYPE=INS;END=18;COUNT=13;KMERS=CGCAT,GCATa,CATab,ATabc,TabcT,abcTT,bcTTA,cTTAG,TTAGA",
		}
		result := strings.Fields(string(out))
		if !reflect.DeepEqual(result, correct) {
			t.Errorf("\nCORRECT:\n%s\nRESULT\n%s", correct, result)
		}
	})
	/// Compound variants ---------------------------------------
	/// TODO change tests because new aligner is case insensitive
	// NOTE: technically for a really complex variant (may never occur),
	// there could be multiple ways to align. I'm just testing the
	// basic functionality of this component, but it'll give the
	// right alignment most of the time in regions that aren't too complex.
	t.Run("COMPOUND-ADJ-SNPs@:7-10", func(t *testing.T) {
		out, _ := exec.Command(
			"bcftools", "view", "-i", "ID=\"6\"", "-H", path).CombinedOutput()
		Check(err)
		correct := []string{
			"contig", "7", "6", "TG", "tg", ".", ".",
			"TYPE=COMPOUND;END=10;COUNT=55;KMERS=CGATA,GATAt,ATAtg,TAtgG,AtgGC,tgGCG,gGCGC,GCGCG",
		}
		result := strings.Fields(string(out))
		if !reflect.DeepEqual(result, correct) {
			t.Errorf("\nCORRECT:\n%s\nRESULT\n%s", correct, result)
		}
	})
	t.Run("COMPOUND-SNP-DEL@:7-10", func(t *testing.T) {
		out, _ := exec.Command(
			"bcftools", "view", "-i", "ID=\"7\"", "-H", path).CombinedOutput()
		Check(err)
		correct := []string{
			"contig", "7", "7", "TG", "t-", ".", ".",
			"TYPE=COMPOUND;END=10;COUNT=100;KMERS=CGATA,GATAt,ATAtG,TAtGC,AtGCG,tGCGC,GCGCG",
		}
		result := strings.Fields(string(out))
		if !reflect.DeepEqual(result, correct) {
			t.Errorf("\nCORRECT:\n%s\nRESULT\n%s", correct, result)
		}
	})
	t.Run("COMPOUND-SNP-INS@:7-10", func(t *testing.T) {
		out, _ := exec.Command(
			"bcftools", "view", "-i", "ID=\"8\"", "-H", path).CombinedOutput()
		Check(err)
		correct := []string{
			"contig", "7", "8", "T-", "tg", ".", ".",
			"TYPE=COMPOUND;END=9;COUNT=1001;KMERS=CGATA,GATAt,ATAtg,TAtgG,AtgGG,tgGGC,gGGCG,GGCGC",
		}
		result := strings.Fields(string(out))
		if !reflect.DeepEqual(result, correct) {
			t.Errorf("\nCORRECT:\n%s\nRESULT\n%s", correct, result)
		}
	})
	// TODO add test where len(merged_deviants) < k



}

// ============================================================================
/// Benchmark on a large set of variants
// ============================================================================

func BenchmarkVariantClassification(b *testing.B) {
	test_fasta, _ := filepath.Abs("../../data/NC_045512.2.fasta")
	test_variants, _ := filepath.Abs(
		"../../data/variants7M_3500Mline_7.64Mgenomes_min100.tsv")
	path, _ := filepath.Abs("test_data/benchmark.vcf")
	out, err := os.Create(path)
	defer out.Close()
	Check(err)
	runtime.GOMAXPROCS(1)
	fmt.Println("running")
	classify_variants.GetVariants(test_variants, test_fasta, 14, out)
}

