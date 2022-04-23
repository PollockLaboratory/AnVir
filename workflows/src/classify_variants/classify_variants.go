package classify_variants

import (
	"bufio"
	// "fmt"
	"path/filepath"
	"os"
	"strconv"
	"strings"
	// TODO use: https://github.com/pkg/errors
	arg "github.com/alexflint/go-arg"

	"annotation/fastaseq"
	"annotation/utils"
	"annotation/vcf"
)

// used when parsing variants
const N_HEADER = 2
const ID = 0
const COUNT = 2
const SEQ = 5

type cliargs struct {
	Reference string `arg:"--reference,required,help:Reference fasta."`
	Variants  string `arg:"--variants,required,help:Variants table."`
	Outfile   string `arg:"--outfile,required,help:Output vcf"`
	K         int    `arg:"-k,required,help:kmer length"`
}
func (c cliargs) Description() string {
	return "Classify variants provided in {variants} with respect to the {reference}.  Output vcf is writen to stdout."
}


type Variant struct {
    // genomic interval of variant
	start int
	end int

	// snp, del, ins, compound
	variant_type string

    // ref seq at the variant genomic position
	ref_allele string

	// alt seq at the genomic position
	// snp: base, del: DEL, ins: inserted sequence,
	// compound: aligned seq to the ref
	alt_allele string    
}

// Take the deviant sequences and merge into a single string.
func MergeDeviants(variant_seq []string, k int) string{
	var sb strings.Builder

	n := len(variant_seq)

	// seed with dev[0], then append the last
	// char from subsequent deviants
	sb.WriteString(variant_seq[1])
	for _, s := range variant_seq[2:n-1] {
		sb.WriteByte(byte(s[k-1]))
	}
	return sb.String()
}

// Given an array containing the anchor sequences (at 0th and last indices)
// and the set of deviant sequences, determine the variant type/genomic position.
// In some cases there could be multiple possible variants return if the anchor
// sequences align to multiple places in the reference in a valid way.
// TODO filter out things that are clearly too long (ie spurious mapping to ref)
func ClassifyVariant(variant_seq []string, k int,
		windowed_ref *fastaseq.WindowedReference,
		contiguous_ref *fastaseq.ContiguousReference) []Variant{
	
	// get anchor sequences
	n := len(variant_seq)
	pre_anchor := variant_seq[0]
	post_anchor := variant_seq[n-1]

	// get all possible "alginments" of the pre/post anchors
	pre_intervals := windowed_ref.Query(pre_anchor)
	post_intervals := windowed_ref.Query(post_anchor)
	interval_pairs := utils.IntervalCartesionProduct(
		pre_intervals, post_intervals)

	// safe to assume at most 3 possibilities since chances
	// of both pre/post intervals aligning to more than 1
	// spot is pretty low.
	variants := make([]Variant, 0, 3)

	// for each pair classify the variant
	for _, anchors := range interval_pairs {
		ref_distance := anchors.Snd.Start - anchors.Fst.End - 1
		n_deviants := n - 2

		// TODO add to end of chain probably spurious mapping of an anchor seq
		// if ref_distance > 100 {
		// 	continue
		if ref_distance == 1 && n_deviants == k {
			/// Simple SNP
			variants = append(variants, Variant{
				start: anchors.Fst.End + 1, // start == end 
				end: anchors.Fst.End + 1,
				variant_type: "SNP",
				ref_allele: contiguous_ref.Query(
					anchors.Fst.End + 1,
					anchors.Fst.End + 1),
				alt_allele: string(variant_seq[1][k-1]),
			})
		} else if merged_deviants := MergeDeviants(variant_seq, k);
		n_deviants == k - 1 && len(merged_deviants) == 2*k-2 {
			// Simple DEL ------------------------------------------------------
			// TODO if I prove that the above 2 properties are equivalient,
			// then I can remove one of those
			variants = append(variants, Variant{
				start: anchors.Fst.End + 1,
				end: anchors.Snd.Start - 1,
				variant_type: "DEL",
				ref_allele: contiguous_ref.Query(
					anchors.Fst.End + 1, anchors.Snd.Start - 1),
				alt_allele: "DEL",
			})
		} else if ref_distance == 0 {
			/// DEL of repeated sequence
			// for example
			//	ref: GGCTGAATAATACGTG
			//  alt: GGCTG---AATACGTG
			//  The location of the DEL is ambiguous so, we report two DELs
			if n_deviants < k {
				// The max overlapping suffix/prefix of the pre/post anchor
				// sequences will give us the repeat sequence that was deleted
				// the deletion either occured in the suffix of the pre anchor
				// or the prefix of the post anchor.
				del_seq := utils.SuffixPrefixOverlap(pre_anchor, post_anchor)
				variants = append(variants, Variant{
					start: anchors.Fst.End - len(del_seq) + 1,
					end: anchors.Fst.End,
					variant_type: "DEL_REPEAT", // should this be its own type?
					ref_allele: del_seq,
					alt_allele: "DEL",
				})
				variants = append(variants, Variant{
					start: anchors.Snd.Start,
					end: anchors.Snd.Start + len(del_seq) - 1,
					variant_type: "DEL_REPEAT", // should this be its own type?
					ref_allele: del_seq,
					alt_allele: "DEL",
				})
			} else {
				/// Simple INS
				// fmt.Println("===================================================================")
				// fmt.Println("DEBUG INFO")
				// fmt.Println("===================================================================")
				// fmt.Printf("ID: %d\n", ID)
				// fmt.Printf("interval_pairs: %+v\n", interval_pairs)
				// fmt.Printf("anchors.Fst: %+v\n", anchors.Fst)
				// i := 0
				// for ;i < n; i++ {
				// 	fmt.Printf("%d:\t%s%s\n", i, strings.Repeat(" ", i), variant_seq[i])
				// }
				// fmt.Printf("anchors.Snd: %+v\n", anchors.Snd)
				// fmt.Printf("merged_deviants: %s\n", merged_deviants)
				// fmt.Printf("len(merged_deviants): %d\n", len(merged_deviants))
				// fmt.Printf("slice range: [%d, %d)\n", k-1, len(merged_deviants)-k+1)
				// fmt.Println("===================================================================")
				variants = append(variants, Variant{
					start: anchors.Fst.End, // 1 before the ins
					end: anchors.Snd.Start, // 1 after the ins
					variant_type: "INS",
					ref_allele: "INS",
					alt_allele: merged_deviants[k-1:len(merged_deviants)-k+1],
				})
			}
		} else if ref_distance < 100 {
			// catch all case for any type of compound variant --------------------
			// take the merged_deviant sequence whose genomic interval will be
			// (anchors.Fst.Start+1) to (anchors.Snd.End-1) (1-based closed)
			// and its associated reference sequence and perform global alignment
			// to catch snps/dels/ins variants that occur in this interval
			ref_seq := contiguous_ref.Query(anchors.Fst.Start + 1, anchors.Snd.End - 1)
			ref_align, alt_align, err := utils.AlignSequences(ref_seq, merged_deviants)
			utils.Check(err)
			variants = append(variants, Variant{
				start: anchors.Fst.End,
				end:   anchors.Snd.Start,
				variant_type: "COMPOUND",
				ref_allele: ref_align[k-1:len(ref_align)-k+1],
				alt_allele: alt_align[k-1:len(alt_align)-k+1],
			})
		}
	}
	return variants
}

func GetVariants(variants_file string, ref_fasta string, k int, out *os.File) {

	// load the reference into windowed and contiguous query structures
	windowed_ref := fastaseq.LoadWindowedReference(ref_fasta, k)
	contiguous_ref := fastaseq.LoadContiguousReference(ref_fasta)

	// Write vcf header to stdout
	vcf.VcfHeader().
		AddReference(ref_fasta).
		AddContig(contiguous_ref.Contig, contiguous_ref.Length()).
		AddInfo("TYPE", "1", "String", "Variant type.").
		AddInfo("END", "1", "Integer", "End position (closed interval)").
		AddInfo("COUNT", "1", "Integer", "Number of occurrences.").
		AddInfo("KMERS", ".", "String", "List of deviant kmer sequences bookended by the prev/next anchor sequences").
		Write(out)

	f, err := os.Open(variants_file)
	utils.Check(err)
	defer f.Close()
	scanner := bufio.NewScanner(f)

	// skip header lines
	for i := 0; i <= N_HEADER; i++ {
		scanner.Scan()
	}
 
	// Parse the variants file, classify, write to vcf (stdout)
	for i:= 0; scanner.Scan(); i++ {
		// get relevant info from the variants table (tab separated)
		fields := strings.Fields(scanner.Text())
		variantID := fields[ID]
		count := fields[COUNT]
		variant_seq := fields[SEQ:]
		// fmt.Printf("%s\n", variantID)

		// get possible variants from this set of deviants
		variants := ClassifyVariant(variant_seq , k, windowed_ref, contiguous_ref)
		for _, v := range variants {
			vcf.VcfRecord().
				SetChrom(contiguous_ref.Contig).
				SetPos(v.start).
				SetID(variantID).
				SetRef(v.ref_allele).
				SetAlt(v.alt_allele).
				SetQual(".").SetFilter(".").
				AddInfo("TYPE", []string{v.variant_type}).
				AddInfo("END", []string{strconv.Itoa(v.end)}).
				AddInfo("COUNT", []string{count}).
				AddInfo("KMERS", variant_seq).
				Write(out)
		}
	}
}

func Main() {
	cli := cliargs{}
	arg.MustParse(&cli)

	outpath, err := filepath.Abs(cli.Outfile)
	utils.Check(err)

	varpath, err := filepath.Abs(cli.Variants)
	utils.Check(err)

	refpath, err := filepath.Abs(cli.Reference)
	utils.Check(err)

	out, err := os.Create(outpath)
	GetVariants(varpath, refpath, cli.K, out)
}
