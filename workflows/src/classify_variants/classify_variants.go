package classify_variants

import (
	"bufio"
	"fmt"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"

	arg "github.com/alexflint/go-arg"

	"annotation/fastaseq"
	. "annotation/utils"
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
	Threads   int    `arg:"--threads,help:n concurrent threads."`
	K         int    `arg:"--k,required,help:kmer length"`
}
func (c cliargs) Description() string {
	return "Classify variants provided in {variants} with respect to the {reference}."
}


type Variant struct {
	// basic metadata
	id string    
	count string // its read as string so why bother right? :)

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

// Use this function to print some basic debug_info about the variant
// to diagnose issues that come up during testing
func debug_info(id string, anchors Pair[Interval, Interval],
		contiguous_ref *fastaseq.ContiguousReference,
		variant_seq []string, k int) {

	merged_deviants := MergeDeviants(variant_seq, k)
	len_merged := len(merged_deviants)

	// ref seq between the anchors
	ref_seq := contiguous_ref.Query(anchors.Fst.End + 1, anchors.Snd.Start - 1)

	fmt.Println("=============================================================")
	fmt.Println("DEBUG INFO")
	fmt.Println("=============================================================")
	fmt.Printf("ID: %s\n", id)
	fmt.Printf("anchors.Fst: %+v\n", anchors.Fst)
	i := 0
	for ;i < len(variant_seq); i++ {
		fmt.Printf("%d:\t%s%s\n", i, strings.Repeat(" ", i), variant_seq[i])
	}
	fmt.Printf("anchors.Snd: %+v\n", anchors.Snd)
	fmt.Println("-------------------------------------------------------------")
	fmt.Println("len ref_seq: ", len(ref_seq))
	fmt.Println("n:", len_merged)
	fmt.Println("k-1:", k-1)
	fmt.Println("n-k+1:", len_merged-k+1)
	fmt.Println("ref interval: ", anchors.Fst.End + 1, anchors.Snd.Start - 1)
	fmt.Println("ref with anchors: ",
		contiguous_ref.Query(anchors.Fst.Start, anchors.Snd.End))
	fmt.Println("merged deviants:   ", merged_deviants)
	alt_seq := merged_deviants[k-1:len_merged-k+1]
	fmt.Println("alt_seq: ", alt_seq)
	fmt.Println("=============================================================")
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
func ClassifyVariant(id string, count string, variant_seq []string, k int,
		windowed_ref *fastaseq.WindowedReference,
		contiguous_ref *fastaseq.ContiguousReference) []Variant{
	
	// get anchor sequences
	n := len(variant_seq)
	pre_anchor := variant_seq[0]
	post_anchor := variant_seq[n-1]

	// get all possible "alginments" of the pre/post anchors
	pre_intervals := windowed_ref.Query(pre_anchor)
	post_intervals := windowed_ref.Query(post_anchor)
	interval_pairs := IntervalCartesionProduct(
		pre_intervals, post_intervals)

	// safe to assume at most 3 possibilities since chances
	// of both pre/post intervals aligning to more than 1
	// spot is pretty low.
	variants := make([]Variant, 0, 3)

	// for each pair classify the variant
	for _, anchors := range interval_pairs {

		// debug_info(id, anchors, variant_seq, k)

		ref_distance := anchors.Snd.Start - anchors.Fst.End - 1
		n_deviants := n - 2

		// TODO add to end of chain probably spurious mapping of an anchor seq
		// if ref_distance > 100 {
		// 	continue
		if ref_distance == 1 && n_deviants == k {
			/// Simple SNP
			variants = append(variants, Variant{
				id: id, count: count,
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
				id: id, count: count,
				start: anchors.Fst.End + 1,
				end: anchors.Snd.Start - 1,
				variant_type: "DEL",
				ref_allele: contiguous_ref.Query(
					anchors.Fst.End + 1, anchors.Snd.Start - 1),
				alt_allele: "DEL",
			})
		} else if n_deviants < k && ref_distance <= 0 {
			/// DEL of repeated sequence
			// TODO test further
			// for example:
			// ref: GGCTGAATAATACGTG
			// alt: GGCTG---AATACGTG
			// The location of the DEL is ambiguous so, we report two DELs
			// The max overlapping suffix/prefix of the pre/post anchor
			// sequences will give us the repeat sequence that was deleted
			// the deletion either occured in the suffix of the pre anchor
			// or the prefix of the post anchor.
			del_seq := SuffixPrefixOverlap(pre_anchor, post_anchor)
			variants = append(variants, Variant{
				id: id, count: count,
				start: anchors.Fst.End - len(del_seq) + 1,
				end: anchors.Fst.End,
				variant_type: "DEL_REPEAT", // should this be its own type?
				ref_allele: del_seq,
				alt_allele: "DEL",
			})
			variants = append(variants, Variant{
				id: id, count: count,
				start: anchors.Snd.Start,
				end: anchors.Snd.Start + len(del_seq) - 1,
				variant_type: "DEL_REPEAT", // should this be its own type?
				ref_allele: del_seq,
				alt_allele: "DEL",
			})
			
		} else if ref_distance == 0 {
				/// Simple INS
				variants = append(variants, Variant{
					id: id, count: count,
					start: anchors.Fst.End, // 1 before the ins
					end: anchors.Snd.Start, // 1 after the ins
					variant_type: "INS",
					ref_allele: "INS",
					alt_allele: merged_deviants[k-1:len(merged_deviants)-k+1],
				})
		} else if ref_distance < 100 {
			// catch all case for any type of compound variant --------------------
			// TODO add more variety of tests
			len_merged := len(merged_deviants)
			ref_seq := contiguous_ref.Query(anchors.Fst.End + 1, anchors.Snd.Start - 1)
			alt_seq := merged_deviants[k-1:len_merged-k+1]
			ref_align, alt_align, err := AlignSequences(ref_seq, alt_seq, true)
			Check(err)
			variants = append(variants, Variant{
				id: id, count: count,
				start: anchors.Fst.End,
				end:   anchors.Snd.Start,
				variant_type: "COMPOUND",
				ref_allele: ref_align,
				alt_allele: alt_align,
			})
		}
	}
	return variants
}

func GetVariants(variants_file string, ref_fasta string, k int, out *os.File) {

	var wg sync.WaitGroup
	var mu sync.Mutex // for concurrent writes to output

	// load the reference into windowed and contiguous query structures
	windowed_ref := fastaseq.LoadWindowedReference(ref_fasta, k)
	contiguous_ref := fastaseq.LoadContiguousReference(ref_fasta)

	// Write vcf header to stdout
	vcf.VcfHeader().
		AddReference(ref_fasta).
		AddContig(contiguous_ref.Contig, contiguous_ref.Length()).
		AddInfo("VARTYPE", "1", "String", "Variant type.").
		AddInfo("END", "1", "Integer", "End position (closed interval)").
		AddInfo("COUNT", "1", "Integer", "Number of occurrences.").
		AddInfo("KMERS", ".", "String", "List of deviant kmer sequences bookended by the prev/next anchor sequences").
		Write(out)

	f, err := os.Open(variants_file)
	Check(err)
	defer f.Close()
	scanner := bufio.NewScanner(f)

	// skip header lines
	for i := 0; i <= N_HEADER; i++ {
		scanner.Scan()
	}
 
	// Parse the variants file, classify, write to vcf (stdout)
	for scanner.Scan() {
		text := scanner.Text()
		wg.Add(1)
		go func(text string) {
			defer wg.Done()
			// get relevant info from the variants table (tab separated)
			fields := strings.Fields(text)
			variantID := fields[ID]
			count := fields[COUNT]
			variant_seq := fields[SEQ:]
			// fmt.Printf("%s\n", variantID)

			// get possible variants from this set of deviants
			// TODO send the ID/count into the func and add fields to Variant struct
			variants := ClassifyVariant(variantID, count, variant_seq, k,
				windowed_ref, contiguous_ref)
			mu.Lock()
			for _, v := range variants {
				vcf.VcfRecord().
					SetChrom(contiguous_ref.Contig).
					SetPos(v.start).
					SetID(v.id).
					SetRef(v.ref_allele).
					SetAlt(v.alt_allele).
					SetQual(".").SetFilter(".").
					AddInfo("VARTYPE", v.variant_type).
					AddInfo("END", strconv.Itoa(v.end)).
					AddInfo("COUNT", v.count).
					AddInfo("KMERS", variant_seq...).
					Write(out)
			}
			mu.Unlock()
		}(text)
	}
	wg.Wait()
}

func Main() {
	cli := cliargs{Threads: 1}
	arg.MustParse(&cli)

	runtime.GOMAXPROCS(cli.Threads)

	outpath, err := filepath.Abs(cli.Outfile)
	Check(err)

	varpath, err := filepath.Abs(cli.Variants)
	Check(err)

	refpath, err := filepath.Abs(cli.Reference)
	Check(err)

	out, err := os.Create(outpath)
	GetVariants(varpath, refpath, cli.K, out)
}
