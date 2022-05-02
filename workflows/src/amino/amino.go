// This package will determine amino acid changes for a single variant.
// See tests for conceptual examples.
package amino

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"annotation/fastaseq"
	. "annotation/utils"
	"annotation/vcf"
)

type Change struct {
	From string
	To string

	// codon number of where change took place
	// these use 1-based indexing.  Don't you
	// just love mixing 1 and 0-based indexing?!
	At int
}
func (c Change)str() string {
	return fmt.Sprintf("%d%s>%s", c.At, c.From, c.To)
}

// return comma separated string of formatted change strings
func formatChanges(changes []Change) string{
	change_strings := make([]string, 0, len(changes))
	for i := range changes {
		change_strings = append(change_strings, changes[i].str())
	}
	return strings.Join(change_strings, ",")
}

// Get the 0-based index position wrt the start pos,
// Assuming start and pos have same indexing offset.
func RelativeCoords(start int, pos int) int {
	return pos - start
}

// given 1-based start of gene and genomic pos query,
// get the 0-based codon index.
func CodonPos(start int, pos int) int {
	return RelativeCoords(start, pos) / 3
}

// given gene start & genomic interval (1-based),
// get the range of codon indices (0-based).
func CodonRange(gstart int, gend int, istart int, iend int) (int, int) {
	// the end codon pos will be CodonRange(gstart, gend)
	return Max(CodonPos(gstart, istart), 0),
		Min(CodonPos(gstart, iend), CodonPos(gstart, gend))
}

// Convert codon index to start/end of query.
//// inputs:
// cindex: codon index used as offset from start of gene
// gstart: start coordinate of gene (1-based)
//// outputs:
// int: start coord (1-based) in reference genome
func CodonIndex2Genomic(cindex int, gstart int) int {
	return gstart + cindex*3
}

// given ref genome, gene start (1based), codon index range (0 based... yep)
// get the codons' genomic sequence
func CodonSeq(ref *fastaseq.ContiguousReference,
		gstart int, cstart int, cend int) string {
	// 1
	// |||___|||___|||___
	//  0  1  2  3  4  5
	// say cstart = 2; cend = 4. we need 3 codons
	// so 4-2+1 = 3
	var seq strings.Builder
	// n := cend-cstart+1
	for i := cstart; i <= cend; i++ {
		// convert codon index to start/end of query
		// gstart + 3*cstart , gstart + 
		gpos := CodonIndex2Genomic(i, gstart)
		seq.WriteString(ref.Query(gpos, gpos+2))
	}
	return seq.String()
}


func GetCodonTable(codons_file string) map[string]byte {
	f, err := os.Open(codons_file)
	Check(err)
	defer f.Close()

	// codon -> amino acid
	table := make(map[string]byte)
	scanner := bufio.NewScanner(f)
	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())
		table[fields[0]] = fields[1][0]
	}

	return table
}

// Load BED with format Chrom  Start  End Gene;
// BED intervals are 0-based half open intervals,
// so to convert to 1-based closed to match a VCF,
// we add 1 to the start and leave the end unchanged.
// Returns map[gene] -> Interval
func GetGeneIntervals(genes_bed string) map[string]Interval{
	path, err := filepath.Abs(genes_bed)
	Check(err)

	f, err := os.Open(path)
	Check(err)
	defer f.Close()

	gz, err := gzip.NewReader(f)
	Check(err)

	gene_intervals := make(map[string]Interval, 12)

	scanner := bufio.NewScanner(gz)
	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())
		gene := fields[3]
		start , _ := strconv.Atoi(fields[1])
		end , _ := strconv.Atoi(fields[2])

		// Bed = 0-based, vcf = 1 based.
		// Make the intervals 1 based for consistency
		// between variants and the list of genes.
		gene_intervals[gene] = Interval{Start:start+1, End:end}
	}
	return gene_intervals
}



// apply a SNP to a codon's sequence given the variant's
// relative position within the codon.
func ApplySNP(codon_seq string, pos int, alt_seq byte) string {
	var sb strings.Builder
	for i := range codon_seq {
		if i == pos {
			sb.WriteByte(alt_seq)
		} else {
			sb.WriteByte(codon_seq[i])
		}
	}
	return sb.String()
}

// apply a deletion to a provided codon_seq, given the relative start
// point and the length of the variant.  Using closed interval.
func ApplyDEL(codon_seq string, pos int, length int) string {
	var sb strings.Builder
	for i := 0; i < pos; i++ {
		sb.WriteByte(codon_seq[i])
	}
	for i := pos; i < pos + length; i++ {}
	for i := pos+length; i < len(codon_seq); i++ {
		sb.WriteByte(codon_seq[i])
	}
	return sb.String()
}

// apply an insertion to a provided codon_seq, given the relative start
// point and inserted sequence. Insertion begins AFTER pos.
func ApplyINS(codon_seq string, pos int, inserted_seq string) string {
	var sb strings.Builder
	var i int
	for i = 0; i <= pos; i++ {
		sb.WriteByte(codon_seq[i])
	}
	for j := range inserted_seq {
		sb.WriteByte(inserted_seq[j])
	}
	for ; i< len(codon_seq); i++ {
		sb.WriteByte(codon_seq[i])
	}
	return sb.String()
}

// apply a compound variant to a provided codon_seq, given the relative start
// point and ref/alt allele sequences. Variant begins AFTER pos.
func ApplyCompound(codon_seq string, pos int, ref_seq string, alt_seq string) string {
	var sb strings.Builder
	var i int
	for i = 0; i <= pos; i++ {
		sb.WriteByte(codon_seq[i])
	}
	for j := 0; j < len(ref_seq); j++ {
		if ref_seq[j] == alt_seq[j] { // matches
			sb.WriteByte(ref_seq[j])
			i++
		} else if ref_seq[j] == '-' { // inserted base
			sb.WriteByte(alt_seq[j])
		} else if alt_seq[j] == '-' { // deleted base
			i++
		} else { // polymorphism
			sb.WriteByte(alt_seq[j])
			i++
		}
	}
	for ; i < len(codon_seq); i++ {
		sb.WriteByte(codon_seq[i])
	}
	return sb.String()
}

func DNA2AminoAcid(dna string, codon_table map[string]byte) (string, bool) {
	var sb strings.Builder
	for i := 0; i + 3 <= len(dna); i+=3 {
		if aa, ok := codon_table[dna[i:i+3]]; ok {
			sb.WriteByte(aa)
		} else {
			return "", false
		}
	}
	return sb.String(), true
}

// Get the list of amino acid changes (from, to, at)
// Note: this will only get local changes to the span of ref/alt
// aa sequences. 
/// Inputs:
// ref_aa, alt_aa: ref/variant aa sequences
// cstart: codon start pos (1-based)
/// Returns:
// slice of changes: each change object consists of
// from: orig aa, to: changed aa, at: codon number (1-based)
func GetChanges(ref_aa , alt_aa string, cstart int) []Change {
	n_ref := len(ref_aa)
	n_alt := len(alt_aa)
	changes := make([]Change, 0, Max(n_ref, n_alt))
	// fmt.Printf("ref: %s\nalt: %s\ncstart: %d\n", ref_aa, alt_aa, cstart)
	if n_ref == 0 { // purely novel inserted aa's
		for i := range alt_aa {
			changes = append(changes,
				Change{From: "ins", To: string(alt_aa[i]), At: cstart})
		}
	} else if n_alt == 0 { // all ref aa's deleted
		for i := range ref_aa {
			changes = append(changes,
				Change{From: string(ref_aa[i]), To: "del", At: cstart + i})
		}
	} else if n_ref == 1 && n_alt == 1 {
		if ref_aa != alt_aa {
			changes = append(changes, Change{
				From: ref_aa, To: alt_aa, At: cstart})
		}
	} else {
		// align ref/alt aa
		// TODO test the codon position
		ref_align, alt_align, err := AlignSequences(ref_aa, alt_aa, false)
		Check(err)
		offset := 0 // wrt starting codon index
		for i := range ref_align {
			if ref_align[i] == '-' { // inserted AA
				// codon = before the current offset
				// offset remains the same for next iter
				changes = append(changes,
					Change{From: "ins", To: string(alt_align[i]),
						At: cstart + offset - 1})
			} else if alt_align[i] == '-' { // deleted AA
				changes = append(changes,
					Change{From: string(ref_align[i]),
						To: "del", At: cstart + offset})
				offset++
			} else { 
				if ref_align[i] != alt_align[i] { // substitution
					changes = append(changes,
						Change{From: string(ref_align[i]),
							To: string(alt_align[i]), At: cstart + offset})
				}
				offset++
			}
		}
	}
	return changes
}

// given and vcf record, gene start, get the amino acid changes
func AminoAcidChanges(ref *fastaseq.ContiguousReference, record *vcf.Record,
		gene_intervals map[string]Interval,
		codon_table map[string]byte) (string, string) {

	// get relevant information from vcf line
	pos := record.Pos
	ref_seq := record.Ref
	alt_seq := record.Alt
	end , _:= strconv.Atoi(record.Info["END"][0])

	gene := record.Info["GENE"][0]
	gstart := gene_intervals[gene].Start
	gend := gene_intervals[gene].End

	// Get codon index range
	// TODO right now do limited range of codons.
	// However, keep in mind that a frameshift will
	// expand this range to the end of the gene.
	// Luckily, none of the logic after this ought 
	// to change if I decide to go that route.
	cstart, cend := CodonRange(gstart, gend, pos, end)

	// get refseq spanning the codons
	ref_codons := CodonSeq(ref, gstart, cstart , cend)

	// get position of variant wrt to the first codon
	var_codon_pos := RelativeCoords(gstart, pos) % 3
	
	// applied variants to the codon seq
	// record.Write(os.Stdout)
	var alt_codons string
	// fmt.Println("-------------------------------------------")
	// fmt.Println(record.Info["VARTYPE"])
	// fmt.Println("-------------------------------------------")
	switch record.Info["VARTYPE"][0] {
	case "SNP":
		alt_codons = ApplySNP(ref_codons, var_codon_pos, alt_seq[0])
	case "DEL", "DEL_REPEAT":
		alt_codons = ApplyDEL(ref_codons, var_codon_pos, end-pos+1)
	case "INS":
		alt_codons = ApplyINS(ref_codons, var_codon_pos, alt_seq)
	case "COMPOUND":
		alt_codons = ApplyCompound(ref_codons, var_codon_pos, ref_seq, alt_seq)
	}
	var frameshift string
	if len(alt_codons) % 3 != 0 {
		// TODO someday find the start of shifted reading frame
		frameshift = "true" 
	} else {
		frameshift = "false"
	}

	// get amino acid sequence of ref/alt
	ref_aa , _:= DNA2AminoAcid(ref_codons, codon_table)
	alt_aa , ok := DNA2AminoAcid(alt_codons, codon_table)
	if ok {
		// get list of changes
		changes := GetChanges(ref_aa, alt_aa, cstart)
		return formatChanges(changes), frameshift
	} else {
		return "AMBIGUOUS", "."
	}
}


func AnnotateChanges(ref_fasta string, vcf_file string,
		genes_bed string, codons_file string, outfile string) {

	codon_table := GetCodonTable(codons_file)
	gene_intervals := GetGeneIntervals(genes_bed)

	ref_path, err := filepath.Abs(ref_fasta)
	Check(err)
	ref := fastaseq.LoadContiguousReference(ref_path)

	vcf_path, err := filepath.Abs(vcf_file)
	Check(err)
	v, err := os.Open(vcf_path)
	Check(err)
	defer v.Close()

	out_path, err := filepath.Abs(outfile)
	Check(err)
	out, err := os.Create(out_path)
	Check(err)
	defer out.Close()

	// capture/update header
	vcf.ParseVCFHeader(v).
		AddInfo("AACHANGES", ".", "String",
			"Changes to amino acid sequence within codons spanned by variant.").
		AddInfo("FRAMESHIFT", ".", "String",
			"true/false - does the variant cause a frameshift mutation?").
		Write(out)

	// Annotate the variants with AA changes
	scanner := bufio.NewScanner(v)
	for scanner.Scan() {
		if line := scanner.Text(); line[0] != '#' {
			//fmt.Fprintf(os.Stderr, "Warning", a ...any)
			record, err := vcf.ParseVCFRecord(scanner.Text())
			if err != nil {
				fmt.Fprintf(os.Stderr, "**Warning**:%s\n**SKIPPING**\n\n", err)
				continue
			}
			// record.Write(os.Stdout)
			
			if _, ok := record.Info["GENE"]; !ok {
				// didn't intersect with gene
				record.AddInfo("GENE", ".")
				record.AddInfo("AACHANGES", ".")
				record.AddInfo("FRAMESHIFT", ".")
			} else {
				change_string, frameshift := AminoAcidChanges(
					ref, record, gene_intervals, codon_table)
				record.AddInfo("AACHANGES", change_string).
					AddInfo("FRAMESHIFT", frameshift).
					Write(out)
			}
		}
	}
}
