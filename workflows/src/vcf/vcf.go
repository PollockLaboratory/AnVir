// ============================================================================
// Simple structs for writing variants to vcf format.
// For our purposes we don't need any samples/FORMAT fields,
// so just the basic fields plus INFO
// ============================================================================
package vcf

import (
	"fmt"
	"os"
	"strings"
)


// ============================================================================
/// VCF header
// ============================================================================

type ContigHeader struct {
	ID string
	Length int
}

type InfoHeader struct {
	ID string
	Number string // handle int or char here
	Type string
	Description string
}

type Header struct {
	Contigs []ContigHeader
	Info []InfoHeader
	Reference string // name or path of reference genome
}

func VcfHeader() *Header{
	return &Header{
		Contigs: make([]ContigHeader, 0, 1),
		Info: make([]InfoHeader, 0, 3),
	}
}

func (hd *Header)AddReference(ref string) *Header{
	hd.Reference = ref
	return hd
}

func (hd *Header)AddContig(contig string, length int) *Header {
	hd.Contigs = append(hd.Contigs,
		ContigHeader{ID:contig, Length:length})
	return hd
}

func (hd *Header)AddInfo(id string, number string,
		datatype string, desc string) *Header {
	hd.Info = append(hd.Info, InfoHeader{
		ID: id, Number: number, Type: datatype, Description: desc,
	})
	return hd
}

func (hd *Header)Write(f *os.File) {
	f.WriteString("##fileformat=VCFv4.3\n")
	fmt.Fprintf(f, "##reference=%s\n", hd.Reference)
	for _, i := range hd.Info {
		fmt.Fprintf(f, "##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\n",
			i.ID, i.Number, i.Type, i.Description)
	}
	for _, c := range hd.Contigs {
		fmt.Fprintf(f, "##contig=<ID=%s,length=%d>\n", c.ID, c.Length)
	}
	f.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
}

// ============================================================================
/// VCF record
// ============================================================================

type InfoRecord struct {
	Name string
	Values []string // just keep them as strings (ie no type checking)
}

// return name=value(,) formatted info field
func (i *InfoRecord)str() string {
	return fmt.Sprintf("%s=%s", i.Name, strings.Join(i.Values, ","))
}

// return ';' delimited list of info fields
func format_info(fields []InfoRecord) string {
	s := make([]string, 0, len(fields))
	for _, f := range fields {
		s = append(s, f.str())
	}
	return strings.Join(s, ";")
}

type Record struct {
	Chrom string
	Pos int
	ID string
	Ref string
	Alt string
	Qual string
	Filter string
	Info []InfoRecord
}
func VcfRecord() *Record {
	return &Record{Info: make([]InfoRecord, 0, 3)}
}
func (r *Record)SetChrom(c string) *Record {
	r.Chrom = c
	return r
}
func (r *Record)SetPos(pos int) *Record {
	r.Pos = pos
	return r
}
func (r *Record)SetID(id string) *Record {
	r.ID = id
	return r
}
func (r *Record)SetRef(ref string) *Record {
	r.Ref = ref
	return r
}
func (r *Record)SetAlt(alt string) *Record {
	r.Alt = alt
	return r
}
func (r *Record)SetQual(qual string) *Record {
	r.Qual = qual
	return r
}
func (r *Record)SetFilter(filter string) *Record {
	r.Filter = filter
	return r
}
func (r *Record)AddInfo(name string, values []string) *Record {
	r.Info = append(r.Info, InfoRecord{Name: name, Values: values})
	return r
}
func (r *Record)Write(f *os.File) {
	fmt.Fprintf(f, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t",
		r.Chrom, r.Pos, r.ID, r.Ref,
		r.Alt, r.Qual, r.Filter)
	fmt.Fprintf(f, "%s\n", format_info((r.Info)))
}
