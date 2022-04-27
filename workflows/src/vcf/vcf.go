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
	"strconv"
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

// type InfoRecord struct {
// 	Name string
// 	Values []string // just keep them as strings (ie no type checking)
// }



type Record struct {
	Chrom string
	Pos int
	ID string
	Ref string
	Alt string
	Qual string
	Filter string
	infokeys []string // interlal list of keys
	Info map[string][]string
}
func VcfRecord() *Record {
	// return &Record{Info: make([]InfoRecord, 0, 3)}
	return &Record{
		Info: make(map[string][]string, 3),
		infokeys: make([]string, 0, 3),
	}
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

// TODO change to variadic function args instead of list of string
func (r *Record)AddInfo(name string, values ...string) *Record {
	// r.Info = append(r.Info, InfoRecord{Name: name, Values: values})
	r.Info[name] = values
	r.infokeys = append(r.infokeys, name)
	return r
}

// parse an info string and add the components to record
func (r *Record)AddInfoFromString(info string) *Record{
	// split into sep info fields
	fields := strings.Split(info, ";")

	// get key,value pairs and add to the record
	for i := range fields {
		info := strings.Split(fields[i], "=")
		values := strings.Split(info[1], ",")
		r.AddInfo(info[0], values...)
	}
	return r
}

// return ';' delimited list of info fields
// func format_info(fields []InfoRecord) string {
func (r *Record)format_info() string {
	s := make([]string, 0, len(r.Info))
	for _, k := range r.infokeys {
		s = append(s, fmt.Sprintf("%s=%s", k, strings.Join(r.Info[k], ",")))
	}
	return strings.Join(s, ";")
}

func (r *Record)Write(f *os.File) {
	fmt.Fprintf(f, "%s\t%d\t%s\t%s\t%s\t%s\t%s\t",
		r.Chrom, r.Pos, r.ID, r.Ref,
		r.Alt, r.Qual, r.Filter)
	fmt.Fprintf(f, "%s\n", r.format_info())
}

// ============================================================================
/// Functions
// ============================================================================

// read line in vcf and put in data structure
func ParseVCFRecord(line string) *Record {
	fields := strings.Fields(line)
	pos, _ := strconv.Atoi(fields[1])
	return VcfRecord().
		SetChrom(fields[0]).
		SetPos(pos).
		SetID(fields[2]).
		SetRef(fields[3]).
		SetAlt(fields[4]).
		SetQual(".").
		SetFilter(".").
		AddInfoFromString(fields[7])
}
