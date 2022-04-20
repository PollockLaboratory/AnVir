package seqmer

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"

	bitsy "github.com/yourbasic/bit"

	globals "AnVir/globals"
)

//import bitsy "github.com/yourbasic/bit"

//
// ## sequence-related structures and functions
//

var nucs = []byte{'G', 'A', 'C', 'T'}
var cnucs = []byte{'C', 'T', 'G', 'A'}

// sequences holds sequence info
type Sequences struct {
	Name      string
	seqfile   string
	seqmap    map[string]string
	seqfilter map[string]*filter
	Filterdir string
	minlength int
	record    bool
	dofilter  bool
	linelimit int
	linemin   int
	count     int
	isalign   bool
	firstlen  int
	totallen  int
	// could add length map
	outfile    string
	filetype   string
	entrystart string
}

// Init creates new parameter structure of hash types
func (seqs *Sequences) Init(seqfile string, minlen int, llimit int, minline int, dorecord bool, filetype string, dofilter bool) {
	seqs.Name = seqfile
	seqs.seqfile = seqfile
	seqs.seqmap = make(map[string]string)
	seqs.seqfilter = make(map[string]*filter)
	seqs.Filterdir = "none"

	seqs.minlength = minlen
	seqs.linelimit = llimit
	seqs.linemin = minline
	seqs.record = dorecord
	seqs.count = 0
	seqs.isalign = false
	seqs.dofilter = dofilter
	fmt.Println("input and recorded filter states", dofilter, seqs.dofilter)
	seqs.firstlen = 0
	seqs.totallen = 0
	seqs.outfile = "sequences.txt"
	seqs.filetype = filetype
	if filetype == "fastq" {
		seqs.entrystart = "@"
	} else {
		seqs.entrystart = ">"
	}
}

// Print prints out sequences
func (seqs *Sequences) Print(printmode string, headers bool) {
	fmt.Println("Printing sequences to ", seqs.outfile)
	fkout, _ := os.Create(seqs.outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fullmode := "long"
	if headers {
		fmt.Fprintf(kwriter, "%s\t%s\t", "Name", "Length")
		if printmode == fullmode {
			fmt.Fprintf(kwriter, "%s\t", "sequence")
		}
		fmt.Fprintf(kwriter, "\n")
	}
	for name, seq := range seqs.seqmap {
		seqlen := len(seq)
		fmt.Fprintf(kwriter, "%s\t%d", name, seqlen)
		if printmode == fullmode {
			fmt.Fprintf(kwriter, "%s", seq)
		}
		fmt.Fprintf(kwriter, "\n")
	}
}

// Print prints out sequences, parsing genbank name
func (seqs *Sequences) Printparse(printmode string, headers bool) {
	fmt.Println("Printing sequences to ", seqs.outfile)
	fkout, _ := os.Create(seqs.outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	parsename := "shortparsename"
	const firstsplitter = "|"
	const splitter = "="
	if printmode == parsename {
		if headers {
			for name, _ := range seqs.seqmap {
				fmt.Fprintf(kwriter, "%s", name)
				fmt.Fprintf(kwriter, "\n")
			}
			fmt.Fprintf(kwriter, "\n")
		}
		count := 0
		for name, _ := range seqs.seqmap {
			fields := strings.Fields(name)
			count++
			fmt.Fprintf(kwriter, "%d\t", count)
			var namemap map[string]string
			namemap = make(map[string]string)
			lastkey := "empty"
			namemap[lastkey] = lastkey
			// for i, j := 0, 1; i < 10; i, j = i+1, j+1
			bitcount := 0
			for _, bit := range fields {
				bitcount++
				bit := strings.Trim(bit, "[]")
				bitpair := strings.SplitN(bit, splitter, 2)
				if bitcount == 1 {
					bitpair = strings.SplitN(bit, firstsplitter, 2)
				}
				if len(bitpair) > 1 {
					namemap[bitpair[0]] = bitpair[1]
					lastkey = bitpair[0]
				} else {
					bitcount--
					namemap[lastkey] = namemap[lastkey] + " " + bit
				}
			}
			for key, value := range namemap {
				if key != "empty" {
					fmt.Fprintf(kwriter, " %s = %s ", key, value)
				}
			}
			fmt.Fprintf(kwriter, "\t%d\n", len(namemap))
		}
	}
}

// filter holds info about read filters
type filter struct {
	name      string
	readlen   int
	pstart    int
	pend      int
	Wuleft    int
	Wuright   int
	fp        int
	lp        int
	direction string

	olist []string // stable ordered list pointing to kmer info
	klist []*oligo // stable ordered list pointing to kmer info
} //

// reads holds info about the reads (but not the sequences)
type reads struct {
	name   string
	length int
	start  int
	end    int
} //

// primer holds start end and direction of primer
type primer struct {
	pstart int
	pend   int
	dir    string
} //

// recordseq adds seq, name, and count total information to seqs
// this is not used much because files too big
func (seqs *Sequences) recordseq(seq string, name string) {
	seqlen := len(seq)
	if (seqlen > seqs.minlength) && seqs.record {
		blurb := "Error 73, sequence name already exists => "
		if _, exists := seqs.seqmap[name]; exists {
			panic(blurb + name)
		}
		// a fancier version would add a number to the name or something instead of panicking
		seqs.seqmap[name] = seq
		seqs.count += 1
		seqs.totallen += len(seq)
		fmt.Println("Recorded sequence ", name, "length", seqlen)
	}
}

// Printseq prints sequence information in controllable fashion
// not used much because files so big
func (seqs *Sequences) Printseq(printmode string) {
	var name, seq string
	var count, lcount int

	// reading stuff
	fmt.Println("File to open is ", seqs.seqfile)
	fpin, err := os.Open(seqs.seqfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)

	// read, record, count kmers
	noisy := true
	for scanner.Scan() {
		lcount += 1
		line := scanner.Text()              // should not include eol
		trimline := strings.TrimSpace(line) // trim off leading and lagging whitespace
		if strings.HasPrefix(line, ">") {
			fmt.Println("header line", line)
			seqs.recordseq(seq, name)
			seq = ""
			name = strings.TrimPrefix(trimline, ">")
			count += 1
			if noisy {
				fmt.Println("New seq", name, "number", count)
			}
		} else {
			seq = seq + trimline
			keepcompany(lcount, len(seq), 10000, 2000, 100000)
		}
	}
	seqs.recordseq(seq, name)
	fmt.Println("Lines counted1\n", count, lcount)
	fmt.Println("Lines counted\n", count, lcount)
}

//
// //  variant-related structures and functions // //
//

// variant holds info for an individual variant
type variant struct {
	name        string   // might be same as the last reference kmer prior
	ID          int      // increment, unique, fixed if read in directly from variant file
	rank        int      // increment, unique, points to spot in list
	index       int      // index for Variants.varcount and .varlist, run-dependent
	prev        string   // the last reference kmer prior
	next        string   // the next reference kmer posterior
	midmer      string   // the middle deviant kmer
	deviants    []string // the deviant kmers, in order
	lastrefpos  int
	nesxtrefpos int
	start       int      // the start of the variant position
	end         int      // the end of the variant position
	length      int      // the length of the variant (width of non-reference)
	varmer      string   // the composition of the variant
	vartype     string   // eg snp, insert, del, compound
	parent      *variant // if it is really an error on something else
} // will make global seqs

// Init variant creates new if nil, then creates deviants string slice
func (v *variant) Init() {
	if v == nil {
		v = new(variant)
	}
	v.deviants = make([]string, 0)
}

// Clear variant empties existing values for prev, next, midmer, deviants; does nothing if nil
// I think I don't want to use this
func (vinfo *variant) Clear() {
	if vinfo != nil {
		vinfo.prev = ""
		vinfo.next = ""
		vinfo.midmer = ""
		for i := range vinfo.deviants {
			vinfo.deviants[i] = ""
		}
	}
}

// Print an individual variant
// this is hard coded to not print anything with Ns
func (vinfo *variant) Print(varcount int, minprint int, vwriter *bufio.Writer) {
	deviant_count := len(vinfo.deviants)
	if varcount >= minprint {
		hasNs := false
		for dev := range vinfo.deviants {
			if strings.Contains(vinfo.deviants[dev], "N") {
				hasNs = true
			}
		}
		if !hasNs {
			hexstring := fmt.Sprintf("%x", vinfo.ID)
			fmt.Fprintf(vwriter, "%d\t%d\t%d\t%s\t%s\t%d\t", vinfo.ID, vinfo.ID, varcount, vinfo.name, hexstring, deviant_count)
			fmt.Fprintf(vwriter, "%s\t", vinfo.prev)
			for dev := range vinfo.deviants {
				fmt.Fprintf(vwriter, "%s\t", vinfo.deviants[dev])
			}
			fmt.Fprintf(vwriter, "%s\t", vinfo.next)
			fmt.Fprintf(vwriter, "\n")
		}
	}
}

// Sync syncs the tags from reference to focal variant
func (v *variant) Sync(refvar *variant) {
	v.next = refvar.next
	v.prev = refvar.prev
	v.deviants = refvar.deviants
}

// multimaps is a structure that exist so that matches can be ordered in series
// leading to efficient identification of variant matches
type multimap struct {
	key1 map[string]*multimap2
}
type multimap2 struct {
	key2 map[string]*variant
}
type multimap3 struct {
	key3 *variant
}

// Init creates new parameter structure of hash types
func (m *multimap) Init(key2 string, key3 string) {
	m.key1 = make(map[string]*multimap2)
	m.key1[key2] = new(multimap2)
	m.key1[key2].Init(key3)
}

// Init creates new parameter structure of hash types
func (m *multimap2) Init(key3 string) {
	m.key2 = make(map[string]*variant)
	m.key2[key3] = new(variant)
}

// Variants holds map of all variants, tracks total
type Variants struct {
	name       string
	currentvar *variant
	varlist    []*variant
	varcount   []int
	matches    map[string]*multimap
	total      int
	klen       int
	minprint   int
	Infile     string
	Outfile    string
	free       bool
	haps       *Haplotypes
} // will make global variants

// Init creates new parameter structure of hash types
func (vars *Variants) Init(klen int, voutfile string, vminprint int) {
	fmt.Println("Running vinit ", voutfile)
	vars.name = "generic_variant_set"
	vars.currentvar = new(variant)
	vars.currentvar.Init()
	vars.varlist = make([]*variant, 0)
	vars.matches = make(map[string]*multimap)
	vars.varcount = make([]int, 0)
	vars.total = 0
	vars.klen = klen // not strictly necessary, but helpful to have on hand
	vars.minprint = vminprint
	vars.Outfile = voutfile
	vars.free = true
}

// clearCurrent makes the values in currentvar blank or Inits if nil
func (vars *Variants) clearCurrent() {
	vars.currentvar = nil
	vars.currentvar = new(variant)
	vars.currentvar.Init()
}

// Addhaps adds a haplotype to vars
func (vars *Variants) Addhaps(haps *Haplotypes) {
	vars.haps = haps
}

// Print outputs variant info
func (vars *Variants) Print() {
	fmt.Println("Opening Variant  Output File", vars.Outfile)
	fvout, _ := os.Create(vars.Outfile)
	defer fvout.Close()
	vwriter := bufio.NewWriter(fvout)
	defer vwriter.Flush() // need this to get output

	var varcount int
	var vinfo *variant
	minprint := vars.minprint
	fmt.Fprintln(vwriter, "Variant dataset", vars.name, "total variants", vars.total, "min to print", vars.minprint, "k", vars.klen)
	fmt.Fprintln(vwriter, "VariantID\torigID\tcount\tname\thexID\tdevnum\tprev\tdeviants\tnext")
	hexcount := 0
	fmt.Println("length of variant list", len(vars.varlist))
	for i := 0; i < (len(vars.varlist)); i++ {
		varcount = vars.varcount[i]
		vinfo = vars.varlist[i]
		deviant_count := len(vinfo.deviants)

		if varcount >= minprint {
			hasNs := false
			for dev := range vinfo.deviants {
				if strings.Contains(vinfo.deviants[dev], "N") {
					hasNs = true
				}
			}
			if !hasNs {
				hexcount++
				hexstring := fmt.Sprintf("%x", hexcount)
				fmt.Fprintf(vwriter, "%d\t%d\t%d\t%s\t%s\t%d\t", hexcount, vinfo.ID, varcount, vinfo.name, hexstring, deviant_count)
				fmt.Fprintf(vwriter, "%s\t", vinfo.prev)
				for dev := range vinfo.deviants {
					fmt.Fprintf(vwriter, "%s\t", vinfo.deviants[dev])
				}
				fmt.Fprintf(vwriter, "%s\t", vinfo.next)
				fmt.Fprintf(vwriter, "\n")
			}
		}
	}
}

// Print outputs variant info using the variant match structure
func (vars *Variants) MatchPrint() {
	matchoutfile := vars.Outfile + "Match.xls"
	fmt.Println("Opening Variant  Output File", matchoutfile)
	fvout, _ := os.Create(matchoutfile)
	defer fvout.Close()
	vwriter := bufio.NewWriter(fvout)
	defer vwriter.Flush() // need this to get output

	var varcount int
	var vinfo *variant
	fmt.Fprintln(vwriter, "Variant dataset", vars.name, "total variants", vars.total, "min to print", vars.minprint, "k", vars.klen)
	fmt.Fprintln(vwriter, "VariantID\tID2\tcount\tname\tID\tdevnum\tprev\tdeviants\tnext")
	fmt.Println("length of variant list", len(vars.varlist))
	for p := range vars.matches {
		for n := range vars.matches[p].key1 {
			for m := range vars.matches[p].key1[n].key2 {
				vinfo = vars.matches[p].key1[n].key2[m]
				varcount = 100000 // fix this hack 4/6/22
				vinfo.Print(varcount, vars.minprint, vwriter)
			}
		}
	}
}

// Read reads in variant info, will set up variant lists and matches structure for fast lookup
func (vars *Variants) Read(varfile string, mincount int) {
	fmt.Println("File to open is ", varfile)
	fpin, err := os.Open(varfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)
	vars.Infile = varfile //

	const splitter = "\t"
	var linecount int
	var prev, ID, count int
	// var prev, next, ID, count int
	// I am not reading next, for now; I should change it so next is after prev in the file
	for scanner.Scan() {
		line := scanner.Text() // should not include eol
		tokens := strings.Split(line, splitter)
		if linecount < 10 {
			fmt.Println("scanning", varfile)
			fmt.Println("line", linecount, line)
			fmt.Println("length of variant list X", len(vars.varlist))
		}
		if linecount == 0 {
		} else if linecount == 1 { // get header location information
			prev = Index(tokens, "prev")
			ID = Index(tokens, "VariantID")
			count = Index(tokens, "count")
			fmt.Println("location of prev ID count", prev, ID, count)
		} else { // read in variant and add to currentvar
			elements := len(tokens)
			vars.addref(tokens[prev]) // things are being added to currentvar
			thiscount, _ := strconv.Atoi(tokens[count])
			thisID, _ := strconv.Atoi(tokens[ID])
			vars.currentvar.ID = thisID
			enddeviants := elements
			for i := prev + 1; i < elements-1; i++ { // this is a hack in case there are blank elements
				if tokens[i] == "" {
					enddeviants = i
					break
				}
			}
			for i := prev + 1; i < (enddeviants - 1); i++ {
				vars.addnonref(tokens[i]) // add all elements except last as deviants, to currentvar
			}
			pretotal := vars.total
			vars.addref(tokens[enddeviants-1]) // add last element as "next"
			if vars.total > pretotal {
				if linecount < 4300 {
					fmt.Println("info A ", linecount, thiscount, thisID, elements, vars.total)
					fmt.Println("info Y ", vars.varcount[vars.total-1])
					fmt.Println("info Z ", tokens[enddeviants-1], enddeviants)
				}
				vars.varcount[vars.total-1] = thiscount // edited to increase based on thisID, vars.total
				vinfo := vars.varlist[vars.total-1]
				fmt.Println("check vinfo", vinfo.prev, "next", vinfo.next)
			} else {
				// fmt.Println("variant was not added, matched earlier", thisID, vars.total)
			}
		} // end read in variant line
		linecount++
	}
	vars.free = false
	fmt.Println("\nleaving after reading lines, freedom of vars is now", vars.free, "lines", linecount, "vars", vars.total)
}

//
// //  haplotype-related structures and functions // //
//

// haplo holds info for a particular haplotype
// note that we will probably want to incorporate this into bitwise addition
// I downloaded yourbasic/bit which operates quickly on sets of integers so might try it
type haplo struct {
	ID            int   // increment, unique
	variants      []int // the  variants that compose the haplotype, in order
	bitset        *bitsy.Set
	bitstring     string
	hapcount      int               // the number of variants so don't have to recalculate all the time
	children      map[string]*haplo // hash of bitsets[]*haplo // map of one-step children
	descends      map[string]*haplo // map of all descendants via series of children
	status        string            // orig, stub, flub
	lineage_count int               // the number of haplotypes for which this is a subset
	descend_count int               // the total hapcounts in those descendant lineages
} //

// Init creates new parameter structure of hash types
func (h *haplo) Init() {
	h.variants = make([]int, 0)
	h.bitset = bitsy.New()
	h.children = make(map[string]*haplo)
	h.descends = make(map[string]*haplo)
}

// VarsToBits makes bitset to correspond to variant IDs
func (h *haplo) varsToBits() {
	if h.bitset == nil {
		h.bitset = bitsy.New()
	} // should have been initialized already though?
	for vID := range h.variants {
		h.bitset.Add(h.variants[vID])
	}
	h.bitstring = h.bitset.String()
}

// Haplotypes holds map of all haplotypes, as list of haplo
type Haplotypes struct {
	name       string
	currenthap *haplo
	haplist    []*haplo          // currently a list, could be a hash with bitstring key
	hapset     map[string]*haplo // hash of bitsets
	total      int               // total number of haplotypes
	hapcount   []int
	subcount   []int
	klen       int
	minprint   int // we might have a separate minimum haplotype count to print
	Outfile    string
	Infile     string
} // will make global haplos

// Init creates new parameter structure of hash types
func (haps *Haplotypes) Init(klen int, houtfile string, hminprint int) {
	fmt.Println("Running hinit ", houtfile)
	haps.name = "generic_haplotype_set"
	haps.currenthap = new(haplo)
	haps.currenthap.Init()
	haps.haplist = make([]*haplo, 0)
	haps.hapset = make(map[string]*haplo)
	haps.total = 0
	haps.hapcount = make([]int, 0)
	haps.klen = klen // not strictly necessary, but helpful to have on hand
	haps.minprint = hminprint
	haps.Outfile = houtfile
}

// Print outputs haplotype info
func (haps *Haplotypes) Print(mode int) {
	fmt.Println("Opening Haplotype  Output File", haps.Outfile)
	fhout, _ := os.Create(haps.Outfile)
	defer fhout.Close()
	hwriter := bufio.NewWriter(fhout)
	defer hwriter.Flush() // need this to get output

	var hinfo *haplo
	fmt.Println("I am about to print haplotypes stuff.", haps.hapcount, len(haps.haplist))
	fmt.Fprintln(hwriter, "Haplotype dataset", haps.name, "total variants", haps.total, "min to print", haps.minprint, "k", haps.klen)
	fmt.Fprint(hwriter, "ID\tcount\thaplo\t")
	if mode > 1 {
		fmt.Fprint(hwriter, "status\tlineages\tdescends\t")
	}
	fmt.Fprintln(hwriter, "variants\t")
	for i := 0; i < len(haps.haplist); i++ {
		hinfo = haps.haplist[i]
		hcount := hinfo.hapcount
		if hcount >= 0 {
			fmt.Fprintf(hwriter, "%d\t%d\t%s\t", hinfo.ID, hcount, hinfo.bitstring)
			if mode > 1 {
				fmt.Fprintf(hwriter, "%s\t%d\t%d\t", hinfo.status, hinfo.lineage_count, hinfo.descend_count)
			}
			for j := 0; j < len(hinfo.variants); j++ {
				fmt.Fprintf(hwriter, "%d\t", hinfo.variants[j])
			}
			fmt.Fprintf(hwriter, "\n")
		}
	}
}

// EdgePrint outputs ancestral edge info
func (haps *Haplotypes) EdgePrint(edgefile string) {
	fmt.Println("Opening Haplotype  Output File", edgefile)
	fhout, _ := os.Create(edgefile)
	defer fhout.Close()
	hwriter := bufio.NewWriter(fhout)
	defer hwriter.Flush() // need this to get output

	var hinfo *haplo
	fmt.Println("I am about to print edges.", haps.hapcount, len(haps.haplist))
	fmt.Fprintln(hwriter, "Haplotype dataset", haps.name, "total variants", haps.total, "min to print", haps.minprint, "k", haps.klen)
	fmt.Fprintln(hwriter, "ID\tcount\tbits\tstatus\tchildcount\tdescendcount\tID\tchild\tposterior\tID\tancestor\tposterior")
	for i := 0; i < len(haps.haplist); i++ {
		hinfo = haps.haplist[i]
		hcount := hinfo.hapcount
		if hcount >= 0 { // temporary until hcount is sensible
			fmt.Fprintf(hwriter, "%d\t%d\t%s\t%s\t", hinfo.ID, hcount, hinfo.bitstring, hinfo.status)
			childcount := len(hinfo.children)
			descendcount := len(hinfo.descends)
			fmt.Fprintf(hwriter, "%d\t%d\t", childcount, descendcount)
			for j := range hinfo.children {
				posterior := 1.0 / float64(childcount)
				fmt.Fprintf(hwriter, "%d\t%s\t%f\t", hinfo.children[j].ID, hinfo.children[j].bitstring, posterior)
			}
			fmt.Fprintf(hwriter, "X\t")
			for j := range hinfo.descends {
				posterior := 1.0 / float64(descendcount)
				fmt.Fprintf(hwriter, "%d\t%s\t%f\t", hinfo.descends[j].ID, hinfo.descends[j].bitstring, posterior)
			}
			fmt.Fprintf(hwriter, "\n")
		}
	}
}

// Read reads in haplotype info from standard file
// starting with variants.Read()
func (haps *Haplotypes) Read(hapfile string) {
	fmt.Println("File to open is ", hapfile)
	fpin, err := os.Open(hapfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)
	haps.Infile = hapfile //

	const splitter = "\t"
	var linecount int
	var ID, count, hapstr int // header indices
	for scanner.Scan() {
		line := scanner.Text() // should not include eol
		tokens := strings.Split(line, splitter)
		if linecount < 10 {
			fmt.Println("scanning", hapfile)
			fmt.Println("line", linecount, line)
			fmt.Println("length of haplotype list", len(haps.haplist))
		}
		if linecount == 0 {
		} else if linecount == 1 {
			hapstr = Index(tokens, "haplo")
			ID = Index(tokens, "ID")
			count = Index(tokens, "count")
		} else {
			elements := len(tokens)
			thiscount, _ := strconv.Atoi(tokens[count])
			thisID, _ := strconv.Atoi(tokens[ID])
			thishstring := tokens[hapstr]
			hinfo := new(haplo)
			hinfo.Init()
			hinfo.ID = thisID
			hinfo.hapcount = thiscount
			haps.haplist = append(haps.haplist, hinfo)
			haps.total++
			if linecount < 10 {
				fmt.Println("line info", thiscount, thisID, len(haps.haplist), thishstring)
			}

			for i := hapstr + 1; i < elements; i++ {
				thisvar, _ := strconv.Atoi(tokens[i])
				if thisvar != 0 {
					hinfo.variants = append(hinfo.variants, thisvar)
				}
			}
			hinfo.varsToBits()
			hinfo.bitstring = hinfo.bitset.String()
			strmatch(hinfo.bitstring, thishstring, "bitstrings don't match") // hard exit if not
			haps.hapset[thishstring] = hinfo
		} // end read in variant line
		linecount++
	}
	fmt.Println("\nleaving after reading lines", linecount, " total haps", haps.total)
}

// addnew adds new haplotype
func (haps *Haplotypes) addnew(bitset *bitsy.Set, bitstr string, status string) {
	if haps == nil {
		fmt.Println("tried to add to nil haplotypes, this is bad ")
	} else if haps.hapset[bitstr] != nil {
		fmt.Println("tried to overwrite existing haplotype in addnew; this is bad ")
	} else {
		newhap := new(haplo)
		newhap.Init()
		newhap.bitset = bitset // create haplotype bitset
		newhap.bitstring = bitstr
		newhap.status = status
		haps.haplist = append(haps.haplist, newhap) // add currhap to list of haplotypes
		haps.hapset[bitstr] = newhap                // add currhap to bitset hash
		haps.hapset[bitstr].hapcount = haps.hapset[bitstr].hapcount + 1
		haps.currenthap = new(haplo) // this may be a holdover
		haps.currenthap.Init()
		haps.total++
	}
}

func (haps *Haplotypes) addhap(hapbits *bitsy.Set, bitstr string) {
	if haps.hapset[bitstr] != nil {
		oldhap := haps.hapset[bitstr]
		if len(oldhap.variants) > 1 && oldhap.status != "orig" {
			fmt.Println("existing haplotype found", oldhap.ID, oldhap.hapcount, oldhap.bitstring, bitstr, oldhap.status)
		}
	} else {
		haps.addnew(hapbits, bitstr, "stub")
	}
}

// AllCombos goes through all the variants j greater than input ultimate variant
// first merging the previous bitset with jth variant, adding the haplotype created
// then recursively calling itself with the combined bitset and j
func (haps *Haplotypes) AllCombos(hinfo *haplo, varlist []int, parent *haplo, prevbitset *bitsy.Set, ultimate int) {
	varlen := len(hinfo.variants)
	for j := ultimate + 1; j < varlen; j++ {
		//fmt.Println("input is ", prevbitset, ultimate)
		// create combined bitset, make corresponding string, and add haplotype
		combobit := new(bitsy.Set)
		combobit.Set(prevbitset)
		combobit.Add(hinfo.variants[j])
		bitstr := combobit.String()
		haps.addhap(combobit, bitstr)

		// add it as a child to parental hinfo
		child := haps.hapset[bitstr]
		child.lineage_count += 1
		hinfo.descends[bitstr] = child
		parentstring := "nothing"
		newvarlist := make([]int, 0)
		newvarlist = append(varlist, hinfo.variants[j])
		for variant := range newvarlist {
			parentbit := new(bitsy.Set)
			parentbit.Set(combobit)
			parentbit.Delete(newvarlist[variant])
			parentstring = parentbit.String()
			haps.addhap(parentbit, parentstring)
			parent = haps.hapset[parentstring]
			parent.children[bitstr] = child
		}
		haps.AllCombos(hinfo, newvarlist, child, combobit, j)
	}
}

func (hinfo *haplo) testprint(marker string, loc int) {
	fmt.Fprint(os.Stdout, "marker\tstatus\tlineages\tdescends\t")
	fmt.Fprintf(os.Stdout, "%s\t%d\t%s\t%d\t%d\t\n", marker, loc, hinfo.status, hinfo.lineage_count, hinfo.descend_count)
}

// Combos scans the original haplotypes, setting initial combo parameters
// then calls the recursive AllCombos function to find all the variant sub combinations
func (haps *Haplotypes) Combos(hapmincombo int) {
	var hinfo *haplo
	fmt.Println("I am about to analyze combinations.", haps.hapcount, len(haps.haplist))
	for i := 0; i < len(haps.haplist); i++ {
		hinfo = haps.haplist[i]
		hinfo.status = "orig"
		hinfo.lineage_count = 1
		hinfo.descend_count = haps.haplist[i].hapcount
	}
	numhaps := len(haps.haplist) // we add to haplist, but we only want to scroll base haplist
	for i := 0; i < numhaps; i++ {
		hinfo = haps.haplist[i]
		hcount := hinfo.hapcount
		if hcount > hapmincombo {
			combobit := new(bitsy.Set)
			bitstr := combobit.String()
			child := haps.hapset[bitstr]
			varlist := make([]int, 0)
			haps.AllCombos(hinfo, varlist, child, combobit, -1)
		}
	}
	fmt.Println("old and new humnaps", numhaps, len(haps.haplist))
}

//
// //  oligo-related structures and functions // //
//

// oligo holds info for individual kmer
type oligo struct {
	name    string
	revcomp string
	poses   []int
	kcount  int
} // will make global seqs

// Init creates new parameter structure of hash types
func (k *oligo) Init(kmer string, kcount int) {
	k.name = kmer
	k.revcomp = rc(kmer)
	k.kcount = kcount
	k.poses = make([]int, 0)
}

// oligos holds kmer map of all kmers, tracks length and total kmers
type Oligos struct {
	name       string
	kmap       map[string]*oligo
	rmap       map[string]*oligo
	kcount     map[string]int
	locprimer  map[int]int
	readcounts map[int]int
	qmatches   map[string]*QSeqMatches
	primers    map[int]*primer
	total      int
	klen       int
	minprint   int
	Outfile    string
	kfile      string
	printNs    bool
	remnant    string
} // will make global kmers

// Init creates new parameter structure of hash types
func (kmers *Oligos) Init(klen int, koutfile string, doNs bool, kminprint int) {
	fmt.Println("Running kinit ", koutfile)
	kmers.name = "generic_seqfile"
	kmers.kmap = make(map[string]*oligo)
	kmers.rmap = make(map[string]*oligo)
	kmers.kcount = make(map[string]int)
	kmers.locprimer = make(map[int]int)
	kmers.readcounts = make(map[int]int)
	kmers.qmatches = make(map[string]*QSeqMatches)
	kmers.primers = make(map[int]*primer)
	kmers.total = 0
	kmers.minprint = kminprint
	kmers.klen = klen
	kmers.Outfile = koutfile
	kmers.printNs = doNs
}

func (kmers *Oligos) Clearqmatches() {
	kmers.qmatches = nil
	kmers.qmatches = make(map[string]*QSeqMatches)
}

// Holds info for matches to reference sequence
type QSeqMatches struct {
	ID       int
	seqname  string
	seqlen   int
	matches  []*oligo
	matchpos []int
} // will be held in sequence info

// Getoutfile outputs kmer counts and positions
func (kmers *Oligos) Getoutfile(infile string) {
	base := strings.SplitN(infile, ".", 2)
	strlen := strconv.Itoa(kmers.klen)
	kmers.Outfile = kmers.Outfile + strlen + "_" + base[0] + ".kcounts"
}

// GetoutfileMulti outputs kmer counts and positions
func (kmers *Oligos) GetoutfileMulti(kcountpre string, basename string, directory string) {
	strlen := strconv.Itoa(kmers.klen)
	kmers.Outfile = directory + kcountpre + strlen + "_" + basename + ".kcounts"
}

// Kprint outputs kmer counts
func (kmers *Oligos) Kprint() {
	fmt.Println("Opening Kmer Count Output File", kmers.Outfile)
	fkout, _ := os.Create(kmers.Outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fmt.Fprintln(kwriter, "kmer\tcount")
	fmt.Println("size of kmers.kcount", len(kmers.kcount))
	fmt.Println("minprint, printNs", kmers.minprint, kmers.printNs)
	for kmer, kcount := range kmers.kcount {
		if (kcount >= kmers.minprint) && (kmers.printNs || (!strings.Contains(kmer, "N"))) {
			fmt.Fprintf(kwriter, "%s\t%d\n", kmer, kcount)
		}
	}
}

// Kposprint2 outputs kmer positions from reference qmers
func (qmers *Oligos) Kposprint2(refmers *Oligos, kmax int) {
	fmt.Println("Opening Kmer Position Output File", qmers.Outfile)
	fkout, _ := os.Create(qmers.Outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fmt.Fprintln(kwriter, "kmer\tcount\tkmerloc\treadcount\tpositions")
	qmin := qmers.minprint
	for qmer, kcount := range qmers.kcount { // fix this so works with structure
		if (kcount >= qmin) && (qmers.printNs || (!strings.Contains(qmer, "N"))) {
			if refmers.kmap[qmer] != nil {
				poses := refmers.kmap[qmer].poses
				if (poses != nil) && (kcount < kmax) && (kcount >= qmin) {
					kmerloc := poses[0]
					readcount := qmers.readcounts[kmerloc]
					fmt.Fprintf(kwriter, "%s\t%d\t%d\t%d", qmer, kcount, kmerloc, readcount)
					for i := range poses {
						fmt.Fprintf(kwriter, "\t%d", poses[i])
					}
					fmt.Fprintf(kwriter, "\n")
				}
			}
		}
	}
}

// Kposprint outputs kmer positions
func (kmers *Oligos) Kposprint(outfile string, kmin int, kmax int) {
	fmt.Println("Opening Kmer Position Output File", outfile)
	fkout, _ := os.Create(outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fmt.Fprintln(kwriter, "kmer\tcount\tpositions")
	for kmer, kcount := range kmers.kcount { // fix this so works with structure
		if (kcount >= kmin) && (kmers.printNs || (!strings.Contains(kmer, "N"))) {
			if kmers.kmap[kmer] != nil {
				poses := kmers.kmap[kmer].poses
				if (poses != nil) && (kcount < kmax) && (kcount >= kmin) {
					fmt.Fprintf(kwriter, "%s\t%d", kmer, kcount)
					for i := range poses {
						fmt.Fprintf(kwriter, "\t%d", poses[i])
					}
					fmt.Fprintf(kwriter, "\n")
				}
			}
		}
	}
}

// Readk inputs kmer counts
func (kmers *Oligos) Readk(kcountfile string) {
	// reading stuff
	kmers.kfile = kcountfile
	fmt.Println("File to open for Readk() is ", kmers.kfile)
	fpin, err := os.Open(kmers.kfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)

	// read, record, count kmers
	const splitter = "\t"
	const bipart = 2
	var linecount int
	for scanner.Scan() {
		if linecount > 0 {
			line := scanner.Text() // should not include eol
			tokens := strings.Split(line, splitter)
			kmer := tokens[0]
			count, _ := strconv.Atoi(tokens[1])
			kmers.kcount[kmer] = count
			kmers.kmap[kmer] = new(oligo)
			kinfo := kmers.kmap[kmer]
			kinfo.name = kmer
			kinfo.kcount = count
			kinfo.revcomp = rc(kmer)
			kmers.rmap[kinfo.revcomp] = kinfo
			kinfo.poses = make([]int, 0) // imagining option to max pos at 10
		}
		linecount++
	}
	fmt.Println("line count in Readk is ", linecount)
}

// ReadPrimers gets primer locations
func (kmers *Oligos) ReadPrimers(primerfile string, direction string) {
	fmt.Println("File to open for ReadPrimers() is ", primerfile)
	fpin, err := os.Open(primerfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)
	// read, record, primer locations
	const splitter = "\t"
	for scanner.Scan() {
		line := scanner.Text() // should not include eol
		tokens := strings.Split(line, splitter)
		pstart, _ := strconv.Atoi(tokens[0])
		pend, _ := strconv.Atoi(tokens[1])
		// this ought to safety check if already exists
		if direction == "right" {
			temp := pstart
			pstart = pend
			pend = temp
		}
		kmers.primers[pstart] = new(primer)
		primer := kmers.primers[pstart]
		primer.pstart = pstart
		primer.pend = pend
		primer.dir = direction
	}
}

// PrintPrimers prints out list of primers
func (kmers *Oligos) PrintPrimers(outfile string) {
	fmt.Println("Opening primer Output File", outfile)
	fkout, _ := os.Create(outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fmt.Fprintln(kwriter, "firstpos\tendpos\tdirection")
	for _, primer := range kmers.primers {
		fmt.Fprintf(kwriter, "%d\t%d\t%s\n", primer.pstart, primer.pend, primer.dir)
	}
}

// QueryNotRefRC puts qmer counts in qnotkmers if not in kmers
func (qnotkmers *Oligos) QueryNotRefRC(kmers *Oligos, qmers *Oligos, kqmers *Oligos, dorevcomp bool) {
	for qmer, _ := range qmers.kcount {
		origimer := qmer
		orcmer := rc(qmer)           // original but rc
		qmer = strings.ToUpper(qmer) // counting lower case will only happen when it hits
		rcmer := rc(qmer)
		// kinfo :=
		if _, ok := kmers.kcount[qmer]; ok { // match to kmers always UC
			kqmers.kcount[qmer] += qmers.kcount[origimer]
		} else if _, ok := kmers.kcount[rcmer]; ok && dorevcomp {
			kqmers.kcount[qmer] += qmers.kcount[origimer]
		} else {
			if dorevcomp {
				fmt.Println("neither primer nor revcomp found", origimer, orcmer)
			} else {
				// fmt.Println("Original primer not found (no revcomp)", origimer)
			}
		}
	}
	for qmer, _ := range qmers.kcount {
		origimer := qmer
		qmer = strings.ToUpper(qmer) // counting lower case will only happen when it hits
		rcmer := rc(qmer)
		if _, ok := kqmers.kcount[qmer]; !ok {
			if dorevcomp {
				if _, ok := kqmers.kcount[rcmer]; !ok {
					qnotkmers.kcount[qmer] += qmers.kcount[origimer]
				}
			} else {
				qnotkmers.kcount[qmer] += qmers.kcount[origimer]
			}
		}
	}
}

// QueryNotRef puts qmer counts in qnotkmers if not in kmers
func (qnotkmers *Oligos) QueryNotRef(kmers *Oligos, qmers *Oligos, kqmers *Oligos, direction string) {
	kqcount := 0
	qnkcount := 0
	for qmer, _ := range qmers.kcount {
		origimer := qmer
		if direction == "right" {
			qmer = rc(qmer) // original but rc
			fmt.Println("flipping directions", origimer, qmer, direction, kmers.kcount[origimer], kmers.kcount[qmer])
		} else {
			fmt.Println("not flipping", origimer, qmer, direction, kmers.kcount[origimer], kmers.kcount[qmer])
		}
		qmer = strings.ToUpper(qmer)         // counting lower case will only happen when it hits
		if _, ok := kmers.kcount[qmer]; ok { // match to kmers always UC
			kqmers.kcount[origimer] += qmers.kcount[origimer]
			kqcount++
		} else {
			fmt.Println("putting in qnotk", origimer, qmer, direction, kmers.kcount[origimer], kmers.kcount[qmer])
			qnotkmers.kcount[origimer] += qmers.kcount[origimer]
			qnkcount++
		}
	}
	fmt.Println("total each type", kqcount, qnkcount)
}

// QueryCounts adds all kmers with qmer substrings to kqmers
func (kqmers *Oligos) QueryCounts(kmers *Oligos, qmers *Oligos, qkmers *Oligos) {
	for kmer, kcount := range kmers.kcount {
		for i := 0; i < (kmers.klen - qmers.klen + 1); i++ {
			submer := kmer[i : i+qmers.klen]
			if _, ok := qmers.kcount[submer]; ok {
				kqmers.kcount[kmer] = kcount
				qkmers.kcount[submer]++
				if submer == "ACCTTACGAAAC" {
					fmt.Println("testmer", kmer, submer)
				}
			}
		}
	}
}

// kadd adds 1 to the kinfo for given token
// see if this works, but maybe add simple map to kmers, then create kinfo if count is 1
func (kmers *Oligos) kadd(kmer string, kpos int) {
	var kinfo *oligo
	if _, ok := kmers.kmap[kmer]; !ok {
		kmers.kmap[kmer] = new(oligo)
	}
	kinfo = kmers.kmap[kmer] // have to reassign because out of loop
	kinfo.kcount += 1
}

//
// //  intermingled seq and kmer management // //
//

// seqKprint outputs kmer counts
func (seqs *Sequences) SeqKprint(outfile string, kmers *Oligos) {
	fmt.Println("Opening Kmer Count Output File", outfile)
	fkout, _ := os.Create(outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	for seqname, filter := range seqs.seqfilter {
		klist := filter.klist
		//		fmt.Println("checking ", seqname, filter, klist)
		if klist != nil {
			fmt.Fprintln(kwriter, "kmer\tcount\tpos\tseqname")
			for pos, kinfo := range klist {
				poses := kinfo.poses
				if poses != nil {
					fmt.Fprintf(kwriter, "%s\t%d\t%d\t%s\n", kinfo.name, kinfo.kcount, pos, seqname)
				}
			}
		}
	}
}

// SeqPrimerPrint outputs primer ranges, barebones
func (seqs *Sequences) SeqPrimerPrint(appendage string, kmers *Oligos, refmers *Oligos, orient string) {
	outfile := orient + appendage
	fmt.Println("Opening Kmer Count Output File", outfile)
	fkout, _ := os.Create(outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	for _, filter := range seqs.seqfilter {
		minpos := 100000
		maxpos := 0
		klist := filter.klist
		olist := filter.olist
		if (klist != nil) && (orient == olist[0]) {
			for _, kinfo := range klist {
				poses := kinfo.poses
				if poses != nil {
					if kinfo.poses[0] > maxpos {
						maxpos = kinfo.poses[0]
					}
					if kinfo.poses[0] < minpos {
						minpos = kinfo.poses[0]
					}
				}
			}
			fmt.Fprintf(kwriter, "%d\t%d\n", minpos, maxpos+kmers.klen)
		}
	}
}

// QeqQKprint outputs kmer counts
func (seqs *Sequences) SeqQKprint(outfile string, kmers *Oligos, refmers *Oligos) {
	fmt.Println("Opening Kmer Count Output File", outfile)
	fkout, _ := os.Create(outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fmt.Fprintln(kwriter, "kmer\tcount\trefpos\tpos\torient\tseqname")
	for seqname, filter := range seqs.seqfilter {
		klist := filter.klist
		olist := filter.olist
		if klist != nil {
			for pos, kinfo := range klist {
				poses := kinfo.poses
				if poses != nil {
					fmt.Fprintf(kwriter, "%s\t%d\t%d\t%d\t%s\t%s\n", kinfo.name, kinfo.kcount, kinfo.poses[0], pos, olist[0], seqname)
				}
			}
		}
	}
}

//
// functions to manage building up variants //
//

//
// some of these are basically duplicated but done so to avoid messing thing up between
// tagvars and haploscan (which call VarFind and HaploBuilder, respectively)
//

// addref2 adds ref kmer prev or next according to whether prev already exists
// if next then closes current variant
func (vars *Variants) addref2(kmer string) {
	vinfo := vars.currentvar
	if vinfo.prev == "" || len(vinfo.deviants) == 0 {
		vinfo.prev = kmer
	} else {
		vinfo.next = kmer
		vars.closecurrent2()
	}
}

// addnonref adds kmer to list of non reference (deviant) kmers for current variant
func (vars *Variants) addnonref2(kmer string) {
	vinfo := vars.currentvar
	vinfo.deviants = append(vinfo.deviants, kmer)
}

// closecurrent2 closes out the current variant and either saves it to variant list or
// adds it to the previously existing tag path (prev,next, and midmer)
func (vars *Variants) closecurrent2() {
	vinfo := vars.currentvar
	oldinfo := vars.getVarMatch(vinfo) // if new info, returns vinfo or nil if not free to add
	if vars.free {
		if oldinfo == vinfo { // this is a pointer comparison, only true if vinfo is new
			vars.AddNewVariant(vinfo)
		} else { // sets count for this variant at 1; updated later by Read() if needed
			vars.varcount[oldinfo.ID-1]++
		}
	}
	vars.clearCurrent() // set currentvar to nil (hopefully garbage collect) and Init
}

//
// end new copies to manage FindVar to avoid messing up other functions that rely on them
//

// addref adds ref kmer prev or next according to whether prev already exists
// if next then closes current variant
func (vars *Variants) addref(kmer string) {
	vinfo := vars.currentvar
	if vinfo.prev == "" || len(vinfo.deviants) == 0 {
		vinfo.prev = kmer
	} else {
		vinfo.next = kmer
		vars.closecurrent()
	}
}

// addnonref adds kmer to list of non reference (deviant) kmers for current variant
func (vars *Variants) addnonref(kmer string) {
	vinfo := vars.currentvar
	vinfo.deviants = append(vinfo.deviants, kmer)
}

func getMidmer(vinfo *variant) string {
	var midmer string // just choosing a deviant kmer from the middle of the list
	devlength := len(vinfo.deviants)
	if devlength > 0 {
		midpoint := len(vinfo.deviants) / 2
		midmer = vinfo.deviants[midpoint]
	} else {
		midmer = "blank" // this is a bit of a hack; why length 0?
	}
	return midmer
}

// getVarMatch matches input variant info with multimap structure to ID matches
// Behavior:
// if free and no match, create and return a new multimap structure to point to vinfo
// if not free and no match, it will return nil, no new multimap structure
// Otherwise, return full multimap pointer, which will either be old variant or vinfo
// this does not add the variant, nor count it or any somesuch
func (vars *Variants) getVarMatch(vinfo *variant) *variant {
	midmer := getMidmer(vinfo)
	if vars.matches[vinfo.prev] == nil {
		if vars.free {
			vars.matches[vinfo.prev] = new(multimap)
			vars.matches[vinfo.prev].Init(vinfo.next, midmer)
			vars.matches[vinfo.prev].key1[vinfo.next].key2[midmer] = vinfo
			return vars.matches[vinfo.prev].key1[vinfo.next].key2[midmer]
		} else {
			return nil
		}
	} else if vars.matches[vinfo.prev].key1[vinfo.next] == nil {
		if vars.free {
			vars.matches[vinfo.prev].key1[vinfo.next] = new(multimap2)
			vars.matches[vinfo.prev].key1[vinfo.next].Init(midmer)
			vars.matches[vinfo.prev].key1[vinfo.next].key2[midmer] = vinfo // move these to bottom
			return vars.matches[vinfo.prev].key1[vinfo.next].key2[midmer]
		} else {
			return nil
		}
	} else if vars.matches[vinfo.prev].key1[vinfo.next].key2[midmer] == nil {
		if vars.free {
			vars.matches[vinfo.prev].key1[vinfo.next].key2[midmer] = vinfo
			return vars.matches[vinfo.prev].key1[vinfo.next].key2[midmer]
		} else {
			return nil
		}
	} // else newvar == false
	return vars.matches[vinfo.prev].key1[vinfo.next].key2[midmer]
}

func (vars *Variants) Increment(vinfo *variant) {
	//this should check if we care
	ID := vinfo.ID
	vars.varcount[ID-1]++
}

// closecurrent tries to find old info for currentvar
// if !vars.free and oldinfo is found, adds the variant to haplotypes
// 		and clears currentvar
func (vars *Variants) closecurrent() {
	vinfo := vars.currentvar
	oldinfo := vars.getVarMatch(vinfo)    // if new info, returns vinfo or nil if not free to add
	if !vars.free && (oldinfo != vinfo) { // build haplotype if recognized variant
		vars.AddVarToHaps(oldinfo) // only adds if not nil
	}
	vars.clearCurrent() // set currentvar to nil (hopefully garbage collect) and Init
}

// closecurrentdenovo closes out the current variant and either saves it to variant list or
// adds it to the previously existing tag path (prev,next, and midmer)
func (vars *Variants) closecurrentdenovo() {
	vinfo := vars.currentvar
	oldinfo := vars.getVarMatch(vinfo) // if new info, returns vinfo or nil if not free to add
	if vars.free && (oldinfo == vinfo) {
		vars.AddNewVariant(vinfo) // currentvar/vinfo is the new variant; adds ID, etc/ currentvar clear
	} // sets count for this variant at 1; updated later by Read() if needed
	if !vars.free && (oldinfo != vinfo) { // build haplotype if recognized variant
		vars.AddVarToHaps(oldinfo) // only adds if not nil
	}
	vars.clearCurrent() // set currentvar to nil (hopefully garbage collect) and Init
}

// AddNewVariant takes current variant, adds to vars memory, and creates new currentvar
// if not variants not free to change, does nothing but clear out currentvar
func (vars *Variants) AddNewVariant(vinfo *variant) {
	if vars.free {
		vars.varlist = append(vars.varlist, vinfo) // add variant info, not ID, to list of variants
		vars.varcount = append(vars.varcount, 1)
		vars.total = len(vars.varlist) // total number of variants
		vinfo.rank = vars.total        // the ID for this variant ID is the new total
		vinfo.ID = vars.total          // the ID for this variant ID is the new total
	}
	vars.currentvar = new(variant) // clean it out
}

// AddVarToHaps takes input variant, and adds to haplotype list
func (vars *Variants) AddVarToHaps(vinfo *variant) {
	if vars.haps == nil {
		fmt.Println("tried to add to nil haplotypes, this is a problem. ")
	} else if vinfo == nil {
		//fmt.Println("vinfo is nil ")
		// this happens all the time so no warning!
	} else {
		hinfo := vars.haps.currenthap
		hinfo.variants = append(hinfo.variants, vinfo.ID)
	}
}

// close current haplotype
func (vars *Variants) closecurrenthap() {
	haps := vars.haps
	currhap := haps.currenthap
	if haps == nil && !vars.free {
		fmt.Println("tried to add to close haplotypes, this probably shouldn't have happened ")
	} else {
		//numvars := len(currhap.variants)
		currhap.varsToBits() // create haplotype bitset
		//fmt.Println("just made it a bitstring ", currhap.bitstring, numvars)
		currbits := currhap.bitstring
		prevhap := haps.hapset[currbits]
		if prevhap == nil {
			haps.haplist = append(haps.haplist, currhap)     // add currhap to list of haplotypes
			haps.hapset[currbits] = currhap                  // add currhap to bitset hash
			haps.hapset[currbits].ID = len(haps.haplist) - 1 // Added Jan 30 2020 but how was this ever working?
		}
		haps.hapset[currbits].hapcount = haps.hapset[currbits].hapcount + 1
		haps.currenthap = new(haplo)
		haps.currenthap.Init()
		haps.total++
	}
}

// Countnonref counts non reference kmers in string, adds to stored counts in kmers
func (kmers *Oligos) Countnonref(seqs *Sequences, seq string, name string, refmers *Oligos, vars *Variants) {
	var kmer string
	for i := 0; i < (len(seq) - kmers.klen + 1); i++ {
		kmer = seq[i : i+kmers.klen]
		if refmers.kcount[kmer] > 0 { // no filter in place on refmer; these are mostly 1
			vars.addref(kmer)
		} else {
			vars.addnonref(kmer)
			if kmers.kmap[kmer] == nil {
				kmers.kmap[kmer] = new(oligo)
				kmers.kmap[kmer].Init(kmer, kmers.kcount[kmer])
			}
			kmers.kcount[kmer]++
			kmers.total++
		}
	}

	// create remnant in case it needs to be pre-pendend to next sequence fragment
	nextpos := len(seq) - kmers.klen + 1
	kmers.remnant = seq[nextpos : nextpos+kmers.klen-1]
}

// Findnonref finds stretches of non reference kmers in string
// stretches are only recorded if pass criteria
// adds nonref counts to storage in kmers
func (kmers *Oligos) Findnonref(seqs *Sequences, seq string, name string, refmers *Oligos, vars *Variants) {
	var kmer string
	for i := 0; i < (len(seq) - kmers.klen + 1); i++ {
		kmer = seq[i : i+kmers.klen]
		if refmers.kcount[kmer] > 0 { // no filter in place on refmer; these are mostly 1
			vars.addref2(kmer) // either start or end tentative stretch
		} else {
			vars.addnonref2(kmer) // continue adding to tentative stretch
			kmers.record(kmer)
		}
	}
	kmers.addremnant(seq)
}

// addremnant adds remnant to kmers in case it needs to be pre-pendend to next sequence fragment
func (kmers *Oligos) addremnant(seq string) {
	nextpos := len(seq) - kmers.klen + 1
	kmers.remnant = seq[nextpos : nextpos+kmers.klen-1]
}

// record adds kmer to the kmers record
func (kmers *Oligos) record(kmer string) {
	if kmers.kmap[kmer] == nil {
		kmers.kmap[kmer] = new(oligo)
		kmers.kmap[kmer].Init(kmer, kmers.kcount[kmer])
	}
	kmers.kcount[kmer]++
	kmers.total++
}

// Countref counts kmers in string, adds to stored counts in kmers
func (kmers *Oligos) Countref(seqs *Sequences, seq string, name string) {
	var kmer string
	for i := 0; i < (len(seq) - kmers.klen + 1); i++ {
		kmer = seq[i : i+kmers.klen]
		if kmers.kmap[kmer] == nil {
			kmers.kmap[kmer] = new(oligo)
			kmers.kmap[kmer].Init(kmer, kmers.kcount[kmer])
		}
		kmers.kcount[kmer]++
		kmers.total++
	}
	fmt.Println("size of kmers.kcount", len(kmers.kcount), kmers.total)
	// create remnant in case it needs to be pre-pendend to next sequence fragment
	nextpos := len(seq) - kmers.klen + 1
	kmers.remnant = seq[nextpos : nextpos+kmers.klen-1]
}

//
// large database sequence readers and kmerizers //
//

func checkfilter(seqfilter map[string]*filter, name string) bool {
	var passfilter bool
	if seqfilter[name] != nil {
		passfilter = true
	} else {
		passfilter = false
	}
	return passfilter
}

// Kmerize reads fasta or fastq file, turns into kmers and counts them
// filtering is implemented but would require the seqs.seqfilter list to be non-nil, needs re-testing
func (seqs *Sequences) Kmerize(kmers *Oligos) {
	var name, seq string
	var count, lcount int

	// reading stuff and setup
	fmt.Println("File to open is ", seqs.seqfile)
	fpin, err := os.Open(seqs.seqfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)
	entrycount := 0 // track line in each entry
	entrylimit := 0 // default every line in entry is counted
	// hackish fastq, just reading the first line after name line
	if seqs.filetype == "fastq" {
		entrylimit = 1
	}
	fmt.Println("In Kmerizer, dofilter ", seqs.dofilter)
	fmt.Println("File type is ", seqs.filetype, "and entry limit is", entrylimit)
	passfilter := false // flag to see if name is in filter list

	// read, record, count kmers
	for scanner.Scan() {
		lcount += 1
		line := scanner.Text()              // should not include eol
		trimline := strings.TrimSpace(line) // trim off leading and lagging whitespace
		nofilter := !seqs.dofilter
		if strings.HasPrefix(line, seqs.entrystart) {
			name = strings.TrimPrefix(trimline, seqs.entrystart)
			count += 1
			entrycount = 1
			kmers.remnant = ""
			//fmt.Println("New seq", name, "number", count)
			if seqs.dofilter {
				passfilter = checkfilter(seqs.seqfilter, name)
			}
		} else {
			if (nofilter || passfilter) && (lcount < seqs.linelimit) && (entrylimit < 1 || entrycount <= entrylimit) {
				seq = kmers.remnant + trimline
				entrycount++
				if (lcount <= seqs.linelimit) && (lcount > seqs.linemin) {
					kmers.Countref(seqs, seq, name)
				}
				keepcompany(lcount, len(seq), 50000, 10000, 50000)
			}
		}
	}
	fmt.Println("Seqs and Lines counted\n", count, lcount)
}

// VarFind reads fasta or fastq query file and kmerizes
// then compares query kmers to reference kmers (Findnonref)
// stretches of non-reference kmers are recorded as variants
// if they pass whatever limits are in place
func (seqs *Sequences) VarFind(kmers *Oligos, refmers *Oligos, vars *Variants) {
	var name, seq string
	var count, lcount int

	// reading stuff
	fmt.Println("File to open is ", seqs.seqfile)
	fpin, err := os.Open(seqs.seqfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)
	entrycount := 0               // track line in each entry
	entrylimit := 0               // default every line in entry is counted
	if seqs.filetype == "fastq" { // hackish fastq, just read first line after name line
		entrylimit = 1
	}
	fmt.Println("File type is ", seqs.filetype, "and entry limit is", entrylimit)
	fmt.Println("In VarFinder, not filtering based on sequence name")

	// read, record, count kmers
	for scanner.Scan() {
		if lcount < seqs.linelimit {
			lcount += 1
			line := scanner.Text()              // should not include eol
			trimline := strings.TrimSpace(line) // trim off leading and lagging whitespace
			if strings.HasPrefix(line, seqs.entrystart) {
				name = strings.TrimPrefix(trimline, seqs.entrystart)
				count += 1
				entrycount = 1
				kmers.remnant = ""
				vars.closecurrent() // if there was a current variant, close it off
				if (count % 1000) == 0 {
					fmt.Println("Doing seq", name, "number", count)
				}
			} else {
				if entrylimit < 1 || entrycount <= entrylimit {
					seq = kmers.remnant + trimline
					entrycount++
					if lcount > seqs.linemin {
						kmers.Findnonref(seqs, seq, name, refmers, vars)
					}
				}
			}
		} else {
			break
		}
	}
	vars.closecurrent() // otherwise last variant left hanging
	fmt.Println("Seqs and Lines counted\n", count, lcount)
}

// HapBuilder reads fasta or fastq file, finds known non-ref variants
// and adds the to a growing haplotype which is closed at the end of the sequence
func (seqs *Sequences) HapBuilder(kmers *Oligos, refmers *Oligos, vars *Variants) {
	var name, seq string
	var count, lcount int

	// reading stuff
	fmt.Println("File to open is ", seqs.seqfile)
	fpin, err := os.Open(seqs.seqfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)
	entrycount := 0               // track line in each entry
	entrylimit := 0               // default every line in entry is counted
	if seqs.filetype == "fastq" { // hackish fastq, just read first line after name line
		entrylimit = 1
	}
	fmt.Println("File type is ", seqs.filetype, "and entry limit is", entrylimit)
	fmt.Println("In HapBuilder, not filtering ")

	// read, record, count kmers
	for scanner.Scan() {
		if lcount < seqs.linelimit {
			lcount += 1
			line := scanner.Text()              // should not include eol
			trimline := strings.TrimSpace(line) // trim off leading and lagging whitespace
			if strings.HasPrefix(line, seqs.entrystart) {
				name = strings.TrimPrefix(trimline, seqs.entrystart)
				count += 1
				entrycount = 1
				kmers.remnant = ""
				vars.closecurrent()    // if there was a current variant, close it off
				vars.closecurrenthap() // if there was a current haplotype, close it off
				if (count % 1000) == 0 {
					fmt.Println("Doing seq", name, "number", count)
				}
			} else {
				if entrylimit < 1 || entrycount <= entrylimit {
					seq = kmers.remnant + trimline
					entrycount++
					if lcount > seqs.linemin {
						kmers.Countnonref(seqs, seq, name, refmers, vars)
					}
				}
			}
		} else {
			break
		}
	}
	vars.closecurrent()    // otherwise last variant left hanging
	vars.closecurrenthap() // otherwise last haplotype left hanging
	fmt.Println("Seqs and Lines counted\n", count, lcount)
}

//
// consider whether these are relevant to seqmerAV //
//

// addqmatch adds the qmatch information for a sequence to kmers qmatches struct
func (kmers *Oligos) addqmatch(name string, seqnum int, seqlen int) {
	if _, ok := kmers.qmatches[name]; !ok {
		kmers.qmatches[name] = new(QSeqMatches)
		kmers.qmatches[name].matches = make([]*oligo, 0)
		kmers.qmatches[name].matchpos = make([]int, 0)
	}
	qmatch := kmers.qmatches[name]
	qmatch.ID = seqnum
	qmatch.seqname = name
	qmatch.seqlen = seqlen
} // this is still lazily redoing the ID and seqname

// addmatch adds the match information to kmers qmatches struct
// if no hit in kmap, looks for match in rmap
func (kmers *Oligos) addmatch(kmer string, qseqname string, readpos int) int {
	qmatch := kmers.qmatches[qseqname]
	kinfo := kmers.kmap[kmer]
	var nohits int
	if kinfo != nil {
		//fmt.Printf("a")
		qmatch.matches = append(qmatch.matches, kinfo)
		qmatch.matchpos = append(qmatch.matchpos, readpos)
	} else {
		kinfo = kmers.rmap[kmer]
		//fmt.Printf("b")
		if kinfo != nil {
			//fmt.Printf("c")
			qmatch.matches = append(qmatch.matches, kinfo)
			qmatch.matchpos = append(qmatch.matchpos, -readpos)
		}
	}
	if kinfo != nil {
		nohits = 0
	} else {
		nohits = 1
		//fmt.Printf("%s %s ", kmer, rc(kmer))
		if kmers.kmap[kmer] == nil {
			//fmt.Printf("nil ")
		} else {
			//fmt.Printf("w ")
		}
	}
	//mt.Printf("%d ", nohits)
	return nohits
}

// matchmer find kmers in string that match stored kmer file
func (kmers *Oligos) matchmer(seq string, name string, seqnum int) {
	var token string
	seqlen := len(seq)
	kmers.addqmatch(name, seqnum, seqlen)
	//fmt.Println("New sequence ", seqnum)
	nohits := 0
	for i := 0; i < (seqlen - kmers.klen + 1); i++ {
		token = seq[i : i+kmers.klen]
		if kmers.kmap[token] == nil {
			//fmt.Printf("n ")
		} else {
			//fmt.Printf("g ")
		}
		nohits += kmers.addmatch(token, name, i+1)
	}
	//fmt.Println("number no hits in sequence", seqlen-nohits-13, nohits, seqlen, kmers.klen)
	if nohits > 0 {
		//fmt.Println(seq)
		//fmt.Println(rc(seq))
		for i := 0; i < (seqlen - kmers.klen + 1); i++ {
			token = seq[i : i+kmers.klen]
			// fmt.Printf("%s ", token)
			if kmers.kmap[token] == nil {
				//fmt.Printf("nil ")
			} else {
				//fmt.Printf("x ")
			}
		}
		//fmt.Println()
	}
	// kmers.printmatch(name)
	nextpos := len(seq) - kmers.klen + 1
	kmers.remnant = seq[nextpos : nextpos+kmers.klen-1]
}

// Kmatchprint outputs kmer positions
func (kmers *Oligos) Kmatchprint(outfile string) {
	fmt.Println("Opening Kmer Matches Output File", outfile)
	fkout, _ := os.Create(outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	fmt.Fprintf(kwriter, "hits\tmean\tvariance\n")
	fmt.Fprintf(kwriter, "SeqID\tSeqnum\tPosition\tDensity\tKmer\n")
	density := 1.02 // this is fixed for the moment, test phase
	for name := range kmers.qmatches {
		var mean float64
		var total, hits int
		var variance, sumsqr float64

		qmatch := kmers.qmatches[name]
		seqname := qmatch.seqname
		ID := qmatch.ID
		for match := range qmatch.matches { // matches []*oligo
			kinfo := qmatch.matches[match]
			poses2 := kmers.kmap[kinfo.name].poses
			hits += counthits(poses2)
			total += getsum(poses2)
		}
		mean = float64(total) / float64(hits)
		for match := range qmatch.matches { // matches []*oligo
			kinfo := qmatch.matches[match]
			poses2 := kmers.kmap[kinfo.name].poses
			sumsqr += sumdiffsqr(poses2, mean)
		}
		if (mean - 1) > 0.0 {
			variance = sumsqr / (float64(hits) - 1.0)
		} else {
			variance = -1.0
		}
		for match := range qmatch.matches { // matches []*oligo
			kinfo := qmatch.matches[match]
			poses2 := kmers.kmap[kinfo.name].poses
			if poses2 != nil {
				for i := range poses2 {
					fmt.Fprintf(kwriter, "%d\t%f\t%f\t", hits, mean, variance)
					fmt.Fprintf(kwriter, "%s\t%d\t%d\t%f\t%s\n", seqname, ID, poses2[i], density, kinfo.name)
				}
			}
		}
	}
}

//
// Basic Tools
//

//Index is from go by example to find the index positio of a string match in a list
func Index(vs []string, t string) int {
	for i, v := range vs {
		if v == t {
			return i
		}
	}
	return -1
}

// strmatch checks that strings match or exits
func strmatch(query string, match string, words string) {
	if query != match {
		fmt.Println("Exiting string comparison failure upon read", query, match, words)
		os.Exit(123)
	}
}

// rc returns reverse complement
func rc(kmer string) string {
	kbits := []byte(kmer)
	rcbits := []byte(kmer)
	klen := len(kmer)
	for i := range kbits {
		for n := range nucs {
			if nucs[n] == kbits[i] {
				rcbits[klen-i-1] = cnucs[n]
			}
		}
	}
	return string(rcbits)
}

func counthits2(intslice []int) int {
	var hits int
	if intslice != nil {
		for i := range intslice {
			if intslice[i] != 0 {
				hits += 1
			}
		}
	}
	return hits
}

func counthits(intslice []int) int {
	var hits int
	if intslice != nil {
		for i := range intslice {
			if intslice[i] > 0 {
				hits += 1
			}
		}
	}
	return hits
}

func getsum(intslice []int) int {
	var sum int
	if intslice != nil {
		for i := range intslice {
			if intslice[i] > 0 {
				sum += intslice[i]
			}
		}
	}
	return (sum)
}

func sumdiffsqr(intslice []int, mean float64) float64 {
	var sumsqr, realval, diff float64
	if intslice != nil {
		for i := range intslice {
			if intslice[i] > 0 {
				realval = float64(intslice[i])
				diff = realval - mean
				sumsqr += diff * diff
			}
		}
	}
	return (sumsqr)
}

// keepcompany outputs count and seqlen at intervals
// we could control this through globals, but...
func keepcompany(count int, seqlen int, early int, earlyinterval int, interval int) {
	if count < early {
		if (count % earlyinterval) == 0 {
			fmt.Println(count, seqlen)
		}
	}
	if (count % interval) == 0 {
		fmt.Println(count, seqlen)
	}
}

//
// //  BED-related structures and functions, eg City // //
//

// City holds map of all hotels, which contain beds
type City struct {
	numhotels   int
	hotels      map[string]*bedlist
	name        string // should be the name of the sequence file
	outfile     string
	openhotel   string
	chrompos    int
	newbedstart int
	newbedlast  int
	bedcount    int
	minbed      int
	mink        int
	bedlim      int
	maxgap      int
} //

// Init creates new parameter structure of hash types
func (city *City) Init(cityoutfile string, minbedlength int, bedlimit int, maxgap int, minkcount int) {
	fmt.Println("in cinit\n")
	city.numhotels = 0
	city.hotels = make(map[string]*bedlist)
	city.name = "city.txt"
	city.outfile = cityoutfile
	city.openhotel = ""
	city.chrompos = 0
	city.newbedstart = 0
	city.newbedlast = 0
	city.bedcount = 0
	city.mink = minkcount // kmer count to extend bed
	city.minbed = minbedlength
	city.bedlim = bedlimit // keep number of beds under control
	city.maxgap = maxgap
}

// Print creates new parameter structure of hash types
func (city *City) Print() {
	fmt.Println("Printing city", city.name)
	fmt.Println("Opening City Output File", city.outfile)
	fkout, _ := os.Create(city.outfile)
	defer fkout.Close()
	kwriter := bufio.NewWriter(fkout)
	defer kwriter.Flush() // need this to get output

	//fmt.Fprintf(kwriter,"%s\t%s\t%s\t%s\n", "Name", "Start", "Stop", "Length")
	for name, hotel := range city.hotels {
		if true {
			for _, bed := range hotel.beds {
				fmt.Fprintf(kwriter, "%s\t%d\t%d\t%d\n", name, bed.start, bed.stop, bed.len)
			}
		} else {
			// fmt.Fprintf(kwriter,"%s\t%s\thas\t%d\thotels\n", name,hotel.name,hotel.bedcount)
		}
	}
	//fmt.Fprintln(kwriter,"Done writing out city", city.name)
}

// bedlist holds map of beds in hotel
type bedlist struct {
	name     string // should be chromosome
	beds     []*beditem
	bedcount int
} //

type beditem struct {
	start int
	stop  int
	len   int
}

// addbed creates new bed in hotel, gives it start and stop, increments hotel bedcount
func (city *City) addbed(bedlength int) {
	bedcount := city.hotels[city.openhotel].bedcount
	city.hotels[city.openhotel].beds = append(city.hotels[city.openhotel].beds, new(beditem)) // create new bed

	city.hotels[city.openhotel].beds[bedcount].start = city.newbedstart
	city.hotels[city.openhotel].beds[bedcount].stop = city.newbedlast
	city.hotels[city.openhotel].beds[bedcount].len = bedlength

	hotel := city.hotels[city.openhotel]
	hotel.bedcount += 1
	city.bedcount += 1
	if hotel.bedcount < 10 {
		fmt.Println("bedcounts", hotel.bedcount, city.bedcount, hotel.name)
	}
}

// Starthotel creates new hotel in city with name (e.g. chromosome name), opens it
func (city *City) Starthotel(name string) {
	//wtf
	if _, ok := city.hotels[name]; !ok {
		city.numhotels += 1
		city.hotels[name] = new(bedlist)
		city.hotels[name].name = name
		city.openhotel = name
		city.chrompos = 0
		city.newbedstart = 0
		city.newbedlast = 0
	}
}

// Closehotels closes hotels in city
func (city *City) Closehotels() {
	city.openhotel = ""
	city.chrompos = 0
	city.newbedstart = 0
	city.newbedlast = 0
}

// closebed saves the bed information to the current hotel
// and closes the city info on this bed if long enough
func (city *City) closebed() {
	bedlength := city.newbedlast - city.newbedstart
	//fmt.Println("Bed Close bl,cp",bedlength,city.chrompos,city.newbedlast,city.newbedstart)
	if bedlength > city.minbed {
		fmt.Println("Adding bed", bedlength, city.chrompos, city.newbedlast, city.newbedstart)
		city.addbed(bedlength)
	}
	city.newbedlast = 0
}

// openbed opens a new bed in the city; an old bed should not be live
func (city *City) openbed(kstart int, kend int) {
	//fmt.Println("a2bed open start cp",kstart,city.chrompos,kend )
	if 0 == kstart {
		city.newbedstart = city.chrompos + kstart + 1 // adding an offset here in bed files; check this
		city.newbedlast = city.chrompos + kend
	}
}

// Read reads in a city of hotels with beds
func (city *City) Read(bedfile string) {
	// reading stuff
	fmt.Println("File to open is ", bedfile)
	fpin, err := os.Open(bedfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)
	city.name = bedfile // a little inconsistent; when read in it is the seqfile

	// read, record, count kmers
	const splitter = "\t"
	const bipart = 2
	var linecount int
	for scanner.Scan() {
		fmt.Println("scanning", bedfile)
		line := scanner.Text() // should not include eol
		tokens := strings.Split(line, splitter)
		name := tokens[0]
		city.Starthotel(name)                         // should check existence and do nothing if already exists
		city.newbedstart, _ = strconv.Atoi(tokens[1]) // should use error to check
		city.newbedlast, _ = strconv.Atoi(tokens[2])
		bedlength, _ := strconv.Atoi(tokens[3])
		city.addbed(bedlength)
		linecount++
	}
	fmt.Println("\nleaving after reading lines", linecount)
}

// addtobeds updates an ongoing repeat stretch or creates new one if needed
func (city *City) addtobeds(kstart int, klen int) {
	if city.newbedlast > 0 { // a bed is open
		newgap := city.chrompos + kstart - city.newbedlast
		if newgap <= city.maxgap {
			city.newbedlast = city.chrompos + kstart + klen
		} else {
			city.closebed()
			city.openbed(kstart, kstart+klen)
		}
	} else {
		city.openbed(kstart, kstart+klen)
	}
}

// Maidservice check for stretches of kmers in seq, tracks beds in hotel
func (kmers *Oligos) Maidservice(seq string, city *City) {
	var token string
	for i := 0; i < (len(seq) - kmers.klen + 1); i++ {
		token = seq[i : i+kmers.klen]
		if kmers.kcount[token] >= city.mink {
			city.addtobeds(i, kmers.klen)
		}
	}
	nextpos := len(seq) - kmers.klen + 1
	kmers.remnant = seq[nextpos : nextpos+kmers.klen-1]
}

// MakeBeds reads fasta file and tracks kmers, adding to beds in chromosomes (hotel)
//  from a genome (city)
func (seqs *Sequences) MakeBeds(kmers *Oligos, city *City) {
	var name, seq string
	var count, lcount int

	// reading stuff
	fmt.Println("File to open is ", seqs.seqfile)
	fpin, err := os.Open(seqs.seqfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)
	city.name = seqs.seqfile

	// read, record, count kmers
	for scanner.Scan() {
		lcount += 1
		line := scanner.Text()              // should not include eol
		trimline := strings.TrimSpace(line) // trim off leading and lagging whitespace
		if strings.HasPrefix(line, ">") {
			name = strings.TrimPrefix(trimline, ">")
			count += 1
			kmers.remnant = ""
			fmt.Println("New seq", name, "number", count)
			city.Starthotel(name)
		} else {
			seq = kmers.remnant + trimline
			if (lcount <= seqs.linelimit) && (lcount > seqs.linemin) && (city.bedcount < city.bedlim) {
				kmers.Maidservice(seq, city)
			}
			city.chrompos += len(seq) - len(kmers.remnant)
		}
	}
	city.Closehotels()
	fmt.Println("Lines counted\n", count, lcount)
}

// Travelers counts kmers in seq if in bed at right hotel
// parses city, so not order dependent
func (kmers *Oligos) travelers(seq string, city *City, name string) {
	//var token string
	//hotel := city.[name] // find the hotel
	// for each bed in hotel, check if starts within the sequence
	// if so, get appropriate sub-sequence and send to countmers
	//kmers.Countmers(seq, seqs.record)

	nextpos := len(seq) - kmers.klen + 1
	kmers.remnant = seq[nextpos : nextpos+kmers.klen-1]

	// we must have updated chrompos and remnant
	// chrompos should probably be a property of seqs, but hey; think about it
}

// Pantopod reads fasta file annotated according to a city of hotels
// and records kmers, so sort of mode of kcount
func (seqs *Sequences) Pantopod(kmers *Oligos, city *City) {
	var name, seq string
	var count, lcount int

	// reading stuff
	fmt.Println("File to open is ", seqs.seqfile)
	fpin, err := os.Open(seqs.seqfile)
	globals.Check(err)
	defer fpin.Close()
	scanner := bufio.NewScanner(fpin)

	// read, record, count kmers
	for scanner.Scan() {
		lcount += 1
		line := scanner.Text()              // should not include eol
		trimline := strings.TrimSpace(line) // trim off leading and lagging whitespace
		if strings.HasPrefix(line, ">") {
			fmt.Println("header line", line)
			name = strings.TrimPrefix(trimline, ">")
			count += 1
			kmers.remnant = ""
			fmt.Println("New seq", name, "number", count)
		} else {
			seq = kmers.remnant + trimline
			if (count <= seqs.linelimit) && (lcount > seqs.linemin) {
				kmers.travelers(seq, city, name)
				city.chrompos += len(seq) - len(kmers.remnant)
			}
			if lcount < 100000 {
				if (lcount % 10000) == 0 {
					fmt.Println(lcount, len(seq))
				}
			}
			if (lcount % 100000) == 0 {
				fmt.Println(lcount, len(seq))
			}
		}
	}
	fmt.Println("Lines counted\n", count, lcount)
}
