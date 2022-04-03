package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"time"
	//"strings"
	//"strconv"

	globals "AnVir/globals"
	seqmer "AnVir/seqmer"
)

// program holds information about this program
type program struct {
	id       int
	name     string
	nickname string
	authors  string
	began    string
	modified string
	uses     string
	runrec   string
	computer string
}

var prog program // need a variable to hold the program structure

// these are the files that hold parameters and control the program settings
const factfile = "factory"
const controlfile = "control"
const modefile = "mode"

// Set adds and prints program information
func (p *program) Set(writer io.Writer) {
	p.name = "tagvars" // read in reference kmer file, count non-reference kmers from seqs file
	// and get unique IDs and sites for variants along with kmer IDs
	fmt.Fprintln(writer, "\n\tRunning Program", p.name)
	p.authors = "David Pollock"
	p.began = "November 2, 2021"
	p.modified = "November 24, 2021"
	fmt.Fprintln(writer, "\tAuthors:", p.authors, "Last Modified:", p.modified, "\n")
}

func main() {
	fmt.Println("Setting Append File")
	fappend, _ := os.OpenFile("access.log", os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	defer fappend.Close()
	writer := bufio.NewWriter(fappend)
	defer writer.Flush() // need this to get output

	// set up the program and globals by reading controls
	prog.Set(writer)
	var globs = globals.New()
	globs.ProgSetUp(controlfile, modefile) // should change through factory and command line
	fmt.Fprintln(writer, "The program was run on", globs.Runstart)
	globs.Print(os.Stdout, "\nStatus after Setup")

	// main program here
	fmt.Println("Starting main program\n")

	refmers := new(seqmer.Oligos) // create global;
	refmers.Init(globs.Geti("klen"), globs.Getf("kcountfile"), globs.Getb("printNs"), globs.Geti("kminprint"))
	refmers.Readk(globs.Getf("kinfile")) // read kmer and counts

	seqs := new(seqmer.Sequences) // create global;
	seqs.Init(globs.Getf("seqfile"), globs.Geti("minseqlen"), globs.Geti("linelimit"), globs.Geti("minline"), globs.Getb("recordseq"), globs.Gets("filetype"), globs.Getb("dofilter"))
	qnkmers := makeKmers(seqs, globs, "qnotk_")

	vars := new(seqmer.Variants) // create global
	vars.Init(globs.Geti("klen"), globs.Getf("varfile"), globs.Geti("kminprint"))
	vars.Read(globs.Getf("varinfile"))
	vars.Print()
	vars.MatchPrint()

	haps := new(seqmer.Haplotypes) // create global
	haps.Init(globs.Geti("klen"), globs.Getf("hapfile"), globs.Geti("kminprint"))
	vars.Addhaps(haps)                        // add haplotype link to vars
	seqs.Tardigrading(qnkmers, refmers, vars) // yet another version
	haps.Print()

	// end main code

	globs.Print(writer, "\nStatus after Program Completion")
	now := time.Now()
	globs.Delta()
	fmt.Fprintln(writer, "\nThe program finished", now)
}

func makeKmers(seqs *seqmer.Sequences, globs *globals.Params, kmerID string) *seqmer.Oligos {
	qbasename := seqs.Name
	qnkmers := new(seqmer.Oligos)
	qnkmers.Init(globs.Geti("klen"), globs.Getf("kcountfile"), globs.Getb("printNs"), globs.Geti("kminprint"))
	qnkmers.Getoutfile(kmerID + globs.Gets("dorevcomp") + "_" + qbasename)
	return qnkmers
}
