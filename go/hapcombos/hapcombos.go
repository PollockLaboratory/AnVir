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
	p.name = "hapcombos" // read in haplotype file and output the unique subset combinations
	fmt.Fprintln(writer, "\n\tRunning Program", p.name)
	p.authors = "David Pollock"
	p.began = "December 14, 2021"
	p.modified = "December 14, 2021"
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

	haps := new(seqmer.Haplotypes) // create global
	haps.Init(globs.Geti("klen"), globs.Getf("hapfile"), globs.Geti("kminprint"))
	haps.Read(globs.Getf("hapinfile")) // read in from standard file
	//haps.Print(1)                      // print out to check mode 1 is simple mode
	haps.Combos(globs.Geti("hapmincombo"))
	haps.Print(2) // print out to check mode 1 is simple mode
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

//Code 	Meaning
//1 	General error
//2 	Misuse of shell builtins
//127 	Command not found
//128 	Fatal error signal
//130 	Ctrl+C termination

// Exit successfully
// os.Exit(0)// no error

// status code should be in the range [0, 125]
