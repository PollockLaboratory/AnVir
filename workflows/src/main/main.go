package main

import (
	"fmt"
	"os"

	// "github.com/valyala/fasttemplate"

	"annotation/amino"
	"annotation/classify_variants"
	"annotation/queryposition"
	"annotation/querywindow"
)

var subprograms = map[string]func(){
	"classify": classify_variants.Main,
	"amino": amino.Main,
	"querywindow": querywindow.Main,
	"queryposition": queryposition.Main,
	// add more as we get more pieces
}

func printUsage() {
	useage := `
anvir -- usage:
Subprograms:
	classify:    classify variants from raw deviant/anchor sequences
	amino:       annotate variants in genes with amino acid changes that span the variant
    querywindow: sequence query reference to get genomic position of sequence
    queryposition: given genomic position, get sequence (1-based closed interval)
`
	fmt.Print(useage)
	os.Exit(1)
}

func main() {
	if len(os.Args) < 2 {
		printUsage()
	} else if p, ok := subprograms[os.Args[1]]; !ok {
		printUsage()
	} else {
		// remove the subprogram from args list. ie remove os.Args[1]
		os.Args = append(os.Args[:1], os.Args[2:]...)
		p()
	}
}
