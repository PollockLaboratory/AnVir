package main

import (
	"fmt"
	"os"

	// "github.com/valyala/fasttemplate"

	"annotation/classify_variants"
)

var subprograms = map[string]func(){
	"classify": classify_variants.Main,
	// add more as we get more pieces
}

func printUsage() {
	useage := `
anvir -- usage:
Subprograms:
	classify
	amino -- NOT IMPLEMENTED YET!
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
