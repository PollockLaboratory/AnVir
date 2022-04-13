package main
/// TODO list
/*
   TODO load reference as a map keyed by sliding window of sequence
   map[string][int]
 */

import (
	"bufio"
	"fmt"
	"flag"
	"os"
	"strings"
	"strconv"
)

// quick error check pattern
func check(e error) {
	if e != nil {
		panic(e)
	}
}

// simple wrapper for shell wc -l program
func line_count(filename string) int {
	f, err := os.Open(filename)
	check(err)
	defer f.Close()
	scanner := bufio.NewScanner(f)
	var i int
	for i = 0; scanner.Scan(); i++ {}
	return i
}

type Variant struct {
	variantID string     // unique id of variant
	count int            // number of occurrences(?) of variant
	variant_seq []string // prev, deviants, next seqs stored in order
	variant_type string  // snp, del, ins, compound (TODO need to think how to rep a compound var - split?)
	start int            // genomic interval of variant
	end int
}

/*
   /// TODO Parse the variants file into a datastructure
   // useing array of structs
   // to test print out the contents of the array
   // assuming the variants slice has the necessary capacity
*/
func parse_variants(variants_file string, variants []Variant) {
	f, err := os.Open(variants_file)
	check(err)
	defer f.Close()
	scanner := bufio.NewScanner(f)

	// skip first two lines
	for i := 0; i < 2; i++ {
		scanner.Scan()
	}
	// load relevant fields into data structure
	for i:= 0; scanner.Scan(); i++ {
		fields := strings.Fields(scanner.Text())
		variants[i].variantID = fields[0]
		variants[i].count, _ = strconv.Atoi(fields[1])
		variants[i].variant_seq = fields[5:]
	}
}

func main() {
	var filename *string = flag.String("variants", "", "Path to variants file")
	flag.Parse()
	if *filename == "" {
		panic("Filname not provided")
	}
	variants := make([]Variant, line_count(*filename) - 2) // don't count header
	parse_variants(*filename, variants)

	for _, v := range variants {
		fmt.Println(v.variantID)
		fmt.Println(v.count)
		fmt.Println(v.variant_seq)
		fmt.Println("------------------------------------------------------------------------------")
	}
}
