package main

import (
	"bufio"
	"fmt"
	"flag"
	"os"
	"strings"
)

// quick error check pattern
func check(e error) {
	if e != nil {
		panic(e)
	}
}

/*
   /// TODO Parse the variants file into a datastructure
   TODO add struct for variant data
   TODO add return type

   /// Format of variants file:
   * first two rows are header
   * columns (tab separated)
	* 0: variantID
	* 1: count
	* 2: name
	* 3: ID (what's the difference?)
	* 4: devnum
	* 5-end: prev,deviants,next
*/
func parse_variants(variants_file string) {
	// open file
	// skip first two lines
	// for each line
	//    split by tab
	f, err := os.Open(variants_file)
	check(err)
	defer f.Close()
	
	scanner := bufio.NewScanner(f)

	// skip first two lines
	for i := 0; i < 2; i++ {
		scanner.Scan()
	}
	for scanner.Scan() {
		line := scanner.Text()
		fmt.Println(line)
		var x []string = strings.Fields(line)
		for _, s := range x {
			fmt.Println(s)
		}
		break
	}
}


func main() {
	var filename *string = flag.String("variants", "", "Path to variants file")
	flag.Parse()
	if *filename == "" {
		panic("Filname not provided")
	}
	parse_variants(*filename)
}
