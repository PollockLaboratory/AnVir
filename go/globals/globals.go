// Package globals contains setup functions for factory, globals, program
package globals

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
	"time"
)

// factory may have other than general parameter file
type factory struct {
	invariant string
	variant   string
	modefile  string
	control   string
}

// generic parameter file, maps of various types
type Params struct {
	Runstart time.Time
	Runend   time.Time
	ints     map[string]int
	bools    map[string]bool
	strs     map[string]string
	files    map[string]string
	nums     map[string]float64
	raw      map[string]string
	fixed    map[string]bool
}

// ProgSetUp reads various control files (factory, mode, control) to set up parameters
func (p *Params) ProgSetUp(controlfile string, modefile string) {
	print := fmt.Println // I just wanted to remember how to do this
	now := time.Now()    // this time stuff should be somewhere else?
	p.strs["progstart"] = now.Format("Mon Jan _2 15:04:05 2006")
	p.Runstart = now
	print("program started", p.strs["progstart"])
	// p.readFact(factfile)  // this is not set up yet
	p.readParams(modefile)
	p.readParams(controlfile) // the main difference is the order. Control file modifies mode settings
	// in future, read factory, then mode, then controlfile
	// * means no more changes by later files
}

// Delta reads various control files (factory, mode, control) to set up parameters
func (p *Params) Delta() string {
	print := fmt.Println // I just wanted to remember how to do this
	now := time.Now()    // this time stuff should be somewhere else?
	p.strs["progend"] = now.Format("Mon Jan _2 15:04:05 2006")
	p.Runend = now
	//delta := now - p.Runstart

	t1 := p.Runstart
	// t2 := t1.Add(time.Second * 341)

	fmt.Println(t1)
	fmt.Println(now)

	elapsed := now.Sub(t1)
	fmt.Println("Elapsed time is", int64(elapsed/(time.Second)), "seconds")
	fmt.Println("Elapsed time 2 is", int64(elapsed/(time.Millisecond*1000)), "seconds")
	fmt.Println("Elapsed time is", int64(elapsed/(time.Millisecond)), "milliseconds")
	fmt.Println("Elapsed time is", int64(elapsed/(time.Millisecond/1000)), "microseconds")
	fmt.Printf("difference = %v\n", elapsed)

	print("program ended", p.strs["progend"])
	//print("time differential",delta)

	return "done"
}

const hash = "#"
const star = "*"
const equals = "="
const bipart = 2

// readParams reads user control file
func (p *Params) readParams(controlfile string) {
	fmt.Println("\nReading control file **", controlfile, "**\n")
	fpin, err := os.Open(controlfile)
	Check(err)
	defer fpin.Close()

	scanner := bufio.NewScanner(fpin)
	for scanner.Scan() {
		line := scanner.Text()
		if len(line) > 0 {
			var name, value, comment, invariant = trisplit(line, equals, hash, star)
			if goodentry(name, value) {
				if !(true == p.fixed[name]) {
					p.raw[name] = value
					if invariant {
						p.fixed[name] = true
					}
				} else { // tried to change invariant
					fmt.Println("\nWarning: tried to change invariant parameter", name, "from", p.raw[name], "to", value, ". \n\tThe invariant value", p.raw[name], "was retained.\n")
				}
			}
			fmt.Println(name, equals, value, comment)
		}
	}
	fmt.Println("\nDone reading control file **", controlfile, "**\n")
}

// readFact reads factory setting file *not yet written properly*
func (p *Params) readFact(factoryfile string) {
	fmt.Println("\nReading factory file **", factoryfile, "**\n")
	fpin, err := os.Open(factoryfile)
	Check(err)
	defer fpin.Close()

	scanner := bufio.NewScanner(fpin)
	for scanner.Scan() {
		line := scanner.Text()
		if len(line) > 0 {
			// in previous version, included file type and bounds
			// as well as invariance control
			// would need to write that
			// on idea would be to have a special optional line with file type and bounds
			// would also then need a bounds checker routine
		}
	}
	fmt.Println("\nDone reading factory file **", factoryfile, "**\n")
}

// goodentry checks if key and value both exist (might change to != nil)
func goodentry(key, valstring string) bool {
	if len(key) > 0 && len(valstring) > 0 {
		return true
	} else {
		return false
	}
}

// trisplit takes line and splits in 3 part based on input separators (second has alternate)
func trisplit(line string, sep1 string, sep2 string, altsep2 string) (string, string, string, bool) {
	var invariant bool
	var name, value, comment string

	tokens := strings.SplitN(line, sep2, bipart) //  split on second separator
	if len(tokens) > 1 {
		comment = tokens[1] // if split, part 2 is a comment
	} else {
		tokens = strings.SplitN(line, altsep2, bipart)
		if len(tokens) > 1 {
			invariant = true
		}
	} // if not, try alt separator 2
	if len(tokens) > 1 {
		comment = tokens[1]
	} // if split, param cannot be further modified
	units := strings.SplitN(tokens[0], sep1, bipart) // split into key (name) value pair; still strings
	name = strings.TrimSpace(units[0])
	// fmt.Println("test trisplit",units,"2",name,"3",tokens)
	if (len(units) > 0) && (len(name) > 0) {
		value = strings.TrimSpace(units[1])
	}

	return name, value, comment, invariant
}

// New creates new parameter structure of hash types
func New() *Params {
	var i = make(map[string]int)
	var b = make(map[string]bool)
	var s = make(map[string]string)
	var f = make(map[string]string)
	var n = make(map[string]float64)
	var r = make(map[string]string)
	var c = make(map[string]bool) // boolean fixed
	return &Params{ints: i, bools: b, strs: s, files: f, nums: n, raw: r, fixed: c}
}

// Print parameters by key type in key = value format
func (p *Params) Print(writer io.Writer, blurb string) {
	if blurb != "" {
		fmt.Fprintln(writer, blurb)
	}
	fmt.Fprintln(writer, "\n\tinteger parameters")
	for key, value := range p.ints {
		fmt.Fprintln(writer, key, "=", value)
	}
	fmt.Fprintln(writer, "\n\boolean parameters")
	for key, value := range p.bools {
		fmt.Fprintln(writer, key, "=", value)
	}
	fmt.Fprintln(writer, "\n\tstring parameters")
	for key, value := range p.strs {
		fmt.Fprintln(writer, key, "=", value)
	}
	fmt.Fprintln(writer, "\n\tfile parameters")
	for key, value := range p.files {
		fmt.Fprintln(writer, key, "=", value)
	}
	fmt.Fprintln(writer, "\n\tnumber parameters")
	for key, value := range p.nums {
		fmt.Fprintln(writer, key, "=", value)
	}
	fmt.Fprintln(writer, "\n\traw parameters")
	for key, value := range p.raw {
		fmt.Fprintln(writer, key, "=", value)
	}
	fmt.Fprintln(writer, "\n\tfixed (invariant, constant) parameters")
	for key, value := range p.fixed {
		if value {
			fmt.Fprintln(writer, key, "is invariant")
		}
	}
}

// Geti and relatives (Getn, Gets, Getf, Getr) put raw params into right map
// This could probably be written more efficiently a 1 prog but likely harder to read
var noparam = 89

func (p *Params) Geti(str string) int {
	if reti, ok := p.ints[str]; ok {
		return reti
	} else if rawstr, ok := p.raw[str]; ok { // check in raw hash
		ret, err := strconv.Atoi(rawstr) // convert to integer
		if err != nil {
			fmt.Println(err)
		}
		p.ints[str] = ret // add to ints map
		return ret
	} else {
		hardExit(noparam, str)
	} // wasn't found anywhere
	return 0
}

func (p *Params) Getn(str string) float64 {
	if retn, ok := p.nums[str]; ok {
		return retn
	} else if rawstr, ok := p.raw[str]; ok { // check in raw hash
		ret, err := strconv.ParseFloat(rawstr, 64) // convert to float
		if err != nil {
			fmt.Println(err)
		}
		p.nums[str] = ret // add to nums map
		return ret
	} else {
		hardExit(noparam, str)
	} // wasn't found anywhere
	return 0
}

func (p *Params) Gets(str string) string {
	if retn, ok := p.strs[str]; ok {
		return retn
	} else if retn, ok := p.raw[str]; ok { // check in raw hash
		p.strs[str] = retn // add to string map
		return retn
	} else {
		hardExit(noparam, str)
	} // wasn't found anywhere
	return ""
}

func (p *Params) Getf(str string) string {
	fmt.Println("\nGetting file", str)
	if retn, ok := p.files[str]; ok {
		fmt.Println("\nfirst", str, p.raw)
		return retn
	} else if retn, ok := p.raw[str]; ok { // check in raw hash
		fmt.Println("\nraw", str)
		p.files[str] = retn // add to files map
		return retn
	} else {
		fmt.Println("\nfailed file", str, p.raw, p.raw[str])
		hardExit(noparam, str)
	} // wasn't found anywhere
	return ""
}

func (p *Params) Getb(str string) bool {
	if retb, ok := p.bools[str]; ok {
		return retb
	} else if rawstr, ok := p.raw[str]; ok { // check in raw hash
		fmt.Println("Getting a boolean from string\n")
		retb, err := strconv.ParseBool(rawstr) // convert to boolean
		// 1, t, T, true, TRUE , True all accepted
		if err != nil {
			fmt.Println(err)
		}
		p.bools[str] = retb // add to bools map
		return retb
	} else {
		hardExit(noparam, str)
	} // wasn't found anywhere
	return false
}

func (p *Params) getr(str string) string { // this is lower case because we shouldn't export it'
	if retn, ok := p.raw[str]; ok {
		return retn
	} else {
		hardExit(noparam, str)
	} // wasn't found anywhere
	return ""
}

// A data structure to hold a key/value pair.
type Pair struct {
	Key   string
	Value int
}

// deleted old way of doing it to not distract; check that this still works!
// A slice of Pairs that implements sort.Interface to sort by Value.
type ByValue []Pair

func (p ByValue) Len() int           { return len(p) }
func (p ByValue) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }
func (p ByValue) Less(i, j int) bool { return p[i].Value < p[j].Value }

// A slice of Pairs that implements sort.Interface to sort by Value.
type ByKey []Pair

func (p ByKey) Len() int           { return len(p) }
func (p ByKey) Swap(i, j int)      { p[i], p[j] = p[j], p[i] }
func (p ByKey) Less(i, j int) bool { return p[i].Value < p[j].Value }

// SortbyVals turns a map into a PairList, then sorts and returns sorted keys.
func SortbyVals(m map[string]int) []string {
	p := make(ByValue, len(m))
	i := 0
	for k, v := range m {
		p[i] = Pair{k, v}
		i++
	} // increment i to point to pairs
	sort.Slice(p, func(i, j int) bool { return p[i].Value > p[j].Value }) // sort by Value
	var sortkeys []string                                                 // to hold sorted keys
	for i := range p {
		sortkeys = append(sortkeys, p[i].Key)
	} //  append each Key in order to sortkeys
	return sortkeys
}

// hardExit prints blurb and flag, exits program
func hardExit(flag int, param string) {
	fmt.Println("\n\nBad parameter", param, "!! Hard Exit!", flag, "\n")
	os.Exit(flag)
}

// Check freaks out if there was an error
func Check(e error) {
	if e != nil {
		panic(e)
	}
}
