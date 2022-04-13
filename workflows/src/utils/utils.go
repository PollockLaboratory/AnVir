package utils

// quick error checking pattern
func Check(e error) {
	if e != nil {
		panic(e)
	}
}
