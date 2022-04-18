package vcf_test

import (
	"bytes"
	"os"
	"io"
	"io/ioutil"
	"testing"

	"annotation/vcf"
	. "annotation/utils"
)

// wrapper function to capture stdout into a string
// Not thread safe, but in this use case it's fine.
func captureStdout(f func()) string {
  old := os.Stdout
  r, w, _ := os.Pipe()
  os.Stdout = w

  f()

  w.Close()
  os.Stdout = old

  var buf bytes.Buffer
  io.Copy(&buf, r)
  return buf.String()
}

func generateTestOutput() {
	vcf.VcfHeader().
		AddReference("blah").
		AddInfo("TEST1", ".", "String", "test1 field").
		AddInfo("TEST2", "1", "Integer", "test2 field").
		AddContig("testchr", 10).
		Write(os.Stdout)
	vcf.VcfRecord().
		SetChrom("testchr").SetPos(1).SetID("1").
		SetRef("A").SetAlt("T").
		SetQual(".").SetFilter(".").
		AddInfo("TEST1", []string{"hello", "world"}).
		AddInfo("TEST2", []string{"2"}).
		Write(os.Stdout)
}

func TestVcfCreation(t *testing.T) {
	correct, err := ioutil.ReadFile("test.vcf")
	Check(err)

	result := captureStdout(func () {
		generateTestOutput()
	})
	if string(correct) != result {
		t.Fatalf("RESULT:\n%s\nCORRECT:\n%s\n", result, correct)
	}
}
