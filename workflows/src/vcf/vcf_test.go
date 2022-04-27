package vcf_test

import (
	"bytes"
	"os"
	"io"
	"io/ioutil"
	"reflect"
	"testing"

	"annotation/vcf"
	. "annotation/utils"
)

func compare[T any](result T, correct T, t *testing.T) {
	if !reflect.DeepEqual(result, correct) {
		t.Errorf("RESULT:\n%+v\nCORRECT:\n%+v\n", result, correct)
	}
}

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
		AddInfo("TEST1", "hello", "world").
		AddInfo("TEST2", "2").
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

func TestParseVCFRecord(t *testing.T) {
	line := "NC_045512.2 2944	3	G	A	.	.	TYPE=SNP;END=2944;COUNT=19745;KMERS=TTACTTACACCACT,TACTTACACCACTA,ACTTACACCACTAG,CTTACACCACTAGG,TTACACCACTAGGC,TACACCACTAGGCA,ACACCACTAGGCAT,CACCACTAGGCATT,ACCACTAGGCATTG,CCACTAGGCATTGA,CACTAGGCATTGAT,ACTAGGCATTGATT,CTAGGCATTGATTT,TAGGCATTGATTTA,AGGCATTGATTTAG,GGCATTGATTTAGA"
	result := vcf.ParseVCFRecord(line)
	correct := vcf.VcfRecord().
		SetChrom("NC_045512.2").
		SetPos(2944).
		SetID("3").
		SetRef("G").
		SetAlt("A").
		SetQual(".").
		SetFilter(".").
		AddInfo("TYPE", "SNP").
		AddInfo("END", "2944").
		AddInfo("COUNT", "19745").
		AddInfo("KMERS", 
			"TTACTTACACCACT", "TACTTACACCACTA", "ACTTACACCACTAG",
			"CTTACACCACTAGG", "TTACACCACTAGGC", "TACACCACTAGGCA",
			"ACACCACTAGGCAT", "CACCACTAGGCATT", "ACCACTAGGCATTG",
			"CCACTAGGCATTGA", "CACTAGGCATTGAT", "ACTAGGCATTGATT",
			"CTAGGCATTGATTT", "TAGGCATTGATTTA", "AGGCATTGATTTAG",
			"GGCATTGATTTAGA")
	compare(*result, *correct, t)
}
