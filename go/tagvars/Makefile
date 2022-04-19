# An example of how to build and run TagVars.

outputs/variants.txt : tagvars inputs/covid_seqs_Jan22.50k.fasta inputs/ref_Wuhan_Oct20.kcounts
	./tagvars

tagvars : tagvars.go
	go build tagvars.go

.PHONY : clean
clean :
	rm tagvars
	rm outputs/*
