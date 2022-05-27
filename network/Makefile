# Simple demonstration of script to build network figure from AnVir output.

OUTPUT=haplo_graph.example.png

${OUTPUT} : resources/haplotype_network_example.1000cmin.xls
	python3 scripts/edges_to_figure.py $< $@

.PHONY: clean
clean :
	rm -f ${OUTPUT}
