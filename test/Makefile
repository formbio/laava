# Local tests: Run the core 3 analysis steps on a downsampled copy of the public PacBio scAAV data

small_inputs = small
out_dir = build
out_small_pfx = $(out_dir)/small


.PHONY: fast slow clean
fast: $(out_small_pfx)_AAV_report.html

slow: $(out_small_pfx).tagged.bam $(out_small_pfx).flipflop_assignments.tsv $(out_small_pfx)_AAV_report.html

clean:
	rm -vf build/*
	rm -vf $(out_small_pfx)_AAV_report.html $(out_small_pfx)_AAV_report.pdf

diffcheck:
	diff build-snapshot/small.per_read.tsv build/small.per_read.tsv && echo "OK"
	diff build-snapshot/small.flipflop.tsv build/small.flipflop.tsv && echo "OK"

$(out_small_pfx)_AAV_report.html: ../src/create_report.R $(out_small_pfx).Rdata
	$^

$(out_small_pfx).Rdata: ../src/calculate_rdata.R ../src/report.Rmd $(out_small_pfx).annotation.txt $(out_small_pfx).flipflop_assignments.tsv
	$< $(out_small_pfx) $(small_inputs)/annotation.txt pAV-CMV-GFP "scaav" $(out_small_pfx).flipflop_assignments.tsv

$(out_small_pfx).flipflop_assignments.tsv: ../src/get_flipflop_config.py $(out_small_pfx).tagged.bam $(out_small_pfx).per_read.tsv
	$^ -o $(out_small_pfx)

$(out_small_pfx).annotation.txt: ../src/prepare_annotation.py $(small_inputs)/annotation.bed $(small_inputs)/reference_names.tsv
	mkdir -p $(out_dir)
	$^ -o $@

$(out_small_pfx).tagged.bam: ../src/summarize_AAV_alignment.py $(small_inputs)/scaav-cba-egfp.bam $(small_inputs)/annotation.txt
	$^ $(out_small_pfx) --cpus 4
