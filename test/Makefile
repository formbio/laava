# Local tests: Run the core analysis steps on downsampled copies of the public PacBio scAAV data

in_dir = samples
out_dir = build
NCPU = 4

.PHONY: all sc ss clean
all: sc ss

clean:
	rm -vf build/*

sc: $(out_dir)/sc_AAV_report.html

ss: $(out_dir)/ss_AAV_report.html


diffcheck:
	diff build-snapshot/sc.per_read.tsv build/sc.per_read.tsv && echo "OK"
	diff build-snapshot/ss.per_read.tsv build/ss.per_read.tsv && echo "OK"
	diff build-snapshot/ss.flipflop.tsv build/ss.flipflop.tsv && echo "OK"

$(out_dir)/sc_AAV_report.html $(out_dir)/ss_AAV_report.html: \
	%_AAV_report.html: \
	../src/create_report.R %.Rdata
	$^


$(out_dir)/sc.Rdata: ../src/calculate_rdata.R ../src/report.Rmd $(out_dir)/sc.annotation.txt $(out_dir)/sc.tagged.bam
	$< $(out_dir)/sc $(out_dir)/sc.annotation.txt pAV-CMV-GFP "scaav"

$(out_dir)/ss.Rdata: ../src/calculate_rdata.R ../src/report.Rmd $(out_dir)/ss.annotation.txt $(out_dir)/ss.flipflop_assignments.tsv $(out_dir)/ss.tagged.bam
	$< $(out_dir)/ss $(out_dir)/ss.annotation.txt pAV-CMV-GFP "scaav" $(out_dir)/ss.flipflop_assignments.tsv


$(out_dir)/ss.flipflop_assignments.tsv: ../src/get_flipflop_config.py $(out_dir)/ss.tagged.bam $(out_dir)/ss.per_read.tsv
	$^ -o $(out_dir)/ss


$(out_dir)/sc.tagged.bam: ../src/summarize_AAV_alignment.py $(in_dir)/sc.subsample005.bam $(out_dir)/sc.annotation.txt
	$^ $(out_dir)/sc --cpus $(NCPU)

$(out_dir)/ss.tagged.bam: ../src/summarize_AAV_alignment.py $(in_dir)/ss.subsample005.bam $(out_dir)/ss.annotation.txt
	$^ $(out_dir)/ss --cpus $(NCPU)


$(out_dir)/sc.annotation.txt $(out_dir)/ss.annotation.txt: \
	$(out_dir)/%.annotation.txt: \
	../src/prepare_annotation.py $(in_dir)/%.annotation.bed $(out_dir)/%.reference_names.tsv
	$^ -o $@


$(out_dir)/sc.reference_names.tsv $(out_dir)/ss.reference_names.tsv: \
	$(out_dir)/%.reference_names.tsv: \
	samples/%.construct.fasta fasta/packaging.fa fasta/hg38.chr19trunc-chrM.fa
	./get_reference_names.sh $^ "pRep2Cap9" $@

