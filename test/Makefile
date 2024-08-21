# Local tests: Run the core analysis steps on downsampled copies of the public PacBio scAAV data

in_dir = samples
out_dir = build
NCPU = 4

.PHONY: all sc ss clean
all: sc ss

clean:
	rm -frv build/*

sc: $(out_dir)/sc_AAV_report.html

ss: $(out_dir)/ss_AAV_report.html


test: $(out_dir)/sc.Rdata $(out_dir)/ss.Rdata
	./compare_tsvs.sh build-snapshot build \
		sc.alignments.tsv sc.per_read.tsv ss.alignments.tsv ss.per_read.tsv ss.flipflop.tsv \
		&& echo "OK" || echo "FAIL"

$(out_dir)/sc_AAV_report.html $(out_dir)/ss_AAV_report.html: \
	%_AAV_report.html: \
	../src/create_report.R %.Rdata
	$^


$(out_dir)/sc.Rdata: ../src/calculate_rdata.R ../src/report.Rmd $(out_dir)/sc.annotation.txt $(out_dir)/sc.tagged.bam
	$< $(out_dir)/sc $(out_dir)/sc.annotation.txt sc.subsample005 sc.subsample005 "scaav"

$(out_dir)/ss.Rdata: ../src/calculate_rdata.R ../src/report.Rmd $(out_dir)/ss.annotation.txt $(out_dir)/ss.flipflop_assignments.tsv $(out_dir)/ss.tagged.bam
	$< $(out_dir)/ss $(out_dir)/ss.annotation.txt ss.subsample005 "ssaav" $(out_dir)/ss.flipflop_assignments.tsv


$(out_dir)/ss.flipflop_assignments.tsv: ../src/get_flipflop_config.py $(out_dir)/ss.tagged.bam $(out_dir)/ss.per_read.tsv
	$^ -o $(out_dir)/ss


$(out_dir)/sc.tagged.bam: ../src/summarize_alignment.py $(in_dir)/sc.subsample005.bam $(out_dir)/sc.annotation.txt
	$^ $(out_dir)/sc --sample-id sc.subsample005 --vector-type scaav --cpus $(NCPU)

$(out_dir)/ss.tagged.bam: ../src/summarize_alignment.py $(in_dir)/ss.subsample005.bam $(out_dir)/ss.annotation.txt
	$^ $(out_dir)/ss --sample-id ss.subsample005 --vector-type ssaav --cpus $(NCPU)


$(out_dir)/sc.annotation.txt: ../src/prepare_annotation.py $(in_dir)/sc.annotation.bed $(out_dir)/sc.reference_names.tsv
	$^ ITR-L ITR-R -o $@

$(out_dir)/ss.annotation.txt: ../src/prepare_annotation.py $(in_dir)/ss.annotation.bed $(out_dir)/ss.reference_names.tsv
	$^ ITR -o $@


$(out_dir)/sc.reference_names.tsv $(out_dir)/ss.reference_names.tsv: \
	$(out_dir)/%.reference_names.tsv: \
	samples/%.construct.fasta fasta/packaging.fa fasta/hg38.chr19trunc-chrM.fa
	../src/get_reference_names.py -o $@ \
		$< --packaging $(word 2,$^) --host $(word 3,$^) \
		--repcap-name pRep2Cap9 --helper-name pHelper --lambda-name Lambda

