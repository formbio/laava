# Local tests: Run the core analysis steps on downsampled copies of the public PacBio example data

in_dir = samples
out_dir = build
NCPU = 4

.PHONY: all clean sc ss test
all: sc ss

clean:
	rm -frv build/*

sc: $(out_dir)/sc_AAV_report.html

ss: $(out_dir)/ss_AAV_report.html

test: sc ss
	pytest test_outputs.py


$(out_dir)/sc_AAV_report.html: ../src/create_report.R ../src/report.Rmd $(out_dir)/sc.annotation.txt $(out_dir)/sc.metadata.tsv $(out_dir)/sc.tagged.bam
	$< $(out_dir)/sc sc $(word 3,$^)

$(out_dir)/ss_AAV_report.html: ../src/create_report.R ../src/report.Rmd $(out_dir)/ss.annotation.txt $(out_dir)/ss.metadata.tsv $(out_dir)/ss.flipflop.tsv.gz
	$< $(out_dir)/ss ss $(word 3,$^)


$(out_dir)/sc.metadata.tsv: ../src/write_sample_metadata.py $(in_dir)/sc.subsample005.bam
	mkdir -p $(out_dir)
	$< sc.subsample005 "Self-comp Example" $(word 2,$^) -o $@

$(out_dir)/ss.metadata.tsv: ../src/write_sample_metadata.py $(in_dir)/ss.subsample005.bam
	mkdir -p $(out_dir)
	$< ss.subsample005 "Single-stranded Example" $(word 2,$^) -o $@


$(out_dir)/ss.flipflop.tsv.gz: ../src/get_flipflop_config.py $(out_dir)/ss.tagged.bam $(out_dir)/ss.per_read.tsv.gz
	$^ -o $(out_dir)/ss


$(out_dir)/sc.tagged.bam: ../src/summarize_alignment.py $(in_dir)/sc.subsample005.bam $(out_dir)/sc.annotation.txt
	$^ $(out_dir)/sc --sample-id sc.subsample005 --vector-type sc --cpus $(NCPU)

$(out_dir)/ss.tagged.bam: ../src/summarize_alignment.py $(in_dir)/ss.subsample005.bam $(out_dir)/ss.annotation.txt
	$^ $(out_dir)/ss --sample-id ss.subsample005 --vector-type ss --cpus $(NCPU)


$(out_dir)/sc.annotation.txt: ../src/prepare_annotation.py $(in_dir)/sc.annotation.bed $(out_dir)/sc.reference_names.tsv
	$^ ITR-L ITR-R -o $@

$(out_dir)/ss.annotation.txt: ../src/prepare_annotation.py $(in_dir)/ss.annotation.bed $(out_dir)/ss.reference_names.tsv
	$^ ITR -o $@


$(out_dir)/sc.reference_names.tsv $(out_dir)/ss.reference_names.tsv: \
	$(out_dir)/%.reference_names.tsv: \
	samples/%.construct.fasta fasta/packaging.fa fasta/hg38.chr19trunc-chrM.fa
	mkdir -p $(out_dir)
	../src/get_reference_names.py -o $@ \
		$< --packaging $(word 2,$^) --host $(word 3,$^) \
		--repcap-name pRep2Cap9 --helper-name pHelper --lambda-name Lambda

