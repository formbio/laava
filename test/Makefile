# Local tests: Run the core analysis steps on downsampled copies of the public PacBio example data

in_dir = samples
out_dir = build
export PATH := $(CURDIR)/../src:$(PATH)


.PHONY: all clean sc ss test test-local sc-local ss-local
all: sc ss

clean:
	rm -frv build/*

test: sc ss
	pytest test_outputs.py

# Local testing without Docker (using conda)
test-local: sc-local ss-local
	conda run -n laava python -m pytest test_outputs.py test_integration_duplicate_removal.py test_unit_laava_modules.py test_output_gold_verification.py -v

sc-local: $(out_dir)/sc.reference_names.tsv
	mkdir -p $(out_dir)
	conda run -n laava bash ../src/make_report.sh \
		sc \
		"Self-comp Example" \
		test \
		$(out_dir)/sc.reference_names.tsv \
		$(in_dir)/sc.subsample005.bam \
		$(in_dir)/sc.annotation.bed \
		wtITR \
		"" \
		mITR \
		sc \
		200 \
		100 \
		100 \
		0.8 \
		"" \
		"" \
		$(out_dir)

ss-local: $(out_dir)/ss.reference_names.tsv
	mkdir -p $(out_dir)
	conda run -n laava bash ../src/make_report.sh \
		ss \
		"Single-stranded Example" \
		test \
		$(out_dir)/ss.reference_names.tsv \
		$(in_dir)/ss.subsample005.bam \
		$(in_dir)/ss.annotation.bed \
		ITR \
		"" \
		"" \
		ss \
		200 \
		100 \
		100 \
		0.8 \
		AAV2 \
		"" \
		$(out_dir)

tc-gia-012-local: $(out_dir)/tc-gia-012.reference_names.tsv
	mkdir -p $(out_dir)
	conda run -n laava bash ../src/make_report.sh \
		tc-gia-012 \
		"TC-GIA-012" \
		test \
		$(out_dir)/tc-gia-012.reference_names.tsv \
		tc-gia-012_aligned.bam \
		$(in_dir)/TC-GIA-012/scAAV_CBA_eGFP_reverse_complement.bed \
		ITR-L \
		"" \
		mITR \
		sc \
		200 \
		100 \
		100 \
		0.8 \
		AAV2 \
		$(in_dir)/TC-GIA-012/AAV2.flipflop.fasta \
		$(out_dir)

sc: $(out_dir)/sc_AAV_report.html

ss: $(out_dir)/ss_AAV_report.html


$(out_dir)/sc_AAV_report.html: \
	../src/make_report.sh $(out_dir)/sc.reference_names.tsv \
		$(in_dir)/sc.subsample005.bam $(in_dir)/sc.annotation.bed \
		../src/aggregate_tables.py ../src/create_report.R \
		../src/emit_target_coords.py ../src/get_flipflop_config.py \
		../src/prepare_annotation.py ../src/summarize_alignment.py \
		../src/write_sample_metadata.py
	mkdir -p $(out_dir)
	$< \
		sc \
		"Self-comp Example" \
		test \
		$(out_dir)/sc.reference_names.tsv \
		$(in_dir)/sc.subsample005.bam \
		$(in_dir)/sc.annotation.bed \
		wtITR \
		"" \
		mITR \
		sc \
		200 \
		100 \
		100 \
		0.8 \
		"" \
		"" \
		$(out_dir)


$(out_dir)/ss_AAV_report.html: \
	../src/make_report.sh $(out_dir)/ss.reference_names.tsv \
		$(in_dir)/ss.subsample005.bam $(in_dir)/ss.annotation.bed \
		../src/aggregate_tables.py ../src/create_report.R \
		../src/emit_target_coords.py ../src/get_flipflop_config.py \
		../src/prepare_annotation.py ../src/summarize_alignment.py \
		../src/write_sample_metadata.py
	mkdir -p $(out_dir)
	$< \
		ss \
		"Single-stranded Example" \
		test \
		$(out_dir)/ss.reference_names.tsv \
		$(in_dir)/ss.subsample005.bam \
		$(in_dir)/ss.annotation.bed \
		ITR \
		"" \
		"" \
		ss \
		200 \
		100 \
		100 \
		0.8 \
		AAV2 \
		"" \
		$(out_dir)

$(out_dir)/sc.reference_names.tsv $(out_dir)/ss.reference_names.tsv: \
	$(out_dir)/%.reference_names.tsv: \
	$(in_dir)/%.construct.fasta fasta/packaging.fa fasta/hg38.chr19trunc-chrM.fa
	mkdir -p $(out_dir)
	conda run -n laava python ../src/get_reference_names.py -o $@ \
		$< --packaging $(word 2,$^) --host $(word 3,$^) \
		--repcap-name pRep2Cap9 --helper-name pHelper --lambda-name Lambda
