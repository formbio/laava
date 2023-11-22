# Slow tests: Run the core 3 analysis steps on 2 public datasets

sc_inputs = CBA-eGFP
ss_inputs = CMV-GFP
out_sc_pfx = build/scAAV
out_ss_pfx = build/ssAAV


all:\
	$(out_sc_pfx).per_read.csv \
	$(out_ss_pfx).per_read.csv \
	$(out_sc_pfx).flipflop_assignments.txt \
	$(out_ss_pfx).flipflop_assignments.txt \
	$(sc_inputs)/scAAV_bc1009_AAV_report.pdf \
	$(ss_inputs)/2022ssAAV_AAV_report.pdf


$(ss_inputs)/2022ssaav_aav_report.pdf: ../plotaavreport.r
	Rscript $< $(ss_inputs)/2022ssaav $(ss_inputs)/annotation.txt

$(out_ss_pfx).flipflop_assignments.txt: ../get_flipflop_config.py $(ss_inputs)/2022ssaav.tagged.sorted.bam $(ss_inputs)/2022ssaav.per_read.csv
	python3 $^ -o $(out_ss_pfx)

$(out_ss_pfx).per_read.csv: ../summarize_aav_alignment.py $(ss_inputs)/mapped.sort_by_read.sam $(ss_inputs)/annotation.txt
	python3 $^ $(out_ss_pfx) --cpus 4


$(sc_inputs)/scAAV_bc1009_AAV_report.pdf: ../plotaavreport.r
	Rscript $< $(sc_inputs)/scAAV_bc1009 $(sc_inputs)/annotation.txt

$(out_sc_pfx).flipflop_assignments.txt: ../get_flipflop_config.py $(sc_inputs)/scAAV_bc1009.tagged.sorted.bam $(sc_inputs)/scAAV_bc1009.per_read.csv
	python3 $^ -o $(out_sc_pfx)

$(out_sc_pfx).per_read.csv: ../summarize_aav_alignment.py $(sc_inputs)/mapped.sort_by_read.sam $(sc_inputs)/annotation.txt
	python3 $^ $(out_sc_pfx) --cpus 4


# Download example input data from PacBio
# Upstream source: https://downloads.pacbcloud.com/public/dataset/AAV/

ss_base_url = https://downloads.pacbcloud.com/public/dataset/AAV/2022-ssAAV-pAV-CMF-GFP/

$(ss_inputs)/annotation.txt:
	wget $(ss_base_url)/1-mapped/annotation.txt -O $@

$(ss_inputs)/mapped.sort_by_read.sam:
	wget $(ss_base_url)/1-mapped/mapped.sort_by_read.sam -O $@

$(ss_inputs)/2022ssAAV.tagged.sorted.bam:
	wget $(ss_base_url)/1-mapped/2022ssAAV.tagged.sorted.bam -O $@
	wget $(ss_base_url)/1-mapped/2022ssAAV.tagged.sorted.bam.bai -O $@.bai

$(ss_inputs)/2022ssAAV.nonmatch_stat.csv.gz:
	wget $(ss_base_url)/1-mapped/2022ssAAV.nonmatch_stat.csv.gz -O $@

$(ss_inputs)/2022ssAAV.per_read.csv:
	wget $(ss_base_url)/1-mapped/2022ssAAV.per_read.csv -O $@

$(ss_inputs)/2022ssAAV.summary.csv:
	wget $(ss_base_url)/1-mapped/2022ssAAV.summary.csv -O $@

$(ss_inputs)/pAV-CMV-GFP.fasta:
	wget $(ss_base_url)/pAV-CMV-GFP.fasta -O $@


sc_base_url = https://downloads.pacbcloud.com/public/dataset/AAV/2021-scAAV-CBA-eGFP/

$(sc_inputs)/annotation.txt:
	wget $(sc_base_url)/2-mapped/annotation.txt -O $@

$(sc_inputs)/mapped.sort_by_read.sam:
	wget $(sc_base_url)/2-mapped/demux.bc1009--bc1009.mapped.sort_by_read.sam -O $@

$(sc_inputs)/scAAV_bc1009.tagged.sorted.bam:
	wget $(sc_base_url)/2-mapped/scAAV_bc1009.tagged.bam -O $@
	wget $(sc_base_url)/2-mapped/scAAV_bc1009.tagged.bam.bai -O $@.bai

$(sc_inputs)/scAAV_bc1009.nonmatch_stat.csv.gz:
	wget $(sc_base_url)/2-mapped/scAAV_bc1009.nonmatch_stat.csv.gz -O $@

$(sc_inputs)/scAAV_bc1009.per_read.csv:
	wget $(sc_base_url)/2-mapped/scAAV_bc1009.per_read.csv -O $@

$(sc_inputs)/scAAV_bc1009.summary.csv:
	wget $(sc_base_url)/2-mapped/scAAV_bc1009.summary.csv -O $@

$(sc_inputs)/scAAV_CBA_eGFP.fasta:
	wget $(sc_base_url)/scAAV_CBA_eGFP.fasta -O $@
