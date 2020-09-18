# ATAC-Seq footprint pipeline. Inspired by HINT-ATAC algorithm, authored by Zhijian Li in the Costa Lab.
# Thanks to Rowan Callahan, Jake Van Campen, and Joey Estabrook, for their code.

# How to use
# 1. Symlink your raw sequencing files into samples/raw.
# 2. Merge your sequencing files by lane, if applicable, before starting the pipeline.

# Assumptions
# Raw files in the format: "{cond}{replicate}_{R1|R2}.fastq.gz" e.g. MOLMD1_R1.fastq.gz where MOLMD is cond and 1 is replicate.
# Replicate ID (if a number) must not be next to another number. So be extra careful if time-series data.
# Same number of replicates for all conditions.
# Paired-end sequencing.
# You've installed RGT to your user directory and configured genomic data for your organism.

import os
import glob
import plotly as plt
import pandas as pd
from datetime import date

SAMPLE, DIR, = glob_wildcards("samples/raw/{sample}_{dir,R1|R2}.fastq.gz")
SAMPLE = sorted(set(SAMPLE))
DIR = sorted(set(DIR))
# SAMPLE wildcard aggregates condition and replicates.

configfile: "config.yaml"
localrules: filter_peaks, consensus_intervals, rm_blacklist, fix_motifs, sort_stats

leuk_cells = ["jurkat-hg38", "k-562-hg38", "kasumi-1-hg38", "MOLM-13-hg38"]

rule all:
	input:
		# trim reads, fastq_screen, fastqc -----------------------------
		expand(["samples/trim/{sample}_R1_val_1.fq.gz",
				"samples/trim/{sample}_R2_val_2.fq.gz",
				"samples/trim/{sample}_R1.fastq.gz_trimming_report.txt",
				"samples/trim/{sample}_R2.fastq.gz_trimming_report.txt"],
				sample = SAMPLE),
		expand(["samples/fastq_screen/{sample}_R1_val_1_screen.{ext}",
				"samples/fastq_screen/{sample}_R2_val_2_screen.{ext}"], 
				sample = SAMPLE, ext = ["html", "txt", "png"]),
		expand(["samples/fastqc/{sample}_R1_val_1_fastqc.{ext}",
				"samples/fastqc/{sample}_R2_val_2_fastqc.{ext}"], 
				sample = SAMPLE, ext = ["html", "zip"]),
		# align reads, get fragment length distribution -----------------
		"samples/align/fragment_length/fragment_length_dist.html",
		# remove mitochondrial reads, deduplicate, quality reads --------
		expand("samples/align/quality_align/{sample}_dedup_rmchrM_quality.bam", sample = SAMPLE),
		# merge bam files -----------------------------------------------
		expand("samples/align/merged/{cond}_merged.{ext}", cond = list(config["CONDITIONS"]), ext = ["bam", "bam.bai"] ),
		# bigwig tracks, call + filter ATAC-Seq peaks, take consensus peaks
		expand("samples/align/bigwig/{cond}.bw", cond = list(config["CONDITIONS"]) ),
		expand("data/macs/filtered/{cond}_filtered.broadPeak", cond = list(config["CONDITIONS"]) ),
		"data/macs/consensus_intervals.bed",
		# footprint + motif match + tracks ------------------------------
		expand("data/differential/match/{cond}_mpbs.bed", cond = list(config["CONDITIONS"]) ),
		expand("data/differential/match/fixed/{cond}_mpbs.bed", cond = list(config["CONDITIONS"]) ),
		expand("data/differential/tracks/{cond}.wig", cond = list(config["CONDITIONS"]) ),
		# differential chromatin oppenness testing, sort results --------
		expand("data/differential/diff_pseudocounts/{control}_{case}/{diff_outputs}",
			control = config["CONTROL"], 
			case =  list( set(config["CONDITIONS"]) - set(config["CONTROL"]) ),
			diff_outputs = ["differential_factor.txt", "differential_statistics.txt", "differential_statistics.pdf", "differential_statistics_sorted_filtered.txt"]),
		# downstream analysis 
		expand("data/differential/chipAtlas/{cond}_{chips}.bed", cond = list(config["CONDITIONS"]), chips = leuk_cells)

# snakemake -j 64 --use-conda --rerun-incomplete --latency-wait 60 --keep-going --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} -N {cluster.N} -o {cluster.o} -e {cluster.e} -t {cluster.t} -J {cluster.J} -c {threads} --mem={cluster.mem}" -s Snakefile

# snakemake -j 64 --use-conda --rerun-incomplete --latency-wait 60 --keep-going --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -N {cluster.N}  -t {cluster.t} -J {cluster.J} -c {cluster.c} --mem={cluster.mem}" -s Snakefile

# quality control -------------------------------------------------------------------------
if config["SEQUENCER"] == "not_illumina":
	adapters = "-a 'file:{}'".format(config["ADAPTER"])
elif config["SEQUENCER"] == "illumina":
	adapters = "--illumina"
else:
	print("Must define sequencer for read trimming: ['not_illumina', 'illumina']. Your input = {}".format(config["SEQUENCER"]))

rule trim_galore_pe:
	input:
		r1 = "samples/raw/{sample}_R1.fastq.gz",
		r2 = "samples/raw/{sample}_R2.fastq.gz"
	output:
		"samples/trim/{sample}_R1.fastq.gz_trimming_report.txt",
		"samples/trim/{sample}_R2.fastq.gz_trimming_report.txt",
		r1 = "samples/trim/{sample}_R1_val_1.fq.gz",
		r2 = "samples/trim/{sample}_R2_val_2.fq.gz"
	message: " -- Trimming {wildcards.sample} -- "
	log:
		"logs/trim_galore/{sample}.log"
	conda:
		"envs/trim_galore.yaml"
	params:
		conf = adapters
	threads: 4
	shell:
		"(trim_galore -j {threads} {params.conf} -q 30 --paired --trim1 --paired -o samples/trim {input.r1} {input.r2}) > {log} 2>&1"

rule fastq_screen:
	input:
		r1 = rules.trim_galore_pe.output.r1,
		r2 = rules.trim_galore_pe.output.r2
	output:
		"samples/fastq_screen/{sample}_R1_val_1_screen.txt",
		"samples/fastq_screen/{sample}_R1_val_1_screen.png",
		"samples/fastq_screen/{sample}_R1_val_1_screen.html",
		"samples/fastq_screen/{sample}_R2_val_2_screen.txt",
		"samples/fastq_screen/{sample}_R2_val_2_screen.png",
		"samples/fastq_screen/{sample}_R2_val_2_screen.html"
	threads: 8
	params:
		conf = config["FASTQC_SCREEN_CONF"]
	conda:
		"envs/fastq_screen.yaml"
	message:
		" -- Fastq-screening {wildcards.sample} -- "
	log:
		"logs/fastq_screen/{sample}.log"
	shell:
		"(fastq_screen --aligner bowtie2 --conf {params.conf} --outdir samples/fastq_screen --threads {threads} {input.r1} {input.r2}) 2>{log}"

# fastq_screen somehow not working with scheduler.

rule fastqc:
	input:
		r1 = rules.trim_galore_pe.output.r1,
		r2 = rules.trim_galore_pe.output.r2
	output:
		"samples/fastqc/{sample}_R1_val_1_fastqc.html",
		"samples/fastqc/{sample}_R2_val_2_fastqc.html",
		"samples/fastqc/{sample}_R1_val_1_fastqc.zip",
		"samples/fastqc/{sample}_R2_val_2_fastqc.zip",
	conda:
		"envs/fastqc.yaml"
	message:
		" -- FastQC'ing {wildcards.sample} -- "
	log:
		"logs/fastqc/{sample}.log"
	threads: 4
	shell:
		"fastqc --outdir samples/fastqc --extract -f fastq {input.r1} {input.r2} -t {threads}"

# rule multiqc:
# 	input:
# 		expand("samples/fastqc/{sample}.html", sample=SAMPLE)
# 	output:
# 		"data/qc/multiqc_report.html"
# 	wrapper:
# 		"0.31.1/bio/multiqc"

# align reads -------------------------------------------------------------------
rule bowtie2:
	input:
		fw = rules.trim_galore_pe.output.r1,
		rv = rules.trim_galore_pe.output.r2
	output:
		"samples/align/{sample}.bam"
	message: "  -- Aligning {wildcards.sample} --  "
	params:
		genome_index = config["GENOME_INDEX"]
	threads: 12
	conda:
		"envs/align.yaml"
	log: "logs/bowtie2/{sample}.err"
	shell:
		"(bowtie2 \
		--very-sensitive -X 2000 --no-mixed --no-discordant \
		--mm -p {threads} \
		-x {params.genome_index} \
		-1 {input.fw} -2 {input.rv} | samtools view -bS - > {output}) 2> {log}"
# stdout is the SAM file. stderr is log file.
# -X 2000 = limit fragment length < 2000bp
# -no-mixed = 
# -no-discordant = 
# pipe output to bam format

# qc metrics -------------------------------------------------------------------

# QC metrics: fragment length distribution, read mapping proportion per chromosome
rule fragment_length:
	input:
		rules.bowtie2.output
	output:
		"samples/align/fragment_length/{sample}.txt"
	conda:
		"envs/align.yaml"
	threads: 4
	shell:
		"samtools view {input} | cut -f9 | awk '$1 > 0' | sort --parallel={threads} -S 50% | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > {output}"
# output fragment length info per sample.
# cut -f9 = fragment length in bam file
# awk '$1 > 0' | sort | uniq -c = get positive fragment length, sort, count the lengths
# sort -b -k2,2n = col1 is fragment size, col2 is counts. sort by counts numerically.
# sed = remove tabs

rule fragment_length_plot:
	input:
		expand("samples/align/fragment_length/{sample}.txt", sample = SAMPLE)
	output:
		"samples/align/fragment_length/fragment_length_dist.html"
	run:
		pd.options.plotting.backend = "plotly"
		df = pd.concat([pd.read_csv(f, index_col = 1, sep = " ", names = [f.split('/')[3][:-4]]) for f in input], axis = 1)
		print("Fragment Length Distribution Table")
		print(df)
		fragment_length_dist = df.plot()
		fragment_length_dist.update_layout( 
			title='Fragment Length Distribution', 
			xaxis_title='Fragment Length (bp)', 
			yaxis_title='Counts', 
			legend_title_text='Samples')
		fragment_length_dist.write_html(str(output))
# read in info as list of df. Concat everything into one df (row = length, col = sample, elements = counts of frag length), plot and label with plotly. 

# process alignment ------------------------------------------------------------------------

# sort, index, deduplicate PCR reads, remove mitochondrial reads, filter for high-quality alignments
rule sort_bams1:
	input:
		rules.bowtie2.output
	output:
		temp("samples/align/{sample}_sorted.bam")
	conda:
		"envs/align.yaml"
	threads: 8
	shell:
		"samtools sort -@ {threads} -m '4G' {input} > {output}; samtools index -b {output}"

rule dedup_collate:
	input:
	  rules.sort_bams1.output
	output:
		temp("samples/align/{sample}_collate.bam")
	conda:
		"envs/align.yaml"
	log:
		"data/logs/samtools_collate_{sample}.log"
	message: " -- Deduplicating {wildcards.sample} -- "
	shell:
		"(samtools collate -o {output} {input}) > {log}"

rule dedup_fixmate:
	input:
		rules.dedup_collate.output
	output:
		temp("samples/align/{sample}_collate_fixmate.bam")
	conda:
		"envs/align.yaml"
	log:
		"data/logs/samtools_fixmate_{sample}.log"
	shell:
		"(samtools fixmate -m {input} {output}) > {log}"

rule sort_bams2:
	input:
		rules.dedup_fixmate.output
	output:
		temp("samples/align/{sample}_collate_fixmate_sort.bam")
	conda:
		"envs/align.yaml"
	threads: 4
	shell:
		"samtools sort -@ {threads} -m 4G {input} > {output}"

rule dedup_markdup:
	input:
		rules.sort_bams2.output
	output:
		"samples/align/{sample}_dedup.bam" # temp
	conda:
		"envs/align.yaml"
	log:
		"data/logs/samtools_markdup_{sample}.log"
	shell:
		"samtools markdup {input} {output} > {log} 2>&1"

rule remove_mito:
	input:
		rules.dedup_markdup.output
	output:
		temp("samples/align/{sample}_dedup_rmChrM.bam")
	conda: 
		"envs/align.yaml"
	message: " -- Removing mitochondrial reads from {wildcards.sample} -- "
	threads: 8
	shell:
		"samtools view -h {input} | grep -v 'chrM' | samtools sort -@ {threads} -m '4G' - | samtools view -bS > {output}"

rule quality_alignment:
	input:
		rules.remove_mito.output
	output:
		"samples/align/quality_align/{sample}_dedup_rmchrM_quality.bam"
	conda: 
		"envs/align.yaml"
	message: " -- Selecting high-quality mappings from {wildcards.sample} -- "
	threads: 8
	shell:
		"samtools view -@ {threads} -m '4G' -q 30 -bh {input} > {output}; samtools index {output}"

# merge replicates -----------------------------------------------------------------

# {sample} wildcard will decompose into {condition} and {replicate}.

rule merge_bams:
	input:
		expand("samples/align/quality_align/{{cond}}_{rep}_dedup_rmchrM_quality.bam", rep = list(range(1, config["REPLICATES"] + 1)) )
	output:
		bam = "samples/align/merged/{cond}_merged.bam",
		index = "samples/align/merged/{cond}_merged.bam.bai"
	conda:
		"envs/align.yaml"
	log:
		"logs/merge_bams/{cond}.txt"
	threads: 8
	message: " -- Merging bam files into {wildcards.cond}. Double check the merge in logs/merge_bams/{wildcards.cond}.txt -- "
	shell:
		"""
		samtools merge -@ {threads} -f {output} {input};
		sleep 10;
		samtools index {output};
		echo $(date) 'samtools merge -@ {threads} -f {output} {input}' > {log}
		"""
# if uneven replicates, maybe add: expand("samples/align/quality_align/{{cond,MOLMC}}{rep}_dedup_rmchrM_quality.bam", rep = list(range(1, 3)) ) 
# where the {{cond,MOLMC}} uses regex to describe the condition that had less / more than ideal replicates. 
# Then adjust rule all to match the expected output files.

# call peaks + consensus -----------------------------------------------------------
rule bigwig:
	input:
		bam = rules.merge_bams.output.bam,
		index = rules.merge_bams.output.index
	output:
		"samples/align/bigwig/{cond}.bw"
	conda:
		"envs/deeptools.yaml"
	threads: 16
	log:
		"logs/bigwig/{cond}.log"
	message: " -- Making bigwig tracks for {wildcards.cond} -- "
	shell:
		"bamCoverage -b {input.bam} -o {output} -of bigwig -p {threads} 2> {log}"

if config["ORGANISM"] == "hg38":
	genome_size = 2.9e9
elif config["ORGANISM"] == "mm10":
	genome_size = 2.6e9
else: 
	print("Must define organism for peak calling: ['hg38', 'mm10']. Your input = {}".format(config["ORGANISM"]))

rule call_peaks:
	input:
		rules.merge_bams.output
	output:
		"data/macs/{cond}_peaks.broadPeak"
	params:
		genome = genome_size
	conda:
		"envs/macs.yaml"
	log:
		"logs/call_peaks/{cond}.log"
	message: " -- Calling peaks for {wildcards.cond} -- "
	shell:
		"macs2 callpeak --treatment {input} --name {wildcards.cond} \
		--format BAMPE --gsize {params.genome} \
		--outdir data/macs --broad 2> {log}"

rule filter_peaks:
	input:
		rules.call_peaks.output
	output:
		"data/macs/filtered/{cond}_filtered.broadPeak"
	shell:
		"awk '{{if ($9 >= 3) {{print}} }}' {input} > {output}"
# each peak has a q-val, and we want to keep peaks with q-val >= 3.
# for context, peaks with q-val of 1 look indistinguishable from background. 3 is OK. 5 looks peaky.
# this is equivalent to filtering with q-val >= 0.50

rule consensus_intervals:
	input:
		expand("data/macs/filtered/{cond}_filtered.broadPeak", cond = config["CONDITIONS"])
	output:
		"data/macs/consensus_intervals.bed",
	params:
		intersect_count = config["INTERSECT_COUNT"]
	conda:
		"envs/bedtools.yaml"
	log: 
		"logs/consensus_intervals/consensus_intervals.log"
	shell:
		"cat {input} | sort -k1,1 -k2,2n | bedtools merge -c 1 -o count -d 1000 | awk '{{ if ($4 >= {params.intersect_count}) {{print}} }}' > {output} 2> {log}"
# open all files, sort by chr and star pos, merge and count intersected intervals within 1000 bp of each other, and filter for intervals with {INTERSECT_COUNT} number of counts (across conditions).
# we want to merge peaks together because they allow HINT-ATAC to find footprints there. likely will not make data noisier bc rgt-hint cross-ref with bam files.

rule rm_blacklist:
	input:
		regions = rules.consensus_intervals.output,
		bl_file = config["BLACKLIST_FILE"]
	output:
		"data/macs/consensus_intervals_rm_blacklist.bed"
	conda:
		"envs/bedtools.yaml"
	log:
		"logs/rm_blacklist/rm_blacklist.txt"
	shell:
		"bedtools intersect -v -a {input.regions} -b {input.bl_file} > {output}"

# footprint + motif match + tracks ------------------------------------------------

# MAKE SURE YOU CONFIGURED RGT GENOMIC DATA FOR YOUR ORGANISMS BEFORE STARTING THE FOLLOWING SECTION.
rule footprinting:
	input:
		bam = rules.merge_bams.output.bam,
		regions = rules.rm_blacklist.output
	output:
		"data/differential/footprints/{cond}.bed"
	conda:
		"envs/rgt.yaml"
	params: 
		organism = config["ORGANISM"]
	message: " -- Finding footprints for {wildcards.cond} -- "
	threads: 8
	log:
		"logs/footprinting/footprint_{cond}.log"
	shell:
		"rgt-hint footprinting \
		--atac-seq \
		--paired-end \
		--organism={params.organism} \
		--output-location=data/differential/footprints \
		--output-prefix={wildcards.cond} \
		{input.bam} {input.regions} 2> {log}"
# find footprint patterns (high shoulder peaks + notch in the middle) in provided open chromatin regions.

rule diff_motifs:
	input:
		expand("data/differential/footprints/{cond}.bed", cond = list(config["CONDITIONS"]) )
	output:
		expand("data/differential/match/{cond}_mpbs.bed", cond = list(config["CONDITIONS"]) )
	conda:
		"envs/rgt.yaml"
	params: 
		organism = config["ORGANISM"]
	message: " -- Finding motifs in footprints -- "
	log:
		"logs/diff_motifs/diff_motifs.log"
	shell:
		"rgt-motifanalysis matching \
		--organism={params.organism} \
		--output-location=data/differential/match \
		--input-files {input} 2> {log}"
# find motifs in said footprints via web scrape of public motif databases.

rule fix_motifs:
	input:
		"data/differential/match/{cond}_mpbs.bed"
	output:
		"data/differential/match/fixed/{cond}_mpbs.bed"
	shell:
		r"sed 's/[\t]*$//' {input} | sed -E 's/MA[0-9.]+([A-Za-z0-9]+)/\1/' > {output}"
# remove trailing tab from the file. It's an issue for tools like bedtools
# strip the 'MA[numeric].1' prefix from gene names.

rule footprint_tracks:
	input:
		footprint = rules.footprinting.output,
		bam = rules.merge_bams.output.bam
	output:
		"data/differential/tracks/{cond}.wig"
	conda:
		"envs/rgt.yaml"
	params: 
		organism = config["ORGANISM"]
	message: " -- Exporting footprint tracks for {wildcards.cond} -- "
	log:
		"logs/footprint_tracks/tracks_{cond}.log"
	shell:
		"rgt-hint tracks \
		--bc --norm \
		--organism={params.organism} \
		--output-location=data/differential/tracks \
		--output-prefix={wildcards.cond} \
		{input.bam} {input.footprint} 2> {log}"
# export tracks of footprints. 1 = very likely a footprint, 0 = not a footprint.

# differential footprint ------------------------------------------------------------
rule differential:
	input:
		controlBam = "samples/align/merged/{control}_merged.bam",
		controlMotif = "data/differential/match/{control}_mpbs.bed",
		caseBam = "samples/align/merged/{case}_merged.bam",
		caseMotif = "data/differential/match/{case}_mpbs.bed",
	output:
		"data/differential/diff_pseudocounts/{control}_{case}/differential_factor.txt",
		"data/differential/diff_pseudocounts/{control}_{case}/differential_statistics.txt",
		"data/differential/diff_pseudocounts/{control}_{case}/differential_statistics.pdf"
	message: " -- Calculating differential chromatin oppenness between {wildcards.control} and {wildcards.case} -- "
	wildcard_constraints:
		control = "".join(config["CONTROL"])
	params:
		outdir = "data/differential/diff_pseudocounts/{control}_{case}",
		organism = config["ORGANISM"]
	threads: 24
	conda: 
		"envs/rgt.yaml"
	log:
		"logs/differential/{control}_{case}.log"
	shell:
		"rgt-hint differential \
		--organism={params.organism} \
		--bc \
		--nc {threads} \
		--conditions={wildcards.control},{wildcards.case} \
		--output-location={params.outdir} \
		--mpbs-files={input.controlMotif},{input.caseMotif} \
		--reads-files={input.controlBam},{input.caseBam}"

rule sort_stats:
	input:
		"data/differential/diff_pseudocounts/{control}_{case}/differential_statistics.txt"
	output:
		"data/differential/diff_pseudocounts/{control}_{case}/differential_statistics_sorted_filtered.txt"
	log:
		"logs/sort_stats/sort_stats_{control}_{case}.log"
	wildcard_constraints:
		control = "".join(config["CONTROL"])
	shell:
		"awk '{{ if ($2 > 1000 && $9 < 0.05) {{print}} }}' {input} | sort -gk9 > {output} 2> {log}"

# downstream analysis --------------------------------------------------------------------

# see if public ChIP-Atlas resource contains any predicted TF footprints.
rule chip_screen:
	input:
		motifs = rules.fix_motifs.output,
		chip = config["CHIP_PREFIX"] + "{chip}.bed.gz"
	output:
		"data/differential/chipAtlas/{cond}_{chip}.bed"
	conda:
		"envs/bedtools.yaml"
	shell:
		"bedtools intersect -wa -wb -f 0.50 -a {input.motifs} -b {input.chip} | awk '{{print tolower($0)}}' | awk '$4 ~ $10' > {output}"
# -wa -wb = write both beds to the output. 
# -f 0.50 = overlap of interval a to interval b must be 50% or higher.
# lower-case all content. $4 is putative TF motif (gene name), $10 is gene name of chip-seq data. query if putative motif is 'in' chip-seq data.

# snakemake -j 4 --use-conda --rerun-incomplete --latency-wait 60 --keep-going --cluster-config cluster.yaml --cluster "sbatch -p {cluster.partition} -N {cluster.N}  -t {cluster.t} -o {cluster.o} -e {cluster.e} -J {cluster.J} -c {threads} --mem={cluster.mem}" -s Snakefile 

# development notes --------------------------------------------------------------------

# Test dataset for development: Dan's drug combination project with GSK + ORY + Quizartinib.
# for i in $(find /home/groups/MaxsonLab/input-data/ATAC/F20FTSUSAT0390_HUMiinR/Clean -name "*24*.gz" | grep -E "24[D|C|X]" | sort); do ln -s $i ${i:82:11}${i:98:3}R${i:101:100}; done
# ls -1 *R1*.gz | awk -F "_" '{print $2}' | sort | uniq > ID
# for i in `cat ./ID`; do zcat L*_${i}_*_R1.fq.gz | gzip >/${i}_R1.fastq.gz; done
# for i in `cat ./ID`; do zcat L*_${i}_*_R2.fq.gz | gzip >/${i}_R2.fastq.gz; done

# in cda7, i renamed samples with a more defined "_" separator. e.g. MOLM24C_1_R1|R2.fastq.gz
# for i in $(find /home/groups/MaxsonLab/kongg/cda6/samples/raw/ -maxdepth 1 -name "*.gz"); do ln -s $i .; done
# for i in $(ls *.gz); do mv $i ${i:0:7}_${i:7:1}_${i:9:2}.fastq.gz; done