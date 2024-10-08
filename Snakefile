## Load config values
configfile: "config/samples.yaml"
configfile: "config/parameters.yaml"

# path to the reference genome fasta
genome_path = config["genome_path"]
transcriptome_path = config["transcriptome_path"]

# samples
samples_dict = config["samples"]
sample_ids = samples_dict.keys()

### -------------------------------------------------------------------
### Target rule
### -------------------------------------------------------------------
rule all:
	input:
		expand("output/{sample}/star/nmappedreads.txt", sample=sample_ids),
		expand("output/{sample}/expressed_hpv_integration_sites.bed", sample=sample_ids),
        expand("output/{sample}/star/hpvRNAAligned.out.sorted.bam", sample=sample_ids)

### -------------------------------------------------------------------
### Rules
### -------------------------------------------------------------------

rule make_fastq:
    input:
        rna = lambda w: config["samples"][w.sample]["rna"]
    output:
        r1 = "output/{sample}/fastq/R1.fastq",
        r2 = "output/{sample}/fastq/R2.fastq"
    conda: "config/conda.yaml"
    threads: 12
    shell:
        """
        samtools fastq -n {input.rna} -1 {output.r1} -2 {output.r2}
        """

rule star_align:
    input:
        r1 = "output/{sample}/fastq/R1.fastq",
        r2 = "output/{sample}/fastq/R2.fastq",
        ref = transcriptome_path
    output:
        "output/{sample}/star/hpvRNAAligned.out.bam",
        "output/{sample}/star/hpvRNAChimeric.out.junction",
        "output/{sample}/star/hpvRNALog.final.out"
    conda: "config/conda.yaml"
    threads: 12
    shell:
        """
        STAR --runThreadN {threads} --genomeDir {input.ref} --readFilesIn {input.r1} {input.r2} \
        --outSAMtype BAM Unsorted --outFileNamePrefix output/{wildcards.sample}/star/hpvRNA \
        --outSAMstrandField intronMotif --chimSegmentMin 12 --chimJunctionOverhangMin 8 \
        --chimOutJunctionFormat 1 --alignSJDBoverhangMin 10 --alignMatesGapMax 100000 \
        --alignIntronMax 100000 --limitSjdbInsertNsj 1500000 --alignSJstitchMismatchNmax 5 -1 5 5 \
        --outSAMattrRGline ID:GRPundef --chimMultimapScoreRange 3 --chimScoreJunctionNonGTAG 4 \
        --chimMultimapNmax 20 --chimNonchimScoreDropMin 10 --peOverlapNbasesMin 12 --peOverlapMMp 0.1 \
        --alignInsertionFlush Right --alignSplicedMateMapLminOverLmate 0 --alignSplicedMateMapLmin 30 \
        --chimOutType Junctions WithinBAM SoftClip 
        """

rule sort:
    input:
        bam = "output/{sample}/star/hpvRNAAligned.out.bam"
    output:
        "output/{sample}/star/hpvRNAAligned.out.sorted.bam"
    conda: "config/conda.yaml"
    threads: 12
    shell:
        """
        sambamba sort {input.bam} -t {threads}
        """

rule num_mapped_reads:
    input:
        "output/{sample}/star/hpvRNALog.final.out"
    output:
        "output/{sample}/star/nmappedreads.txt"
    conda: "config/conda.yaml"
    threads: 12
    shell:
        """
        grep "Uniquely mapped reads number" {input} | awk '{{print $NF}}' > {output}
        """

rule filter_juntion:
    input:
        junction = "output/{sample}/star/hpvRNAChimeric.out.junction"
    output:
        "output/{sample}/star/hpvRNAChimeric.out.filtered.junction"
    conda: "config/conda.yaml"
    shell:
        """
        grep 'HPV' {input.junction} > {output}
        """

rule find_RNA_breakpoints:
    input:
        junction = "output/{sample}/star/hpvRNAChimeric.out.filtered.junction"
    output:
        "output/{sample}/hpvRNAbreakpoints.bed"
    conda: "config/conda.yaml"
    shell:
        """
        scripts/findHPVRNABreakpoints.py {input.junction} output/{wildcards.sample} > {output}
        """

rule sort_bed:
    input:
        "output/{sample}/hpvRNAbreakpoints.bed"
    output:
        "output/{sample}/hpvRNAbreakpoints.sorted.bed"
    conda: "config/conda.yaml"
    shell:
        """
        bedtools sort -i {input} > {output}
        """

rule sort_bed2:
    input:
        hpv = lambda w: config["samples"][w.sample]["ont-dna-sites"]
    output:
        "output/{sample}/hpv_integration_sites.sorted.bed"
    conda: "config/conda.yaml"
    shell:
        """
        bedtools sort -i {input.hpv} > {output}
        """

rule bedtools_dist:
    input:
        a="output/{sample}/hpvRNAbreakpoints.sorted.bed",
        b="output/{sample}/hpv_integration_sites.sorted.bed"
    output:
        "output/{sample}/expressed_hpv_integration_sites.bed"
    conda: "config/conda.yaml"
    shell:
        """
        bedtools closest -D "a" -a {input.a} -b {input.b} | awk '$8 !~ /\./ && sqrt($12*$12) < 100000' > {output}
        """
