#!/usr/bin/env python
import glob
import os
if not os.path.exists("loglsf"):
    os.makedirs("loglsf")

configfile: "config.yaml"
FASTQDIR = config["resources"]["fastqdir"]
SAMPLES, = glob_wildcards(FASTQDIR + "/{sample}_R1.fastq.gz")


rule all:
    input:
        "qc/multiqc.html",
         expand("aln/{sample}_Aligned.sortedByCoord.out.bam", sample= SAMPLES),
         expand("kallisto/{sample}", sample= SAMPLES)

rule fastp:
    input: 
        R1 = FASTQDIR + "/{sample}_R1.fastq.gz",
        R2 = FASTQDIR + "/{sample}_R2.fastq.gz"
    output:
        R1_qc = "qc/{sample}_R1.fastq.gz",
        R2_qc = "qc/{sample}_R2.fastq.gz",
        JSON = "qc/{sample}_fastp.json"
    params:
        fastp = config["tools"]["FASTP"],
        qualified_quality_phred = 15,
        unqualified_percent_limit = 40,
        nova_seq = "-g"
    shell:
        """
        {params.fastp} -i {input.R1} -o {output.R1_qc} -I {input.R2} -O {output.R2_qc} \
                -j {output.JSON} -q {params.qualified_quality_phred} \
                -u {params.unqualified_percent_limit} \
                {params.nova_seq} >/dev/null
        """

rule multiqc:
    input:
        lambda wildcards: expand("qc/{sample}_fastp.json", sample=SAMPLES) 
    output:
        "qc/multiqc.html"
    envmodules:
        "gcc/4.8.5",
        "multiqc/1.8"
    shell:
        "multiqc `pwd`/qc"


rule alignment:
    input:
        R1= "qc/{sample}_R1.fastq.gz",
        R2= "qc/{sample}_R2.fastq.gz"
    output:
        "aln/{sample}_Aligned.sortedByCoord.out.bam"
    params:
        genomedir=config["resources"]["genomedir"]
        prefix="aln/{sample}_"
    threads:
        8
    envmodules:
        "gcc/4.8.5",
        "star/2.7.3a",
        "samtools/1.8"
    shell:
        """
        STAR --genomeDir {params.genomedir} --readFilesIn {input.R1} {input.R2} \
                --readFilesCommand zcat \
                --outSAMtype BAM SortedByCoordinate \
                --quantMode GeneCounts \
                --outFileNamePrefix {params.prefix} \
                --runThreadN {threads} && samtools index -@ {threads} {output}
        """

rule kallisto:
    input:
        R1= "qc/{sample}_R1.fastq.gz",
        R2= "qc/{sample}_R2.fastq.gz"
    output:
        directory("kallisto/{sample}")
    params:
        cDNAidx=config["resources"]["cDNAidx"]
    threads:
        6
    envmodules:
        "gcc/4.8.2",
        "kallisto/0.43.0"
    shell:
        """
        kallisto quant --fusion -i ${params.cDNAidx} -o {output} {input.R1} ${input.R2} \
                -t {threads} -b 100 > /dev/null 2>&1"
        """
