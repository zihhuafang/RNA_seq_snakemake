# RNA_seq_snakemake

This is the snakemake pipeline for transcriptome profiling. 
It consists three steps:
1. QC rna sequencing reads
2. Align reads passed QC using STAR 
3. Quantify the abundance of transcripts using kallisto
