configfile: "config.yaml"

rule all:
    input:
        expand("outputs/multiqc_report_trimmed/{sample}_trimmed_multiqc_report.html", sample=config["samples"]),
        expand("outputs/multiqc_report_raw/{sample}_multiqc_report.html", sample=config["samples"]),
        expand("outputs/STAR/{samplename}/Aligned.sortedByCoord.out.bam", samplename=config["samplename"])
        
rule raw_fastqc:
    input:
        "data/samples/{sample}.fastq"
    output:
        "outputs/fastqc_raw/{sample}_fastqc.html",
        "outputs/fastqc_raw/{sample}_fastqc.zip"
    shell:
        "fastqc {input} -o outputs/fastqc_raw/"

rule raw_multiqc:
    input:
        "outputs/fastqc_raw/{sample}_fastqc.html",
        "outputs/fastqc_raw/{sample}_fastqc.zip"
    output:
        "outputs/multiqc_report_raw/{sample}_multiqc_report.html"
    shell:
        "multiqc ./outputs/fastqc_raw/ -n {output}"

rule trim_bbduk:
    input:
        read="data/samples/{sample}.fastq",
        adapters="data/adapters.fa"
    output:
        out="outputs/trimmed_reads/{sample}_trimmed.fastq"
    shell:
        "./bin/bbmap/bbduk.sh in={input.read} out={output.out} ref={input.adapters} ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10"

        
rule trimmed_fastqc:
    input:
        "outputs/trimmed_reads/{sample}_trimmed.fastq"
    output:
        "outputs/fastqc_trimmed/{sample}_trimmed_fastqc.html",
        "outputs/fastqc_trimmed/{sample}_trimmed_fastqc.zip"
    shell:
        "fastqc {input} -o outputs/fastqc_trimmed/"

rule trimmed_multiqc:
    input:
        "outputs/fastqc_trimmed/{sample}_trimmed_fastqc.html",
        "outputs/fastqc_trimmed/{sample}_trimmed_fastqc.zip"
    output:
        "outputs/multiqc_report_trimmed/{sample}_trimmed_multiqc_report.html"
    shell:
        "multiqc ./outputs/fastqc_trimmed/ -n {output}"
        
rule star_indices:
    input:
        fa="data/chr19_20Mb.fa",
        gtf="data/chr19_20Mb.gtf"
    output:
        "outputs/STAR/chr19_20Mb/Genome",
        directory("outputs/STAR/chr19_20Mb/")
    shell:
        "mkdir -p ./outputs/STAR/chr19_20Mb && STAR --runThreadN 1 --runMode genomeGenerate --genomeDir outputs/STAR/chr19_20Mb/ --genomeFastaFiles {input.fa} --sjdbGTFfile {input.gtf} --genomeSAindexNbases 11"

rule star_map:
    input:
        gendir="outputs/STAR/chr19_20Mb/",
        read1="outputs/trimmed_reads/{samplename}1_001_trimmed.fastq",
        read2="outputs/trimmed_reads/{samplename}2_001_trimmed.fastq"
    output:
        "outputs/STAR/{samplename}/Aligned.sortedByCoord.out.bam"
    shell:
        "mkdir -p outputs/STAR/{wildcards.samplename} && STAR --runThreadN 1 --genomeDir {input.gendir} --readFilesIn {input.read1} {input.read2} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix outputs/STAR/{wildcards.samplename}/"
