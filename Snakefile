configfile: "config.yaml"

rule all:
    input:
        expand("outputs/multiqc_report_raw/{filename}{num}_001_multiqc_report.html", filename=[filename for method in config["method"].keys() for filename in config["method"][method]["files"]], num=config["num"]),
        expand("outputs/trimmed_reads/{filename}2_001_trimmed.fastq", filename=[filename for method in config["method"].keys() for filename in config["method"][method]["files"]]),
        expand("outputs/multiqc_report_trimmed/{filename}{num}_001_trimmed_multiqc_report.html", filename=[filename for method in config["method"].keys() for filename in config["method"][method]["files"]], num=config["num"]),
        expand("outputs/STAR/{filename}/Aligned.sortedByCoord.out.bam.bai", filename=[filename for method in config["method"].keys() for filename in config["method"][method]["files"]]),
        expand("outputs/STAR/{filename}/counts_2.txt", filename=[filename for method in config["method"].keys() for filename in config["method"][method]["files"]]),
        [f"outputs/STAR/all/{x}/vulcano-plot.eps" for x in config["method"]],
        "outputs/STAR/all/pca-plot.eps"
        
rule raw_fastqc:
    input:
        "data/samples/{filename}{num}_001.fastq"
    output:
        "outputs/fastqc_raw/{filename}{num}_001_fastqc.html",
        "outputs/fastqc_raw/{filename}{num}_001_fastqc.zip"
    shell:
        "fastqc {input} -o outputs/fastqc_raw/"

rule raw_multiqc:
    input:
        "outputs/fastqc_raw/{filename}{num}_001_fastqc.html",
        "outputs/fastqc_raw/{filename}{num}_001_fastqc.zip"
    output:
        "outputs/multiqc_report_raw/{filename}{num}_001_multiqc_report.html"
    shell:
        "multiqc ./outputs/fastqc_raw/ -n {output}"

rule trim_bbduk:
    input:
        read1="data/samples/{filename}1_001.fastq",
        read2="data/samples/{filename}2_001.fastq",
        adapters="data/adapters.fa"
    output:
        out1="outputs/trimmed_reads/{filename}1_001_trimmed.fastq",
        out2="outputs/trimmed_reads/{filename}2_001_trimmed.fastq"
    shell:
        "./bin/bbmap/bbduk.sh -Xmx2g in1={input.read1} in2={input.read2} out1={output.out1} out2={output.out2} ref={input.adapters} ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10"

        
rule trimmed_fastqc:
    input:
        "outputs/trimmed_reads/{filename}{num}_001_trimmed.fastq"
    output:
        "outputs/fastqc_trimmed/{filename}{num}_001_trimmed_fastqc.html",
        "outputs/fastqc_trimmed/{filename}{num}_001_trimmed_fastqc.zip"
    shell:
        "fastqc {input} -o outputs/fastqc_trimmed/"

rule trimmed_multiqc:
    input:
        "outputs/fastqc_trimmed/{filename}{num}_001_trimmed_fastqc.html",
        "outputs/fastqc_trimmed/{filename}{num}_001_trimmed_fastqc.zip"
    output:
        "outputs/multiqc_report_trimmed/{filename}{num}_001_trimmed_multiqc_report.html"
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
        read1="outputs/trimmed_reads/{filename}1_001_trimmed.fastq",
        read2="outputs/trimmed_reads/{filename}2_001_trimmed.fastq"
    output:
        "outputs/STAR/{filename}/Aligned.sortedByCoord.out.bam"
    shell:
        "mkdir -p outputs/STAR/{wildcards.filename} && STAR --runThreadN 1 --genomeDir {input.gendir} --readFilesIn {input.read1} {input.read2} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix outputs/STAR/{wildcards.filename}/"
        
rule bam_index:
    input:
        "outputs/STAR/{filename}/Aligned.sortedByCoord.out.bam"
    output:
        "outputs/STAR/{filename}/Aligned.sortedByCoord.out.bam.bai"
    shell:
        "samtools index {input}"
        
rule bam_sort_name:
    input:
        "outputs/STAR/{filename}/Aligned.sortedByCoord.out.bam"
    output:
        "outputs/STAR/{filename}/Aligned.sortedByCoord.out.sortedbyname.bam"
    shell:
        "samtools sort -n -o {output} {input}"

rule feature_counts:
    input:
        bam="outputs/STAR/{filename}/Aligned.sortedByCoord.out.sortedbyname.bam",
        gtf="data/chr19_20Mb.gtf"
    output:
        out1="outputs/STAR/{filename}/counts_1.txt",
        out2="outputs/STAR/{filename}/counts_2.txt"
    shell:
        "featureCounts -p -t exon -g gene_id -a {input.gtf} -o {output.out1} -s 1 {input.bam} && featureCounts -p -t exon -g gene_id -a {input.gtf} -o {output.out2} -s 2 {input.bam}"


rule feature_counts_per_sample:
    input:
        bam=lambda wildcards: [f"outputs/STAR/{name}/Aligned.sortedByCoord.out.sortedbyname.bam" for name in config["method"][wildcards.x]["files"]],
        gtf="data/chr19_20Mb.gtf"
    params:
        s=lambda wildcards: [x for x in config["method"][wildcards.x]["s"]],
    output:
        outA="outputs/STAR/all/{x}/counts.txt"
    shell:
        "mkdir -p outputs/STAR/all/{wildcards.x}/ && featureCounts -p -t exon -g gene_id -a {input.gtf} -o {output.outA} -s {params.s} {input.bam}"

rule de_analysis:
    input:
        counts = "outputs/STAR/all/{x}/counts.txt"
    output:
        plot="outputs/STAR/all/{x}/vulcano-plot.eps",
        pval="outputs/STAR/all/{x}/DE.csv",
        cts="outputs/STAR/all/{x}/cts.csv",
        coldata="outputs/STAR/all/{x}/coldata.csv"
    shell:
        "./bin/DE.R -i {input.counts} -o {output.pval} -v {output.plot}"

rule pca_plot:
    input:
        expand("outputs/STAR/all/{x}/cts.csv", x=config["method"])
    output:
        "outputs/STAR/all/pca-plot.eps"
    shell:
        "./bin/PCA.R {input} {output}"
