def downloadSRA(SRA_accession):
    SAMPLE=[]
    with open(SRA_accession) as f:
        for line in f:
            line = line.rstrip('\n')
            SAMPLE.append(line)

    return SAMPLE

def write_URL_files(SRA_accession):
    website_for ='ftp://ftp-trace.ncbi.nlm.nih.gov/sra' + \
        '/sra-instant/reads/ByStudy/sra/SRP/SRP014/SRP014652/'
    with open(SRA_accession, "r") as f:
        for line in f:
            l = line.strip()
            download_site = website_for + l +"/" +l+ ".sra"
            with open(l + "temp.txt", "w") as fw:
                fw.write(download_site)



write_URL_files("SRR_accession.txt")
SAMPLES = downloadSRA("SRR_accession.txt")


rule all:
    input:
        expand('./trimmed_FASTQ/FASTQC/{sample}.trimmed_fastqc.html', sample=SAMPLES),
        expand('./trimmed_FASTQ/FASTQC/{sample}.trimmed_fastqc.zip', sample=SAMPLES),
        expand('./trimmed_FASTQ/{sample}.bam', sample=SAMPLES),
        expand('./trimmed_FASTQ/hisat2_summary_{sample}.txt', sample=SAMPLES),
        "./featureCount_result/Round8.featureCount",
        "./featureCount_result/FC_summary_Round8.txt",
        expand('./kallisto_result/{sample}',sample=SAMPLES),


rule wget_SRA:
    input:
        "./{sample}temp.txt"
    output:
        temp("./{sample}.sra")
    shell:"""
        wget -i {input}
        """


rule fastq_dump_SRA:
    input: 
        "./{sample}.sra"
    output:
        temp("./{sample}.fastq.gz")
    shell:"""
        fastq-dump --gzip {input}
        """


rule CUTADAPT_SRA:
    input:
        "./{sample}.fastq.gz",
    output:
        temp("./trimmed_FASTQ/{sample}.trimmed.fastq.gz")
    threads: 2
    shell:
       'cutadapt -a AGATCGGAAGAGC -m 80 -q 20 -o {output} {input}'
       #'cutadapt -a AGATCGGAAGAGC -m 80 -q 20 -o {output} {input} 1>info.txt'
       
rule FASTQC_SRA:
    input:
        expand("./trimmed_FASTQ/{sample}.trimmed.fastq.gz",sample=SAMPLES)
    output:
        expand('./trimmed_FASTQ/FASTQC/{sample}.trimmed_fastqc.html', sample=SAMPLES),
        expand('./trimmed_FASTQ/FASTQC/{sample}.trimmed_fastqc.zip', sample=SAMPLES),

    shell:"""
        fastqc -t 8 {input} -o ./trimmed_FASTQ/FASTQC 
        """

# This is for single-end files.

rule KALLISTO_ALN:
    input:
        "./trimmed_FASTQ/{sample}.trimmed.fastq.gz"
    output:
        "./kallisto_result/{sample}"
    shell:
        " kallisto quant -i ~/GenomeFiles/AGPv3/ZmaAGPv3.22_cdna_kallisto.idx -t 2 -o {output} -b 20 --single -l 180 -s 20 {input}"
        

rule HISAT_ALN:
    input:
        "./trimmed_FASTQ/{sample}.trimmed.fastq.gz"
    output:
        "./trimmed_FASTQ/{sample}.bam",
        "./trimmed_FASTQ/hisat2_summary_{sample}.txt"
    shell:
        "hisat2 -x ~/GenomeFiles/AGPv3/Sequence/hisat2_zm3_index/genome_snp_tran_ercc -p 4 -t -U {input} 2>{output[1]} | samtools view -Sbo {output[0]} -  "

rule FEATURE_COUNT:
    input:
        expand("./trimmed_FASTQ/{sample}.bam",sample=SAMPLES)
    output:
        "./featureCount_result/Round8.featureCount",
        "./featureCount_result/FC_summary_Round8.txt"
    shell:
        "featureCounts -a ~/GenomeFiles/Ensembl/Zmav3.gtf -o {output[0]} {input} 2> {output[1]}"
     
# cufflink process is too slow, going to discard it. But can try kallisto since it's ultra fast.


