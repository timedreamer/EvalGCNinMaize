def downloadSRA(SRA_accession):
    SAMPLE=[]
    with open(SRA_accession) as f:
        for line in f:
            line = line.rstrip('\n')
            SAMPLE.append(line)

    return SAMPLE

# def write_URL_files(SRA_accession):
    # website_for ='ftp://ftp-trace.ncbi.nlm.nih.gov/sra' + \
        # '/sra-instant/reads/ByStudy/sra/SRP/SRP072/SRP072207/'
    # with open(SRA_accession, "r") as f:
        # for line in f:
            # l = line.strip()
            # download_site = website_for + l +"/" +l+ ".sra"
            # with open(l + "temp.txt", "w") as fw:
                # fw.write(download_site)

def write_URL_files(SRA_accession):
    website_for ='ftp://ftp-trace.ncbi.nlm.nih.gov/sra' +\
                  '/sra-instant/reads/ByRun/sra/'
    with open(SRA_accession,"r") as f:
        for line in f:
            l = line.strip()
            threeChar=l[0:3]
            sixChar = l[0:6]
            download_site = website_for + threeChar + "/" + sixChar + "/"\
                            + l + "/" + l + ".sra"
            with open(l + "temp.txt","w") as fw:
                fw.write(download_site)

write_URL_files("SRR_accession.txt")
SAMPLES = downloadSRA("SRR_accession.txt")



rule all:
    input:
        expand('./trimmed_FASTQ/FASTQC/{sample}_1.trimmed_fastqc.html', sample=SAMPLES),
        expand('./trimmed_FASTQ/FASTQC/{sample}_1.trimmed_fastqc.zip', sample=SAMPLES),
        expand('./trimmed_FASTQ/FASTQC/{sample}_2.trimmed_fastqc.html', sample=SAMPLES),
        expand('./trimmed_FASTQ/FASTQC/{sample}_2.trimmed_fastqc.zip', sample=SAMPLES),
        expand('./trimmed_FASTQ/{sample}_PE.bam', sample=SAMPLES),
        expand('./trimmed_FASTQ/hisat2_summary_{sample}_PE.txt', sample=SAMPLES),
        "./featureCount_result/Round14.featureCount",
        "./featureCount_result/FC_summary_Round14.txt",
        expand('./kallisto_result/{sample}_PE',sample=SAMPLES),


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
        temp("./{sample}_1.fastq.gz"),
        temp("./{sample}_2.fastq.gz"),
    shell:"""
        fastq-dump --split-files --gzip {input}
        """


rule CUTADAPT_SRA:
    input:
        "./{sample}_1.fastq.gz",
        "./{sample}_2.fastq.gz"
    output:
        temp("./trimmed_FASTQ/{sample}_1.trimmed.fastq.gz"),
        temp("./trimmed_FASTQ/{sample}_2.trimmed.fastq.gz"),
    threads: 4
    shell:
       'cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 80 -q 20 -o {output[0]} -p {output[1]} {input[0]} {input[1]}'
       
rule FASTQC_SRA:
    input:
        expand("./trimmed_FASTQ/{sample}_1.trimmed.fastq.gz",sample=SAMPLES),
        expand("./trimmed_FASTQ/{sample}_2.trimmed.fastq.gz",sample=SAMPLES),
    output:
        expand('./trimmed_FASTQ/FASTQC/{sample}_1.trimmed_fastqc.html', sample=SAMPLES),
        expand('./trimmed_FASTQ/FASTQC/{sample}_1.trimmed_fastqc.zip', sample=SAMPLES),
        expand('./trimmed_FASTQ/FASTQC/{sample}_2.trimmed_fastqc.html', sample=SAMPLES),
        expand('./trimmed_FASTQ/FASTQC/{sample}_2.trimmed_fastqc.zip', sample=SAMPLES)
    shell:
        "fastqc -t 8 {input} -o ./trimmed_FASTQ/FASTQC "

        
#This is for paired-end
rule KALLISTO_ALN:
    input:
        "./trimmed_FASTQ/{sample}_1.trimmed.fastq.gz",
        "./trimmed_FASTQ/{sample}_2.trimmed.fastq.gz",
    output:
        "./kallisto_result/{sample}_PE"
    shell:
        " kallisto quant -i ~/GenomeFiles/AGPv3/ZmaAGPv3.22_cdna_kallisto.idx -t 2 -o {output} -b 20 {input[0]} {input[1]}"
        
        
rule HISAT_ALN:
    input:
        "./trimmed_FASTQ/{sample}_1.trimmed.fastq.gz",
        "./trimmed_FASTQ/{sample}_2.trimmed.fastq.gz",
    output:
        "./trimmed_FASTQ/{sample}_PE.bam",
        "./trimmed_FASTQ/hisat2_summary_{sample}_PE.txt"
        
    shell:
        "hisat2 -x ~/GenomeFiles/AGPv3/Sequence/hisat2_zm3_index/genome_snp_tran_ercc -p 4 -t -1 {input[0]} -2 {input[1]} 2>{output[1]} |samtools view -Sbo {output[0]} - "

rule FEATURE_COUNT:
    input:
        expand("./trimmed_FASTQ/{sample}_PE.bam",sample=SAMPLES)
    output:
        "./featureCount_result/Round14.featureCount",
        "./featureCount_result/FC_summary_Round14.txt"
    shell:
        "featureCounts -p -a ~/GenomeFiles/Ensembl/Zmav3.gtf -o {output[0]} {input} 2> {output[1]}"
     
        
