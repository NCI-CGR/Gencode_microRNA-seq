### This pipeline is for general QC and alignment of miRNA-seq data

## vim: ft=python
import sys
import os
import glob
import itertools

shell.prefix("set -eo pipefail; ")
configfile: "config.json"
os.makedirs("log", exist_ok=True)
# localrules: all

# define wildcards
# def parse_sampleID(fname):
#     base=os.path.basename(fname)
#     return fname.split('/')[-1].split('_')[0]

def parse_sampleID(fname):
    base=os.path.basename(fname)
    sample = base.rsplit('.fastq.gz')[0].rsplit('.fq.gz')[0]
    return sample

# file = sorted(glob.glob(os.path.join(config['fastq_dir'],'*.fastq.gz')), key=parse_sampleID)
file = sorted(glob.glob('merged_fastq/*.fastq.gz'), key=parse_sampleID)
sample_files_dict = {}
for key, value in itertools.groupby(file, parse_sampleID):
    sample_files_dict[key] = list(value)
       
rule all:
    input:
          expand("star_align/{sample}/{sample}Aligned.sortedByCoord.out.bam",sample=sample_files_dict.keys()),
          "pretrim_qc/preQC_multiqc_report.html",
          "posttrim_qc/postQC_multiqc_report.html",
          "star_align/log/star_align_multiqc_report.html",
          "reads_count/reads_count.csv"      


rule cutadapt:
    input: 
          fastq="merged_fastq/{sample}.fastq.gz",
          adapters=config["adapters"]

    output:
          trimmed="trimmed/{sample}.trim.fastq.gz",
          too_short="trimmed/{sample}_too_short.fastq.gz",
          json_report="trimmed/{sample}_report.json",
          log_report="trimmed/{sample}.log"
    
    singularity:
        "docker://quay.io/biocontainers/cutadapt:4.4--py39hf95cd2a_1"
    threads: 10
    shell:
          """
          cutadapt \
            --adapter file:{input.adapters} \
            --minimum-length 15 \
            --overlap 5 \
            --too-short-output={output.too_short} \
            --quality-cutoff 10,10 \
            -o {output.trimmed} \
            --cores {threads} \
            --json {output.json_report} \
            {input.fastq} \
            > {output.log_report} \
            2>log/{wildcards.sample}_cutadapt.err
          """
      
rule pretrim_qc:
    input: 
          "merged_fastq/{sample}.fastq.gz" 
    output: 
          "pretrim_qc/{sample}_fastqc.zip",
          "pretrim_qc/{sample}_fastqc.html"
    threads: 10
    singularity:
            "docker://cgrlab/qctools:v2.0"
    shell:
          """
          fastqc {input} -o pretrim_qc -f fastq --noextract 2>log/{wildcards.sample}_preqc.err
          """


rule posttrim_qc:
    input: 
          "trimmed/{sample}.trim.fastq.gz"
    output: 
          "posttrim_qc/{sample}.trim_fastqc.zip",
          "posttrim_qc/{sample}.trim_fastqc.html"
    threads: 10
    singularity:
            "docker://cgrlab/qctools:v2.0"
    shell:
          """
          fastqc {input} -o posttrim_qc -f fastq --noextract 2>log/{wildcards.sample}_postqc.err
          """



rule star_index:
    input:
            fasta=config['genome_fasta'],
            annot=config['genome_annotation']
    output:
            "star_index/complete.txt"
    threads: 24
    resources:
            mem_mb=500000
    singularity:
            "docker://quay.io/biocontainers/star:2.7.10b--h6b7c446_1"
    shell:
            """
            MAX_RAM=$(({resources.mem_mb} * 1000000))
            STAR --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir star_index \
            --sjdbGTFfile {input.annot} \
            --sjdbOverhang 1 \
            --limitGenomeGenerateRAM $MAX_RAM \
            --genomeFastaFiles {input.fasta} 2>log/star_index.err 
            touch {output}
            """

rule star_align:
    input:
          mirna_annot=config["miRNA_Annotation"],
          trimmed_fastq="trimmed/{sample}.trim.fastq.gz",
          star_index_complete="star_index/complete.txt"
    output:
          "star_align/{sample}/{sample}Aligned.sortedByCoord.out.bam",
          "star_align/{sample}/{sample}Log.final.out",
          "star_align/{sample}/{sample}ReadsPerGene.out.tab"
    threads: 24
    singularity:
            "docker://quay.io/biocontainers/star:2.7.10b--h6b7c446_1"
    params:
          index="star_index",
          mirna_annot = config["miRNA_Annotation"]
    shell:
          """
          STAR --runThreadN {threads} \
          --genomeDir {params.index} \
          --readFilesIn {input.trimmed_fastq} \
          --outFileNamePrefix star_align/{wildcards.sample}/{wildcards.sample} \
          --readFilesCommand zcat \
          --sjdbGTFfile {input.mirna_annot} \
          --alignEndsType EndToEnd \
          --outFilterMismatchNmax 1 \
          --outFilterMultimapScoreRange 0 \
          --quantMode TranscriptomeSAM GeneCounts \
          --outReadsUnmapped Fastx \
          --outSAMtype BAM SortedByCoordinate \
          --outFilterMultimapNmax 10 \
          --outSAMunmapped Within \
          --outFilterScoreMinOverLread 0 \
          --outFilterMatchNminOverLread 0 \
          --outFilterMatchNmin 16 \
          --alignSJDBoverhangMin 1000 \
          --alignIntronMax 1 \
          --outWigType wiggle \
          --outWigStrand Stranded \
          --outWigNorm RPM 2>log/{wildcards.sample}_star_align.err 
          """

rule multiqc:
    input:
          expand("pretrim_qc/{sample}_fastqc.html",sample=sample_files_dict.keys()),
          expand("posttrim_qc/{sample}.trim_fastqc.html",sample=sample_files_dict.keys()),
          expand("star_align/{sample}/{sample}Log.final.out",sample=sample_files_dict.keys())
          expand("trimmed/{sample}.log", sample=sample_files_dict.keys())
    output:
          "pretrim_qc/preQC_multiqc_report.html",
          "posttrim_qc/postQC_multiqc_report.html",
          "star_align/log/star_align_multiqc_report.html"
    threads: 8
    singularity:
            "docker://cgrlab/qctools:v2.0"
    shell:
          """
          multiqc pretrim_qc/. --title preQC -o pretrim_qc 2>log/multiqc_preqc.err
          multiqc trimmed/. --title Cutadapt -o trimmed 2>log/multiqc_cutadapt.err
          multiqc posttrim_qc/. --title postQC -o posttrim_qc/ 2>log/multiqc_postqc.err
          mkdir -p star_align/log
          cp star_align/*/*Log.final.out star_align/log
          multiqc star_align/log/. --title star_align -o star_align/log 2>log/multiqc_star.err
          """

## library is sense stranded
rule merge:
    input:  
          expand("star_align/{sample}/{sample}ReadsPerGene.out.tab",sample=sample_files_dict.keys())
    output:
          "reads_count/reads_count.csv"
    singularity:
      "docker://cgrlab/gencode-microrna-seq:v1.0"
    shell:
          """
          Rscript merge.R 2>log/merge_count.err
          """

