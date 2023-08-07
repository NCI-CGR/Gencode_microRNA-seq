# Gencode_microRNA-seq
microRNA-seq workflow utilizing Cutadapt and STAR to generate a Sample-Gene read count matrix

#### Authors: Ben Jordan (ben.jordan@nih.gov) and Komal Jain (komal.jain@nih.gov)

## I. Description
Major steps in the workflow include:
1) Trimming of adapters and low-quality reads using Cutadapt and retention of reads with a minimum length of 15 nt
2) Generating QC reports using FastQC and aggregating results using MultiQC
3) Aligning trimmed reads to GRCh38 human reference genome (https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz)) using STAR, and quantifying microRNA reads according to GENCODE V29 (ENCFF470CZH.gtf, https://www.encodeproject.org/files/ENCFF470CZH/) genome annotation file downloaded from Gencode
Using the trimmed reads from Step 1, mapping to the miRBase gtf by extracting the hairpin/stem-loop sequences from hsa.gff3 file downloaded from here (https://www.mirbase.org/ftp.shtml).  The gff3 file was converted to gtf using gffread software, then a custom script was used to bring the gtf file to the format expected by STAR mapper
5) Merging reads-count tables of all samples for both miRBase and Encode gtf resulting in two count matrices

![DAG](dag.jpeg)
## II. Dependencies
* [Snakemake](https://snakemake.readthedocs.io/en/stable/)

**Docker Containers**
* [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [MultiQC](https://multiqc.info)
* [STAR](https://github.com/alexdobin/STAR)
* [R](https://www.r-project.org)
## III. Input
* merged fastq files stored in directory: `merged_fastq`
* `adapters.fa` with adapter sequences
* reference genome sequence and annotation files
* microRNA annotation file
* `config.json` with updated file locations
## IV. Output
* Trimmed reads in directory: `trimmed/`
* QC reports of pre-trimmed reads in direcotry: `pretrim_qc/`
* QC reports of post-trimmed reads in direcotry: `posttrim_qc/`
* STAR index in directory: `star_index/`
* STAR alignent results and statistics reports in directory : `star_align/`
* merged reads-count table: `reads_count/reads_count.csv`

## V. Running the workflow
1. Update `config.json` file with the following:
    - reference fasta
    - reference annotation
    - adapters fasta file (optional)
    - miRNA annotations (optional)

2. Run singularity

Example: running on CCAD2 cluster
```
module load slurm
module load singularity
module load conda

conda activate snakemake

mkdir -p singularity_cache

snakemake \
    --use-singularity \
    --singularity-prefix singularity_cache \
    --keep-going \
    --local-cores $SLURM_CPUS_PER_TASK \
    --jobs 10 \
    --slurm \
    --max-jobs-per-second 1 \
    --max-status-checks-per-second 0.01 \
    --latency-wait 120 all
```
## VI. Working directory structure
```bash
.
├── adapters.fa
├── log
│   └── log files
├── merged_fastq
│   └── {sample}.fastq.gz
├── merge.R
├── posttrim_qc
│   ├── posttrim_qc_multiqc_report.html
│   ├── {sample}.trim_fastqc.html
│   └── {sample}.trim_fastqc.zip
├── pretrim_qc
│   ├── pretrim_qc_multiqc_report.html
│   ├── {sample}_fastqc.html
│   └── {sample}_fastqc.zip
├── reads_count
│   └── reads_count.csv
├── sample_names.txt
├── Snakefile
├── star_align
│   ├── log
│   │   ├── {sample}Log.final.out
│   │   └── star_align_multiqc_report.html
│   └── {sample}
│       ├── {sample}Aligned.sortedByCoord.out.bam
│       ├── {sample}Log.final.out
│       ├── {sample}ReadsPerGene.out.tab
│       └── other star output files
├── star_index
│   ├── ENCFF628BVT.gtf
│   ├── gencode.v24.primary_assembly.annotation-tRNAs-ERCC_phiX.gtf
│   └── other index files
└── trimmed
    ├── {sample}_too_long.fastq.gz
    ├── {sample}_too_short.fastq.gz
    ├── {sample}.trim.fastq.gz
    └── {sample}.trim.log
```
