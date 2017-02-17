*Steps to process RNAseq

Assumptions:
 
--Output_directory

----Fasta_files (original files)

----Bam_sam (sorted and indexed)

----Expression_Values

#Steps to run expression analysis

1. FASTQC quality check of raw reads
    1. FASTQC will provide information including the quality of reads along the length of a read, the rate of duplication within the reads, and biases of sequences present in the raw reads.
2. Go from fasta/fastq to trimmed fasta/fastq (using FlexBar) 
    1. Flexbar performs several pre-processing steps prior to alignment including filtering out reads with uncalled bases, adaptor and barcode detection and removal, read separation based on barcode identification, and quality-based trimming from the right side of the read (as read quality tends to deteriorate toward the end/right side of a read). 
2. Trimmed fasta(fastq) to SAM file (aligned to genome using HISAT2)
    1. HISAT2 requires minimal memory requirements and alignment of reads to genome is rapid (typically within thirty minutes, depending on the provided parameters).
3. SAM to BAM (using samtools)
4. Quality control on aligned BAM file using SAMSTAT
    1. SAMSTAT output provides the number of reads that aligned to the genome (with varying degrees of confidence, indicated by a MAPQ score) and unmapped reads. Information on the rates of errors and mismatches among mapped reads is also included.
5. BAM to sorted BAM (using samtools)
    1. This step organizes the BAM file based on the position of the reads on the chromosome
6. Index sorted BAM for viewing (using samtools)
7. For each sample do steps 1 to 5
8. Get values of expression per gene (using featureCounts across all samples)
9. Calculate significantly changed genes (using R packages such as edgeR, DESeq)

#Steps to do INDEL/SNPS

1. Taked trimmed fasta(fastq) from flexbar (step 1 under Expression analysis)
2. Align to genome (using BWA)
    1. BWA is used in place of HISAT2 as it is better at mapping low-diverging sequences and is more sensitive for picking up SNPs and Indels
3. Add readgroup to SAM file using Packard
    1. Could be optional, based on protocol from people in CQS
4. Convert to BAM (using samtools)
5. Sort and index using SAMTOOLS
6. Quality control using SAMSTAT
7. gatk_intervals
    1. Converts seq data to intervals making it easier to analyze
8.  gatk_realignment(sample_base)
9.  gatk_recalibration(sample_base)
10. gatk_realign_recal(sample_base)
11. mark_dup(sample_base)
12. dup_index(sample_base)
13. Final entity search
    1. This is where you do it for SNP or INDEL

#Packages required
Installation
------------
[Flexbar] (https://github.com/seqan/flexbar)

[hisat2] (https://ccb.jhu.edu/software/hisat2/index.shtml)

[samtools] (http://www.htslib.org/)

[featurecount] (http://subread.sourceforge.net/)

[samstat] (http://samstat.sourceforge.net/)

#Files Required for INDEL/SNP Analysis
[Indel VCF File] (https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz)

[SNP VCF File] (https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz)