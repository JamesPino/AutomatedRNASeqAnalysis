*Steps to process RNAseq

Assumptions:
 
--Output_directory

----Fasta_files (orignal files)

----Bam_sam (sorted and indexed)

----Expression_Values

#Steps to run expression analysis

1. Go from fasta(or fastq) to trimmed fasta (using FlexBar) 
    1. (FASTQC to quality check, quality of sequence reads prior to aligning)
2. fasta trimmed to SAM (aligned to genome using HISAT2)
3. SAM to BAM (using SAMTOOLS)
    1. Quality control on aligned BAM using SAMSTAT
Tells percent of reads aligned and quality of alignment
4. BAM to sorted BAM (using SAMTOOLS)
5. Index sorted BAM for viewing (using SAMTOOLS)
6. For each sample do steps 1 to 5
7. Get values of expression per gene (using FeatureCount across all samples)
8. Calculate significantly changed genes (using R packages such as edgeR, DESeq)

#Steps to do INDEL/SNPS

1. Taked trimmed fasta(fastq) from flexbar (step 1 under Expression analysis)

2. Align to genome (using BWA)
    1. Ouput is a sam file
3. Add readgroup to SAM file using Packard
    1. Could be optional, based on protocol from people in CQS
4. Convert to BAM (using SAMTOOL)
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
