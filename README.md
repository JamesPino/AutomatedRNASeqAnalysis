Steps to process RNAseq

Assumptions:
 
--Output_directory

----Fasta_files (orignal files)

----Bam_sam (sorted and indexed)

----Expression_Values

#Steps to run

1. Go from Fasta to trimmed fasta (using FlexBar) 
..1. (FASTQC to quality check, quality of sequence reads prior to aligning)
2. Trimmed to SAM (aligned to genome using HISAT2)
3. SAM to BAM (using SAMTOOLS)
3.1 Quality control on aligned BAM using SAMSTAT
Tells percent of reads aligned and quality of alignment
4. BAM to sorted BAM (using SAMTOOLS)
5. Index sorted BAM for viewing (using SAMTOOLS)
6. Do 1 to 5 per sample (output featurecounts per sample)
7. Get values of expression per gene (using FeatureCount across all samples)
8. Calculate significantly changed genes (using R packages such as edgeR, DESeq)


#Packages required
Installation
------------
[Flexbar] (https://github.com/seqan/flexbar)

[hisat2] (https://ccb.jhu.edu/software/hisat2/index.shtml)

[samtools] (http://www.htslib.org/)

[featurecount] (http://subread.sourceforge.net/)

[samstat] (http://samstat.sourceforge.net/)