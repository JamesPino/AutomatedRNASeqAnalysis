cortez_config = dict(flexbar='/usr/local/bin/flexbar_v2.5_macosx/flexbar',
                     hisat2='/usr/local/bin/hisat2-2.0.4/hisat2',
                     samtools='/usr/local/bin/samtools',
                     featurecounts='/usr/local/bin/featureCounts',
                     samstat='/usr/local/bin/samstat',
                     gatk='/Users/temporary/Sources/GenomeAnalysisTK.jar',
                     bwa='/usr/local/bin/bwa-0.7.15/bwa',
                     picard='/Users/temporary/Sources/picard.jar', )

human_cortez_reference_config = dict(
    adaptors='/Users/temporary/genomes/adapters/illumina_truseq.fasta',
    transcripts='/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf',
    hisat2_index='/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/HISAT2Index.hisat_genome',
    bwa_index='/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome',
    reference_genome='/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa',
    indel_vcf_file='/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/Mills_and_1000G_gold_standard.indels.hg38.vcf',
    snp_vcf_file='/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/1000G_phase1.snps.high_confidence.hg38.vcf', )

mouse_cortez_reference_config = dict(
    adaptors='/Users/temporary/genomes/adapters/illumina_truseq.fasta',
    transcripts='/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf',
    hisat2_index='/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/HISAT2_index/mm10.genome',
    bwa_index='/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/version0.7.15/genome_indel',
    reference_genome='/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome_indel.fa',
    indel_vcf_file='/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/version0.7.15/mgp.v5.merged.indels.dbSNP142.normed_with_chr.vcf')