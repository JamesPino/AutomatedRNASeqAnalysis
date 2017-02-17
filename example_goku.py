from RNASeqAnalyzer import RNASeqAnalyzer


if __name__ == '__main__':
    goku_ref_config = dict(
            adaptors='/mnt/d/RNA_Seq_sources/UCSC/hg38/illumina_truseq.fasta',
            transcripts='/mnt/d/RNA_Seq_sources/UCSC/hg38/Annotation/Genes/genes.gtf',
            hisat2_index='/mnt/d/RNA_Seq_sources/UCSC/hg38/Sequence/HISAT2Index/genome',
            reference_genome='/mnt/d/RNA_Seq_sources/UCSC/hg38/Sequence/BWAIndex/genome.fa', )

    goku_config = dict(
            featurecounts='/home/pinojc/Sources/subread-1.5.1-source/bin/featureCounts',
            samstat='samstat', samtools='samtools',
            hisat2='/home/pinojc/Sources/hisat2-2.0.5/hisat2',
            flexbar='/home/pinojc/Sources/flexbar_v2.5_linux64/flexbar')

    x = RNASeqAnalyzer('3612-DC-1',
                       '/mnt/d/LP_data/Test',
                       4,
                       ref_config=goku_ref_config,
                       executables_config=goku_config,
                       write_bash=False)
    # x.gene_exp()
    # x.flexbar_trim()
    x.hisat2_alignment()
    with open('run_1.sh', 'w') as f:
        f.write(x._bash_file)
    # print(x._bash_file)