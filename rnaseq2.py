#!/usr/bin/python

import os
import gzip
import subprocess

# General information - modify appropriately for each experiment
#sample_base = "3612-DC-1"
species = 'human'  # mouse or human
read_type = 'PE'  # Valid options: "SE" or "PE"
read_length = 150
sample_suffix = 'fastq'
n_cpus = 8
pc = 'goku'


if pc == 'lisa':
    # Setting up programs
    flexbar = '/usr/local/bin/flexbar_v2.5_macosx'
    hisat2 = '/usr/local/bin/hisat2-2.0.4'
    samtools = '/usr/local/bin/samtools'
    featurecounts = '/usr/local/bin/featureCounts'
    samstat = '/usr/local/bin/samstat'

    #experiment specific information
    output_directory = "/Users/temporary/projects/20160816_rnaseq"
    fasta_directory = "/Users/temporary/projects/20160816_rnaseq/fasta"

    # Reference files
    if species == 'mouse':
        adaptors = '/Users/temporary/genomes/adapters/illumina_truseq.fasta'
        transcripts = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
        hisat2_index = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/HISAT2_index/mm10.genome'
        reference_genome = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome_indel.fa'

    elif species == 'human':
        adaptors = '/Users/temporary/genomes/adapters/illumina_truseq.fasta'
        transcripts = '/Users/temporary/genomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf'
        hisat2_index = '/Users/temporary/genomes/Homo_sapiens/UCSC/hg19/Sequence/HISAT2_index/hg19.genome'
        reference_genome = '/Users/temporary/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa'

    else:
        print("Error - invalid species. Valid arguments are 'human' or 'mouse'")
        quit()

elif pc == 'goku':
    #Setting up programs
    flexbar = '/usr/bin/flexbar'
    hisat2 = '/home/pinojc/RNASeq_sources/Software/hisat2-2.0.4/hisat2'
    samtools = '/home/pinojc/RNASeq_sources/Software/samtools-1.3.1/samtools'
    samstat = '/usr/local/bin/samstat'
    featurecounts = '/home/pinojc/RNASeq_sources/Software/subread-1.5.1-source/bin/featureCounts'
    output_directory = "/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/LP_data"
    fasta_directory = "/home/pinojc/LisaData/E65_rna_fasta"

    #experiment specific information
    adaptors = '/home/pinojc/RNASeq_sources/illumina_truseq.fasta'
    output_directory = "/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/LP_data"
    fasta_directory = "/home/pinojc/LisaData/E65_rna_fasta"

    if species == 'mouse':
        transcripts = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
        hisat2_index = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/HISAT2_index/mm10.genome'
        reference_genome = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome_indel.fa'

    elif species == 'human':
        transcripts = '/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf'
        hisat2_index = '/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Sequence/HISAT2Index/genome'
        reference_genome = '/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa'

    else:
        print(
            "Error - invalid species. Valid arguments are 'human' or 'mouse'")
        quit()


# Adaptor Trimming via Flexbar

# Flexbar options
# '-a' = adaptor sequences (fasta format)
# '-n' = number of threads
# '-u' = max uncalled bases for each read to pass filtering
# '-m' = min read length to remain after filtering/trimming
# '-t' = prefix for output file names
# '-ao' = adapter min overlap
# '-ae' = adapter trim end

def flexbar_trim():
    path_to_executable = flexbar
    suffix_for_output = '-t {}/{}-trimmed'.format(fasta_directory,
                                                    sample_base)
    adaptor_trim_end = '-ae ANY'
    adaptor_overlap = '-ao 5'
    path_to_adaptors = "-a {}".format(adaptors)
    number_max_uncalled_bases_pass = '-u {}'.format(read_length)
    main_read_length_to_remain = '-m 18'
    threads = '-n {}'.format(n_cpus)

    if read_type not in ('SE', 'PE'):
        print("Error. Invalid read type. Valid options are SE or PE Aborting.")
        quit()

    elif read_type == 'SE':
        path_to_reads = '-r {}/{}.{}'.format(fasta_directory, sample_base,
                                             sample_suffix)

        command = [path_to_executable, path_to_reads, path_to_adaptors, threads, suffix_for_output, adaptor_trim_end, adaptor_overlap,
                   number_max_uncalled_bases_pass, main_read_length_to_remain,
                   parallel_cores]

    elif read_type == 'PE':

        reads = ' -r {}/{}_1.{}'.format(fasta_directory, sample_base,
                                       sample_suffix)
        paired_reads = '-p {}/{}_2.{}'.format(fasta_directory, sample_base,
                                              sample_suffix)

        command = [path_to_executable, reads, paired_reads, path_to_adaptors, threads, suffix_for_output, adaptor_overlap, adaptor_trim_end,
                   number_max_uncalled_bases_pass, main_read_length_to_remain]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                      stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()

# HISAT Alignment
def hisat2_alignment(sample_base):
    print("Starting to do {}".format(sample_base))
    if not os.path.exists('BAM_files'):
        os.mkdir('BAM_files')

    path_to_executable = hisat2
    output_name = '-S BAM_files/{}.sam'.format(sample_base)
    threads = '-p {}'.format(n_cpus)
    indices = '-x {}'.format(hisat2_index)
    if read_type == 'SE':
        unpaired_sample = '-U {}/{}-trimmed.fastq'.format(fasta_directory,
                                                         sample_base)
        command = [path_to_executable, threads, indices, unpaired_sample,
                   output_name]

    elif read_type == 'PE':
        left_fasta_sample = '-1 {}/{}-trimmed_1_1.fastq'.format(fasta_directory,
                                                           sample_base)
        right_fasta_sample = '-2 {}/{}-trimmed_1_2.fastq'.format(fasta_directory,
                                                            sample_base)
        command = [path_to_executable, threads, indices,
                   left_fasta_sample, right_fasta_sample, output_name]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Done running {}".format(sample_base))



def bam_sort(sample_base):
    path_to_executable = '{} sort'.format(samtools)
    path_to_samples = './BAM_files/{}.bam'.format(sample_base)
    output_filename = '-o ./BAM_files/{}.sorted.bam'.format(sample_base)
    threads = '--threads {}'.format(n_cpus)
    command = [path_to_executable, threads, output_filename, path_to_samples]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Done running {}".format(sample_base))

# Conversion from SAM to BAM and indexing
def sam_to_bam(sample_base):
    path_to_executable = '{} view'.format(samtools)
    path_to_samples = '-S -b BAM_files/{}.sam'.format(sample_base)
    output_filename = '-o BAM_files/{}.bam'.format(sample_base)
    threads = '--threads {}'.format(n_cpus)
    command = [path_to_executable, path_to_samples,threads, output_filename]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Done running {}".format(sample_base))
    #os.remove('BAM_files/{}.sam'.format(sample_base))


# SAMSTAT Quality Check
def samstat_analysis(sample_base):
    if not os.path.exists('SAMSTAT_analysis'):
        os.mkdir('SAMSTAT_analysis')

    path_to_executable = samstat
    path_to_samples = './BAM_files/{}.sorted.bam'.format(sample_base)
    command = [path_to_executable, path_to_samples]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print output.strip()
        rc = process.poll()
    print("Done running {}".format(sample_base))
    os.rename('./BAM_files/{}.sorted.bam.samstat.html'.format(sample_base),
              './SAMSTAT_analysis/{}.sorted.bam.samstat.html'.format(sample_base))


# FeatureCounts - align reads to genes
# FeatureCounts options
# -a annotation file
# -o name of output file with read counts
# -s <int> perform strand specific read counting, possible values: 0 (unstranded), 1 (stranded), and 2 (reversely stranded)
# -p fragments counted instead of reads (paired end specific)
# -B count only read pairs that have both ends successfully aligned only
# -C don't count reads that have two ends mapping to different chromosomes or mapping to same chromosome but on different strands
# -Q minimum mapping quality score a read must satisfy in order to be counted
# -t feature type in GTF file, default is "exon"
# -g specify attribute in GTF file, default is gene_id
def featurecounts_analysis(sample_base):
    if read_type == 'SE':
        path_to_executable = featurecounts
        annotation_file = "-a 'transcripts'"
        output_name = "-o gene_counts.txt"
        important_options = "-T 'n_cpus'"
        gtf_feature = '-t exon'
        gtf_attibute = '-g gene_id'
        quality_score = '-Q 30'
        input_files = '/BAM_files/{}.sorted.bam'.format(sample_base)
        command = [path_to_executable, important_options, gtf_feature,
                   gtf_attibute, quality_score, annotation_file,
                   output_name, input_files]
        print(command)
        subprocess.Popen(command,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    elif read_type == 'PE':
        path_to_executable = featurecounts
        annotation_file = "-a {}".format(transcripts)
        output_name = "-o gene_counts.txt"
        important_options = "-p -B -C -T {}".format(n_cpus)
        gtf_feature = '-t exon'
        gtf_attibute = '-g gene_id'
        quality_score = '-Q 30'
        out_string = ''
        for i in range(1,5):
            input_files = ' ./BAM_files/{0}-{1}.sorted.bam'.format(sample_base, i)
            out_string+= input_files


        command = [path_to_executable, important_options, gtf_feature,
                   gtf_attibute, quality_score, annotation_file, output_name,
                   out_string]
        call_code = ' '.join(command)
        print(call_code)
        process = subprocess.Popen([call_code], shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                print output.strip()
            rc = process.poll()
        print("Done running {}".format(sample_base))

        quit()


os.chdir(output_directory)
wd = os.getcwd()
print(wd)

# flexbar_trim()
# quit()
# hisat2_alignment("3612-DC-1")
# hisat2_alignment("3612-DC-2")
# hisat2_alignment("3612-DC-3")
# hisat2_alignment("3612-DC-4")
# quit()

# sam_to_bam("3612-DC-1")
# sam_to_bam("3612-DC-2")
# sam_to_bam("3612-DC-3")
# sam_to_bam("3612-DC-4")

# bam_sort("3612-DC-1")
# bam_sort("3612-DC-2")
# bam_sort("3612-DC-3")
# bam_sort("3612-DC-4")
# quit()

# bam_index()
# quit()

# samstat_analysis("3612-DC-1")
# samstat_analysis("3612-DC-2")
# samstat_analysis("3612-DC-3")
# samstat_analysis("3612-DC-4")
# quit()
sample_base = "3612-DC"
featurecounts_analysis(sample_base)
