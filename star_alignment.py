#!/usr/bin/python
import os
import subprocess

number_of_samples = '4'  # Number of sequencing samples
species = 'human'  # Valid options: mouse or human
read_type = 'PE'  # Valid options: "SE" or "PE"
read_length = 150
sample_suffix = 'fastq'
compression_suffix = 'gz'
n_cpus = 8
pc = 'puma'

if pc == 'puma':
    # Setting up programs
    STAR = '/home/pinojc/Sources/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR'
    STAR_fusion = '/home/pinojc/Sources/STAR-Fusion_v1.1.0/STAR-Fusion'
    CTAT_lib = '/home/pinojc/Sources/GRCh38_gencode_v26_CTAT_lib_July192017'

    # experiment specific information
    output_directory = "/home/pinojc/lisa_rnaseq"
    fasta_directory = "/home/pinojc/lisa_rnaseq/fastq"

    if species == 'human':
        STAR_index = '/home/pinojc/RNA_Seq/data/GRCh38/star_indices_overhang100'
        # reference_genome = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome_indel.fa'


def STAR_alignment(sample_base):
    print("Starting aligning {}".format(sample_base))
    # if not os.path.exists('{}/BAM_files'.format(output_directory)):
    #     os.mkdir('{}/BAM_files'.format(output_directory))

    path_to_executable = STAR
    threads = '--runThreadN {}'.format(n_cpus)
    indices = '--genomeDir {}'.format(STAR_index)
    fusion_options = '--twopassMode Basic --outReadsUnmapped None --chimSegmentMin 12 --chimJunctionOverhangMin 12 ---align SJDBoverhangMin 10 --alignMatesGapMax 100000 --alignIntronMax 100000 --chimSegmentReadGapMax parameter 3 --alignSJstitchMismatchNmax 5 -1 5 5 --limitBAMsortRAM 31532137230 --outSAMtype BAM SortedByCoordinate'

    if read_type == 'SE':
        reads = '--readFilesIn {}/{}-trimmed_1.{}'.format(fasta_directory,
                                                          sample_base, sample_suffix)
    elif read_type == 'PE':
        reads = '--readFilesIn {0}/{1}-trimmed_1.{2} {0}/{1}-trimmed_2.{2}'.format(fasta_directory, sample_base,
                                                                                   sample_suffix)

    command = [path_to_executable, indices, reads, fusion_options, threads]
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
    print("Done aligning {}".format(sample_base))


def STAR_fusion_search(sample_base):
    print("Starting fusion search for {}".format(sample_base))
    path_to_executable = STAR_fusion
    threads = '--runThreadN {}'.format(n_cpus)
    CTAT = '--genome_lib_dir {}'.format(CTAT_lib)
    output = '--output_dir {}'.format(output_directory)

    if read_type == 'SE':
        reads = '--left_fq {}/{}-trimmed_1.{}'.format(fasta_directory,
                                                      sample_base, sample_suffix)
    elif read_type == 'PE':
        reads = '--left_fq {0}/{1}-trimmed_1.{2} --right_fq {0}/{1}-trimmed_2.{2}'.format(fasta_directory, sample_base,
                                                                                          sample_suffix)

    command = [path_to_executable, CTAT, reads, output]
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
    print("Done with fusion search for {}".format(sample_base))


def STAR_a(sample):
    for i in range(1, 5):
        # STAR_alignment('{}-{}'.format(sample, i))
        STAR_fusion_search('{}-{}'.format(sample, i))


STAR_a('3612-DC')
