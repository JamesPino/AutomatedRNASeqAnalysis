#!/usr/bin/python

import os
import subprocess

# General information - modify appropriately for each experiment
# sample_base = "3612-DC-1"  # Naming scheme of samples without the sample number
number_of_samples = '4'  # Number of sequencing samples
species = 'human'  # Valid options: mouse or human
read_type = 'PE'  # Valid options: "SE" or "PE"
read_length = 150
sample_suffix = 'fastq'
compression_suffix = 'gz'
n_cpus = 8
pc = 'goku'
# experiment_name = '20160816_rnaseq'  # Name of experiment, also name of the output folder for all files
entity_searched = 'snp'  # Valid options include 'indel', 'snp', or 'none'
confidence_threshold = 20

elif pc == 'goku':  # Setting up programs
flexbar = '/usr/bin/flexbar'
hisat2 = '/home/pinojc/RNASeq_sources/Software/hisat2-2.0.4/hisat2'
samtools = '/home/pinojc/RNASeq_sources/Software/samtools-1.3.1/samtools'
samstat = '/usr/local/bin/samstat'
featurecounts = '/home/pinojc/RNASeq_sources/Software/subread-1.5.1-source/bin/featureCounts'
output_directory = "/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/LP_data"
fasta_directory = "/home/pinojc/LisaData/E65_rna_fasta"
gatk = '/home/pinojc/RNASeq_sources/Software/GenomeAnalysisTK.jar'
picard = '/home/pinojc/RNASeq_sources/Software/picard.jar'
STAR = '/home/pinojc/RNASeq_sources/Software/STAR'

# experiment specific information
adaptors = '/home/pinojc/RNASeq_sources/illumina_truseq.fasta'

if species == 'mouse':
    transcripts = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
    hisat2_index = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/HISAT2_index/mm10.genome'
    reference_genome = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome_indel.fa'

elif species == 'human':
    transcripts = '/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf'
    hisat2_index = '/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Sequence/HISAT2Index/genome'
    reference_genome = '/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa'
    snp_vcf_file = '/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/1000G_phase1.snps.high_confidence.hg38.vcf'
    indel_vcf_file = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/Mills_and_1000G_gold_standard.indels.hg38.vcf'
    genome_directory = '/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Sequence/STARIndex'

else:
    print(
        "Error - invalid species. Valid arguments are 'human' or 'mouse'")
    quit()  # Adaptor Trimming via Flexbar
def flexbar_trim(sample_base):
    # Flexbar options
    # '-a' = adaptor sequences (fasta format)
    # '-n' = number of threads
    # '-u' = max uncalled bases for each read to pass filtering
    # '-m' = min read length to remain after filtering/trimming
    # '-t' = prefix for output file names
    # '-ao' = adapter min overlap
    # '-ae' = adapter trim end
    print("Start trimming {}".format(sample_base))
    path_to_executable = flexbar
    suffix_for_output = '-t {}/{}-trimmed'.format(fasta_directory, sample_base)
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
        reads = '-r {}/{}.{}.{}'.format(fasta_directory, sample_base, sample_suffix, compression_suffix)

    elif read_type == 'PE':
        reads = ' -r {}/{}_1.{}.{} -p {}/{}_2.{}.{}'.format(fasta_directory, sample_base, sample_suffix,
                                                            compression_suffix, fasta_directory, sample_base,
                                                            sample_suffix, compression_suffix)

    command = [path_to_executable, reads, path_to_adaptors, threads, suffix_for_output, adaptor_overlap,
               adaptor_trim_end, number_max_uncalled_bases_pass, main_read_length_to_remain]
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
    print("Done trimming {}".format(sample_base))


# # STAR 2-pass alignemnt Alignment
# def star_file_read_in(sample_base):
#     path_to_executable = STAR
#     genome_dir = '--genomeDir {}'.format(genome_directory)
#     if read_type == 'SE':
#         reads = '--readFilesIn {}/{}-trimmed.{}'.format(fasta_directory,
#                                              sample_base, sample_suffix)
#     elif read_type == 'PE':
#         reads = '--readFilesIn {}/{}-trimmed_1.{} {}/{}-trimmed_2.{}'.format(fasta_directory, sample_base, sample_suffix,
#                                                                      fasta_directory, sample_base, sample_suffix)
#     threads = '--runThreadN {}'.format(n_cpus)
#     command = [path_to_executable, genome_dir, reads, threads]
#     call_code = ' '.join(command)
#     print(call_code)
#     process = subprocess.Popen([call_code], shell=True,
#                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#     while True:
#         output = process.stdout.readline()
#         if output == '' and process.poll() is not None:
#             break
#         if output:
#             print output.strip()
#         rc = process.poll()
#     print("Done aligning {}".format(sample_base))
#
#
# def star_initial_alignment(sample_base):
#     print("Starting aligning {}".format(sample_base))
#     if not os.path.exists('BAM_files'):
#         os.mkdir('BAM_files')
#
#     path_to_executable = STAR
#     run_mode = '--runMode genomeGenerate'
#     genome_dir = '--genomeDir {}'.format(genome_directory)
#     ref_genome = '--genomeFastaFiles {}'.format(reference_genome)
#     options = '--sjdbFileChrStartEnd {}/BAM_files --sjdbOverhang 75 --runThreadN {}'.format(output_directory, n_cpus)
#
#     output_name = '-S ./BAM_files/{}.sam'.format(sample_base)
#     threads = '-p {}'.format(n_cpus)
#     indices = '-x {}'.format(hisat2_index)
#
#     command = [path_to_executable, threads, indices, reads, output_name]
#     call_code = ' '.join(command)
#     print(call_code)
#     process = subprocess.Popen([call_code], shell=True,
#                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#     while True:
#         output = process.stdout.readline()
#         if output == '' and process.poll() is not None:
#             break
#         if output:
#             print output.strip()
#         rc = process.poll()
#     print("Done aligning {}".format(sample_base))
#
#
# def star_second_pass(sample_base):
#     print("Starting aligning {}".format(sample_base))
#     if not os.path.exists('BAM_files'):
#         os.mkdir('BAM_files')
#
#     path_to_executable = STAR
#     run_mode = '--runMode genomeGenerate'
#     genome_dir = '--genomeDir {}'.format(genome_directory)
#     ref_genome = '--genomeFastaFiles {}'.format(reference_genome)
#     options = '--sjdbFileChrStartEnd {}/BAM_files --sjdbOverhang 75 --runThreadN {}'.format(output_directory, n_cpus)
#
#     output_name = '-S ./BAM_files/{}.sam'.format(sample_base)
#     threads = '-p {}'.format(n_cpus)
#     indices = '-x {}'.format(hisat2_index)
#
#     command = [path_to_executable, threads, indices, reads, output_name]
#     call_code = ' '.join(command)
#     print(call_code)
#     process = subprocess.Popen([call_code], shell=True,
#                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#     while True:
#         output = process.stdout.readline()
#         if output == '' and process.poll() is not None:
#             break
#         if output:
#             print output.strip()
#         rc = process.poll()
#     print("Done aligning {}".format(sample_base))


def bwa_alignment(sample_base):
    # global_output = ''
    print("Starting alignment for {}".format(sample_base))
    if not os.path.exists('BWA_BAM_files'):
        os.mkdir('BWA_BAM_files')

    path_to_executable = '{} mem'.format(bwa)
    first_pair_reads = "{}/{}-trimmed_1.{}".format(fasta_directory, sample_base, sample_suffix)
    second_pair_reads = "{}/{}-trimmed_2.{}".format(fasta_directory, sample_base, sample_suffix)
    important_options = '-T 15 -M -t {}'.format(n_cpus)
    path_to_reference = reference_genome
    export_to_file = '> ./BWA_BAM_files/{}.sam'.format(sample_base)
    command = [path_to_executable, important_options, path_to_reference, first_pair_reads, second_pair_reads,
               export_to_file]
    call_code = ' '.join(command)
    print(call_code)
    process = subprocess.Popen([call_code], shell=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
            # if output:
            #         # global_output += output
            # print output.strip()
        rc = process.poll()
    print("Done aligning {}".format(sample_base))
    # with open('{}.sam'.format(sample_base)) as f:
    #     f.write(global_output)


# Conversion from SAM to BAM and sorting
def sam_to_bam(sample_base):
    print("Start sam to bam conversion {}".format(sample_base))
    path_to_executable = '{} view'.format(samtools)
    path_to_samples = '-S -b ./BAM_files/{}.sam'.format(sample_base)
    output_filename = '-o ./BAM_files/{}.bam'.format(sample_base)
    threads = '--threads {}'.format(n_cpus)
    command = [path_to_executable, path_to_samples, threads, output_filename]
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
    print("Done sam to bam conversion {}".format(sample_base))
    # os.remove('./BAM_files/{}.sam'.format(sample_base))


def bam_sort(sample_base):
    print("Start sorting {}".format(sample_base))
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
    print("Done sorting {}".format(sample_base))
    os.remove('./BAM_files/{}.bam'.format(sample_base))


def bwa_read_group(sample_base):
    print("Starting read group addition for {}".format(sample_base))
    path_to_executable = 'java -jar {}'.format(picard)
    picard_program = "AddOrReplaceReadGroups"
    input_files = 'I=./BWA_BAM_files/{}.sorted.bam'.format(sample_base)
    output_files = 'O=./BWA_BAM_files/{}.rg.sorted.bam'.format(sample_base)
    necessary_parameters = "RGID={} RGLB={} RGPL=ILLUMINA RGPU=ILLUMINA RGSM={}".format(sample_base, sample_base,
                                                                                        sample_base)
    command = [path_to_executable, picard_program, input_files, output_files, necessary_parameters]
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
    print("Done adding read group {}".format(sample_base))
    # os.remove('./BWA_BAM_files/{}.sam'.format(sample_base))


def mark_dup(sample_base):
    print("Starting mark dup for {}".format(sample_base))
    path_to_executable = "java -jar {}".format(picard)
    picard_program = 'MarkDuplicates'
    input_files = 'I=./BWA_BAM_files/{}.rg.sorted.bam'.format(sample_base)
    output_file = 'O=./BWA_BAM_files/{}.dedeuped.bam'
    options = 'CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=ouput.{}.metrics'.format(entity_searched)
    command = [path_to_executable, picard_program, input_files, output_file, options]
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
    print("Done with mark dup for {}".format(sample_base))


def gatk_splitntrim(sample_base):
    print("Starting split for {}".format(sample_base))
    path_to_executable = "java -jar {}".format(gatk)
    gatk_program = '-T SplitNCigarReads'
    path_to_reference = "-R {}".format(reference_genome)
    input_files = '-I ./BWA_BAM_files/{}.sorted.bam'.format(sample_base)
    output_file = '-o ./BWA_BAM_files/{}.split.bam'.format(sample_base)
    important_options = '-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS'
    command = [path_to_executable, gatk_program, path_to_reference, input_files, output_file, important_options]
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
    print("Done split for {}".format(sample_base))


def variant_calling(sample_base):
    print("Starting variant calling for {}".format(sample_base))
    path_to_executable = "java -jar {}".format(gatk)
    gatk_program = '-T HaplotypeCaller'
    path_to_reference = "-R {}".format(reference_genome)
    input_files = '-I ./BWA_BAM_files/{}.split.bam'.format(sample_base)
    options = '-dontUseSoftClippedBases -stand_all_conf {}'.format(confidence_threshold)
    output_file = './BWA_BAM_files/{}.{}.vcf'.format(sample_base, entity_searched)
    command = [path_to_executable, gatk_program, path_to_reference, input_files, output_file, important_options]
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
    print("Done split for {}".format(sample_base))
