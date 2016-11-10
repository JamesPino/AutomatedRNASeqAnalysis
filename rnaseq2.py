#!/usr/bin/python

###########################################################################################
# RNAseq analysis for Cortez Lab
# Authors: Lisa Poole and James Pino

# Requirements for this experiment - sequencing data files must be in a folder entitled
# "fasta" within a folder corresponding to the experiment name; modify general information
# in the first part of the script (starting with sample base).
##########################################################################################

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
pc = 'cortez_mac'
# experiment_name = '20160816_rnaseq'  # Name of experiment, also name of the output folder for all files
entity_searched = 'indel'  # Valid options include 'indel', 'snp', or 'none'

if pc == 'cortez_mac':
    # Setting up programs
    flexbar = '/usr/local/bin/flexbar_v2.5_macosx/flexbar'
    hisat2 = '/usr/local/bin/hisat2-2.0.4/hisat2'
    samtools = '/usr/local/bin/samtools'
    featurecounts = '/usr/local/bin/featureCounts'
    samstat = '/usr/local/bin/samstat'
    gatk = '/Users/temporary/Sources/GenomeAnalysisTK.jar'
    bwa = '/usr/local/bin/bwa-0.7.15/bwa'
    picard = '/Users/temporary/Sources/picard.jar'

    # experiment specific information
    output_directory = "/Users/temporary/projects/E65_RPE_RNA_seq"
    fasta_directory = "/Users/temporary/projects/E65_RPE_RNA_seq/fastq"

    # Reference files
    if species == 'mouse':
        adaptors = '/Users/temporary/genomes/adapters/illumina_truseq.fasta'
        transcripts = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
        hisat2_index = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/HISAT2_index/mm10.genome'
        bwa_index = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/version0.7.15/genome_indel'
        reference_genome = reference_genome = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome_indel.fa'
        indel_vcf_file = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/BWAIndex/version0.7.15/mgp.v5.merged.indels.dbSNP142.normed_with_chr.vcf'

    elif species == 'human':
        adaptors = '/Users/temporary/genomes/adapters/illumina_truseq.fasta'
        transcripts = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf'
        hisat2_index = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/HISAT2Index.hisat_genome'
        bwa_index = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome'
        reference_genome = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa'
        indel_vcf_file = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/Mills_and_1000G_gold_standard.indels.hg38.vcf'
        snp_vcf_file = '/Users/temporary/genomes/Homo_sapiens_hg38/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/1000G_phase1.snps.high_confidence.hg38.vcf'

    else:
        print("Error - invalid species. Valid arguments are 'human' or 'mouse'")
        quit()

elif pc == 'goku':
    # Setting up programs
    flexbar = '/usr/bin/flexbar'
    hisat2 = '/home/pinojc/RNASeq_sources/Software/hisat2-2.0.4/hisat2'
    samtools = '/home/pinojc/RNASeq_sources/Software/samtools-1.3.1/samtools'
    samstat = '/usr/local/bin/samstat'
    featurecounts = '/home/pinojc/RNASeq_sources/Software/subread-1.5.1-source/bin/featureCounts'
    output_directory = "/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/LP_data"
    fasta_directory = "/home/pinojc/LisaData/E65_rna_fasta"
    samtools = '/home/pinojc/RNASeq_sources/Software/samtools-1.3.1/samtools'
    gatk = '/home/pinojc/RNASeq_sources/Software/GenomeAnalysisTK.jar'
    bwa = '/home/pinojc/RNASeq_sources/Software/bwa.kit/bwa'
    picard = '/home/pinojc/RNASeq_sources/Software/picard.jar'
    samstat = '/usr/local/bin/samstat'

    # experiment specific information
    adaptors = '/home/pinojc/RNASeq_sources/illumina_truseq.fasta'
    output_directory = "/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/LP_data/"
    fasta_directory = "/home/pinojc/LisaData/E65_rna_fasta"

    if species == 'mouse':
        transcripts = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf'
        hisat2_index = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/HISAT2_index/mm10.genome'
        reference_genome = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome_indel.fa'

    elif species == 'human':
        transcripts = '/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf'
        hisat2_index = '/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Sequence/HISAT2Index/genome'
        bwa_index = '/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome'
        reference_genome = '/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa'

    else:
        print(
            "Error - invalid species. Valid arguments are 'human' or 'mouse'")
        quit()


# Adaptor Trimming via Flexbar
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


# HISAT Alignment
def hisat2_alignment(sample_base):
    print("Starting aligning {}".format(sample_base))
    if not os.path.exists('BAM_files'):
        os.mkdir('BAM_files')

    path_to_executable = hisat2
    output_name = '-S ./BAM_files/{}.sam'.format(sample_base)
    threads = '-p {}'.format(n_cpus)
    indices = '-x {}'.format(hisat2_index)

    if read_type == 'SE':
        reads = '-U {}/{}-trimmed.{}'.format(fasta_directory,
                                             sample_base, sample_suffix)
    elif read_type == 'PE':
        reads = '-1 {}/{}-trimmed_1.{} -2 {}/{}-trimmed_2.{}'.format(fasta_directory, sample_base, sample_suffix,
                                                                     fasta_directory, sample_base, sample_suffix)

    command = [path_to_executable, threads, indices, reads, output_name]
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
    os.remove('./BAM_files/{}.sam'.format(sample_base))


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


def bam_index(sample_base):
    print("Start indexing {}".format(sample_base))
    path_to_executable = '{} index'.format(samtools)
    path_to_samples = './BAM_files/{}.sorted.bam'.format(sample_base)
    threads = '--threads {}'.format(n_cpus)
    command = [path_to_executable, threads, path_to_samples]
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
    print("Done indexing {}".format(sample_base))


# SAMSTAT Quality Check
def samstat_analysis(sample_base):
    print("Start SAMSTAT check {}".format(sample_base))
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
    print("Done SAMSTAT check {}".format(sample_base))
    os.rename('./BAM_files/{}.sorted.bam.samstat.html'.format(sample_base),
              './SAMSTAT_analysis/{}.hisat2.sorted.bam.samstat.html'.format(sample_base))


# FeatureCounts - align reads to genes
def featurecounts_analysis(sample_base):
    # FeatureCounts options
    # -a annotation file
    # -o name of output file with read counts
    # -p fragments counted instead of reads (paired end specific)
    # -B count only read pairs that have both ends successfully aligned only
    # -C don't count reads that have two ends mapping to different chromosomes or mapping to same chromosome but on different strands
    # -Q minimum mapping quality score a read must satisfy in order to be counted
    # -t feature type in GTF file, default is "exon"
    # -g specify attribute in GTF file, default is gene_id

    path_to_executable = featurecounts
    annotation_file = "-a {}".format(transcripts)
    output_name = "-o gene_counts.txt"
    gtf_feature = '-t exon'
    gtf_attibute = '-g gene_id'
    quality_score = '-Q 30'
    out_string = ''
    for i in range(1, 5):
        input_files = ' ./BAM_files/{0}-{1}.sorted.bam'.format(sample_base, i)
        out_string += input_files
    if read_type == 'SE':
        important_options = "-T {}".format(n_cpus)
    elif read_type == 'PE':
        important_options = "-p -B -C -T {}".format(n_cpus)

    command = [path_to_executable, important_options, gtf_feature, gtf_attibute, quality_score, annotation_file,
               output_name, out_string]
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
    print("Done running {}".format(sample_base))


# INDEL/SNP Search Options - via BWA
def bwa_alignment(sample_base):
    global_output = ''
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


def bwa_read_group(sample_base):
    print("Starting read group addition for {}".format(sample_base))
    path_to_executable = 'java -jar {}'.format(picard)
    picard_program = "AddOrReplaceReadGroups"
    input_files = 'I=./BWA_BAM_files/{}.sam'.format(sample_base)
    output_files = 'O=./BWA_BAM_files/{}.rg.sam'.format(sample_base)
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
    os.remove('./BWA_BAM_files/{}.sam'.format(sample_base))


def bwa_sam_to_bam(sample_base):
    print("Starting sam to bam conversion for {}".format(sample_base))
    path_to_executable = '{} view'.format(samtools)
    path_to_samples = '-S -b ./BWA_BAM_files/{}.rg.sam'.format(sample_base)
    output_filename = '-o ./BWA_BAM_files/{}.bam'.format(sample_base)
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
    print("Done converting sam to bam conversion for {}".format(sample_base))
    os.remove('./BWA_BAM_files/{}.rg.sam'.format(sample_base))


def bwa_bam_sort(sample_base):
    print("Start sorting {}".format(sample_base))
    path_to_executable = '{} sort'.format(samtools)
    path_to_samples = './BWA_BAM_files/{}.bam'.format(sample_base)
    output_filename = '-o ./BWA_BAM_files/{}.sorted.bam'.format(sample_base)
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
    os.remove('./BWA_BAM_files/{}.bam'.format(sample_base))


def bwa_samstat_analysis(sample_base):
    print("Starting samstat analysis for {}".format(sample_base))

    path_to_executable = samstat
    path_to_samples = './BWA_BAM_filess/{}.sorted.bam'.format(sample_base)
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
    print("Done with samstat analysis for {}".format(sample_base))
    os.rename('./BWA_BAM_files/{}.sorted.bam.samstat.html'.format(sample_base),
              './SAMSTAT_analysis/{}.bwa.sorted.bam.samstat.html'.format(sample_base))


# Start of steps involving selection of 'indel' or 'snp'
def gatk_intervals(sample_base):
    print("Starting gatk intervals for {}".format(sample_base))
    path_to_executable = "java -jar {}".format(gatk)
    gatk_program = '-T RealignerTargetCreator'
    path_to_reference = "-R {}".format(reference_genome)
    input_files = '-I ./BWA_BAM_files/{}.sorted.bam'.format(sample_base)

    if entity_searched == 'indel':
        output_file = '-o ./BWA_BAM_files/{}.indel.intervals'.format(sample_base)
        path_to_vcf = "--known {}".format(indel_vcf_file)
    elif entity_searched == 'snp':
        output_file = '-o ./BWA_BAM_files/{}.snp.intervals'.format(sample_base)
        path_to_vcf = "--known {}".format(snp_vcf_file)

    command = [path_to_executable, gatk_program, path_to_reference, input_files, path_to_vcf, output_file]
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
    print("Done with gatk intervals for {}".format(sample_base))


def gatk_realignment(sample_base):
    print("Starting gatk realignment for {}".format(sample_base))
    path_to_executable = "java -jar {}".format(gatk)
    gatk_program = '-T IndelRealigner'
    path_to_reference = "-R {}".format(reference_genome)
    options = '--maxReadsForRealignment 999999 --maxReadsInMemory 999999'

    if entity_searched == 'indel':
        input_files = '-I ./BWA_BAM_files/{}.sorted.bam'.format(sample_base)
        intervals = '-targetIntervals ./BWA_BAM_files/{}.indel.intervals'.format(sample_base)
        output_file = '-o ./BWA_BAM_files/{}.indel.realigned.bam'.format(sample_base)
    elif entity_searched == 'snp':
        input_files = '-I ./BWA_BAM_files/{}.sorted.bam'.format(sample_base)
        intervals = '-targetIntervals ./BWA_BAM_files/{}.snp.intervals'.format(sample_base)
        output_file = '-o ./BWA_BAM_files/{}.snp.realigned.bam'.format(sample_base)

    command = [path_to_executable, gatk_program, path_to_reference, input_files, intervals, options, output_file]
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
    print("Done with gatk realignment for {}".format(sample_base))


def gatk_recalibration(sample_base):
    print("Starting gatk indel recalibration for {}".format(sample_base))
    path_to_executable = "java -jar {}".format(gatk)
    gatk_program = '-T BaseRecalibrator'
    path_to_reference = "-R {}".format(reference_genome)
    options = '-l INFO'

    if entity_searched == 'indel':
        input_files = '-I ./BWA_BAM_files/{}.indel.realigned.bam'.format(sample_base)
        output_file = '-o ./BWA_BAM_files/{}.indel.recal.table'.format(sample_base)
        path_to_vcf = "-knownSites {}".format(indel_vcf_file)
    elif entity_searched == 'snp':
        input_files = '-I BWA_BAM_files/{}.snp.realigned.bam'.format(sample_base)
        output_file = '-o ./BWA_BAM_files/{}.snp.recal.table'.format(sample_base)
        path_to_vcf = "--known {}".format(snp_vcf_file)

    command = [path_to_executable, gatk_program, path_to_reference, input_files, options, path_to_vcf, output_file]
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
    print("Done with gatk indel recalibration for {}".format(sample_base))


def gatk_realign_recal(sample_base):
    print("Starting gatk realign recal for {}".format(sample_base))
    path_to_executable = "java -jar {}".format(gatk)
    gatk_program = '-T PrintReads'
    path_to_reference = "-R {}".format(reference_genome)

    if entity_searched == 'indel':
        input_files = '-I ./BWA_BAM_files/{}.indel.realigned.bam'.format(sample_base)
        options = '-BQSR BWA_BAM_files/{}.indel.recal.table'.format(sample_base)
        output_file = '-o ./BWA_BAM_files/{}.indel.realigned.recal.bam'.format(sample_base)
    elif entity_searched == 'snp':
        input_files = '-I ./BWA_BAM_files/{}.snp.realigned.bam'.format(sample_base)
        options = '-BQSR ./BWA_BAM_files/{}.snp.recal.table'.format(sample_base)
        output_file = '-o ./BWA_BAM_files/{}.snp.realigned.recal.bam'.format(sample_base)

    command = [path_to_executable, gatk_program, path_to_reference, input_files, options, output_file]
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
    print("Done with gatk realign recal for {}".format(sample_base))
    # os.remove('./BWA_BAM_files/{}.indel.realigned.bam'.format(sample_base))


def mark_dup(sample_base):
    print("Starting mark dup for {}".format(sample_base))
    path_to_executable = "java -jar {}".format(picard)
    picard_program = 'MarkDuplicates'

    if entity_searched == 'indel':
        input_files = 'I=./BWA_BAM_files/{}.indel.realigned.recal.bam'.format(sample_base)
        output_file = 'O=./BWA_BAM_files/{}.indel.realigned.recal.dupmarked.bam'.format(sample_base)
        options = 'M=./BWA_BAM_files/{}-indel-marked_dup_metrics.txt'.format(sample_base)
    elif entity_searched == 'snp':
        input_files = 'I=./BWA_BAM_files/{}.snp.realigned.recal.bam'.format(sample_base)
        output_file = 'O=./BWA_BAM_files/{}.snp.realigned.recal.dupmarked.bam'.format(sample_base)
        options = 'M=./BWA_BAM_files/{}-snp-marked_dup_metrics.txt'.format(sample_base)

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


def dup_index(sample_base):
    print("Starting indel dup index for {}".format(sample_base))
    path_to_executable = '{} index'.format(samtools)

    if entity_searched == 'indel':
        path_to_samples = './BWA_BAM_files/{}.indel.realigned.recal.dupmarked.bam'.format(sample_base)
    elif entity_searched == 'snp':
        path_to_samples = './BWA_BAM_files/{}.snp.realigned.recal.dupmarked.bam'.format(sample_base)

    command = [path_to_executable, path_to_samples]
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
    print("Done with dup index for {}".format(sample_base))


def final_entity_search(sample_base):
    print("Starting final search for {}".format(sample_base))
    path_to_executable = "java -jar {}".format(gatk)
    gatk_program = '-T UnifiedGenotyper -l INFO'
    path_to_reference = "-R {}".format(reference_genome)
    options = '-A Coverage -A AlleleBalance -G Standard -stand_call_conf 50.0 -stand_emit_conf 10.0 -mbq 20 -deletions 0.05 -dcov 1000'

    if entity_searched == 'indel':
        search_item = '-glm INDEL'
        input_files = '-I ./BWA_BAM_files/{}.indel.realigned.recal.dupmarked.bam'.format(sample_base)
        path_to_vcf = "-D {}".format(indel_vcf_file)
        output_file = '--out {}.indel.vcf -metrics {}.indel.outmetrics.txt'.format(sample_base, sample_base)
    elif entity_searched == 'snp':
        search_item = '-glm SNP'
        input_files = '-I ./BWA_BAM_files/{}.snp.realigned.recal.dupmarked.bam'.format(sample_base)
        path_to_vcf = "-D {}".format(snp_vcf_file)
        output_file = '--out {}.snps.vcf -metrics {}.snp.outmetrics.txt'.format(sample_base, sample_base)

    command = [path_to_executable, gatk_program, path_to_reference, input_files, path_to_vcf, output_file, options,
               search_item]
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
    print("Done with final search for {}".format(sample_base))


def RNAseq_analysis(sample_base):
    flexbar_trim(sample_base)
    hisat2_alignment(sample_base)
    sam_to_bam(sample_base)
    bam_sort(sample_base)
    bam_index(sample_base)
    samstat_analysis(sample_base)
    if entity_searched == 'snp':
        bwa_genome_search(sample_base)
    elif entity_searched == 'indel':
        bwa_genome_search(sample_base)
    elif entity_searched == 'none':
        quit()


def setup_bwa(sample_base):
    bwa_alignment(sample_base)
    bwa_read_group(sample_base)
    bwa_sam_to_bam(sample_base)
    bwa_bam_sort(sample_base)
    bwa_samstat_analysis(sample_base)


def bwa_genome_search(sample_base):
    # if cont:
    setup_bwa(sample_base)
    gatk_intervals(sample_base)
    gatk_realignment(sample_base)
    gatk_recalibration(sample_base)
    gatk_realign_recal(sample_base)
    mark_dup(sample_base)
    dup_index(sample_base)
    final_entity_search(sample_base)
    # if entity_searched == 'indel':
    #     os.remove('./BWA_BAM_files/{}.indel.realigned.recal.bam'.format(sample_base))
    # elif entity_searched == 'snp':
    #     os.remove('./BWA_BAM_files/{}.snp.realigned.recal.bam'.format(sample_base))
    return False


def gene_exp(sample_base):
    """ Runs basic rna analysis


    :param sample_base:
    :return:
    """
    flexbar_trim(sample_base)
    hisat2_alignment(sample_base)
    sam_to_bam(sample_base)
    bam_sort(sample_base)
    bam_index(sample_base)
    samstat_analysis(sample_base)


def RNAseq_analysis(sample_base):
    gene_exp(sample_base)
    cont = True

    if entity_searched == 'snp':
        cont = bwa_genome_search(sample_base, cont)
    elif entity_searched == 'indel':
        cont = bwa_genome_search(sample_base, cont)
    elif entity_searched == 'none':
        quit()

    # def do_snp(sample_base):
    # does only snp

    # def do_indel(sample_base):
    # does only indel

    # def protocol_both(sample_base):
    #     do_snp()
    #     do_indel()

    # def protocol_1(sample_base):
    #     gene_exp()
    #     protocol_both()

    # def protocol_2(sample_base):
    # already did gene_exp, no need to do again
    # protocol_both()


os.chdir(output_directory)
wd = os.getcwd()
print(wd)

# RNAseq_analysis('3612-DC-1')
# RNAseq_analysis('3612-DC-2')
# RNAseq_analysis('3612-DC-3')
# RNAseq_analysis('3612-DC-4')
# featurecounts_analysis()

# bwa_genome_search('3612-DC-1')
bwa_genome_search('3612-DC-2')
bwa_genome_search('3612-DC-3')
bwa_genome_search('3612-DC-4')
