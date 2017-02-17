#!/usr/bin/python

"""
RNAseq analysis for Cortez Lab
Authors: Lisa Poole and James Pino

Requirements for this experiment - sequencing data files must be in a folder
entitled "fasta" within a folder corresponding to the experiment name;
modify general information in the first part of the script
(starting with sample base).
"""
import time
import os
import subprocess
import warnings

# General information - modify appropriately for each experiment
# Number of sequencing samples
number_of_samples = '4'
# Valid options: mouse or human
species = 'human'
# Valid options: "SE" or "PE"
read_type = 'PE'
read_length = 150
sample_suffix = 'fastq'
compression_suffix = 'gz'
n_cpus = 8

source_config = dict(flexbar=None, hisat2=None, samtools=None,
                     featurecounts=None, samstat=None, )

reference_config = dict(adaptors=None, transcripts=None, hisat2_index=None,
                        reference_genome=None, )

goku_config = dict(flexbar='/usr/bin/flexbar',
                   hisat2='/home/pinojc/RNASeq_sources/Software/hisat2-2.0.4/hisat2',
                   samtools='/home/pinojc/RNASeq_sources/Software/samtools-1.3.1/samtools',
                   samstat='/usr/local/bin/samstat',
                   featurecounts='/home/pinojc/RNASeq_sources/Software/subread-1.5.1-source/bin/featureCounts',
                   fasta_directory="/home/pinojc/LisaData/E65_rna_fasta",
                   gatk='/home/pinojc/RNASeq_sources/Software/GATK/GenomeAnalysisTK.jar',
                   bwa='/home/pinojc/RNASeq_sources/Software/bwa.kit/bwa',
                   picard='/home/pinojc/RNASeq_sources/Software/picard.jar', )


class RNASeqAnalyzer(object):
    def __init__(self, sample_base, output_directory, n_cpu, ref_config,
                 executables_config, write_bash=False):
        self.sample_base = sample_base
        self.n_cpu = n_cpu
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        self.fasta_dir = os.path.join(output_directory, 'fasta')
        if not os.path.exists(self.fasta_dir):
            warnings.warn("Must have directory within provided output "
                          "directory {} named 'fasta' that contains "
                          "fasta(or fastq) files with sample base name")
        os.chdir(output_directory)
        self.all_output = ''

        for i in reference_config:
            if i not in ref_config:
                print("Please provide path to {} in ref_config")
        for i in ref_config:
            if not os.path.exists(ref_config[i]):
                err = "{} path of  {} does not exist".format(i, ref_config[i])
                warnings.warn(err)
        self.adaptors = ref_config['adaptors']
        self.transcripts = ref_config['transcripts']
        self.hisat2_index = ref_config['hisat2_index']
        self.reference_genome = ref_config['reference_genome']
        self.write_bash = write_bash
        self._bash_file = ''

        for i in source_config:
            if i not in executables_config:
                print("Please provide path to {} in executables_config".format(
                        i))
        for i in executables_config:
            which(executables_config[i])

        self._exe = executables_config
        self.flexbar = self._exe['flexbar']
        self.hisat2 = self._exe['hisat2']
        self.samstat = self._exe['samstat']
        self.samtools = self._exe['samtools']
        self.featurecounts = self._exe['featurecounts']

    def setup(self):
        """
        Create directory of all outputs

        Returns
        -------

        """

    def flexbar_trim(self):
        """
        Removed adapters for better and quicker alignment

        input :
            zipped fasta(fastq) file
            If single-ended, only 1, if paired-end there will be 2
            Labeled with sample_base_1, sample_base_2

        Returns
        -------

        """
        # Flexbar options
        # '-a' = adaptor sequences (fasta format)
        # '-n' = number of threads
        # '-u' = max uncalled bases for each read to pass filtering
        # '-m' = min read length to remain after filtering/trimming
        # '-t' = prefix for output file names
        # '-ao' = adapter min overlap
        # '-ae' = adapter trim end
        print("Start trimming {}".format(self.sample_base))

        warnings.warn("Assuming the use of Illumina sequence adapters!\n"
                      "Please verify")

        path_to_executable = self._exe['flexbar']

        suffix_for_output = '-t {}/{}-trimmed'.format(self.fasta_dir,
                                                      self.sample_base)
        adaptor_trim_end = '-ae ANY'
        adaptor_overlap = '-ao 5'
        path_to_adaptors = "-a {}".format(self.adaptors)
        number_max_uncalled_bases_pass = '-u {}'.format(read_length)
        main_read_length_to_remain = '-m 18'
        threads = '-n {}'.format(self.n_cpu)

        if read_type not in ('SE', 'PE'):
            print("Error. Invalid read type. Valid options are "
                  "SE or PE Aborting.")
            quit()

        elif read_type == 'SE':
            reads = '-r {}/{}.{}.{}'.format(self.fasta_dir, self.sample_base,
                                            sample_suffix, compression_suffix)

        elif read_type == 'PE':
            reads = ' -r {0}/{1}_1.{2}.{3} -p {0}/{1}_2.{2}.{3}'.format(
                    self.fasta_dir, self.sample_base, sample_suffix,
                    compression_suffix)

        command = [path_to_executable, reads, path_to_adaptors, threads,
                   suffix_for_output, adaptor_overlap, adaptor_trim_end,
                   number_max_uncalled_bases_pass, main_read_length_to_remain]
        if self.write_bash:
            self._bash_file += ' '.join(i for i in command)
            self._bash_file += '\n'
        else:
            self._run(command)
        print("Done trimming {}".format(self.sample_base))

    # HISAT Alignment
    def hisat2_alignment(self):
        """
        Aligns trimmed fasta(fastq) to genome

        #TODO check to see if genome fasta exists
        #TODO make sure HISAT2 index exists for that genome
        Can be done with HISAT2 index command
        Returns
        -------
        outputs a SAM format

        """
        print("Starting aligning {}".format(self.sample_base))
        path_to_executable = self._exe['hisat2']
        if not os.path.exists("BAM_files"):
            print("Creating BAM_files directory")
            os.mkdir("BAM_files")
        output_name = '-S ./BAM_files/{}.sam'.format(self.sample_base)
        threads = '-p {}'.format(self.n_cpu)

        indices = '-x {}'.format(self.hisat2_index)

        if read_type == 'SE':
            reads = '-U {0}/{1}-trimmed.{2}'.format(self.fasta_dir,
                                                    self.sample_base,
                                                    sample_suffix)
        elif read_type == 'PE':
            reads = '-1 {0}/{1}-trimmed_1.{2} ' \
                    '-2 {0}/{1}-trimmed_2.{2}'.format(self.fasta_dir,
                                                        self.sample_base,
                                                        sample_suffix, )

        command = [path_to_executable, threads, indices, reads, output_name]
        if self.write_bash:
            self._bash_file += ' '.join(i for i in command)
            self._bash_file += '\n'
        else:
            self._run(command)
        print("Done aligning {}".format(self.sample_base))

    def sam_to_bam(self):
        """ Conversion from SAM to BAM and sorting

        Returns
        -------

        """
        print("Start sam to bam conversion {}".format(self.sample_base))
        path_to_executable = '{} view'.format(self.samtools)
        path_to_samples = '-S -b ./BAM_files/{}.sam'.format(self.sample_base)
        output_filename = '-o ./BAM_files/{}.bam'.format(self.sample_base)
        threads = '--threads {}'.format(self.n_cpu)
        command = [path_to_executable, path_to_samples, threads,
                   output_filename]
        if self.write_bash:
            self._bash_file += ' '.join(i for i in command)
            self._bash_file += '\n'
        else:
            self._run(command)

    def bam_index(self):
        """ Indexs SAM file to BAM


        Returns
        -------

        """
        print("Start indexing {}".format(self.sample_base))
        # samtools = '/home/pinojc/RNASeq_sources/Software/./sambamba_v0.6.5'
        # TODO see why I chose to use SAMBA? Benchmark
        path_to_executable = '{} index'.format(self.samtools)
        path_to_samples = './BAM_files/{}.sorted.bam'.format(self.sample_base)
        threads = '-p --nthreads={}'.format(self.n_cpu)
        threads = ''
        command = [path_to_executable, threads, path_to_samples]
        if self.write_bash:
            self._bash_file += ' '.join(i for i in command)
            self._bash_file += '\n'
        else:
            self._run(command)
        print("Done indexing {}".format(self.sample_base))

    # SAMSTAT Quality Check
    def samstat_analysis(self):
        """ quality control check on alignment
        Open file and see what is print what is good and bad

        Returns
        -------

        """
        print("Start SAMSTAT check {}".format(self.sample_base))
        if not os.path.exists('quality_control'):
            os.mkdir('quality_control')

        path_to_executable = self.samstat
        path_to_samples = './BAM_files/{}.sorted.bam'.format(self.sample_base)
        command = [path_to_executable, path_to_samples]
        if self.write_bash:
            self._bash_file += ' '.join(i for i in command)
            self._bash_file += '\n'
        else:
            self._run(command)
            os.rename('BAM_files/{}.sorted.bam.samstat.html'.format(
                    self.sample_base),
                    'quality_control/{}.hisat2.sorted.bam.samstat.html'.format(
                            self.sample_base))
        print("Done SAMSTAT check {}".format(self.sample_base))

    def bam_sort(self, directory):
        """ sorts bam file and collapses duplicates

        append .sorted to  "sample_base_name.bam"
            like "sample_base_name.sorted.bam"

        Parameters
        ----------
        directory : str
            output directory

        Returns
        -------

        """
        print("Start sorting {}".format(self.sample_base))
        path_to_executable = '{} sort'.format(self.samtools)
        path_to_samples = ' {}/{}.bam'.format(directory, self.sample_base)
        output_filename = '-o {}/{}.sorted.bam'.format(directory,
                                                       self.sample_base)
        threads = '--threads {}'.format(self.n_cpu)
        command = [path_to_executable, threads, output_filename,
                   path_to_samples]
        if self.write_bash:
            self._bash_file += ' '.join(i for i in command)
            self._bash_file += '\n'
        else:
            self._run(command)
        print("Done sorting {}".format(self.sample_base))

    # FeatureCounts - align reads to genes
    def featurecounts_analysis(self, list_of_samples):
        """
        Takes all sample sorted bam files and compares all genes across samples
        Need GTF file for gene ids ( downloaded with genome)
        #TODO check for GTF file
        Parameters
        -------
        list_of_samples: list_like
            list of all samples
        Returns
        -------

        """

        # FeatureCounts options
        # -a annotation file
        # -o name of output file with read counts
        # -p fragments counted instead of reads (paired end specific)
        # -B count only read pairs that have both ends successfully aligned only
        # -C don't count reads that have two ends mapping to different chromosomes or mapping to same chromosome but on different strands
        # -Q minimum mapping quality score a read must satisfy in order to be counted
        # -t feature type in GTF file, default is "exon"
        # -g specify attribute in GTF file, default is gene_id

        path_to_executable = self.featurecounts
        annotation_file = "-a {}".format(self.transcripts)
        output_name = "-o gene_counts.txt"
        gtf_feature = '-t exon'
        gtf_attibute = '-g gene_id'
        quality_score = '-Q 30'
        out_string = ''
        if list_of_samples is None:
            print("Need list of samples")
            quit()

        for i in list_of_samples:
            input_files = ' ./BAM_files/{0}.sorted.bam'.format(i)
            out_string += input_files
        if read_type == 'SE':
            important_options = "-T {}".format(n_cpus)
        elif read_type == 'PE':
            important_options = "-p -B -C -T {}".format(n_cpus)

        command = [path_to_executable, important_options, gtf_feature,
                   gtf_attibute, quality_score, annotation_file, output_name,
                   out_string]
        if self.write_bash:
            self._bash_file += ' '.join(i for i in command)
            self._bash_file += '\n'
        else:
            self._run(command)
        print("Done running {}".format(self.sample_base))

    def _run(self, command):
        st = time.time()
        call_code = ' '.join(command)
        print("Running command {}".format(call_code))
        process = subprocess.Popen([call_code], shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                print(output.strip())
                self.all_output += output.strip()
            rc = process.poll()
        print('Finished - time taken = {}'.format(time.time() - st))

    def gene_exp(self):
        """ Runs basic rna analysis


        :param sample_base:
        :return:
        """
        self.flexbar_trim()
        self.hisat2_alignment()
        self.sam_to_bam()
        self.bam_sort('BAM_files')
        self.bam_index()
        self.samstat_analysis()

    def write_log_file(self, file_name):
        with open(file_name, 'w') as f:
            f.write(self.all_output)


def which(program):
    import os

    def _is_exe(filepath):
        return os.path.isfile(filepath) and os.access(filepath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if _is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if _is_exe(exe_file):
                return exe_file
    warnings.warn("{} isn't executable!\n"
                  "Make sure it is installed.".format(program))
    return None


if __name__ == "__main__":

    output_directory = "/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/LP_data",

    goku_ref_config = dict(
            adaptors='/Users/temporary/genomes/adapters/illumina_truseq.fasta',
            transcripts='/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf',
            hisat2_index='/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Sequence/HISAT2Index/genome',
            reference_genome='/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/genome.fa', )

    vegeta_config = dict(
            featurecounts='/home/pinojc/Sources/subread-1.5.1-source/bin/featureCounts',
            samstat='samstat', samtools='samtools',
            hisat2='/home/pinojc/Sources/hisat2-2.0.5/hisat2',
            flexbar='/home/pinojc/Sources/flexbar_v2.5_linux64/flexbar')

    x = RNASeqAnalyzer('test', 'test', 1, ref_config=goku_ref_config,
                       executables_config=vegeta_config, write_bash=True)
    x.gene_exp()
    print(x._bash_file)
