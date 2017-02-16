from RNASeqAnalyzer import RNASeqAnalyzer
import os
entity_searched = 'snp'  # Valid options include 'indel', 'snp', or 'none'
# INDEL/SNP Search Options - via BWA
class SNPIndel(RNASeqAnalyzer):
    def __init__(self):
        self.bwa_output = 'BWA_BAM_files'
    def bwa_alignment(self):
        print("Starting alignment for {}".format(self.sample_base))
        if not os.path.exists('BWA_BAM_files'):
            os.mkdir('BWA_BAM_files')

        path_to_executable = '{} mem'.format(bwa)
        first_pair_reads = "{}/{}-trimmed_1.{}".format(fasta_directory,
                                                       self.sample_base,
                                                       sample_suffix)
        second_pair_reads = "{}/{}-trimmed_2.{}".format(fasta_directory,
                                                        self.sample_base,
                                                        sample_suffix)
        important_options = '-T 15 -M -t {}'.format(self.n_cpu)
        path_to_reference = reference_genome
        export_to_file = '> ./BWA_BAM_files/{}.sam'.format(self.sample_base)
        command = [path_to_executable, important_options, path_to_reference,
                   first_pair_reads, second_pair_reads,
                   export_to_file]
        self._run(command)
        print("Done aligning {}".format(self.sample_base))
        # with open('{}.sam'.format(sample_base)) as f:
        #     f.write(global_output)

    def bwa_read_group(self):
        print("Starting read group addition for {}".format(self.sample_base))
        path_to_executable = 'java -jar {}'.format(picard)
        picard_program = "AddOrReplaceReadGroups"
        input_files = 'I=./BWA_BAM_files/{}.sam'.format(self.sample_base)
        output_files = 'O=./BWA_BAM_files/{}.rg.sam'.format(self.sample_base)
        necessary_parameters = "RGID={0} RGLB={0} RGPL=ILLUMINA" \
                               "RGPU=ILLUMINA RGSM={0}".format(
            self.sample_base)

        command = [path_to_executable, picard_program, input_files,
                   output_files, necessary_parameters]
        self._run(command)
        print("Done adding read group {}".format(self.sample_base))
        os.remove('BWA_BAM_files/{}.sam'.format(self.sample_base))



    def bwa_sam_to_bam(self):
        print("Starting sam to bam conversion for {}".format(self.sample_base))
        path_to_executable = '{} view'.format(samtools)
        path_to_samples = '-S -b BWA_BAM_files/{}.rg.sam'.format(
            self.sample_base)
        output_filename = '-o BWA_BAM_files/{}.bam'.format(self.sample_base)
        threads = '--threads {}'.format(n_cpus)
        command = [path_to_executable, path_to_samples, threads,
                   output_filename]
        self._run(command)
        print("Done converting sam to bam conversion "
              "for {}".format(self.sample_base))
        os.remove('BWA_BAM_files/{}.rg.sam'.format(self.sample_base))



    def bwa_index(self):
        print("Starting index for {}".format(self.sample_base))
        samtools = '/home/pinojc/RNASeq_sources/Software/./sambamba_v0.6.5'
        path_to_executable = '{} index'.format(samtools)
        path_to_samples = './BWA_BAM_files/{}.sorted.bam'.format(
            self.sample_base)
        threads = '-p --nthreads={}'.format(self.n_cpu)
        command = [path_to_executable, path_to_samples]
        command = [path_to_executable, threads, path_to_samples]
        self._run(command)
        print("Done with index for {}".format(self.sample_base))

    def bwa_samstat_analysis(self):
        print("Starting samstat analysis for {}".format(self.sample_base))

        path_to_executable = samstat
        path_to_samples = './BWA_BAM_files/{}.sorted.bam'.format(
            self.sample_base)
        command = [path_to_executable, path_to_samples]
        self._run(command)
        print("Done with samstat analysis for {}".format(self.sample_base))
        os.rename('./BWA_BAM_files/{}.sorted.bam.samstat.html'.format(
            self.sample_base),
                  './SAMSTAT_analysis/{}.bwa.sorted.bam.samstat.html'.format(
                          self.sample_base))

    # Start of steps involving selection of 'indel' or 'snp'
    def gatk_intervals(self):
        print("Starting gatk intervals for {}".format(self.sample_base))
        path_to_executable = "java -jar {}".format(gatk)
        gatk_program = '-T RealignerTargetCreator'
        path_to_reference = "-R {}".format(reference_genome)
        nt = '-nt {}'.format(self.n_cpu)
        input_files = '-I ./BWA_BAM_files/{}.sorted.bam'.format(
            self.sample_base)

        if entity_searched == 'indel':
            output_file = '-o ./BWA_BAM_files/{}.indel.intervals'.format(
                    self.sample_base)
            path_to_vcf = "--known {}".format(indel_vcf_file)
        elif entity_searched == 'snp':
            output_file = '-o ./BWA_BAM_files/{}.snp.intervals'.format(
                self.sample_base)
            path_to_vcf = "--known {}".format(snp_vcf_file)

        command = [path_to_executable, gatk_program,nt,  path_to_reference,
                   input_files, path_to_vcf, output_file]
        self._run(command)
        print("Done with gatk intervals for {}".format(self.sample_base))

    def gatk_realignment(self):
        print("Starting gatk realignment for {}".format(self.sample_base))
        path_to_executable = "java -jar {}".format(gatk)
        gatk_program = '-T IndelRealigner'
        # nt = '-nt {}'.format(self.n_cpu)
        path_to_reference = "-R {}".format(reference_genome)
        options = '--maxReadsForRealignment 999999 --maxReadsInMemory 999999'

        if entity_searched == 'indel':
            input_files = '-I ./BWA_BAM_files/{}.sorted.bam'.format(
                self.sample_base)
            intervals = '-targetIntervals ./BWA_BAM_files/{}.indel.intervals'.format(
                    self.sample_base)
            output_file = '-o ./BWA_BAM_files/{}.indel.realigned.bam'.format(
                    self.sample_base)
        elif entity_searched == 'snp':
            input_files = '-I ./BWA_BAM_files/{}.sorted.bam'.format(
                self.sample_base)
            intervals = '-targetIntervals ./BWA_BAM_files/{}.snp.intervals'.format(
                    self.sample_base)
            output_file = '-o ./BWA_BAM_files/{}.snp.realigned.bam'.format(
                    self.sample_base)

        command = [path_to_executable, gatk_program,  path_to_reference,
                   input_files, intervals, options, output_file]
        self._run(command)
        print("Done with gatk realignment for {}".format(self.sample_base))

    def gatk_recalibration(self):
        print(
        "Starting gatk indel recalibration for {}".format(self.sample_base))
        path_to_executable = "java -jar {}".format(gatk)
        gatk_program = '-T BaseRecalibrator'
        path_to_reference = "-R {}".format(reference_genome)
        options = '-l INFO'

        if entity_searched == 'indel':
            input_files = '-I ./BWA_BAM_files/{}.indel.realigned.bam'.format(
                    self.sample_base)
            output_file = '-o ./BWA_BAM_files/{}.indel.recal.table'.format(
                    self.sample_base)
            path_to_vcf = "-knownSites {}".format(indel_vcf_file)
        elif entity_searched == 'snp':
            input_files = '-I BWA_BAM_files/{}.snp.realigned.bam'.format(
                    self.sample_base)
            output_file = '-o ./BWA_BAM_files/{}.snp.recal.table'.format(
                    self.sample_base)
            path_to_vcf = "-knownSites {}".format(snp_vcf_file)

        command = [path_to_executable, gatk_program, path_to_reference,
                   input_files, options, path_to_vcf, output_file]
        self._run(command)
        print(
        "Done with gatk indel recalibration for {}".format(self.sample_base))

    def gatk_realign_recal(self):
        print("Starting gatk realign recal for {}".format(self.sample_base))
        path_to_executable = "java -jar {}".format(gatk)
        gatk_program = '-T PrintReads'
        path_to_reference = "-R {}".format(reference_genome)

        if entity_searched == 'indel':
            input_files = '-I ./BWA_BAM_files/{}.indel.realigned.bam'.format(
                    self.sample_base)
            options = '-BQSR BWA_BAM_files/{}.indel.recal.table'.format(
                    self.sample_base)
            output_file = '-o ./BWA_BAM_files/{}.indel.realigned.recal.bam'.format(
                    self.sample_base)
        elif entity_searched == 'snp':
            input_files = '-I ./BWA_BAM_files/{}.snp.realigned.bam'.format(
                    self.sample_base)
            options = '-BQSR ./BWA_BAM_files/{}.snp.recal.table'.format(
                    self.sample_base)
            output_file = '-o ./BWA_BAM_files/{}.snp.realigned.recal.bam'.format(
                    self.sample_base)

        command = [path_to_executable, gatk_program, path_to_reference,
                   input_files, options, output_file]
        self._run(command)
        print("Done with gatk realign recal for {}".format(self.sample_base))
        # os.remove('./BWA_BAM_files/{}.indel.realigned.bam'.format(sample_base))

    def mark_dup(self):
        print("Starting mark dup for {}".format(self.sample_base))
        path_to_executable = "java -jar {}".format(picard)
        picard_program = 'MarkDuplicates'

        if entity_searched == 'indel':
            input_files = 'I=./BWA_BAM_files/{}.indel.realigned.recal.bam'.format(
                    self.sample_base)
            output_file = 'O=./BWA_BAM_files/{}.indel.realigned.recal.dupmarked.bam'.format(
                    self.sample_base)
            options = 'M=./BWA_BAM_files/{}-indel-marked_dup_metrics.txt'.format(
                    self.sample_base)
        elif entity_searched == 'snp':
            input_files = 'I=./BWA_BAM_files/{}.snp.realigned.recal.bam'.format(
                    self.sample_base)
            output_file = 'O=./BWA_BAM_files/{}.snp.realigned.recal.dupmarked.bam'.format(
                    self.sample_base)
            options = 'M=./BWA_BAM_files/{}-snp-marked_dup_metrics.txt'.format(
                    self.sample_base)

        command = [path_to_executable, picard_program, input_files,
                   output_file, options]
        self._run(command)
        print("Done with mark dup for {}".format(self.sample_base))

    def dup_index(self):
        print("Starting indel dup index for {}".format(self.sample_base))
        path_to_executable = '{} index'.format(samtools)

        if entity_searched == 'indel':
            path_to_samples = 'BWA_BAM_files/{}.indel.realigned.recal.dupmarked.bam'.format(
                    self.sample_base)
        elif entity_searched == 'snp':
            path_to_samples = 'BWA_BAM_files/{}.snp.realigned.recal.dupmarked.bam'.format(
                    self.sample_base)

        command = [path_to_executable, path_to_samples]
        self._run(command)
        print("Done with dup index for {}".format(self.sample_base))

    def final_entity_search(self):
        print("Starting final search for {}".format(self.sample_base))
        path_to_executable = "java -jar {}".format(gatk)
        gatk_program = '-T UnifiedGenotyper -l INFO'
        path_to_reference = "-R {}".format(reference_genome)
        options = '-A Coverage -A AlleleBalance -G Standard ' \
                  '-stand_call_conf 50.0 -stand_emit_conf 10.0 ' \
                  '-mbq 20 -deletions 0.05 -dcov 1000'

        if entity_searched == 'indel':
            search_item = '-glm INDEL'
            input_files = '-I ./BWA_BAM_files/{}.indel.realigned.recal.dupmarked.bam'.format(
                    self.sample_base)
            path_to_vcf = "-D {}".format(indel_vcf_file)
            output_file = '--out {0}.indel.vcf -metrics {0}.indel.outmetrics.txt'.format(
                    self.sample_base, )
        elif entity_searched == 'snp':
            search_item = '-glm SNP'
            input_files = '-I ./BWA_BAM_files/{}.snp.realigned.recal.dupmarked.bam'.format(
                    self.sample_base)
            path_to_vcf = "-D {}".format(snp_vcf_file)
            output_file = '--out {0}.snps.vcf -metrics {0}.snp.outmetrics.txt'.format(
                    self.sample_base)

        command = [path_to_executable, gatk_program, path_to_reference,
                   input_files, path_to_vcf, output_file, options,
                   search_item]
        self._run(command)
        print("Done with final search for {}".format(self.sample_base))

        def setup_bwa(self):
            # self.bwa_alignment()
            # self.bwa_read_group()
            # self.bwa_sam_to_bam()
            # self.bam_sort('BWA_BAM_files')
            # self.bwa_samstat_analysis()
            self.bwa_index()

        def bwa_genome_search(self):
            # if cont:
            # self.setup_bwa()
            # self.gatk_intervals()
            # self.gatk_realignment()
            self.gatk_recalibration()
            self.gatk_realign_recal()
            self.mark_dup()
            self.dup_index()
            self.final_entity_search()
            # if entity_searched == 'indel':
            #     os.remove('./BWA_BAM_files/{}.indel.realigned.recal.bam'.format(sample_base))
            # elif entity_searched == 'snp':
            #     os.remove('./BWA_BAM_files/{}.snp.realigned.recal.bam'.format(sample_base))
            return False