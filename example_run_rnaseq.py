import subprocess
import os

# setting up an output directory and file saving

run_name = 'run_sample_xyz'
output_directory = 'sample1/{}'.format(run_name)
read_type = 'SE'
n_cpus = 16

# paths to programs used
gtk_path = 'GenomeAnalysisTK.jar'
flexbar_path = "/usr/local/bin/flexbar_v2.5_macosx/flexbar"
fastqc_path = "/usr/local/bin/FastQC/fastqc"
star_path = "/home/pinojc/Sources/STAR-STAR_2.4.2a/bin/Linux_x86_64/STAR"
samstat = "/usr/local/bin/samstat"
cufflinks = "/home/pinojc/Sources/cufflinks-2.2.1.Linux_x86_64/cufflinks"
feature_counts = "/usr/local/bin/featureCounts"
samtools = "/home/pinojc/Sources/samtools-1.3.1/samtools"

# static variables
# set your fasta files here
genome_fasta_path = '/Users/temporary/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa'
index_path = '/home/pinojc/RNA_Seq/STARIndex_mm10_vM10'
fastq_path = '/home/pinojc/TEST_RNA_SEQ/fastq'
adapters = "/Users/temporary/genomes/adapters/illumina_truseq.fasta"
transcripts = "/home/pinojc/RNA_Seq/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"
star_index = "/home/pinojc/RNA_Seq/STARIndex_mm10_vM10"
ref_genome = "/home/pinojc/RNA_Seq/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
# you do these for all the files that won't change.




if not os.path.exists(output_directory):
    os.makedirs(output_directory)
os.chdir(output_directory)

def run_star(sample):
    local_out_dir = os.path.join('alignments','star',sample)
    if not os.path.exists(local_out_dir):
        os.makedirs(local_out_dir)
    num_cores = '--runThreadN {}'.format(n_cpus)
    genome_path = '--genomeDir {}'.format(star_index)
    if read_type == 'SE':
        read_files = '--readFilesIn {0}/{1}_1.fastq'.format(fastq_path, sample)
    elif read_type == 'PE':
        read_files ='{0}/{1}_1.fastq {0}/{1}_2.fastq'.format(fastq_path, sample)
    else:
        print('Must provide a read type (SE of PE)')
        read_files =None
        quit()
    sjdb_gtf_files = '--sjdbGTFfile {}'.format(transcripts)
    output_name = '--outFileNamePrefix {}/{}'.format(local_out_dir,sample)
    # .. continue all the commands you would use for each
    bam_output = '--outSAMtype BAM SortedByCoordinate'
    remove_noncan = '--outFilterIntronMotifs RemoveNoncanonical'
    working_d = os.getcwd()
    print(working_d)
    command = [star_path,num_cores, genome_path,bam_output,remove_noncan,
               read_files,sjdb_gtf_files, output_name ]
    print(' '.join(i for i in command))
    p = subprocess.Popen(command,)
                         #stdout=subprocess.STDOUT,
                         #stderr=subprocess.STDOUT)


def run_gtk():
    # This {} symbol replaces what you have in the .format() command
    path_to_exectuable = 'java -jar {} -T RealignerTargetCreator'.format(gtk_path)
    r_option_path = '-R {}'.format(genome_fasta_path)
    # This assumes the bam files are placed into a directory inside of the output
    bam_location = '-I {}/bam_files/{}.bam'.format(output_directory, run_name)
    output_name = '{}_gtk_output.intervals'
    vcf_location = '--known /Users/temporary/new.vcf'
    command = [path_to_exectuable, r_option_path, vcf_location,bam_location,
               output_name,
               ]
    print(command)
    p = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)



run_star('sample-84')

# run_gtk()