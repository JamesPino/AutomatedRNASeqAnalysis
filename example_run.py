from RNASeq_pipeline import RNASeq_pipeline

output_directory = "/media/pinojc/68f7ba6a-cdf6-4761-930d-9c1bb724e40d/home/LP_data/"
rna_seq = RNASeq_pipeline('3612-DC-2', n_cpu=8, output_directory=output_directory)
rna_seq.bwa_genome_search()
