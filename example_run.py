from RNASeq_pipeline import RNASeq_pipeline

rna_seq = RNASeq_pipeline('3612-DC-1', n_cpu=4)
rna_seq.hisat2_alignment()
rna_seq.write_log_file('output_run_1.txt')
