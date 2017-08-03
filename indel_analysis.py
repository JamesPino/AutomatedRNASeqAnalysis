import pybedtools

WT3 = pybedtools.BedTool(
    '/Users/lisapoole/Desktop/65_RNAseq_on_SMARCAL1-knockout_human_cells/attachments/RNA_analysis/genome_analysis_bwa/final_indel_files/3612-DC-1.indel.vcf')
WT4 = pybedtools.BedTool(
    '/Users/lisapoole/Desktop/65_RNAseq_on_SMARCAL1-knockout_human_cells/attachments/RNA_analysis/genome_analysis_bwa/final_indel_files/3612-DC-2.indel.vcf')
KO28 = pybedtools.BedTool(
    '/Users/lisapoole/Desktop/65_RNAseq_on_SMARCAL1-knockout_human_cells/attachments/RNA_analysis/genome_analysis_bwa/final_indel_files/3612-DC-3.indel.vcf')
KO66 = pybedtools.BedTool(
    '/Users/lisapoole/Desktop/65_RNAseq_on_SMARCAL1-knockout_human_cells/attachments/RNA_analysis/genome_analysis_bwa/final_indel_files/3612-DC-4.indel.vcf')
WT = pybedtools.BedTool('/Users/lisapoole/Desktop/intersection-of-WT3-WT4.vcf')
KO = pybedtools.BedTool('/Users/lisapoole/Desktop/intersection-of-KO28-KO66.vcf')

WT3_and_WT4 = WT3.intersect(WT4)
KO28_and_KO66 = KO28.intersect(KO66)
all = WT.intersect(KO)

# e = WT3_and_WT4.saveas('intersection-of-WT3-WT4.vcf', trackline='track name="WT3 and WT4"')
e = all.saveas('intersection-of-all.vcf', trackline='track name="WT and KO"')
