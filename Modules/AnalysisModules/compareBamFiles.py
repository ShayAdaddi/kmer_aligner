import pysam
import os
import pandas as pd
from Modules import analysis_module
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_theme(style="darkgrid")
import collections

# M	BAM_CMATCH	0
# I	BAM_CINS	1
# D	BAM_CDEL	2
# N	BAM_CREF_SKIP	3
# S	BAM_CSOFT_CLIP	4
# H	BAM_CHARD_CLIP	5
# P	BAM_CPAD	6
# =	BAM_CEQUAL	7
# X	BAM_CDIFF	8
# B	BAM_CBACK	9

def get_cigar_operation(dict,op_id):
    return dict[op_id]

def cigar_counter(cigar,counter):
    for tup in cigar:
        if tup[0] == 2:
            counter.update({'reads_del': 1})
            counter.update({'{}'.format(tup[1]): 1})
        if tup[0] == 4:
            counter.update({'reads_soft_clip': 1})
            counter.update({'bp_accum_soft_clip': tup[1]})
        if tup[0] == 5:
            counter.update({'reads_hard_clip': 1})
            counter.update({'bp_accum_hard_clip': tup[1]})

os.chdir('../../')
root_dir = os.getcwd() + '/Data/SRA_SRR1561378/mouse_genome_data/'
print(root_dir)

# bam1_path = root_dir + 'RAG1_TX_sorted.bam'
# bam2_path = root_dir + 'RAG1_TX_new_sorted.bam'

bam1_path = root_dir + 'res_sorted.bam'
bam2_path = root_dir + 'res_new_sorted.bam'
# bam3_path = root_dir + '/virtual_genome/res_vg.bam'
# bam4_path = root_dir + '/virtual_genome/res_vg_new.bam'

bamfile1 = pysam.AlignmentFile(bam1_path, "rb", ignore_truncation=True)
bamfile2 = pysam.AlignmentFile(bam2_path, "rb", ignore_truncation=True)
# bamfile3 = pysam.AlignmentFile(bam2_path, "rb", ignore_truncation=True)
# bamfile4 = pysam.AlignmentFile(bam2_path, "rb", ignore_truncation=True)

pairs = [(0,'MATCH'),(1,'INS'),(2,'DEL'),(3,'REF_SKIP'),(4,'SOFT_CLIP'),(5,'HARD_CLIP'),(6,'PAD'),(7,'EQUAL'),(8,'DIFF'),(9,'BACK')]
bam_cigar_dict = dict(pairs)

c1 = collections.Counter()
c2 = collections.Counter()
# c3 = collections.Counter()
# c4 = collections.Counter()

# for b1,b2,b3,b4 in zip(bamfile1.fetch(until_eof=True),bamfile2.fetch(until_eof=True),bamfile3.fetch(until_eof=True),bamfile4.fetch(until_eof=True)):
for b1,b2 in zip(bamfile1.fetch(until_eof=True),bamfile2.fetch(until_eof=True)):
    cigar_counter(b1.cigar,c1)
    cigar_counter(b2.cigar,c2)
    # cigar_counter(b3.cigar,c3)
    # cigar_counter(b4.cigar,c4)

print('------------Summary-----------')
print('files:')
print(bam1_path)

bam1_df = pd.DataFrame.from_dict(c1.items())
bam1_df.to_excel(root_dir + "compareBamFiles_res.xlsx")
bam2_df = pd.DataFrame.from_dict(c2.items())
bam2_df.to_excel(root_dir + "compareBamFiles_res_new.xlsx")
# bam3_df = pd.DataFrame.from_dict(c3.items())
# bam3_df.to_excel(root_dir + "res_vg.xlsx")
# bam4_df = pd.DataFrame.from_dict(c4.items())
# bam4_df.to_excel(root_dir + "res_new_vg.xlsx")
print("done")




