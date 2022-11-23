import pysam
import os
import sys
import pandas as pd
from Modules import analysis_module
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

sns.set_theme(style="darkgrid")
from pyfaidx import Fasta

# genome_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/hg38.fa'
genome_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/mm9_mouse_genome.fa'
gRNA1 = 'GCTCCCTACGCATGCGTCCCAGG'
gRNA2 = 'AATGGCTCAGGTTTGTCGCGCGG'

def sw_local_alignment(read, reference_seq,print_score_th):
    local_align_result = pairwise2.align.localms(read, reference_seq, 10, -4, -6, -1,
                                                 one_alignment_only=True)  # match, mismatch, open, extend (10, -4, -6, -1)
    align_score = local_align_result[0][2]
    # align1, align2, score, local_align_start, local_align_end = local_align_result[0]
    if align_score > print_score_th:
        print(format_alignment(*local_align_result[0]))
    return local_align_result, align_score

os.chdir('../../')
root_dir = os.getcwd() + '/Data/SRA_SRR1561378/mouse_genome_data/'
print(root_dir)

sys.stdout = open(root_dir + 'ogeen_interesting_deletions_to_align_guide_res2.txt', 'w')
# local align thresholds
genome = Fasta(genome_path, sequence_always_upper=True)

locations_to_align = pd.read_csv(root_dir + 'ogeen_interesting_deletions_to_align_guide.txt', sep='\t')
locations_to_align.columns = ['chr','start','end']
# regions = [('chr11',3092029,3092129),('chr11',70458457,70458557),('chr12',56593531,56593631),('chr13',90970620,90970720),('chr7',78188120,78188220),('chr7',151470810,151470910)]
regions = []
for index, row in locations_to_align.iterrows():
    location = (row.chr,row.start,row.end)
    regions.append(location)
    print(location)


# print(regions[0][0])
print_score_th = 180
for region in regions:

    reference_seq = genome[region[0]][region[1]:region[2]].seq
    print('----- gRNA1 for region {} -----'.format(region))
    local_align_result, align_score = sw_local_alignment(gRNA1, reference_seq,print_score_th)
    print('----- gRNA2 for region {} -----'.format(region))
    local_align_result, align_score = sw_local_alignment(gRNA2, reference_seq,print_score_th)
#
#
sys.stdout.close()
