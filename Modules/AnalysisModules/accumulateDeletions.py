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
VEGFA_guide = 'GGTGAGTGAGTGTGTGCGTGTGG'

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

def sw_local_alignment(read, reference_seq):
    local_align_result = pairwise2.align.localms(read, reference_seq, 10, -4, -10, -1,
                                                 one_alignment_only=True)  # match, mismatch, open, extend (10, -4, -6, -1)
    align_score = local_align_result[0][2]
    # align1, align2, score, local_align_start, local_align_end = local_align_result[0]
    print(format_alignment(*local_align_result[0]))
    return local_align_result, align_score


def get_cigar_operation(dict, op_id):
    return dict[op_id]


def get_deletions_from_cigar(cigar, min_del_len, max_del_len):
    del_list = []
    for tup in cigar:
        if tup[0] != 2:
            continue
        if min_del_len < tup[1] < max_del_len:
            del_list.append(tup[1])
    return del_list


os.chdir('../../')
root_dir = os.getcwd() + '/Data/SRA_SRR1561378/mouse_genome_data/'
print(root_dir)

sys.stdout = open(root_dir + 'res_new_sorted_del_analysis.txt', 'w')

bam_path = root_dir + 'res_new_sorted.bam'
bamfile = pysam.AlignmentFile(bam_path, "rb", ignore_truncation=True)

pairs = [(0, 'MATCH'), (1, 'INS'), (2, 'DEL'), (3, 'REF_SKIP'), (4, 'SOFT_CLIP'), (5, 'HARD_CLIP'), (6, 'PAD'),
         (7, 'EQUAL'), (8, 'DIFF'), (9, 'BACK')]
bam_cigar_dict = dict(pairs)

del_dict = {}
mapq_th = 30
min_del = 30
max_del = 90
stop = 10000
i = 0
for b in bamfile.fetch(until_eof=True):
    # if i == stop:
    #     break
    if i % 100000 == 0:
        print("processing read {}".format(i))
    i += 1

    d = 0
    # check mapq
    if b.mapping_quality >= mapq_th:
        # get list of del size
        del_list = get_deletions_from_cigar(b.cigar, min_del, max_del)
        if del_list:
            # take max_del in case more than 1
            d = max(del_list)
            if (b.reference_name, d) in del_dict.keys():
                del_dict[b.reference_name, d].append(b)
            else:
                del_dict[b.reference_name, d] = []
                del_dict[b.reference_name, d].append(b)

print('------------Summary-----------')
print('Processed {}'.format(bam_path))
print('Results: \n')

# local align thresholds
genome = Fasta(genome_path, sequence_always_upper=True)
num_of_read_th = 4
val_range_th = 200
expand = 500

for dd in del_dict.items():
    val_ranges = []
    print('{},deletion_length_{}'.format(dd[0][0], dd[0][1]))
    for r in dd[1]:
        start = r.reference_start
        print(start)
        val_ranges.append(start)
    max_val = max(val_ranges)
    min_val = min(val_ranges)
    val_range = max_val - min_val
    dd[1].append(val_range)
    dd[1].append(len(dd[1]))
    print("value_range = {}, num_of_reads = {}".format(val_range, len(dd[1])-2))

    # Perform local align
    if dd[1][-2] < val_range_th and dd[1][-1] >= num_of_read_th:
        for r in dd[1]:
            try:
                start = r.reference_start
                reference_seq = genome[r.reference_name][start - expand:start + expand].seq
                local_align_result, align_score = sw_local_alignment(r.seq, reference_seq)
                # VEGFA guide alignment
                # print('-----VEGFA guide alignment to genome-----')
                # local_align_result, align_score = sw_local_alignment(VEGFA_guide, reference_seq)
                # print('----- VEGFA guide alignment to read -----')
                # local_align_result, align_score = sw_local_alignment(VEGFA_guide, r.seq)
            except:
                break

sys.stdout.close()
