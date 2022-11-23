import ast

import khmer
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
from pyfaidx import Fasta

sns.set()
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import itertools
import collections

from Modules import kmer_operations
from datetime import datetime

def generateKmerPermutationsDictionary(K):
  alphabet = 'ACGTN'
  permutations_dict = {}

  for idx, output in enumerate(itertools.product(alphabet, repeat=K)):
      permutation = ''.join(output)
      permutations_dict[permutation] = format(idx, '019b')

  return permutations_dict


def kmers_binary_representation(read,region):
    K = 8
    permutation_dict = generateKmerPermutationsDictionary(K)
    alphabet = ['A', 'C', 'G', 'T', 'N']
    cg = khmer.Countgraph(K, len(alphabet) ** K, 1)
    read_kmers = cg.get_kmers(read)
    read_bin_kmers = [permutation_dict[x] for x in read_kmers]
    region_kmers = cg.get_kmers(region)
    region_bin_kmers = [permutation_dict[x] for x in region_kmers]

    region_kmers_counter = collections.Counter()
    for kmer in region_bin_kmers:
        region_kmers_counter[kmer] += 1

    read_kmers_counter = collections.Counter()
    for kmer in read_bin_kmers:
        read_kmers_counter[kmer] += 1

    plot_kmer_binary_representation_dict(region_kmers_counter)
    # plot_kmer_binary_representation_dict(read_kmers_counter)

    ################################################
    region_kmers_full = permutation_dict.copy()
    region_kmers_full.update({}.fromkeys(region_kmers_full,0))
    for kmer in region_bin_kmers:
        region_kmers_full[kmer]=1

    # plot_kmer_binary_representation_dict(region_kmers_full)
    return

def plot_kmer_binary_representation_dict(kmers_counter):
    plt.figure(figsize=(6, 24))
    Y = list(kmers_counter.keys())
    X = list(kmers_counter.values())
    plt.scatter(X, Y, alpha=0.7, c=X, cmap='winter')
    plt.title("Sequence K-Mer Expression")
    plt.xlabel("K-Mer Counts")
    plt.ylabel("K-Mer Binary Representation ")
    plt.xticks(np.arange(min(X), max(X) + 1, 1))
    plt.show()


def read_kmer_analysis(read):
    regions_kmers = []
    with open('Data/regionsKmersRaw_k8.txt', 'r') as f:
        for region_idx, region in enumerate(f.readlines()):
            region_kmers = ast.literal_eval(region)
            regions_kmers.append(set(region_kmers))
    alphabet = ['A', 'C', 'G', 'T', 'N']
    K= 8
    cg = khmer.Countgraph(K, len(alphabet) ** K, 1)
    jac_sim_res = np.zeros(len(regions_kmers))
    l = set(cg.get_kmers(read))
    for region_idx, region in enumerate(regions_kmers):
        intersection = len(l.intersection(region))
        union = len(l.union(region))
        if union == 0:
            print('Divide by zero')
            continue
        # print("Intersection = {} Union = {} #kmers = {} ".format(intersection,union,len(region)))
        jac_sim_res[region_idx] = intersection / union
    return jac_sim_res

def generate_genome_object(genome_path):
    genome = Fasta(genome_path, sequence_always_upper=True)
    return genome

def check_unmapped(bam_path):
    bamfile = pysam.AlignmentFile(bam_path, "rb", ignore_truncation=True)
    stat_list = bamfile.get_index_statistics()
    # 0-Chr / 1-mapped / 2-unmapped / 3-total
    # Check unmapped reads
    for item in stat_list:
        if item[2] > 0:
            print(item)

def comapre_alignments(bam_path1,bam_path2):
    bamfile1 = pysam.AlignmentFile(bam_path1, "rb", ignore_truncation=True)
    bamfile2 = pysam.AlignmentFile(bam_path2, "rb", ignore_truncation=True)

    for read in bamfile1.fetch(until_eof=True):
        print(read.query_name)
        break

def check_overlap(read_start,read_end,regions_df):
    # Assuming no overlaps between regions
    min_overlap = 1
    region_index_low = regions_df.loc[
        ((read_end - regions_df['chromStart']) > min_overlap) & (read_start < regions_df['chromStart'])
        ].index.tolist()

    region_index_high = regions_df.loc[
        ((regions_df['chromEnd'] - read_start) > min_overlap) & (read_end > regions_df['chromEnd'])
        ].index.tolist()

    if (not region_index_low) & (not region_index_high):
        return [0,-1]

    else:
        if region_index_low:
            region_start = regions_df.loc[region_index_low[0]].at['chromStart']
            ov = read_end - region_start
            return ov.item(), region_index_low[0]
        if region_index_high:
            region_end = regions_df.loc[region_index_high[0]].at['chromEnd']
            ov = region_end - read_start
            return ov.item(), region_index_high[0]

def plot_kmer_read_cap_heatmap(similarity_mtrx):
    plt.figure(figsize=(18, 8))
    print('Plotting Read to Region Similarity')
    X = list(range(len(similarity_mtrx)))
    ax = sns.histplot(x=X,y=similarity_mtrx,bins=len(similarity_mtrx))
    print(similarity_mtrx.shape)
    plt.xlabel('Region Index', fontsize=18)
    plt.ylabel('Similarity', fontsize=18)
    plt.title('Similarity Score For Each Region', fontsize=20)
    plt.xticks(X)
    plt.xticks(rotation=45)
    plt.show()

def plot_similarity_histogram(similarity_mtrx,width):

    plt.figure(figsize=(18, 8))
    X = list(range(len(similarity_mtrx)))
    plt.bar(X, similarity_mtrx, width, color='b',alpha=0.7)
    plt.title("Jaccard Similarity")
    plt.xticks(rotation=90)
    plt.xlabel("Regions")
    plt.ylabel("Jaccard Similarity")
    # plt.xticks(np.arange(min(X), max(X) + 1, 4))
    plt.show()

def local_alignment(read,region):
    local_align_result = pairwise2.align.localms(read, region, 10, -10, -10, -9, one_alignment_only=True)
    align_score = local_align_result[0][2]
    return local_align_result,align_score


def generate_regions_sequences_list(path):
    regions_raw = []
    with open(path, 'r') as f:
        for region_idx, region in enumerate(f.readlines()):
            regions_raw.append(region)
    return regions_raw

def generate_regions_kmers_list(path):
    regions_kmers_list = []
    with open(path, 'r') as f:
        for region_idx, region in enumerate(f.readlines()):
            region_kmers = ast.literal_eval(region)
            regions_kmers_list.append(region_kmers)
    return regions_kmers_list


def calculate_regions_average_length(regions_list):
    len_sum = 0
    for i,region in enumerate(regions_list):
        len_sum += len(region)

    total_avg = len_sum/(i+1)
    print("Average length of region = {}".format(total_avg))
    return total_avg


def generate_bed_df(bed_path):
    df_bed = pd.read_csv(bed_path, sep='\t', comment='t', header=None)
    header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
    df_bed.columns = header[:len(df_bed.columns)]
    return df_bed


def unmatched_reads_meta_df(reads_meta_df_path):
    df = pd.read_pickle(reads_meta_df_path)
    df = df[df['label_region'] != -1]
    df_unmatched = df[df['label_region'] != df['predict_region']]
    predictions = df_unmatched['predict_confidence'].tolist()
    return df_unmatched,predictions


def extract_x_maximums_from_list(l,num_of_max):
    ind_max = l.argsort()[-num_of_max:][::-1]
    print('{} maximal values are {} at locations {}'.format(num_of_max,l[ind_max], ind_max))
    return ind_max, l[ind_max]


def two_step_alignment(reads_meta_df,bamfile,regions_raw_list):
    # Random choice
    # random state allow keeping the random choice
    matched_df = reads_meta_df[reads_meta_df['predict_region'] == reads_meta_df['label_region']]
    matched_df_avg_high = matched_df[matched_df['predict_confidence'] > 0.05]
    unmatched_df = reads_meta_df[reads_meta_df['predict_region'] != reads_meta_df['label_region']]

    # Sample from max
    # ind_max, val = extract_x_maximums_from_list(matched_df_avg_high['predict_confidence'].to_numpy(), 1)
    ind_max, val = extract_x_maximums_from_list(reads_meta_df['predict_confidence_second_max'].to_numpy(), 3)
    read_sampled = reads_meta_df.iloc[ind_max[2]]
    # Sample from random choice of matched
    # read_sampled = matched_df_avg_high.sample(n=1, random_state=1)

    # read_chr = read_sampled.chr.item()
    # read_start = read_sampled.start.item()
    # read_end = read_sampled.end.item()

    read_chr = read_sampled.chr
    read_start = read_sampled.start
    read_end = read_sampled.end

    for read in bamfile.fetch(read_chr,read_start,read_end):
        similarity_array = read_kmer_analysis(read.seq)
        predicted_region = np.argmax(similarity_array)
        print("Predicted region={}".format(predicted_region))
        extract_x_maximums_from_list(similarity_array, 1)
        # plot_similarity_histogram(similarity_array, 1)
        local_alignment_region = predicted_region
        local_align_result, align_score = local_alignment(read.seq, regions_raw_list[local_alignment_region])
        print(format_alignment(*local_align_result[0]))
        kmers_binary_representation(read.seq, regions_raw_list[local_alignment_region].rstrip("\n"))

        break

    return


def kmer_local_alignment(df,bed_df,regions_kmers_list,k):
    ko = kmer_operations.KmerOperationsClass()
    alphabet = ['A', 'C', 'G', 'T', 'N']
    cg = khmer.Countgraph(k, len(alphabet) ** k, 1)

    df['first_match_kmer_pos'] = -1
    df['updated_start'] = -1
    df['start_gap'] = -1
    df['last_match_kmer_pos'] = -1
    df['updated_end'] = -1
    df['end_gap'] = -1
    df['sequence_length_gap'] = -1

    for index, read in df.iterrows():
        read_kmers = cg.get_kmers(read.sequence)
        region_kmers = regions_kmers_list[read['predict_region']]

        [element_first,indx_first,element_last,indx_last] = ko.kmer_local_alignment2(read_kmers,region_kmers)
        if (indx_first != None) and (indx_last != None) :
            updated_start = bed_df.loc[read['predict_region']]['chromStart'] + indx_first
            updated_end = bed_df.loc[read['predict_region']]['chromStart'] + indx_last
            df.at[index, 'first_match_kmer_pos'] = indx_first
            df.at[index, 'updated_start'] = updated_start
            df.at[index, 'start_gap'] = abs(read.start - updated_start)
            df.at[index, 'last_match_kmer_pos'] = indx_last
            df.at[index, 'updated_end'] = updated_end
            df.at[index, 'end_gap'] = abs(read.end - updated_end)
            df.at[index, 'sequence_length_gap'] = abs((updated_end - updated_start) - (read.end - read.start))

    return df


def main():

    # Paths
    bam_path1 = 'Data/Rag1_chr11_0_10000000.bam'
    bamfile1 = pysam.AlignmentFile(bam_path1, "rb", ignore_truncation=True)
    bam_path2 = 'Data/bwa_virtualReferenceGenome_Rag1_chr11_0_10000000.sorted.sorted.bam'
    bamfile2 = pysam.AlignmentFile(bam_path2, "rb", ignore_truncation=True)
    bam_path3 = 'Data/Rag1.sorted.dedup.bam'
    bamfile3 = pysam.AlignmentFile(bam_path3, "rb", ignore_truncation=True)

    regions_kmers_path = 'Data/regionsKmersRaw_k8.txt'
    raw_regions_path = 'Data/regionsRawSequence_k8.txt'
    bedfile_path = 'Data/RAG1.bed'
    reads_meta_df_path = 'Data/reads_meta_df_29-08-2021_18-28.pkl'
    reads_meta_df = pd.read_pickle(reads_meta_df_path)

    df_unmatched = unmatched_reads_meta_df(reads_meta_df_path)

    df_bed = generate_bed_df(bedfile_path)

    regions_raw = generate_regions_sequences_list(raw_regions_path)
    regions_kmers_list = generate_regions_kmers_list(regions_kmers_path)

    average_region_length = calculate_regions_average_length(regions_raw)

    # two_step_alignment(reads_meta_df,bamfile3,regions_raw)
    K=8
    reads_meta_df_roi = reads_meta_df[reads_meta_df['label_region'] == reads_meta_df['predict_region']]
    print(len(reads_meta_df_roi))
    updated_df = kmer_local_alignment(reads_meta_df_roi,df_bed,regions_kmers_list,K)
    updated_df.to_pickle(reads_meta_df_path.strip('.pkl') + '_updated.pkl')
    return



if __name__ == '__main__':
    main()

# for read in bamfile3.fetch('chr11'):
#
#     print("Ground Truth region = {}".format(unmatched.loc[59].at['label_region']))
#     similarity_array = read_kmer_analysis(read.seq)
#     predicted_region = np.argmax(similarity_array)
#     ind_max = similarity_array.argsort()[-num_of_max:][::-1]
#     print(similarity_array[ind_max])
#     plot_similarity_histogram(similarity_array, 1)
#     print("Predicted region={}".format(predicted_region))
#     local_alignment_region = predicted_region #61
#     local_align_result, align_score = local_alignment(read.seq, regions_raw[local_alignment_region])
#     print(format_alignment(*local_align_result[0]))
#     print("DONE!")
#
# check_unmapped(bam_path1)
# bamfile = pysam.AlignmentFile(bam_path1, "rb", ignore_truncation=True)
# for read in bamfile.fetch('chr11'):
#     if read.is_unmapped:
#         # overlap,region = check_overlap(read.reference_start,read.reference_end,df_bed)
#         # print("Correct region = {}".format(region))
#         similarity_array = read_kmer_analysis(read.seq)
#         predicted_region = np.argmax(similarity_array)
#         if predicted_region == 341:
#             plot_similarity_histogram(similarity_array,1)
#             print("Predicted region={}".format(predicted_region))
#             local_align_result,align_score = local_alignment(read.seq,regions_raw[predicted_region])
#             print(format_alignment(*local_align_result[0]))
#

