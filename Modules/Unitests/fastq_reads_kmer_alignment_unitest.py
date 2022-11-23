from Modules import kmer_operations, services
import kmer_aligner
import pandas as pd
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import os
import sys

os.chdir('../../')
root_dir = os.getcwd() + '/Data/'
print(root_dir)
data_set_dir = 'data_set_03-12-2021_17-36/'
if os.path.exists(root_dir + data_set_dir + 'kmer_output.txt'):
    os.remove(root_dir + data_set_dir + 'kmer_output.txt')

# sys.stdout = open(root_dir + data_set_dir + 'kmer_output.txt','w')


def sw_local_alignment(read, region):
    local_align_result = pairwise2.align.localms(read, region, 10, -10, -10, -1/10000000, one_alignment_only=True) # match, mismatch, open, extend
    align_score = local_align_result[0][2]
    return local_align_result, align_score


def is_classified_correctly(label, ground_truth_df, indx):
    if label == ground_truth_df.loc[indx].region_label:
        return True
    else:
        return False


def main():
    print("main start")

    # Inputs:
    fastq_path = root_dir + data_set_dir + 'r2_synthetic_data_03-12-2021_17-36.fq'
    print("Working on fastq path {}".format(fastq_path))
    ground_truth_df_path = root_dir + data_set_dir + 'df_synthetic_data_03-12-2021_17-36.pkl'
    regions_kmers_path = root_dir + 'regionsKmersRaw_k8.txt'
    regions_sequences_path = root_dir + 'regionsRawSequence_k8.txt'
    bwa_miss_df_path = root_dir + data_set_dir + 'bwa_wrong_region.pkl'

    # Initialize objects
    ka_obj = kmer_aligner.KmerAlignerClass(8, 0.05, regions_kmers_path, regions_sequences_path)
    ko_obj = kmer_operations.KmerOperationsClass()
    se_obj = services.Services()

    # Read ground truth dataframe
    ground_truth_df = pd.read_pickle(ground_truth_df_path)
    regions_kmers_list = ko_obj.generate_regions_kmers_list(regions_kmers_path)
    regions_raw = se_obj.generate_regions_sequences_list(regions_sequences_path)

    # Read BWA miss step1 pikle df
    bwa_region_miss_df = pd.read_pickle(bwa_miss_df_path)

    correct_step_1_classification = 0
    correct_step_1_vs_BWA = 0
    for idx, read in enumerate(ka_obj.get_fastq_lines2(fastq_path)):
        # read classification to region
        read_name = read[0]
        read_name = read_name[1:].split('/')[0]
        read = read[1].replace("\n", "")
        max_similarity, label, running_time = ka_obj.label_read_two_sided(read)
        if is_classified_correctly(label, ground_truth_df, idx):
            correct_step_1_classification +=1
            if read_name in bwa_region_miss_df.read_name:
                correct_step_1_vs_BWA +=1
                print("V")



        # kmer local alignment
        # region_kmers = regions_kmers_list[label]
        # read_kmers = ko_obj.getKmersList(read, 8)
        # ko_obj.kmer_local_alignment2(read_kmers, region_kmers, 1, 2, idx)
        # SW Local alignment for reference
        # region_seq = regions_raw[label]
        # local_align_result, align_score = sw_local_alignment(read, region_seq)
        # align1, align2, score, local_align_start, local_align_end = local_align_result[0]
        # print(format_alignment(*local_align_result[0]))

    # Outuput
    print(" ---------------------------------------------------------------------------------------------- ")
    print("Total reads                                                          = {}".format(idx))
    print(
        "Region label match                                                     = {}".format(correct_step_1_classification / idx))
    print("Total Step 1 missed by BWA                                           = {}".format(len(bwa_region_miss_df)))
    print("Total Step 1 matches of BWA misses by Kmer                           = {}".format(correct_step_1_vs_BWA))
    # print("Total valid reads (remove missed regions)                           = {}".format(regions_match))
    # print("Deletions match to ground truth                                     = {}".format(
    #     bwa_deletions_match / regions_match))
    return


if __name__ == '__main__':
    main()
    sys.stdout.close()
