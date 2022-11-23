import argparse
from pathlib import Path
import ast
import khmer
import time
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import offline

sns.set()

NUM_OF_READS_IN_FASTQ = 3690324
SUBSAMPLE_READS = 100


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '\n': ''}
    new = ''.join([complement[base] for base in seq[::-1]])
    return new


def plot_kmer_read_cap_heatmap(similarity_mtrx):
    plt.figure(figsize=(18, 8))
    print('Plotting heatmap')
    ax = sns.heatmap(similarity_mtrx, cmap='Wistia')
    print(similarity_mtrx.shape)
    plt.xlabel('Iceberg index', fontsize=18)
    plt.ylabel('Read index', fontsize=18)
    plt.title('Read to Capture Similarity Confidence', fontsize=20)
    plt.show()


class KmerAlignerClass:

    def __init__(self, k, jaccard_threshold, regions_kmers_path, regions_sequences_path):
        self.k = k
        self.jaccard_threshold = jaccard_threshold
        # khmer
        self.alphabet = ['A', 'C', 'G', 'T', 'N']
        self.cg = khmer.Countgraph(self.k, len(self.alphabet) ** self.k, 1)
        # load regions
        regions = []
        with open(regions_kmers_path, 'r') as f:
            for iceberg_idx, iceberg in enumerate(f.readlines()):
                iceberg_indices = ast.literal_eval(iceberg)
                regions.append(set(iceberg_indices))
        self.regions = regions
        regions_seq = []
        with open(regions_sequences_path, 'r') as f:
            for line in f.readlines():
                regions_seq.append(line)
        self.regions_seq = regions_seq

    def get_fastq_lines(self, file):

        # if re.search('.gz$', file):
        #     fastq = gzip.open(file, 'rb')
        # else:
        fastq = open(file, 'r')
        with fastq as f:
            while True:
                # tic0 = time.perf_counter()
                l1 = f.readline()
                if not l1:
                    break
                l2 = f.readline()  # read string
                l3 = f.readline()
                l4 = f.readline()
                # tic1 = time.perf_counter()
                # print("get_fastq_lines Timing {} (sec)".format(tic1-tic0))
                yield l2

    def label_read_two_sided(self, read):
        # calc read kmers
        max_similarity = 0
        label = -1
        tic0 = time.perf_counter()
        read_kmers_set = set(self.cg.get_kmers(read))
        read_kmers_set_rc = set(self.cg.get_kmers(reverse_complement(read)))
        tic1 = time.perf_counter()
        for iceberg_indx, iceberg_kmers in enumerate(self.regions):
            sim_org = len(read_kmers_set.intersection(iceberg_kmers))
            sim_rc = len(read_kmers_set_rc.intersection(iceberg_kmers))
            if sim_org > max_similarity or sim_rc > max_similarity:
                max_similarity = max(sim_org, sim_rc)
                label = iceberg_indx

        return max_similarity, label, (tic1 - tic0)

    def label_read_debug(self, read):
        # calc read kmers
        max_similarity = 0
        label = -1
        read_kmers_set = set(self.cg.get_kmers(read))
        res = np.zeros(len(self.regions), dtype=float)
        for iceberg_indx, iceberg_kmers in enumerate(self.regions):
            res[iceberg_indx] = len(read_kmers_set.intersection(iceberg_kmers)) / len(
                read_kmers_set.union(iceberg_kmers))
        return res

    def label_read_two_sided_debug(self, read):
        # calc read kmers
        read_kmers_set = set(self.cg.get_kmers(read))
        read_kmers_set_rc = set(self.cg.get_kmers(reverse_complement(read)))
        res = np.zeros(len(self.regions), dtype=float)
        for iceberg_indx, iceberg_kmers in enumerate(self.regions):
            res[iceberg_indx] = max(len(read_kmers_set.intersection(iceberg_kmers)),
                                    len(read_kmers_set_rc.intersection(iceberg_kmers)))
        return res

    # def local_alignment(self, read, region):
    #     tic0 = time.perf_counter()
    #     res = pairwise2.align.globalms(read, region, 10, -10, -10, -9, penalize_end_gaps=(False, True),
    #                                    one_alignment_only=True)
    #     # local_align_start = region_start_location + get_offset(res[0][0])
    #     if (res[0][2]):
    #         local_align_score = res[0][2]
    #     tic1 = time.perf_counter()
    #     print("local alignment timing is {} ms ".format((tic1 - tic0) * 1000))
    #     return local_align_score

    def kmer_calc_pair_end_reads(self, r1, files_dir):

        avg_time = 0
        avg_time_kmers_only = 0
        abs_read_idx = 0
        # jac_sim_res = np.zeros((int(NUM_OF_READS_IN_FASTQ/SUBSAMPLE_READS)+1, len(self.regions)), dtype=float)
        jac_sim_res = []
        # local_align_scores = []
        for read_idx, read_r1 in enumerate(self.get_fastq_lines(str(files_dir / r1))):
            if read_idx % SUBSAMPLE_READS != 0:
                continue
            tic0 = time.perf_counter()
            [max_similarity1, iceberg_indx1, kmers_timing] = self.label_read_two_sided(read_r1)
            tic1 = time.perf_counter()
            jac_sim_res.append([max_similarity1, iceberg_indx1])
            # local_alignment_score = self.local_alignment(read_r1, self.regions_seq[iceberg_indx1])
            # local_align_scores.append(local_alignment_score)

            # jac_sim_res[abs_read_idx, iceberg_indx1] = max_similarity1
            # res = self.label_read_two_sided_debug(read_r1)
            # jac_sim_res[abs_read_idx,:] = res
            # print('LABEL1 - region label is {} with score of {}'.format(iceberg_indx1,max_similarity1))
            abs_read_idx += 1
            avg_time += (tic1 - tic0)
            avg_time_kmers_only += kmers_timing

        return jac_sim_res, avg_time / (NUM_OF_READS_IN_FASTQ / SUBSAMPLE_READS), avg_time_kmers_only / (
                    NUM_OF_READS_IN_FASTQ / SUBSAMPLE_READS)  # , local_align_scores


def main(args):
    try:

        offline_obj = offline.OfflineClass()
        regions_kmers_path = args['files_dir'] + args['path_to_regions_kmer']
        regions_sequences_path = args['files_dir'] + args['path_to_regions_seq']
        print("Start Init")
        tic0 = time.perf_counter()
        kc = KmerAlignerClass(args['k'], args['jaccard_threshold'], regions_kmers_path, regions_sequences_path)
        tic1 = time.perf_counter()
        print("End init - Timing {}(sec)".format((tic1 - tic0)))
        # res_mat, local_align_scores = kc.kmer_calc_pair_end_reads(args['read1_1'], args['curr_in_dir'])

        print("Start compute")
        tic2 = time.perf_counter()
        res_mat, duration, duration_kmers_only = kc.kmer_calc_pair_end_reads(args['read1_1'], args['curr_in_dir'])
        tic3 = time.perf_counter()
        print(
            "End compute for {} reads - Timing {}(sec), avg per read of label_read_two_sided {}, kmers only {} (sec)".format(
                (NUM_OF_READS_IN_FASTQ / SUBSAMPLE_READS), (tic3 - tic2), duration, duration_kmers_only))
        # print(res_mat)
        # Show results
        # plot_kmer_read_cap_heatmap(res_mat)

    except Exception as e:
        raise e


def parse_args():
    """
    parse_args:
    Parses steps arguments (used only in the --help command).
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--files_dir', help="ff")
    parser.add_argument('--read1_1', help="ff")
    parser.add_argument('--read2_1', help="ff")
    # parser.add_argument('index1_1', help="ff")
    # parser.add_argument('index2_1', help="ff")
    parser.add_argument('--out_dir', help="ff")
    parser.add_argument('--tag1', default='', help="ff")
    parser.add_argument('--tag2', default='', help="ff")
    parser.add_argument('--read1_2', default='', help="ff", metavar='')
    parser.add_argument('--read2_2', default='', help="ff", metavar='')
    # parser.add_argument('--index1_2', default='', help="ff", metavar='')
    # parser.add_argument('--index2_2', default='', help="ff", metavar='')
    parser.add_argument('--k', default=6, help="ff", type=int, metavar='')
    parser.add_argument('--jaccard_threshold', default=0.05, type=float, help="ff", metavar='')
    parser.add_argument('--path_to_regions_kmer', default="", help="ff", type=str, metavar='')
    parser.add_argument('--path_to_regions_seq', default="", help="ff", type=str, metavar='')

    args, rest = parser.parse_known_args()
    return dict([(k, v) for k, v in vars(args).items()])


if __name__ == '__main__':
    args = parse_args()
    args['curr_in_dir'] = Path(args['files_dir'])
    print(args)
    main(args)
