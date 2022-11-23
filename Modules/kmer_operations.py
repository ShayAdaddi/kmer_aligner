import khmer
import numpy as np
import collections
import matplotlib.pyplot as plt
import itertools
import ast
from datetime import datetime

def reject_outliers(data, m, max_dist):
    d = np.array(data,dtype=int)
    dist = min(m * np.std(d), max_dist)
    d_after_process = d[abs(d - np.median(d)) < dist]
    if d_after_process.size == 0:
        d_after_process = np.median(d)
    return d_after_process


def apply_sliding_window(all_pos,read_length):
    all_pos_set = set(all_pos)
    window_size = read_length*0.7
    median = np.median(all_pos)
    max_intersection = 0
    max_i = -1
    for i in range(int(median-window_size), int(median+window_size)):
        i_start = max(0,i-read_length/2)
        i_end = i+read_length/2
        current_window = range(int(i_start), int(i_end))
        intersected = len(set(current_window).intersection(all_pos_set))
        if intersected > max_intersection:
            max_intersection = intersected
            max_i = i
            indx_first = i_start
            indx_last = i_end

    return indx_first, indx_last


class KmerOperationsClass:
    def __init__(self):
        return

    def kmer_local_alignment_debug(self, sequence, reference, k):
        alphabet = ['A', 'C', 'G', 'T', 'N']
        cg = khmer.Countgraph(k, len(alphabet) ** k, 1)

        seq_kmers = cg.get_kmers(sequence)
        ref_kmers = cg.get_kmers(reference)
        seq_kmers_set = set(cg.get_kmers(sequence))
        ref_kmers_set = set(cg.get_kmers(reference))

        intersection = len(seq_kmers_set.intersection(ref_kmers_set))
        union = len(seq_kmers_set.union(ref_kmers_set))
        if union == 0:
            print('Divide by zero')
            return

        for element_first in seq_kmers:
            try:
                indx_first = ref_kmers.index(element_first)
                break
            except ValueError as e:
                print(f'Element "{element_first}" not found in the list: ', e)
                continue

        for element_last in seq_kmers[::-1]:
            try:
                indx_last = ref_kmers.index(element_last)
                break
            except ValueError as e:
                print(f'Element "{element_last}" not found in the list: ', e)
                continue

        ############################################## DEBUG print ###################################################
        print("Sequence:                                            {}".format(sequence))
        print("Reference                                            {}".format(reference))
        print("First kmer matched                                   [{}] position {}".format(element_first, indx_first))
        print("First kmer matched                                   {}".format(' ' * indx_first + u'\u2191'))
        print("Last kmer matched                                    [{}] position {}".format(element_last, indx_last))
        print("Last kmer matched                                    {}".format(' ' * indx_last + u'\u2191'))
        print("Union={} Intersection={} Jaccard={}".format(union, intersection, intersection / union))
        print("Length gap between sequence length to anchor points={}".format(
            abs(indx_last - indx_first - len(sequence))))

        return

    def kmer_local_alignment(self, sequence_kmer_list, reference_kmer_list):

        element_first, indx_first, element_last, indx_last = None, None, None, None

        for element_first in sequence_kmer_list:
            try:
                indx_first = reference_kmer_list.index(element_first)
                break
            except ValueError as e:
                # print(f'Element "{element_first}" not found in the list')
                continue

        for element_last in sequence_kmer_list[::-1]:
            try:
                indx_last = reference_kmer_list.index(element_last)
                break
            except ValueError as e:
                # print(f'Element "{element_last}" not found in the list')
                continue

        return [element_first, indx_first, element_last, indx_last]

    def kmer_local_alignment2(self, sequence_kmer_list, reference_kmer_list,debug,method,read_index):
        intersected_kmers = set(sequence_kmer_list).intersection(set(reference_kmer_list))
        d = {}  # Final dictionary
        all_pos = []
        all_pos_updated = []
        for element in intersected_kmers:
            indexes = []
            for i, j in enumerate(reference_kmer_list):
                if j == element:
                    indexes.append(i)
            d[element] = indexes
            if indexes is not None:
                all_pos.extend(indexes)
                # print(f"Item {element} is found at index {indexes}")

        avg_pos = sum(all_pos) / len(all_pos)
        for element in d.keys():
            tmp = d[element]
            val = min(d[element], key=lambda x: abs(x - avg_pos))
            d[element] = val
            if tmp[0] != val:
                print("Updated kmer pos {} to {}".format(tmp[0],val))

        d_l = list(d.keys())
        element_first = d_l[0]
        indx_first = d[d_l[0]]
        element_last = d_l[-1]
        indx_last = d[d_l[-1]]

        read_length = 140
        if method == 1:
            kmers_median = np.median(all_pos)
            indx_first = kmers_median - read_length/2
            indx_last = kmers_median + read_length/2

        elif method == 2:
            [indx_first, indx_last ] = apply_sliding_window(all_pos,read_length)

        elif method == 3:
            all_pos_updated = reject_outliers(all_pos,3,100)
            all_pos_updated_sorted = np.sort(all_pos_updated)
            indx_first = all_pos_updated_sorted[0]
            indx_last = all_pos_updated_sorted[-1]
            element_first = ""
            element_last = ""
        ### DEBUG ####
        if debug:
            print("The length of the kmer local alignment is {}-{}={}".format(indx_last,indx_first,indx_last-indx_first))
            arr = np.zeros(len(reference_kmer_list))
            arr_updated = np.zeros(len(reference_kmer_list))
            for item in all_pos:
                arr[item] = 1
            for item in all_pos_updated:
                arr_updated[item] = 1
            plt.plot(arr,'o')
            plt.plot(arr_updated, 'x')
            # Plot the align segment
            point1 = [indx_first, 1]
            point2 = [indx_last, 1]
            x_values = [point1[0], point2[0]]
            y_values = [point1[1], point2[1]]

            plt.plot(x_values, y_values,'r','-')
            plt.title('Kmer matches positions')
            plt.xlabel('Position in region')
            plt.ylabel('1=Match 0=No match')
            plt.show()
            # today = datetime.today()
            # timestamp = today.strftime("%d-%m-%Y_%H-%M")
            # filename = '/export/home/fpd/PycharmProjects/kmer_aligner/Outputs/deletion_kmer_local_alignment_read_{}_method_{}_{}'.format(read_index,method,timestamp)
            # plt.savefig(filename)
            # plt.clf()

        return [element_first, indx_first, element_last, indx_last]

    def kmer_to_List(kmer_dict):
        # TODO: Add amount of kmer repetitions
        l = []
        for k, v in kmer_dict.items():
            l.append(k)

        return l

    def kmer_to_sorted_unique_idx_list(kmer_dict, idx_dict):
        # TODO: Add amount of kmer repetitions
        idx_list = []

        for k, v in kmer_dict.items():
            idx_list.append(idx_dict[k])
        arr = np.sort(np.unique(np.array(idx_list)))

        return list(arr)

    def kmer_counter_calculation(seq, k):
        c = collections.Counter()
        # Calculate how many kmers of length k there are
        num_kmers = len(seq) - k + 1
        for i in range(num_kmers):
            c.update({seq[i:i + k]: 1})

        return c

    def generate_kmer_permutations_dictionary(K):
        alphabet = 'ACGTN'
        permutations_dict = {}
        idx_dict = {}
        for idx, output in enumerate(itertools.product(alphabet, repeat=K)):
            permutation = ''.join(output)
            permutations_dict[idx] = permutation
            idx_dict[permutation] = idx

        return permutations_dict, idx_dict

    def getKmersList(self, sequence, k):
        alphabet = ['A', 'C', 'G', 'T', 'N']
        cg = khmer.Countgraph(k, len(alphabet) ** k, 1)
        seq_kmers = cg.get_kmers(sequence)
        return seq_kmers


    def generate_regions_kmers_list(self,path):
        regions_kmers_list = []
        with open(path, 'r') as f:
            for region_idx, region in enumerate(f.readlines()):
                region_kmers = ast.literal_eval(region)
                regions_kmers_list.append(region_kmers)
        return regions_kmers_list