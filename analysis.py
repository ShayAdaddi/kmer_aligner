import argparse
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import collections
from numpy.distutils.command.config import config
import seaborn as sns

sns.set()


def parse_args():
    """
    parse_args:
    Parses steps arguments (used only in the --help command).
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--path_to_reads_df', default="", help="Full path to the reads metadata df", type=str,
                        metavar='')
    parser.add_argument('--bed_filename', default='Data/RAG1.bed', help="Filename of assay bed file", type=str,
                        metavar='')
    args, rest = parser.parse_known_args()
    return dict([(k, v) for k, v in vars(args).items()])


def check_overlap(read_start, read_end, regions_df):
    # Assuming no overlaps between regions
    min_overlap = 1
    region_index_low = regions_df.loc[
        ((read_end - regions_df['chromStart']) > min_overlap) & (read_start <= regions_df['chromStart'])
        ].index.tolist()

    region_index_high = regions_df.loc[
        ((regions_df['chromEnd'] - read_start) > min_overlap) & (read_end >= regions_df['chromEnd'])
        ].index.tolist()

    if (not region_index_low) & (not region_index_high):
        return [0, -1]

    else:
        if region_index_low:
            region_start = regions_df.loc[region_index_low[0]].at['chromStart']
            ov = read_end - region_start
            return ov.item(), region_index_low
        if region_index_high:
            region_end = regions_df.loc[region_index_high[0]].at['chromEnd']
            ov = region_end - read_start
            return ov.item(), region_index_high


def confusion_matrix(df):
    y_real = df['label_region']
    y_pred = df['predict_region']

    # Calculate Metrics
    total_classified = len(y_pred)
    match_mask = (y_real == y_pred)
    match = match_mask.sum()
    non_overlap_mask = (y_real == -1)
    non_overlap = (y_real == -1).sum()
    false_match = total_classified - match - non_overlap
    valid_classification = total_classified - non_overlap

    # Output
    print(" ---------------------------------------------------------------------------------------------- ")
    print("Total reads                      = {}".format(total_classified))
    print("Not overlapped reads             = {}".format(non_overlap))
    print("Total valid reads                = {}".format(valid_classification))
    print("Accurate classification          = {}".format(match / valid_classification))
    print("False classification             = {}".format(false_match / valid_classification))
    print("Matches AVG Jaccard Similarity   = {}".format(df[match_mask]['predict_confidence'].mean()))
    print("Matches Jaccard Similarity STD   = {}".format(df[match_mask]['predict_confidence'].std()))
    print("False AVG Jaccard Similarity     = {}".format(df[non_overlap_mask]['predict_confidence'].mean()))
    print("False Jaccard Similarity STD     = {}".format(df[non_overlap_mask]['predict_confidence'].std()))
    print(" ---------------------------------------------------------------------------------------------- ")

    return


def non_overlap_histogram(bed_df, df):
    # extract non overlap reads
    y_real = df['label_region']
    non_overlap_df = df[(y_real == -1)]
    abs_index = 0
    counts = {}
    for index, read in non_overlap_df.iterrows():
        if abs_index % 10000 == 0:
            print("Running non overlap read #{}".format(index))
        bed_df_chrom = bed_df[bed_df['chrom'] == read.chr]
        overlap, region = check_overlap(read.start, read.end, bed_df_chrom)
        if overlap not in counts:
            counts[overlap] = 0
        counts[overlap] += 1
        abs_index += 1

    # today = datetime.today()
    # timestamp = today.strftime("%d-%m-%Y_%H-%M")
    # filename = "Data/non_overlap_dict_" + timestamp + '.pickle'
    # file_to_write = open(filename, "wb")
    # pickle.dump(counts, file_to_write)
    plot_overlapped_histogram(counts)

    return


def accumulated_dict(d):
    overlaps = []
    values = []
    od = collections.OrderedDict(sorted(d.items()))
    for k in od.keys():
        overlaps.append(k)
        values.append(od[k])

    summed = np.cumsum(values)
    accum_dict = dict(zip(overlaps, summed))
    return accum_dict


def plot_overlapped_histogram(dict, islog=True):
    plt.figure(figsize=(14, 7))
    plt.style.use('seaborn-whitegrid')
    plt.bar(dict.keys(), dict.values(), log=islog, facecolor='#2ab0ff', edgecolor='#169acf', linewidth=0.5, alpha=0.7)
    plt.title("Non Overlapped Reads Histogram")
    plt.xticks(rotation=90)
    plt.xlabel("Overlap")
    plt.ylabel("Read Count")
    plt.show()


def plot_local_align_attrib(df, attrib):
    X = list(range(len(df)))
    plt.scatter(x=X,
                y=df[attrib],
                c='DarkBlue',
                label=attrib,
                alpha=0.7)

    y_avg = df[attrib].mean()
    y_std = df[attrib].std()
    y_median = df[attrib].median()

    y_avg_vec = [y_avg] * len(X)
    y_avg_std_top_vec = [(y_avg + y_std)] * len(X)
    y_median_vec = [y_median] * len(X)

    plt.plot(X, y_avg_vec, color='red', lw=1, ls='--', label="{} mean".format(attrib), alpha=0.7)
    plt.text(0, y_avg_vec[0], str(y_avg_vec[0]))

    plt.plot(X, y_median_vec, lw=1, ls='dashed', label="{} median".format(attrib), alpha=0.7)
    plt.text(0, y_median_vec[0], str(y_median_vec[0]))

    plt.plot(X, y_avg_std_top_vec, color='red', lw=1, ls='--', label="{} std".format(attrib), alpha=0.7)
    plt.text(0, y_avg_std_top_vec[0], str(y_avg_std_top_vec[0]))

    plt.title("Local align {} plot".format(attrib))
    plt.xlabel("Reads")
    plt.ylabel("{}".format(attrib))
    plt.legend()
    plt.show()
    return


def plot_local_align_attrib_histogram(df, attrib):
    sns.displot(df, x=attrib, log=True)
    plt.title("Local align {} log distribution plot".format(attrib))
    plt.show()
    return


def plot_local_align_multi_attrib(df, attrib_list):

    X = list(range(len(df)))
    for attr in attrib_list:
        plt.scatter(x=X,
                    y=df[attr],
                    label=attr,
                    alpha=0.7)
        y_avg = df[attr].mean()
        y_std = df[attr].std()
        y_median = df[attr].median()

        y_avg_vec = [y_avg] * len(X)
        y_avg_std_top_vec = [(y_avg + y_std)] * len(X)
        y_median_vec = [y_median] * len(X)

        plt.plot(X, y_avg_vec, lw=1, ls='--', label="{} mean".format(attr), alpha=0.7)
        plt.text(0, y_avg_vec[0], str(y_avg_vec[0]))

        plt.plot(X, y_median_vec, lw=1, ls='dashed', label="{} median".format(attr), alpha=0.7)
        plt.text(0, y_median_vec[0], str(y_median_vec[0]))

        plt.plot(X, y_avg_std_top_vec, lw=1, ls='dotted', label="{} std".format(attr), alpha=0.7)
        plt.text(0, y_avg_std_top_vec[0], str(y_avg_std_top_vec[0]))

    plt.title("Local align attributes plot")
    plt.xlabel("Reads")
    plt.ylabel("Attributes counts")
    plt.legend()
    plt.show()
    return


def main(args):
    df_path = args['path_to_reads_df']
    df = pd.read_pickle(df_path)

    # Generate Bedfile DF
    df_bed = pd.read_csv(args['bed_filename'], sep='\t', comment='t', header=None)
    header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
    df_bed.columns = header[:len(df_bed.columns)]

    confusion_matrix(df)

    # plot_local_align_attrib(df, 'start_gap')
    # plot_local_align_attrib(df, 'end_gap')
    plot_local_align_attrib(df, 'sequence_length_gap')

    # plot_local_align_multi_attrib(df, ['sequence_length_gap'])
    # plot_local_align_multi_attrib(df, ['start_gap', 'end_gap', 'sequence_length_gap'])

    plot_local_align_attrib_histogram(df,'sequence_length_gap')

    # print(df[df.predict_confidence_second_max == df.predict_confidence_second_max.max()])
    # non_overlap_histogram(df_bed, df)
    # file_to_read = open("Data/non_overlap_dict_09-08-2021_11-03.pickle", "rb")
    # non_overlapped_dict = pickle.load(file_to_read)
    # plot_overlapped_histogram(non_overlapped_dict)
    # accum = accumulated_dict(non_overlapped_dict)
    # plot_overlapped_histogram(accum,islog=False)

    print("done")


if __name__ == '__main__':
    args = parse_args()
    main(args)
