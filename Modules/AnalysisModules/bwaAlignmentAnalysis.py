import pysam
import os
import pandas as pd
from Modules import analysis_module
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_theme(style="darkgrid")


def generate_bed_file_df(bed_file_path):
    # Generate Bedfile DF
    df = pd.read_csv(bed_file_path, sep='\t', comment='t', header=None)
    header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
    df.columns = header[:len(df.columns)]
    return df


os.chdir('../../')
root_dir = os.getcwd() + '/Data/data_set_03-12-2021_17-36/'
print(root_dir)

bedfile_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/RAG1_-Capture-seq-List-CRISPR-regions-hg38-final.bed'
bed_df = generate_bed_file_df(bedfile_path)

# # Load ground truth dataframe
excel_path = root_dir + 'df_synthetic_data.xlsx'

# Load bwa result
bam_path = root_dir + 'bwa1.bam'

bamfile = pysam.AlignmentFile(bam_path, "rb", ignore_truncation=True)

gt_df = pd.DataFrame(columns=['index', 'read_name', 'chr', 'read_start', 'read_end', 'region_label', 'del_len'])
bwa_df = pd.DataFrame(
    columns=['index', 'read_name', 'chr', 'read_start', 'read_end', 'region_label', 'del_len', 'read_len_no_clipping'])

gt_df_pkl_path = root_dir + 'gt_df.pkl'
bwa_df_pkl_path = root_dir + 'bwa_df.pkl'

analysis_obj = analysis_module.AnalysisClass()
#
# for indx, read in enumerate(bamfile.fetch(until_eof=True)):
#
#     gt_read_name = read.qname
#     gt_read_name_parts = gt_read_name.split(':')
#     del_start_pos = int(gt_read_name_parts[3].split('=')[1])
#     del_end_pos = int(gt_read_name_parts[4].split('=')[1])
#     del_len = 0
#
#     if read.is_read1:
#         read_end = int(gt_read_name_parts[5].split('=')[1])
#         read_start = read_end - 140
#         # check if read start is inside deletions segment
#         if read_start < del_end_pos:
#             read_start = read_start - del_len
#             del_len = del_end_pos - del_start_pos
#
#     if read.is_read2:
#         read_start = int(gt_read_name_parts[2].split('=')[1])
#         read_end = read_start + 140
#         # check if read end is inside deletions segment
#         if read_end > del_start_pos:
#             read_end = read_end + del_len
#             del_len = del_end_pos - del_start_pos
#
#     gt_read_info = pd.Series({'index': indx, 'read_name': gt_read_name, 'chr': gt_read_name_parts[1],
#                               'read_start': read_start, 'read_end': read_end,
#                               'region_label': int(gt_read_name_parts[0].split('_')[1]),
#                               'del_len': del_len})
#     gt_df.loc[indx] = gt_read_info
#
#     del_len_list_from_cigar = [int(item[1]) for item in read.cigartuples if 2 == item[0]]
#     if del_len_list_from_cigar:
#         del_len = del_len_list_from_cigar[0]
#     else:  # list is empty
#         del_len = -1
#
#     read_len_no_clipping = sum([int(item[1]) for item in read.cigartuples if (item[0] != 4 and item[0] != 5)])
#     label = analysis_obj.check_label(read.reference_start, read.reference_end, bed_df)
#     bwa_read_info = pd.Series({'index': indx, 'read_name': read.qname, 'chr': read.reference_name,
#                                'read_start': read.reference_start, 'read_end': read.reference_end,
#                                'region_label': label, 'del_len': del_len, 'read_len_no_clipping': read_len_no_clipping})
#     bwa_df.loc[indx] = bwa_read_info
#
# gt_df.to_pickle(gt_df_pkl_path)
# bwa_df.to_pickle(bwa_df_pkl_path)
# #
gt_df = pd.read_pickle(gt_df_pkl_path)
bwa_df = pd.read_pickle(bwa_df_pkl_path)

print("Size of ground truth dataframe: {}".format(len(gt_df)))
print("Size of bwa dataframe: {}".format(len(bwa_df)))

total_reads = len(gt_df)
# comparing region labels
diff_df = gt_df['region_label'].compare(bwa_df['region_label'])
bwa_wrong_region = bwa_df.loc[diff_df.index]
bwa_wrong_region.to_pickle(root_dir + 'bwa_wrong_region.pkl')

regions_nonmatch = len(diff_df)
regions_match = total_reads - regions_nonmatch

headers = ['index', 'gt_start', 'gt_end', 'bwa_start', 'bwa_end', 'gt_end_start', 'bwa_end_start', 'diff_start',
           'diff_end', 'gt_del_len', 'bwa_del_len','gt_len_with_del', 'bwa_read_len_no_clipping']
data = [gt_df['index'], gt_df['read_start'], gt_df['read_end'], bwa_df['read_start'], bwa_df['read_end'],
        abs(gt_df['read_end'] - gt_df['read_start']), abs(bwa_df['read_end'] - bwa_df['read_start']),
        abs(gt_df['read_start'] - bwa_df['read_start']), abs(gt_df['read_end'] - bwa_df['read_end']),
        gt_df['del_len'], bwa_df['del_len'], gt_df['del_len'] + 140 ,bwa_df['read_len_no_clipping']]

combined_df = pd.concat(data, axis=1, keys=headers)
combined_df['gt_len'] = 140


# remove non match regions
combined_df = combined_df.drop(diff_df.index)

combined_df_bwa_del_len = combined_df[combined_df['bwa_del_len'] > 0]
combined_df['del_diff'] = combined_df['bwa_del_len'] - combined_df['gt_del_len']
combined_df_bwa_len_match = combined_df_bwa_del_len[
    combined_df_bwa_del_len['bwa_del_len'] == combined_df_bwa_del_len['gt_del_len']]
bwa_deletions_match = len(combined_df_bwa_len_match)

combined_df_bwa_read_len_match = combined_df[combined_df['bwa_read_len_no_clipping'] == 140]
bwa_read_len_match = len(combined_df_bwa_read_len_match)

# Output BWA Comparison
print(" ---------------------------------------------------------------------------------------------- ")
print("Total reads                                                         = {}".format(total_reads))
print("Region label match                                                  = {}".format(regions_match / total_reads))
print("Total valid reads (remove missed regions)                           = {}".format(regions_match))
print("Deletions match to ground truth                                     = {}".format(
    bwa_deletions_match / regions_match))

# Plots


# ax = combined_df_bwa_del_len.iloc[0:50].plot(x='index', y=['gt_del_len', 'bwa_del_len'], kind='bar')  # .iloc[0:10]
# ax.set_xlabel("Reads index")
# ax.set_ylabel("Deletion count (bps)")
# ax.set_title("BWA Deletions comparison to ground truth for {}".format(root_dir.split('/')[-2]))
# plt.show()
# # Partial plot due to large data on bar plot
# ax = combined_df.iloc[0:50].plot(x='index', y=['gt_len_with_del', 'bwa_read_len_no_clipping'], kind='bar')  # .iloc[0:10] # stacked=True
# ax.set_xlabel("Reads index")
# ax.set_ylabel("Read length (bps)")
# ax.set_title("BWA read length comparison to ground truth for {}".format(root_dir.split('/')[-2]))
# plt.show()

# Scatter plots
# sns.scatterplot(x="gt_del_len", y="bwa_del_len", data=combined_df);
# plt.show()
#
# sns.scatterplot(x="gt_len_with_del", y="bwa_read_len_no_clipping", data=combined_df);
# plt.show()

# Deletions histogram
# fig, axes = plt.subplots(1, 2)
# fig.suptitle('Deletions distribution for {}'.format(root_dir.split('/')[-2]))
# g0 = sns.histplot(ax=axes[0], data=combined_df, x="gt_del_len")
# axes[0].set_title('Ground truth deletions distribution')
# g0.set(ylim=(0, 500))
# g1 = sns.histplot(ax=axes[1], data=combined_df, x="bwa_del_len")
# axes[1].set_title('BWA result deletions distribution')
# g1.set(ylim=(0, 500))
# plt.show()

# Deletions diff histogram
sns.histplot(data=combined_df, x="del_diff", binwidth=1).set(title='Deletions diff distribution for {}'.format(root_dir.split('/')[-2]))
plt.show()


# Output Kmers
# print(" ---------------------------------------------------------------------------------------------- ")
# print("Total reads                      = {}".format(total_classified))
# print("Not overlapped reads             = {}".format(non_overlap))
# print("Total valid reads                = {}".format(valid_classification))
# print("Accurate classification          = {}".format(match / valid_classification))
# print("False classification             = {}".format(false_match / valid_classification))
# print("Matches AVG Jaccard Similarity   = {}".format(df[match_mask]['predict_confidence'].mean()))
# print("Matches Jaccard Similarity STD   = {}".format(df[match_mask]['predict_confidence'].std()))
# print("False AVG Jaccard Similarity     = {}".format(df[non_overlap_mask]['predict_confidence'].mean()))
# print("False Jaccard Similarity STD     = {}".format(df[non_overlap_mask]['predict_confidence'].std()))
# print(" ---------------------------------------------------------------------------------------------- ")

print("Done")
