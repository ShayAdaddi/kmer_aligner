import random
import numpy as np
import pysam
import pandas as pd
from datetime import datetime
from Modules import analysis_module
import offline


def fetch_randomly_from_bedfile(chr_num,bed_df):
    # Random choice from entry of specific chr
    # If chr number does not exists in the bed file return -1
    # header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
    chr_regions_random_sample = bed_df[bed_df.chrom == "chr{}".format(chr_num)].sample()
    start = chr_regions_random_sample.chromStart
    end = chr_regions_random_sample.chromEnd
    return start,end

def generate_deletion_data():
    # Inputs:
    # Probability for deletion (0-1)
    # Deletion start location (rand inside read)
    # Deletion length (ran in range(a,b))

    # Define hyper parameters
    num_of_reads = 10
    deletion_prob = 1
    start_low, start_high = 20, 50
    del_len_low, del_len_high = 50, 80
    offset_low , offset_high = 1, 1000
    do_deletion_arr = [True,False]
    chr_low , chr_high = 1,22

    path_to_bamfile = '../../Data/Rag1.sorted.dedup.bam'
    bamfile = pysam.AlignmentFile(path_to_bamfile, "rb", ignore_truncation=True)
    bedfile_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/RAG1.bed'
    offline_obj = offline.OfflineClass()
    bed_df = offline_obj.generate_bed_file_df(bedfile_path)

    reads_meta_df = pd.DataFrame(columns=['index', 'chr_read', 'read_start', 'read_end', 'read', 'read_name',
                                          'deletion', 'deletion_start', 'deletion_length', 'new_sequence',
                                          'region_label','read_start_pos_in_region','read_end_pos_in_region'])

    an_obj = analysis_module.AnalysisClass()

    i = 0
    while i < num_of_reads:
        if i % 100 == 0: print("Working on read {}".format(i))
        deletion_start = random.randrange(start_low, start_high)
        deletion_length = random.randrange(del_len_low, del_len_high)
        do_deletion = np.random.choice(do_deletion_arr, 1, p=[deletion_prob,1-deletion_prob])
        chr_num = random.randrange(chr_low, chr_high)

        [start,end] = fetch_randomly_from_bedfile(chr_num,bed_df)
        offset = random.randrange(offset_low,offset_high)
        for read in bamfile.fetch('chr{}'.format(chr_num),start+offset,end+offset):
            if len(read.seq) == 140:
                break
        if do_deletion:
            new_seq = read.seq[:deletion_start] + read.seq[deletion_start+deletion_length:]
        else:
            new_seq = read.seq
        # get ground truth region
        read_region_label = an_obj.check_label(read.reference_start, read.reference_end, bed_df, full=True)

        if read_region_label == -1:
            continue

        region_start_pos = bed_df.loc[read_region_label].chromStart
        reads_meta_df.loc[i] = pd.Series(
            {'index': i, 'chr_read': read.reference_name, 'read_start': read.reference_start,
             'read_end': read.reference_end, 'read': read.seq, 'read_name': read.query_name,
             'deletion': do_deletion, 'deletion_start': deletion_start, 'deletion_length': deletion_length,
             'new_sequence': new_seq, 'region_label': read_region_label,
             'read_start_pos_in_region': read.reference_start - region_start_pos,
             'read_end_pos_in_region': read.reference_end - region_start_pos})

        i += 1

    # Write to file
    output_data_folder = "/export/home/fpd/PycharmProjects/kmer_aligner/Data/"
    today = datetime.today()
    timestamp = today.strftime("%d-%m-%Y_%H-%M")

    with open(output_data_folder + 'rag1_deletion_' + timestamp + '.fq', 'w+') as f:
        for index, row in reads_meta_df.iterrows():

            f.writelines("@{}/1\n".format(row.read_name))
            f.writelines("{}\n".format(row.new_sequence))
            f.writelines('+\n')
            f.writelines('B' * 140)
            f.writelines('\n')

    reads_meta_df.to_pickle(output_data_folder + 'rag1_deletion_df_' + timestamp + '.pkl')
    print('Done')


def main():
    print("Start")
    generate_deletion_data()
    print("Finished")


if __name__ == '__main__':
    main()
    print("Main Done")
