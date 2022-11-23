import random
import pysam
import pandas as pd
from datetime import datetime
from Bio.Seq import Seq
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
sns.set
from Modules import services, analysis_module
import kmer_aligner, offline


def plot_local_align_attrib_histogram(df, attrib):
    sns.displot(df, x=attrib, log=True)
    plt.title("{} log distribution plot".format(attrib))
    plt.show()
    return


def generate_translocation_data():
    # Read will be generated from two reads that located in different chromosomes according to raffle of split point
    # READ1 ACGTACGTTACGAGATCTAGATCGG
    #       -------------------------*--------------------------------
    # READ2                           ATTCTTTTACGAGAGGGACTAGACGGGACTAT

    # * = Split point (raffled for each new read generated)

    # Define hyper parameters
    generated_reads_number = 1
    split_point_low = 40
    split_point_high = 100
    chr_low = 1
    chr_high = 22
    offset_low, offset_high = 1, 1000

    path_to_bamfile = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/Rag1.sorted.dedup.bam'
    bamfile = pysam.AlignmentFile(path_to_bamfile, "rb", ignore_truncation=True)

    reads_meta_df = pd.DataFrame(columns=['index', 'chr_read1', 'read1_start', 'read1_end', 'read1', 'read1_name',
                                          'region_label1','chr_read2', 'read2_start', 'read2_end','read2',
                                          'read2_name','region_label2','split_point', 'sequence', 'mate_sequence'])

    an_obj = analysis_module.AnalysisClass()
    bedfile_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/RAG1.bed'
    offline_obj = offline.OfflineClass()
    bed_df = offline_obj.generate_bed_file_df(bedfile_path)
    services_obj = services.Services()

    i = 0
    while i < generated_reads_number:
        if i % 100 == 0: print("Working on read {}".format(i))
        read1_chr_num = random.randrange(chr_low, chr_high)
        read2_chr_num = random.randrange(chr_low, chr_high)
        split_point = random.randrange(split_point_low, split_point_high)

        [start, end] = services_obj.fetch_randomly_from_bedfile(read1_chr_num, bed_df)
        offset = random.randrange(offset_low, offset_high)
        for read1 in bamfile.fetch('chr{}'.format(read1_chr_num), start + offset, end + offset):
            if len(read1.seq) == 140:
                break
        read1_region_label = an_obj.check_label(read1.reference_start, read1.reference_end, bed_df, full=True)
        if read1_region_label == -1:
            continue

        [start, end] = services_obj.fetch_randomly_from_bedfile(read2_chr_num, bed_df)
        offset = random.randrange(offset_low, offset_high)
        for read2 in bamfile.fetch('chr{}'.format(read2_chr_num), start + offset, end + offset):
            if len(read2.seq) == 140:
                break
        read2_region_label = an_obj.check_label(read2.reference_start, read2.reference_end, bed_df, full=True)
        if read2_region_label == -1:
            continue

        translocation_seq = read1.seq[:split_point] + read2.seq[split_point:]
        seq2 = Seq(translocation_seq)
        mate = seq2.reverse_complement()
        reads_meta_df.loc[i] = pd.Series(
            {'index': i, 'chr_read1': read1.reference_name, 'read1_start': read1.reference_start,
             'read1_end': read1.reference_end, 'read1': read1.seq, 'read1_name': read1.query_name,
             'region_label1': read1_region_label,'chr_read2': read2.reference_name,
             'read2_start': read2.reference_start, 'read2_end': read2.reference_end,'read2': read2.seq,
             'read2_name':read2.query_name, 'region_label2': read2_region_label,'split_point': split_point,
             'sequence': translocation_seq, 'mate_sequence': mate})

        i += 1

    # Write to file
    today = datetime.today()
    timestamp = today.strftime("%d-%m-%Y_%H-%M")

    with open('Data/' + 'rag1_translocation_R1_' + timestamp + '.fq', 'w+') as f1, \
            open('Data/' + 'rag1_translocation_R2_' + timestamp + '.fq', 'w+') as f2:
        for index, row in reads_meta_df.iterrows():
            f1.writelines("@{}/1\n".format(row.read1_name))
            f1.writelines("{}\n".format(row.sequence))
            f1.writelines('+\n')
            f1.writelines('B' * 140)
            f1.writelines('\n')

            f2.writelines("@{}/2\n".format(row.read1_name))
            f2.writelines("{}\n".format(row.mate_sequence))
            f2.writelines('+\n')
            f2.writelines('B' * 140)
            f2.writelines('\n')

    reads_meta_df.to_pickle('Data/' + 'rag1_translocation_df_' + timestamp + '.pkl')
    print('Done')


def classify_reads():
    # ---------------------- Classify -------------- #
    services_obj = services.Services()
    regions_kmers_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/regionsKmersRaw_k8.txt'
    regions_sequences_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/regionsRawSequence_k8.txt'
    fastq_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/rag1_translocation_R1_23-10-2021_22-28.fq'

    k, jaccard_th = 8, 0.02
    kmer_aligner_obj = kmer_aligner.KmerAlignerClass(k, jaccard_th, regions_kmers_path, regions_sequences_path)
    jac_sim_res = []
    for read_idx, read_r1 in enumerate(services_obj.get_fastq_lines(fastq_path)):
        [max_similarity, region_indx, kmers_timing] = kmer_aligner_obj.label_read_two_sided(read_r1)
        jac_sim_res.append([max_similarity, region_indx])

    print("Done")
    return jac_sim_res


def main():
    print("Start")
    generate_translocation_data()
    # result = classify_reads()
    # today = datetime.today()
    # timestamp = today.strftime("%d-%m-%Y_%H-%M")
    # with open('Data/classify_reads_result_' + timestamp + '.pkl', 'wb') as f:
    #     pickle.dump(result, f)


if __name__ == '__main__':
    main()
    print("main Done")
