import argparse
import numpy as np
import pandas as pd
from pyfaidx import Fasta
import pysam
import ast
import khmer
from datetime import datetime
import pickle
from Modules import visualizations
from Modules import analysis_module

def parse_args():
    """
    parse_args:
    Parses steps arguments (used only in the --help command).
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--path_to_base_dir', default="", help="Full path to the folder with the input files", type=str,
                        metavar='')
    parser.add_argument('--bed_filename', default='RAG1_-Capture-seq-List-CRISPR-regions-hg38-final.bed',
                        help="file name of the assay bed file", type=str, metavar='')
    parser.add_argument('--bam_filename', default='Rag1_chr11_0_10000000.bam', help="file name of the bam file",
                        type=str, metavar='')
    parser.add_argument('--regions_kmers_filename', default="",
                        help="file name of the regions kmers calculated offline", type=str, metavar='')
    parser.add_argument('--regions_sequences_filename', default="",
                        help="file name of the regions sequences calculated offline", type=str, metavar='')
    parser.add_argument('--genome_filename', default="", help="file name of reference genome", type=str, metavar='')
    # parser.add_argument('--output_dir', help="ff")
    parser.add_argument('--k', default=8, help="what [k] to use for kmer calculation", type=int, metavar='')
    parser.add_argument('--jaccard_threshold', default=0.05, type=float, help="jaccard similarity threshold",
                        metavar='')
    parser.add_argument('--subsample_reads', default=100, help="sample every X reads", type=int, metavar='')
    parser.add_argument('--reads_limit', default=1000000, help="define maximum amount of reads for memory constraints",
                        type=int, metavar='')
    parser.add_argument('--save', default=False, help="Save the results to pickle ?", type=bool, metavar='')

    args, rest = parser.parse_known_args()
    return dict([(k, v) for k, v in vars(args).items()])


def main(args):
    try:

        print("Start Kmer Aligner Unit test")

        # Generate Bedfile DF
        df = pd.read_csv(args['path_to_base_dir'] + args['bed_filename'], sep='\t', comment='t', header=None)
        header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
        df.columns = header[:len(df.columns)]
        # Load Genome
        genome = Fasta(args['path_to_base_dir'] + args['genome_filename'], sequence_always_upper=True)
        bamfile = pysam.AlignmentFile(args['path_to_base_dir'] + args['bam_filename'], "rb", ignore_truncation=True)
        # Extract regions kmers file into list
        regions_kmers = []
        with open(args['path_to_base_dir'] + args['regions_kmers_filename'], 'r') as f:
            for region_idx, region in enumerate(f.readlines()):
                region_kmers = ast.literal_eval(region)
                regions_kmers.append(set(region_kmers))

        alphabet = ['A', 'C', 'G', 'T', 'N']
        cg = khmer.Countgraph(args['k'], len(alphabet) ** args['k'], 1)
        # Generate results array
        num_of_reads_to_process = int(args['reads_limit'] / args['subsample_reads'])
        jac_sim_res = np.zeros((num_of_reads_to_process, len(regions_kmers)))
        print("Allocated results array for [{},{}] results".format(num_of_reads_to_process, len(regions_kmers)))
        a_obj = analysis_module.AnalysisClass()
        # Generate reads meta dataframe
        reads_meta_df = pd.DataFrame(columns=['index', 'chr', 'name', 'start', 'end',
                                              'label_region', 'predict_region', 'predict_confidence'
            , 'cigar', 'is_unmapped', 'second_max', 'predict_confidence_second_max', 'sequence'])
        read_idx = 0
        abs_read_idx = 0
        for read in bamfile.fetch(until_eof=True):
            if abs_read_idx == num_of_reads_to_process:
                break
            if read_idx % args['subsample_reads'] == 0:
                print("compute read #{} out of {} - completion percentage = {} %".format(abs_read_idx,
                                                                                         num_of_reads_to_process,
                                                                                         int((
                                                                                                     abs_read_idx + 1) * 100 / num_of_reads_to_process)))
                # if read.is_unmapped:
                # print("compute read #{}".format(abs_read_idx))
                # metadata preparations
                read_region_label = a_obj.check_label(read.reference_start, read.reference_end, df, full=True)
                # read_region_label = -1
                reads_meta_df.loc[abs_read_idx] = pd.Series(
                    {'index': abs_read_idx, 'chr': read.reference_name, 'name': read.qname,
                     'start': read.reference_start, 'end': read.reference_end,
                     'label_region': read_region_label, 'predict_region': -1,
                     'predict_confidence': -1, 'cigar': read.cigarstring, 'is_unmapped': read.is_unmapped,
                     'second_max': -1, 'predict_confidence_second_max': -1, 'sequence': read.seq})
                l = set(cg.get_kmers(read.seq))
                for region_idx, region in enumerate(regions_kmers):
                    intersect_kmers = l.intersection(region)
                    intersection = len(intersect_kmers)
                    union = len(l.union(region))
                    if union == 0:
                        print('Divide by zero')
                        continue

                    jac_sim_res[abs_read_idx, region_idx] = intersection / union

                abs_read_idx += 1
            read_idx += 1

        reads_label_region_mapping = np.argmax(jac_sim_res, axis=1)
        reads_meta_df['predict_region'] = reads_label_region_mapping
        reads_meta_df['predict_confidence'] = np.max(jac_sim_res, axis=1)
        jac_sim_res[:, reads_label_region_mapping] = 0
        reads_label_region_mapping_second_max = np.argmax(jac_sim_res, axis=1)
        reads_meta_df['second_max'] = reads_label_region_mapping_second_max
        reads_meta_df['predict_confidence_second_max'] = np.max(jac_sim_res, axis=1)

        # v_obj = visualizations.Visualizations()
        # v_obj.plot_kmer_read_cap_heatmap(jac_sim_res)
        # a_obj.calculate_confusion_matrix(reads_meta_df)
        # Export reads metadata dataframe to pickle
        if (args['save']):
            today = datetime.today()
            timestamp = today.strftime("%d-%m-%Y_%H-%M")
            reads_meta_df.to_pickle(args['path_to_base_dir'] + 'reads_meta_df_' + timestamp + '.pkl')

        # DEBUG
        # valid_df = reads_meta_df[reads_meta_df.label_region != -1]
        # valid_df = valid_df[valid_df.is_unmapped == False]
        # unmatch_df = valid_df[valid_df.predict_region != valid_df.label_region]
        # unmatch_df.to_pickle('Data/debug_df.pkl')


    except Exception as e:
        raise e


if __name__ == '__main__':
    args = parse_args()
    main(args)
