import collections
import os
import re
import gzip
import subprocess
import argparse
import logging
from pathlib import Path
import directories_operations
from loggers_factory import set_logging_config
import sys
import pandas as pd
#import khmer
tag_suffix = '.tag'
no_tag_suffix = '.no-tag'
trimmed_suffix = '.trimmed'
untrimmed_suffix = '.untrimmed'

tag_dir = 'TAG-RESULTS'


class tag_handler:

    def __init__(self, tag1, tag2, k,jaccard_threshold):
        self.tag1 = tag1
        self.tag2 = tag2
        self.k = k
        self.jaccard_threshold = jaccard_threshold

        #alphabet_length = 5  # ['A', 'C', 'G', 'T', 'N']
        #self.kmers_count_graph = khmer.Countgraph(self.k, alphabet_length ** self.k, 1)

        self.tag1_kmer_dict = self.naiveKmerCountsCalculation(tag1)
        self.tag2_kmer_dict = self.naiveKmerCountsCalculation(tag2)


    # iterate over fastq file.
    def get_fastq_lines(self, file):
        """
        get_fastq_lines:
        Yield one read from the given fastq file (for iterate the reads in the file).

        :param file: fastq file name.

        return: four lines which represent one read.
        """
        if re.search('.gz$', file):
            fastq = gzip.open(file, 'rb')
        else:
            fastq = open(file, 'r')
        with fastq as f:
            while True:
                l1 = f.readline()
                if not l1:
                    break
                l2 = f.readline()
                l3 = f.readline()
                l4 = f.readline()
                yield [l1, l2, l3, l4]

    def naiveKmerCountsCalculation(self, sequence):
        """
        naiveKmerCountsCalculation:
        Creates a kmers dictionary for the sequence
        (with kmers as keys and amount of appearance in the sequence as the values).

        return: The kmers dictionary.
        """
        k_mer_counts_dict = {}
        # Loop over the string (read/tag)
        for i in range(len(sequence) - self.k):
            kmer = sequence[i:i + self.k]
            # Add the kmer to the dictionary if it's not there or sum +1 to its value
            if kmer not in k_mer_counts_dict:
                k_mer_counts_dict[kmer] = 0
            k_mer_counts_dict[kmer] += 1
        # sort the kmer dictionary by amount
        #k_mer_counts_dict = collections.OrderedDict(sorted(k_mer_counts_dict.items()))
        return k_mer_counts_dict

    # def naiveKmerCountsCalculation(self, sequence):
    #     """
    #     naiveKmerCountsCalculation:
    #     Creates a kmers dictionary for the sequence
    #     (with kmers as keys and amount of appearance in the sequence as the values).
    #
    #     return: The kmers dictionary.
    #     """
    #     kmers = self.kmers_count_graph.get_kmers(sequence)
    #
    #     k_mer_counts_dict = {}
    #     # Loop over the string (read/tag)
    #     for kmer in kmers:
    #         # Add the kmer to the dictionary if it's not there or sum +1 to its value
    #         if kmer not in k_mer_counts_dict:
    #             k_mer_counts_dict[kmer] = 0
    #         k_mer_counts_dict[kmer] += 1
    #     # sort the kmer dictionary by amount
    #     #k_mer_counts_dict = collections.OrderedDict(sorted(k_mer_counts_dict.items()))
    #     return k_mer_counts_dict

    def jaccard_similarity(self, kmer_dict):
        """
        jaccard_similarity:
        Calculate the Jaccard similarity of the two kmers dictionary keys
        (where the kmers are keys and their amount of appearance in the sequence are the values).

        :param kmer_dict_a: kmers dictionary.
        :param kmer_dict_b: kmers dictionary.

        :return: float which represent the Jaccard similarity result.
        """

        a = set(kmer_dict.keys())
        t1 = set(self.tag1_kmer_dict.keys())
        t2 = set(self.tag2_kmer_dict.keys())
        direction = ''

        # intersection_t1 = len(a.intersection(t1))
        # intersection_t2 = len(a.intersection(t2))
        jaccard_t1 = len(a.intersection(t1)) / len(a.union(t1))
        jaccard_t2 = len(a.intersection(t2)) / len(a.union(t2))

        if jaccard_t1 >= jaccard_t2:
            direction = '+'
            jaccard = jaccard_t1
        else:
            direction = '-'
            jaccard = jaccard_t2

        return jaccard, direction

    def log_subprocess(self, proc):
        """
        log_subprocess:
        Adds the process console outputs to the log.

        :param proc: The process which created using subprocess.Popen
        """
        output = '\n'
        for line in proc.stdout:
            output += line
        proc.wait()
        logging.info(output)

    def generate_kmers_df(self, read_r1, read_r2, r1_curr_read_kmer_dict, r2_curr_read_kmer_dict, r1_res, r2_res, direction_r1,direction_r2, kmer_df):
        kmer_df.append({'read r1':read_r1[1], 'read r1 id':read_r1[0], 'jaccard results r1': r1_res, 'r1 direction':direction_r1,
                        'read r2':read_r2[1], 'read r2 id':read_r2[0], 'jaccard results r2': r2_res, 'r2 direction':direction_r2})
        # print('aaaaaaa dddfffff ilai ilai',read_r1)
        return kmer_df

    def classify_tag_pair_end_reads(self, r1:str, r2, files_dir, out_dir)->str:
        """
        classify_tag_pair_end_reads:
        Classifies the reads in to two files - reads with tag and reads without tag based on the Jaccard similarity on
        the kmers of the read and the tags.

        :param r1: The fastq file name of r1.
        :param r2: The fastq file name of r2.
        :param i_r1:The index fastq file name of r1
        :param i_r2:The index fastq file name of r2.
        :param files_dir: Directory where the four files are stored.
        :param out_dir:Directory where the output fastq files will be saved.

        :return: Four files names - two files for r1 (tag and without tag) and two files for r2 (tag and without tag).
        """
        logging.info(f'classify reads in {r1} and {r2} to reads with tag and reads without tag using Jacard similarity on the read and the tag kmers (k={self.k})...')

        r1_tag = open(out_dir / r1.replace('.fastq', tag_suffix + '.fastq'), 'w')
        r1_no_tag = open(out_dir / r1.replace('.fastq', no_tag_suffix + '.fastq'), 'w')
        r2_tag = open(out_dir / r2.replace('.fastq', tag_suffix + '.fastq'), 'w')
        r2_no_tag = open(out_dir / r2.replace('.fastq', no_tag_suffix + '.fastq'), 'w')

        kmer_df = []
        for read_r1, read_r2 in zip(self.get_fastq_lines(str(files_dir / r1)), self.get_fastq_lines(str(files_dir / r2))):
            # generating the kmer dict for the read
            r1_curr_read_kmer_dict = self.naiveKmerCountsCalculation(read_r1[1])
            r1_res, direction_r1 = self.jaccard_similarity(r1_curr_read_kmer_dict)

            r2_curr_read_kmer_dict = self.naiveKmerCountsCalculation(read_r2[1])
            r2_res, direction_r2 = self.jaccard_similarity(r2_curr_read_kmer_dict)

            kmer_df = self.generate_kmers_df(read_r1, read_r2, r1_curr_read_kmer_dict, r2_curr_read_kmer_dict, r1_res, r2_res, direction_r1, direction_r2, kmer_df)
            if r1_res > self.jaccard_threshold or r2_res > self.jaccard_threshold:
                for line_r1, line_r2 in zip(read_r1, read_r2):
                    r1_tag.write(line_r1)
                    r2_tag.write(line_r2)
            else:
                for line_r1, line_r2 in zip(read_r1, read_r2):
                    r1_no_tag.write(line_r1)
                    r2_no_tag.write(line_r2)

        r1_tag.close()
        r1_no_tag.close()
        r2_tag.close()
        r2_no_tag.close()
        kmer_df = pd.DataFrame(kmer_df)
        kmer_df.to_csv(out_dir / r1.replace('.fastq','-kmer-validation.csv'))
        return r1.replace('.fastq', tag_suffix + '.fastq'), \
               r1.replace('.fastq', no_tag_suffix + '.fastq'), \
               r2.replace('.fastq', tag_suffix + '.fastq'), \
               r2.replace('.fastq', no_tag_suffix + '.fastq')

    def trim_tag_pair_end_reads(self, r1, r2, files_dir, out_dir):
        '''
        trim_tag_pair_end_reads:
        Trim the tags from the reads (in a pair-end reads experiment).
        '''
        try:
            logging.info('trim tag only in reads that classified as reads with tag...')
            out_r = r1.replace('.fastq', trimmed_suffix + '+-' + '.fastq')
            out_l = r2.replace('.fastq', trimmed_suffix + '+-' + '.fastq')
            r1_out_untrimmed = r1.replace('.fastq', untrimmed_suffix + '.fastq')
            r2_out_untrimmed = r2.replace('.fastq', untrimmed_suffix + '.fastq')
            logging.info('cutadapt summary:')
            #here cutadapt split the outputs to files.
            # cmd = f'~/.local/bin/cutadapt -b {self.tag1} -B {self.tag2} --untrimmed-output {str(out_dir / r1_out_untrimmed)} --untrimmed-paired-output {str(out_dir / r2_out_untrimmed)} ' \
            #       f'-o {str(out_dir / out_r)} -p {str(out_dir / out_l)} ' \
            #       f'{str(files_dir / r1)} {str(files_dir / r2)}'

            # FR paired-end reads.
            cmd = f'~/.local/bin/cutadapt -b "{self.tag1};min_overlap={self.k}" -B "{self.tag2};min_overlap=8" ' \
                  f'-o {str(out_dir / out_r)} -p {str(out_dir / out_l)} ' \
                  f'{str(files_dir / r1)} {str(files_dir / r2)}'
            proc = subprocess.Popen(cmd, shell=True, env=os.environ.copy(), stdout=subprocess.PIPE,
                                    universal_newlines=True, stderr=subprocess.STDOUT)
            self.log_subprocess(proc)

            # RF paired-end reads.
            out_r_final = r1.replace('.fastq', trimmed_suffix + '-+' + '.fastq')
            out_l_final = r2.replace('.fastq', trimmed_suffix + '-+' + '.fastq')

            cmd = f'~/.local/bin/cutadapt -b "{self.tag2};min_overlap={self.k}" -B "{self.tag1};min_overlap={self.k}" ' \
                  f'-o {str(out_dir / out_r_final)} -p {str(out_dir / out_l_final)} ' \
                  f'{str(out_dir / out_r)} {str(out_dir / out_l)}'
            proc = subprocess.Popen(cmd, shell=True, env=os.environ.copy(), stdout=subprocess.PIPE,
                                    universal_newlines=True, stderr=subprocess.STDOUT)
            self.log_subprocess(proc)
            return out_r_final, out_l_final

        except Exception as e:
            logging.error(f'couldnt run cutadapt step (read trimming) on files: {r1} and {r2} please check files paths and make sure cutadapt is install correctly (by typing  ~/.local/bin/cutadapt --version in conlsole).', exc_info=False)
            # sys.exit()
            raise e

def validate_arguments(args):
    """
    validate_arguments:
    Validates the step arguments from the given arguments.

    :param args: The iceberg pipeline arguments.
    """
    try:
        directories_operations.validate_path_existence(args['curr_in_dir'])
        if args['tag1'] == '' or args['tag2'] == '':
            raise ValueError()
        if args['read1_1'] != '':
            for file in [args['read1_1'], args['read2_1']]:
                if file == '':
                    raise ValueError()
                directories_operations.validate_path_existence(Path(args['curr_in_dir'], file))
        if args['read1_2'] != '':
            for file in [args['read1_2'], args['read2_2']]:
                if file == '':
                    raise ValueError()
                directories_operations.validate_path_existence(Path(args['curr_in_dir'], file))
        if args['read1_1'] == '' or args['read1_1'] == '' and args['read1_2'] == '':
            raise ValueError()

    except ValueError as e:
        logging.error('not all necessary arguments for tag step where giving.')
        raise e

def prepare_step(args):
    """
    prepare_step:
    Prepares the step by setting the logging configurations, setting the out directory and
    validates the step arguments.

    :param args: The iceberg pipeline arguments.
    """
    set_logging_config('tag', Path(args['out_dir'], 'DEBUG', 'LOGS', 'tag.log'))
    args['curr_out_dir'] = directories_operations.set_curr_out_dir(args['out_dir'], tag_dir)
    validate_arguments(args)

def main(args):
    try:
        th = tag_handler(args['tag1'], args['tag2'], args['kmer'], args['jaccard_threshold'])
        r1_tag_1, r1_no_tag_1, r2_tag_1, r2_no_tag_1 = th.classify_tag_pair_end_reads(args['read1_1'], args['read2_1'],
                                                                                     # args['index1_1'], args['index2_1'],
                                                                                      args['curr_in_dir'], args['curr_out_dir'])

        r1_tag_trimmed_1, r2_tag_trimmed_1 = th.trim_tag_pair_end_reads(r1_tag_1, r2_tag_1,
                                                                        #args['index1_1'], args['index2_1'],
                                                                        args['curr_out_dir'], args['curr_out_dir'])
        # r1_tag_trimmed_1, r2_tag_trimmed_1 = th.trim_tag_pair_end_reads(args['read1_1'], args['read2_1'],
        #                                                                 args['index1_1'], args['index2_1'],
        #                                                                 args['curr_in_dir'], args['curr_out_dir'])
        args['read1_1'] = r1_tag_trimmed_1
        args['read2_1'] = r2_tag_trimmed_1

        if args['read1_2'] != '':
            r1_tag_2, r1_no_tag_2, r2_tag_2, r2_no_tag_2 = th.classify_tag_pair_end_reads(args['read1_2'], args['read2_2'],
                                                                                          #args['index1_2'], args['index2_2'],
                                                                                          args['curr_in_dir'], args['curr_out_dir'])

            r1_tag_trimmed_2, r2_tag_trimmed_2 = th.trim_tag_pair_end_reads(r1_tag_2, r2_tag_2,
                                                                            #args['index1_2'], args['index2_2'],
                                                                            args['curr_out_dir'], args['curr_out_dir'])
            # r1_tag_trimmed_2, r2_tag_trimmed_2 = th.trim_tag_pair_end_reads(args['read1_2'], args['read2_2'],
            #                                                                 args['index1_2'], args['index2_2'],
            #                                                                 args['curr_in_dir'], args['curr_out_dir'])
            args['read1_2'] = r1_tag_trimmed_2
            args['read2_2'] = r2_tag_trimmed_2
        args['curr_in_dir'] = directories_operations.set_curr_in_dir(args['out_dir'], tag_dir)
        logging.info('tag step done.')
    except Exception as e:
        logging.error('error occurred while running tag step, please read the following for more information.',
                      exc_info=False)
        # sys.exit()
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
    #parser.add_argument('index1_1', help="ff")
    #parser.add_argument('index2_1', help="ff")
    parser.add_argument('--out_dir', help="ff")
    parser.add_argument('--tag1', default='', help="ff")
    parser.add_argument('--tag2', default='', help="ff")
    parser.add_argument('--read1_2', default='', help="ff", metavar='')
    parser.add_argument('--read2_2', default='', help="ff", metavar='')
    #parser.add_argument('--index1_2', default='', help="ff", metavar='')
    #parser.add_argument('--index2_2', default='', help="ff", metavar='')
    parser.add_argument('--kmer', default=8, help="ff", type=int, metavar='')
    parser.add_argument('--jaccard_threshold', default=0.05, type=float, help="ff", metavar='')

    args, rest = parser.parse_known_args()
    return dict([(k, v) for k, v in vars(args).items()])


if __name__ == '__main__':
    args = parse_args()
    print(args)
    print(args['files_dir'])
    args['curr_in_dir'] = Path(args['files_dir'])
    directories_operations.set_curr_out_dir(Path(args['out_dir'], 'DEBUG'), 'LOGS')
    prepare_step(args)

    # print(args['curr_in_dir'])
    main(args)
