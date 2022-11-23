import sys

from pyfaidx import Fasta
from Classes import randomReadsGeneratorClass
import pandas as pd
import pickle
from datetime import datetime
import os
os.chdir('../../../')
root_dir = os. getcwd() + '/Data/'

# Hyper parameters and Macros
genome_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/hg38.fa'
genome = Fasta(genome_path,sequence_always_upper=True)
bedfile_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/RAG1_-Capture-seq-List-CRISPR-regions-hg38-final.bed'
read_size = 140
num_of_reads_per_region = 10
today = datetime.today()
timestamp = today.strftime("%d-%m-%Y_%H-%M")
r1_output_name = 'r1_synthetic_data_' + timestamp + '.fq'
r2_output_name = 'r2_synthetic_data_' + timestamp + '.fq'
excel_output_name = 'df_synthetic_data_' + timestamp + '.xlsx'
pickle_output_name = 'df_synthetic_data_' + timestamp + '.pkl'
folder_name = 'data_set_' + timestamp +'/'
# make dataset dir
os.mkdir(root_dir + folder_name)
sys.stdout = open(root_dir + folder_name + 'readme.txt','w')

def main():
    re_obj = randomReadsGeneratorClass.RandomReadsGeneratorClass(bedfile_path,genome,read_size)
    r1,r2,info_df = re_obj.rand_reads_for_all_regions2()
    info_df.to_excel(root_dir + folder_name + excel_output_name)
    info_df.to_pickle(root_dir + folder_name + pickle_output_name)
    re_obj.write_fastq_file2(root_dir + folder_name + r1_output_name, r1)
    re_obj.write_fastq_file2(root_dir + folder_name + r2_output_name, r2)

if __name__ == '__main__':
    main()
    sys.stdout.close()