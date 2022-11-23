import pickle
import pandas as pd
from Modules import analysis_module

df_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/rag1_translocation_df_23-10-2021_22-42.pkl'
result_list_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/classify_reads_result_23-10-2021_22-43.pkl'

df_bed = pd.read_csv('/export/home/fpd/PycharmProjects/kmer_aligner/Data/RAG1.bed', sep='\t', comment='t', header=None)
header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
df_bed.columns = header[:len(df_bed.columns)]

# analysis_obj = analysis_module.AnalysisClass()
# region_to_chrom_map = analysis_obj.get_chromosomes_from_regions_map(df_bed)
with open(result_list_path, 'rb') as f:
    result_list = pickle.load(f)

gt_df = pd.read_pickle(df_path)
gt_df['predicted'] = result_list

print("IM HERE")
print("Done")