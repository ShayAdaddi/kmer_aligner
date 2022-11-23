import pandas as pd
import pysam

ROOT_DIR = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/SRA_SRR1561378/mouse_mock/'
bedfile_path = ROOT_DIR + 'bait-100-50-mouse-selection-stringent.bed'
bwa_filename = "ogeen_res_new_fixmate_sorted_dedup.bam"
bwa_path = ROOT_DIR + bwa_filename

def generate_bed_file_df(bedfile_path):
    # Generate Bedfile DF
    df = pd.read_csv(bedfile_path, sep='\t', comment='t', header=None)
    header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand','readCoverage']
    df.columns = header[:len(df.columns)]
    return df

def main():
    try:
        bedfile_df = generate_bed_file_df(bedfile_path)
        bedfile_df['read_count'] = 0
        print(bedfile_df.head())

        # go over all reads
        bwa_res = pysam.AlignmentFile(bwa_path, "rb")
        read_index=0
        for al in bwa_res.fetch():
            if read_index % 1000 == 0:
                print("processing read {}".format(read_index))
            read_start = al.reference_start
            read_end = al.reference_end
            chr = al.reference_name
            logical_index = (bedfile_df['chrom']==chr) & ((read_start>=bedfile_df['chromStart']) & (read_start<=bedfile_df['chromEnd'])) | ((read_end>=bedfile_df['chromStart']) & (read_end<=bedfile_df['chromEnd']))
            if logical_index.any():
                for row,index in bedfile_df.iterrows():
                    if row.read_counts > 0:
                        print("{}".format(row.chrom))
                break
            read_index += 1

    except Exception as e:
        raise e


if __name__ == '__main__':
    main()
