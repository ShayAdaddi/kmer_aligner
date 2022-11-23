import pandas as pd

ROOT_DIR = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/SRA_SRR1561378/mouse_genome_data/'
bedfile_path = ROOT_DIR + 'res_sorted_bed_coverage.bed'
del_obesrvation_path = ROOT_DIR + 'ogeen_del_new_mapq30_min30_max99.txt'
del_header = ['chrom', 'location', 'mapq', 'cigar']
trans_observation_path = ROOT_DIR + 'Ogeen_translocation.bed'
trans_header = ['chrom', 'start', 'end', 'str', 'b', 'c']

def generate_bed_file_df(bedfile_path):
    # Generate Bedfile DF
    df = pd.read_csv(bedfile_path, sep='\t', comment='t', header=None)
    header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand','readCoverage']
    df.columns = header[:len(df.columns)]
    return df

def load_observation(obsfile_path,obs_header):
    df_obs = pd.read_csv(obsfile_path, sep=' ', header=None)
    df_obs.columns = obs_header[:len(df_obs.columns)]
    return df_obs

def load_observation_tab(obsfile_path,obs_header):
    df_obs = pd.read_csv(obsfile_path, sep='\t', comment='t', header=None)
    df_obs.columns = obs_header[:len(df_obs.columns)]
    return df_obs

def main():
    try:
        bedfile_df = generate_bed_file_df(bedfile_path)

        bedfile_df['del_intersect'] = False
        del_df = load_observation(del_obesrvation_path,del_header)

        # Going over all observations
        for index, row in del_df.iterrows():
            if index % 10000 == 0:
                print('Working on observation # {}'.format(index))
            obs_chr, obs_loc = row.chrom, row.location
            # check if obs_loc in range of relevant chrom in bed file (chrom==obs_chr)
            # chr_bedfile = bedfile_df[bedfile_df['chrom']== obs_chr]
            # Example: obs_chr, obs_loc = 'chr10' , 5734116
            mask = (bedfile_df['chrom']==obs_chr) & (bedfile_df['chromEnd'] > obs_loc) & (bedfile_df['chromStart'] < obs_loc)
            bedfile_df.loc[mask, 'del_intersect'] = True
            #DEBUG
            if mask.any():
                print("{} location={} cigar={}".format(row.chrom,row.location,row.cigar))

        bedfile_df['trans_intersect'] = False
        trans_df = load_observation_tab(trans_observation_path, trans_header)
        for index, row in trans_df.iterrows():
            if index % 10000 == 0:
                print('Working on observation # {}'.format(index))
            obs_chr, obs_loc = row.chrom, row.start
            # check if obs_loc in range of relevant chrom in bed file (chrom==obs_chr)
            # chr_bedfile = bedfile_df[bedfile_df['chrom']== obs_chr]
            # Example: obs_chr, obs_loc = 'chr10' , 5734116
            mask = (bedfile_df['chrom']==obs_chr) & (bedfile_df['chromEnd'] > obs_loc) & (bedfile_df['chromStart'] < obs_loc)
            bedfile_df.loc[mask, 'trans_intersect'] = True

        # Export
        # bedfile_df.to_excel(ROOT_DIR + "ogeen_obs_res.xlsx")

    except Exception as e:
        raise e


if __name__ == '__main__':
    main()
