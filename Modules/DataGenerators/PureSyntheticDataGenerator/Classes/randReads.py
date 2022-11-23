import pickle
import pandas as pd
import numpy as np
from scipy import stats

class randReads(object):
    def __init__(self):
        self.min_reads = 300
        self.iceberg_df = None
        self.sam_df = None
        self.sonication_mean = 170
        self.sonication_std = 25

    def set_sonication_mean(self, mean):
        self.sonication_mean = mean

    def get_sonication_mean(self):
        return self.sonication_mean

    def set_sonication_std(self, std):
        self.sonication_std = std

    def get_sonication_std(self):
        return self.sonication_std

    def load_iceberg_dataframe_from_pickle(self, pkl_path):
        '''
        :param pkl_path: Path to dam DF
        :return: None
        fills self.iceberg_df with the DF from the pickle
        '''
        self.iceberg_df = pickle.load(open(pkl_path, 'rb'))
        # Type check
        assert type(self.iceberg_df) == pd.core.frame.DataFrame


    def load_sam_dataframe_from_pickle(self, pkl_path):
        '''
        :param pkl_path: Path to dam DF
        :return: None
        fills self.dam_df with the DF from the pickle
        '''
        self.sam_df = pickle.load(open(pkl_path, 'rb'))
        # Type check
        assert type(self.sam_df) == pd.core.frame.DataFrame

    def rand_iceberg(self):
        '''
        :return: DataFrame series (DF row) of Iceberg that fulfills the min num of counts
        '''
        return self.iceberg_df[self.iceberg_df['Count'] > self.min_reads].sample(1).iloc[0]

    def rand_reads_from_iceberg(self, iceberg_row, size=1):
        read_right = np.random.choice([r[0]-iceberg_row['Peak'] for r in iceberg_row['Reads']], size = size)
        sonication = stats.norm.rvs(self.sonication_mean, self.sonication_std, size=size).astype(int)
        read_left = np.add(read_right, sonication)
        return list(zip(read_right, read_left))