import pandas as pd
import random
from Bio.Seq import Seq
import numpy as np


class Read:

    def __init__(self, name, sequence, strand, quality):
        self.name = name
        self.sequence = sequence
        self.strand = strand
        self.quality = quality

    def get_name(self):
        return self.name

    def get_sequence(self):
        return self.sequence

    def get_strand(self):
        return self.strand

    def get_quality(self):
        return self.quality


class RandomReadsGeneratorClass:

    def __init__(self, bed_file_path, genome, read_size):
        self.bed_df = self.generate_bed_file_df(bed_file_path)
        self.genome = genome
        self.read_size = read_size
        self.del_start_low = 20
        self.del_start_high = 50
        self.del_len_low = 70
        self.del_len_high = 90
        self.mu = 300
        self.sigma = 50
        self.max_sonication_size = 450

        # Print to readme
        print("Test Parameters:\ndel_start_low = {}\ndel_start_high = {}\ndel_len_low = {}\ndel_len_high = {}\nmu = {}\nsigma = {}\nmax_sonication_size = {}\n ".format(self.del_start_low,self.del_start_high,self.del_len_low,self.del_len_high,self.mu,self.sigma,self.max_sonication_size))

    def get_rand_del_start(self):
        return random.randrange(self.del_start_low, self.del_start_high)

    def get_rand_del_length(self):
        return random.randrange(self.del_len_low, self.del_len_high)

    def get_sonication_segment_size(self):
        return int(min(np.random.normal(self.mu, self.sigma, 1), self.max_sonication_size))

    def get_bed_df(self):
        return self.bed_df

    def generate_bed_file_df(self, bed_file_path):
        # Generate Bedfile DF
        df = pd.read_csv(bed_file_path, sep='\t', comment='t', header=None)
        header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
        df.columns = header[:len(df.columns)]
        return df

    def extract_region(self, chr, start, end):
        return self.genome[chr][start:end]

    def rand_read_from_region(self, region, do_deletion):
        read_pos = random.randrange(self.read_size, len(region.seq) - self.read_size)
        read = region[read_pos:read_pos + self.read_size]
        read_rc = read.reverse.complement
        if do_deletion:
            deletion_start = random.randrange(self.del_start_low, self.del_start_high)
            deletion_length = random.randrange(self.del_len_low, self.del_len_high)
            new_read = read.seq[:deletion_start] + read.seq[deletion_start + deletion_length:]
            seq = Seq(new_read)
            new_read_rc = seq.reverse_complement()
            read_info = pd.Series(
                {'read_org': read.seq, 'new_read': new_read, 'read_org_rc': read_rc.seq,
                 'new_read_rc': new_read_rc, 'read_pos': read_pos, 'del_start_pos': read_pos + deletion_start,
                 'del_len': deletion_length})
            return new_read, new_read_rc, read_info

        else:
            read_info = pd.Series(
                {'read_org': read.seq, 'new_read': read.seq, 'read_org_rc': read_rc.seq,
                 'new_read_rc': read_rc.seq, 'read_pos': read_pos, 'del_start_pos': -1, 'del_len': -1})
            return read.seq, read_rc.seq, read_info

    def rand_read_from_region2(self, region, do_deletion):

        # r1 - right to left
        # r2 - left to right

        region_del_pos_start = random.randrange(3 * self.read_size, len(region) - 3 * self.read_size)
        if do_deletion:
            rand_del_length = self.get_rand_del_length()
        else:
            rand_del_length = 0

        new_region = region.seq[:region_del_pos_start] + region.seq[(region_del_pos_start + rand_del_length):]  # region with dels
        sonication_size = self.get_sonication_segment_size()
        sonication_del_pos = random.randrange(0, sonication_size)

        r2_start = region_del_pos_start - sonication_del_pos
        r2 = new_region[r2_start:r2_start + self.read_size]

        r1_start = region_del_pos_start + sonication_size - sonication_del_pos - self.read_size
        r1 = new_region[r1_start:r1_start + self.read_size]
        seq = Seq(r1)
        r1_rc = seq.reverse_complement()

        # Extract info

        read_info = pd.Series(
            {'region_del_pos_start': region_del_pos_start, 'del_length': rand_del_length,
             'sonication_size': int(sonication_size), 'sonication_del_pos': sonication_del_pos})

        return r1_rc, r2, read_info

    def rand_reads_for_all_regions2(self):

        r1 = []
        r2 = []
        reads_df = pd.DataFrame(columns=['index', 'region_label', 'region_start', 'region_end', 'region_del_pos_start',
                                         'del_length', 'sonication_size', 'sonication_del_pos'])

        for indx, row in self.bed_df.iterrows():
            region = self.extract_region(row.chrom, row.chromStart, row.chromEnd)  # Extract from Genome
            print("Working on region {}".format(row.name))

            read_r1, read_r2, read_info = self.rand_read_from_region2(region, True)
            series = pd.concat([pd.Series({'index': indx, 'region_label': row.name, 'region_start': row.chromStart, 'region_end': row.chromEnd}),
                                read_info])
            reads_df.loc[indx] = series

            # Extract info in region coordinates
            del_start_in_org_region = row.chromStart + read_info['region_del_pos_start']
            del_end_in_org_region = del_start_in_org_region + read_info['del_length']
            sonic_start_in_org_region = del_start_in_org_region - read_info['sonication_del_pos']
            sonic_end_in_region = del_end_in_org_region + read_info['sonication_size'] - read_info['sonication_del_pos']

            # Print information
            print("region label = {} chr = {} region_start = {} region_end = {}".format(row.name,row.chrom,row.chromStart,row.chromEnd))
            print("del_len = {} sonication_size = {}".format(read_info['del_length'],read_info['sonication_size']))
            print("sonication_start = {} del_start = {} del_end = {} sonication_end = {}".format(sonic_start_in_org_region,del_start_in_org_region,del_end_in_org_region,sonic_end_in_region))
            print("r1 = {}".format(read_r1))
            print("r2 = {}".format(read_r2))


            # Generate reads
            read_r1_obj = Read("@region_{}:{}:sonic_start_pos={}:del_start_pos={}:del_end_pos={}:sonic_end_pos={}:read{}/{}".format(row.name,row.chrom,sonic_start_in_org_region,del_start_in_org_region,del_end_in_org_region,sonic_end_in_region,indx,1), read_r1, '+', 'B' * len(read_r1))
            # read_r1_obj = Read("@READ_{}/{}".format(indx,1), read_r1, '+', 'B' * len(read_r1))

            print('----------------------')
            # print("del_length={}, sonication_length={}".format(read_info['del_length'],read_info['sonication_size']))
            # print(read_r1_obj.name)
            read_r2_obj = Read("@region_{}:{}:sonic_start_pos={}:del_start_pos={}:del_end_pos={}:sonic_end_pos={}:read{}/{}".format(row.name,row.chrom,sonic_start_in_org_region,del_start_in_org_region,del_end_in_org_region,sonic_end_in_region,indx,2), read_r2, '+', 'B' * len(read_r2))
            # read_r2_obj = Read("@READ_{}/{}".format(indx,2), read_r2, '+', 'B' * len(read_r2))

            # print(read_r2_obj.name)
            print('----------------------')
            r1.append(read_r1_obj)
            r2.append(read_r2_obj)

        return r1, r2, reads_df

    def rand_reads_for_all_regions(self, number_of_reads_per_region):
        r1 = []
        r2 = []
        reads_df = pd.DataFrame(columns=['index', 'read_name', 'chr', 'region_label', 'region_start', 'read_org',
                                         'new_read', 'read_org_rc', 'new_read_rc', 'read_pos', 'del_start_pos',
                                         'del_len'])
        for indx, row in self.bed_df.iterrows():
            region = self.extract_region(row.chrom, row.chromStart, row.chromEnd)
            print("Working on region {}".format(row.name))
            i = 0
            while i < number_of_reads_per_region:
                # if i % 100 == 0: print("Working on read {}".format(i))

                read, read_rc, read_info = self.rand_read_from_region(region, True)

                r1.append(read)
                r2.append(read_rc)

                series = pd.concat([pd.Series({'index': indx * number_of_reads_per_region + i,
                                               'read_name': "@/{}{}".format(1, read[0:20]), 'chr': row.chrom,
                                               'region_label': row.name, 'region_start': row.chromStart}),
                                    read_info])
                reads_df.loc[indx * number_of_reads_per_region + i] = series

                i += 1

        return r1, r2, reads_df

    def write_fastq_file(self, path, reads_list, lib_number):

        with open(path, 'w+') as f:
            for item in reads_list:
                f.writelines("@/{}{}\n".format(lib_number, item[0:20]))
                f.writelines("{}\n".format(item))
                f.writelines('+\n')
                f.writelines('B' * len(item))
                f.writelines('\n')
        return

    def write_fastq_file2(self, path, reads_obj_list):

        with open(path, 'w+') as f:
            for item in reads_obj_list:
                f.writelines("{}\n".format(item.get_name()))
                f.writelines("{}\n".format(item.get_sequence()))
                f.writelines('{}\n'.format(item.get_strand()))
                f.writelines('{}\n'.format(item.get_quality()))
        return
