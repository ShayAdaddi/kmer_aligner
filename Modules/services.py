class Services:
    def __init__(self):
        return

    def get_fastq_lines(self, file):
        fastq = open(file, 'r')
        with fastq as f:
            while True:
                # tic0 = time.perf_counter()
                l1 = f.readline() # read name
                if not l1:
                    break
                l2 = f.readline()  # read sequence
                l3 = f.readline()
                l4 = f.readline()  # read quality
                # tic1 = time.perf_counter()
                # print("get_fastq_lines Timing {} (sec)".format(tic1-tic0))
                yield l2


    def reverse_complement(self,seq):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '\n': ''}
        new = ''.join([complement[base] for base in seq[::-1]])
        return new


    def fetch_randomly_from_bedfile(self, chr_num, bed_df):
        # Random choice from entry of specific chr
        # If chr number does not exists in the bed file return -1
        # header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
        chr_regions_random_sample = bed_df[bed_df.chrom == "chr{}".format(chr_num)].sample()
        start = chr_regions_random_sample.chromStart
        end = chr_regions_random_sample.chromEnd
        return start, end


    def generate_regions_sequences_list(self,path):
        regions_raw = []
        with open(path, 'r') as f:
            for region_idx, region in enumerate(f.readlines()):
                regions_raw.append(region)
        return regions_raw