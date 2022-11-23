import pandas as pd
from pyfaidx import Fasta
import khmer
from Modules import kmer_operations


class OfflineClass:

    def init(self):

        # DEFAULTS
        ko = kmer_operations.KmerOperationsClass()
        self.k = 4
        self.alphabet = ['A', 'C', 'G', 'T', 'N']
        self.cg = khmer.Countgraph(self.k, len(self.alphabet) ** self.k, 1)
        # INITIALIZATION
        self.output_folder = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/'
        self.df = self.generate_bed_file_df('/export/home/fpd/PycharmProjects/kmer_aligner/Data/RAG1.bed')
        self.full_genome = self.generate_genome_object('/export/home/fpd/PycharmProjects/kmer_aligner/Data/hg38.fa')
        self.LUT = ko.generate_kmer_permutations_dictionary(self.k)
        self.regions_extractor_and_files_generator(self.df, self.full_genome)
        self.virtual_genome_generator(self.df, self.full_genome)

        return

    def generate_genome_object(self, genome_path):
        genome = Fasta(genome_path, sequence_always_upper=True)
        return genome

    def generate_bed_file_df(self, bedfile_path):
        # Generate Bedfile DF
        df = pd.read_csv(bedfile_path, sep='\t', comment='t', header=None)
        header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
        df.columns = header[:len(df.columns)]
        return df

    def virtual_genome_generator(self, bed_df, genome):
        chunks_size = 50  # Size of each row in genome
        buffer_size = 0  # Amount of A's to add between regions
        expansion_size = 200  # size of region expansion to build virtual genome
        with open(self.output_folder + 'virtualReferenceGenome.fa', 'w') as f1:
            for index, row in bed_df.iterrows():
                print('Working on Region {}'.format(index))
                ref_seq = genome[row.chrom][row.chromStart - expansion_size:row.chromEnd + expansion_size].seq
                chunks = [ref_seq[i:i + chunks_size] for i in range(0, len(ref_seq), chunks_size)]
                f1.write(">region_" + str(index) + "\n")
                if buffer_size != 0: f1.write(buffer_size * 'A' + '\n')  # padding
                for c in chunks:
                    f1.write(c + "\n")
                f1.write("\n\n")
            f1.close()

        return

    def regions_extractor_and_files_generator(self, bed_df, genome):
        with open(self.output_folder + 'regionsKmersOptimized_k4.txt', 'w') as f1, \
                open(self.output_folder + 'regionsKmersRaw_k4.txt', 'w') as f2, \
                open(self.output_folder + 'regionsRawSequence_k4.txt', 'w') as f3:
            for index, row in bed_df.iterrows():
                print('Working on Iceberg {}'.format(index))
                ref_seq = genome[row.chrom][row.chromStart:row.chromEnd].seq
                # RAW SEQUENCE
                f3.write(str(ref_seq))
                f3.write("\n")
                # RAW KMERS
                kmers_dict = self.ko.kmer_counter_calculation(ref_seq, self.k)
                l = self.ko.kmerToList(kmers_dict)
                f2.write(str(l))
                f2.write("\n")
                # OPTIMIZED KMERS (ONLY UNIQUE)
                l2 = self.ko.kmer_to_sorted_unique_idx_list(kmers_dict, self.LUT[1])
                f1.write(str(l2))
                f1.write("\n")

            f1.close()
            f2.close()
            f3.close()

        return


def main():
    try:
        offline_obj = OfflineClass()
        offline_obj.init()

    except Exception as e:
        raise e


if __name__ == '__main__':
    main()
