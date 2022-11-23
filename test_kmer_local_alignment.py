from Modules import kmer_operations

##### Kmer local alignment tests ######
K=8

ko = kmer_operations.KmerOperationsClass()
a = 'CCCCGGGGTTGGCCAC'
b = 'AAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCGGGGGGGGGGGGGGACGTNAATTGGCCACGTNAGCTN'
ko.kmer_local_alignment_debug(a,b,K)

a_kmers = ko.getKmersList(a,K)
b_kmers = ko.getKmersList(b,K)
ko.kmer_local_alignment2(a_kmers,b_kmers)