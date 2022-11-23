import numpy as np
import pylab as pl
from matplotlib import collections  as mc
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import gridspec

# file_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/comparison_all_datasets_07092022.csv'
# df = pd.read_csv(file_path, sep=',')
#
# i=1
# all_chr_list_x = []
# all_chr_list_y = []
# for i in range(1,20):
#     chrom_str='chr'+str(i)
#     print("Working on {}".format(chrom_str))
#     df_chr = df[df['chrom']==chrom_str]
#     chr_segments_list_x = []
#     chr_segments_list_y = []
#     for index, row in df_chr.iterrows():
#         chr_segments_list_x.append((i,i))
#         chr_segments_list_y.append((row.chromStart,row.chromEnd))
#     all_chr_list_x.append(chr_segments_list_x)
#     all_chr_list_y.append(chr_segments_list_y)
#     i+=1
#
#
# fig=plt.figure()
# gs = gridspec.GridSpec(4, 1)
# fig.suptitle('Chromosome 1')
# chr = 1
# ax1 = plt.subplot(gs[0])
# ax1.plot((1,10000),(1,1),linewidth=7.0,color='r', label='chr1')
# plt.setp(ax1.get_xticklabels(), visible=False)
# ax2 = plt.subplot(gs[1],sharex=ax1)
# ax2.plot((3,5000),(1,1),(7098,8088),(1,1))
# plt.setp(ax2.get_xticklabels(), visible=False)
# ax3 = plt.subplot(gs[2],sharex=ax1)
# ax3.scatter(5000,2)
# plt.setp(ax3.get_xticklabels(), visible=False)
# ax4 = plt.subplot(gs[3],sharex=ax1)
# ax4.bar(4980,50,width=30)
# plt.subplots_adjust(hspace=.0)


# for chr in range(0,len(all_chr_list_x)):
#     plt.plot(all_chr_list_x[chr],all_chr_list_y[chr],linewidth=7.0,color='r', label='chromosome')
#
# plt.xticks(np.arange(1, 20, step=1))
#
# plt.show()


















"""
Demonstrates plotting chromosome ideograms and genes (or any features, really)
using matplotlib.
1) Assumes a file from UCSC's Table Browser from the "cytoBandIdeo" table,
saved as "ideogram.txt". Lines look like this::
    #chrom  chromStart  chromEnd  name    gieStain
    chr1    0           2300000   p36.33  gneg
    chr1    2300000     5300000   p36.32  gpos25
    chr1    5300000     7100000   p36.31  gneg
2) Assumes another file, "ucsc_genes.txt", which is a BED format file
   downloaded from UCSC's Table Browser. This script will work with any
   BED-format file.
"""

from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas


# Here's the function that we'll call for each dataframe (once for chromosome
# ideograms, once for genes).  The rest of this script will be prepping data
# for input to this function
#
def chromosome_collections(df, y_positions, height,  **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        print(chrom)
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']


# Height of each ideogram
chrom_height = 1

# Spacing between consecutive ideograms
chrom_spacing = 1

# Height of the gene track. Should be smaller than `chrom_spacing` in order to
# fit correctly
gene_height = 0.4

# Padding between the top of a gene track and its corresponding ideogram
gene_padding = 0.1

# Width, height (in inches)
figsize = (6, 8)

# Decide which chromosomes to use
chromosome_list = ['chr{}'.format(i) for i in range(1, 20)]

# Keep track of the y positions for ideograms and genes for each chromosome,
# and the center of each ideogram (which is where we'll put the ytick labels)
ybase = 0
chrom_ybase = {}
gene_ybase = {}
chrom_centers = {}

# Iterate in reverse so that items in the beginning of `chromosome_list` will
# appear at the top of the plot
for chrom in chromosome_list[::-1]:
    chrom_ybase[chrom] = ybase
    chrom_centers[chrom] = ybase + chrom_height / 2.
    gene_ybase[chrom] = ybase - gene_height - gene_padding
    ybase += chrom_height + chrom_spacing

# Read in ideogram.txt, downloaded from UCSC Table Browser
ideo = pandas.read_table(
    '/export/home/fpd/PycharmProjects/kmer_aligner/Data/SRA_SRR1561378/mouse_genome_data/mouse_mm9_table.txt',
    skiprows=1,
    names=['chrom', 'start', 'end', 'name', 'gieStain']
)

# Filter out chromosomes not in our list
ideo = ideo[ideo.chrom.apply(lambda x: x in chromosome_list)]

# Add a new column for width
ideo['width'] = ideo.end - ideo.start

# Colors for different chromosome stains
color_lookup = {
                  'gneg': (.7, .7, .7),
                'gpos25': (.6, .6, .6),
                'gpos33': (.5, .5, .5),
                'gpos50': (.4, .4, .4),
                'gpos66': (.3, .3, .3),
                'gpos75': (.2, .2, .2),
               'gpos100': (.0, .0, .0),
                  'acen': (.8, .4, .4),
                  'gvar': (.8, .8, .8),
                 'stalk': (.9, .9, .9),
               }

# Add a new column for colors
ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x])


# Same thing for genes
genes = pandas.read_table(
    '/export/home/fpd/PycharmProjects/kmer_aligner/Data/SRA_SRR1561378/mouse_genome_data/mouse_mm9_gene_table.txt',
    names=['chrom', 'start', 'end', 'name'],
    usecols=range(4))
genes = genes[genes.chrom.apply(lambda x: x in chromosome_list)]
genes['width'] = genes.end - genes.start
genes['colors'] = '#2243a8'


fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111)

# Now all we have to do is call our function for the ideogram data...
print("adding ideograms...")
for collection in chromosome_collections(ideo, chrom_ybase, chrom_height):
    ax.add_collection(collection)

# ...and the gene data
print("adding genes...")
for collection in chromosome_collections(
    genes, gene_ybase, gene_height, alpha=0.5, linewidths=0
):
    ax.add_collection(collection)

print("done")
# Axes tweaking
ax.set_yticks([chrom_centers[i] for i in chromosome_list])
ax.set_yticklabels(chromosome_list)
ax.axis('tight')
plt.show()