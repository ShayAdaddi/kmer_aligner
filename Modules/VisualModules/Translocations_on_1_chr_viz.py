from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas as pd
import matplotlib.patches as mpatches
from matplotlib.pyplot import cm
import numpy as np

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
height = 0.9
spacing = 2

file_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/SRA_SRR1561378/mouse_genome_data/Ogeen_trans_chr1_list.csv' #Ogeen
# file_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/SRA_ERP001058_Liver_cancer_human/Li_trans_chr1_list.csv' #Li
df_trans = pd.read_csv(file_path, sep=',')

def ideograms(fn):
    last_chrom = None
    fin = open(fn)
    fin.readline()
    xranges, colors = [], []
    ymin = 0

    for line in fin:
        chrom, start, stop, label, stain = line.strip().split('\t')
        start = int(start)
        stop = int(stop)
        width = stop - start
        if chrom == last_chrom or (last_chrom is None):
            xranges.append((start, width))
            colors.append(color_lookup[stain])
            last_chrom = chrom
            continue

        ymin += height + spacing
        yrange = (ymin, height)
        yield xranges, yrange, colors, last_chrom
        xranges, colors = [], []
        xranges.append((start, width))
        colors.append(color_lookup[stain])
        last_chrom = chrom

    # last one
    ymin += height + spacing
    yrange = (ymin, height)
    yield xranges, yrange, colors, last_chrom

fig = plt.figure()
ax = fig.add_subplot(111)
d = {}
yticks = []
yticklabels = []

# ideogram.txt downloaded from UCSC's table browser
for xranges, yrange, colors, label in ideograms('/export/home/fpd/PycharmProjects/kmer_aligner/Data/SRA_SRR1561378/mouse_genome_data/mouse_mm9_table_chr1.txt'): #Ogeen
# for xranges, yrange, colors, label in ideograms('/export/home/fpd/PycharmProjects/kmer_aligner/Data/SRA_ERP001058_Liver_cancer_human/human_hg38_table_chr1.txt'):
    coll = BrokenBarHCollection(xranges, yrange, facecolors=colors)
    ax.add_collection(coll)
    center = yrange[0] + yrange[1]/2.
    yticks.append(center)
    yticklabels.append(label)
    d[label] = xranges

legend_patches_handles = []
n=6 #Ogeen
# n=10 # Li
color = iter(cm.Paired(np.linspace(0, 1, n)))
for index, row in df_trans.iterrows():
    c = next(color)
    x_1 = [(row.start, row.end - row.start)]
    y_1 = (3.9+index*0.2, 0.1)
    print(y_1)
    text_offset = 100000
    text_angle = 15
    ax.broken_barh(x_1, y_1, linewidths=(5,), color=c, alpha=0.7)
    ax.text(row.start, 4+index*0.2, str(row.start), bbox=dict(facecolor='white', alpha=0.1),fontsize=8,rotation=text_angle)
    if (row.end-text_offset) - row.start > 100000:
        ax.text(row.end - text_offset, 4 + index * 0.2, str(row.end), bbox=dict(facecolor='white', alpha=0.1),
                fontsize=8, rotation=text_angle)
    else:
        ax.text(row.end - text_offset, 3.85 + index * 0.2, str(row.end), bbox=dict(facecolor='white', alpha=0.1),
                fontsize=8, rotation=text_angle)
    legend_patches_handles.append(mpatches.Patch(color=c, label=row.read_count))

fig.suptitle("O'Geen et al. Potential Translocation/Long Deletions Events", fontsize=16)
# fig.suptitle("Li Y et al. Potential Translocation/Long Deletions Events", fontsize=16)
ax.axis('tight')
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels)
ax.set_xticks([])

plt.legend(handles=legend_patches_handles,title="read count")

plt.show()