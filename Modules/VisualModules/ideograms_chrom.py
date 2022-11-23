from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas as pd
import matplotlib.patches as mpatches

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
####
file_path = '/export/home/fpd/PycharmProjects/kmer_aligner/Data/comparison_all_datasets_07092022.csv'

def bedfile_xcoverage_by_labels(file_path):

    df = pd.read_csv(file_path, sep=',')
    d_probes = {}
    d_coverage_res = {}
    d_coverage_res_new = {}
    d_sv_observed = {}
    for i in range(1,20):
        chrom_str='chr'+str(i)
        print("Working on {}".format(chrom_str))
        df_chr = df[df['chrom']==chrom_str]
        chr_probe_list_x = []
        chr_coverage_res_list_x = []
        chr_coverage_res_new_list_x = []
        chr_observed_svs_res_new_list_x = []
        for index, row in df_chr.iterrows():
            width = row.chromEnd-row.chromStart
            # probes
            chr_probe_list_x.append((row.chromStart,width))
            # res coverage
            if row.res_readCoverage > 0:
                chr_coverage_res_list_x.append((row.chromStart,width))
            # res new coverage
            if row.res_new_readCoverage > 0:
                chr_coverage_res_new_list_x.append((row.chromStart,width))
            # res new observed svs
            if row.res_new_readCoverage > 0 and row.res_new_del_intersect == True:
                chr_observed_svs_res_new_list_x.append(row.chromStart+width/2)

        d_probes[chrom_str]=chr_probe_list_x
        d_coverage_res[chrom_str] = chr_coverage_res_list_x
        d_coverage_res_new[chrom_str] = chr_coverage_res_new_list_x
        d_sv_observed[chrom_str] = chr_observed_svs_res_new_list_x

    #add chrom X
    chrom_str = 'chrX'
    for index, row in df_chr.iterrows():
        width = row.chromEnd-row.chromStart
        # probes
        chr_probe_list_x.append((row.chromStart,width))
        # res coverage
        if row.res_readCoverage > 0:
            chr_coverage_res_list_x.append((row.chromStart, width))
        # res new coverage
        if row.res_readCoverage > 0:
            chr_coverage_res_new_list_x.append((row.chromStart, width))
        # res new observed svs
        if row.res_new_readCoverage > 0 and row.res_new_del_intersect == True:
            chr_observed_svs_res_new_list_x.append(row.chromStart + width / 2)

    d_probes[chrom_str] = chr_probe_list_x
    d_coverage_res[chrom_str] = chr_coverage_res_list_x
    d_coverage_res_new[chrom_str] = chr_coverage_res_new_list_x
    d_sv_observed[chrom_str] = chr_observed_svs_res_new_list_x

    return d_probes, d_coverage_res, d_coverage_res_new, d_sv_observed


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


probes_x_dict, coverage_res_dict, coverage_res_new_dict, observed_svs_res_new_dict = bedfile_xcoverage_by_labels(file_path)


fig = plt.figure()
ax = fig.add_subplot(111)
d = {}
yticks = []
yticklabels = []

# ideogram.txt downloaded from UCSC's table browser
y=0
for xranges, yrange, colors, label in ideograms('/export/home/fpd/PycharmProjects/kmer_aligner/Data/SRA_SRR1561378/mouse_genome_data/mouse_mm9_table_ordered.txt'):
    coll = BrokenBarHCollection(xranges, yrange, facecolors=colors)
    ax.add_collection(coll)
    center = yrange[0] + yrange[1]/2.
    yticks.append(center)
    yticklabels.append(label)
    d[label] = xranges

    # Add Probes
    print(label)
    if label in probes_x_dict:
        # probes
        chr_x_probes = probes_x_dict[label]
        probes = BrokenBarHCollection(chr_x_probes, (yrange[0]+yrange[1],0.3),linewidths=(5,),color='springgreen',
                                      alpha=0.7)
        ax.add_collection(probes)
        # res coverage
        # chr_x_res_coverage = coverage_res_dict[label]
        # res_coverage = BrokenBarHCollection(chr_x_res_coverage, (yrange[0]+yrange[1]+0.2 ,0.2), linewidths=(8,),color='green',alpha=0.7)
        # ax.add_collection(res_coverage)
        #res new coverage
        chr_x_res_new_coverage = coverage_res_new_dict[label]
        res_new_coverage = BrokenBarHCollection(chr_x_res_new_coverage, (yrange[0]+yrange[1]+0.3, 0.3),
                                                linewidths=(5,),color='cornflowerblue', alpha=0.7,
                                                )
        ax.add_collection(res_new_coverage)
        # res new observed sv's
        chr_x_res_new_observed_svs = observed_svs_res_new_dict[label]
        svs_y_values = [yrange[0]+yrange[1]+0.6] * len(chr_x_res_new_observed_svs)
        ax.scatter(chr_x_res_new_observed_svs,svs_y_values,marker='x',color='red',linewidths=(1,))

    # Add Translocation marks
    trans_chr1_x_vals = [124187037,159768112,168036729,83379183,4313197,181813381]
    trans_chr1_y_vals = [62.3,62.3,62.3,62.3,62.3,62.3]
    ax.scatter(trans_chr1_x_vals, trans_chr1_y_vals, marker='x', color='purple', linewidths=(1,))

fig.suptitle("O'Geen et al. CaptureSeq Experiment Probes, Coverage & SV Observations", fontsize=16)
# fig.suptitle("O'Geen et al. CaptureSeq Experiment Probes", fontsize=16)
ax.axis('tight')
ax.set_yticks(yticks)
ax.set_yticklabels(yticklabels)
ax.set_xticks([])
legend_handles = []
legend_handles.append(mpatches.Patch(color='springgreen', label='Capture Probs'))
legend_handles.append(mpatches.Patch(color='cornflowerblue', label='Reads Coverage (Suggested Configuration)'))
legend_handles.append(mpatches.Patch(color='red', label='Novel Observed SVs'))
legend_handles.append(mpatches.Patch(color='purple', label='SVs as reported in Section 4.3'))
ax.legend(handles=legend_handles,loc='center right')

plt.show()