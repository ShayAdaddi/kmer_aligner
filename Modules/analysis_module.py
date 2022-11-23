import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
import pandas as pd


class AnalysisClass:
    def __init__(self):
        return

    def calculate_confusion_matrix(self,df):
        plt.figure(figsize=(18, 8))
        print('Plotting Confusion Matrix')
        confusion_matrix = pd.crosstab(df['label_region'], df['predict_region'], rownames=['Actual'],
                                       colnames=['Predicted'])
        # confusion_matrix_object = ConfusionMatrix(df['label_region'], df['predict_region'])
        # confusion_matrix_object.print_stats()
        sns.heatmap(confusion_matrix, annot=False, cmap='Wistia', xticklabels=1, yticklabels=1)
        plt.title('Predicted Region Confusion Matrix', fontsize=20)
        plt.show()

    def check_label(self,read_start, read_end, regions_df, full=True):
        # Assuming no overlaps between regions
        min_overlap = 10
        if full:
            region_index = regions_df.loc[
                (regions_df['chromStart'] <= read_start) & (regions_df['chromEnd'] >= read_end)
                ].index.tolist()
        if not full:
            region_index = regions_df.loc[
                (regions_df['chromStart'] <= read_start) & (regions_df['chromEnd'] >= read_end) |
                ((read_end - regions_df['chromStart']) > min_overlap) & (read_start < regions_df['chromStart']) |
                ((regions_df['chromEnd'] - read_start) > min_overlap) & (read_end > regions_df['chromEnd'])
                ].index.tolist()

        if not region_index:
            return -1
        else:
            return region_index[0]

    def get_chromosomes_from_regions_map(self, regions_df):
        d = {}
        for index, row in regions_df.iterrows():
            d[index] = row.chrom
        return d

    def plot_local_align_attrib(self,df, attrib):
        X = list(range(len(df)))
        plt.scatter(x=X,
                    y=df[attrib],
                    c='DarkBlue',
                    label=attrib,
                    alpha=0.7)

        y_avg = df[attrib].mean()
        y_std = df[attrib].std()
        y_median = df[attrib].median()

        y_avg_vec = [y_avg] * len(X)
        y_avg_std_top_vec = [(y_avg + y_std)] * len(X)
        y_median_vec = [y_median] * len(X)

        plt.plot(X, y_avg_vec, color='red', lw=1, ls='--', label="{} mean".format(attrib), alpha=0.7)
        plt.text(0, y_avg_vec[0], str(y_avg_vec[0]))

        plt.plot(X, y_median_vec, lw=1, ls='dashed', label="{} median".format(attrib), alpha=0.7)
        plt.text(0, y_median_vec[0], str(y_median_vec[0]))

        plt.plot(X, y_avg_std_top_vec, color='red', lw=1, ls='--', label="{} std".format(attrib), alpha=0.7)
        plt.text(0, y_avg_std_top_vec[0], str(y_avg_std_top_vec[0]))

        plt.title("Local align {} plot".format(attrib))
        plt.xlabel("Reads")
        plt.ylabel("{}".format(attrib))
        plt.legend()
        plt.show()
        return

    def plot_local_align_attrib_histogram(self,df, attrib):
        sns.displot(df, x=attrib, log=True)
        plt.title("Local align {} log distribution plot".format(attrib))
        plt.show()
        return

    def plot_local_align_multi_attrib(self,df, attrib_list):

        X = list(range(len(df)))
        for attr in attrib_list:
            plt.scatter(x=X,
                        y=df[attr],
                        label=attr,
                        alpha=0.7)
            y_avg = df[attr].mean()
            y_std = df[attr].std()
            y_median = df[attr].median()

            y_avg_vec = [y_avg] * len(X)
            y_avg_std_top_vec = [(y_avg + y_std)] * len(X)
            y_median_vec = [y_median] * len(X)

            plt.plot(X, y_avg_vec, lw=1, ls='--', label="{} mean".format(attr), alpha=0.7)
            plt.text(0, y_avg_vec[0], str(y_avg_vec[0]))

            plt.plot(X, y_median_vec, lw=1, ls='dashed', label="{} median".format(attr), alpha=0.7)
            plt.text(0, y_median_vec[0], str(y_median_vec[0]))

            plt.plot(X, y_avg_std_top_vec, lw=1, ls='dotted', label="{} std".format(attr), alpha=0.7)
            plt.text(0, y_avg_std_top_vec[0], str(y_avg_std_top_vec[0]))

        plt.title("Local align attributes plot")
        plt.xlabel("Reads")
        plt.ylabel("Attributes counts")
        plt.legend()
        plt.show()
        return

    def plot_local_align_multi_attrib_correlation(self,df, attrib_list):

        plt.scatter(x=df[attrib_list[0]],
                    y=df[attrib_list[1]],
                    alpha=0.7)

        plt.title("Local align attributes plot")
        plt.xlabel("{}".format(attrib_list[0]))
        plt.ylabel("{}".format(attrib_list[1]))
        plt.legend()
        plt.show()
        return