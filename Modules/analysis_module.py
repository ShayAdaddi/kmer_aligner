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

    def check_label(self,read_start, read_end, regions_df, full=False):
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