import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


class VisualizationsClass:
    def __init__(self):
        return

    def plot_kmer_read_cap_heatmap(self,similarity_mtrx):
        plt.figure(figsize=(18, 8))
        print('Plotting Read to Region Similarity')
        ax = sns.heatmap(similarity_mtrx, cmap='Wistia')
        print(similarity_mtrx.shape)
        plt.xlabel('Region Index', fontsize=18)
        plt.ylabel('Read Index', fontsize=18)
        plt.title('Read to Capture Similarity Confidence', fontsize=20)
        plt.show()