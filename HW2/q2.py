import sys
import numpy as np
import matplotlib.pyplot as plt

def Dispersion(data: np.ndarray):
    """
    Calculates the dispersion of gene data.

    Input:
        data: The gene data. Columns correspond to different samples. Each row is a different gene.

    Output:
        The dispersion of each gene accross all the samples. 
    """
    std = np.std(data, axis = 1)
    mean = np.mean(data, axis = 1)
    return std/mean    


def Q2A(data: np.ndarray, labels: np.ndarray):
    """
    Calculates dispersion of data, then plots the log2-log2 dispersion vs the mean of genes
    accross different samples.

    Input:
        data: The gene data.
    """
    dexamethasone_data = data[:, np.asarray(labels == 1).nonzero()[0]]
    ethanol_data = data[:, np.asarray(labels == 2).nonzero()[0]]

    
    save_path = ".\\images\\mean_dispersion.png"

    fig = plt.figure()
    axs = fig.subplots()
    axs.set_xlabel("$Log_2$(Gene Mean)")
    axs.set_ylabel("$Log_2$(Gene Dispersion)")
    axs.set_title("Mean vs Dispersion")
    
    x = np.log2(Dispersion(dexamethasone_data))
    y = np.log2(np.mean(dexamethasone_data, axis = 1))
    line1 = axs.scatter(x, y, s = 5, color = 'royalblue')

    x = np.log2(Dispersion(ethanol_data))
    y = np.log2(np.mean(ethanol_data, axis = 1))
    line2 = axs.scatter(x, y, s = 5, color = 'forestgreen')

    axs.legend((line1, line2),("Dexamethasone", "Ethanol"),loc = "lower left")

    fig.savefig(save_path)
    plt.close(fig)
    # save_path = ".\\images\\ethanol_dispersion.png"

    # fig = plt.figure()
    # axs = fig.subplots()
    # axs.set_xlabel("$Log_2$(Gene Mean)")
    # axs.set_ylabel("$Log_2$(Gene Dispersion)")
    # axs.set_title("Ethanol Treated")
    # axs.scatter(x, y, s = 5)
    # fig.savefig(save_path)
    # plt.close(fig)

def Q2B(data: np.ndarray, labels: np.ndarray, genes: np.ndarray):
    """
    Computes the log2 fold change between the two categories. 
    Then plots the log2-log2 fold change vs mean count.

    Input:
        data: normalized data by size factor
        labels: the treatment label
        genes: gene names
    """
    dexamethasone_data = data[:, np.asarray(labels == 1).nonzero()[0]]
    ethanol_data = data[:, np.asarray(labels == 2).nonzero()[0]]

    fold_change = np.mean(dexamethasone_data, axis = 1)/np.mean(ethanol_data, axis = 1)
    global_mean = np.mean(data, axis = 1)

    x = np.log2(global_mean)
    y = np.log2(fold_change)
    save_path = ".\\images\\Q2b.png"
    fig = plt.figure()
    axs = fig.subplots()
    axs.scatter(x, y, s = 5)
    axs.set_title("Fold Change vs Mean Expression") 
    axs.set_ylabel("$log_2$(Fold Change)")
    axs.set_xlabel("$log2$(Mean)")
    fig.savefig(save_path)
    plt.close(fig)

    sorted_labeled_fc = sorted(zip(fold_change, gene_names),
                                key = lambda t: t[0], reverse = True)
    
    with open(".\\images\\q2b.txt", "w") as f:
        s1, s2 = "Gene Name", "log2 Fold Change"
        f.write(f"{s1:<20} {s2:<20}") 
        f.write("\n----------------------------------------------------------")   
        for i in range(10):
            f.write(f"\n{sorted_labeled_fc[i][1]: <20} {sorted_labeled_fc[i][0]: 6.6f}")
        f.write("\n----------------------------------------------------------")   
        for i in range(10):
            f.write(f"\n{sorted_labeled_fc[-i][1]: <20} {sorted_labeled_fc[-i][0]: 6.6f}")


if __name__ == "__main__":
    labels_path = "C:\\users\\tzhan\\pythonws\\02710_spring2023\\hw2\\labels.txt"
    # labels_path = "labels.txt"
    condition_labels = np.loadtxt(labels_path)
    counts_path = "C:\\users\\tzhan\\pythonws\\02710_spring2023\\hw2\\size_count_normalized_counts.txt"
    # counts_path = "size_count_normalized_counts.txt"
    counts_data = np.loadtxt(counts_path)
    genes_path = "C:\\users\\tzhan\\pythonws\\02710_spring2023\\hw2\\GeneNames.txt"
    # genes_path = "gene_names.txt"
    gene_names = np.loadtxt(genes_path, dtype = str)
    
    Q2A(counts_data, condition_labels)  
    Q2B(counts_data, condition_labels, gene_names) 