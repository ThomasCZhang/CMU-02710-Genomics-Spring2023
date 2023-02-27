# normalization.py
# HW2, Computational Genomics, Spring 2022
# andrewid: tczhang

# WARNING: Do not change the file name; Autograder expects it.

import sys
import numpy as np
import matplotlib.pyplot as plt

PER_MILLION = 1/1000000
PER_KILOBASE = 1/1000

# Do not change this function signature
def rpkm(raw_counts, gene_lengths):
    """Find the normalized counts for raw_counts
    Returns: a matrix of same size as raw_counts
    """
    r = raw_counts
    l = gene_lengths[:, np.newaxis]
    m = np.sum(raw_counts, axis = 0)[np.newaxis,:]
    var_rpkm = r/(m*PER_MILLION)/(l*PER_KILOBASE)
    return var_rpkm
    
# Do not change this function signature
def tpm(raw_counts, gene_lengths):
    """Find the normalized counts for raw_counts
    Returns: a matrix of same size as raw_counts
    """
    r = raw_counts
    l = gene_lengths[:, np.newaxis]
    rpk = r/(l*PER_KILOBASE)
    var_tpm = rpk/(np.sum(rpk, axis = 0)*PER_MILLION)
    return var_tpm
   
# define any helper function here    


# Do not change this function signature
def size_factor(raw_counts) -> np.ndarray:
    """Find the normalized counts for raw_counts
    Returns: a matrix of same size as raw_counts
    """
    geo_means = np.exp(np.mean(np.log(raw_counts), axis = 1))[:, np.newaxis]
    s = np.median(raw_counts/geo_means, axis = 0)
    sf = raw_counts/s
    return sf

def MakeBoxPlot(data: np.ndarray, PlotTitle: str, ylabel:str, savePath: str):
    fig = plt.figure()
    axs = fig.subplots()
    bp = axs.boxplot(data, flierprops = {'marker': '.', 'markersize': 5,
                                          'markerfacecolor': (0.2, 0.6, 1),
                                            'markeredgewidth': 0})
    fig.suptitle(PlotTitle)
    axs.set_ylabel(ylabel)
    axs.set_xlabel("Sample")
    fig.savefig(savePath)
    plt.close(fig)

if __name__=="__main__":
    raw_counts=np.loadtxt(sys.argv[1])
    gene_lengths=np.loadtxt(sys.argv[2])
    
    rpkm1=rpkm(raw_counts, gene_lengths)
    tpm1=tpm(raw_counts, gene_lengths)
    size_factor1=size_factor(raw_counts)

    save_path = "size_count_normalized.txt"
    np.savetxt(save_path, size_factor1)

    # TODO: write plotting code here
    plot_title = "Raw Data"
    y_label = "$log_2$(Raw Data Counts)"
    save_path = ".\\images\\rawdata.png"
    MakeBoxPlot(np.log2(raw_counts), plot_title, y_label, save_path)

    plot_title = "RPKM"
    y_label = "$log_2$(RPKM Normalized Counts)"
    save_path = ".\\images\\rpkm.png"
    MakeBoxPlot(np.log2(rpkm1), plot_title, y_label, save_path)

    plot_title = "TPM"
    y_label = "$log_2$(TPM Normalized Counts)"
    save_path = ".\\images\\tpm.png"
    MakeBoxPlot(np.log2(tpm1), plot_title, y_label, save_path)

    plot_title = "Size Factor"
    y_label = "$log_2$(Size Factor Normalized Counts)"
    save_path = ".\\images\\size_factor.png"
    MakeBoxPlot(np.log2(size_factor1), plot_title, y_label, save_path)
    

