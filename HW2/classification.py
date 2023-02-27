# classification.py
# HW2, Computational Genomics, Spring 2022
# andrewid: tczhang

# WARNING: Do not change the file name; Autograder expects it.

import sys

import numpy as np
from scipy.sparse import csc_matrix, save_npz, load_npz

from matplotlib import pyplot as plt
import pandas as pd
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.svm import SVC

import torch
from torch import nn
from torch.utils.data import Dataset, DataLoader

def get_top_gene_filter(data, n_keep = 2000):
    """Select top n_keep most dispersed genes.

    Args:
        data (n x m matrix): input gene expression data of shape num_cells x num_genes
        n_keep (int): number of genes to be kepted after filtration; default 2000

    Returns:
        filter (array of length n_keep): an array of column indices that can be used as an
            index to keep only certain genes in data. Each element of filter is the column
            index of a highly-dispersed gene in data.
    """
    # Remove columns that only consists of 0's
    col_sums = np.sum(np.abs(data), axis = 0)
    zero_cols = np.nonzero(col_sums == 0)[0]
    data = np.delete(data, zero_cols, axis = 1)

    # Calculated Dispersion
    std = np.std(data, axis = 0)
    mean = np.mean(data, axis = 0)
    dispersion = std/mean

    idx_dispersion = zip(range(len(dispersion)), dispersion)
    idx_dispersion = sorted(idx_dispersion, key = lambda x: x[1], reverse=True)
    top_idx = [0 for x in range(n_keep)]
    for x in range(n_keep):
        top_idx[x] = idx_dispersion[x][0]
    return top_idx

def reduce_dimensionality_pca(filtered_train_gene_expression, filtered_test_gene_expression, n_components = 20):
    """Train a PCA model and use it to reduce the training and testing data.
    
    Args:
        filtered_train_gene_expression (n_train x num_top_genes matrix): input filtered training expression data 
        filtered_test_gene_expression (n_test x num_top_genes matrix): input filtered test expression data 
        
    Return:
        (reduced_train_data, reduced_test_data): a tuple of
            1. The filtered training data transformed to the PC space.
            2. The filtered test data transformed to the PC space.
    """
    combined_data = np.concatenate([filtered_train_gene_expression, filtered_test_gene_expression], axis = 0)
    pca = PCA(n_components)    
    pca = pca.fit(combined_data)
    train_pca = pca.transform(filtered_train_gene_expression)
    test_pca = pca.transform(filtered_test_gene_expression)
    return train_pca, test_pca

def plot_transformed_cells(reduced_train_data, train_labels):
    """Plot the PCA-reduced training data using just the first 2 principal components.
    
    Args:
        reduced_train_data (n_train x num_components matrix): reduced training expression data
        train_labels (array of length n_train): array of cell type labels for training data
        
    Return:
        None

    """
    x = reduced_train_data[:, 0]
    y = reduced_train_data[:, 1]
    plot_df = pd.DataFrame({"x": x, "y": y, "label": train_labels})
    axs = sns.scatterplot(data = plot_df, x="x", y = "y", hue = 'label', s = 10)
    axs.set_xlabel("PCA Feature 1")
    axs.set_ylabel("PCA Feature 2")
    axs.set_title("Training Data PCA")
    fig = axs.get_figure()
    fig.savefig(".\\images\\pca.png")

    
def train_and_evaluate_svm_classifier(reduced_train_data, reduced_test_data, train_labels, test_labels):
    """Train and evaluate a simple SVM-based classification pipeline.
    
    Before passing the data to the SVM module, this function scales the data such that the mean
    is 0 and the variance is 1.
    
    Args:
        reduced_train_data (n_train x num_components matrix): reduced training expression data
        train_labels (array of length n_train): array of cell type labels for training data
        
    Return:
        (classifier, score): a tuple consisting of
            1. classifier: the trained classifier
            2. The score (accuracy) of the classifier on the test data.

    """
    pipeline = make_pipeline(SVC())
    pipeline.fit(reduced_train_data, train_labels)
    test_accuracy = pipeline.score(reduced_test_data, test_labels)
    return pipeline, test_accuracy


        
if __name__ == "__main__":

    train_gene_expression = load_npz(sys.argv[1]).toarray()
    test_gene_expression = load_npz(sys.argv[2]).toarray()
    train_labels = np.load(sys.argv[3])
    test_labels = np.load(sys.argv[4])
    
    # top_gene_filter = get_top_gene_filter(train_gene_expression)
    # filtered_test_gene_expression = test_gene_expression[:, top_gene_filter]
    # filtered_train_gene_expression = train_gene_expression[:, top_gene_filter]

    # pca_train_data, pca_test_data = reduce_dimensionality_pca(filtered_train_gene_expression,
    #                                                            filtered_test_gene_expression)
    # plot_transformed_cells(pca_train_data, train_labels)
    # classifier, accuracy = train_and_evaluate_svm_classifier(pca_train_data, pca_test_data, train_labels, test_labels)
 
    mode = sys.argv[5]
    if mode == "svm_pipeline":
        top_gene_filter = get_top_gene_filter(train_gene_expression)
        filtered_test_gene_expression = test_gene_expression[:, top_gene_filter]
        filtered_train_gene_expression = train_gene_expression[:, top_gene_filter]

        pca_train_data, pca_test_data = reduce_dimensionality_pca(filtered_train_gene_expression,
                                                                filtered_test_gene_expression)

        classifier, test_accuracy = train_and_evaluate_svm_classifier(pca_train_data, pca_test_data,
                                                                  train_labels, test_labels)

        print(f"Training Accuracy: {classifier.score(pca_train_data, train_labels): .10f} \n"+\
              f"Testing Accuracy: {test_accuracy: .10f}")
