# classification.py
# HW2, Computational Genomics, Spring 2022
# andrewid: 

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

    pass

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
    
    pass

def plot_transformed_cells(reduced_train_data, train_labels):
    """Plot the PCA-reduced training data using just the first 2 principal components.
    
    Args:
        reduced_train_data (n_train x num_components matrix): reduced training expression data
        train_labels (array of length n_train): array of cell type labels for training data
        
    Return:
        None

    """
    
    pass
    
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
    
    pass
        
if __name__ == "__main__":
    train_gene_expression = load_npz(sys.argv[1]).toarray()
    test_gene_expression = load_npz(sys.argv[2]).toarray()
    train_labels = np.load(sys.argv[3])
    test_labels = np.load(sys.argv[4])
    
    top_gene_filter = get_top_gene_filter(train_gene_expression)
    filtered_test_gene_expression = test_gene_expression[:, top_gene_filter]
    filtered_train_gene_expression = train_gene_expression[:, top_gene_filter]
        
    mode = sys.argv[5]
    if mode == "svm_pipeline":
        # TODO: Implement the pipeline here
        pass
