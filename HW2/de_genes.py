# de_genes.py
# HW2, Computational Genomics, Spring 2022
# andrewid: tczhang

# WARNING: Do not change the file name; Autograder expects it.

import sys
import numpy as np
import matplotlib.pyplot as plt


# Do not change this function signature

def bh(genes, pvals, alpha):
    """(list, list, float) -> numpy array
    applies benjamini-hochberg procedure
    
    Parameters
    ----------
    genes: name of genes 
    pvalues: corresponding pvals
    alpha: desired false discovery rate
    
    Returns
    -------
    array containing gene names of significant genes.
    gene names do not need to be in any specific order.
    """
    n = len(genes)
    ordered_pval = sorted(zip(genes, pvals), key = lambda t: t[1])
    rank = 0
    while rank < n and ordered_pval[rank][1] < alpha*(rank+1)/n:
        rank += 1
    significant_genes = np.asarray([x[0] for x in ordered_pval[0:rank]])
    return significant_genes

# define any helper function here       
    
if __name__=="__main__":
    # Here is a free test case
    genes=['a', 'b', 'c']
    input1 = [0.01, 0.04, 0.1]
    print(bh(genes, input1, 0.05))
