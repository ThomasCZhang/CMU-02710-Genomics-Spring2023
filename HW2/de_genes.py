# de_genes.py
# HW2, Computational Genomics, Spring 2022
# andrewid:

# WARNING: Do not change the file name; Autograder expects it.

import sys
import numpy as np


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

    pass

# define any helper function here    

if __name__=="__main__":
    # Here is a free test case
    genes=['a', 'b', 'c']
    input1 = [0.01, 0.04, 0.1]
    print(bh(genes, input1, 0.05))
