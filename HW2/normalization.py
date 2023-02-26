# normalization.py
# HW2, Computational Genomics, Spring 2022
# andrewid: 

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
    
    pass
    
# Do not change this function signature
def tpm(raw_counts, gene_lengths):
    """Find the normalized counts for raw_counts
    Returns: a matrix of same size as raw_counts
    """

    pass
   
# define any helper function here    


# Do not change this function signature


def size_factor(raw_counts):
    """Find the normalized counts for raw_counts
    Returns: a matrix of same size as raw_counts
    """
    pass
    

if __name__=="__main__":
    raw_counts=np.loadtxt(sys.argv[1])
    gene_lengths=np.loadtxt(sys.argv[2])
    
    rpkm1=rpkm(raw_counts, gene_lengths)
    tpm1=tpm(raw_counts, gene_lengths)
    size_factor1=size_factor(raw_counts)

    # TODO: write plotting code here
    pass
