# NW.py
# HW1, Computational Genomics, Spring 2022
# andrewid: tczhang

# WARNING: Do not change the file name; Autograder expects it.

import sys
import numpy as np

def ReadFASTA(filename):
    fp=open(filename, 'r')
    Sequences={}
    tmpname=""
    tmpseq=""
    for line in fp:
        if line[0]==">":
            if len(tmpseq)!=0:
                Sequences[tmpname]=tmpseq
            tmpname=line.strip().split()[0][1:]
            tmpseq=""
        else:
            tmpseq+=line.strip()
    Sequences[tmpname]=tmpseq
    fp.close()
    return Sequences

# You may define any helper functions for Needleman-Wunsch algorithm here
def Backtrack(seq1: str, seq2: str, path_matrix: np.ndarray) -> tuple[str, str]:
    """
    Returns the "path" taken by Needleman Wunsch algorithm for aligning two strings.

    Input:
        seq1, seq2: The two strings being aligned

        path_matrix: A 3 dimensional numpy ndarray that contains info on the path taken by Needleman Wunsh
            The third dimension holds info on whether a letter was taken from seq1, seq2 or from both.
    
    Output:
        Two strings corresponding to the input sequences. "-" represents a gap.
    """
    final_seq1 = ""
    final_seq2 = ""
    
    row = len(seq1)
    col = len(seq2)
    
    while (row != 0) or (col != 0):
        if path_matrix[row, col, 2] == 1:
            row -= 1
            col -= 1
            final_seq1 = seq1[row] + final_seq1
            final_seq2 = seq2[col] + final_seq2
        elif path_matrix[row, col, 1] == 1:
            col -= 1
            final_seq1 = "-" + final_seq1
            final_seq2 = seq2[col] + final_seq2
        elif path_matrix[row, col, 0] == 1:
            row -= 1
            final_seq1 = seq1[row] + final_seq1
            final_seq2 = "-" + final_seq2
    
    return final_seq1, final_seq2

# Do not change this function signature
def needleman_wunsch(seq1: str, seq2: str) -> tuple[int, str, str]:
    """Find the global alignment for seq1 and seq2
    Returns: 3 items as so:
    the alignment score, alignment in seq1 (str), alignment in seq2 (str)
    """
    score_matrix = np.zeros((len(seq1)+1, len(seq2)+1))  # Matrix that stores the scores of the subproblems that have been solved
    path_matrix = np.zeros((len(seq1)+1, len(seq2)+1, 3))  # Matrix that scores the "path" has been taken so far.
        # [i,j,0] take from str1,
        # [i,j,1] take from str2,
        # [i,j,2] take from both strings.

    # Defining the scoring rules
    char_match = 1
    char_mismatch = -2
    gap = -1

    for i in range(1, score_matrix.shape[0]):
        score_matrix[i,0] = score_matrix[i-1, 0] + gap
        path_matrix[i, 0, 0] = 1
    
    for i in range(1, score_matrix.shape[1]):
        score_matrix[0, i] = score_matrix[0, i-1] + gap
        path_matrix[0, i, 1] = 1

    for i in range(1, score_matrix.shape[0]):
        for j in range(1, score_matrix.shape[1]):
            match_score = char_mismatch
            if seq1[i-1] == seq2[j-1]:
                match_score = char_match

            seq1_score = score_matrix[i-1,j] + gap              # Score for taking letter from seq1.
            seq2_score = score_matrix[i, j-1] + gap             # Score for taking letter from seq2.
            both_score = score_matrix[i-1, j-1] + match_score   # Score for taking letters from both seq1 and seq2.

            score_matrix[i,j] = max(seq1_score, seq2_score, both_score)
            if score_matrix[i,j] == seq1_score:
                path_matrix[i,j,0] = 1
            if score_matrix[i,j] == seq2_score:
                path_matrix[i,j,1] = 1
            if score_matrix[i,j] == both_score:
                path_matrix[i,j,2] = 1

    final_seq1, final_seq2 = Backtrack(seq1, seq2, path_matrix)
    return (score_matrix[len(seq1), len(seq2)], final_seq1, final_seq2)

if __name__=="__main__":
    Sequences=ReadFASTA(sys.argv[1])
    assert len(Sequences.keys())==2, "fasta file contains more than 2 sequences."
    seq1=Sequences[list(Sequences.keys())[0]]
    seq2=Sequences[list(Sequences.keys())[1]]

    
    score, align1, align2 = needleman_wunsch(seq1, seq2)

    print('Score: ', score)
    print('Seq1: ', align1)
    print('Seq2: ', align2)