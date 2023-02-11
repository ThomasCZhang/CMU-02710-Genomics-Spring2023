# NW.py
# HW1, Computational Genomics, Spring 2022
# andrewid:

# WARNING: Do not change the file name; Autograder expects it.

import sys

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

# Do not change this function signature
def needleman_wunsch(seq1, seq2):
    """Find the global alignment for seq1 and seq2
    Returns: 3 items as so:
    the alignment score, alignment in seq1 (str), alignment in seq2 (str)
    """
    raise NotImplementedError

if __name__=="__main__":
    # path = ".\\HW1\\test.txt"
    # Sequences = ReadFASTA(path)

    Sequences=ReadFASTA(sys.argv[1])
    assert len(Sequences.keys())==2, "fasta file contains more than 2 sequences."
    seq1=Sequences[list(Sequences.keys())[0]]
    seq2=Sequences[list(Sequences.keys())[1]]

    
    score, align1, align2 = needleman_wunsch(seq1, seq2)

    print('Score: ', score)
    print('Seq1: ', align1)
    print('Seq2: ', align2)

def needleman_wunsch(seq1: str, seq2: str) -> tuple[int, str, str]:
    pass