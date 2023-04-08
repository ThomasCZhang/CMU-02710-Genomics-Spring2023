import os, math
from glob import glob
import numpy as np

def main():
    dirname = os.path.dirname(__file__)
    filepath = os.path.join(dirname, "provided_data", "werepyre-test-annotation.fasta")
    true_state_path = ReadFasta(filepath)[0]
    
    filepath = os.path.join(dirname, "provided_data", "werepyre-test.fasta")
    sequences = ReadFasta(filepath)[0]
    emission_vals, state_vals, transition_mm, emission_mm = GetInitialHMM()
    state_path = Viterbi(sequences, emission_vals, state_vals, transition_mm, emission_mm)

    num_diff = CountLetterDifferences(true_state_path, state_path)    
    percent_diff = 1-num_diff/len(sequences)
    answerpath = os.path.join(dirname, "answer1.txt")
    with open(answerpath, 'w') as f:
        for i, char in enumerate(state_path):
            if i%60==0 and i != 0:
                f.write("\\\\\n")
            f.write(char)
        f.write("\n")
        f.write(f"1 - {num_diff}/{len(sequences)} = {percent_diff}")

    emission_vals, state_vals, transition_mm, emission_mm = GetTrueHMM()
    state_path = Viterbi(sequences, emission_vals, state_vals, transition_mm, emission_mm)

    num_diff = CountLetterDifferences(true_state_path, state_path)    
    percent_diff = 1-num_diff/len(sequences)
    answerpath = os.path.join(dirname, "answer2.txt")
    with open(answerpath, 'w') as f:
        for i, char in enumerate(state_path):
            if i%60==0 and i != 0:
                f.write("\\\\\n")
            f.write(char)
        f.write("\n")
        f.write(f"1 - {num_diff}/{len(sequences)} = {percent_diff}")


    

def ReadFasta(filename):
    with open(filename) as f:
        sequences = []
        for line in f:
            if line[0] == ">":
                sequence = ""
                sequences.append(sequence)
            else:
                sequences[-1] = "".join([sequences[-1], line.strip()])

    return sequences

def GetInitialHMM():
    emission_vals = ["A", "C", "T", "G"]
    state_vals = ["V", "W"]
    transition_mm = {"V": {"V": 0.75, "W": 0.25}, "W": {"V": 0.25, "W": 0.75}}
    emission_mm = {"V": {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}, "W": {"A": 0.2, "C": 0.3, "G": 0.3, "T": 0.2}}
    for key in transition_mm:
        for k in transition_mm[key]:
            transition_mm[key][k] = np.log(transition_mm[key][k])
        for k in emission_mm[key]:
            emission_mm[key][k] = np.log(emission_mm[key][k])

    return emission_vals, state_vals, transition_mm, emission_mm

def GetTrueHMM():
    emission_vals = ["A", "C", "T", "G"]
    state_vals = ["V", "W"]
    transition_mm = {"V": {"V": 0.99, "W": 0.010}, "W": {"V": 0.014, "W": 0.986}}
    emission_mm = {"V": {"A": 0.35, "C": 0.15, "G": 0.15, "T": 0.35}, "W": {"A": 0.175, "C": 0.325, "G": 0.325, "T": 0.175}}
    for key in transition_mm:
        for k in transition_mm[key]:
            transition_mm[key][k] = np.log(transition_mm[key][k])
        for k in emission_mm[key]:
            emission_mm[key][k] = np.log(emission_mm[key][k])

    return emission_vals, state_vals, transition_mm, emission_mm

def Viterbi(output:str, emission_vals: list[str], state_vals: list[str], transition_mm: dict[str, dict[str, float]],
            emission_mm: dict[str, dict[str, float]]) -> str:
    """
    Viterbi: Finds the most likely sequence of states that create a specific outcome based on a state transition
    matrix and an emission matrix.

    Input:
        output: The given output.
        emission_vals: The possible emission values.
        state_vals: The possible states.
        transition_mm: The markov model for state transitions.
        emission_mm: The model for emissions based on state.
    
    Output:
        The most likely sequence of states.
    """

    score_matrix = [[0 for y in output] for x in state_vals]
    for x in range(len(state_vals)):
        score_matrix[x][0] = emission_mm[state_vals[x]][output[0]]+np.log(0.5) # Assume equal starting probabilities

    prior_dict = {}
    for x in state_vals:
        prior_dict[x] = False
    path_matrix = [[prior_dict.copy() for y in output] for x in state_vals]

    for x in range(1, len(output)):
        for y in range(len(state_vals)):
            current_state = state_vals[y]
            scores = {}
            for z in range(len(state_vals)):
                prior_state = state_vals[z]
                scores[prior_state] = score_matrix[z][x-1]+\
                    transition_mm[prior_state][current_state]+\
                    emission_mm[current_state][output[x]]
            score_matrix[y][x] = max(scores.values())
            
            for key in scores:
                if score_matrix[y][x] == scores[key]:
                    path_matrix[y][x][key] = True
    
    final_state = 0
    final_state_score = score_matrix[final_state][-1]
    for x in range(1, len(state_vals)):
        if score_matrix[x][-1] > final_state_score:
            final_state_score = score_matrix[x][-1]
            final_state = x

    state_path = Backtrack(final_state, state_vals, path_matrix)
    return state_path

def Backtrack(final_state: int, state_vals: list[str], path_matrix: list[list[dict[str, bool]]]):
    """
    Backtrack: Backtrack through a viterbi graph.
    Input:
        final_state: The state to start backtracking from.
        state_vals: a list of the possible states.
        path_matrix: A matrix that holds the path used to solve a viterbi graph.
    Output:
        The path taken as a string.
    """
    state_path = state_vals[final_state]
    current_state = final_state
    i = len(path_matrix[0])-1
    while i > 0:
        for idx, key in enumerate(state_vals):
            if path_matrix[current_state][i][key] == True:
                state_path = "".join([key, state_path])
                current_state = idx
                break
        i -= 1
    return state_path

def CountLetterDifferences(s1: str, s2: str) -> int:
    if len(s1) != len(s2):
        raise Exception("Length of sequences are not the same.")

    count = 0
    for v, e in zip(s1, s2):
        if v != e:
            count += 1
    return count

if __name__ == "__main__":
    main()