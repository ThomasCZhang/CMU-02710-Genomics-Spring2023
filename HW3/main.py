import os
import numpy as np

def main():
    dirname = os.path.dirname(__file__)

    filepath = os.path.join(dirname, "provided_data", "werepyre-test-annotation.fasta")
    true_state_path = ReadFasta(filepath)[0]
    filepath = os.path.join(dirname, "provided_data", "werepyre-test.fasta")
    test_seq = ReadFasta(filepath)[0]
    filepath = os.path.join(dirname, "provided_data", "werepyre-train.fasta")
    train_seq = ReadFasta(filepath)

    # Viterbi Training on Training Data
    emission_vals, state_vals, transition_mm, emission_mm = GetInitialHMM()
    _, transition_mm, emission_mm = ViterbiTraining(train_seq, emission_vals, state_vals, transition_mm, emission_mm)
    state_path = ViterbiDecoding(test_seq, emission_vals, state_vals, transition_mm, emission_mm)
    answerpath = os.path.join(dirname, "answer1.txt")
    WriteAnswer(answerpath, state_path, true_state_path, test_seq)

    # Viterbi decoding on just the test data with true parameters.
    emission_vals, state_vals, transition_mm, emission_mm = GetTrueHMM()
    state_path = ViterbiDecoding(test_seq, emission_vals, state_vals, transition_mm, emission_mm)
    answerpath = os.path.join(dirname, "answer2.txt")
    WriteAnswer(answerpath, state_path, true_state_path, test_seq)

def WriteAnswer(answerpath, state_path, true_state_path, test_seq):
    num_diff = CountLetterDifferences(true_state_path, state_path)    
    percent_diff = 1-num_diff/len(test_seq)
    with open(answerpath, 'w') as f:
        for i, char in enumerate(state_path):
            if i%60==0 and i != 0:
                f.write("\\\\\n")
            f.write(char)
        f.write("\n")
        f.write(f"1 - {num_diff}/{len(test_seq)} = {percent_diff: 0.6f} = {100*percent_diff: 6.4f}%")

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

def ViterbiTraining(sequences: list[str] , emission_vals: list[str], state_vals: list[str], transition_mm: dict[str, dict[str, float]],
    emission_mm: dict[str, dict[str, float]]) -> str:
    """
    Performs viterbi training determine state sequences.
    Input:
        seq: the input sequences.
        emission_vals: the list of possible emitted characters.
        state_vals: the list of possible state characters.
        transition_mm: A dictionary holding all the transition probabilities.
        emission_mm: A dictionary holding all the emission probabilities.
    Output:
        List of predicted state sequences for each input sequences.
        The transition probabilities.
        The emission probabilities.
    """
    current_path = ["" for _ in sequences]
    new_paths = [None for _ in sequences]
    while current_path != new_paths:
        current_path = new_paths
        for i, seq in enumerate(sequences):
            state_path = ViterbiDecoding(seq, emission_vals, state_vals, transition_mm, emission_mm)
            new_paths[i]  = state_path
        emission_mm = UpdateEmission(sequences, new_paths, state_vals, emission_vals)
        transition_mm = UpdateTransition(new_paths, state_vals)

    return current_path, transition_mm, emission_mm

def UpdateEmission(sequences: list[str], state_paths: list[str], state_vals: list[str], emission_vals: list[str]):
    """
    Creates an emission matrix based on a sequence of emitted characters and a sequence of states.
    Input:
        seq: the sequence of emitted characters
        state_path: the sequence of states
        state_vals: the possible states.
    Output:
        The updated emission probabiltiies.
    """
    emission_mm = {}
    for s in state_vals:
        emission_mm[s] = {}
        for e in emission_vals:
            emission_mm[s][e] = 0
    
    for seq, path in zip(sequences, state_paths):
        for char, state in zip(seq, path):
            emission_mm[state][char] += 1

    for s in state_vals:
        total = 0
        for e in emission_vals:
            total += emission_mm[s][e] 
        for e in emission_vals:
            emission_mm[s][e] /= total
    return emission_mm

def UpdateTransition(state_paths: list[str], state_vals: list[str]):
    """
    Creates an transition matrix based on a sequence of emitted characters and a sequence of states.
    Input:
        seq: the sequence of emitted characters
        state_path: the sequence of states
        state_vals: the possible states.
    Output:
        the updated transition probabilties.
    """
    transition_mm = {}
    for a in state_vals:
        transition_mm[a] = {}
        for b in state_vals:
            transition_mm[a][b] = 0
    
    for path in state_paths:
        for i in range(len(path)-1):
            curr_state = path[i]
            next_state = path[i+1]
            transition_mm[curr_state][next_state] += 1

    for a in state_vals:
        total = 0
        for b in state_vals:
            total += transition_mm[a][b]
        for b in state_vals:
            transition_mm[a][b] /= total
    return transition_mm


def ViterbiDecoding(
    seq: str , emission_vals: list[str], state_vals: list[str], transition_mm: dict[str, dict[str, float]],
    emission_mm: dict[str, dict[str, float]]) -> str:
    """
    Uses the Viterbi algorithm to determine the sequence of states for an input sequence.
    Input:
        seq: The input sequence.
        state_0: The initial probabilities of each state.
        transitions: transition probabilities
        emissions: emission probabilities
    Output:
        Sequence of states.
    """
        
    score_matrix = np.zeros([len(state_vals), len(seq)])
    path_matrix = np.zeros(score_matrix.shape) - 1  # Initialize a matrix of -1. 0 = vampire, 1 = werewolf
    for row, state in enumerate(state_vals):
        score_matrix[row, 0] = np.log(0.5) + emission_mm[state][seq[0]]

    for col in range(1, score_matrix.shape[1]):
        for row in range(score_matrix.shape[0]):
            new_state = state_vals[row]
            best_score = -np.inf
            best_prev_state = -1
            for i, old_state in enumerate(state_vals):  # Number of states
                score = score_matrix[i, col - 1] + transition_mm[old_state][new_state] + emission_mm[new_state][seq[col]]
                if score > best_score:
                    best_score = score
                    best_prev_state = i

            score_matrix[row, col] = best_score
            path_matrix[row, col] = best_prev_state

    start_idx = int(score_matrix[:, -1].argmax())
    path = Backtrack(path_matrix, state_vals, start_idx)
    return path

def Backtrack(path_matrix: np.ndarray, state_vals: list[str], start_idx: int) -> str:
    """
    Generates the path of states given the path matrix.
    Input:
        path_matrix: matrix containing path data.
        start_idx: the starting row in the last column of path_matrix.
    Output:
        the path as a string.
    """
    current_idx = start_idx
    path = ""
    for i in range(path_matrix.shape[1] - 1, -1, -1):
        path = "".join([ state_vals[current_idx], path])
        current_idx = int(path_matrix[current_idx, i])
    return path

def CountLetterDifferences(s1: str, s2: str) -> int:
    """
    Count the number of different letters between two strings element wise.
    Input:
        s1, s2: the input strings
    Output:
        The number of differing characters.
    """
    if len(s1) != len(s2):
        raise Exception("Length of sequences are not the same.")

    count = 0
    for v, e in zip(s1, s2):
        if v != e:
            count += 1
    return count

if __name__ == "__main__":
    main()
