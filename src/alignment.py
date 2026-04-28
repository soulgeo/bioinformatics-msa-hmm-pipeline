from itertools import product

import numpy as np

ALPHA = 1  # ID: p19153 ends in an odd number

SAME = 1
DIFF = -ALPHA / 2
DELTA = -ALPHA


def sim(char1, char2):
    return SAME if char1 == char2 else DIFF


def align_pair(seq1, seq2):
    # Calculate the similarity matrix
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    sim_matrix = np.empty((rows, cols))

    for i in range(rows):
        sim_matrix[i, 0] = DELTA * i
    for j in range(cols):
        sim_matrix[0, j] = DELTA * j

    for i, j in product(range(1, rows), range(1, cols)):
        match = sim_matrix[i - 1, j - 1] + sim(seq1[i - 1], seq2[j - 1])
        delete = sim_matrix[i - 1, j] + DELTA
        insert = sim_matrix[i, j - 1] + DELTA
        sim_matrix[i, j] = max(match, delete, insert)

    # Using the similarity matrix, find an optimal alignment
    alignment1 = []
    alignment2 = []

    i, j = len(seq1), len(seq2)

    while i > 0 or j > 0:
        if (
            i > 0
            and j > 0
            and sim_matrix[i, j]
            == sim_matrix[i - 1, j - 1] + sim(seq1[i - 1], seq2[j - 1])
        ):
            alignment1.append(seq1[i - 1])
            alignment2.append(seq2[j - 1])
            i -= 1
            j -= 1

        elif i > 0 and sim_matrix[i, j] == sim_matrix[i - 1, j] + DELTA:
            alignment1.append(seq1[i - 1])
            alignment2.append("-")
            i -= 1

        else:
            alignment1.append("-")
            alignment2.append(seq2[j - 1])
            j -= 1

    output1 = ''.join(reversed(alignment1))
    output2 = ''.join(reversed(alignment2))
    return output1, output2


def align_multiple_sequences(sequences):
    if len(sequences) < 2:
        return sequences

    if len(sequences) == 2:
        return align_pair(sequences[0], sequences[1])

    c, other_sequences = sequences[0], sequences[1:]

    aligned = [c]
    for seq in other_sequences:
        aligned_c, aligned_seq = align_pair(c, seq)

        new_gap_indices = [i for i, char in enumerate(aligned_c) if char == "-"]

        current_gap_indices = [idx for idx, char in enumerate(aligned[0]) if char == "-"]

        # Step A: Insert old gaps in new alignment
        temp_new_seq = list(aligned_seq)
        for gap in sorted(current_gap_indices, reverse=True):
            temp_new_seq.insert(gap, "-")

        # Step B: Insert new gaps in old alignments
        for i in range(len(aligned)):
            temp_seq = list(aligned[i])
            for gap in sorted(new_gap_indices, reverse=True):
                temp_seq.insert(gap, "-")
            aligned[i] = "".join(temp_seq)

        aligned.append("".join(temp_new_seq))

    return aligned
