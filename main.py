import json
import os

from src.alignment import msa
from src.generate import generate_datasets
from src.hmm import (
    calculate_emissions_matrix,
    calculate_transitions_matrix,
    create_states_sequence,
)
from src.viterbi import viterbi
from src.retrain import retrain_emissions_matrix, retrain_transitions_matrix

A_COUNT = 20
B_COUNT = 140
C_COUNT = 40


def create_datasets():
    if not os.path.exists("datasets"):
        os.makedirs("datasets")

    datasets = generate_datasets(A_COUNT, B_COUNT, C_COUNT)

    for i, dataset in enumerate(datasets):
        letter = chr(ord('`') + i + 1)
        filename = f"datasets/dataset_{letter}.txt"
        with open(filename, "w") as f:
            output = "\n".join(dataset)
            f.write(output)


def main():
    expected_files = [
        "datasets/dataset_a.txt",
        "datasets/dataset_b.txt",
        "datasets/dataset_c.txt",
    ]
    if not all(os.path.exists(f) for f in expected_files):
        create_datasets()

    dataset_a = []

    file_a = open("datasets/dataset_a.txt", "r")
    for line in file_a:
        dataset_a.append(line.strip())
    file_a.close()
    
    dataset_b = []

    file_b = open("datasets/dataset_b.txt", "r")
    for line in file_b:
        dataset_b.append(line.strip())
    file_b.close()

    msa_a = msa(dataset_a)

    # 1. Build the Initial Draft HMM (from Dataset A)
    states_seq_a = create_states_sequence(msa_a)
    num_matches = states_seq_a.count("M") # Save this architectural constant!
    
    initial_emissions = calculate_emissions_matrix(msa_a, states_seq_a)
    initial_transitions = calculate_transitions_matrix(msa_a, states_seq_a)

    # 2. Viterbi Training (Guesstimating the paths for Dataset B)
    dataset_b_paths = []
    dataset_b_clean = []
    
    for seq in dataset_b:
        # Crucial: Remove all gaps from Dataset B before running Viterbi
        clean_seq = seq.replace("-", "") 
        dataset_b_clean.append(clean_seq)
        
        path = viterbi(clean_seq, initial_emissions, initial_transitions)
        dataset_b_paths.append(path)

    # 3. Generating the Trained HMM!
    trained_emissions = retrain_emissions_matrix(dataset_b_clean, dataset_b_paths, num_matches)
    trained_transitions = retrain_transitions_matrix(dataset_b_paths, num_matches)

    print("HMM Successfully Trained!")
    print(json.dumps(trained_emissions, indent=2))
    print(json.dumps(trained_transitions, indent=2))

if __name__ == "__main__":
    main()
