import json
import os

from src.alignment import msa
from src.generate import generate_datasets
from src.hmm import (
    calculate_emissions_matrix,
    calculate_transitions_matrix,
    create_states_sequence,
)

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

    msa_a = msa(dataset_a)
    for seq in msa_a:
        print(seq)

    states_seq = create_states_sequence(msa_a)
    print(''.join(states_seq), '\n')

    emissions = calculate_emissions_matrix(msa_a, states_seq)
    print(json.dumps(emissions, indent=2))

    transitions = calculate_transitions_matrix(msa_a, states_seq)
    print(json.dumps(transitions, indent=2))


if __name__ == "__main__":
    main()
