import os

from src.training.alignment import msa
from src.generation.generate import generate_datasets, generate_random_gene
from src.training.hmm import (
    calculate_emissions_matrix,
    calculate_transitions_matrix,
    create_states_sequence,
)
from src.training.viterbi import forward_algorithm, viterbi 
from src.training.retrain import retrain_emissions_matrix, retrain_transitions_matrix

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

    datasets = {}
    for char in ['a', 'b', 'c']:
        datasets[char] = []
        file = open(f"datasets/dataset_{char}.txt", "r")
        for line in file:
            datasets[char].append(line.strip())
        file.close()

    msa_a = msa(datasets['a'])

    # 1. Build the Initial Draft HMM (from Dataset A)
    states_seq_a = create_states_sequence(msa_a)
    num_matches = states_seq_a.count("M")
    
    initial_emissions = calculate_emissions_matrix(msa_a, states_seq_a)
    initial_transitions = calculate_transitions_matrix(msa_a, states_seq_a)

    # 2. Viterbi Training
    dataset_b_paths = []
    
    for gene in datasets['b']:
        path = viterbi(gene, initial_emissions, initial_transitions)
        dataset_b_paths.append(path)

    # 3. Generating the Trained HMM
    trained_emissions = retrain_emissions_matrix(datasets['b'], dataset_b_paths, num_matches)
    trained_transitions = retrain_transitions_matrix(dataset_b_paths, num_matches)

    dataset_c_scores = []
    for gene in datasets['c']:
        score = forward_algorithm(gene, trained_emissions, trained_transitions)
        dataset_c_scores.append(score)

    print(dataset_c_scores)

    random_genes = [generate_random_gene() for _ in range(40)]
    random_genes_scores = []
    for gene in random_genes:
        score = forward_algorithm(gene, trained_emissions, trained_transitions)
        random_genes_scores.append(score)

    print(random_genes_scores)


if __name__ == "__main__":
    main()
