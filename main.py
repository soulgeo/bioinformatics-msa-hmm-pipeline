import os

from src.training.hmm import HMM
from src.training.alignment import msa
from src.generation.generate import generate_datasets, generate_random_gene

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
        with open(f"datasets/dataset_{char}.txt", "r") as f:
            for line in f:
                datasets[char].append(line.strip())

    msa_a = msa(datasets['a'])

    # Build the Initial Draft HMM (from Dataset A)
    hmm = HMM.from_msa(msa_a)

    # Viterbi Training and Generating the Trained HMM
    hmm.train(datasets['b'])

    dataset_c_scores = []
    for gene in datasets['c']:
        score = hmm.forward(gene)
        dataset_c_scores.append(score)

    print("Dataset C scores:", dataset_c_scores)
    print("")

    random_genes = [generate_random_gene() for _ in range(40)]
    random_genes_scores = []
    for gene in random_genes:
        score = hmm.forward(gene)
        random_genes_scores.append(score)

    print("Random genes scores:", random_genes_scores)


if __name__ == "__main__":
    main()
