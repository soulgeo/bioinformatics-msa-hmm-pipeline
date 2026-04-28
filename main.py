import random
from collections import Counter

from src.generate import generate_string

A_COUNT = 20
B_COUNT = 140
C_COUNT = 40
GENE_COUNT = A_COUNT + B_COUNT + C_COUNT


def remove_sampled_elements(original_list, sampled_list):
    counts = Counter(original_list)
    to_remove = Counter(sampled_list)

    remaining_counts = counts - to_remove

    return list(remaining_counts.elements())


def main():
    genes = []
    for _ in range(GENE_COUNT):
        genes += [generate_string()]

    dataset_a = random.sample(genes, A_COUNT)
    genes = remove_sampled_elements(genes, dataset_a)

    dataset_b = random.sample(genes, B_COUNT)
    genes = remove_sampled_elements(genes, dataset_b)

    dataset_c = random.sample(genes, C_COUNT)

    print(f"Dataset A: {dataset_a}")
    print(f"Dataset B: {dataset_b}")
    print(f"Dataset C: {dataset_c}")


if __name__ == "__main__":
    main()
