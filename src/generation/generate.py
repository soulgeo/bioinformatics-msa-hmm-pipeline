from collections import Counter
import random

ALPHABET = {"A", "C", "G", "T"}
PATTERNS = ["AGATTA", "ACTTTGCA", "AACTGGCAA", "CATTCAGT"]

A_MIN, A_MAX = 1, 3
B_MIN, B_MAX = 0, 2
C_MIN, C_MAX = 1, 3

RANDOM_GENE_LENGTH = 60
VARIANCE = 3


def generate_gene_from_patterns():
    gene = []

    # a) Random symbols at the beginning
    a_count = random.randint(A_MIN, A_MAX)
    for _ in range(a_count):
        gene.append(random.choice(list(ALPHABET)))

    # b) Process patterns with mutations
    a_count = random.randint(A_MIN, A_MAX)
    for pattern in PATTERNS:
        mutated = pattern.split()
        b_count = random.randint(B_MIN, B_MAX)
        for _ in range(b_count):
            index = random.randint(0, len(mutated) - 1)
            mutated[index] = random.choice(list(ALPHABET | {""}))

        gene.append(mutated)

    # c) Random symbols at the end
    c_count = random.randint(C_MIN, C_MAX)
    for _ in range(c_count):
        gene.append(random.choice(list(ALPHABET)))

    return ''.join(gene)


def remove_sampled_genes(original_list, sampled_list):
    counts = Counter(original_list)
    to_remove = Counter(sampled_list)

    remaining_counts = counts - to_remove

    return list(remaining_counts.elements())


def generate_datasets(*args: int):
    genes = []
    for _ in range(sum(args)):
        genes.append(generate_gene_from_patterns())

    datasets = []
    for i, arg in enumerate(args):
        datasets.append(random.sample(genes, arg))
        genes = remove_sampled_genes(genes, datasets[i])

    return datasets


def generate_random_gene():
    gene = []
    variance = random.randint(-VARIANCE, VARIANCE)
    for _ in range(RANDOM_GENE_LENGTH + variance):
        gene.append(random.choice(list(ALPHABET)))

    return ''.join(gene)
