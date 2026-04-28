from collections import Counter
import random

ALPHABET = {"A", "C", "G", "T"}
PATTERNS = ["AGATTA", "ACTTTGCA", "AACTGGCAA", "CATTCAGT"]

A_MIN, A_MAX = 1, 3
B_MIN, B_MAX = 0, 2
C_MIN, C_MAX = 1, 3


def generate_gene():
    output = ""

    # a) random symbols at the beginning
    a_count = random.randint(A_MIN, A_MAX)
    for _ in range(a_count):
        output += random.choice(list(ALPHABET))

    # b) process patterns with mutations
    a_count = random.randint(A_MIN, A_MAX)
    for pattern in PATTERNS:
        mutated = pattern
        b_count = random.randint(B_MIN, B_MAX)
        for _ in range(b_count):
            index = random.randint(0, len(mutated) - 1)
            mutated = (
                mutated[:index]
                + random.choice(list(ALPHABET | {""}))
                + mutated[index + 1 :]
            )
        output += mutated

    # c) random symbols at the end
    c_count = random.randint(C_MIN, C_MAX)
    for _ in range(c_count):
        output += random.choice(list(ALPHABET))

    return output


def remove_sampled_genes(original, sampled):
    counts = Counter(original)
    to_remove = Counter(sampled)

    remaining_counts = counts - to_remove

    return list(remaining_counts.elements())


def generate_datasets(*args: int):
    genes = []
    for _ in range(sum(args)):
        genes += [generate_gene()]

    datasets = []
    for i, val in enumerate(args):
        datasets += [random.sample(genes, val)]
        genes = remove_sampled_genes(genes, datasets[i])

    return datasets
