from random import choice, randint

ALPHABET = {"A", "C", "G", "T"}
PATTERNS = ["AGATTA", "ACTTTGCA", "AACTGGCAA", "CATTCAGT"]

A_MIN, A_MAX = 1, 3
B_MIN, B_MAX = 0, 2
C_MIN, C_MAX = 1, 3


def generate_string():
    output = ""

    a_count = randint(A_MIN, A_MAX)
    for _ in range(a_count):
        output += choice(list(ALPHABET))

    for pattern in PATTERNS:
        mutated = pattern
        b_count = randint(B_MIN, B_MAX)
        for _ in range(b_count):
            index = randint(0, len(mutated) - 1)
            mutated = (
                mutated[:index]
                + choice(list(ALPHABET | {""}))
                + mutated[index + 1 :]
            )
        output += mutated

    c_count = randint(C_MIN, C_MAX)
    for _ in range(c_count):
        output += choice(list(ALPHABET))

    return output
