from random import choice, randint

ALPHABET = {"A", "C", "G", "T"}
PATTERNS = ["AGAΤΤA", "ACTTTGCA", "AACTGGCAA", "CATTCAGT"]

a_min, a_max = 1, 3
b_min, b_max = 0, 2
c_min, c_max = 1, 2


def generate_string():
    output = ""

    a_count = randint(a_min, a_max)
    for _ in range(a_count):
        output += choice(list(ALPHABET))

    for pattern in PATTERNS:
        mutated = pattern
        b_count = randint(b_min, b_max)
        for _ in range(b_count):
            index = randint(0, len(mutated) - 1)
            mutated = (
                mutated[:index]
                + choice(list(ALPHABET | {""}))
                + mutated[index + 1 :]
            )
        output += mutated

    c_count = randint(c_min, c_max)
    for _ in range(c_count):
        output += choice(list(ALPHABET))

    return output
