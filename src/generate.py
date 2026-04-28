from random import choice, randint

alphabet = {"A", "C", "G", "T"}
patterns = ["ΑGAΤΤΑ", "ΑCTTTGCΑ", "ΑACTGGCAA", "CATTCAGT"]

a_min, a_max = 1, 3
b_min, b_max = 1, 2
c_min, c_max = 1, 2


def generate_string():
    out = ""

    a_count = randint(a_min, a_max)
    for _ in range(a_count):
        out += choice(list(alphabet))

    for pattern in patterns:
        mutated = pattern
        b_count = randint(b_min, b_max)
        for _ in range(b_count):
            index = randint(0, len(mutated) - 1)
            mutated = (
                mutated[:index]
                + choice(list(alphabet | {""}))
                + mutated[index + 1 :]
            )
        out += mutated

    c_count = randint(c_min, c_max)
    for _ in range(c_count):
        out += choice(list(alphabet))

    return out
