import os
from src.alignment import msa
from src.generate import generate_datasets
from src.hmm import create_states_sequence

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

    aligned_dataset_a = msa(dataset_a)
    for seq in aligned_dataset_a:
        print(seq)

    states = create_states_sequence(aligned_dataset_a)
    print(''.join(states))


if __name__ == "__main__":
    main()
