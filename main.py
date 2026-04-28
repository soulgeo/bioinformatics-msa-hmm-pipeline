from src.generate import generate_datasets

A_COUNT = 20
B_COUNT = 140
C_COUNT = 40


def create_datasets():
    datasets = generate_datasets(A_COUNT, B_COUNT, C_COUNT)
    for i, dataset in enumerate(datasets):
        letter = chr(ord('`') + i + 1)
        with open(f"datasets/dataset_{letter}.txt", "w") as f:
            output = "\n".join(dataset)
            f.write(output)


def main():
    create_datasets()


if __name__ == "__main__":
    main()
