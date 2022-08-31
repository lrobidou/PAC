# generate data for unit tests of pac
import random as rd

def rd_kmer(k):
    return "".join(rd.choice("ACTG") for _ in range(k))


def random_kmc_files(filename, nb_file, nb_line, k, max_value):
    """Write i files in tests/unit/data/{filename}_{i}.txt and a file tests/unit/data/{filename}_fof.txt"""
    for i in range(nb_file):
        with open(f"tests/unit/data/{filename}_{i}.txt", 'w') as fichier:
            for _ in range(nb_line):
                fichier.write(f"{rd_kmer(k)}\t{rd.randint(1, max_value)}\n")
    with open(f"tests/unit/data/{filename}_fof.txt", "w") as fichier:
        for i in range(nb_file):
            fichier.write(f"tests/unit/data/{filename}_{i}.txt\n")

    

def main():
    random_kmc_files("10000_31-mers", 3, 10000, 31, 20)


if __name__ == '__main__':
    main()