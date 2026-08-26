import sys
import csv
from tabulate import tabulate


def main():
    samples_fp = sys.argv[1]
    metadata_fp = sys.argv[2]

    metadata = {}

    for row in csv.reader(open(metadata_fp, newline="", encoding="utf-8")):
        metadata[row[0]] = row[3]

    samples = []
    for line in open(samples_fp):
        sample = line.strip()
        samples.append([sample, metadata[sample]])
    print(len(samples))

    latex = tabulate(
        samples,
        tablefmt="latex",
        showindex=True,  # range(1, len(samples) + 1),
    )
    print(latex)


if __name__ == "__main__":
    main()
