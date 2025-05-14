#!/usr/bin/env python
import argparse
import os


class GlobalAlign:

    def __init__(
        self,
        match: int = 1,
        mismatch: int = -1,
        gap: int = -1,
        similarity_matrix_path: str = None,
    ):
        """sets initial values"""
        self.match = match
        self.mismatch = mismatch
        self.indel = gap
        self.matrix = self.parse_matrix(similarity_matrix_path)
        self.f_matrix = None

    def parse_matrix(self, similarity_matrix_path):
        """Parse similarity_matrix_path and return matrix"""
        matrix = None
        if similarity_matrix_path is not None and os.path.exists(
            similarity_matrix_path
        ):
            matrix = {}
            with open(similarity_matrix_path, mode="r", encoding="utf-8") as reader:
                lines = reader.readlines()
                for l in lines:
                    if l[0] == "#":  # skip comments
                        continue

                    if l[0] == " ":
                        columns = l.split()
                    else:
                        scores = l.split()
                        if len(scores) == 0:
                            continue
                        letter = scores[0]
                        for index, item in enumerate(columns):
                            matrix[(letter, item)] = scores[index + 1]
        return matrix

    def update_matrix(self, similarity_matrix_path):
        """Updates the GlobalAlign similarity_matrix"""
        self.matrix = self.parse_matrix(similarity_matrix_path)

    def similarity_matrix(self, a, b):
        """Give the results of a similarity matrix for 2 given elements.
        Defaults to basic similarity matrix"""
        if a == "U":
            a = "T"
        if b == "U":
            b = "T"
        if self.matrix is None:  # basic similarity matrix
            if a == b:
                return self.match
            return self.mismatch
        if a == " " or b == " ":
            return self.mismatch
        return int(self.matrix.get((a, b), self.mismatch))

    def get_f_matrix(self, sequence1: str, sequence2: str):
        """Calculate F matrix"""
        rows, cols = (len(sequence2) + 1, len(sequence1) + 1)
        f_matrix = [
            [0 for i in range(cols)] for j in range(rows)
        ]  # fill matrix with 0s
        for i in range(0, len(sequence2) + 1):  # assign gap penalty
            f_matrix[i][0] = self.indel * i
        for j in range(0, len(sequence1) + 1):  # assign gap penalty
            f_matrix[0][j] = self.indel * j
        for i in range(1, len(sequence2) + 1):
            for j in range(1, len(sequence1) + 1):
                # match or mismatch
                match = f_matrix[i - 1][j - 1] + self.similarity_matrix(
                    sequence2[i - 1], sequence1[j - 1]
                )
                # gap insertion or deletion
                deletion = f_matrix[i - 1][j] + self.indel
                insertion = f_matrix[i][j - 1] + self.indel
                f_matrix[i][j] = max(match, deletion, insertion)
        return f_matrix

    def align_sequences(self, sequence1: str, sequence2: str):
        """Align 2 sequences and set the alignments as the results"""
        self.f_matrix = self.get_f_matrix(sequence1, sequence2)
        alignments = self.find_alignment(sequence1, sequence2)
        score = self.f_matrix[len(sequence2)][len(sequence1)]
        return (alignments, score)

    def find_alignment(self, sequence1: str, sequence2: str):
        """function to find alignment of 2 sequences"""
        i = len(sequence2)
        j = len(sequence1)
        stack = [("", "", i, j)]
        results = []

        while stack:
            alignment_a, alignment_b, i, j = stack.pop()

            if i == 0 and j == 0:
                results.append((alignment_b, alignment_a))
                continue

            current_score = self.f_matrix[i][j]

            if (
                i > 0
                and j > 0
                and current_score
                == self.f_matrix[i - 1][j - 1]
                + self.similarity_matrix(sequence2[i - 1], sequence1[j - 1])
            ):  # match or mismatch
                stack.append(
                    (
                        sequence2[i - 1] + alignment_a,
                        sequence1[j - 1] + alignment_b,
                        i - 1,
                        j - 1,
                    )
                )
            # if there is a deletion then match alignmentA with a gap
            if i > 0 and current_score == self.f_matrix[i - 1][j] + self.indel:
                stack.append(
                    (sequence2[i - 1] + alignment_a, "-" + alignment_b, i - 1, j)
                )
            # if there is a insertion then match alignmentB with a gap
            if j > 0 and current_score == self.f_matrix[i][j - 1] + self.indel:
                stack.append(
                    ("-" + alignment_a, sequence1[j - 1] + alignment_b, i, j - 1)
                )
        return results


def read_sequence(file_path):
    """Read a FASTQ sequnce from a text file"""
    with open(file_path, mode="r", encoding="utf-8") as reader:
        sequence = reader.readline().strip()
        return sequence


def main():
    parser = argparse.ArgumentParser(
        prog="globalalign", description="find the global alignment of two sequences"
    )
    parser.add_argument(
        "sequence1", type=str, help="the first sequence"
    )  # file location
    parser.add_argument("sequence2", type=str, help="the second sequence")
    parser.add_argument("--match", type=int, default=1)
    parser.add_argument("--mismatch", type=int, default=-1)
    parser.add_argument("--gap", type=int, default=-1, help="Gap penalty")
    parser.add_argument("--similarity_matrix", type=str)

    args = parser.parse_args()
    filename1 = args.sequence1
    filename2 = args.sequence2
    match = args.match
    mismatch = args.mismatch
    indel = args.gap
    similarity_matrix = args.similarity_matrix

    alligner = GlobalAlign(match, mismatch, indel, similarity_matrix)
    sequence1 = read_sequence(filename1)
    sequence2 = read_sequence(filename2)
    alignments = set(alligner.align_sequences(sequence1, sequence2))
    print("sequence1")
    print(sequence1)
    print("sequence2")
    print(sequence2)
    for i, alignment in enumerate(alignments):
        print(f"Best Alignment {i+1}")
        print(alignment[0])
        print(alignment[1])


if __name__ == "__main__":
    main()
