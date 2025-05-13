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
        self.results = None
        self.sequence1 = None
        self.sequence2 = None

    def parse_matrix(self, similarity_matrix_path):
        """Parse similarity_matrix_path and return matrix"""
        matrix = None
        if similarity_matrix_path is not None and os.path.exists(
            similarity_matrix_path
        ):
            matrix = {}
            # dirname = os.path.dirname(__file__)
            # filenameS = os.path.join(dirname, 'matrices/',similarity_matrix_path)
            with open(similarity_matrix_path, "r") as reader:
                lines = reader.readlines()
                for l in lines:
                    if l[0] == "#":
                        continue
                    elif l[0] == " ":
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
        self.matrix = self.parse_matrix(similarity_matrix_path)

    def similarityMatrix(self, a, b):
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
        return int(self.matrix.get((a, b)))

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
                match = f_matrix[i - 1][j - 1] + self.similarityMatrix(
                    sequence2[i - 1], sequence1[j - 1]
                )
                # gap insertion or deletion
                deletion = f_matrix[i - 1][j] + self.indel
                insertion = f_matrix[i][j - 1] + self.indel
                f_matrix[i][j] = max(match, deletion, insertion)
        return f_matrix

    def align_sequences(self, sequence1: str, sequence2: str):
        self.f_matrix = self.get_f_matrix(sequence1, sequence2)
        self.results = []
        self.find_allignment(
            "", "", len(sequence2), len(sequence1), 0, sequence1, sequence2
        )
        allignments = self.results
        self.results = None
        return allignments

    def find_allignment(
        self, AllA, AllB, i_x, j_x, score, sequence1: str, sequence2: str
    ):
        if i_x == 0 and j_x == 0:  # write alignment
            self.results.append((AllB, AllA))
            return

        # traceback
        nextScores = [
            self.f_matrix[i_x - 1][j_x] + self.indel,
            self.f_matrix[i_x][j_x - 1] + self.indel,
            self.f_matrix[i_x - 1][j_x - 1]
            + self.similarityMatrix(sequence2[i_x - 1], sequence1[j_x - 1]),
        ]
        mNext = max(nextScores)
        newScore = mNext + score

        # if there is a match or mismatch
        if (
            i_x > 0
            and j_x > 0
            and self.f_matrix[i_x][j_x]
            == self.f_matrix[i_x - 1][j_x - 1]
            + self.similarityMatrix(sequence2[i_x - 1], sequence1[j_x - 1])
        ):
            AllAn = sequence2[i_x - 1] + AllA
            AllBn = sequence1[j_x - 1] + AllB
            self.find_allignment(
                AllAn, AllBn, i_x - 1, j_x - 1, newScore, sequence1, sequence2
            )
        # if there is a deletion then match allignmentA with a gap
        if (
            i_x > 0
            and self.f_matrix[i_x][j_x] == self.f_matrix[i_x - 1][j_x] + self.indel
        ):
            AllAn = sequence2[i_x - 1] + AllA
            AllBn = "-" + AllB
            self.find_allignment(
                AllAn, AllBn, i_x - 1, j_x, newScore, sequence1, sequence2
            )
        # if there is a insertion then match allignmentB with a gap
        if (
            j_x > 0
            and self.f_matrix[i_x][j_x] == self.f_matrix[i_x][j_x - 1] + self.indel
        ):
            AllAn = "-" + AllA
            AllBn = sequence1[j_x - 1] + AllB
            self.find_allignment(
                AllAn, AllBn, i_x, j_x - 1, newScore, sequence1, sequence2
            )


def read_sequence(file_path):
    with open(file_path, "r") as reader:
        sequence = reader.readline()
        if sequence[len(sequence) - 1] == "\n":
            sequence = sequence[:-1]
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
        print(f"Best Allignment {i+1}")
        print(alignment[0])
        print(alignment[1])


if __name__ == "__main__":
    main()
