#!/usr/bin/env python
import sys
import argparse
def main():
    parser = argparse.ArgumentParser(prog='globalalign',description='find the global alignment of two sequences')
    parser.add_argument('sequence1', type = str, help = 'the first sequence') #file location
    parser.add_argument('sequence2', type = str, help = 'the second sequence')
    parser.add_argument('--match', type = int, default = 1)
    parser.add_argument('--mismatch', type = int, default = -1)
    parser.add_argument('--gap', type = int, default = -1, help = 'Gap penalty')
    parser.add_argument('--gapExtension', type = int, default = -1, help = 'Gap Extension penalty')
    parser.add_argument('--similarity_matrix', type = str)

    args = parser.parse_args()
    filename1 = args.sequence1
    filename2 = args.sequence2
    match = args.match
    mismatch = args.mismatch
    indel = args.gap
    extension = args.gapExtension
 
    def similarityMatrix(a,b): #basic similarity matrix
        if(a==b):
            return match
        else:
            return mismatch

    with open(filename1, 'r') as reader:
        sequence1 = reader.readline()
        sequenceLength1 = len(sequence1)
        if(sequence1[len(sequence1)-1] == '\n'):
            sequenceLength1 = sequenceLength1 - 1
            sequence1 = sequence1[:sequenceLength1]
        print("Sequence 1")
        print(sequence1)

    with open(filename2, 'r') as reader:
        sequence2 = reader.readline()
        sequenceLength2 = len(sequence2)
        if(sequence2[len(sequence2)-1] == '\n'):
            sequenceLength2 = sequenceLength2 - 1
            sequence2 = sequence2[:sequenceLength2]
        print("Sequence 2")
        print(sequence2)

    rows, cols = (sequenceLength2 + 1, sequenceLength1 + 1)
    Fmatrix = [[0 for i in range(cols)] for j in range(rows)]
    Mmatrix = [[0 for i in range(cols)] for j in range(rows)]
    Ixmatrix = [[0 for i in range(cols)] for j in range(rows)]
    Iymatrix = [[0 for i in range(cols)] for j in range(rows)]
    if (indel != extension):
        for i in range(0, sequenceLength2 + 1):
            Fmatrix[i][0] = indel + (i-1)*extension
            Mmatrix[i][0] = float("-inf")
            Ixmatrix[i][0] = float("-inf")
            Iymatrix[i][0] = indel + i * extension
        for j in range(0, sequenceLength1 + 1):
            Fmatrix[0][j] = indel + (j-1)*extension
            Mmatrix[0][j] = float("-inf")
            Ixmatrix[0][j] = indel + j * extension
            Iymatrix[0][j] = float("-inf")
        for i in range(1, sequenceLength2 + 1):
            for j in range(1, sequenceLength1 + 1):
                    M = Mmatrix[i-1][j-1] + similarityMatrix(sequence2[i-1], sequence1[j-1])
                    D = Ixmatrix[i-1][j-1] + similarityMatrix(sequence2[i-1], sequence1[j-1])
                    I = Iymatrix[i-1][j-1] + similarityMatrix(sequence2[i-1], sequence1[j-1])
                    Mmatrix[i][j] = max(M,D,I)
                    Ixmatrix[i][j] = max(Mmatrix[i][j-1] + indel + extension, Ixmatrix[i][j-1] + extension, Iymatrix[i][j-1] + indel + extension)
                    Iymatrix[i][j] = max(Mmatrix[i-1][j] + indel + extension, Iymatrix[i-1][j] + extension, Ixmatrix[i-1][j] + indel + extension)
                    Fmatrix[i][j] = max(Mmatrix[i][j], Ixmatrix[i][j], Iymatrix[i][j])
    else:
        for i in range(0, sequenceLength2 + 1):
            Fmatrix[i][0] = indel * i
        for j in range(0, sequenceLength1 + 1):
            Fmatrix[0][j] = indel * j
        for i in range(1, sequenceLength2 + 1):
            for j in range(1, sequenceLength1 + 1):
                    M = Fmatrix[i-1][j-1] + similarityMatrix(sequence2[i-1], sequence1[j-1])
                    D = Fmatrix[i-1][j] + indel
                    I = Fmatrix[i][j-1] + indel
                    Fmatrix[i][j] = max(M,D,I)

    AllignmentA = ""
    AllignmentB = ""
    i = sequenceLength2
    j = sequenceLength1
    #find optimal allignment, there may be more than 1 but this only finds 1
    while (i > 0 or j > 0):
        #if there is a match
        if (i > 0 and j > 0 and Fmatrix[i][j] == (Fmatrix[i-1][j-1] + similarityMatrix(sequence2[i-1],sequence1[j-1]))):
            AllignmentA = sequence2[i - 1] + AllignmentA
            AllignmentB = sequence1[j - 1] + AllignmentB
            i = i - 1
            j = j - 1
        #if there is a deletion then match allignmentA with a gap
        elif (i > 0 and Fmatrix[i][j] == Fmatrix[i-1][j] + indel):
            AllignmentA = sequence2[i - 1] + AllignmentA
            AllignmentB = "-" + AllignmentB
            i = i - 1
        #if there is a insertion then match allignmentB with a gap
        else:
            AllignmentA = "-" + AllignmentA
            AllignmentB = sequence1[j - 1] + AllignmentB
            j = j - 1
    print("Global Allignment")
    print(AllignmentB)
    print(AllignmentA)
    with open ("Output.txt", "w") as output_file:
        output_file.writelines([AllignmentB,"\n", AllignmentA])
if __name__ == '__main__':
    main()