#!/usr/bin/env python
import sys

def main():
    if len(sys.argv) < 3: #check to make sure 2 files are specified
        print("Need to specify file names")
        return
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]
    match = 1
    mismatch = -1
    indel = -1
    if len(sys.argv) == 6: #alter scoring scheme if user specifies
        match = int(sys.argv[3])
        mismatch = int(sys.argv[4])
        indel = int(sys.argv[5])
 
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
if __name__ == '__main__':
    main()