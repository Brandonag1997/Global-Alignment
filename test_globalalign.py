import pytest
from globalalign import *


test_alligner = GlobalAlign(1, -1, -1)


def test_align_sequences1():
    sequence1 = "GCATGCU"
    sequence2 = "GATTACA"
    alignments = set(test_alligner.align_sequences(sequence1, sequence2))
    assert alignments == {
        ("GCA-TGCU", "G-ATTACA"),
        ("GCAT-GCU", "G-ATTACA"),
        ("GCATG-CU", "G-ATTACA"),
    }


def test_align_sequences2():
    sequence1 = "GCATGCG"
    sequence2 = "GATTACA"
    alignments = set(test_alligner.align_sequences(sequence1, sequence2))
    assert alignments == {
        ("GCAT-GCG", "G-ATTACA"),
        ("GCATG-CG", "G-ATTACA"),
        ("GCA-TGCG", "G-ATTACA"),
    }
