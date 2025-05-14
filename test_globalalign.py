import pytest
from globalalign import *


test_aligner = GlobalAlign(1, -1, -1)


def test_align_sequences1():
    """test align_sequences gives optimal alignments"""
    sequence1 = "GCATGCU"
    sequence2 = "GATTACA"
    alignments = set(test_aligner.align_sequences(sequence1, sequence2)[0])
    assert alignments == {
        ("GCA-TGCU", "G-ATTACA"),
        ("GCAT-GCU", "G-ATTACA"),
        ("GCATG-CU", "G-ATTACA"),
    }


def test_align_sequences2():
    """test align_sequences gives optimal alignments"""
    sequence1 = "GCATGCG"
    sequence2 = "GATTACA"
    alignments = set(test_aligner.align_sequences(sequence1, sequence2)[0])
    assert alignments == {
        ("GCAT-GCG", "G-ATTACA"),
        ("GCATG-CG", "G-ATTACA"),
        ("GCA-TGCG", "G-ATTACA"),
    }


def test_align_sequences3():
    """test align_sequences gives correct alignment score"""
    sequence1 = "GCATGCU"
    sequence2 = "GATTACA"
    score = test_aligner.align_sequences(sequence1, sequence2)[1]
    assert score == 0


def test_align_sequences4():
    """test align_sequences gives correct alignment score"""
    sequence1 = "GCATGCG"
    sequence2 = "GATTACA"
    score = test_aligner.align_sequences(sequence1, sequence2)[1]
    assert score == 0
