from pMSMSearch import msm_dist
import sys, csv, math


def test_pruned_msm_simple():
    sequences = []
    query = []
    with open(".\\tests\\unit\\data_testing.tsv") as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for i, sequence in enumerate(tsv_file):
            for j, value in enumerate(sequence):
                sequence[j] = float(value)
            if i == 0:
                query = sequence
            else:
                sequences.append(sequence)
    pred_msm_dist = msm_dist(sequence[1:], query[1:], len(sequence)-1, len(query)-1)
    assert pred_msm_dist == 4.4
