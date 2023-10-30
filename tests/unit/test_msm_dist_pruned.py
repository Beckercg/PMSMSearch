from pMSMSearch import msm_dist_pruned, msm_dist
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
    pred_pruned_msm_dist = msm_dist_pruned(sequence, query, len(sequence)-1, len(query)-1, math.inf)
    pred_msm_dist = msm_dist(sequence[1:], query[1:], len(sequence)-1, len(query)-1)
    assert pred_pruned_msm_dist == pred_msm_dist
