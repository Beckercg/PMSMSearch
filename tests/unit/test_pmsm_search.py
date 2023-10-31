from pMSMSearch import msm_dist_pruned, msm_dist
import sys, csv, math


def run_msm(queryNumber, sequenceNumber):
    sequence = []
    query = []
    with open(".\\tests\\unit\\data_testing.tsv") as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for i, line in enumerate(tsv_file):
            for j, value in enumerate(line):
                line[j] = float(value)
            if i == queryNumber:
                query = line
            if i == sequenceNumber:
                sequence = line
        pred_pruned_msm_dist = msm_dist_pruned(sequence, query, len(sequence)-1, len(query)-1, math.inf)
        pred_msm_dist = msm_dist(sequence[1:], query[1:], len(sequence)-1, len(query)-1)
        return [pred_pruned_msm_dist, pred_msm_dist]

def test_pruned_msm_simple_same():
    assert run_msm(0,0)[0] == 0
def test_pruned_msm_simple_calculated():
    assert run_msm(0,1)[0] == 4.4
def test_pruned_msm_simple_compare():
    assert run_msm(0,1)[0] == run_msm(0,1)[1]