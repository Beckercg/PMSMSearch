from Archiv.pMSMSearch import msm_dist
import csv


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
    pred_msm_dist = msm_dist(sequence[1:], query[1:], len(sequence)-1, len(query)-1)
    return pred_msm_dist

def test_msm_simple_same():
    assert run_msm(0,0) == 0.0
def test_pruned_msm_simple_calculated():
    assert run_msm(0,1) == 4.4