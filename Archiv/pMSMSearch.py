import sys, os,\
    math, time, csv, pdb
import numpy as np

C = 0.1
UCRARCHIVE_PATH = ".\\UCRArchive\\"

def fill_with_inf(start, end, array):
    for i in range(start, end-1):
        array[i] = math.inf
    return array

def get_lower_bound(m, n, xCoord, yCoord):
    return abs(m-n-xCoord+yCoord) * C

def cost(newPoint, x, y):
    if (x <= newPoint and newPoint <= y) or (y <= newPoint and newPoint <= x):
        return C
    else:
        return C + min(abs(newPoint - x), abs(newPoint - y))

def msm_dist_pruned(ts1, ts2, m, n, bsf):
    tmpArray = np.full(n+1, float('inf'))
    tmp = 0
    sc = 1
    ec = 1
    ts1[0], ts2[0]= math.inf, math.inf
    i=1
    while i < m+1:
        smallerFound = False
        xi = ts1[i]
        ec_next = i
        j=sc
        while j < n+1:
            yj = ts2[j]
            d1 = tmp + abs(xi - yj)
            d2 = tmpArray[j] + cost(xi, ts1[i - 1], yj)
            d3 = tmpArray[j - 1] + cost(yj, xi, ts2[j - 1])
            tmp = tmpArray[j]
            tmpArray[j] = min(d1, min(d2, d3))
            lb = get_lower_bound(m, n, i, j)
            if (tmpArray[j] + lb) > bsf:
                if not smallerFound:
                    sc = j + 1
                if j > ec:
                    tmpArray = fill_with_inf(j + 1, n + 1, tmpArray)
                    break
            else:
                smallerFound = True
                ec_next = j + 1
            j += 1
        tmpArray = fill_with_inf(1, sc, tmpArray)
        tmp = float('inf')
        ec = ec_next
        i += 1
    return tmpArray[n]

def msm_dist(sequence, query, m, n):
    i=1
    j=1
    costArray = np.full((m,n), float('inf'))
    costArray[0,0] = abs(sequence[0] - query[0])
    while i < m:
        costArray[i][0] = costArray[i-1][0] + cost(query[i], query[i - 1], sequence[0])
        i += 1
    while j < n:
        costArray[0][j] = costArray[0][j-1] + cost(sequence[j], query[0], sequence[j - 1])
        j += 1
    i=1
    j=1
    while i < m:
        j=1
        while j < n:

            d1 = costArray[i-1][j-1] + abs(query[i] - sequence[j])
            d2 = costArray[i-1][j] + cost(query[i], query[i - 1], sequence[j])
            d3 = costArray[i][j-1] + cost(sequence[j], query[i], sequence[j - 1])
            costArray[i][j] = min(d1, min(d2, d3))
            j += 1
        i += 1
    return costArray[m-1][n-1]

def main(*args):
    class_train = []
    class_test = []
    true_pred = 0
    query_count = 0

    if os.path.getsize(UCRARCHIVE_PATH + args[1] + "\\" + args[1] + "_TRAIN.tsv") == 0:
        raise Exception("Sequence Path is empty")
    if os.path.getsize(UCRARCHIVE_PATH + args[1] + "\\" + args[1] + "_TEST.tsv") == 0:
        raise Exception("Query Path is empty")
    start_time = time.time()
    with open(UCRARCHIVE_PATH + args[1] + "\\" + args[1] + "_TEST.tsv") as file:
        for value in file:
            query_count += 1
            query = value.strip().split('\t')
            #Convert string to float
            for i, value in enumerate(query):
                query[i] = float(value)
            q_class = query[0]
            class_test.append(query[0])
            print(f"Class: {query[0]}, Query: {query[1:]}")
            #initialize bsf and pred_class before every query
            bsf = math.inf
            pred_class = math.inf
            with open(UCRARCHIVE_PATH + args[1] + "\\" + args[1] + "_TRAIN.tsv") as file:
                tsv_file = csv.reader(file, delimiter="\t")
                for idx, sequence in enumerate(tsv_file):
                    #Convert string to float
                    for j, value in enumerate(sequence):
                        sequence[j] = float(value)
                    s_class = sequence[0]
                    class_train.append(query[0])
                    pred_dist = msm_dist_pruned(sequence, query, len(query)-1, len(sequence)-1, bsf)
                    #pred_dist = msm_dist(sequence[1:], query[1:], len(query)-1, len(sequence)-1)
                    if bsf > pred_dist:
                        print("Better found at: ", idx, " with class: ", s_class, " class should be: ", q_class)
                        print("pred_dist: ", pred_dist, " bsf: ", bsf)
                        bsf = pred_dist
                        pred_class = s_class
                if pred_class == q_class:
                    true_pred += 1
                print(f"1NN class: {pred_class}")
                print(f"Query class: {q_class}")
        acc = true_pred/query_count
        print(f"TP: {true_pred}")
        print(f"Accuracy: {acc}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python test.py <dataset>")
        sys.exit(1)
    main(*sys.argv[0:])