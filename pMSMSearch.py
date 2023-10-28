import sys, os,\
    math, time, csv, pdb
import numpy as np

C = 0.1
UCRARCHIVE_PATH = ".\\UCRArchive_2018\\UCRArchive_2018\\"

def fill_with_inf(start, end, array):
    for i in range(start, end-1):
        array[i] = math.inf
    return array

def get_lower_bound(m, n, xCoord, yCoord):
    return abs(m-n-xCoord+yCoord) * C

def Cost(newp, x, y):
    if newp <= min(x,y) or newp >= max(x,y):
        return C + min(abs(newp - x), abs(newp - y))
    return C

def msm_dist_pruned(ts1, ts2, m, n, bsf):
    tmp_array = np.full(n, float('inf'))
    tmp = 0
    sc = 1
    ec = 1
    ts1[0], ts2[0]= math.inf, math.inf
    for i in range(1, m):
        smaller_found = False
        xi = ts1[i]
        ec_next = i
        for j in range(sc, n):
            yj = ts2[j]
            d1 = tmp + abs(xi - yj)
            d2 = tmp_array[j] + Cost(xi, ts1[i - 1], yj)
            d3 = tmp_array[j - 1] + Cost(yj, xi, ts2[j - 1])
            tmp = tmp_array[j]
            tmp_array[j] = min(d1, min(d2, d3))
            lb = get_lower_bound(m, n, i, j)
            if (tmp_array[j] + lb) > bsf:
                if not smaller_found:
                    sc = j + 1
                if j > ec:
                    tmp_array = fill_with_inf(j + 1, n + 1, tmp_array)
                    break
            else:
                smaller_found = True
                ec_next = j + 1
        tmp_array = fill_with_inf(1, sc, tmp_array)
        tmp = float('inf')
        ec = ec_next
    return tmp_array[n-1]

def msm_dist(sequence, query, m, n):
    cost_array = np.full((m,n), float('inf'))
    cost_array[0,0] = abs(sequence[0] - query[0])
    for i in range(1, m-1):
        cost_array[i][0] = cost_array[i-1][0] + Cost(sequence[i], sequence[i-1], query[0])
    for j in range(1, n-1):
        cost_array[j][0] = cost_array[j-1][0] + Cost(query[i], sequence[0], query[j-1])

    for i in range(1, m-1 ):
        for j in range(1, n-10):
            d1 = cost_array[i-1][j-1] + abs(sequence[i] - query[j])
            d2 = cost_array[i-1][j] + Cost(sequence[i], sequence[i - 1], query[j])
            d3 = cost_array[i][j-1] + Cost(query[j], sequence[i], query[j-1])
            cost_array[j] = min(d1, min(d2, d3))
    return cost_array[m-1][n-1]

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
        query_count += 1
        for value in file:
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
                for sequence in tsv_file:
                    #Convert string to float
                    for i, value in enumerate(sequence):
                        sequence[i] = float(value)
                    s_class = sequence[0]
                    class_train.append(query[0])
                    pred_dist = msm_dist_pruned(sequence, query, len(query), len(sequence), bsf)
                    #pred_dist = msm_dist(sequence[1:], query[1:], len(query)-1, len(sequence)-1)
                    if bsf > pred_dist:
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