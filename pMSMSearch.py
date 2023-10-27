import sys, os,\
    math, time, csv
import numpy as np

C = 0.1

def fill_with_inf(start, end, array):
    for i in range(start, end):
        array[i] = math.inf
    return array

def get_lower_bound(m, n, xCoord, yCoord):
    return abs(m-n-xCoord+yCoord) * C

def Cost(newp, x, y):
    if newp < min(x,y) or newp > max(x,y):
        return C + min(abs(newp - x), abs(newp - y))
    return C

def msm_dist_pruned(ts1, ts2, n, m, bsf):
    tmp_array = np.full(n + 1, float('inf'))
    tmp = 0
    sc = 1
    ec = 1
    for i in range(1, m + 1):
        smaller_found = False
        xi = ts1[i]
        ec_next = i
        for j in range(sc, len(tmp_array)):
            yj = ts2[j]
            d1 = tmp + abs(xi - yj)
            d2 = tmp_array[j] + Cost(xi, ts1[i - 1], yj)
            d3 = tmp_array[j - 1] + Cost(yj, xi, ts2[j - 1])
            tmp = tmp_array[j]
            tmp_array[j] = min(d1, d2, d3)
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
    return [tmp_array[-1]]

def main(*args):
    m, r = -1, -1
    print(args)
    m = int(args[3])
    q = []
    bsf = math.inf
    bc = math.inf
    print("Read files")
    if os.path.getsize(args[1]) == 0:
        raise Exception("Sequence Path is empty")
    if os.path.getsize(args[2]) == 0:
        raise Exception("Query Path is empty")
    with open(args[2], 'r') as file:
        for value in file:
            query = value.strip().split('\t')
            for i, value in enumerate(query):
                query[i] = float(value)
    print(query)
    with open(args[1], 'r') as file:
        lines = file.readline()

        tsv_file = csv.reader(file, delimiter="\t")
        print ("tsv: ", tsv_file)
        for sequence in tsv_file:
            print("sequence=", sequence)
            #Convert string to float
            for i, value in enumerate(sequence):
                sequence[i] = float(value)
            print(sequence)
            print('Sequence read')
            tempdist = msm_dist_pruned(sequence[1:], query[1:], len(query)-2, len(sequence)-2, bsf)
            print(tempdist[0])
            if bsf > tempdist[0]:
                bsf = tempdist
            bc = sequence[0]

    start_time = time.time()
    print(f"1NN class: {bc}")
    print(f"Query class: {query[0]}")
    for i, arg in enumerate(args, 1):
        print(f"Argument {i}: {arg}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python test.py <filepath> <querypath> <querysize>")
        sys.exit(1)

    main(*sys.argv[0:])