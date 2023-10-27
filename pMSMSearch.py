import sys, os,\
    math, time
import numpy as np

def msm_dist_pruned(self):
    tmp_array = np.full(self.n + 1, float('inf'))
    tmp = 0
    sc = 1
    ec = 1
    for i in range(1, self.m + 1):
        smaller_found = False
        xi = self.ts1[i]
        ec_next = i
        for j in range(sc, len(tmp_array)):
            yj = self.ts2[j]
            d1 = tmp + abs(xi - yj)
            d2 = tmp_array[j] + self.C(xi, self.ts1[i - 1], yj)
            d3 = tmp_array[j - 1] + self.C(yj, xi, self.ts2[j - 1])
            tmp = tmp_array[j]
            tmp_array[j] = min(d1, d2, d3)
            lb = self.get_lower_bound(i, j)
            if (tmp_array[j] + lb) > self.upper_bound:
                if not smaller_found:
                    sc = j + 1
                if j > ec:
                    tmp_array = self.fill_with_inf(j + 1, self.n + 1, tmp_array)
                    break
            else:
                smaller_found = True
                ec_next = j + 1
        tmp_array = self.fill_with_inf(1, sc, tmp_array)
        tmp = float('inf')
        ec = ec_next
    return [tmp_array[-1]]

def main(*args):
    m, r = -1, -1
    print(args)
    m = int(args[3])
    q = []
    print("Read files")
    if os.path.getsize(args[1]) == 0:
        raise Exception("Sequence Path is empty")
    if os.path.getsize(args[2]) == 0:
        raise Exception("Query Path is empty")
    with open(args[1], 'r') as file:
        fp = file.read()
    with open(args[2], 'r') as file:
        for line in file:
            values = file[0].strip().split('\t')
            if len(values) >= 2:
                q.append(values[1:])


    start_time = time.time()

    for i, arg in enumerate(args, 1):
        print(f"Argument {i}: {arg}")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python test.py <filepath> <querypath> <querysize>")
        sys.exit(1)

    main(*sys.argv[0:])