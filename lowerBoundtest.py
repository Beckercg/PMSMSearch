import sys, os, \
    math, time, csv
import numpy as np

UCRARCHIVE_PATH = ".\\UCRArchive_2018\\"


def main(queryNumber, sequenceNumber):
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
        n = len(query)
        if n > len(sequence):
            n = len(sequence)
        j = 1
        if query[j] < sequence[j]:
            q_gr_s = True
        else:
            q_gr_s = False

        max_is_q = query[0]
        max_is_s = sequence[0]
        min_is_q = query[0]
        min_is_s = sequence[0]
        j += 1
        glb = 0
        while j < n:
            if query[j] < sequence[j]:
                # calculate lowerbound for last intersection
                if not q_gr_s:
                    glb += abs(max_is_q - max_is_s) + abs(min_is_q - min_is_s)
                    print("switch at: ", j)
                    print("before: ", query[j - 1], sequence[j - 1])
                    print("now: ", query[j], sequence[j])
                q_gr_s = True
            else:
                if q_gr_s:
                    glb += abs(max_is_q - max_is_s) + abs(min_is_q - min_is_s)
                    print("switch at: ", j)
                    print("before: ", query[j - 1], sequence[j - 1])
                    print("now: ", query[j], sequence[j])
                q_gr_s = False
            if query[j] > max_is_q:
                max_is_q = query[j]
            if query[j] < min_is_q:
                min_is_q = query[j]
            if sequence[j] > max_is_s:
                max_is_s = sequence[j]
            if sequence[j] < min_is_s:
                min_is_s = sequence[j]
            j += 1
        print("GLB IS: ", glb)
        print("Maxquery: ", max(query), "Maxsequence", max(sequence))
        print("Minquery: ", min(query), "Minsequence", min(sequence))
        print("Firstquery: ", query[1], "Firstsequence", sequence[1])
        print("GLB: ", abs(max(query) - max(sequence)) + abs(min(query) - min(sequence)) + abs(query[1] - sequence[1]))

        return 1


if __name__ == "__main__":
    main(4, 5)