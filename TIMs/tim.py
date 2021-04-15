import csv
import numpy as np
import bottleneck as bn  # faster nan calcs
import heapq
from collections import defaultdict
import sys


def distance(values):
    """
    calculates the absolute difference between the tissue methylation
    and the median for a CpG

    values: array of methylation proportions for all tissues

    return: distance between each tissue methylation prop and median, CpG
    median methylation
    """

    median = bn.nanmedian(values)
    dist = np.abs(values - median)
    return dist, median


def add_to_heap(heap, n, dist, median, cpg, percent, num):
    """
    keep track of the top CpGs using a max heap of size n

    heap: max heap keeping track of top tims for a specific tissues
    n: size of max heap
    dist: distance between tissue methylation and median methylation
    median: median methylation for CpG
    cpg: cpg positional information (tuple of chrom, start, end)
    percent: tissue methylation percent
    num: tissue number

    return: none
    """

    if len(heap) < n:  # only keep track of the n CpGs with greatest distance
        heapq.heappush(
            heap, (dist, cpg, median, percent, num)
        )  # note- works because it only checks value of dist, not the rest of the tuple
    else:
        heapq.heappushpop(heap, (dist, cpg, median, percent, num))


def get_cpgs(heap_list):
    """
    converts heap to a default dict in order to print

    heap_list: list of max heaps- size is the number of tissues

    return: dictionary of CpGs and TIM values
    """

    cpgs = defaultdict(list)

    for heap in heap_list:
        for value in heap:  # iterate over list of heaps

            # dictionary is CpG position info: number of tissue, distance between
            # tissue meth and median, methylation percent for tissue, median CpG meth
            cpgs[tuple(value[1])].append((value[4], value[0], value[3], value[2]))

    return cpgs


if __name__ == "__main__":

    # user parameters
    input_file = sys.argv[1]  # input bedfile of WGBS data
    output_file = sys.argv[2]  # path to output file
    num_values = int(sys.argv[3])  # number of values to keep as tims
    tissues = int(sys.argv[4])  # total number of tissues to calc tims for
    depth_filter = int(sys.argv[5])  # depth filter
    nan_filter = int(sys.argv[6])  # number of nans to allow in TIM calc

    # list of heaps needed to count max distances for each tissue
    distance_heaps = [[] for i in range(tissues)]

    with open(input_file, "r") as f:
        bed = csv.reader(f, delimiter="\t")

        for line in bed:
            cpg = line[0:3]  # positional information
            meth = np.asarray(line[3::2], dtype="float")
            depth = np.asarray(line[4::2], dtype="float")

            median_depth = bn.nanmedian(depth)

            np.seterr(
                divide="ignore", invalid="ignore"
            )  # ignore Nans when there are no counts
            percents = meth / depth

            nan_count = np.count_nonzero(np.isnan(percents))

            if (
                median_depth >= depth_filter and nan_count < nan_filter
            ):  # data must pass basic quality to be a tim
                dist, median = distance(percents)

                for i, col_dist in enumerate(dist):
                    if not np.isnan(col_dist):
                        add_to_heap(
                            distance_heaps[i],
                            num_values,
                            col_dist,
                            median,
                            cpg,
                            percents[i],
                            i,
                        )

    cpgs = get_cpgs(distance_heaps)  # get TIM heap for printing

    # write the output file as tab delimited file with columns being:
    # chrom, start, end, tissue # for tim, absolute difference, methylation prop for tissue, median methylation
    # for all other tissues
    with open(output_file, "w") as o:
        dist_out = csv.writer(o, delimiter="\t")

        # header
        dist_out.writerow(
            [
                "chrom",
                "start",
                "end",
                "tissue number",
                "difference",
                "tissue methylation",
                "other tissue methylation",
            ]
        )

        for c in cpgs:

            if len(cpgs[c]) == 1:  # exclude CpGs where it is a TIM for > 1 tissue

                position = [c[0], c[1], c[2]]
                entry = [cpgs[c][0][0], cpgs[c][0][1], cpgs[c][0][2], cpgs[c][0][3]]

                dist_out.writerow(position + entry)
