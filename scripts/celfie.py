#!/usr/bin/env python

import argparse
import os
import sys

import bottleneck as bn  # substantially speeds up calculations with nan's
import numpy as np
import pandas as pd
import pickle as pkl  # to save output

np.seterr(divide="ignore", invalid="ignore")

################  support functions   ################


def add_pseudocounts(value, array, meth, meth_depths):
    """finds values of gamma where logll cannot be computed, adds pseudo-counts to make
    computation possible

    value: checks for a value that will prevent computation; either 0 or 1
    array: gamma array to check for inproper value
    meth: np array of methylation counts
    meth_depths: np array of total number of reads (meth counts + unmethylated counts)
    """

    axis0, axis1 = np.where(
        array == value  # find indices where value isn't able to be computed
    )

    meth[axis0, axis1] += 1  # add one read to methylated counts
    meth_depths[axis0, axis1] += 2  # adds two reads to total counts


def check_gamma(array):
    """checks for values of gamma where log likelihood cannot be computed, returns
    true if can be computed

    array: np array to check
    """

    return (0 in array) or (1 in array)


########  expectation-maximization algorithm  ########


def expectation(gamma, alpha):
    """calculates the components needed for loglikelihood for each iteration of gamma and alpha

    gamma: np matrix of the estimated 'true' methylation proportions
    alpha: np matrix of estimated mixing proportions
    """

    alpha = alpha.T[:, np.newaxis, :]
    gamma = gamma[..., np.newaxis]

    p0 = (1.0 - gamma) * alpha
    p1 = gamma * alpha

    p0 /= np.nansum(p0, axis=0)[np.newaxis, ...]
    p1 /= np.nansum(p1, axis=0)[np.newaxis, ...]

    return p0, p1


def log_likelihood(p0, p1, x_depths, x, y_depths, y, gamma, alpha):
    """calculates the log likelihood P(X, Z, Y | alpha, gamma)

    p0: probability that read is methylated
    p1: probability read is unmethylated
    x_depths: input read depths
    x: input methylated reads
    y_depths: reference matrix read depths
    y: reference methylated counts
    gamma: estimated true methylation proportions
    alpha: estimated mixing proportions
    """

    tissues, sites, individuals = p0.shape[0], p0.shape[1], p0.shape[2]

    # Reshape arrays for faster computation
    alpha = alpha.T[:, np.newaxis, :]
    gamma = gamma[..., np.newaxis]

    y = y[..., np.newaxis]
    y_depths = y_depths[..., np.newaxis]

    x = x.T[np.newaxis, ...]
    x_depths = x_depths.T[np.newaxis, ...]

    ll = 0
    ll += np.sum((y + p1 * x) * np.log(gamma))
    ll += np.sum((y_depths - y + p0 * (x_depths - x)) * np.log(1.0 - gamma))
    ll += np.sum((p1 * x + (x_depths - x) * p0) * np.log(alpha))

    return ll


def maximization(p0, p1, x, x_depths, y, y_depths):

    """maximizes log-likelihood, calculated in the expectation step
    calculates new alpha and gamma given these new parameters

    p0: probability that read is methylated
    p1: probability read is unmethylated
    x_depths: input read depths
    x: input methylated reads
    y_depths: reference matrix read depths
    y: reference methylated counts
    """

    individuals = p0.shape[2]

    # initialize vector
    ones_vector = np.ones(shape=(y.shape[0]))
    new_alpha = np.zeros((x.shape[0], y.shape[0]))

    # in case of overflow or error, transform nans to 0 and inf to large float
    p0 = np.nan_to_num(p0)
    p1 = np.nan_to_num(p1)
    x = np.nan_to_num(x)
    x_depths = np.nan_to_num(x_depths)

    # break up calculation into two terms
    term0 = 0
    term1 = 0

    for n in range(individuals):

        new_alpha[n, :] = np.dot(p1[:, :, n], x[n, :]) + np.matmul(
            p0[:, :, n], (x_depths[n, :] - x[n, :])
        )

        term1 += p1[:, :, n] * (np.outer(ones_vector, x[n, :]))
        term0 += p0[:, :, n] * (np.outer(ones_vector, x_depths[n, :] - x[n, :]))

    gamma = (term1 + y) / (term0 + term1 + y_depths)  # calculate new gamma

    # check if gamma goes out of bounds, if so add psuedocounts to misbehaving y values
    if check_gamma(gamma):
        add_pseudocounts(1, gamma, y, y_depths)
        add_pseudocounts(0, gamma, y, y_depths)
        gamma = (term1 + y) / (term0 + term1 + y_depths)  # recalculate gamma

    # return alpha to be normalized to sum to 1
    normalized_new_alpha = new_alpha / np.sum(new_alpha, axis=1)[:, np.newaxis]
    return normalized_new_alpha, gamma


########################  run em  ########################


def em(x, x_depths, y, y_depths, num_iterations, convergence_criteria):
    """take in the input cfdna matrices and the reference data and
    runs the EM for the specified number of iterations, or stops once the
    convergence_criteria is reached

    x: methylated cfDNA read counts
    x_depths: depth of cfDNA
    y: methylated reference counts
    y_depths: depth of cfDNA
    convergence_criteria: difference between alpha + gamma before stopping

    """

    # randomly intialize alpha for each iteration
    alpha = np.random.uniform(size=(x.shape[0], y.shape[0]))
    alpha /= np.sum(alpha, axis=1)[:, np.newaxis]  # make alpha sum to 1

    # begin by checking for instances where there are no counts for y or y_depths
    add_pseudocounts(1, np.nan_to_num(y / y_depths), y, y_depths)
    add_pseudocounts(0, np.nan_to_num(y / y_depths), y, y_depths)

    # intialize gamma to reference values
    gamma = y / y_depths

    # perform EM for a given number of iterations
    for i in range(num_iterations):

        p0, p1 = expectation(gamma, alpha)
        a, g = maximization(p0, p1, x, x_depths, y, y_depths)

        # check convergence of alpha and gamma
        alpha_diff = np.mean(abs(a - alpha)) / np.mean(abs(alpha))
        gamma_diff = np.mean(abs(g - gamma)) / np.mean(abs(gamma))

        if (
            alpha_diff + gamma_diff < convergence_criteria
        ):  # if convergence criteria, break
            break

        else:  # set current evaluation of alpha and gamma
            alpha = a
            gamma = g

    ll = log_likelihood(
        p0, p1, x_depths, x, y_depths, y, gamma, alpha
    )  # print ll for random restarts

    return alpha, gamma, ll


################## read in data #######################


def define_arrays(sample, num_samples, num_unk):
    """
    takes input data matrix- cfDNA and reference, and creates the arrays to run in EM. Adds
    specified number of unknowns to estimate


    sample: pandas dataframe of data (samples and reference). Assumes there is 3 columns (chrom, start, end)
    before the samples and before the reference
    num_samples: number of samples to deconvolve
    num_unk: number of unknowns to estimate
    """

    test = sample.iloc[:, 3 : (num_samples * 2) + 3].values.T
    train = sample.iloc[:, (num_samples * 2) + 3 + 3 :].values.T

    x = test[::2, :]
    x_depths = test[1::2, :]

    y = train[::2, :]
    y_depths = train[1::2, :]

    # add one unknown component
    unknown = np.zeros((num_unk, y_depths.shape[1]))
    y_depths_unknown = np.append(y_depths, unknown, axis=0)
    y_unknown = np.append(y, unknown, axis=0)

    return (
        np.nan_to_num(x),
        np.nan_to_num(x_depths),
        np.nan_to_num(y_unknown),
        np.nan_to_num(y_depths_unknown),
    )


def parse_header_names(header):

    parsed_header = []

    for i in range(0, len(header), 2):
        parsed_header.append(header[i].split("_")[0])

    return parsed_header


def get_header(sample, num_samples, num_unk):
    """
    gets the tissue and sample names to be used in generating an interpretable output file

    sample: dataframe of input data- with header
    num_samples: number of cfDNA samples
    num_unk: number of unknowns to be estimated
    """

    header = list(sample)

    samples = parse_header_names(
        header[3 : (num_samples * 2) + 3]
    )  # samples are first part of header
    tissues = parse_header_names(
        header[(num_samples * 2) + 3 + 3 :]
    )  # tissues are second part of header

    unknowns = ["unknown" + str(i) for i in range(1, num_unk + 1)]

    return samples, tissues + unknowns


def write_output(output_file, output_matrix, header, index):
    """
    write estimated methylation proportions and tissue proportions as txt file

    output_file: outputfile name
    output_matrix: celfie estimate
    header: tissue names
    index: either number of cpgs or number of samples, depending on type of output
    written
    """

    output = pd.DataFrame(output_matrix)
    output.columns = header
    output.insert(
        0, "", index
    )  # insert either the sample names or cpg numbers as first col

    output.to_csv(output_file, sep="\t", index=False)  # save as text file


################## run #######################

if __name__ == "__main__":

    # read command line input parameters
    parser = argparse.ArgumentParser(
        description="CelFiE - Cell-free DNA decomposition. CelFie estimated the cell type of origin proportions of a cell-free DNA sample."
    )
    parser.add_argument("input_path", help="the path to the input file")
    parser.add_argument("output_directory", help="the path to the output directory")
    parser.add_argument("num_samples", type=int, help="Number of cfdna samples")
    parser.add_argument(
        "-m",
        "--max_iterations",
        default=1000,
        type=int,
        help="How long the EM should iterate before stopping, unless convergence criteria is met. Default 1000.",
    )
    parser.add_argument(
        "-u",
        "--unknowns",
        default=1,
        type=int,
        help="Number of unknown categories to be estimated along with the reference data. Default 1.",
    )
    parser.add_argument(
        "-p",
        "--parallel_job_id",
        default=1,
        type=int,
        help="Replicate number in a simulation experiment. Default 1. ",
    )
    parser.add_argument(
        "-c",
        "--convergence",
        default=0.0001,
        type=float,
        help="Convergence criteria for EM. Default 0.0001.",
    )
    parser.add_argument(
        "-r",
        "--random_restarts",
        default=10,
        type=int,
        help="CelFiE will perform several random restarts and select the one with the highest log-likelihood. Default 10.",
    )
    args = parser.parse_args()

    # make output directory if it does not exist
    if not os.path.exists(args.output_directory) and args.parallel_job_id == 1:
        os.makedirs(args.output_directory)
        print("made " + args.output_directory + "/")
        print()
    else:
        print("writing to " + args.output_directory + "/")

    data_df = pd.read_csv(
        args.input_path, delimiter="\t"
    )  # read input samples/reference data

    print(f"finished reading {args.input_path}")
    print()

    output_alpha_file = f"{args.output_directory}/{args.parallel_job_id}_tissue_proportions.txt"
    output_gamma_file = f"{args.output_directory}/{args.parallel_job_id}_methylation_proportions.txt"

    print(f"beginning generation of {args.output_directory}")
    print()

    # make input arrays and add the specified number of unknowns
    x, x_depths, y, y_depths = define_arrays(data_df, int(args.num_samples), int(args.unknowns))

    # get header for output files
    samples, tissues = get_header(data_df, args.num_samples, args.unknowns)

    # Run EM with the specified iterations and convergence criteria
    random_restarts = []

    for i in range(args.random_restarts):
        alpha, gamma, ll = em(
            x, x_depths, y, y_depths, args.max_iterations, args.convergence
        )
        random_restarts.append((ll, alpha, gamma))

    ll_max, alpha_max, gamma_max = max(
        random_restarts
    )  # pick best random restart per replicate

    # write estimates as text files
    write_output(output_alpha_file, alpha_max, tissues, samples)
    write_output(
        output_gamma_file, gamma_max.T, tissues, list(range(len(gamma_max[1])))
    )
