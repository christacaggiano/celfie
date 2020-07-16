#!/usr/bin/env python

########  imports  ########
import numpy as np 
import pandas as pd 
import pickle as pkl # to save output 
import bottleneck as bn # substantially speeds up calculations with nan's
import os 
import sys 
np.seterr(divide='ignore', invalid='ignore')

################  support functions   ################


def add_pseudocounts(value, array, meth, meth_depths):
    """ finds values of gamma where logll cannot be computed, adds pseudo-counts to make 
    computation possible 

    value: checks for a value that will prevent computation; either 0 or 1
    array: gamma array to check for inproper value 
    meth: np array of methylation counts 
    meth_depths: np array of total number of reads (meth counts + unmethylated counts) """

    # find indices where value isn't able to be computed
    axis0, axis1 = np.where(array == value)

    # add one read to methylated counts 
    meth[axis0, axis1] += 1

    # add two reads to total counts
    meth_depths[axis0, axis1] += 2


def check_gamma(array): 
    """ checks for values of gamma where log likelihood cannot be computed, returns 
    true if can be computed 
    
    array: np array to check """ 

    return (0 in array) or (1 in array) 


########  expectation-maximization algorithm  ########


def expectation(gamma, alpha):
    """ calculates the components needed for loglikelihood for each iteration of gamma and alpha 

    gamma: np matrix of the estimated 'true' methylation proportions 
    alpha: np matrix of estimated mixing proportions """
    

    individuals, tissues = alpha.shape
    sites = gamma.shape[1]

    alpha = alpha.T[:, np.newaxis, :]
    gamma = gamma[..., np.newaxis]

    p0 = (1. - gamma) * alpha
    p1 = gamma * alpha

    p0 /= np.nansum(p0, axis=0)[np.newaxis, ...]
    p1 /= np.nansum(p1, axis=0)[np.newaxis, ...]
            
    return p0, p1


def log_likelihood(p0, p1, x_depths, x, y_depths, y, gamma, alpha):
    """ calculates the log likelihood P(X, Z, Y | alpha, gamma)

    p0: probability that read is methylated 
    p1: probability read is unmethylated 
    x_depths: input read depths
    x: input methylated reads
    y_depths: reference matrix read depths
    y: reference methylated counts
    gamma: estimated true methylation proportions 
    alpha: estimated mixing proportions """

    # Reshape arrays for faster computation
    alpha = alpha.T[:, np.newaxis, :]
    gamma = gamma[..., np.newaxis]
    y = y[..., np.newaxis]
    y_depths = y_depths[..., np.newaxis]
    x = x.T[np.newaxis, ...]
    x_depths = x_depths.T[np.newaxis, ...]

    ll = 0
    ll += np.sum((y + p1 * x) * np.log(gamma))
    ll += np.sum((y_depths - y + p0 * (x_depths - x)) * np.log(1. - gamma))
    ll += np.sum((p1 * x + (x_depths - x) * p0) * np.log(alpha))

    return ll 


def maximization(p0, p1, x, x_depths, y, y_depths): 

    """ maximizes log-likelihood, calculated in the expectation step
    calculates new alpha and gamma given these new parameters 

    p0: probability that read is methylated 
    p1: probability read is unmethylated 
    x_depths: input read depths
    x: input methylated reads
    y_depths: reference matrix read depths
    y: reference methylated counts """


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
       
        new_alpha[n, :] = np.dot(p1[:, :, n], x[n,:]) + np.matmul(p0[:, :, n], (x_depths[n, :]-x[n, :]))
    
        term1 +=  p1[:, :, n] * (np.outer(ones_vector, x[n, :]))
        term0 +=  p0[:, :, n] * (np.outer(ones_vector, x_depths[n, :]-x[n, :]))
    
    gamma = (term1 + y) / (term0 + term1 + y_depths)  # calculate new gamma 

    # check if gamma goes out of bounds, if so add psuedocounts to misbehaving y values
    if check_gamma(gamma): 
        add_pseudocounts(1, gamma, y, y_depths)
        add_pseudocounts(0, gamma, y, y_depths)
        gamma = (term1 + y) / (term0 + term1 + y_depths)   # recalculate gamma
    
    # return alpha to be normalized to sum to 1
    normalized_new_alpha = new_alpha / np.sum(new_alpha, axis=1)[:, np.newaxis]
    return normalized_new_alpha, gamma 
 

 ########################  run em  ########################
   
 
def em(x, x_depths, y, y_depths, num_iterations, convergence_criteriav): 
    
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

        if alpha_diff + gamma_diff < convergence_criteriav:  # if convergence criteria, break 
            break 

        else:  # set current evaluation of alpha and gamma 
            alpha = a 
            gamma = g 
    

    ll = log_likelihood(p0, p1, x_depths, x, y_depths, y, gamma, alpha)  # print ll for random restarts
    
    return alpha, gamma, ll

################## read in data #######################

def define_arrays(sample, num_samples, num_unk):
    """
    takes input data matrix- cfDNA and reference, and creates the arrays to run in EM. Adds 1 unknown
    sample: pandas dataframe of data (samples and reference). Assumes there is 3 columns (chrom, start, end)
    before the samples and before the reference 
    """

    # Split dataframe in two parts
    split = (num_samples * 2) + 3
    test_df = sample.iloc[:, :split]
    train_df = sample.iloc[:, split:]

    # Get genomic coordinates of test samples
    get_idx = lambda df: df.iloc[:, 0] + '-' + df.iloc[:, 1].astype(str) + '-' + df.iloc[:, 2].astype(str)
    test_idx = get_idx(test_df)
    test_df = test_df.set_index(test_idx, drop=False)

    # Get genomic coordinates of reference samples
    train_idx = get_idx(train_df)
    train_df = train_df.set_index(train_idx, drop=False)

    # Re-order reference samples by coordinates
    train_df = train_df.loc[test_idx, :]

    # Convert dataframes to arrays
    test = test_df.iloc[:, 3:].values.T
    train = train_df.iloc[:, 3:].values.T

    x = test[::2, :]
    x_depths = test[1::2, :]

    y = train[::2, :]
    y_depths = train[1::2, :]

    # add one unknown component 
    unknown = np.zeros((num_unk, y_depths.shape[1]))
    y_depths_unknown = np.append(y_depths, unknown, axis=0)
    y_unknown = np.append(y, unknown, axis=0)

    return np.nan_to_num(x), np.nan_to_num(x_depths), np.nan_to_num(y_unknown), np.nan_to_num(y_depths_unknown)

################## run #######################



if __name__=="__main__": 
    
    # read command line input parameters 
    data = sys.argv[1]
    output_dir = sys.argv[2]
    num_samples = sys.argv[3]
    iterations = sys.argv[4]
    num_unk = sys.argv[5]
    iteration_number = sys.argv[6]
    convergence_criteria = sys.argv[7]
    num_random_restart = sys.argv[8]

    # make output directory if it does not exist
    if not os.path.exists(output_dir) and int(iteration_number)==1:
        os.makedirs(output_dir)
        print("made " + output_dir + "/")
        print()
    else: 
        print("writing to " + output_dir + "/")

    data_df = pd.read_csv(data, header=None, delimiter="\t")  # read input samples/reference data
    
    print("finished reading " + str(data))
    print()

    output_alpha_file = output_dir + "/" + iteration_number + "_alpha.pkl"
    output_gamma_file = output_dir + "/" + iteration_number + "_gamma.pkl"

    print("beginning generation of " + output_alpha_file)
    print()

    # make input arrays and add the specified number of unknowns 
    x, x_depths, y, y_depths = define_arrays(data_df, int(num_samples), int(num_unk))

    # Run EM with the specified iterations and convergence criteria
    random_restarts = []

    for i in range(int(num_random_restart)): 
        alpha, gamma, ll = em(x, x_depths, y, y_depths, int(iterations), float(convergence_criteria))
        random_restarts.append((ll, alpha, gamma))

    ll_max, alpha_max, gamma_max = max(random_restarts)  # pick best random restart per replicate

    # write estimates as pickle files 
    with open(output_dir + "/" + str(iteration_number) + "_alpha.pkl", "wb") as f:
        pkl.dump(alpha_max, f)

    with open(output_dir + "/" + str(iteration_number) + "_gamma.pkl", "wb") as f:
        pkl.dump(gamma_max, f)

