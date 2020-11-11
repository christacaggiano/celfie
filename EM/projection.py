
### Imports ###

import numpy as np 
import pandas as pd 
import random 
import pickle as pkl
import sys 
import os
from scipy.optimize import minimize 
from scipy.stats import binom


##### Functions #####

def generate_alpha(pkl_file):

    """
    reads in pickle file dictating tissue proportions 
    :param str pickle file of the tissue proportions (people X tissues)
    :return: values loaded as a numpy array 

    """

    t = pkl.load(open(pkl_file, "rb"))
    return t.T


def generate_gamma(i, j): 
    """
    randomly generates a matrix of methylation values between 0 and 1 
    :param int tissue: number of tissues 
    :param int cpg: number of cpgs
    :return: cpg x tissue matrix of random methylation values 
    """
    
    gamma = np.zeros((i, j))
    
    for n in range (j): 
        gamma[:, n] = np.random.uniform(0, 1, size=i)  # draws from a uniform distribution 
    
    return gamma


def generate_beta(alpha, gamma, x_depths): 
    """
    generates the cfDNA reads based on the generated tissue proportions, the true 
    methylation values in the reference and the depths 
    :param array alpha: tissue props 
    :param array gamma: methylation values for the reference 
    :param array x_depths: simulated read depths for each CpG in each individual in cfDNA input 
    :return: methylation reads for the cfDNA input (number of sites X number of individuals)

    """
    
    total_indiv = alpha.shape[0]
    i, total_cpg = gamma.shape
    
    beta = np.zeros((total_indiv, total_cpg))
    
    for n in range(total_indiv): 
        for j in range(total_cpg): 
           
            depth = x_depths[n, j]  # depth at a paticular cpg and person 
            gamma_cpg = gamma[:, j]  # "true" methylation value in the reference

            mix = np.random.choice(i, depth, replace=True, p=alpha[n, ])  # assign reads based on the tissue proportions for that individual 
            probability = gamma_cpg[mix]
            
            beta[n, j] = np.sum(np.random.binomial(1, probability, size=depth))  # the beta is the sum of all the individual reads coming from the tissues contributing to that cpg in that individual 
            
    return beta
            


def generate_counts(count, probability): 
    """
    generate the methylation read counts for the reference data 
    :param array count: read depths 
    :param array probability: probability of each cpg being methylated for each tissue 
    :return: array of methylated reads

    """
    
    return np.random.binomial(count, probability, size=probability.shape)



def generate_depths(depth, input_shape):
    """
    :param int depth: read depth
    :param tuple input_shape: number of tissues X number of CpGs
    :return: array of ints where each number represents the number of total reads in a tissue at a cpg
    """ 
   
    return np.random.poisson(depth, input_shape)


def generate_replicate(i, j, n, depth, gamma_depth, pkl_file, rep_num): 
    
    # method is only run for one individual at a time, select the person 
    # to run it for 
    alpha_int = generate_alpha(pkl_file)[rep_num].reshape(1, -1)  # reshape to be compatible dimensions
    gamma_int = generate_gamma(i, j)

    Y_int_depths = generate_depths(gamma_depth, (i, j))
    Y_int = generate_counts(Y_int_depths, gamma_int)
    
    X_int_depths = generate_depths(depth, (n, j))
    X_int = generate_beta(alpha_int, gamma_int, X_int_depths) 

    # add one unknown 
    Y_int_depths[0] = 0 
    Y_int[0] = 0 

    return Y_int_depths, Y_int, X_int_depths, X_int, alpha_int, gamma_int


def compute_projection(proportions, z=1):
    """
    performs a projection onto the L1 ball, as in Duchi et al 
    :param array proportions: estimate of a the proportions 
    :param int z: set to 1- indicates projection to probablistic simplex
    :return: projected estimate of the cell type proportions 
    """
    
    proportions_sorted = np.sort(proportions)[::-1]
    cumulative_sum = np.cumsum(proportions_sorted)

    for j, summed_val in enumerate(cumulative_sum): 
        value = proportions_sorted[j] - ((summed_val - z)/(j+1))

        if value > 0: 
            p = j

    theta = (cumulative_sum[p] - z)/(p+1)
    
    return (proportions - theta).clip(min=0)  


def counts_log_likelihood(alpha_est, X, X_depth, gamma):
    """
    calculates a binomial log likelihood 
    :param array alpha_est: estimate of the cell type proportions 
    :param array X: methylation counts for cfDNA input 
    :param array X_depth: total depths for cfDNA input
    :param array gamma: reference methylation proportions
    """

    alpha_est = compute_projection(alpha_est).flatten()  # compute the projection of the estimates
    b = np.matmul(alpha_est, gamma)  # the probability a cfDNA cpg comes from a reference tissue is the weighted average of the estimates of the tissues contributing to that person 
    ll = np.sum(binom.logpmf(X, X_depth, b, loc=0))  # log likelihood

    return -ll  # optimize negative ll 



def perform_optimization(alpha_est, gamma, X, X_depth):
    """
    use scipy optimize to perform minimization using BFGS
    :param individuals: number of individuals (always pretty much 1)
    :param tissues: number of tissues
    :param proportions_est:
    :return: truth, guess
    """

    bounds = tuple((0, 1) for x in range(np.shape(alpha_est)[1]))  # constrains that values in array must be prop

    # perform minimization using scipy optimization, BFGS technique. Max iterations==10,000
    alpha_est_guess = minimize(counts_log_likelihood, alpha_est, args=(X, X_depth, gamma),
                          bounds=bounds, method="L-BFGS-B", options={'maxiter': 1000, 'ftol': 1e-08})
    
    return alpha_est_guess["x"]


if __name__ == "__main__": 

    output_dir = sys.argv[1]  # output directory
    rep_num = int(sys.argv[2])  # individual to run least squares regression for
    i = int(sys.argv[3])  # number of tissues
    j = int(sys.argv[4])  # number of CpGs
    n = int(sys.argv[5])  # number of individuals 
    depth = int(sys.argv[6])  # depth in cfDNA
    gamma_depth = int(sys.argv[7])  # depth in reference
    pkl_file = str(sys.argv[8])  # tissue proportions to simulate from 


    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    Y_depths, Y, X_depths, X, alpha_true, gamma_int = generate_replicate(i, j, n, depth, gamma_depth, pkl_file)

    gamma = Y/Y_depths  # let the reference be the proportion of methylated to unmethlated reads

    # initial guess for the tissue estimates
    alpha_int = np.zeros((1, i)) + 0.5
    alpha_int = alpha_int/alpha_int.sum()

    alpha = perform_optimization(alpha_int, np.nan_to_num(gamma), X, X_depths)
    alpha = compute_projection(alpha)


    # print the true and estimate of tissue proportions
    with open(output_dir + "/" + str(rep_num) + "_alpha_est.pkl", "wb") as f:
        pkl.dump(alpha, f)

    with open(output_dir + "/" + str(rep_num)  + "_alpha_true.pkl", "wb") as f:
        pkl.dump(alpha_true, f)













