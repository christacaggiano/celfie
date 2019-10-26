#!/usr/bin/env python

# ## Imports


import numpy as np 
import pandas as pd 
import pickle as pkl
import bottleneck as bn
import sys 
import os 



# ## Functions

def generate_proportion_fixed(individuals, tissue, percentage, column):

    rows = []
    for i in range(individuals):
        vals = np.random.multinomial(100-percentage, np.ones(tissue-1)/(tissue-1), size=1)[0].tolist()
        vals.insert(column, percentage)
        t = np.array(vals)
        return(np.array(t/t.sum()))
#         rows.append([x/100 for x in vals])  # proportions must sum to 1

    return rows


def generate_alpha(tissue, people, precision=10):
    
    a = np.zeros((7,  tissue))
    cols = [0, 0.01, 0.1, 1,  5, 50, 100]

    for i in range(a.shape[0]): 
        
        a[i, :] = np.array(generate_proportion_fixed(1, tissue, cols[i], 0))

    return a 
    

    

def generate_gamma(i, j, t): 
    """
    randomly generates a matrix of methylation values between 0 and 1 
    :param int tissue: number of tissues 
    :param int cpg: number of cpgs
    :return: cpg x tissue matrix of random methylation values 
    """
    
    gamma = np.zeros((i, j))
    
    for n in range (j): 
        gamma[:, n] = np.random.uniform(0, 1, size=i)
    
    return gamma



def generate_beta(alpha, gamma, x_depths): 
    
    total_indiv = alpha.shape[0]
    i, total_cpg = gamma.shape
    
    beta = np.zeros((total_indiv, total_cpg))
    
    for n in range(total_indiv): 
        for j in range(total_cpg): 
           
            depth = x_depths[n, j]
            gamma_cpg = gamma[:, j]

            mix = np.random.choice(i, depth, replace=True, p=alpha[n, :])
            probability = gamma_cpg[mix]
            
            beta[n, j] = np.sum(np.random.binomial(1, probability, size=depth))
            
    return beta
            


def generate_counts(count, probability): 
    
    return np.random.binomial(count, probability, size=probability.shape)



def generate_depths(depth, shape): 
    
    return np.random.poisson(depth, shape)
    

# ## EM 

def expectation(gamma, alpha): 
    
    individuals, tissues = alpha.shape
    sites = gamma.shape[1]

    p0 = np.zeros((tissues, sites, individuals))
    p1 = np.zeros((tissues, sites, individuals))

    for n in range(individuals):
        
        for j in range(sites):
            p0_j = (1-gamma[:, j])*alpha[n, : ]
            p0_j = p0_j/(np.sum(p0_j))
            
            p1_j = (gamma[:, j])*alpha[n, : ]
            p1_j = p1_j/(np.sum(p1_j))
            
            p0[:, j, n] = p0_j 
            p1[:, j, n] = p1_j 
            
    return p0, p1

def add_pseduocounts(value, array, meth, meth_depths):
    """ finds values of gamma where logll cannot be computed, adds pseudo-counts to make 
    computation possible 

    value: checks for a value that will prevent computation; either 0 or 1
    array: gamma array to check for inproper value 
    meth: np array of methylation counts 
    meth_depths: np array of total number of reads (meth counts + unmethylated counts) """

    axis0, axis1 = np.where(array==value)  # find indices where value isn't able to be computed 

    for i, j in zip(axis0, axis1): 
        meth[i, j] += 1  # add one read to methylated counts 
        meth_depths[i, j] += 2  # adds two reads to total counts


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

    # pre-defines probability matrices  
    p0 = np.zeros((tissues, sites, individuals))
    p1 = np.zeros((tissues, sites, individuals))

    for n in range(individuals):
        for j in range(sites):
            p0_j = (1-gamma[:, j])*alpha[n, : ]
            p0_j = p0_j/(bn.nansum(p0_j))
            
            p1_j = (gamma[:, j])*alpha[n, : ]
            p1_j = p1_j/(bn.nansum(p1_j))
            
            p0[:, j, n] = p0_j 
            p1[:, j, n] = p1_j 
            
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

    tissues, sites, individuals = p0.shape[0], p0.shape[1], p0.shape[2]

    # for easier calculation, split the ll into two terms and calculate iteratively 
    term1 = 0
    term2 = 0 
    ll = 0 

    for n in range(individuals): 
        for j in range(sites): 
            for i in range(tissues):

                # calculates the q-function; see writeup 
                term1 += np.nan_to_num(y[i, j] + (p1[i, j, n]*x[n, j]))*np.log(gamma[i,j]) + (y_depths[i, j] - y[i, j] + p0[i, j, n]*(x_depths[n, j] - x[n, j]))*np.log(1-gamma[i,j])
                term2 += np.nan_to_num(((p1[i, j, n]*x[n, j]) + ((x_depths[n, j] - x[n, j])*p0[i, j, n])) * np.log(alpha[n, i]))
                
                ll += term1 + term2
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
    
        term1 +=  p1[:, :, n]*(np.outer(ones_vector, x[n,:]))
        term0 +=  p0[:, :, n]*(np.outer(ones_vector, x_depths[n,:]-x[n,:]))
    
    gamma = (term1 + y) / (term0 + term1 + y_depths)  # calculate new gamma 

    # check if gamma goes out of bounds, if so add psuedocounts to misbehaving y values
    if check_gamma(gamma): 
        add_pseduocounts(1, gamma, y, y_depths)
        add_pseduocounts(0, gamma, y, y_depths)
        gamma = (t1 + y) / (t0 + t1 + y_depths)  # recalculate gamma
    
    # return alpha to be normalized to sum to 1 
    return np.array([row/row.sum() for row in new_alpha]), gamma 
 

 ########################  run em  ########################
   
 
def em(x, x_depths, y, y_depths, num_iterations): 
    
    # randomly intialize alpha for each iteration 
    alpha = np.random.uniform(size=(x.shape[0], y.shape[0]))
    alpha = np.array([row/row.sum() for row in alpha])  # make alpha sum to 1 

    # begin by checking for instances where there are no counts for y or y_depths
    add_pseduocounts(1, np.nan_to_num(y/y_depths), y, y_depths)
    add_pseduocounts(0, np.nan_to_num(y/y_depths), y, y_depths)
    
    # intialize gamma to reference values 
    gamma = y/y_depths

    # perform EM for a given number of iterations 
    for i in range(num_iterations): 

        p0, p1 = expectation(gamma, alpha) 
        # ll = log_likelihood(p0, p1, x_depths, x, y_depths, y, gamma, alpha)
        a, g = maximization(p0, p1, x, x_depths, y, y_depths)
        
        # check convergence of alpha and gamma
        alpha_diff = np.mean(abs(a-alpha))/np.mean(abs(alpha))
        gamma_diff = np.mean(abs(g-gamma))/np.mean(abs(gamma))

        if alpha_diff + gamma_diff < 0.0001:  # if convergence criteria, break 
            break 

        else:  # set current evaluation of alpha and gamma 
            alpha = a 
            gamma = g 
    

    print(log_likelihood(p0, p1, x_depths, x, y_depths, y, gamma, alpha))  # print ll for random restarts
    
    return alpha, gamma


def generate_em_replicate(i, j, depth, gamma_depth, t): 
    
    alpha_int = generate_alpha(i, n)

    gamma_int = generate_gamma(i, j, t)
    
    Y_int_depths = generate_depths(gamma_depth, (i, j))
    Y_int = generate_counts(Y_int_depths, gamma_int)
    
    X_int_depths = generate_depths(depth, (n, j))

    X_int = generate_beta(alpha_int, gamma_int, X_int_depths) 

    # Y_int_depths[-1] = 0 
    # Y_int[-1] = 0  
   
    return Y_int_depths, Y_int, X_int_depths, X_int, alpha_int, gamma_int

if __name__=="__main__": 

    # np.random.seed(500)

    output_dir = sys.argv[1]
    rep_num = int(sys.argv[2])
    i = int(sys.argv[3])
    j = int(sys.argv[4])
    n = int(sys.argv[5])
    depth = int(sys.argv[6])
    gamma_depth = int(sys.argv[7])
    t = float(sys.argv[8])

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


    Y_depths, Y, X_depths, X, alpha_true, gamma_int = generate_em_replicate(i, j, depth, gamma_depth, t)

    alpha, gamma = em(X, X_depths, Y, Y_depths, 1000)

    with open(output_dir + "/" + str(rep_num) + "_alpha_est.pkl", "wb") as f:
        pkl.dump(alpha, f)

    with open(output_dir + "/" + str(rep_num)  + "_alpha_true.pkl", "wb") as f:
        pkl.dump(alpha_true, f)

    with open(output_dir + "/" + str(rep_num) + "_gamma_est.pkl", "wb") as f:
        pkl.dump(gamma, f)

    with open(output_dir + "/" + str(rep_num)  + "_gamma_true.pkl", "wb") as f:
        pkl.dump(gamma_int, f)




