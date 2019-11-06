##### Imports ########
import numpy as np 
import pandas as pd 
from random import randint
import pickle as pkl
import bottleneck 
import sys 
import time 
   

def preg(sample, sample_size):
    """
    takes input data matrix- cfDNA and reference, and creates the arrays to run in EM 
    sample: pandas dataframe of data 
    sample_size: if the data is being sampled, include what percentage is to be samples 
    """

    test = sample.iloc[:, :30].values.T
    train = sample.iloc[:, 33:].values.T

    x = test[::2, :]
    x_depths = test[1::2, :]

    y = train[::2, :]
    y_depths = train[1::2, :]

    # add one unknown component 
    unknown = np.zeros((1, y_depths.shape[1]))
    y_depths_unknown = np.append(y_depths, unknown, axis=0)
    y_unknown = np.append(y, unknown, axis=0)

    return np.nan_to_num(x), np.nan_to_num(x_depths), np.nan_to_num(y_unknown), np.nan_to_num(y_depths_unknown)

