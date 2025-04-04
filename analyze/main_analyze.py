import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from analyze import *

def main():
    # Example usage: specify the data folder and lists of beta and sigma values
    # Adjust these lists as needed
    folder = "../data/data_local"  # Folder where the observable files are stored
    L = 100
    betas = [0.1, 1.0, 10.0, 100.0]
    sigmas = [0.0, 0.1, 0.5]
    init="random"
    run=0
    # Create time series figures for each combination of beta and sigma
    for sigma in sigmas:
        for beta in betas:
            #pass
            analyze_time_series(folder, L, beta, sigma, init, run)

    # Create mean magnetization plot for each sigma
    analyze_mean_magnetization(folder, L, betas, sigmas, init, run)


if __name__ == '__main__':
    main()
