#!/opt/homebrew/bin/python3
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from analyze import analyze_time_series, analyze_mean_magnetization


def main():
    # Example usage: specify the data folder and lists of beta and sigma values
    # Adjust these lists as needed
    # folder = "../data/data_local"
    if 1:
        folder = "../data/data_local/data_pool"
        L = 100
        Ts = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
        sigmas = [0.0]
        init = "random"
        run = 0
        # Create time series figures for each combination of beta and sigma
        for sigma in sigmas:
            for T in Ts:
                # analyze_time_series(folder, L, T, sigma, init, run)
                pass

        # Create mean magnetization plot for each sigma
        analyze_mean_magnetization(folder, L, Ts, sigmas, init, run)

    if 0:
        folder = "../data/20250404"
        L = 100
        betas = np.concatenate([np.arange(0.3, 1.01, 0.1), np.arange(5.0, 41.1, 1.0)])
        sigmas = np.arange(0.0, 0.51, 0.1)
        init = "random"
        run = 0
        # Create time series figures for each combination of beta and sigma
        for sigma in sigmas:
            for beta in betas:
                pass
                analyze_time_series(folder, L, beta, sigma, init, run)

        # Create mean magnetization plot for each sigma
        analyze_mean_magnetization(folder, L, betas, sigmas, init, run)


if __name__ == "__main__":
    main()
