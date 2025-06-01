#!/opt/homebrew/bin/python3
import numpy as np
from analyze import *


def main():
    # Example usage: specify the data folder and lists of beta and sigma values
    # Adjust these lists as needed
    # folder = "../data/data_local"
    if 0:
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
        folder = "../data/20250407"
        L = 100
        Ts = np.arange(0.1, 2.51, 0.1)
        print("Ts", Ts)
        sigmas = np.arange(0.0, 0.91, 0.1)
        method = "single"
        run = 0
        # Create time series figures for each combination of beta and sigma
        for sigma in sigmas:
            for T in Ts:
                pass
                # analyze_time_series(folder, L, T, sigma, init, run)

        # Create mean magnetization plot for each sigma
        print("Ts1", Ts)
        analyze_mean_magnetization(folder, L, Ts, sigmas, method, run)

    if 0:

        folder = "../data/20250409"
        L = 100
        Ti = 0.05
        Tf = 2.0
        nT = 20
        sigmas = np.arange(0.05, 1.01, 0.05)
        method = "single"
        run = 0
        # Create mean magnetization plot for each sigma
        analyze_parallel_mean_magnetization(folder, L, Ti, Tf, nT, sigmas, method, run)

    if 0:
        folder = "../data/20250429"
        L = 100
        Ti = 0.02
        Tf = 2.50
        nT = 100
        sigmas = np.arange(0.10, 1.01, 0.05)
        method = "single"
        doswap = 0
        Mrun = 50
        analyze_parallel_mean_magnetization_multirun_per_sigma(folder, L, Ti, Tf, nT, sigmas, method, doswap, Mrun)

    if 1:
        #folder = "../data/20250430"
        #folder = "../data/20250519"
        folder = "../data/20250522"
        Ls = np.arange(300, 381, 20)
        #Ti = 0.02
        #Tf = 0.50
        #nT = 20
        Ti = 1.60
        Tf = 2.20
        nT = 20

        sigma = 0.20
        method = "single"
        doswap = 0
        Mrun = 50
        analyze_parallel_mean_magnetization_multirun_per_L(folder, Ls, Ti, Tf, nT, sigma, method, doswap, Mrun)

    if 0:
        folder = "../data/20250430"
        L = 100
        Ti = 0.02
        Tf = 0.50
        nT = 20
        sigma = 0.20
        method = "single"
        doswap = 0
        for run in range(20):
            analyze_parallel_time_series(folder, L, Ti, Tf, nT, sigma, method, doswap, run)

    if 0:
        # histogram of E and m2
        folder = "../data/20250430"
        Ls = np.arange(60, 401, 20)
        #Ls = [60, 120, 240, 360, 400]
        Ti = 0.02
        Tf = 0.50
        nT = 20
        sigma = 0.20
        method = "single"
        doswap = 0
        Mrun = 50
        #analyze_parallel_histgram(folder, Ls, Ti, Tf, nT, sigma, method, doswap, Mrun)
        analyze_parallel_histgram_perT(folder, Ls, Ti, Tf, nT, sigma, method, doswap, Mrun)



#TODO: add a plot for the m vs. steps

if __name__ == "__main__":
    main()
