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

    if 1:
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


if __name__ == "__main__":
    main()
