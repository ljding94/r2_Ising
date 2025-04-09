import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def analyze_time_series(folder, L, T, sigma, method, run):
    """
    Analyze the observable time series for a given combination of T and sigma.
    The expected data file name is assumed to be: obs_T_<T>_sigma_<sigma>.txt
    with columns: m, m2, E
    """
    finfo = f"L{L:.0f}_T{T:.1f}_sigma{sigma:.1f}_{method}_run{run}"
    # Construct filename based on T and sigma
    filename = f"{folder}/obs_{finfo}.txt"
    if not os.path.exists(filename):
        print(f"File {filename} does not exist.")
        return

    # Load data using pandas
    df = pd.read_csv(filename)
    # Assume the file has a header: m,m2,E
    steps = np.arange(len(df))

    # Create a figure with three subplots for m, m2, and E vs MC step
    fig, axs = plt.subplots(2, 1, figsize=(5, 6), sharex=True)

    axs[0].plot(steps, df["m"], "b-")
    axs[0].set_ylabel("m")
    axs[0].set_title(f"T = {T}, sigma = {sigma} : m vs MC step")

    axs[1].plot(steps, df["E"], "g-")
    axs[1].set_ylabel("E")
    axs[1].set_xlabel("MC step")
    axs[1].set_title(f"T = {T}, sigma = {sigma} : E vs MC step")

    plt.tight_layout()
    # Save the figure
    plt.savefig(f"{folder}/time_series_{finfo}.png")
    plt.close()


def analyze_mean_magnetization(folder, L, Ts, sigmas, method, run):
    """
    For multiple sigma values, analyze the mean magnetization for each T and plot |<m>| vs T.
    Assumes data file for each T is named: obs_L{L}_T{T}_sigma{sigma}_init{init}.txt
    """
    plt.figure(figsize=(4, 8))
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)
    print("Ts",Ts)
    for sigma in sigmas:
        T_vals = []
        mean_abs_m = []
        std_abs_m = []  # To store the standard deviation of |m|
        chi_vals = []
        mean_E = []
        std_E = []

        for T in Ts:
            print("Analyzing T =", T, "for sigma =", sigma)  # Debugging statement to track the current T and sigma
            finfo = f"L{L:.0f}_T{T:.1f}_sigma{sigma:.1f}_{method}_run{run}"
            filename = f"{folder}/obs_{finfo}.txt"
            print(f"Analyzing file: {filename}")  # Debugging statement to check the filename being processed
            if not os.path.exists(filename):
                print(f"File {filename} does not exist. Skipping T = {T} for sigma = {sigma}")
                continue
            df = pd.read_csv(filename)
            print(df)
            # Compute the average magnetization over MC steps
            # Skip first half of data for thermalization
            half_index = len(df) // 2
            abs_m = np.abs(df["m"].iloc[half_index:])  # Use the second half of the data for averaging
            abs_m_mean = abs_m.mean()  # This is <|m|>
            abs_m_std = abs_m.std()  # Standard deviation of |m|, can be used for error bars if needed
            m2 = np.power(df["m"], 2)
            m2_mean = m2.iloc[half_index:].mean()
            chi = (m2_mean - abs_m_mean**2) / T
            mean_E.append(df["E"].iloc[half_index:].mean())
            std_E.append(df["E"].iloc[half_index:].std())
            T_vals.append(T)
            mean_abs_m.append(abs_m_mean)
            std_abs_m.append(abs_m_std)  # Store the standard deviation of |m| for potential error bars
            chi_vals.append(chi)
        N = half_index

        # ax1.plot(1.0/np.array(T_vals), mfan_m2_vals, 'o-', label=f'{sigma:.1f}')
        print("np.array(mean_abs_m)**2", np.array(mean_abs_m) ** 2)
        print("chi_vals", chi_vals)
        mean_abs_m = np.array(mean_abs_m)
        std_abs_m = np.array(std_abs_m)
        ax1.errorbar(T_vals, mean_abs_m ** 2, yerr=2*mean_abs_m*std_abs_m/np.sqrt(N/100), ls="-", marker="o", ms="3", label=f"{sigma:.1f}")
        ax2.plot(T_vals, chi_vals, "o--", ms="3", label=f"{sigma:.1f} ")
        ax3.errorbar(T_vals, mean_E, yerr=std_E, ls="-", marker="o", ms="3", label=f"{sigma:.1f}")

    ax1.set_ylim(0, 1.1)
    ax1.set_xlabel(r"$T$")
    ax1.set_ylabel(r"$<|m|>^2$")
    ax1.legend(title=r"$\sigma$", loc="upper right")
    ax2.set_xlabel(r"$T$")
    ax2.set_ylabel(r"$\chi = (<m^2> - <m>^2)/T$")
    ax2.legend(title=r"$\sigma$", loc="upper right")
    ax3.set_xlabel(r"$T$")
    ax3.set_ylabel(r"$<E>$")
    ax3.legend(title=r"$\sigma$", loc="upper right")

    plt.tight_layout()
    output_fig = os.path.join(folder, "mean_magnetization_all.png")
    plt.savefig(output_fig)
    plt.show()  # Show the plot for interactive environments, if needed
    plt.close()
    print(f"Saved mean magnetization figure to {output_fig}")



def analyze_parallel_mean_magnetization(folder, L, Ti, Tf, nT, sigmas, method, run):
    plt.figure(figsize=(4, 8))
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)
    for sigma in sigmas:
        T_vals = np.geomspace(Ti, Tf, nT)
        mean_abs_m = []
        std_abs_m = []  # To store the standard deviation of |m|
        chi_vals = []
        mean_E = []
        std_E = []
        print("Analyzing sigma =", sigma)  # Debugging statement to track the current T and sigma
        finfo = f"L{L:.0f}_Ti{Ti:.2f}_Tf{Tf:.2f}_nT{nT:.0f}_sigma{sigma:.1f}_method_{method}_run_{run}"
        filename = f"{folder}/parallel_obs_{finfo}.txt"
        print(f"Analyzing file: {filename}")  # Debugging statement to check the filename being processed
        if not os.path.exists(filename):
            print(f"File {filename} does not exist.  for sigma = {sigma}")
            continue
        df = pd.read_csv(filename, skiprows=1)
        half_index = len(df) // 2
        for i in range(nT):
            abs_m = np.abs(df[f"m_{i}"].iloc[half_index:])  # Use the second half of the data for averaging
            abs_m_mean = abs_m.mean()  # This is <|m|>
            abs_m_std = abs_m.std()  # Standard deviation of |m|, can be used for error bars if needed
            m2 = np.power(df[f"m_{i}"], 2)
            m2_mean = m2.iloc[half_index:].mean()
            chi = (m2_mean - abs_m_mean**2) / T_vals[i]
            mean_E.append(df[f"m_{i}"].iloc[half_index:].mean())
            std_E.append(df[f"m_{i}"].iloc[half_index:].std())
            mean_abs_m.append(abs_m_mean)
            std_abs_m.append(abs_m_std)  # Store the standard deviation of |m| for potential error bars
            chi_vals.append(chi)
        N = half_index

        # ax1.plot(1.0/np.array(T_vals), mfan_m2_vals, 'o-', label=f'{sigma:.1f}')
        print("np.array(mean_abs_m)**2", np.array(mean_abs_m) ** 2)
        print("chi_vals", chi_vals)
        mean_abs_m = np.array(mean_abs_m)
        std_abs_m = np.array(std_abs_m)
        ax1.errorbar(T_vals, mean_abs_m ** 2, yerr=2*mean_abs_m*std_abs_m/np.sqrt(N/100), ls="-", marker="o", ms="3", label=f"{sigma:.1f}")
        ax2.plot(T_vals, chi_vals, "o--", ms="3", label=f"{sigma:.1f} ")
        ax3.errorbar(T_vals, mean_E, yerr=std_E, ls="-", marker="o", ms="3", label=f"{sigma:.1f}")

    ax1.set_ylim(0, 1.1)
    ax1.set_xlabel(r"$T$")
    ax1.set_ylabel(r"$<|m|>^2$")
    ax1.legend(title=r"$\sigma$", loc="upper right", ncol=2, handlelength=0.5)
    ax2.set_xlabel(r"$T$")
    ax2.set_ylabel(r"$\chi = (<m^2> - <m>^2)/T$")
    ax2.legend(title=r"$\sigma$", loc="upper right", ncol=2, handlelength=0.5)
    ax3.set_xlabel(r"$T$")
    ax3.set_ylabel(r"$<E>$")
    ax3.legend(title=r"$\sigma$", loc="upper right", ncol=2, handlelength=0.5)

    plt.tight_layout()
    output_fig = os.path.join(folder, "mean_magnetization_all.png")
    plt.savefig(output_fig)
    plt.show()  # Show the plot for interactive environments, if needed
    plt.close()
    print(f"Saved mean magnetization figure to {output_fig}")
