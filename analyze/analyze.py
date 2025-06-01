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
    print("Ts", Ts)
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
        ax1.errorbar(T_vals, mean_abs_m**2, yerr=2 * mean_abs_m * std_abs_m / np.sqrt(N / 100), ls="-", marker="o", ms="3", label=f"{sigma:.1f}")
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
    plt.figure(figsize=(6, 12))
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)
    colormap = plt.get_cmap("rainbow")  # you can change this to another colormap if desired
    colors = colormap(np.linspace(0, 1, len(sigmas)))

    for idx, sigma in enumerate(sigmas):
        color = colors[idx]
        T_vals = np.geomspace(Ti, Tf, nT)
        mean_abs_m = []
        std_abs_m = []  # To store the standard deviation of |m|
        chi_vals = []
        mean_E = []
        std_E = []
        print("Analyzing sigma =", sigma)  # Debugging statement to track the current T and sigma
        finfo = f"L{L:.0f}_Ti{Ti:.2f}_Tf{Tf:.2f}_nT{nT:.0f}_sigma{sigma:.2f}_method_{method}_run_{run}"
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
            mean_E.append(df[f"E_{i}"].iloc[half_index:].mean())
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
        ax1.errorbar(T_vals, mean_abs_m**2, yerr=2 * mean_abs_m * std_abs_m / np.sqrt(N), ls="-", marker="o", ms="3", color=color, label=f"{sigma:.2f}")
        ax2.plot(T_vals, chi_vals, "o--", ms="3", color=color, label=f"{sigma:.2f} ")
        ax3.errorbar(T_vals, mean_E, yerr=std_E / np.sqrt(N), ls="-", marker="o", ms="3", color=color, label=f"{sigma:.2f}")

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


def get_multi_run_data(folder, L, Ti, Tf, nT, sigma, method, doswap, Mrun):
    # return mean_m2, std_m2, mean_E, std_E, mean_chi, std_chi versus all T

    # T_vals = np.geomspace(Ti, Tf, nT)
    T_vals = np.linspace(Ti, Tf, nT)
    all_abs_m2_mean = []
    all_E_mean = []
    all_chi_mean = []
    for run in range(Mrun):
        finfo = f"L{L:.0f}_Ti{Ti:.2f}_Tf{Tf:.2f}_nT{nT:.0f}_sigma{sigma:.2f}_method_{method}_doswap{doswap}_run_{run}"
        filename = f"{folder}/parallel_obs_{finfo}.txt"
        print(f"Analyzing file: {filename}")

        if not os.path.exists(filename):
            print(f"File {filename} does not exist.  for sigma = {sigma}")
            continue
        df = pd.read_csv(filename, skiprows=1)
        half_index = len(df) // 2
        mean_abs_m_sq_per_T = []
        mean_E_per_T = []
        mean_chi_per_T = []
        for i in range(nT):
            abs_m = np.abs(df[f"m_{i}"].iloc[half_index:])
            abs_m_mean = abs_m.mean()  # This is <|m|>
            m2 = np.power(df[f"m_{i}"], 2)
            m2_mean = m2.iloc[half_index:].mean()
            chi = (m2_mean - abs_m_mean**2) / T_vals[i]
            E_mean = df[f"E_{i}"].iloc[half_index:].mean()
            mean_abs_m_sq_per_T.append(abs_m_mean**2)
            mean_E_per_T.append(E_mean)
            mean_chi_per_T.append(chi)
        all_abs_m2_mean.append(mean_abs_m_sq_per_T)
        all_E_mean.append(mean_E_per_T)
        all_chi_mean.append(mean_chi_per_T)

    return T_vals, all_abs_m2_mean, all_E_mean, all_chi_mean

def analyze_parallel_time_series(folder, L, Ti, Tf, nT, sigma, method, doswap, run):
    finfo = f"L{L:.0f}_Ti{Ti:.2f}_Tf{Tf:.2f}_nT{nT:.0f}_sigma{sigma:.2f}_method_{method}_doswap{doswap}_run_{run}"
    filename = f"{folder}/parallel_obs_{finfo}.txt"
    if not os.path.exists(filename):
        print(f"File {filename} does not exist.")
        return
    df = pd.read_csv(filename, skiprows=1)

    print("df.columns", df.columns)
    print("df[m_0]", df["m_0"])
    # Assume the file has a header: m,m2,E
    steps = np.arange(len(df))

    # Create a figure with three subplots for m, m2, and E vs MC step
    fig, axs = plt.subplots(nT, 2, figsize=(5, 20), sharex=True)

    T_vals = np.linspace(Ti, Tf, nT)
    colormap = plt.get_cmap("rainbow")  # you can change this to another colormap if desired
    colors = colormap(np.linspace(0, 1, nT))
    for i in range(nT):
        axs[i, 0].plot(steps, df[f"m_{i}"], "-", lw=1, label=f"{T_vals[i]:.2f}", color=colors[i])
        axs[i, 1].plot(steps, df[f"E_{i}"], "-", lw=1, label=f"{T_vals[i]:.2f}", color=colors[i])
        axs[i, 0].legend(title=r"$T$")

    axs[nT - 1, 0].set_ylabel("m")
    axs[nT - 1, 0].set_xlabel("MC step")

    axs[nT - 1, 1].set_ylabel("E")
    axs[nT - 1, 1].set_xlabel("MC step")

    plt.tight_layout(pad=0.1)
    # Save the figure
    plt.savefig(f"{folder}/time_series_{finfo}.png")
    plt.close()


def analyze_parallel_mean_magnetization_multirun_per_sigma(folder, L, Ti, Tf, nT, sigmas, method, doswap, Mrun):
    plt.figure(figsize=(6, 12))
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)
    colormap = plt.get_cmap("rainbow")  # you can change this to another colormap if desired
    colors = colormap(np.linspace(0, 1, len(sigmas)))

    mean_all_abs_m2_mean = []
    std_all_abs_m2_mean = []
    mean_all_E_mean = []
    std_all_E_mean = []
    mean_all_chi_mean = []
    std_all_chi_mean = []
    for idx, sigma in enumerate(sigmas):
        color = colors[idx]
        T_vals, all_abs_m2_mean, all_E_mean, all_chi_mean = get_multi_run_data(folder, L, Ti, Tf, nT, sigma, method, doswap, Mrun)

        mean_all_abs_m2_mean.append(np.mean(all_abs_m2_mean, axis=0))
        std_all_abs_m2_mean.append(np.std(all_abs_m2_mean, axis=0))
        mean_all_E_mean.append(np.mean(all_E_mean, axis=0))
        std_all_E_mean.append(np.std(all_E_mean, axis=0))
        mean_all_chi_mean.append(np.mean(all_chi_mean, axis=0))
        std_all_chi_mean.append(np.std(all_chi_mean, axis=0))

        ax1.errorbar(T_vals, np.mean(all_abs_m2_mean, axis=0), yerr=np.std(all_abs_m2_mean, axis=0) / np.sqrt(Mrun), ls="-", marker="o", ms="3", color=color, label=f"{sigma:.2f}")
        ax2.errorbar(T_vals, np.mean(all_chi_mean, axis=0), yerr=np.std(all_chi_mean, axis=0) / np.sqrt(Mrun), ls="--", marker="o", ms="3", color=color, label=f"{sigma:.2f}")
        ax3.errorbar(T_vals, np.mean(all_E_mean, axis=0), yerr=np.std(all_E_mean, axis=0) / np.sqrt(Mrun), ls="-", marker="o", ms="3", color=color, label=f"{sigma:.2f}")

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

    # save the Tval, mean_abs_m2_mean, std_abs_m2_mean, mean_E_mean, std_E_mean, mean_chi_mean, std_chi_mean for all sigma to a csv file
    outputfile = os.path.join(folder, "mean_magnetization_multirun.csv")
    # Create a DataFrame to save the results
    df = pd.DataFrame()
    # First column is temperature values
    df["T"] = T_vals
    # Add data for each sigma
    for i, sigma in enumerate(sigmas):
        df[f"sigma_{sigma:.2f}_m2_mean"] = mean_all_abs_m2_mean[i]
        df[f"sigma_{sigma:.2f}_m2_std"] = std_all_abs_m2_mean[i]
        df[f"sigma_{sigma:.2f}_E_mean"] = mean_all_E_mean[i]
        df[f"sigma_{sigma:.2f}_E_std"] = std_all_E_mean[i]
        df[f"sigma_{sigma:.2f}_chi_mean"] = mean_all_chi_mean[i]
        df[f"sigma_{sigma:.2f}_chi_std"] = std_all_chi_mean[i]
    # Save to CSV
    df.to_csv(outputfile, index=False)
    print(f"Saved data to {outputfile}")

    plt.tight_layout()
    output_fig = os.path.join(folder, "mean_magnetization_all_multirun.png")
    plt.savefig(output_fig)
    plt.show()  # Show the plot for interactive environments, if needed
    plt.close()
    print(f"Saved mean magnetization figure to {output_fig}")


def analyze_parallel_mean_magnetization_multirun_per_L(folder, Ls, Ti, Tf, nT, sigma, method, doswap, Mrun):
    plt.figure(figsize=(6, 12))
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)
    colormap = plt.get_cmap("rainbow")  # you can change this to another colormap if desired
    colors = colormap(np.linspace(0, 1, len(Ls)))

    mean_all_abs_m2_mean = []
    std_all_abs_m2_mean = []
    mean_all_E_mean = []
    std_all_E_mean = []
    mean_all_chi_mean = []
    std_all_chi_mean = []
    for idx, L in enumerate(Ls):
        color = colors[idx]
        T_vals, all_abs_m2_mean, all_E_mean, all_chi_mean = get_multi_run_data(folder, L, Ti, Tf, nT, sigma, method, doswap, Mrun)

        mean_all_abs_m2_mean.append(np.mean(all_abs_m2_mean, axis=0))
        std_all_abs_m2_mean.append(np.std(all_abs_m2_mean, axis=0))
        mean_all_E_mean.append(np.mean(all_E_mean, axis=0))
        std_all_E_mean.append(np.std(all_E_mean, axis=0))
        mean_all_chi_mean.append(np.mean(all_chi_mean, axis=0))
        std_all_chi_mean.append(np.std(all_chi_mean, axis=0))

        print("np.shape(T_vals)", np.shape(T_vals))
        print("np.shape(all_abs_m2_mean)", np.shape(all_abs_m2_mean))
        ax1.errorbar(T_vals, np.mean(all_abs_m2_mean, axis=0), yerr=np.std(all_abs_m2_mean, axis=0) / np.sqrt(Mrun), ls="-", marker="o", ms="3", color=color, label=f"{L:.0f}")
        ax2.errorbar(T_vals, np.mean(all_chi_mean, axis=0), yerr=np.std(all_chi_mean, axis=0) / np.sqrt(Mrun), ls="--", marker="o", ms="3", color=color, label=f"{L:.0f}")
        ax3.errorbar(T_vals, np.mean(all_E_mean, axis=0), yerr=np.std(all_E_mean, axis=0) / np.sqrt(Mrun), ls="-", marker="o", ms="3", color=color, label=f"{L:.0f}")

    ax1.set_ylim(0, 1.1)
    ax1.set_xlabel(r"$T$")
    ax1.set_ylabel(r"$<|m|>^2$")
    ax1.legend(title=r"$L$", loc="upper right", ncol=2, handlelength=0.5)
    ax2.set_xlabel(r"$T$")
    ax2.set_ylabel(r"$\chi = (<m^2> - <m>^2)/T$")
    ax2.legend(title=r"$L$", loc="upper right", ncol=2, handlelength=0.5)
    ax3.set_xlabel(r"$T$")
    ax3.set_ylabel(r"$<E>$")
    ax3.legend(title=r"$L$", loc="upper right", ncol=2, handlelength=0.5)

    # save the Tval, mean_abs_m2_mean, std_abs_m2_mean, mean_E_mean, std_E_mean, mean_chi_mean, std_chi_mean for all sigma to a csv file
    outputfile = os.path.join(folder, "mean_magnetization_multirun.csv")
    # Create a DataFrame to save the results
    df = pd.DataFrame()
    # First column is temperature values
    df["T"] = T_vals
    # Add data for each sigma
    for i, L in enumerate(Ls):
        df[f"L_{L:.2f}_m2_mean"] = mean_all_abs_m2_mean[i]
        df[f"L_{L:.2f}_m2_std"] = std_all_abs_m2_mean[i]
        df[f"L_{L:.2f}_E_mean"] = mean_all_E_mean[i]
        df[f"L_{L:.2f}_E_std"] = std_all_E_mean[i]
        df[f"L_{L:.2f}_chi_mean"] = mean_all_chi_mean[i]
        df[f"L_{L:.2f}_chi_std"] = std_all_chi_mean[i]
    # Save to CSV
    df.to_csv(outputfile, index=False)
    print(f"Saved data to {outputfile}")

    plt.tight_layout()

    output_fig = os.path.join(folder, "mean_magnetization_all_multirun_per_L.png")
    plt.savefig(output_fig)
    plt.show()  # Show the plot for interactive environments, if needed
    plt.close()
    print(f"Saved mean magnetization figure to {output_fig}")



def get_multi_run_hist(folder, L, Ti, Tf, nT, sigma, method, doswap, Mrun):
    # return mean_m2, std_m2, mean_E, std_E, mean_chi, std_chi versus all T

    # T_vals = np.geomspace(Ti, Tf, nT)
    T_vals = np.linspace(Ti, Tf, nT)
    all_E_per_T = [[] for i in range(len(T_vals))]
    all_m2_per_T = [[] for i in range(len(T_vals))]
    for run in range(Mrun):
        finfo = f"L{L:.0f}_Ti{Ti:.2f}_Tf{Tf:.2f}_nT{nT:.0f}_sigma{sigma:.2f}_method_{method}_doswap{doswap}_run_{run}"
        filename = f"{folder}/parallel_obs_{finfo}.txt"
        print(f"Analyzing file: {filename}")
        if not os.path.exists(filename):
            print(f"File {filename} does not exist.  for sigma = {sigma}")
            continue
        df = pd.read_csv(filename, skiprows=1)
        half_index = len(df) // 2
        for i in range(nT):
            abs_m = np.abs(df[f"m_{i}"].iloc[half_index:])
            abs_m_mean = abs_m.mean()  # This is <|m|>
            m2 = np.power(df[f"m_{i}"], 2)
            all_m2_per_T[i] += list(m2.iloc[half_index:])
            all_E_per_T[i] += list(df[f"E_{i}"].iloc[half_index:])

    return T_vals, all_m2_per_T, all_E_per_T



def analyze_parallel_histgram(folder, Ls, Ti, Tf, nT, sigma, method, doswap, Mrun):

    fig, axs = plt.subplots(nT, 2, figsize=(6, 2 * nT))
    colormap = plt.get_cmap("rainbow")
    colors = colormap(np.linspace(0, 1, len(Ls)))

    for idx, L in enumerate(Ls):

        color = colors[idx]
        T_vals, all_m2_per_T, all_E_per_T = get_multi_run_hist(folder, L, Ti, Tf, nT, sigma, method, doswap, Mrun)
        sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=min(Ls), vmax=max(Ls)))
        for i in range(nT):
            axs[i, 0].hist(all_m2_per_T[i], density=True, histtype="step", lw=1, color=color)
            axs[i, 1].hist(all_E_per_T[i], density=True, histtype="step", lw=1, color=color)
            axs[i, 0].text(0.5, 0.6, f"T = {T_vals[i]:.2f}", transform=axs[i,0].transAxes)
            axs[i, 1].text(0.5,0.6,f"T = {T_vals[i]:.2f}", transform=axs[i,1].transAxes)
            if(idx==0):
                fig.colorbar(sm, ax=axs[i,1], orientation="vertical", label="L", shrink=0.5)
                axs[i, 0].set_xlim(-0.1, 1.1)
                axs[i, 1].set_xlim(-10, 60)
                axs[i, 0].set_xlabel("m2", fontsize=9)
                axs[i, 0].set_ylabel("Density", fontsize=9)
                axs[i, 1].set_xlabel("E", fontsize=9)
                axs[i, 1].set_ylabel("Density", fontsize=9)


    plt.tight_layout()
    plt.savefig(f"{folder}/histogram_Ls.png")
    plt.show()
    plt.close()
    print(f"Saved histogram figure to {folder}/histogram_Ls.png")



def analyze_parallel_histgram_perT(folder, Ls, Ti, Tf, nT, sigma, method, doswap, Mrun):
    nL = len(Ls)
    fig, axs = plt.subplots(nL, 2, figsize=(6, 2 * nL))
    colormap = plt.get_cmap("rainbow")
    colors = colormap(np.linspace(0, 1, nT))

    for idx, L in enumerate(Ls):
        n_E_bin_per_T = []
        E_bin_center_per_T = []
        n_m2_bin_per_T = []
        m2_bin_center_per_T = []
        T_vals, all_m2_per_T, all_E_per_T = get_multi_run_hist(folder, L, Ti, Tf, nT, sigma, method, doswap, Mrun)
        sm = plt.cm.ScalarMappable(cmap=colormap, norm=plt.Normalize(vmin=min(T_vals), vmax=max(T_vals)))
        for i in range(nT):
            color = colors[i]
            n_m2, bins_m2_edges, _ =  axs[idx, 0].hist(all_m2_per_T[i], bins=100, density=True, histtype="step", lw=1, color=color)
            n_m2_bin_per_T.append(n_m2)
            m2_bin_center_per_T.append((bins_m2_edges[:-1] + bins_m2_edges[1:]) / 2)
            print("np.shape(n_m2), np.shape(m2_bin_center_per_T[-1])")
            print(np.shape(n_m2), np.shape(m2_bin_center_per_T[-1]))
            n_E, bins_E_edges, _ =  axs[idx, 1].hist(all_E_per_T[i], bins=100, density=True, histtype="step", lw=1, color=color)
            n_E_bin_per_T.append(n_E)
            E_bin_center_per_T.append((bins_E_edges[:-1] + bins_E_edges[1:]) / 2)
            if(i==0):
                axs[idx, 0].text(0.5, 0.6, f"L = {Ls[idx]:.0f}", transform=axs[idx,0].transAxes)
                axs[idx, 1].text(0.5,0.6,f"L = {Ls[idx]:.0f}", transform=axs[idx,1].transAxes)
                fig.colorbar(sm, ax=axs[idx,1], orientation="vertical", label="T", shrink=0.5)
                axs[idx, 0].set_xlim(-0.1, 1.1)
                axs[idx, 1].set_xlim(-10, 60)
                axs[idx, 0].set_xlabel("m2", fontsize=9)
                axs[idx, 0].set_ylabel("Density", fontsize=9)
                axs[idx, 1].set_xlabel("E", fontsize=9)
                axs[idx, 1].set_ylabel("Density", fontsize=9)

        df = pd.DataFrame()
        df["L"] = [L] * len(m2_bin_center_per_T[0])
        for i in range(nT):
            df[f"T_{T_vals[i]:.2f}_m2_bin_center"] = m2_bin_center_per_T[i]
            df[f"T_{T_vals[i]:.2f}_m2_bin_count"] = n_m2_bin_per_T[i]
            df[f"T_{T_vals[i]:.2f}_E_bin_center"] = E_bin_center_per_T[i]
            df[f"T_{T_vals[i]:.2f}_E_bin_count"] = n_E_bin_per_T[i]
        # Save to CSV
        outputfile = os.path.join(folder, f"histogram_L{L:.0f}.csv")
        df.to_csv(outputfile, index=False)
        print(f"Saved data to {outputfile}")

    plt.tight_layout()
    plt.savefig(f"{folder}/histogram_Ts.png")
    plt.show()
    plt.close()
    print(f"Saved histogram figure to {folder}/histogram_Ts.png")
