import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def analyze_time_series(folder, L, beta, sigma, init, run):
    """
    Analyze the observable time series for a given combination of beta and sigma.
    The expected data file name is assumed to be: obs_beta_<beta>_sigma_<sigma>.txt
    with columns: m, m2, E
    """
    finfo = f"L{L:.0f}_beta{beta:.1f}_sigma{sigma:.1f}_init{init}_run{run}"
    # Construct filename based on beta and sigma
    filename = f"{folder}/obs_{finfo}.txt"
    if not os.path.exists(filename):
        print(f"File {filename} does not exist.")
        return

    # Load data using pandas
    df = pd.read_csv(filename)
    # Assume the file has a header: m,m2,E
    steps = np.arange(len(df))

    # Create a figure with three subplots for m, m2, and E vs MC step
    fig, axs = plt.subplots(3, 1, figsize=(5, 8), sharex=True)

    axs[0].plot(steps, df['m'], 'b-')
    axs[0].set_ylabel('m')
    axs[0].set_title(f'beta = {beta}, sigma = {sigma} : m vs MC step')

    axs[1].plot(steps, df['m2'], 'r-')
    axs[1].set_ylabel('m^2')
    axs[1].set_title(f'beta = {beta}, sigma = {sigma} : m^2 vs MC step')

    axs[2].plot(steps, df['E'], 'g-')
    axs[2].set_ylabel('E')
    axs[2].set_xlabel('MC step')
    axs[2].set_title(f'beta = {beta}, sigma = {sigma} : E vs MC step')

    plt.tight_layout()
    # Save the figure
    plt.savefig(f"{folder}/time_series_{finfo}.png")
    plt.close()


def analyze_mean_magnetization(folder, L, betas, sigmas, init, run):
    """
    For multiple sigma values, analyze the mean magnetization for each beta and plot |<m>| vs beta.
    Assumes data file for each beta is named: obs_L{L}_beta{beta}_sigma{sigma}_init{init}.txt
    """
    plt.figure(figsize=(4,6))
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)

    for sigma in sigmas:
        beta_vals = []
        mean_abs_m = []
        chi_vals = []

        for beta in betas:
            finfo = f"L{L:.0f}_beta{beta:.1f}_sigma{sigma:.1f}_init{init}_run{run}"
            filename = f"{folder}/obs_{finfo}.txt"
            if not os.path.exists(filename):
                print(f"File {filename} does not exist. Skipping beta = {beta} for sigma = {sigma}")
                continue
            df = pd.read_csv(filename)
            # Compute the average magnetization over MC steps
            # Skip first half of data for thermalization
            half_index = len(df) // 2
            m_mean = df['m'].iloc[half_index:].mean()
            m2_mean = df['m2'].iloc[half_index:].mean()
            chi = m2_mean - m_mean**2
            beta_vals.append(beta)
            mean_abs_m.append(np.abs(m_mean))
            chi_vals.append(chi)

        ax1.plot(1.0/np.array(beta_vals), mean_abs_m, 'o-', label=f'{sigma}')
        ax2.plot(1.0/np.array(beta_vals), chi_vals, 'o--', label=f'{sigma} ')

    ax1.set_xlabel(r'$T = 1/\beta$')
    ax1.set_ylabel('|<m>|')
    ax1.legend(title=r'$\sigma$', loc='upper right')
    ax2.set_xlabel(r'$T = 1/\beta$')
    ax2.set_ylabel(r'$\chi = <m^2> - <m>^2$')
    ax2.legend(title=r'$\sigma$', loc='upper right')
    plt.tight_layout()
    output_fig = os.path.join(folder, "mean_magnetization_all.png")
    plt.savefig(output_fig)
    plt.close()
    print(f"Saved mean magnetization figure to {output_fig}")


def main():
    # Example usage: specify the data folder and lists of beta and sigma values
    # Adjust these lists as needed
    folder = 'data'  # Folder where the observable files are stored
    betas = [0.5, 1.0, 1.5, 2.0]  # List of inverse temperature values
    sigmas = [0.1, 0.5]           # List of random field strengths

    # Create time series figures for each combination of beta and sigma
    for sigma in sigmas:
        for beta in betas:
            analyze_time_series(folder, beta, sigma)

    # Create mean magnetization plot for each sigma
    analyze_mean_magnetization(folder, betas, sigmas)


if __name__ == '__main__':
    main()
