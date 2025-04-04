# r2_Ising_RandH

This project simulates a one-dimensional Ising model with long-range 1/r² interactions and a random field using Monte Carlo methods. It supports both single-spin Metropolis updates and a Wolff cluster update algorithm with an additional Metropolis acceptance step to account for the random field contribution.

## Overview

The simulation models an Ising chain of spins with the following features:
	•	Interactions: The interaction strength between spins decays as 1/(r_eff)², where r_eff is computed using periodic boundary conditions.
	•	Random Field: Each spin experiences an independent random field drawn from a normal distribution scaled by a parameter σ.
	•	Updates: Two Monte Carlo update schemes are provided:
	•	Single-spin Metropolis Update (MC_update_single): Attempts to flip a randomly selected spin.
	•	Cluster Update (MC_update_cluster): Uses the Wolff algorithm to grow a cluster based on bond probabilities and then applies a Metropolis criterion to accept the cluster flip.
	•	Observables: At each measurement, the simulation computes the magnetization, its square, and the total energy.

## Files
	•	r2_Ising.h / r2_Ising.cpp: Contains the class definition and implementation for the r2_Ising model.
	•	main.cpp: Provides the main function to set up and run the simulation with user-defined parameters.
	•	Analysis Scripts (optional): Scripts (e.g., in Python) are available to analyze the observable data (m, m², E) as a function of Monte Carlo steps and to calculate statistical averages.

## Requirements
	•	A C++ compiler supporting C++11 (or later).
	•	CMake (optional, if using CMake for build automation).
	•	Python (with numpy, pandas, and matplotlib) for data analysis.

## Build Instructions

You can compile the project using your preferred build system. For example, using g++:

g++ -std=c++11 -O2 -o r2_Ising main.cpp r2_Ising.cpp

Or using CMake, create a CMakeLists.txt file and run:

mkdir build
cd build
cmake ..
make

## Running the Simulation

The executable expects five command-line arguments:
	1.	L: System size (number of spins).
	2.	β (beta): Inverse temperature.
	3.	σ (sigma): Random field strength.
	4.	init: Initialization type (ordered or random).
    5.  run_num: for tracking different realization of h field
	6.	folder: Folder where simulation outputs (observables and configuration files) will be saved.

Example usage:

./r2_Ising 100 1.0 0.5 random 0 output_folder

This command runs a simulation on a system of 100 spins with β = 1.0, σ = 0.5, using a random initial configuration, and saves the output in the output_folder folder.

Output Files
	•	obs_L{L}_beta{beta}_sigma{sigma}_init{init}.txt: CSV file containing observables (m, m², E) measured over the simulation.
	•	config_L{L}_beta{beta}_sigma{sigma}_init{init}.txt: File containing the final spin configuration and corresponding random fields.

## Data Analysis

For analysis, you can use the provided Python scripts. For instance, you might have a script that:
	1.	Plots time series: Displays m, m², and E versus Monte Carlo step for each combination of β and σ.
	2.	Computes Mean Magnetization: Plots |⟨m⟩| versus β for different values of σ on the same graph.

An example analysis function is provided in analyze.py which:
	•	Reads observable data from files.
	•	Generates subplots for the time series of m, m², and E.
	•	Computes and plots |⟨m⟩| vs. β for multiple σ values.

Ensure that your data files follow the naming convention described above for the analysis scripts to work correctly.
