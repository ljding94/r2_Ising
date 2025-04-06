// Copyright[2024] [Lijie Ding]
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <chrono>
#include <filesystem>
#include "r2_Ising.h"

// todo:
// 1. implement cluster update, form cluster using Wolff
// 2. implement simulated tempering,

int main(int argc, char const *argv[])
{
    std::clock_t c_start = std::clock();

    std::cout << "running with input argc:" << argc << "\n";
    for (int i = 0; i < argc; i++)
    {
        std::cout << argv[i] << " ";
    }

    std::string folder;
    std::string finfo;
    int L;
    double beta;
    double T;
    double sigma;
    int run_num;
    std::string init;
    // precision run with specified parameters
    if (argc == 7)
    {
        L = std::atoi(argv[1]);      // system size
        T = std::atof(argv[2]);   // inverse temperature
        beta = 1.0 / T; // convert temperature to beta, beta = 1/T
        sigma = std::atof(argv[3]);  // random field strength
        init = std::string(argv[4]); // initialization type, e.g., "random" or "ordered
        run_num = std::atoi(argv[5]);
        folder = std::string(argv[6]);
        finfo = "L" + std::string(argv[1]) + "_T" + std::string(argv[2]) + "_sigma" + std::string(argv[3]) + "_init" + init + "_run" + std::string(argv[5]); // filename info
    }
    else
    {
        std::cout << "input error\n";
        return 0;
    }
    std::cout << "Running simulation for L = " << L << ", T = " << T << ", sigma = " << sigma << ", init = " << init << "\n";
    int N = 10000;                              // number of sweep
    int M_sweep = 10 * L;                   // update per sweep
    r2_Ising r2_ising_1d(L, beta, sigma, init); // create an instance of the r2_Ising class

    r2_ising_1d.run_simulation(N, M_sweep, folder, finfo); // run the simulation
    std::cout << "Simulation completed. Results saved in folder: " << folder << "\n";

    std::clock_t c_end = std::clock();
    double time_elapsed = static_cast<double>(c_end - c_start) / CLOCKS_PER_SEC;
    std::cout << "\nTime elapsed: " << time_elapsed << " seconds" << std::endl;
    return 0;
}

// run time reference: N=10000 L=100, M=L*L, MC_update_single: 11.8s
// N1000 L=100, M_sweep=L, MC_update_cluster: 20~40s