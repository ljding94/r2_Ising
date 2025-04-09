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
    double Ti, Tf;
    int nT; // number of temperature steps
    double sigma;
    int run_num;
    std::string method;
    // precision run with specified parameters
    if (argc == 9)
    {
        L = std::atoi(argv[1]);        // system size
        Ti = std::atof(argv[2]);       // starting temperature
        Tf = std::atof(argv[3]);       // final temperature
        nT = std::atoi(argv[4]);       // number of temperature steps
        sigma = std::atof(argv[5]);    // random field strength
        method = std::string(argv[6]);
        run_num = std::atoi(argv[7]);
        folder = std::string(argv[8]);
        finfo = "L" + std::string(argv[1]) + "_Ti" + std::string(argv[2]) +
                "_Tf" + std::string(argv[3]) + "_nT" + std::to_string(nT) +
                "_sigma" + std::string(argv[5]) +
                "_method_" + method +
                "_run_" + std::string(argv[7]);
    }
    else
    {
        std::cout << "input error\n";
        return 0;
    }
    std::cout << "Running simulation for L = " << L << ", Ti = " << Ti
              << ", Tf = " << Tf << ", nT = " << nT
              << ", sigma = " << sigma
              << ", method = " << method
              << ", folder = " << folder
              << ", finfo = " << finfo
              << "\n";
    int N = 10000;        // number of sweep
    int M_sweep = 10 * L; // update per sweep
    if (method == "single")
    {
        M_sweep = 100 * L * L;
    }

    r2_Ising r2_ising_1d(L, Ti, Tf, nT, sigma, method); // create an instance of the r2_Ising class

    if (nT == 1)
    {
        r2_ising_1d.run_simulation(N, M_sweep, folder, finfo); // run the simulation
        std::cout << "Simulation completed. Results saved in folder: " << folder << "\n";
    } else {
        r2_ising_1d.run_parallel_simulation(N, M_sweep, folder, finfo);
        std::cout << "Parallel simulation completed. Results saved in folder: " << folder << "\n";
    }

    std::clock_t c_end = std::clock();
    double time_elapsed = static_cast<double>(c_end - c_start) / CLOCKS_PER_SEC;
    std::cout << "\nTime elapsed: " << time_elapsed << " seconds" << std::endl;
    return 0;
}

// run time reference: N=10000 L=100, M=L*L, MC_update_single: 11.8s
// N1000 L=100, M_sweep=L, MC_update_cluster: 20~40s

// parallel running N10000, L50, nT 5, M=10*L*L : 100s