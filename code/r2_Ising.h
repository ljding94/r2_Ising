#ifndef _R2_ISING_H
#define _R2_ISING_H
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

// observable
struct observable
{
    // geometric
    double m;      // magnetization, sum(s_i)/L
    double E; // energy,
};

struct bead
{
    std::vector<double> r{0, 0, 0}; // position
    std::vector<int> neighbor;      // keep track of neighbors
};

class r2_Ising
{
public:
    int L;                             // system size, LxL
    double beta;                     // inverse temperature
    double sigma;                   // random field strength
    std::string method; // initialization type, e.g., "random" or "ordered"
    std::vector<std::vector<double>> Jij; // interaction matrix lookup, Jij[i][j] = J(i,j)
    std::vector<double> hi; // random field, hi[i] = h(i), i.e. the random field at site i, size L

    std::vector<int> Spins;       // all Spins
    double E_sys; // total energy of the system, E = -sum<ij> J(i,j) s_i s_j + sum_i h_i s_i, J(i,j) = 0.5/(r_eff)^2, where r_eff = min(|i - j|, L - |i - j|)

    // randomnumber generators
    std::mt19937 gen;
    std::uniform_real_distribution<> rand_uni; // uniform distribution
    std::normal_distribution<> rand_norm;      // normal distribution

    // initialization
    r2_Ising(int L_, double beta_, double sigma_, std::string init_);

    int MC_update_single();
    int MC_update_cluster(); // Added missing semicolon

    observable measure_observable();

    void run_simulation(int N, int M_sweep, std::string folder, std::string finfo);

    // experiment
    void save_config_to_file(std::string filename);

    void save_observable_to_file(std::string filename, std::vector<observable> obs_ensemble);
};
#endif