#include "r2_Ising.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <string>
#include <queue>
#include <algorithm>

// Constructor: initializes the 1D Ising model with 1/r^2 interactions and a random field.
r2_Ising::r2_Ising(int L_, double Ti_, double Tf_, int nT_, double sigma_, std::string method_)
{
    L = L_; // System size (L)
    T = Ti_;
    beta = 1.0 / T; // Inverse temperature
    sigma = sigma_; // Random field strength
    // Seed the generator
    std::random_device rd;
    std::mt19937 gen_(rd()); // Random number generator, seeded with rd()
    std::uniform_real_distribution<> rand_uni_(0, 1);
    std::normal_distribution<> rand_norm_(0, 1); // normal
    gen = gen_;
    rand_uni = rand_uni_;   // Uniform distribution for random numbers in [0, 1]
    rand_norm = rand_norm_; // Normal distribution for random fields

    // initialize spins
    Spins.resize(L); // Resize the Spins vector to hold L elements
    // Store the methodialization type
    method = method_;

    // Initialize random fields hi from a normal distribution with standard deviation sigma
    hi.resize(L);
    for (int i = 0; i < L; i++)
    {
        hi[i] = sigma * rand_norm(gen);
    }

    for (int i = 0; i < L; i++)
    {
        // initialize spins to aline with the field and add a small bias
        //Spins[i] = (rand_uni(gen) < 0.6)? 1: -1;
        Spins[i] = 1; // initialize all spins to +1
    }


    // Initialize interaction matrix Jij (using periodic boundary conditions)
    // Jij[i][j] = 1/(r_eff)^2, where r_eff = min(|i - j|, L - |i - j|)
    Jij.resize(L, std::vector<double>(L, 0.0));
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < L; j++)
        {
            if (i == j)
                continue; // No self-interaction
            int r = std::abs(i - j);
            int r_eff = std::min(r, L - r);
            if (r_eff > 0) // just to be safe
            {
                Jij[i][j] = 1.0 / (r_eff * r_eff);
            }
        }
    }

    // E = 1/2 * sum_ij J(i,j) (s_i-s_j)^2 - sum_i h_i s_i
    // E = 1/2 * sum_ij J(i,j) (s_i^2 + s_j^2 - 2 s_i s_j) - sum_i h_i s_i
    // E = sum_ij J(i,j) (1 - s_i s_j) - sum_i h_i s_i


    // Sum over pairs i<j for interactions and add the field contributions.
    E_sys = 0.0;
    for (int i = 0; i < L-1; i++)
    {
        for (int j = i + 1; j < L; j++)
        {
            E_sys += Jij[i][j] * (1- Spins[i] * Spins[j]);
        }
        E_sys += - hi[i] * Spins[i];
    }

    // initilize all replicas
    nT = nT_;
    T_replicas.resize(nT);
    beta_replicas.resize(nT);
    Spins_replicas.resize(nT);
    E_sys_replicas.resize(nT);
    T_replicas[0] = Ti_;
    beta_replicas[0] = 1.0 / T_replicas[0];

    double T_step = (Tf_ - Ti_) / (nT - 1);
    for (int i = 1; i < nT; i++)
    {
        T_replicas[i] = T_replicas[i - 1] + T_step;
        beta_replicas[i] = 1.0 / T_replicas[i];
    }
    /*
    double T_step = std::pow((Tf_ / Ti_), 1.0 / (nT - 1));

    for (int i = 1; i < nT; i++)
    {
        T_replicas[i] = T_replicas[i - 1] * T_step;
        beta_replicas[i] = 1.0 / T_replicas[i];
    }
        */

    for (int i = 0; i < nT; i++)
    {
        Spins_replicas[i] = Spins;
        E_sys_replicas[i] = E_sys;
    }
    std::cout << "T_replicas: ";
    for (int i = 0; i < nT; i++)
    {
        std::cout << T_replicas[i] << ",";
    }
    std::cout << "\n";
}

// MC_update: Performs one Metropolis update by attempting to flip a randomly selected spin.
// Energy change is computed as: ΔE = 2 s_i ( h_i + sum_{j≠i} J_{ij} s_j )
int r2_Ising::MC_update_single()
{
    // Select a random site i
    int i = rand_uni(gen) * L; // Randomly select a spin index in the range [0, L-1]
    int s_old = Spins[i];
    int s_new = -s_old;

    // Compute the energy change ΔE if spin i is flipped.
    double deltaE = 2 * s_old * hi[i];
    for (int j = 0; j < L; j++)
    {
        if (j == i)
            continue;
        deltaE += 2 * s_old * Jij[i][j] * Spins[j];
    }

    // Metropolis acceptance criterion
    if (deltaE <= 0 || rand_uni(gen) < std::exp(-beta * deltaE)) // to accelerate, no need to draw random number if deltaE <= 0
    {
        Spins[i] = s_new;
        E_sys += deltaE;
        return 1;
    }
    return 0;
}

// MC_update_cluster: Performs a Wolff cluster update with a Metropolis acceptance step in the presence of a random field
// The algorithm:
// 1. Choose a random seed spin and initialize the cluster with it.
// 2. For each spin in the cluster, consider every other spin. If the other spin has the same orientation,
//    add it to the cluster with probability P = 1 - exp(-2*beta*Jij[i][j]).
// 3. Compute the energy change from flipping the entire cluster, including field contributions and
//    interactions between cluster spins and spins outside the cluster:
//    ΔE = 2 * [sum_{i in cluster} hi[i]*s_i + sum_{i in cluster, j not in cluster} Jij[i][j]*s_i*s_j]
// 4. Accept the flip with the Metropolis criterion using ΔE.
// 5. If accepted, flip all spins in the cluster and update E_sys; return the number of spins flipped; otherwise, return 0.

int r2_Ising::MC_update_cluster()
{
    // Choose a random seed spin
    int seed = static_cast<int>(rand_uni(gen) * L);
    int seed_spin = Spins[seed];

    // Initialize cluster data structures
    std::vector<bool> in_cluster(L, false);
    std::vector<int> cluster;
    std::queue<int> to_check;

    // Add the seed to the cluster
    in_cluster[seed] = true;
    cluster.push_back(seed);
    to_check.push(seed);

    // Grow the cluster using the Wolff algorithm
    double p_add;
    while (!to_check.empty())
    {
        int i = to_check.front();
        to_check.pop();

        // Iterate over all other sites
        for (int j = 0; j < L; j++)
        {
            if (in_cluster[j])
                continue; // already in cluster
            // Only consider spins with the same orientation
            if (Spins[j] != Spins[i])
                continue;

            // Calculate bond probability
            //p_add = 1.0 - std::exp(-2 * beta * Jij[i][j]);
            p_add = 1.0 - std::exp(-2 * beta * Jij[i][j] - 2 * beta * hi[j] * Spins[j]); // H adjusted

            if (rand_uni(gen) < p_add)
            {
                in_cluster[j] = true;
                cluster.push_back(j);
                to_check.push(j);
            }
        }
    }

    // Compute energy change if the cluster is flipped
    // Field contribution: ΔE_field = 2 * sum_{i in cluster} hi[i] * s_i
    // Interaction contribution: ΔE_int = 2 * sum_{i in cluster, j not in cluster} Jij[i][j] * s_i * s_j
    double deltaE = 0.0;
    for (int i : cluster)
    {
        deltaE += 2 * hi[i] * Spins[i];
        for (int j = 0; j < L; j++)
        {
            if (!in_cluster[j])
            {
                deltaE += 2 * Jij[i][j] * Spins[i] * Spins[j];
            }
        }
    }

    // Metropolis acceptance criterion for the cluster flip
    if (deltaE <= 0 || rand_uni(gen) < std::exp(-beta * deltaE))
    {
        // Flip all spins in the cluster
        for (int i : cluster)
        {
            Spins[i] = -Spins[i];
        }
        E_sys += deltaE;
        return static_cast<int>(cluster.size());
    }

    return 0;
}

// measure_observable: Computes observables such as magnetization and energy.
observable r2_Ising::measure_observable()
{
    observable obs;
    int sum_m = 0;
    for (int i = 0; i < L; i++)
    {
        sum_m += Spins[i];
    }
    obs.m = sum_m / double(L); // Magnetization, m = sum(s_i) / L
    obs.E = E_sys;
    return obs;
}

// run_simulation: Runs the simulation over many MC steps, measures observables, and saves data.
void r2_Ising::run_simulation(int N, int M_sweep, std::string folder, std::string finfo)
{
    std::vector<observable> obs_ensemble;
    double acceptance_rate = 0.0;
    double sweep_acc_rate = 0.0;
    for (int n = 0; n < N; n++)
    {
        sweep_acc_rate = 0;
        for (int m = 0; m < M_sweep; m++)
        {
            // Perform a single Metropolis update, accumulate acceptance rate
            // sweep_acc_rate += MC_update_single(); // accumulate the number of accepted updates
            if (method == "single")
            {
                sweep_acc_rate += MC_update_single();
            }
            else if (method == "cluster")
            {
                sweep_acc_rate += MC_update_cluster(); // accumulate the number of accepted
            }
            else
            {
                std::cerr << "Unknown method: " << method << "\n";
            }
        }
        acceptance_rate += sweep_acc_rate / (double)M_sweep; // acceptance rate for this sweep
        obs_ensemble.push_back(measure_observable());

        // Print progress every 10% of the total steps
        if (n % (N / 10) == 0)
        {
            std::cout << "Progress: " << (n * 100.0) / N << "%" << std::endl;
        }
    }
    std::cout << obs_ensemble.size() << " observables to save\n";

    std::cout << "Final acceptance rate: "
              << (double)acceptance_rate / (double)N * 100.0
              << "%\n"; // Print the acceptance rate for debugging

    // Save observables and final configuration to files.
    std::string obs_filename = folder + "/obs_" + finfo + ".txt";
    save_observable_to_file(obs_filename, obs_ensemble);

    std::string config_filename = folder + "/config_" + finfo + ".txt";
    save_config_to_file(config_filename);
}

// save_config_to_file: Writes the current spin configuration to a file.
void r2_Ising::save_config_to_file(std::string filename)
{
    std::ofstream ofs(filename);
    if (!ofs)
    {
        std::cerr << "Error opening file " << filename << " for writing configuration.\n";
        return;
    }
    ofs << "s" << ",h\n";
    for (int i = 0; i < L; i++)
    {
        // Write each spin and its corresponding random field to the file
        ofs << Spins[i] << "," << hi[i] << "\n"; // format: s_i, h_i
    }

    ofs.close();
}

// save_observable_to_file: Writes the observable data (magnetization, m2, energy) to a file.
void r2_Ising::save_observable_to_file(std::string filename, std::vector<observable> obs_ensemble)
{
    std::ofstream ofs(filename);
    if (!ofs)
    {
        std::cerr << "Error opening file " << filename << " for writing observable data.\n";
        return;
    }
    std::cout << obs_ensemble.size() << " observables to save\n";
    ofs << "m,E\n"; // Header for the CSV file
    for (int i = 0; i < obs_ensemble.size(); i++)
    {
        // Write each observable to the file
        ofs << obs_ensemble[i].m << ","
            << obs_ensemble[i].E << "\n"; // format: m, m^2, E
    }
    ofs.close();
}

int r2_Ising::MC_update_single_replica(int rep)
{

    // Select a random site i
    int i = rand_uni(gen) * L; // Randomly select a spin index in the range [0, L-1]
    int s_old = Spins_replicas[rep][i];
    int s_new = -s_old;

    // Compute the energy change ΔE if spin i is flipped.
    double deltaE = 2 * s_old * hi[i];
    for (int j = 0; j < L; j++)
    {
        if (j == i)
            continue;
        deltaE += 2 * s_old * Jij[i][j] * Spins_replicas[rep][j];
    }

    // Metropolis acceptance criterion
    if (deltaE <= 0 || rand_uni(gen) < std::exp(-beta_replicas[rep] * deltaE)) // to accelerate, no need to draw random number if deltaE <= 0
    {
        Spins_replicas[rep][i] = s_new;
        E_sys_replicas[rep] += deltaE;
        return 1;
    }
    return 0;
}

// Returns the number of accepted exchanges.
int r2_Ising::MC_replica_exchange()
{
    int numAccepted = 0;
    // Assuming beta_replicas.size() is the number of replicas.
    int numReplicas = beta_replicas.size();
    // Loop over adjacent pairs of replicas.
    for (int i = 0; i < numReplicas - 1; i++)
    {
        // Calculate the acceptance probability:
        // Δ = (β_i - β_j) * (E_i - E_j)
        double delta = (beta_replicas[i] - beta_replicas[i + 1]) *
                       (E_sys_replicas[i] - E_sys_replicas[i + 1]);
        double swapProb = std::min(1.0, std::exp(delta));
        double r = rand_uni(gen);
        if (r < swapProb)
        {
            // Swap the configurations (Spins and energy) between replica i and i+1.
            std::swap(Spins_replicas[i], Spins_replicas[i + 1]);
            std::swap(E_sys_replicas[i], E_sys_replicas[i + 1]);
            numAccepted++;
        }
    }
    return numAccepted;
}

observable r2_Ising::measure_observable_replica(int rep)
{
    observable obs;
    int sum_m = 0;
    for (int i = 0; i < L; i++)
    {
        sum_m += Spins_replicas[rep][i];
    }
    obs.m = sum_m / double(L); // Magnetization, m = sum(s_i) / L
    obs.E = E_sys_replicas[rep];
    return obs;
}

void r2_Ising::run_parallel_simulation(int N, int M_sweep, std::string folder, std::string finfo, int doswap)
{
    int numReplicas = beta_replicas.size();

    // TODO: need to implement the following, need to check
    // Optionally, prepare observables storage for each replica.
    std::vector<std::vector<observable>> obs_ensemble_replica(numReplicas);

    std::vector<double> acceptance_rates(numReplicas, 0.0);
    // Main loop over Monte Carlo sweeps.
    double acceptance_rate_swap = 0.0;
    for (int sweep = 0; sweep < N; ++sweep)
    {
        // For each replica, perform M_sweep MC update steps.
        for (int rep = 0; rep < numReplicas; ++rep)
        {
            double acceptance_rate_sweep = 0;
            for (int step = 0; step < M_sweep; ++step)
            {

                // Example pseudo-code:
                std::vector<int> currentSpins = Spins_replicas[rep];
                double currentE = E_sys_replicas[rep];
                acceptance_rate_sweep += MC_update_single_replica(rep);
            }
            acceptance_rates[rep] += acceptance_rate_sweep / M_sweep; // Average acceptance rate for this replica
            // Optionally measure and store obs_ensemble_replica for this replica.
            obs_ensemble_replica[rep].push_back(measure_observable_replica(rep));
        }

        // Attempt replica exchanges among all replicas every 10 sweeps
        if(doswap && sweep % 10 == 1 )
        {
            acceptance_rate_swap += MC_replica_exchange();
        }

        // Print progress every 10% of the total steps
        if (sweep % (N / 10) == 0)
        {
            std::cout << "Progress: " << (sweep * 100.0) / N << "%" << std::endl;
        }
    }
    // (numReplicas-1) swap attempt per exchange update
    std::cout << "swap acceprance rate: "
              << (acceptance_rate_swap / N / (numReplicas - 1) * 100.0) << "%\n";
    for (int rep = 0; rep < numReplicas; ++rep)
    {
        std::cout << "Replica " << rep << " final acceptance rate: "
                  << (acceptance_rates[rep] * 10 / N * 100.0) << "%" << std::endl;
    }

    // At the end of the simulation, write the results to files.
    save_parallel_observable_to_file(folder, finfo, obs_ensemble_replica); // example for one replica
}
void r2_Ising::save_parallel_observable_to_file(std::string folder, std::string finfo, std::vector<std::vector<observable>> obs_ensemble_replicas)
{
    std::string filename = folder + "/parallel_obs_" + finfo + ".txt";
    std::ofstream ofs(filename);
    if (!ofs)
    {
        std::cerr << "Error opening file " << filename << " for writing parallel observable data.\n";
        return;
    }

    int numReplicas = obs_ensemble_replicas.size();
    int numSteps = obs_ensemble_replicas[0].size();

    // Write header line with temperatures
    ofs << "T";
    for (int rep = 0; rep < numReplicas; ++rep)
    {
        ofs << "," << T_replicas[rep];
    }
    ofs << "\n";

    // Write column headers
    for (int rep = 0; rep < numReplicas; ++rep)
    {
        ofs << "m_" << rep << ",E_" << rep;
        if (rep < numReplicas - 1)
            ofs << ",";
    }
    ofs << "\n";

    // Write data - each row corresponds to one MC step
    for (int step = 0; step < numSteps; ++step)
    {
        for (int rep = 0; rep < numReplicas; ++rep)
        {
            ofs << obs_ensemble_replicas[rep][step].m << ","
                << obs_ensemble_replicas[rep][step].E;
            if (rep < numReplicas - 1)
                ofs << ",";
        }
        ofs << "\n";
    }

    ofs.close();
    std::cout << "Saved observables for " << numReplicas << " replicas over " << numSteps << " steps to " << filename << std::endl;
}