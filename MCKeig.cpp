#ifndef MCKEig
#define MCKEig

#include "Position.cpp"
#include "Particle.cpp"
#include "Surface.cpp"
#include "Cell.cpp"
#include "Geometry.cpp"
#include "Material.cpp"
#include "Source.cpp"
#include "Input.cpp"

#include <queue>
#include <vector>
#include <random>
#include <algorithm>
#include <memory>

// lowerWeight must be less than 1
void weightWindows(std::vector<std::shared_ptr<Particle>> &bank, double lowerWeight, double upperWeight, double survivalWeight, Rand& rng) {
    std::vector<std::shared_ptr<Particle>> temp_bank;
    int i = 0;
    for (std::shared_ptr<Particle> p:bank) {
        i = i + 1;
        // roulette
        if(p->getWeight() < lowerWeight) {
            if(rng.getRand() < p->getWeight() / survivalWeight){
                p->setWeight(survivalWeight);
                temp_bank.push_back(p);
            }
        }
        // split
        else if (p->getWeight() > upperWeight){
            double num_splits = round(p->getWeight() / survivalWeight);
            p->multiplyWeight(1.0 / num_splits);
            for(size_t i = 0; i < num_splits; i++) {
                temp_bank.push_back(p);
            }
        }
        else {
            temp_bank.push_back(p);
        }
    }
    bank = temp_bank;
}

double sumBankWeights(std::vector<std::shared_ptr<Particle>> &bank) {
    double sum = 0;
    for (std::shared_ptr<Particle> part : bank) {
        sum += part->getWeight();
    }
    return sum;
}

void comb(std::vector<std::shared_ptr<Particle>> &particle_bank, int n_particles) {
    int bank_size = particle_bank.size();
    if (bank_size > n_particles) {
        double bank_weight = sumBankWeights(particle_bank);
        auto rng_shuffler = std::default_random_engine {};
        std::shuffle(std::begin(particle_bank), std::end(particle_bank), rng_shuffler);
        particle_bank = std::vector<std::shared_ptr<Particle>>(particle_bank.begin(), particle_bank.begin() + n_particles);
        double new_total_weight = sumBankWeights(particle_bank);
        for (std::shared_ptr<Particle> p:particle_bank) {
            p->multiplyWeight(bank_weight / new_total_weight);
        }
    }
}

void runTransientFixedSource(Geometry *geo, double deltaT) {
    Rand rng;
    int initial_particles = 10000;
    int delayed_bank_max_size = 10000;
    int cycles = 100;
    int inactive_cycles = 50;

    // create bank for initial guess of fission source, B*psi_n
    std::vector<std::shared_ptr<Particle>> part_bank;
    Position initial_pos = {0,0,0};
    Direction initial_dir = Direction(1,0,0);
    for(size_t i = 0; i < initial_particles; i++) {
        part_bank.push_back(std::make_shared<Particle>(initial_pos, initial_dir, geo));
    }
    std::vector<std::shared_ptr<Particle>> delayed_neutron_bank; // delayed neutron source from the last time step, C*psi_(n-1)
    std::vector<std::shared_ptr<Particle>> new_delayed_neutron_bank; // delayed neutron source created the last time step, C*psi_(n-1)

    // progress through time steps
    for (double t = 0; t < deltaT * 100; t=t + deltaT) {
        double k_eig_sum = 0;
        double lifetime_sum = 0;
        std::vector<double> beta_eff_sum = {0,0,0,0,0,0};
        int adjoint_weighted_cycles = 0;

        double inactive_weight = sumBankWeights(part_bank);

        // fission source iteration
        for (size_t cycle = 0; cycle < cycles; cycle++) {
            std::cout << "--------------Cycle " << cycle << "-----------" << std::endl;
            std::cout << "Num Neutrons: " << part_bank.size() << std::endl;
            std::vector<std::shared_ptr<Particle>> fission_source_bank; // neutrons that will be used in the next fission source iteration
            double source_particles = 0;

            // run particles in current fission source step, add generated neutrons to banks for either next fission source step or time step
            for (size_t i = 0; i < part_bank.size(); i++) {
                std::shared_ptr<Particle> p = part_bank[i];
                while(p->isAlive()) {
                    p->move(fission_source_bank, new_delayed_neutron_bank, cycle >= inactive_cycles, source_particles, deltaT, rng);
                }
            }
            // run fission source from transient bank
            for (size_t i = 0; i < delayed_neutron_bank.size(); i++) {
                std::shared_ptr<Particle> p = delayed_neutron_bank[i];
                while(p->isAlive()) {
                    p->move(fission_source_bank, new_delayed_neutron_bank, cycle >= inactive_cycles, source_particles, deltaT, rng);
                }
            }
            //std::cout << "Source particles: " << source_particles << std::endl;
            //std::cout << "Delayed neutrons in bank: " << delayed_neutron_bank.size() << std::endl;
            double k = 1.0 * source_particles / sumBankWeights(part_bank);
            if (cycle >= inactive_cycles) {
                k_eig_sum += k;
            }
            std::cout << "k = " << k << " for cycle " << cycle << std::endl;

            // TODO: adjust particles' weights by k to maintain weight, then do splitting/routletting if necessary
            double weight_adjustment = initial_particles / sumBankWeights(fission_source_bank);
            part_bank = fission_source_bank; // set new fission source iteration bank

            // renormalize fission source bank weight
            for (std::shared_ptr<Particle> p:part_bank) {
                p->multiplyWeight(weight_adjustment);
            }

            comb(part_bank, initial_particles);

            double total_fission_source_weight = sumBankWeights(part_bank);
            //std::cout << "part_bank weight: " << total_fission_source_weight << "\n" << std::endl;

            // calculate "adjoint weighted" lifetime for 10 generations back
            if (cycle > inactive_cycles + 10) {
                double lifetimeTotal = 0;
                std::vector<double> betaTotals = {0,0,0,0,0,0};
                for (std::shared_ptr<Particle> p:part_bank) {
                    lifetimeTotal += p->getLifetimeContribution(10);
                    for (int del_group = 0; del_group < 6; del_group++) {
                        betaTotals[del_group] += p->getBetaContribution(del_group, 10);
                    }
                }
                //std::cout << "Effective Lifetime: " << lifetimeTotal / total_fission_source_weight << std::endl;
                //std::cout << "Beta Effective: " << betaTotal / total_fission_source_weight << std::endl;
                lifetime_sum += lifetimeTotal / total_fission_source_weight;
                for (int del_group = 0; del_group < 6; del_group++) {
                    beta_eff_sum[del_group] += betaTotals[del_group] / total_fission_source_weight;
                }
                adjoint_weighted_cycles++;
            }

            // print and clear tallies
            if (cycle == cycles - 1) {
                //geo->printFluxes();
            }
            if (cycle == inactive_cycles - 1) {
                geo->clearTallies();
            }
            
        }
        geo->clearTallies();

        delayed_neutron_bank.clear();
        // set up delayed neutron bank for next time step
        for (std::shared_ptr<Particle> p:new_delayed_neutron_bank) {
            delayed_neutron_bank.push_back(std::make_shared<Particle>(*p));
        }

        // adjust weights of delayed neutrons for next time step (term 3 in Evan's dissertation TFS equation)
        for (std::shared_ptr<Particle> p:new_delayed_neutron_bank) {
            p->f1WeightAdjust(deltaT); // hardcode group 1 as delayed neutron group, using 1 group delayed neutrons, use 6 in future
        }
        

        // Comb delayed neutron bank (keep size of delayed neutron vector under max value)
        double total_bank_weight = sumBankWeights(delayed_neutron_bank);
        std::cout << "Delayed bank total weight: " << total_bank_weight << std::endl;
        weightWindows(delayed_neutron_bank, 0.01, 10, 1, rng);
        comb(delayed_neutron_bank, delayed_bank_max_size);
        total_bank_weight = sumBankWeights(delayed_neutron_bank);
        std::cout << "Delayed bank total weight: " << total_bank_weight << std::endl;

        // Weight windows: low-weight neutrons in delayed neutron bank
        // weightWindows(delayed_neutron_bank, 0.01, 10, 1, rng);

        // Print results for this time step
        double avg_k_eig = k_eig_sum / (cycles - inactive_cycles);
        double avg_lifetime = lifetime_sum / adjoint_weighted_cycles;
        std::vector<double> avg_beta_eff = {0,0,0,0,0,0};
        for (int del_group = 0; del_group < 6; del_group++) {
            avg_beta_eff[del_group] = beta_eff_sum[del_group] / adjoint_weighted_cycles;
        }
        std::cout << "Average k-eig from last " << cycles - inactive_cycles << " active cycles: " << avg_k_eig << std::endl;
        std::cout << "Average lifetime from last " << adjoint_weighted_cycles << " active cycles: " << avg_lifetime << std::endl;
        for (int del_group = 0; del_group < 6; del_group++) {
            std::cout << "Average beta-eff (group " << del_group << ") from last " << adjoint_weighted_cycles << " active cycles: " << avg_beta_eff[del_group] << std::endl;
        }
        std::cout << "Total particles: " << part_bank.size() << std::endl;
    }
}

int main(int argc, char *argv[]) {

    // command line input
    std::string simType = "TFS";
    if (argc >= 2) {
        simType = argv[1];
    }
    std::string input_file = "Inputs/concentric-spheres-2.txt";
    if (argc >= 3) {
        input_file = argv[2];
    }

    double deltaT = 100;

    Input input(input_file);
    Geometry *inputTestGeo = input.getGeometry();

    runTransientFixedSource(inputTestGeo, deltaT);

    return 0;
}



#endif