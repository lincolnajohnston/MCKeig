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

void runTransientFixedSource(Geometry *geo, double deltaT) {
    Rand rng;
    int initial_particles = 10000;
    int delayed_bank_max_size = 10000;
    int cycles = 100;
    int inactive_cycles = 50;
    // create particle bank for fission source iteration
    // std::vector<std::shared_ptr<Particle>> part_bank(initial_particles, std::make_shared<Particle>({0,0,0}, {1,0,0}, geo));
    std::vector<std::shared_ptr<Particle>> part_bank;
    //std::cout << part_bank[0].get()->getWeight() << std::endl;
    Position initial_pos = {0,0,0};
    Direction initial_dir = Direction(1,0,0);
    for(size_t i = 0; i < initial_particles; i++) {
        part_bank.push_back(std::make_shared<Particle>(initial_pos, initial_dir, geo));
    }
    std::vector<std::shared_ptr<Particle>> delayed_neutron_bank; // neutrons passed to the next time step through delayed neutron or census

    // progress through time steps
    for (double t = 0; t < deltaT * 20; t=t + deltaT) {
        double k_eig_sum = 0;

        double inactive_weight = sumBankWeights(part_bank);
        double beta_adjoint_sum = 0;
        double adjoint_sum = 0;

        // adjust weights of delayed neutrons for next time step (term 3 in Evan's dissertation TFS equation)
        for (std::shared_ptr<Particle> p:delayed_neutron_bank) {
            p->f1WeightAdjust(deltaT, 1); // hardcode group 1 as delayed neutron group
        }

        // fission source iteration
        for (size_t cycle = 0; cycle < cycles; cycle++) {
            //std::cout << "--------------Cycle " << cycle << "-----------" << std::endl;
            //std::cout << "Num Neutrons: " << part_bank.size() << std::endl;
            std::vector<std::shared_ptr<Particle>> fission_source_bank; // neutrons that will be used in the next fission source iteration
            double source_particles = 0;

            // run particles in current fission source step, add generated neutrons to banks for either next fission source step or time step
            for (size_t i = 0; i < part_bank.size(); i++) {
                std::shared_ptr<Particle> p = part_bank[i];
                while(p->isAlive()) {
                    if (cycle >= inactive_cycles) { // active cycles (add neutrons to delayed bank)
                        p->move(fission_source_bank, delayed_neutron_bank, source_particles, deltaT, rng);
                    }
                    else {  // inactive cycles (don't add neutrons to delayed bank)
                        double trash_variable = 0;
                        p->move(fission_source_bank, fission_source_bank, source_particles, deltaT, rng);
                    }
                }
            }
            //std::cout << "Source particles: " << source_particles << std::endl;
            //std::cout << "Delayed neutrons in bank: " << delayed_neutron_bank.size() << std::endl;
            double k = 1.0 * source_particles / sumBankWeights(part_bank);
            if (cycle >= inactive_cycles) {
                k_eig_sum += k;
            }
            std::cout << "k = " << k << " for this cycle" << std::endl;

            if (cycle >= inactive_cycles) {  // active cycles (particles are lost to delayed bank)
                part_bank = fission_source_bank;
            }
            else {  // inactive cycles (particles not lost to delayed bank, # of source particles kept constant)
                weightWindows(fission_source_bank, 0.01, 10, 1, rng);
                // Repeat particle locations if number of neutrons decreased in this generation
                int n_fission_source_neutrons = fission_source_bank.size();
                if (n_fission_source_neutrons < initial_particles) {
                    int new_particles = part_bank.size() - n_fission_source_neutrons;
                    for (size_t i = 0; i < new_particles; i++) {
                        int random_index = floor(rng.getRand() * n_fission_source_neutrons);
                        fission_source_bank.push_back(std::make_shared<Particle>(*fission_source_bank[random_index]));
                    }
                    part_bank = fission_source_bank;
                    double new_total_weight = sumBankWeights(part_bank);
                    for (std::shared_ptr<Particle> p:part_bank) {
                        p->multiplyWeight(inactive_weight / new_total_weight);
                    }
                }
                // Randomly sample fission neutron locations if neutron population increased this generation
                else if (n_fission_source_neutrons > initial_particles) {
                    double total_fission_source_weight = sumBankWeights(fission_source_bank);
                    auto rng_shuffler = std::default_random_engine {};
                    std::shuffle(std::begin(fission_source_bank), std::end(fission_source_bank), rng_shuffler);
                    part_bank = std::vector<std::shared_ptr<Particle>>(fission_source_bank.begin(), fission_source_bank.begin() + initial_particles);
                    double new_total_weight = sumBankWeights(part_bank);
                    for (std::shared_ptr<Particle> p:part_bank) {
                        p->multiplyWeight(inactive_weight / new_total_weight);
                    }
                }
            }

            double total_fission_source_weight = sumBankWeights(part_bank);
            std::cout << "part_bank weight: " << total_fission_source_weight << "\n" << std::endl;

            if (cycle == cycles - 1) {
                geo->printFluxes();
            }
            if (cycle == inactive_cycles - 1) {
                geo->clearTallies();
            }
            
        }
        double beta_eff;
        double lambda_eff;
        geo->getAverageBetaAndLambda(beta_eff, lambda_eff);
        std::cout << "Beta-eff: " << beta_eff << std::endl;
        std::cout << "Lambda-eff: " << lambda_eff << std::endl;

        geo->clearTallies();

        int n_delayed_neutrons = delayed_neutron_bank.size();
        double total_bank_weight = sumBankWeights(delayed_neutron_bank);
        std::cout << "Delayed bank total weight: " << total_bank_weight << std::endl;
        // Comb delayed neutron bank
        if (n_delayed_neutrons > delayed_bank_max_size) {
            auto rng_shuffler = std::default_random_engine {};
            std::shuffle(std::begin(delayed_neutron_bank), std::end(delayed_neutron_bank), rng_shuffler);
            delayed_neutron_bank = std::vector<std::shared_ptr<Particle>>(delayed_neutron_bank.begin(), delayed_neutron_bank.begin() + delayed_bank_max_size);
            double new_total_weight = sumBankWeights(delayed_neutron_bank);
            for (std::shared_ptr<Particle> p:delayed_neutron_bank) {
                p->multiplyWeight(total_bank_weight / new_total_weight);
            }
        }

        total_bank_weight = sumBankWeights(delayed_neutron_bank);
        std::cout << "Delayed bank total weight: " << total_bank_weight << std::endl;
        // Weight windows: low-weight neutrons in delayed neutron bank
        weightWindows(delayed_neutron_bank, 0.01, 10, 1, rng);


        // add delayed neutrons to fission source iteration bank for next time step
        for (std::shared_ptr<Particle> p:delayed_neutron_bank) {
            part_bank.push_back(std::make_shared<Particle>(*p));
        }

        int n_prompt_neutrons = part_bank.size();
        double total_fission_source_weight = sumBankWeights(part_bank);
        std::cout << "part_bank weight: " << total_fission_source_weight << std::endl;

        // Comb fission source particle neutron bank
        if (n_prompt_neutrons > initial_particles) {
            auto rng_shuffler = std::default_random_engine {};
            std::shuffle(std::begin(part_bank), std::end(part_bank), rng_shuffler);
            part_bank = std::vector<std::shared_ptr<Particle>>(part_bank.begin(), part_bank.begin() + initial_particles);
            double new_total_weight = sumBankWeights(part_bank);
            for (std::shared_ptr<Particle> p:part_bank) {
                p->multiplyWeight(total_fission_source_weight / new_total_weight);
            }
        }

        // Print eigenvalue results for this time step
        double avg_k_eig = k_eig_sum / (cycles - inactive_cycles);
        std::cout << "Average k-eig from last " << cycles - inactive_cycles << " active cycles: " << avg_k_eig << std::endl;
        std::cout << "Total particles: " << part_bank.size() << std::endl;
    }
}

int main(int argc, char *argv[]) {

    // command line input
    std::string simType = "TFS";
    if (argc >= 2) {
        simType = argv[1];
    }
    std::string input_file = "Inputs/concentric-spheres.txt";
    if (argc >= 3) {
        input_file = argv[2];
    }

    double deltaT = 0.01;

    Input input(input_file);
    Geometry *inputTestGeo = input.getGeometry();

    runTransientFixedSource(inputTestGeo, deltaT);

    return 0;
}



#endif