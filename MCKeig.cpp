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

/*void runKEig(Geometry *geo) {
    int initial_particles = 100000;
    int cycles = 100;
    int inactive_cycles = 50;
    double k_eig_sum = 0;
    // create particle bank
    std::vector<Particle> part_bank(initial_particles, Particle({0,0,0}, {1,0,0}, geo));
    for(size_t i = 0; i < initial_particles; i++) {
        part_bank[i] = (Particle({0,0,0}, {1,0,0}, geo));
    }

    for (size_t cycle = 0; cycle < cycles; cycle++) {
        std::cout << "--------------Cycle " << cycle << "-----------" << std::endl;
        std::cout << "Num Neutrons: " << part_bank.size() << std::endl;
        //std::vector<Particle> prompt_neutron_bank(initial_particles, Particle({0,0,0}, {1,0,0}, geo));
        std::vector<Particle> fission_source_bank;
        //std::vector<Particle> delayed_neutron_bank;
        int source_particles = 0;
        // run particles
        for (size_t i = 0; i < initial_particles; i++) {
            Particle p = part_bank[i];
            //std::cout << "New Particle" << std::endl;
            
            while(p.isAlive()) {
                p.move(fission_source_bank, fission_source_bank, source_particles);
            }
            //std::cout << "Killed\n\n" << std::endl;
        }
        double k = 1.0 * source_particles / initial_particles;
        if (cycle >= inactive_cycles) {
            k_eig_sum += k;
        }
        std::cout << "k = " << k << " for this cycle\n" << std::endl;

        // Repeat particle locations if number of neutrons decreased in this generation
        if (source_particles < initial_particles) {
            int new_particles = initial_particles - source_particles;
            for (size_t i = 0; i < new_particles; i++) {
                int random_index = floor(Rand::getRand() * source_particles);
                fission_source_bank.push_back(fission_source_bank[random_index]);
            }
            part_bank = fission_source_bank;
        }
        // Randomly sample fission neutron locations if neutron population increased this generation
        else if (source_particles > initial_particles) {
            for (size_t i = 0; i < initial_particles; i++) {
                int random_index = floor(Rand::getRand() * source_particles);
                part_bank[i] = fission_source_bank[random_index];
            }
        }
        
    }


    double avg_k_eig = k_eig_sum / (cycles - inactive_cycles);
    std::cout << "Average k-eig from last " << cycles - inactive_cycles << " active cycles: " << avg_k_eig << std::endl;

}*/

// lowerWeight must be less than 1
void weightWindows(std::vector<Particle> &bank, double lowerWeight, double upperWeight, double survivalWeight, Rand& rng) {
    std::vector<Particle> temp_bank;
    int i = 0;
    for (Particle p:bank) {
        i = i + 1;
        // roulette
        if(p.getWeight() < lowerWeight) {
            if(rng.getRand2() < p.getWeight() / survivalWeight){
                p.setWeight(survivalWeight);
                temp_bank.push_back(p);
            }
        }
        // split
        else if (p.getWeight() > upperWeight){
            double num_splits = round(p.getWeight() / survivalWeight);
            p.multiplyWeight(1.0 / num_splits);
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

double sumBankWeights(std::vector<Particle> &bank) {
    double sum = 0;
    for (Particle part : bank) {
        sum += part.getWeight();
    }
    return sum;
}

void runTransientFixedSource(Geometry *geo, double deltaT) {
    Rand rng;
    int initial_particles = 100;
    int delayed_bank_max_size = 100;
    int cycles = 100;
    int inactive_cycles = 50;
    // create particle bank for fission source iteration
    std::vector<Particle> part_bank(initial_particles, Particle({0,0,0}, {1,0,0}, geo));
    for(size_t i = 0; i < initial_particles; i++) {
        part_bank[i] = (Particle({0,0,0}, {1,0,0}, geo));
    }
    std::vector<Particle> delayed_neutron_bank; // neutrons passed to the next time step through delayed neutron or census

    // progress through time steps
    for (double t = 0; t < deltaT * 20; t=t + deltaT) {
        double k_eig_sum = 0;

        double inactive_weight = sumBankWeights(part_bank);

        // adjust weights of delayed neutrons for next time step (term 3 in Evan's dissertation TFS equation)
        for (Particle &p:delayed_neutron_bank) {
            p.f1WeightAdjust(deltaT);
        }

        // fission source iteration
        for (size_t cycle = 0; cycle < cycles; cycle++) {
            //std::cout << "--------------Cycle " << cycle << "-----------" << std::endl;
            //std::cout << "Num Neutrons: " << part_bank.size() << std::endl;
            std::vector<Particle> fission_source_bank; // neutrons that will be used in the next fission source iteration
            double source_particles = 0;

            // run particles in current fission source step, add generated neutrons to banks for either next fission source step or time step
            for (size_t i = 0; i < part_bank.size(); i++) {
                Particle p = part_bank[i];
                while(p.isAlive()) {
                    if (cycle >= inactive_cycles) { // active cycles (add neutrons to delayed bank)
                        p.move(fission_source_bank, delayed_neutron_bank, source_particles, deltaT, rng);
                    }
                    else {  // inactive cycles (don't add neutrons to delayed bank)
                        double trash_variable = 0;
                        p.move(fission_source_bank, fission_source_bank, trash_variable, deltaT, rng);
                    }
                }
            }
            //std::cout << "Source particles: " << source_particles << std::endl;
            //std::cout << "Delayed neutrons in bank: " << delayed_neutron_bank.size() << std::endl;
            double k = 1.0 * source_particles / sumBankWeights(part_bank);;
            if (cycle >= inactive_cycles) {
                k_eig_sum += k;
            }
            //std::cout << "k = " << k << " for this cycle\n" << std::endl;

            if (cycle >= inactive_cycles) {  // active cycles (particles are lost to delayed bank)
                part_bank = fission_source_bank;
            }
            else {  // inactive cycles (particles not lost to delayed bank, # of source particles kept constant)
                weightWindows(fission_source_bank, 0.01, 10, 1, rng);
                // Repeat particle locations if number of neutrons decreased in this generation
                int n_fission_source_neutrons = fission_source_bank.size();
                if (n_fission_source_neutrons < initial_particles) {
                    //int new_particles = part_bank.size() - n_fission_source_neutrons;
                    for (size_t i = 0; i < initial_particles; i++) {
                        int random_index = floor(rng.getRand2() * n_fission_source_neutrons);
                        fission_source_bank.push_back(fission_source_bank[random_index]);
                    }
                    part_bank = fission_source_bank;
                    double new_total_weight = sumBankWeights(part_bank);
                    for (Particle &p:part_bank) {
                        p.multiplyWeight(inactive_weight / new_total_weight);
                    }
                }
                // Randomly sample fission neutron locations if neutron population increased this generation
                else if (n_fission_source_neutrons > initial_particles) {
                    double total_fission_source_weight = sumBankWeights(fission_source_bank);
                    auto rng_shuffler = std::default_random_engine {};
                    std::shuffle(std::begin(fission_source_bank), std::end(fission_source_bank), rng_shuffler);
                    part_bank = std::vector<Particle>(fission_source_bank.begin(), fission_source_bank.begin() + initial_particles);
                    double new_total_weight = sumBankWeights(part_bank);
                    for (Particle &p:part_bank) {
                        p.multiplyWeight(inactive_weight / new_total_weight);
                    }
                }
            }

            double total_fission_source_weight = sumBankWeights(part_bank);
            std::cout << "part_bank weight: " << total_fission_source_weight << std::endl;

            //geo->printTallies();
            geo->clearTallies();
            
        }

        int n_delayed_neutrons = delayed_neutron_bank.size();
        double total_bank_weight = sumBankWeights(delayed_neutron_bank);
        std::cout << "Delayed bank total weight: " << total_bank_weight << std::endl;
        // Comb delayed neutron bank
        if (n_delayed_neutrons > delayed_bank_max_size) {
            auto rng_shuffler = std::default_random_engine {};
            std::shuffle(std::begin(delayed_neutron_bank), std::end(delayed_neutron_bank), rng_shuffler);
            delayed_neutron_bank = std::vector<Particle>(delayed_neutron_bank.begin(), delayed_neutron_bank.begin() + delayed_bank_max_size);
            double new_total_weight = sumBankWeights(delayed_neutron_bank);
            for (Particle &p:delayed_neutron_bank) {
                p.multiplyWeight(total_bank_weight / new_total_weight);
            }
        }

        total_bank_weight = sumBankWeights(delayed_neutron_bank);
        std::cout << "Delayed bank total weight: " << total_bank_weight << std::endl;
        // Weight windows: low-weight neutrons in delayed neutron bank
        weightWindows(delayed_neutron_bank, 0.01, 10, 1, rng);


        // add delayed neutrons to fission source iteration bank for next time step
        for (Particle p:delayed_neutron_bank) {
            part_bank.push_back(p);
        }

        int n_prompt_neutrons = part_bank.size();
        double total_fission_source_weight = sumBankWeights(part_bank);
        std::cout << "part_bank weight: " << total_fission_source_weight << std::endl;

        // Comb fission source particle neutron bank
        if (n_prompt_neutrons > initial_particles) {
            auto rng_shuffler = std::default_random_engine {};
            std::shuffle(std::begin(part_bank), std::end(part_bank), rng_shuffler);
            part_bank = std::vector<Particle>(part_bank.begin(), part_bank.begin() + initial_particles);
            double new_total_weight = sumBankWeights(part_bank);
            for (Particle &p:part_bank) {
                p.multiplyWeight(total_fission_source_weight / new_total_weight);
            }
        }

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
    std::string input_file = "input.txt";
    if (argc >= 3) {
        input_file = argv[2];
    }

    double deltaT = 0.01;

    Input input(input_file);
    Geometry *inputTestGeo = input.getGeometry();

    /*if (simType == "keig") {
        std::cout << "Running K-eigenvalue simulation" << std::endl;
        runKEig(geo);
    }
    else if (simType == "TFS") {
        std::cout << "Running Transient Fixed Source simulation" << std::endl;
        runTransientFixedSource(geo, deltaT);
    }*/
    //runTransientFixedSource(geo, deltaT);
    runTransientFixedSource(inputTestGeo, deltaT);

    return 0;
}



#endif