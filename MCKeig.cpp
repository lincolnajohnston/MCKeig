#ifndef MCKEig
#define MCKEig

#include "Position.cpp"
#include "Particle.cpp"
#include "Surface.cpp"
#include "Cell.cpp"
#include "Geometry.cpp"
#include "Material.cpp"
#include "Source.cpp"

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
    int initial_particles = 10000;
    int delayed_bank_max_size = 10000;
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
            int source_particles = 0;

            // run particles in current fission source step, add generated neutrons to banks for either next fission source step or time step
            for (size_t i = 0; i < part_bank.size(); i++) {
                Particle p = part_bank[i];
                while(p.isAlive()) {
                    if (cycle >= inactive_cycles) { // active cycles (add neutrons to delayed bank)
                        p.move(fission_source_bank, delayed_neutron_bank, source_particles, deltaT, rng);
                    }
                    else {  // inactive cycles (don't add neutrons to delayed bank)
                        p.move(fission_source_bank, fission_source_bank, source_particles, deltaT, rng);
                    }
                }
            }
            //std::cout << "Source particles: " << source_particles << std::endl;
            std::cout << "Delayed neutrons in bank: " << delayed_neutron_bank.size() << std::endl;
            double k = 1.0 * source_particles / part_bank.size();
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


Geometry *buildGeometry (std::string simType, double deltaT) {

    Plane *centerline = new Plane({0,0,0}, {1,0,0});

    // define sphere of radius 1 centered at the origin
    double sphereVelocity = 10000; // neutron velocity in cm/s
    double sphereScatterXS = 5;
    double sphereCapXS = 1.37;
    double sphereCensusXS = 1 / (deltaT * sphereVelocity);
    if (simType == "keig") {
        sphereCensusXS = 0;
    }
    double sphereFisXS = 1.63;
    double sphereNu = 2;
    double sphereTotalXS = sphereScatterXS + sphereCapXS + sphereCensusXS + sphereFisXS;
    double sphereBeta = 0.01;
    double sphereDecayConst = 10; // units of decays per second
    Material *sphereMat = new Material(sphereTotalXS, sphereScatterXS, sphereCapXS, sphereCensusXS, sphereFisXS, sphereNu, sphereBeta, sphereDecayConst);
    Sphere *sphere = new Sphere({0,0,0}, 0.7937);
    Cell *sphere_cell = new Cell("inner_sphere_1", {sphere, centerline}, {0, 0}, sphereMat);
    Cell *sphere0_cell = new Cell("inner_sphere_2", {sphere, centerline}, {0, 1}, sphereMat);

    // define sphere of radius 1 centered at the origin
    double sphere2Velocity = 10000; // neutron velocity in cm/s
    double sphere2ScatterXS = 5;
    double sphere2CapXS = 1.37;
    double sphere2CensusXS = 1 / (deltaT * sphereVelocity);
    if (simType == "keig") {
        sphere2CensusXS = 0;
    }
    double sphere2FisXS = 1.63;
    double sphere2Nu = 2;
    double sphere2TotalXS = sphere2ScatterXS + sphere2CapXS + sphere2CensusXS + sphere2FisXS;
    double sphere2Beta = 0.01;
    double sphere2DecayConst = 10; // units of decays per second
    Material *sphere2Mat = new Material(sphere2TotalXS, sphere2ScatterXS, sphere2CapXS, sphere2CensusXS, sphere2FisXS, sphere2Nu, sphere2Beta, sphere2DecayConst);
    Sphere *sphere2 = new Sphere({0,0,0}, 1);
    Cell *sphere2_cell = new Cell("outer_sphere_1", {sphere2, sphere, centerline}, {0,1,0}, sphere2Mat);
    Cell *sphere3_cell = new Cell("outer_sphere_2", {sphere2, sphere, centerline}, {0,1,1}, sphere2Mat);

    // define a cube of side lengths 2, centered at the origin, cell excludes center sphere
    double cubeVelocity = 10000; // neutron velocity in cm/s
    double cubeScatterXS = 0;
    double cubeCapXS = 2;
    double cubeCensusXS = 1 / (deltaT * cubeVelocity);
    if (simType == "keig") {
        cubeCensusXS = 0;
    }
    double cubeFisXS = 2;
    double cubeNu = 2;
    double cubeTotalXS = cubeScatterXS + cubeCapXS + cubeCensusXS + cubeFisXS;
    double cubeBeta = 0.01;
    double cubeDecayConst = 10; // units of decays per second
    Material *cubeMat = new Material(cubeTotalXS, cubeScatterXS, cubeCapXS, cubeCensusXS, cubeFisXS, cubeNu, cubeBeta, cubeDecayConst);
    Plane *wall1 = new Plane({-1,0,0}, {1,0,0});
    Plane *wall2 = new Plane({1,0,0}, {-1,0,0});
    Plane *wall3 = new Plane({0,-1,0}, {0,1,0});
    Plane *wall4 = new Plane({0,1,0}, {0,-1,0});
    Plane *wall5 = new Plane({0,0,-1}, {0,0,1});
    Plane *wall6 = new Plane({0,0,1}, {0,0,-1});
    Cell *box = new Cell("outer_box", {wall1, wall2, wall3, wall4, wall5, wall6, sphere2}, {1, 1, 1, 1, 1, 1, 1}, cubeMat);

    // add cells to geometry
    Geometry *geo = new Geometry({sphere_cell, sphere0_cell, sphere2_cell, sphere3_cell, box});
    return geo;
}

int main(int argc, char *argv[]) {

    // Determine whether to run K-eigenvalue or Transient fixed source simulation based on command line input
    std::string simType = "TFS";
    if (argc >= 2) {
        simType = argv[1];
    }

    double deltaT = 0.01;
    Geometry *geo = buildGeometry(simType, deltaT);

    /*if (simType == "keig") {
        std::cout << "Running K-eigenvalue simulation" << std::endl;
        runKEig(geo);
    }
    else if (simType == "TFS") {
        std::cout << "Running Transient Fixed Source simulation" << std::endl;
        runTransientFixedSource(geo, deltaT);
    }*/
    runTransientFixedSource(geo, deltaT);

    delete geo;
    return 0;
}



#endif