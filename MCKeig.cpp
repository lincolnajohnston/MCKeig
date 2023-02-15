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

void runKEig(Geometry *geo) {
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

}

void runTransientFixedSource(Geometry *geo, double deltaT) {
    int initial_particles = 100000;
    int cycles = 100;
    int inactive_cycles = 50;
    // create particle bank
    std::vector<Particle> part_bank(initial_particles, Particle({0,0,0}, {1,0,0}, geo));
    for(size_t i = 0; i < initial_particles; i++) {
        part_bank[i] = (Particle({0,0,0}, {1,0,0}, geo));
    }

    for (double i = 0; i < deltaT * 3; i=i + deltaT) {
        std::vector<Particle> delayed_neutron_bank;
        double k_eig_sum = 0;
        //initial_particles = part_bank.size();

        for (size_t cycle = 0; cycle < cycles; cycle++) {
            std::cout << "--------------Cycle " << cycle << "-----------" << std::endl;
            std::cout << "Num Neutrons: " << part_bank.size() << std::endl;
            //std::vector<Particle> prompt_neutron_bank(initial_particles, Particle({0,0,0}, {1,0,0}, geo));
            std::vector<Particle> fission_source_bank;
            int source_particles = 0;
            // run particles
            for (size_t i = 0; i < part_bank.size(); i++) {
                Particle p = part_bank[i];
                //std::cout << "New Particle" << std::endl;
                
                while(p.isAlive()) {
                    if (cycle >= inactive_cycles) { // active cycles
                        p.move(fission_source_bank, delayed_neutron_bank, source_particles);
                    }
                    else {  // inactive cycles
                        p.move(fission_source_bank, fission_source_bank, source_particles);
                    }
                }
                //std::cout << "Killed\n\n" << std::endl;
            }
            std::cout << "Source particles: " << source_particles << std::endl;
            std::cout << "Delayed neutrons in bank: " << delayed_neutron_bank.size() << std::endl;
            double k = 1.0 * source_particles / part_bank.size();
            if (cycle >= inactive_cycles) {
                k_eig_sum += k;
            }
            std::cout << "k = " << k << " for this cycle\n" << std::endl;

            if (cycle >= inactive_cycles) {  // active cycles
                part_bank = fission_source_bank;
            }
            else {  // inactive cycles
                // Repeat particle locations if number of neutrons decreased in this generation
                int n_fission_source_neutrons = fission_source_bank.size();
                if (n_fission_source_neutrons < part_bank.size()) {
                    int new_particles = part_bank.size() - n_fission_source_neutrons;
                    for (size_t i = 0; i < new_particles; i++) {
                        int random_index = floor(Rand::getRand() * n_fission_source_neutrons);
                        fission_source_bank.push_back(fission_source_bank[random_index]);
                    }
                    part_bank = fission_source_bank;
                }
                // Randomly sample fission neutron locations if neutron population increased this generation
                else if (n_fission_source_neutrons > part_bank.size()) {
                    for (size_t i = 0; i < part_bank.size(); i++) {
                        int random_index = floor(Rand::getRand() * n_fission_source_neutrons);
                        part_bank[i] = fission_source_bank[random_index];
                    }
                }
            }
            
        }
        //part_bank = delayed_neutron_bank;
        // Repeat particle locations if number of neutrons decreased in this generation (kind of like splitting and rouletting without weight modification)
        int n_delayed_neutrons = delayed_neutron_bank.size();
        if (n_delayed_neutrons < initial_particles) {
            int new_particles = initial_particles - n_delayed_neutrons;
            for (size_t i = 0; i < new_particles; i++) {
                int random_index = floor(Rand::getRand() * n_delayed_neutrons);
                delayed_neutron_bank.push_back(delayed_neutron_bank[random_index]);
            }
            part_bank = delayed_neutron_bank;
        }
        // Randomly sample fission neutron locations if neutron population increased this generation
        else if (n_delayed_neutrons > initial_particles) {
            part_bank.clear();
            for (size_t i = 0; i < initial_particles; i++) {
                int random_index = floor(Rand::getRand() * n_delayed_neutrons);
                part_bank.push_back(delayed_neutron_bank[random_index]);
            }
        }


        double avg_k_eig = k_eig_sum / (cycles - inactive_cycles);
        std::cout << "Average k-eig from last " << cycles - inactive_cycles << " active cycles: " << avg_k_eig << std::endl;
        std::cout << "Total particles: " << part_bank.size() << std::endl;
    }
}


// What TFS should actually do I think: Given some intial source term (neutrons and their distribution), run some inactive cycles just like was done
// in the k-eigenvalue function to get the correct shape of the neutrons, then start running active cycles, during this time, the 
// delayed and censused neutrons will be added to a separate bank (with some probability to weight them into the current time step as well).
// If the reactor isn't prompt supercritical, then over time the number of neutrons in the fission source iteration (for k calculation)
// will decrease to a small value at the end of the active cycles. Then, we move to the next time step, where the delayed and censused
// neutrons will become the source term for the next fission source iteration. The shape should be relatively similar to the previous 
// time step but the inactive cycles will still be repeated (keeping the same initial neutron population but changing distribution).
// In a perfectly critical reactor, the delayed + censused neutron bank for each time step will have the same size equal to the size
// of the bank of neutrons in the initial source term

Geometry *buildGeometry (std::string simType, double deltaT) {

    // define sphere of radius 1 centered at the origin
    double sphereVelocity = 10000; // neutron velocity in cm/s
    double sphereScatterXS = 5;
    double sphereCapXS = 1.33;
    double sphereCensusXS = 1 / (deltaT * sphereVelocity);
    if (simType == "keig") {
        sphereCensusXS = 0;
    }
    double sphereFisXS = 1.67;
    double sphereNu = 2;
    double sphereTotalXS = sphereScatterXS + sphereCapXS + sphereCensusXS + sphereFisXS;
    double sphereBeta = 0.007;
    double sphereDecayConst = 1000; // units of decays per second
    Material *sphereMat = new Material(sphereTotalXS, sphereScatterXS, sphereCapXS, sphereCensusXS, sphereFisXS, sphereNu, sphereBeta, sphereDecayConst);
    Sphere *sphere = new Sphere({0,0,0}, 1);
    Cell *sphere_cell = new Cell("inner_sphere", {sphere}, {0}, sphereMat);

    // define a cube of side lengths 2, centered at the origin, cell excludes center sphere
    double cubeVelocity = 10000; // neutron velocity in cm/s
    double cubeScatterXS = 0;
    double cubeCapXS = 0;
    double cubeCensusXS = 1 / (deltaT * cubeVelocity);
    if (simType == "keig") {
        cubeCensusXS = 0;
    }
    double cubeFisXS = 0;
    double cubeNu = 2;
    double cubeTotalXS = cubeScatterXS + cubeCapXS + cubeCensusXS + cubeFisXS;
    double cubeBeta = 0.007;
    double cubeDecayConst = 0.1; // units of decays per second
    Material *cubeMat = new Material(cubeTotalXS, cubeScatterXS, cubeCapXS, cubeCensusXS, cubeFisXS, cubeNu, cubeBeta, cubeDecayConst);
    Plane *wall1 = new Plane({-1,0,0}, {1,0,0});
    Plane *wall2 = new Plane({1,0,0}, {-1,0,0});
    Plane *wall3 = new Plane({0,-1,0}, {0,1,0});
    Plane *wall4 = new Plane({0,1,0}, {0,-1,0});
    Plane *wall5 = new Plane({0,0,-1}, {0,0,1});
    Plane *wall6 = new Plane({0,0,1}, {0,0,-1});
    Cell *box = new Cell("outer_box", {wall1, wall2, wall3, wall4, wall5, wall6, sphere}, {1, 1, 1, 1, 1, 1, 1}, cubeMat);

    // add cells to geometry
    Geometry *geo = new Geometry({sphere_cell});
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

    if (simType == "keig") {
        std::cout << "Running K-eigenvalue simulation" << std::endl;
        runKEig(geo);
    }
    else if (simType == "TFS") {
        std::cout << "Running Transient Fixed Source simulation" << std::endl;
        runTransientFixedSource(geo, deltaT);
    }
    

    delete geo;
    return 0;
}



#endif