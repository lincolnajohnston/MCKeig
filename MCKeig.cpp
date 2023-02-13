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
    int initial_particles = 1000000;
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
        std::vector<Particle> prompt_neutron_bank;
        std::deque<Particle> delayed_neutron_bank;
        int source_particles = 0;
        // run particles
        for (size_t i = 0; i < initial_particles; i++) {
            Particle p = part_bank[i];
            //std::cout << "New Particle" << std::endl;
            
            while(p.isAlive()) {
                p.move(prompt_neutron_bank, delayed_neutron_bank, source_particles);
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
                prompt_neutron_bank[source_particles + i] = prompt_neutron_bank[random_index];
            }
            part_bank = prompt_neutron_bank;
        }
        // Randomly sample fission neutron locations if neutron population increased this generation
        else if (source_particles > initial_particles) {
            for (size_t i = 0; i < initial_particles; i++) {
                int random_index = floor(Rand::getRand() * source_particles);
                part_bank[i] = prompt_neutron_bank[random_index];
            }
        }
        
    }

    double avg_k_eig = k_eig_sum / (cycles - inactive_cycles);
    std::cout << "Average k-eig from last " << cycles - inactive_cycles << " active cycles: " << avg_k_eig << std::endl;

}

void runTransientFixedSource(Geometry *geo, double deltaT) {
    int initial_particles = 100;  // initial particles spawned at the origin
    int source_rate = 5000;  // particles generated at every time step
    int cycles = 200;
    int inactive_cycles = 50;
    double k_eig_sum = 0;
    // create initial particle bank
    std::vector<Particle> part_bank(initial_particles, Particle({0,0,0}, {1,0,0}, geo));
    for(size_t i = 0; i < initial_particles; i++) {
        part_bank[i] = (Particle({0,0,0}, {1,0,0}, geo));
    }
    std::deque<Particle> delayed_neutron_bank;

    // run timesteps
    int generated_particles = 0;
    for (size_t cycle = 0; cycle < cycles; cycle++) {
        std::cout << "--------------t = " << cycle * deltaT << "-----------" << std::endl;
        std::cout << "Num Neutrons: " << part_bank.size() << "\n\n\n" << std::endl;
        //std::vector<Particle> prompt_neutron_bank(generated_particles, Particle({0,0,0}, {1,0,0}, geo));
        //std::vector<Particle> delayed_neutron_bank(generated_particles, Particle({0,0,0}, {1,0,0}, geo));
        std::vector<Particle> generated_neutron_bank;
        
        // generate fixed source neutrons
        Source *source = new Source(source_rate, geo);
        source->generateParticles(part_bank);

        initial_particles = part_bank.size();
        generated_particles = 0;
        // move particles in this generation
        for (size_t i = 0; i < initial_particles; i++) {
            Particle p = part_bank[i];
            //std::cout << "New Particle" << std::endl;
            
            while(p.isAlive()) {
                p.move(generated_neutron_bank, delayed_neutron_bank, generated_particles);
            }
            //std::cout << "Killed\n\n" << std::endl;
        }
        
        // add delayed neutrons to bank
        for (int i = delayed_neutron_bank.size() - 1; i >= 0; i--) {
            Particle p = delayed_neutron_bank[i];
            
            // move particle from delayed to prompt group with probability based on decay constant of material of cell particle is in
            double decayProb = 1 - exp(-1 * p.getCell()->getDecayConst() * deltaT);
            if (Rand::getRand() < decayProb) {
                generated_neutron_bank.push_back(p);
                delayed_neutron_bank.erase(delayed_neutron_bank.begin() + i);
            }

        }
        std::cout << "Delayed Neutron Length: " << delayed_neutron_bank.size() << std::endl;

        part_bank = generated_neutron_bank;
    }
}

Geometry *buildGeometry (std::string simType, double deltaT) {

    // define sphere of radius 1 centered at the origin
    double sphereVelocity = 10000; // neutron velocity in cm/s
    double sphereScatterXS = 5;
    double sphereCapXS = 1.5;
    double sphereCensusXS = 1 / (deltaT * sphereVelocity);
    if (simType == "keig") {
        sphereCensusXS = 0;
    }
    double sphereFisXS = 1.5;
    double sphereNu = 2;
    double sphereTotalXS = sphereScatterXS + sphereCapXS + sphereCensusXS + sphereFisXS;
    double sphereBeta = 0.07;
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