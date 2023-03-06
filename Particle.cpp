#ifndef PARTICLE
#define PARTICLE

#include "Position.cpp"
#include "Cell.cpp"
#include "Geometry.cpp"
#include <queue>
#include <vector>

class Particle {
    protected:
        Position location;
        Direction direction;
        bool alive;
        Geometry *geo;
        Cell *currentCell;
        double epsilon = pow(10,-10);
        double weight;

    public:

    Particle(Position location, Direction direction, Geometry *geo, double weight = 1):location(location), direction(direction),
         geo(geo), weight(weight) {
        currentCell = geo->cellAtLocation(location);
        alive = true;
    }
    Position getLocation() {
        return location;
    }

    Cell *getCell() {
        return currentCell;
    }

    void setWeight(double newWeight) {
        weight = newWeight;
    }

    void multiplyWeight(double factor) {
        weight = weight * factor;
    }

    void move(std::vector<Particle> &generated_neutron_bank, std::vector<Particle> &delayed_neutron_bank, int &source_particles, double timestep, Rand& rng) {
        double edge_dist = currentCell->distToEdge(location, direction);
        double collision_dist = currentCell->distToNextCollision();

        double move_dist = edge_dist + epsilon;
        if (collision_dist < edge_dist) {
            move_dist = collision_dist;

            location.x += move_dist * direction.getI();
            location.y += move_dist * direction.getJ();
            location.z += move_dist * direction.getK();

            currentCell->tallyTL(move_dist, weight);


            std::string collision_name = currentCell->sample_collision(rng);
            if (collision_name == "scat") {
                direction.isotropicScatter(rng);
            }
            else if (collision_name == "cap") {
                alive = false;
            }
            else if (collision_name == "census") {
                delayed_neutron_bank.push_back(*this);

                source_particles += 1;
                alive = false;
            }
            else if (collision_name == "fis") {
                for (size_t i = 0; i < currentCell->getNu(); i++) {
                    // implicit delayed and prompt neutrons
                    /*direction.isotropicScatter();
                    Particle delayed_particle = *this;
                    direction.isotropicScatter();
                    Particle prompt_particle = *this;
                    delayed_particle.multiplyWeight(currentCell->getBeta());
                    prompt_particle.multiplyWeight(1 - currentCell->getBeta());
                    delayed_neutron_bank.push_back(delayed_particle);
                    generated_neutron_bank.push_back(prompt_particle);*/

                    // analog case
                    direction.isotropicScatter(rng);
                    double ksi = rng.getRand2();
                    if (ksi < currentCell->getBeta()) {
                        Particle fission_source_del_neut = *this;
                        fission_source_del_neut.f3WeightAdjust(timestep);
                        Particle time_bank_del_neut = *this;
                        time_bank_del_neut.f2WeightAdjust(timestep);

                        generated_neutron_bank.push_back(fission_source_del_neut);
                        delayed_neutron_bank.push_back(time_bank_del_neut);
                    }
                    else {
                        generated_neutron_bank.push_back(*this);
                    }
                    
                    source_particles += 1;
                }
                alive = false;
            }
            
        }
        else {
            location.x += move_dist * direction.getI();
            location.y += move_dist * direction.getJ();
            location.z += move_dist * direction.getK();

            currentCell->tallyTL(move_dist, weight);

            currentCell = geo->cellAtLocation(location);

            // if particle leaked, kill it
            if (!currentCell) {
                alive = false;
            }
        }
        
    }

    double getWeight() {
        return weight;
    }

    // adjust weight by factor equal to probability of precursor surviving the length of the time step
    void f1WeightAdjust(double timestep) {
        multiplyWeight(exp(-1 * currentCell->getDecayConst() * timestep));
    }

    // adjust weight by probability of precursor decaying within this time step and being being considered as a delayed neutron for the next time step
    void f2WeightAdjust(double timestep) {
        double lambda = currentCell->getDecayConst();
        multiplyWeight((1 - exp(-1 * lambda * timestep) - lambda * timestep * exp(-1 * lambda * timestep)) / (lambda * timestep));
    }

    // adjust weight by probability of precursor decaying within this time step and being being considered as a delayed neutron for the current fission source iteration
    void f3WeightAdjust(double timestep) {
        double lambda = currentCell->getDecayConst();
        multiplyWeight(exp(-1 * lambda * timestep) * (1 - exp(lambda * timestep) + lambda * timestep * exp(lambda * timestep)) / (lambda * timestep));
    }

    bool isAlive() {
        return alive;
    }

};

#endif