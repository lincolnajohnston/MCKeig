#ifndef PARTICLE
#define PARTICLE

#include "Position.cpp"
#include "Cell.cpp"
#include "Geometry.cpp"
#include <queue>
#include <vector>
#include <memory>

class Particle {
    protected:
        Position location;
        Direction direction;
        bool alive;
        Geometry *geo;
        Cell *currentCell;
        double epsilon = pow(10,-10);
        double weight;
        int group;
        int delayed_group;
        int prompt_group;

    public:

    Particle(Position location, Direction direction, Geometry *geo, double weight = 1, int group = 0):location(location), direction(direction),
         geo(geo), weight(weight), group(group), delayed_group(1), prompt_group(0) {
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

    void move(std::vector<std::shared_ptr<Particle>> &generated_neutron_bank, std::vector<std::shared_ptr<Particle>> &delayed_neutron_bank, double &source_particles, double timestep, Rand& rng) {
        double edge_dist = currentCell->distToEdge(location, direction);
        double collision_dist = currentCell->distToNextCollision(rng, group);

        double move_dist = edge_dist + epsilon;
        if (collision_dist < move_dist) {
            move_dist = collision_dist;

            location.x += move_dist * direction.getI();
            location.y += move_dist * direction.getJ();
            location.z += move_dist * direction.getK();

            currentCell->tallyTL(move_dist, weight, group);

            std::string collision_name = currentCell->sample_collision(rng, group);
            if (collision_name == "scat") {
                direction.isotropicScatter(rng);
                group = currentCell->sample_scatter_group(rng, group);
            }
            else if (collision_name == "cap") {
                alive = false;
            }
            else if (collision_name == "census") {
                delayed_neutron_bank.push_back(std::make_shared<Particle>(*this));

                source_particles += this->getWeight();
                alive = false;
            }
            else if (collision_name == "fis") {
                for (size_t i = 0; i < currentCell->getNu(group); i++) {
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
                    double ksi = rng.getRand();
                    if (ksi < currentCell->getBeta()) { // delayed neutron
                        int del_group = currentCell->sampleDelayedGroup(rng);
                        std::shared_ptr<Particle> fission_source_del_neut = std::make_shared<Particle>(*this);
                        fission_source_del_neut->f3WeightAdjust(timestep,del_group); // hardcode group 1 as delayed neutron group
                        fission_source_del_neut->sampleDelayedNeutronEnergy();

                        std::shared_ptr<Particle> time_bank_del_neut = std::make_shared<Particle>(*this);
                        time_bank_del_neut->f2WeightAdjust(timestep, del_group); // hardcode group 1 as delayed neutron group
                        time_bank_del_neut->sampleDelayedNeutronEnergy();

                        generated_neutron_bank.push_back(fission_source_del_neut);
                        delayed_neutron_bank.push_back(time_bank_del_neut);

                        // tally that delayed fission took place with adjoint weighting, will need to update this when neutron group sampling is more complicated
                    }
                    else { // prompt neutron
                        std::shared_ptr<Particle> prompt_neut = std::make_shared<Particle>(*this);
                        prompt_neut->samplePromptNeutronEnergy();
                        generated_neutron_bank.push_back(prompt_neut);
                        // tally that prompt fission took place with adjoint weighting, will need to update this when neutron group sampling is more complicated
                    }
                    
                    source_particles += this->getWeight();
                }
                alive = false;
            }
            
        }
        else {
            location.x += move_dist * direction.getI();
            location.y += move_dist * direction.getJ();
            location.z += move_dist * direction.getK();

            currentCell->tallyTL(move_dist, weight, group);

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
    void f1WeightAdjust(double timestep, int delayed_group) {
        multiplyWeight(exp(-1 * currentCell->getDecayConst(delayed_group) * timestep));
    }

    // adjust weight by probability of precursor decaying within this time step and being being considered as a delayed neutron for the next time step
    void f2WeightAdjust(double timestep, int delayed_group) {
        double lambda = currentCell->getDecayConst(delayed_group);
        multiplyWeight((1 - exp(-1 * lambda * timestep) - lambda * timestep * exp(-1 * lambda * timestep)) / (lambda * timestep));
    }

    // adjust weight by probability of precursor decaying within this time step and being being considered as a delayed neutron for the current fission source iteration
    void f3WeightAdjust(double timestep, int delayed_group) {
        double lambda = currentCell->getDecayConst(delayed_group);
        multiplyWeight(exp(-1 * lambda * timestep) * (1 - exp(lambda * timestep) + lambda * timestep * exp(lambda * timestep)) / (lambda * timestep));
    }

    bool isAlive() {
        return alive;
    }

    // sample an energy group for this particle if it originates from a delayed neutron event
    void sampleDelayedNeutronEnergy() {
        // TODO: More accurately sample energy group when there are more than 2 groups
        group = delayed_group;
    }

    // sample an energy group for this particle if it originates from a prompt fission
    void samplePromptNeutronEnergy() {
        // TODO: More accurately sample energy group when there are more than 2 groups
        group = prompt_group;
    }

};

#endif