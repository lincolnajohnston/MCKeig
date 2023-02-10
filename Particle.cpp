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

    void multiplyWeight(double factor) {
        weight = weight * factor;
    }

    void move(std::vector<Particle> &generated_neutron_bank, std::vector<Particle> &delayed_neutron_bank, int &source_particles) {
        double edge_dist = currentCell->distToEdge(location, direction);
        double collision_dist = currentCell->distToNextCollision();
        //std::cout << "Dist To Edge: " << edge_dist << std::endl;
        //std::cout << "Dist To Collision: " << collision_dist << std::endl;

        double move_dist = edge_dist + epsilon;
        if (collision_dist < edge_dist) {
            move_dist = collision_dist + epsilon;

            location.x += move_dist * direction.getI();
            location.y += move_dist * direction.getJ();
            location.z += move_dist * direction.getK();


            std::string collision_name = currentCell->sample_collision();
            if (collision_name == "scat") {
                //std::cout << "Particle scattered at " << location.x << ", " << location.y << ", " << location.z << std::endl;
                direction.isotropicScatter();
            }
            else if (collision_name == "cap") {
                //std::cout << "Particle captured at " << location.x << ", " << location.y << ", " << location.z << std::endl;
                alive = false;
            }
            else if (collision_name == "census") {
                //std::cout << "Particle censused at " << location.x << ", " << location.y << ", " << location.z << std::endl;
                delayed_neutron_bank.push_back(*this);

                source_particles += 1;
                alive = false;
            }
            else if (collision_name == "fis") {
                //std::cout << "Particle fissioned at " << location.x << ", " << location.y << ", " << location.z << std::endl;
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
                    direction.isotropicScatter();
                    double ksi = Rand::getRand();
                    if (ksi < currentCell->getBeta()) {
                        delayed_neutron_bank.push_back(*this);
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

            currentCell = geo->cellAtLocation(location);

            if (!currentCell) {
                //std::cout << "Particle leaked at position: " << location.x << ", " << location.y << ", " << location.z << std::endl;
                alive = false;
            }
            else {
                //std::cout << "New Cell: " << currentCell->getName() << std::endl;
            }
        }
        
    }


    

    bool isAlive() {
        return alive;
    }

};

#endif