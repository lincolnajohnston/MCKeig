#ifndef PROGENITOR
#define PROGENITOR

#include "Particle.cpp"
#include <memory>


class Progenitor {
    protected:
        double weight;
        double travel_time; // travel time before fissioning
        int energy_group; // energy group progenitor inhabited immediately following its creation
        Cell *source_cell;
        std::shared_ptr<Progenitor> prog;

    public:

    Progenitor(double weight, double travel_time, int energy_group, Cell* source_cell, std::shared_ptr<Progenitor> prog): weight(weight), 
        travel_time(travel_time), energy_group(energy_group), source_cell(source_cell), prog(prog) { }

    std::shared_ptr<Progenitor> getProg() {
        return prog;
    }

    double getTravelTime() {
        return travel_time;
    }

    int getEnergyGroup() {
        return energy_group;
    }

    Cell *getCell() {
        return source_cell;
    }

    double getWeight() {
        return weight;
    }


};

#endif