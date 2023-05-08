#ifndef PROGENITOR
#define PROGENITOR

#include "Particle.cpp"
#include <memory>


class Progenitor {
    protected:
        double weight;
        double travel_time; // travel time before fissioning
        int delayed_group; // delayed group progenitor created when it fissioned (-1 if prompt)
        int timestep; // timestep that this progenitor neutron was created
        std::shared_ptr<Progenitor> prog;

    public:

    Progenitor(double weight, double travel_time, int delayed_group, int timestep, std::shared_ptr<Progenitor> prog): weight(weight), 
        travel_time(travel_time), delayed_group(delayed_group), timestep(timestep), prog(prog) { }

    std::shared_ptr<Progenitor> getProg() {
        return prog;
    }

    double getTravelTime() {
        return travel_time;
    }

    int getDelayedGroup() {
        return delayed_group;
    }

    int getTimestep() {
        return timestep;
    }


};

#endif