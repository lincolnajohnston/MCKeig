#ifndef SOURCE
#define SOURCE

#include <vector>
#include "Particle.cpp"

class Source {
    private:
        double source_rate;
         Geometry *geo;
    public:
    Source(double source_rate, Geometry *geo):source_rate(source_rate), geo(geo) {

    }

    void generateParticles(std::vector<Particle> &part_bank) {
        for (int i = 0; i < source_rate; i++) {
            part_bank.push_back(Particle({0,0,0}, {1,0,0}, geo));
        }
    }
};

#endif