#ifndef RAND
#define RAND

#include <random>

class Rand {
    public:
        std::mt19937_64 rng;
        std::uniform_real_distribution<double> unif;
        Rand():unif(0,1) {
            rng.seed(12345);
        }

    double getRand() {
        return (double)unif(rng);
    }
};

#endif