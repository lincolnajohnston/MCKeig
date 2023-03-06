#ifndef POSITION
#define POSITION

#include <iostream>
#include <vector>
#include "math.h"
#include "Rand.cpp"

struct Position {
    double x;
    double y;
    double z;

    std::vector<double> getVector() {
        return {x,y,z};
    }

    double getX() {
        return x;
    }

    double getY() {
        return y;
    }

    double getZ() {
        return z;
    }

};

class Direction {
    public:
        double i;
        double j;
        double k;
        double magnitude;
    
    public:
        Direction(double i_in, double j_in, double k_in) {
            double norm_const = sqrt(i_in*i_in + j_in*j_in + k_in*k_in);
            i = i_in / norm_const;
            j = j_in / norm_const;
            k = k_in / norm_const;
        }

    std::vector<double> getDirection() {
        return {i,j,k};
    }

    double getI() {
        return i;
    }

    double getJ() {
        return j;
    }

    double getK() {
        return k;
    }

    void isotropicScatter(Rand &rng) {
        double mu = 2 * rng.getRand2() - 1;
        double phi = 2 * M_PI * rng.getRand2();

        i = mu;
        j = sqrt(1 - pow(mu,2)) * sin(phi);
        k = sqrt(1 - pow(mu,2)) * cos(phi);

    }
};

#endif