#ifndef CELL
#define CELL

#include "Position.cpp"
#include "Surface.cpp"
#include "Rand.cpp"
#include "Material.cpp"
#include <tuple>

class Cell {
    protected:

        std::vector<Surface *> surfaces;
        std::vector<bool> senses;
        std::string name;
        int groups;
        Material *mat;
        std::vector<double> TLtally;
        std::vector<double> adjoint_flux;
        double volume;

        std::vector<double> beta_sums;
        double fission_sum;

        //temporary vectors used in the creation of adjoint flux
        std::vector<double> adjoint_weighted_fission_source;
        std::vector<double> fission_source;

    public:

    Cell(std::string name):name(name) {

    }

    Cell(std::string name, std::vector<Surface *> surfs, std::vector<bool> senses, int groups, Material *mat, double volume = 0):name(name),
        surfaces(surfs), senses(senses), groups(groups), mat(mat), TLtally(groups, 0), adjoint_flux(groups, 0),volume(volume), 
        beta_sums(mat->getNumDelayedGroups(), 0), fission_sum(0), adjoint_weighted_fission_source(groups, 0), fission_source(groups, 0) {
        if (surfs.size() != senses.size()) {
            std::cout << "Error: Senses not defined for each surface" << std::endl;
        }
    }

    void addSurface(Surface  *surf, double sense) {
        surfaces.push_back(surf);
        senses.push_back(sense);
    }

    bool locationInCell(Position p) {
        for (size_t i = 0; i < surfaces.size(); i++) {
            
            if (surfaces[i]->senseOfPosition(p) != senses[i]) {
                return false;
            }
        }
        return true;
    }

    double distToEdge(Position p, Direction d) {
        double min_dist = std::numeric_limits<double>::max();
        for (Surface *surf:surfaces) {
            double cur_dist = surf->distanceToSurface(p,d);
            if (cur_dist < min_dist) {
                min_dist = cur_dist;
            }
        }
        return min_dist;
    }

    std::string getName() {
        return name;
    }

    double getVelocity(int group) {
        return mat->getVelocity(group);
    }

    double getNu(int group) {
        return mat->getNu(group);
    }

    double getBeta() {
        return mat->getBeta();
    }

    double getBeta(int del_group) {
        return mat->getBeta(del_group);
    }

    double getDecayConst(int del_group) {
        return mat->getDecayConst(del_group);
    }

    int sampleDelayedGroup(Rand &rng) {
        return mat->sampleDelayedGroup(rng);
    }

    double distToNextCollision(Rand &rng, int group) {
        double rand_val = rng.getRand();
        return -1 * log(rand_val)/mat->getTotalXS(group);
    }

    std::string sample_collision(Rand &rng, int group) {
        double ksi = rng.getRand();

        if (ksi < mat->getScatterXS(group) / mat->getTotalXS(group)) {
            return "scat";
        }
        else if (ksi < (mat->getScatterXS(group) + mat->getCaptureXS(group)) / mat->getTotalXS(group)) {
            return "cap";
        }
        else if (ksi < (mat->getScatterXS(group) + mat->getCaptureXS(group) + mat->getCensusXS(group)) / mat->getTotalXS(group)) {
            return "census";
        }
        else {
            return "fis";
        }
    }

    int sample_scatter_group(Rand &rng, int in_group) {
        std::vector<double> scatter_mat = mat->getScatterMatrixXS(in_group);
        double total_scatter_xs = mat->getScatterXS(in_group);

        double ksi = rng.getRand() * total_scatter_xs;
        double scatterXSsum = 0;

        for(int i = 0; i < scatter_mat.size(); i++) {
            scatterXSsum += scatter_mat[i];
            if (ksi < scatterXSsum) {
                return i;
            }
        }
    }

    void tallyTL(double dist, double weight, int group) {
        TLtally[group] += dist * weight;
    }

    double getTLTally(int group) {
        return TLtally[group];
    }

    void clearTL() {
        for (int g = 0; g < groups; g++) {
            TLtally[g] = 0;
        }
    }

    // record every time fission takes place, tally if it was a delayed neutron, weighted by importance
    void tallyBeta(double weight, int del_group, int E_group) { // energy group of outgoing neutron after fission used for importance
        if (del_group >= 0) {
            beta_sums[del_group] += weight * adjoint_flux[E_group];
        }
        fission_sum += weight * adjoint_flux[E_group];
    }

    double getVolume() {
        return volume;
    }

    void tallyAdjointFlux(double weight, double current_E_group, int prog_E_group) {
        adjoint_weighted_fission_source[prog_E_group] += weight;
        fission_source[current_E_group] += weight;
    }

    void createAdjointFluxSolution() {
        for(int i = 0; i < adjoint_flux.size(); i++) {
            if (fission_source[i] == 0) { // if no data collected on importance of flux in this cell at this energy, set adjoint to 1
                adjoint_flux[i] = 1;
            }
            else {
                adjoint_flux[i] = adjoint_weighted_fission_source[i] / fission_source[i];
            }
        }
    }

    double getAdjointFlux(int energy_group) {
        return adjoint_weighted_fission_source[energy_group];
    }

    std::tuple<std::vector<double>, double> getBetaEff() {
        return std::make_tuple(beta_sums, fission_sum);
    }


};

#endif