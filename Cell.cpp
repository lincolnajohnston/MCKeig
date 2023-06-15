#ifndef CELL
#define CELL

#include "Position.cpp"
#include "Surface.cpp"
#include "Rand.cpp"
#include "Material.cpp"

class Cell {
    protected:

        std::vector<Surface *> surfaces;
        std::vector<bool> senses;
        std::string name;
        int groups;
        Material *mat;
        std::vector<double> TLtally;
        double volume;

    public:

    Cell(std::string name):name(name) {

    }

    Cell(std::string name, std::vector<Surface *> surfs, std::vector<bool> senses, int groups, Material *mat, double volume = 0):name(name),
        surfaces(surfs), senses(senses), groups(groups), mat(mat), TLtally(groups, 0), volume(volume) {
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

    double getNuSigmaF(int group) {
        return mat->getNu(group) * mat->getFissionXS(group);
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

    double getVolume() {
        return volume;
    }


};

#endif