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
        std::vector<double> adjoint_flux;
        double volume;
        double adjoint_flux_sum;
        double adjoint_beta_sum;
        double adjoint_lambda_sum;

    public:

    Cell(std::string name):name(name) {

    }

    Cell(std::string name, std::vector<Surface *> surfs, std::vector<bool> senses, int groups, Material *mat, std::vector<double> adjoint_flux, double volume = 0):name(name),
        surfaces(surfs), senses(senses), groups(groups), mat(mat), TLtally(groups, 0), adjoint_flux(adjoint_flux), volume(volume) {
        if (surfs.size() != senses.size()) {
            std::cout << "Error: Senses not defined for each surface" << std::endl;
        }
        adjoint_flux_sum = 0;
        adjoint_beta_sum = 0;
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

    double getNu(int group) {
        return mat->getNu(group);
    }

    double getBeta(int group) {
        return mat->getBeta(group);
    }

    double getDecayConst(int group) {
        return mat->getDecayConst(group);
    }

    double distToNextCollision(int group) {
        double rand_val = ((double)rand())/(double)RAND_MAX;
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

    double getAdjointFlux(double g) {
        return adjoint_flux[g];
    }

    void tallyBeta(double beta_tally, double g) {
        adjoint_flux_sum += adjoint_flux[g];
        adjoint_beta_sum += beta_tally * adjoint_flux[g];
        adjoint_lambda_sum += getDecayConst(g) * adjoint_flux[g];
    }

    void resetBetas() {
        adjoint_flux_sum = 0;
        adjoint_beta_sum = 0;
        adjoint_lambda_sum = 0;
    }

    void addBetaTallies(double &adjoint_flux_outer_sum, double &adjoint_beta_outer_sum, double &adjoint_lambda_outer_sum) {
        adjoint_flux_outer_sum += adjoint_flux_sum;
        adjoint_beta_outer_sum += adjoint_beta_sum;
        adjoint_lambda_outer_sum += adjoint_lambda_sum;
    }

};

#endif