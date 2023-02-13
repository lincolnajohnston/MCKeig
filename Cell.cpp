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
        Material *mat;
        std::vector<std::string> xs_names;

    public:

    Cell(std::string name):name(name) {

    }

    Cell(std::string name, std::vector<Surface *> surfs, std::vector<bool> senses, Material *mat):name(name),
        surfaces(surfs), senses(senses), mat(mat), xs_names({"total", "scat", "cap", "fis", "nu"}) {
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

    double getNu() {
        return mat->getNu();
    }

    double getBeta() {
        return mat->getBeta();
    }

    double getDecayConst() {
        return mat->getDecayConst();
    }

    double distToNextCollision() {
        double rand_val = ((double)rand())/(double)RAND_MAX;
        return -1 * log(rand_val)/mat->getTotalXS();
    }

    std::string sample_collision() {
        double ksi = Rand::getRand();

        if (ksi < mat->getScatterXS() / mat->getTotalXS()) {
            return "scat";
        }
        else if (ksi < (mat->getScatterXS() + mat->getCaptureXS()) / mat->getTotalXS()) {
            return "cap";
        }
        else if (ksi < (mat->getScatterXS() + mat->getCaptureXS() + mat->getCensusXS()) / mat->getTotalXS()) {
            return "census";
        }
        else {
            return "fis";
        }
    }

};

#endif