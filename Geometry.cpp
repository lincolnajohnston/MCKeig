#ifndef GEOMETRY
#define GEOMETRY

#include "Cell.cpp"
#include <tuple>

class Geometry {
    private:
        std::vector<Cell *> cells;
        int groups;
        int delayed_groups = 6;  // TODO: get from the input file
    
    public:
        Geometry() {}

        Geometry(std::vector<Cell *> cells, int groups):cells(cells), groups(groups) {}

        void addCell(Cell *newCell) {
            cells.push_back(newCell);
        }

        Cell * cellAtLocation(Position p) {
            for (Cell *cell:cells) {
                if (cell->locationInCell(p)) {
                    return cell;
                }
            }
            return NULL;
        }

        void clearTallies() {
            for (Cell *cell:cells) {
                cell->clearTL();
            }
        }

        void printTallies() {
            for (Cell *cell:cells) {
                for (int g = 0; g < groups; g++) {
                    std::cout << cell->getName() << " group " << g << ": " << cell->getTLTally(g) << std::endl;
                }
            }
        }

        void printFluxes() {
            for (Cell *cell:cells) {
                for (int g = 0; g < groups; g++) {
                    std::cout << cell->getName() << " group " << g << ": " << cell->getTLTally(g)/cell->getVolume() << std::endl;
                }
            }
        }

        int getGroups() {
            return groups;
        }

        void createCellAdjointFluxes() {
            for (Cell *cell:cells) {
                cell->createAdjointFluxSolution();
            }
        }

        std::vector<double> getBetaEff() {
            std::vector<double> beta_sums(delayed_groups, 0);
            double fission_sum = 0;
            for (Cell *cell:cells) {
                for(int g = 0; g < delayed_groups; g++) {
                    beta_sums[g] += std::get<0>(cell->getBetaEff())[g];
                }
                fission_sum += std::get<1>(cell->getBetaEff());
            }

            std::vector<double> betas(delayed_groups, 0);
            for(int g = 0; g < delayed_groups; g++) {
                betas[g] = beta_sums[g] / fission_sum;
            }
            return betas;
        }


};

#endif