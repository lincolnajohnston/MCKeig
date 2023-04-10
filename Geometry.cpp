#ifndef GEOMETRY
#define GEOMETRY

#include "Cell.cpp"

class Geometry {
    private:
        std::vector<Cell *> cells;
        int groups;
    
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
                cell->resetBetas();
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

        void getAverageBetaAndLambda(double &beta_eff, double &lambda_eff) {
            double adjoint_flux_total = 0;
            double adjoint_beta_total = 0;
            double adjoint_lambda_total = 0;
            for (Cell *cell:cells) {
                cell->addBetaTallies(adjoint_flux_total, adjoint_beta_total, adjoint_lambda_total);
            }
            beta_eff = adjoint_beta_total / adjoint_flux_total;
            lambda_eff = adjoint_lambda_total / adjoint_flux_total;
        }

        int getGroups() {
            return groups;
        }


};

#endif