#ifndef GEOMETRY
#define GEOMETRY

#include "Cell.cpp"

class Geometry {
    private:
        std::vector<Cell *> cells;
    
    public:
        Geometry() {}

        Geometry(std::vector<Cell *> cells):cells(cells) {}

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
                std::cout << cell->getName() << ": " << cell->getTLTally() << std::endl;
            }
        }

        void printFluxes() {
            for (Cell *cell:cells) {
                std::cout << cell->getName() << " Flux (scaled by arbitrary constant): " << cell->getTLTally() / cell->getVolume() << std::endl;
            }
        }

        double getAverageBeta() {
            double adjoint_flux_total = 0;
            double adjoint_beta_total = 0;
            for (Cell *cell:cells) {
                cell->addBetaTallies(adjoint_flux_total, adjoint_beta_total);
            }
            return adjoint_beta_total / adjoint_flux_total;
        }


};

#endif