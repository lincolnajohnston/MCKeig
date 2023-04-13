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


};

#endif