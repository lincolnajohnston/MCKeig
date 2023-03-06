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
            }
        }

        void printTallies() {
            for (Cell *cell:cells) {
                std::cout << cell->getName() << ": " << cell->getTLTally() << std::endl;
            }
        }

};

#endif