#ifndef INPUT
#define INPUT

#include <iostream>
#include <fstream>
#include <string>
#include <bits/stdc++.h>
#include "Geometry.cpp"

class Input {
    public:
        std::string input_file_name;
        std::vector<Surface *> surfaces;
        std::vector<std::string> surface_names;
        std::vector<Material *> materials;
        std::vector<std::string> material_names;
        std::vector<Cell *> cells;
        double deltaT;
        std::string simType;

        Input(std::string input_file_name):input_file_name(input_file_name) { }

    // read one singular line containing information about a surface, add it to the vector of surfaces
    void readSurface(std::string line_string) {
        std::stringstream line(line_string);
        std::string surf_name;
        getline(line, surf_name, ' ');
        surface_names.push_back(surf_name);

        std::string surf_type;
        getline(line, surf_type, ' ');
        if (surf_type == "Plane") {
            std::string x0,y0,z0,u,v,w;
            getline(line, x0, ' ');
            getline(line, y0, ' ');
            getline(line, z0, ' ');
            getline(line, u, ' ');
            getline(line, v, ' ');
            getline(line, w, ' ');

            Plane *plane = new Plane({std::stod(x0),std::stod(y0),std::stod(z0)}, {std::stod(u),std::stod(v),std::stod(w)});
            surfaces.push_back(plane);
        }
        else if (surf_type == "Sphere") {
            std::string x0,y0,z0,r;
            getline(line, x0, ' ');
            getline(line, y0, ' ');
            getline(line, z0, ' ');
            getline(line, r, ' ');

            Sphere *sphere = new Sphere({std::stod(x0),std::stod(y0),std::stod(z0)}, std::stod(r));
            surfaces.push_back(sphere);
        }
    }

    // read one singular line containing information about a material, add it to the vector of materials
    void readMaterial(std::string line_string) {
        std::stringstream line(line_string);
        std::string mat_name, mat_file_name;
        getline(line, mat_name, ' ');
        getline(line, mat_file_name, ' ');
        material_names.push_back(mat_name);

        std::ifstream myfile; 
        myfile.open(mat_file_name);
        std::string totalXS, scatterXS, capXS, fisXS, nu, beta, decayConst;
        getline(myfile, totalXS);
        getline(myfile, scatterXS);
        getline(myfile, capXS);
        getline(myfile, fisXS);
        getline(myfile, nu);
        getline(myfile, beta);
        getline(myfile, decayConst);

        double v = 10000;  // HARDCODED HERE, PUT delta_t in material class, energy in particle class

        Material *mat = new Material(std::stod(totalXS), std::stod(scatterXS), std::stod(capXS), 1 / (deltaT * v), std::stod(fisXS), std::stod(nu), std::stod(beta), std::stod(decayConst));
        materials.push_back(mat);
    }

    // read one singular line containing information about a cell, add it to the vector of cells
    void readCell(std::string line_string) {
        std::stringstream line(line_string);
        std::string cell_name;
        getline(line, cell_name, ' ');

        std::string mat_name;
        getline(line, mat_name, ' ');
        int mat_index = -1;
        for(int i = 0; i < material_names.size(); i++) {
            if (material_names[i] == mat_name) {
                mat_index = i;
                break;
            }
        }
        Material *mat = materials[mat_index];

        std::string surf_name;
        std::vector<Surface *> cellSurfaces;
        std::vector<bool> senses;
        int surf_index = -1;
        while(getline(line, surf_name, ' ')) {
            int sense;
            if (surf_name.at(0) == '-') {
                senses.push_back(false);
                surf_name = surf_name.substr(1);
            }
            else {
                senses.push_back(true);
            }
            for(int i = 0; i < surface_names.size(); i++) {
                if (surface_names[i] == surf_name) {
                    surf_index = i;
                    break;
                }
            }
            cellSurfaces.push_back(surfaces[surf_index]);
        }
        cells.push_back(new Cell(cell_name, cellSurfaces, senses, mat));
    }

    // create and return geometry using the input file specified in the creation of the Input object
    Geometry *getGeometry() {
        std::ifstream myfile; 
        myfile.open(input_file_name);
        std::string file_string;
        int input_phase = 0;
        if ( myfile.is_open() ) {
            while(std::getline(myfile, file_string)) {
                std::cout << "line: " << file_string << std::endl;
                if (file_string == "") {
                    continue;
                }
                if (file_string.rfind("!", 0) == 0) {
                    if (file_string.rfind("!Simulation", 0) == 0) {
                        input_phase = 1;
                    }
                    else if (file_string.rfind("!Materials", 0) == 0) {
                        input_phase = 2;
                    }
                    else if (file_string.rfind("!Surfaces", 0) == 0) {
                        input_phase = 3;
                    }
                    else if (file_string.rfind("!Cells", 0) == 0) {
                        input_phase = 4;
                    }
                    else if (file_string.rfind("!Tallies", 0) == 0) {
                        input_phase = 5;
                    }
                    continue;
                }
                switch (input_phase) {
                    case 0:
                        {
                            std::cout << "phase 0" << std::endl;
                            break;
                        }
                    case 1:
                        {
                            std::cout << "Simulation phase" << std::endl;
                            std::cout << file_string << std::endl;
                            std::stringstream line(file_string);
                            getline(line, simType, ' ');
                            std::string delta_t_string;
                            getline(line, delta_t_string, ' ');
                            deltaT = std::stod(delta_t_string);
                            break;
                        }
                    case 2:
                        {
                            std::cout << "Materials phase" << std::endl;
                            readMaterial(file_string);
                            break;
                        }
                    case 3:
                        {
                            std::cout << "Surfaces phase" << std::endl;
                            readSurface(file_string);
                            break;
                        }
                    case 4:
                        {
                            std::cout << "Cells Phase" << std::endl;
                            readCell(file_string);
                            break;
                        }
                    case 5:
                        {
                            std::cout << "Tallies Phase" << std::endl;
                            // TODO: implement different types of tallies and method of tracking/displaying them
                            break;
                        }
                }
            }
        }
        return new Geometry(cells);
    }
};

#endif