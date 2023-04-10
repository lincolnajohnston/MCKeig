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
        std::vector<std::vector<double>> adjoint_flux;
        std::vector<Surface *> surfaces;
        std::vector<std::string> surface_names;
        std::vector<Material *> materials;
        std::vector<std::string> material_names;
        std::vector<Cell *> cells;
        double deltaT;
        int groups;
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
        std::string group_str, totalXS_str, scatterXS_str, capXS_str, fisXS_str, nu_str, beta_str, decayConst_str;

        getline(myfile, group_str);
        int groups = std::stoi(group_str);

        getline(myfile, totalXS_str);
        getline(myfile, scatterXS_str);
        getline(myfile, capXS_str);
        getline(myfile, fisXS_str);
        getline(myfile, nu_str);
        getline(myfile, beta_str);
        getline(myfile, decayConst_str);
        std::cout << totalXS_str << std::endl;
        std::cout << decayConst_str << std::endl;

        std::stringstream totalXS_line(totalXS_str);
        std::stringstream scatterXS_line(scatterXS_str);
        std::stringstream capXS_line(capXS_str);
        std::stringstream fisXS_line(fisXS_str);
        std::stringstream nu_line(nu_str);
        std::stringstream beta_line(beta_str);
        std::stringstream decayConst_line(decayConst_str);

        std::string totalXS_str_temp, scatterXS_str_temp, capXS_str_temp, fisXS_str_temp, nu_str_temp, beta_str_temp, decayConst_str_temp;
        std::vector<double> totalXS, scatterXS, capXS, censusXS, fisXS, nu , beta, decayConst;
        std::vector< std::vector<double> > scatter_matrix_XS;
        scatter_matrix_XS.resize(groups);
        double v = 10000;  // HARDCODED HERE, put delta_t in material class, energy in particle class
        for (int i = 0; i < groups; i++) {
            getline(totalXS_line, totalXS_str_temp, ' ');
            totalXS.push_back(std::stod(totalXS_str_temp));

            getline(scatterXS_line, scatterXS_str_temp, ' ');
            scatterXS.push_back(std::stod(scatterXS_str_temp));

            getline(capXS_line, capXS_str_temp, ' ');
            capXS.push_back(std::stod(capXS_str_temp));

            censusXS.push_back(1/v);

            getline(fisXS_line, fisXS_str_temp, ' ');
            fisXS.push_back(std::stod(fisXS_str_temp));

            getline(nu_line, nu_str_temp, ' ');
            nu.push_back(std::stod(nu_str_temp));

            getline(beta_line, beta_str_temp, ' ');
            beta.push_back(std::stod(beta_str_temp));

            getline(decayConst_line, decayConst_str_temp, ' ');
            decayConst.push_back(std::stod(decayConst_str_temp));

            std::string scatter_matrix_str;
            std::string scatter_matrix_str_temp;
            getline(myfile, scatter_matrix_str);
            std::stringstream scatter_matrix_line(scatter_matrix_str);
            for (int j = 0; j < groups; j++) {
                getline(scatter_matrix_line, scatter_matrix_str_temp, ' ');
                scatter_matrix_XS[i].push_back(std::stod(scatter_matrix_str_temp));
            }
        }


        Material *mat = new Material(groups, totalXS, scatterXS, scatter_matrix_XS, capXS, censusXS, fisXS, nu , beta, decayConst);
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
            if (surf_name.find_first_not_of("0123456789.") == std::string::npos) {
                break;
            }

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
        double volume = std::stod(surf_name);
        cells.push_back(new Cell(cell_name, cellSurfaces, senses, groups, mat, adjoint_flux[cells.size()], volume));
    }

    void readAdjointFlux(std::string adjoint_flux_filename) {
        std::ifstream myfile; 
        myfile.open(adjoint_flux_filename);
        std::string flux_line;
        std::string groups_string;

        getline(myfile, groups_string);
        int groups = std::stoi(groups_string);

        std::string adjoint_flux_str;
        int cell_num = 0;
        while(getline(myfile, flux_line)) {
            std::stringstream line(flux_line);
            adjoint_flux.resize(cell_num + 1);
            for (int j = 0; j < groups; j++) {
                getline(line, adjoint_flux_str, ' ');
                adjoint_flux[cell_num].push_back(std::stod(adjoint_flux_str));
            }
            cell_num++;
        }

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
                if (file_string.rfind("%", 0) == 0) {
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
                            std::string identifier;
                            getline(line, identifier, ' ');
                            if (identifier == "TFS") {
                                simType == identifier;
                                std::string delta_t_string;
                                getline(line, delta_t_string, ' ');
                                deltaT = std::stod(delta_t_string);
                                std::string group_string;
                                getline(line, group_string, ' ');
                                groups = std::stoi(group_string);
                            }
                            else if (identifier == "adjoint_flux") {
                                std::string adjoint_flux_filename;
                                getline(line, adjoint_flux_filename, ' ');
                                readAdjointFlux(adjoint_flux_filename);
                            }
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
        return new Geometry(cells, groups);
    }
};

#endif