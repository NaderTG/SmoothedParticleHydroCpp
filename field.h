#ifndef FIELD_H
#define FIELD_H

#include <iostream>
#include <string>
#include "Particle.h"
#include "Cell.h"
#include "Vec.h"
#include "params.h"
#include <vector>
#include <fstream> // ifstream
#include <sstream> // stringstream

class Field{
public:
    //Implementation of periodic domain
    std::vector<Cell*> cells_domain;
    std::vector<Particle*> particle_list;
    
    int num_cells;
    int num_particles;
    int num_cells_1D;
    double width, height; //Total width of the field
    double dx, dy; //cell width
    double init_denisty, init_mass;
    int num_cells_hor, num_cells_ver; //number of cells of each diractions
    double _search_radius;
    std::string init_filename;
    std::string image_name;
    Field();
    Field(Params);
    
    ~Field(){};
    
    void initField();
    void Update_field();
    void loadParticles();
    
    void calcForce();
    void calcVelocity();
    void calcDensity();
    
};

Field::Field(){
    Params _params("test");
    num_cells = _params.num_cells;
    num_cells_hor = _params.num_cells_hor;
    num_cells_ver = _params.num_cells_ver;
    _search_radius = 0.1;
    
    init_mass = 1.0;
    init_denisty = 1.0;
    image_name = "test.pgm";
    width = _params.width;
    height = _params.height;
    dx = _params.dx;
    dy = _params.dx;
    loadParticles();
    intiField();
    
}
Field::Field(Params _params){
    
    num_cells = _params.num_cells;
    num_cells_hor = _params.num_cells_hor;
    num_cells_ver = _params.num_cells_ver;
    
    _search_radius = _params.search_radius;
    width = _params.width;
    height = _params.height;
    dx = _params.dx;
    dy = _params.dy;
    
    init_denisty = _params.init_denisty;
    init_mass = _params.init_mass;
    
    init_filename = _params.filename;
    image_name = _params.image_name;
    loadParticles();
    initField();
    
    
}

void Field::loadParticles(){
    
    
    //Read PGM file (needs rewriting!!)
    
    
    int xsize, ysize;
    
    double** _pos_list;
    std::ifstream infile(image_name);
    std::stringstream ss;
    std::string inputLine = "";
    // First line : version
    getline(infile,inputLine);
    if(inputLine.compare("P2") != 0) std::cout << "Version error" << std::endl;
    else std::cout << "Version : " << inputLine << std::endl;
    int _read_temp;
    // Second line : comment
    getline(infile,inputLine);
    std::cout << "Comment : " << inputLine << std::endl;
    
    // Continue with a stringstream
    ss << infile.rdbuf();
    // Third line : size
    int numcols, numrows;
    ss >> numcols >> numrows;
    std::cout << numcols << " columns and " << numrows << " rows" << std::endl;
    
    int counter = 0;
    ss >> _read_temp;
    // Following lines : data
    
    for(int row = 0; row < numrows; ++row){
        for (int col = 0; col < numcols; ++col) {
            ss >> _read_temp;
            if(_read_temp == 1){
                particle_list.push_back(new Particle(dx*col, dy*row, counter));
                particle_list.back()->mass = init_mass;
                particle_list.back()->density = init_denisty;
                counter++;
            }
        }
    }
    
    //Might put debug
    infile.close();
    
    
}


void Field::initField(){
    
    double x, y;
    Vec<2> _position(0.0);
    int x_idx, y_idx, idx_temp, idx;
    double dx_half = dx/ 2.0;
    double dy_half = dy / 2.0;
    
    for(int i = 0; i < num_cells_ver ; i++){ //Horizontal
        for(int j = 0; j < num_cells_hor ; j++){ //Vertical
            idx = num_cells_hor*i + j;
            
            _position[0] = i*dx+ dx_half;
            _position[1] = j*dy + dy_half;
            
            cells_domain.push_back(new Cell(_position, idx, 0, 0));
            
            for(int k1= 0; k1 < 3; k1++){
                for(int k2= 0; k2 < 3; k2++){
                    x_idx = (i + k1) % num_cells_hor;
                    y_idx = (j + k2) % num_cells_ver;
                    
                    idx_temp = num_cells_hor*x_idx + y_idx;
                    cells_domain[idx]->neighbour_index.push_back(idx_temp);
                    cells_domain[idx]->cell_width = dx;
                    cells_domain[idx]->cell_height = dy;
                }
            }
            cells_domain[idx]->neighbour_index.erase(cells_domain[idx]->neighbour_index.begin());
            
        }
    }
    
    //Assign particles to cells;
    for(int i = 0; i < num_particles; i++){
        _position = particle_list[i]->position;
        
        x_idx = (int) floor(_position[0]/ dx );
        y_idx = (int) floor(_position[1]/ dy );
        idx_temp = num_cells_hor*x_idx + y_idx;
        particle_list[i]->cell_ID = idx_temp;
        
        particle_list[i]->part_cell_order = cells_domain[idx_temp]->parts_idx.size();
        cells_domain[idx_temp]->parts_idx.push_back(i);
        
    }
    
}

void Field::calcDensity(){
    Vec<2> = _pos;
    Vec<2> = _pos_neighbour;
    int _idx_neighbour;
    for (int i = 0; i <  num_particles; i++){
        _pos = particle_list.at(i)->position;
        //particle_list.at(i)->density = particle_list.at(i)->mass ;
        
        for(int j = 0; j < particle_list.at(i)->neighbour_size(); j++){
            _idx_neighbour = particle_list.at(i)->neighbours_ID.at(j);
            _pos_neighbour = particle_list.at(_idx_neighbour)->position;
            
            // particle_list.at(i)->density += (particle_list.at(i)->mass)*W(_pos, _pos_neighbour);
            
        }
        
    }
    
}

void Field::calcForce(){
    
    calcDensity();
    
    Vec<2> = _pos;
    Vec<2> = _pos_neighbour;
    double _rho_i, _rho_j;
    int _idx_neighbour;
    
    for (int i = 0; i <  num_particles; i++){
        _pos = particle_list.at(i)->position;
       _rho_i = particle_list.at(i)->density ;
        
        for(int j = 0; j < particle_list.at(i)->neighbour_size(); j++){
            _idx_neighbour = particle_list.at(i)->neighbours_ID.at(j);
            _pos_neighbour = particle_list.at(_idx_neighbour)->position;
            
            // particle_list.at(i)->density += particle_list.at(i)->mass*W(_pos, _pos_neighbour);
            
        }
        
    }
    
    
    
}

void Field::Update_field(){
    
    Vec<2> _pos2(0.0); int _x_temp, _y_temp;
    Vec<2> _pos_other(0.0);
    int _idx_part,  _cell_idx , _cell_idx_neighbour;
    
    for(int i = 0; i < num_particles; i++){
        _cell_idx = particle_list.at(i)->cell_ID;
        _pos2 = particle_list.at(i)->position;
        //Clear neighbour list
        particle_list.at(i)->neighbours_ID.clear();
        for(int j = 0; j < cells_domain.at(_cell_idx)->getNumParts(); j++){
            //define dist function!!!!
            _idx_part = cells_domain.at(_cell_idx)->parts_idx.at(j);
            _pos_other =particle_list.at(_idx_part)->position;
            if(dist(_pos2, _pos_other , _search_radius)){
                particle_list.at(i)->neighbours_ID.push_back(_idx_part);
            }
            
        }
        //Loop on the neighbouring cells to find all particles in the radius;
        for(int k = 0; k < cells_domain.at(_cell_idx)->getNumNeighbours(); k++){
            _cell_idx_neighbour = cells_domain.at(_cell_idx)->neighbour_index[k];
            
            for(int j = 0; j < cells_domain.at(_cell_idx_neighbour)->getNumParts(); j++){
                //define dist function!!!!
                _idx_part = cells_domain.at(_cell_idx_neighbour)->parts_idx.at(j);
                _pos_other =particle_list.at(_cell_idx_neighbour)->position;
                if(dist(_pos2, _pos_other , _search_radius)){
                    particle_list.at(i)->neighbours_ID.push_back(_idx_part);
                }
                
            }
            
        }
        
        //erase the particle from the particle list in the cell
        cells_domain.at(_cell_idx)->parts_idx.erase(particle_list.at(i)->part_cell_order);
        
        _pos2 =particle_list.at(i)->position;
        _x_temp = (int) floor(_pos2[0] / dx);
        _y_temp = (int) floor(_pos2[1] / dy);
        _cell_idx = _x_temp*num_cells_hor + _y_temp;
        
        cells_domain.at(_cell_idx)->parts_idx.push_back(i);
        //Add the cell number to the particle
        particle_list.at(i)->cell_ID = _cell_idx;
        
    }
    
    
}
