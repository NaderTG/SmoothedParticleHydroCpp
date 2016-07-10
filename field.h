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
    double w_h, h_h; //cell width
    int num_cells_w, num_cells_h; //number of cells of each diractions
    double _search_radius;
    std::string init_filename;
    Field();
    Field(Params);
    
    ~Field(){};
    
    void initField();
    void Update_field();
    void loadParticles();
    
    void calcForce();
    void calcVelocity();
    
};

Field::Field(){
    Params _params();
    num_cells = _params.num_cells;
    num_cells_w = (int) (sqrt(num_cells));
    num_cells_h = (int) (sqrt(num_cells));
    
    num_particles = _params.num_particles;
    
    init_filename = "default.pgm";
    width = _params.width;
    w_h = _params.dx;
    h_h = _params.dx;
    loadParticles();
    intiField();
    
}
Field::Field(Params _params){
    
    num_cells = _params.num_cells;
    num_particles = _params.num_particles;
    width = _params.width;
    w_h = _params.dx;
    h_h = _params.dx;
    
    init_filename = _params.filename;
    
    loadParticles();
    initField();
    
    
}

void Field::loadParticles(){
    
    
    //Read PGM file (needs rewriting!!)
    
    int xsize, ysize;
    double dx, dy;
    ifstream infile(init_filename);
    stringstream ss;
    string inputLine = "";
    // First line : version
    getline(infile,inputLine);
    if(inputLine.compare("P2") != 0) cerr << "Version error" << endl;
    else cout << "Version : " << inputLine << endl;
    int _read_temp;
    // Second line : comment
    getline(infile,inputLine);
    cout << "Comment : " << inputLine << endl;
    
    // Continue with a stringstream
    ss << infile.rdbuf();
    // Third line : size
    ss >> numcols >> numrows;
    cout << numcols << " columns and " << numrows << " rows" << endl;
    int counter = 0;
    // Following lines : data
    for(row = 0; row < numrows; ++row)
        for (col = 0; col < numcols; ++col) {
            ss >> _read_temp;
            if(_read_temp == 1){
                particle_list.push_back(new Particle(dx*col, dy*row, counter));
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
    double w_half = width/ 2.0;
    double h_half = height / 2.0;
    
    for(int i = 0; i < num_cells_h ; i++){ //Horizontal
        for(int j = 0; j < num_cells_w ; j++){ //Vertical
            idx = num_cells_w*i + j;
            
            _position[0] = i*width + w_half;
            _position[1] = j*height + h_half;
            
            cells_domain.push_back(new Cell(_position, idx, 0, 0));
            
            for(int k1= 0; k1 < 3; k1++){
                for(int k2= 0; k2 < 3; k2++){
                    x_idx = (i + k1)% num_cells_h;
                    y_idx = (j + k2) % num_cells_w;
                    
                    idx_temp = num_cells_w*x_idx + y_idx;
                    cells_domain[idx]->neighbour_index.push_back(idx_temp);
                    
                }
            }
            cells_domain[idx]->neighbour_index.erase(cells_domain[idx]->neighbour_index.begin());
            
        }
    }
    
    //Assign particles to cells;
    for(int i = 0; i < num_particles; i++){
        _position = particle_list[i]->position;
        
        x_idx = (int) floor(_position[0]/ width );
        y_idx = (int) floor(_position[1]/ width );
        idx_temp = num_cells_w*x_idx + y_idx;
        particle_list[i]->cell_ID = idx_temp;
        
        particle_list[i]->part_cell_order = cells_domain[idx_temp]->parts_idx.size();
        cells_domain[idx_temp]->parts_idx.push_back(i);
        
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
        _x_temp = (int) floor(_pos2[0] / width);
        _y_temp = (int) floor(_pos2[1] / height);
        _cell_idx = _x_temp*num_cells_w + _y_temp;
        
        cells_domain.at(_cell_idx)->parts_idx.push_back(i);
        //Add the cell number to the particle
        particle_list.at(i)->cell_ID = _cell_idx;
        
    }
    
    
}
