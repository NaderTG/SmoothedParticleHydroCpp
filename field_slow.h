//
//  field_slow.h
//  SPH3
//
//  Created by Nader on 21/08/2016.
//  Copyright (c) 2016 Nader. All rights reserved.
//

#ifndef SPH3_field_slow_h
#define SPH3_field_slow_h


#include <iostream>
#include <string>
#include "Particle.h"
#include "Cell.h"
#include "Vec.h"
#include "params.h"
#include "Kernel.h"
#include <vector>
#include <fstream> // ifstream
#include <sstream> // stringstream


bool dist(Vec<2>& _pos_a, Vec<2>& _pos_b , double radius){
    bool result = false;
    Vec<2> _r = _pos_a - _pos_b;
    double _r_norm = _r.l2norm();
    
    result = (_r_norm < radius);
    
    return result;
}

class Field_slow{
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
    double pressure_denisty;
    
    double nu;
    double lambda;
    
    std::string init_filename;
    std::string image_name;
    
    ker_spline W_fun;
    
    Field_slow();
    Field_slow(Params&);
    
    ~Field_slow(){};
    
    friend std::ostream& operator<<(std::ostream& output, Field_slow& _A){
        for(int i = 0; i < _A.num_particles; i++){
             output << _A.particle_list.at(i)->position[0] << "\t" << _A.particle_list.at(i)->position[1] << "\t" << _A.particle_list.at(i)->density <<  "\t" << "[" <<  _A.particle_list.at(i)->force[0] << ", " <<  _A.particle_list.at(i)->force[1] << "]\t" <<  _A.particle_list.at(i)->pressure << "\n";
            
//                        output << _A.particle_list.at(i)->position[0] << "\t" << _A.particle_list.at(i)->position[1] << "\t" << _A.particle_list.at(i)->density <<  "\n" ;
            
        }
        return output;
    }
    
    
    void initField();
    void Update_field();
    void loadParticles();
    
    void loadTest();
    void neighbourListing();
    
    void calcForce();
    void calcVelocity(double);
    void calcInitVelocity(double);
    void calcPressure();
    void calcDensity();
    
};

Field_slow::Field_slow(){
    Params _params("test");
    num_cells = _params.num_cells;
    num_cells_hor = _params.num_cells_hor;
    num_cells_ver = _params.num_cells_ver;
    _search_radius = 0.2;
    pressure_denisty = _params.pressure;
    W_fun = *new ker_spline();
    
    init_mass = 1.0;
    init_denisty = 1.0;
    image_name = "toystar.pgm";
    width = _params.width;
    height = _params.height;
    dx = _params.dx;
    dy = _params.dx;
   loadParticles();
    //loadTest();
    initField();
    
}
Field_slow::Field_slow(Params& _params){
    
    num_cells = _params.num_cells;
    num_cells_hor = _params.num_cells_hor;
    num_cells_ver = _params.num_cells_ver;
    
    _search_radius = _params.search_radius;
    width = _params.width;
    height = _params.height;
    dx = _params.dx;
    dy = _params.dy;
    
    init_denisty = _params.init_density ;
    init_mass = _params.init_mass;
    pressure_denisty = _params.pressure;
    nu =  1.0;
    
    W_fun = *new ker_spline();
    
    // W_fun.width = _params.kernel_radius;
    std::cout << "Kernel radius = " << _params.kernel_radius << std::endl;
    W_fun.changeWidth(_params.kernel_radius);
    //    lambda = 0.1; //Needs calculating
    
    double temp = (_params.star_mass*2.0 /(_params.star_rad * _params.star_rad));
    lambda = 2.0*pressure_denisty*(M_1_PI)*(temp*temp) /_params.star_mass ;
    
    std::cout << "Lambda = " << lambda << std::endl;
    init_filename = _params.filename;
    image_name = _params.image_name;
    loadParticles();
    //loadTest();
    initField();
    
    
}

void Field_slow::loadTest(){
    std::cout << "Loading!\n";
    int num_part = 20;
    
    double radius = 0.25; double _x, _y;
    
    double theta = (2.0*M_PI)/ (double)  (num_part -1);
    
    
    //std::cout << "theta = " << theta << std::endl;
    
    for(int i = 0; i <  num_part ; i++){
        radius = (double)(i + 1)/(double) num_part;
        _x = radius*(sin((i)*theta)) ;
        _y = radius*(cos(i*theta)) ;
        //    std::cout << "x = [" << _x << " , " << _y << "]\n";
        particle_list.push_back(new Particle(_x, _y, i));
        particle_list.back()->mass = init_mass;
        particle_list.back()->density = init_denisty;
        
    }
    
    num_particles = (int) particle_list.size();
    
}

void Field_slow::loadParticles(){
    
    
    //Read PGM file (needs rewriting!!)
    
    
    int xsize, ysize;
    double x_pos, y_pos;
    double** _pos_list;
    std::ifstream infile("toystar.pgm");
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
    
    double dx2 = 2.0 / (double) numcols;
    double dy2 = 2.0 / (double) numrows;
    int counter = 0;
    ss >> _read_temp;
    // Following lines : data
    
    for(int row = 0; row < numrows; ++row){
        for (int col = 0; col < numcols; ++col) {
            ss >> _read_temp;
            if(_read_temp > 100){
                x_pos = -(2.0 / 2.0) +dx2*col;
                y_pos = -(2.0 /2.0) + dy2*row;
                particle_list.push_back(new Particle(x_pos, y_pos, counter));
                particle_list.back()->mass = init_mass;
                particle_list.back()->density = init_denisty;
                counter++;
            }
        }
    }
    
    num_particles = (int) particle_list.size();
    
    std::cout << "Num parts = " << num_particles << std::endl;
    //Might put debug
    infile.close();
    
    
}


void Field_slow::initField(){
    
    double x, y;
    Vec<2> _position(0.0);
    int x_idx, y_idx, idx_temp, idx;
    double dx_half = dx / 2.0;
    double dy_half = dy / 2.0;
    
    for(int i = 0; i < num_cells_ver ; i++){ //Horizontal
        for(int j = 0; j < num_cells_hor ; j++){ //Vertical
            idx = num_cells_hor*i + j;
            
            _position[0] =-(width / 2.0)+  i*dx+ dx_half;
            _position[1] = -(height / 2.0) + j*dy + dy_half;
            
            cells_domain.push_back(new Cell(_position, idx,  0));
            
            cells_domain[idx]->cell_width = dx;
            cells_domain[idx]->cell_height = dy;
            
            //Populating neighbour list
            
            for(int k1= -1; k1 < 2; k1++){
                for(int k2= -1; k2 < 2; k2++){
                    
                    if((k1 == 0) && (k2 == 0)){continue;}
                    
                    x_idx = (((i + k1) % num_cells_hor ) +  num_cells_hor) %  num_cells_hor;
                    y_idx = (((j + k2) % num_cells_ver) + num_cells_ver ) % num_cells_ver;
                    
                    idx_temp = num_cells_hor*x_idx + y_idx;
                    cells_domain[idx]->neighbour_index.push_back(idx_temp);
                    
                }
            }
            
            
        }
    }
    
    //Assign particles to cells;
    for(int i = 0; i < num_particles; i++){
        _position = particle_list[i]->position;
        
        x_idx = (int) floor((_position[0] + (width / 2.0)) / dx );
        y_idx = (int) floor((_position[1] + (height / 2.0))/ dy );
        idx_temp = num_cells_hor*x_idx + y_idx;
        particle_list[i]->cell_ID = idx_temp;
        
        particle_list[i]->part_cell_order = (int) cells_domain[idx_temp]->parts_idx.size();
        cells_domain[idx_temp]->parts_idx.push_back(i);
        
    }
    //Create neighbour list
    neighbourListing();
    
    
}


void Field_slow::neighbourListing(){
    Vec<2> _pos_curr(0.0); int _x_temp, _y_temp;
    Vec<2> _pos_other(0.0);
    int _idx_part,  _cell_idx , _cell_idx_neighbour;
    
    
    for(int i = 0; i < num_particles; i++){
        
        _cell_idx = particle_list.at(i)->cell_ID;
        _pos_curr = particle_list.at(i)->position;
        
        //Clear neighbour list
        particle_list.at(i)->neighbours_ID.clear();
        
        //Loop on the particle's home cell
        for(int j = 0; j < cells_domain.at(_cell_idx)->getNumParts() - 1; j++){
            //define dist function!!!!
            _idx_part = cells_domain.at(_cell_idx)->parts_idx.at(j);
            if(_idx_part == i){continue;}
            _pos_other =particle_list.at(_idx_part)->position;
            if(dist(_pos_curr, _pos_other , _search_radius)){
                particle_list.at(i)->neighbours_ID.push_back(_idx_part);
            }
            
        }
        
        //Loop on the neighbouring cells to find all particles in the radius;
        for(int k = 0; k < cells_domain.at(_cell_idx)->getNumNeighbours(); k++){
            _cell_idx_neighbour = cells_domain.at(_cell_idx)->neighbour_index[k];
            
            
            for(int j = 0; j < cells_domain.at(_cell_idx_neighbour)->getNumParts(); j++){
                _idx_part = cells_domain.at(_cell_idx_neighbour)->parts_idx.at(j);
                _pos_other =particle_list.at(_idx_part)->position;
                if(dist(_pos_curr, _pos_other , _search_radius)){
                    particle_list.at(i)->neighbours_ID.push_back(_idx_part);
                }
                
            }
            
        }
        
        
    }
    
    
}

void Field_slow::calcDensity(){
    Vec<2> _pos;
    Vec<2>  _pos_neighbour;
    
    // std::cout << "Width of kernel = " << W_fun.width ;
    double temp_density = 0.0;
    for (int i = 0; i <  num_particles; i++){
        
        _pos = particle_list.at(i)->position;
        particle_list.at(i)->density = particle_list.at(i)->mass*W_fun(_pos, _pos) ;
        
        for(int j = i+1; j < num_particles; j++){
           // _idx_neighbour = particle_list.at(i)->neighbours_ID.at(j);
            _pos_neighbour = particle_list.at(j)->position;
            
            temp_density = (particle_list.at(i)->mass)*W_fun(_pos, _pos_neighbour);
            particle_list.at(i)->density += temp_density;
            particle_list.at(j)->density += temp_density;
            //I need to add the fact that both particles are symmetric.
            
        }
        
        //  std::cout << "Particle [" << i << "] = " <<  particle_list.at(i)->density << std::endl;
    }
    
}

void Field_slow::calcPressure(){
    double _density =0.0;
    for (int i = 0; i <  num_particles; i++){
        _density = particle_list.at(i)->density;
        particle_list.at(i)->pressure = pressure_denisty*(_density)*(_density);
    }
}

void Field_slow::calcForce(){
    
    // calcDensity();
    
    Vec<2>  _pos(0.0);
    Vec<2>  _pos_neighbour(0.0);
    Vec<2>  _force(0.0);
    double _rho_i, _rho_j;
    double _pressure_i, _pressure_j;
    
    double temp;
    
    for (int i = 0; i <  num_particles; i++){
        particle_list.at(i)->resetForce(nu, lambda);
    }
    
    for (int i = 0; i <  num_particles; i++){
               _pos = particle_list.at(i)->position;
        _rho_i = particle_list.at(i)->density ;
        _pressure_i =particle_list.at(i)->pressure;
        
        //particle_list.at(i)->resetForce(nu, lambda);
        
        
        
        for(int j = i+1; j < num_particles; j++){
            
            _pos_neighbour = particle_list.at(j)->position;
            _pressure_j = particle_list.at(j)->pressure;
            _rho_j =particle_list.at(j)->density;
            
            temp =-particle_list.at(i)->mass * (_pressure_i/ (_rho_i*_rho_i) + _pressure_j/ (_rho_j*_rho_j));
            _force = temp*W_fun.gradient(_pos, _pos_neighbour);
             
            particle_list.at(i)->force += _force;
            particle_list.at(j)->force -= _force;
            
        }
        
    }
    
    
    
}

void Field_slow::calcInitVelocity(double _dt){
    
    for(int i = 0; i < num_particles; i++){
        particle_list[i]->corrector(_dt);
    }
    
}

void Field_slow::calcVelocity(double _dt){
    
    for(int i = 0; i < num_particles; i++){
        particle_list[i]->predictor(_dt);
    }
}

void Field_slow::Update_field(){
    
    Vec<2> _position(0.0); int x_idx, y_idx;
    Vec<2> _pos_other(0.0);
    int idx_temp,  _cell_idx , _cell_idx_neighbour;
    
    //Not optimal, need to think of a better way to avoid the 2*Nlog(N)
    //std::cout << "neighbour\n";
    neighbourListing();
    //std::cout << "neighbour1\n";
    int erase_idx = 0;
    
    
    //Reset cell-particle list
    for(int i = 0; i < cells_domain.size(); i++){
        //  std::cout << "Resetting Cell[ " << i<< " ]\n";
        cells_domain.at(i)->reset_cell();
    }
    // std::cout << "neighbour2\n";
    for(int i = 0; i < num_particles; i++){
        //std::cout << "Particle[" << i << "]\n";
        _position = particle_list[i]->position;
        // std::cout << "Particle position[" << i << "]\n";
        x_idx = (int) floor((_position[0] + (width / 2.0)) / dx );
        y_idx = (int) floor((_position[1] + (height / 2.0))/ dy );
        idx_temp = num_cells_hor*x_idx + y_idx;
        //  std::cout << "Particle idx[" << i << "]\n";
        
        if(idx_temp > (int) cells_domain.size()){continue;}
        particle_list[i]->cell_ID = idx_temp;
        //  std::cout << "Particle cell ID[" << i << "] = " << idx_temp << "\n";
        //  std::cout << "Particle cell ID[" << i << "] = " << cells_domain[idx_temp]->getNumParts()<< "\n";
        //  std::cout << "Cell[ " << idx_temp << " ] = ";
        // std::cout <<  cells_domain[idx_temp]->parts_idx.size() << std::endl;
        particle_list[i]->part_cell_order = (int) cells_domain[idx_temp]->parts_idx.size();
        cells_domain[idx_temp]->parts_idx.push_back(i);
        
    }
    
    
}


#endif
