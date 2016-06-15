#ifndef FIELD_H
#define FIELD_H

#include <iostream>
#include "Particle.h"
#include "Cell.h"
#include "Vec.h"
#include "params.h"
#include <vector>

class Field{
public:
    //Implementation of periodic domain
    std::vector<Cell*> cells_domain;
    std::vector<Particle*> particle_list;

    int num_cells;
    int num_particles;

    double width, height; //Total width of the field
    double w_h, h_h; //cell width
    int num_cells_w, num_cells_h; //number of cells of each diractions

    Field();
    Field(Params);
    
    ~Field(){};

    void initField();
    void Update_field();


};

Field::Field(){
    Params _params();
    num_cells = _params.num_cells;
    num_cells_w = (int) (sqrt(num_cells));
    num_cells_h = (int) (sqrt(num_cells));

    num_particles = _params.num_particles;
    width = _params.width;
    w_h = _params.dx;
    h_h = _params.dx;

    intiField();

}
Field::Field(Params _params){

    num_cells = _params.num_cells;
    num_particles = _params.num_particles;
    width = _params.width;
    w_h = _params.dx;
    h_h = _params.dx;
    initField();

}

void Field::initField(){

    //double _x_pos, _y_pos;
    Vec<2> _pos(0.0);
    Vec<2> _vel(0.0);

    for(int i = 0; i < num_particles; i++){
        //position random
        // d = i*w_h; //make it an integrer as h is equal to one
        // _pos = d2hscf(d, L, w_h);
        //Particle Radius
        particle_list.push_back(Particle( _pos, _vel,1.0, 1.0));
        particle_list.at(i)->particle_ID = i;
    }


    for(int i = 0; i < num_cells_w; i++ ){

        for(int j = 0; j < num_cells_h; j++){
            _pos[0] =(double)(i*w_h) + w_h/2.0;
            _pos[1] = (double)(j*h_h) + h_h/2.0;
            //Cell(Vec<2>, int, int, int)
            cells_domain.push_back(new Cell(_pos, (i*num_cells_h + j), (i*num_cells_h + j), 0 ));
            //d = hsfc2pos(_pos, w_h, L);
            //

        }

    }
    
    Vec<2> _pos2(0.0); int _x_temp, _y_temp, _cell_idx;
    //Put Particles in cells
    for(int i = 0; i < num_particles; i++){
     
        _pos2 =particle_list.at(i).position;
        _x_temp = (int) floor(_pos2[0] / h);
        _y_temp = (int) floor(_pos2[1] / h);
        _cell_idx = _x_temp*num_cells_w + _y_temp;
        
        cells_domain.at(_cell_idx)->parts_idx.push_back(i);
        //Add the cell number to the particle
        particle_list.at(i)->cell_ID = _cell_idx;
    }



}

void Field::Update_field(){
    
    Vec<2> _pos2(0.0); int _x_temp, _y_temp, _cell_idx;
    Vec<2> _pos_other(0.0);
    int _idx_part;
    for(int i = 0; i < num_particles; i++){
        _cell_idx = particle_list.at(i)->cell_ID;
        _pos2 = particle_list.at(i)->position;
        //Clear neighbour list
        
        for(int j = 0; j < cells_domain.at(_cell_idx)->getNumParts(); j++){
            //define dist function!!!!
            _idx_part = cells_domain.at(_cell_idx)->parts_idx.at(j);
            _pos_other =particle_list.at(_idx_part)->position;
            if(dist(_pos2, _pos_other , h)){
                particle_list.at(i)->neighbours_ID.push_back(_idx_part);
            }
            
        }
        //Loop on the neighbouring cells to find all particles in the radius;
        for(int j = 0; j < cells_domain.at(_cell_idx)->getNumNeighbours(); j++){
            
            
            
            
        }
        
        
        _pos2 =particle_list.at(i).position;
        _x_temp = (int) floor(_pos2[0] / h);
        _y_temp = (int) floor(_pos2[1] / h);
        _cell_idx = _x_temp*num_cells_w + _y_temp;
        
        cells_domain.at(_cell_idx)->parts_idx.push_back(i);
        //Add the cell number to the particle
        particle_list.at(i)->cell_ID = _cell_idx;
        
    }
    
    
}
