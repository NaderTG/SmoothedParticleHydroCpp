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



}
