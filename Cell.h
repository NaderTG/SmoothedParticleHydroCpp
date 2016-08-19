//
//  Cell.h
//  SPH
//
//  Created by Nader on 06/06/2016.
//  Copyright (c) 2016 Nader. All rights reserved.
//

#ifndef SPH_Cell_h
#define SPH_Cell_h

#include "Vec.h"
#include "Particle.h"
#include <vector>


class Cell{
public:
    int cell_ID;
    
    int num_particles;
    double cell_width;
    double cell_height;
    int cell_type; //0: Fluid, 1:hard surface
    Vec<2> position;
    
    Cell();
    Cell(Vec<2>, int, int);
    ~Cell(){};
    
    void setCellID(int );
    
    void setCellType(int);
    
    void setPosition(Vec<2>);
    
    
    std::vector<int> neighbour_index;
    std::vector<int> parts_idx;
    int getNumParts();
    void setCell_type(int _type);
    void eraseFromList(int _idx);
    int getNumNeighbours();
    bool inCell(Particle*);
    
};



Cell::Cell(){
    //These are place holder values
    position[0] = 0.0; position[1] = 0.0;
    
    cell_ID = 1;
    
    num_particles = 1;
    cell_type = 0;
    
}

Cell::Cell(Vec<2> _pos, int _cell_id, int _type){
    position = _pos;
    cell_ID = _cell_id;
    
    cell_type = _type;
}

void Cell::setCellID(int _id){ cell_ID = _id;}
void Cell::setCell_type(int _type) {cell_type = _type;}
void Cell::setPosition(Vec<2> _pos) { position = _pos;}


int Cell::getNumParts(){
    return (int) parts_idx.size();
}

void Cell::eraseFromList(int _idx){
    if((int) parts_idx.size() > 0){
        parts_idx.erase(parts_idx.begin() + _idx);
    }
}
int Cell::getNumNeighbours(){
    return (int) neighbour_index.size();
    
}

bool Cell::inCell(Particle* _part){
    bool _result = false;
    Vec<2> _pos(0.0);
    
    _pos = (_part->position) - position ;
    
    if(abs(_pos[0]) < (cell_width / 2.0) && abs(_pos[1]) < (cell_height / 2.0)){
        _result =  true;
    }
    
    return _result;
}
#endif
