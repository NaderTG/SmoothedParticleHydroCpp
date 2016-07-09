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
    int cell_hsfc_ID;
    int num_particles;
    double cell_width;
    int cell_type; //0: Fluid, 1:hard surface
    Vec<2> position;
    
    Cell();
    Cell(Vec<2>, int, int, int);
    ~Cell(){};
    
    void setCellID(int );
    void setCellHSFC(int);
    void setCellType(int);
    
    void setPosition(Vec<2>);
    
    vector<Particle> parts_in_cell;
    vector<int> neighbour_index;
    vector<int> parts_idx;
    int getNumParts();
    int getNumNeighbours();
    
};



Cell::Cell(){
    //These are place holder values
    position[0] = 0.0; position[1] = 0.0;
    
    cell_ID = 1;
    cell_hsfc_ID = 1;
    num_particles = 1;
    cell_type = 0;
    
}

Cell::Cell(Vec<2> _pos, int _cell_id, int _cell_hsfc, int _type){
    position = _pos;
    cell_ID = _cell_id;
    cell_hsfc_ID = _cell_hsfc;
    cell_type = _type;
}

void Cell::setCellID(int _id){ cell_ID = _id;}
void Cell::setCellHSFC(int _id){cell_hsfc_ID = _id;}
void Cell::cell_type(int _type) {cell_type = _type;}
void Cell::setPosition(Vec<2> _pos) { position = _pos;}


int Cell::getNumParts(){
    return parts_in_cell.size();
}

int Cell::getNumNeighbours(){
    return neighbour_index.size();
    
}
#endif
