//
//  SPH_system.h
//  SPH
//
//  Created by Nader on 06/06/2016.
//  Copyright (c) 2016 Nader. All rights reserved.
//

#ifndef SPH_SPH_system_h
#define SPH_SPH_system_h

#include <vector>
#include "Vec.h"
#include "Particle.h"

class SPH_system{
    
public:
    int num_particles;
    vector<Particle*> All_particles;            //It might need changing to <Particle>
    double cell_width;
    
    
};


//Place Particle in the correct cell
//Also use it to update particle/cell order
void Particle_cell_order(){
    double _temp_x, _temp_y;
    int _cell_i, _cell_j;
    int _cell_idx, _cell_hsfc;
    for(int i = 0; i < num_particles; i++){
        
        //To find a cell, just divide position by h
        
        //_temp_x = (Particle->position[0])/(cell_width);
        _cell_i = floor(_temp_x);
        _cell_j = floor(_temp_y);
        //_cell_idx = _cell_i*num_cells_x + _cell_j;
        //convert to HSFC
        //Cell.at(_cell_hsfc) //Push an index to a particle
    }
    
    
}

//Construct a Grid and order cells
void Grid_construct(){
    
    double h;
    int num_cells;
    int num_cells_1D;
    
    for(int i=0; i < num_cells_1D; i++){
        for(int j=0; j < num_cells_1D; j++){
            
            Cells->position[0] = h*i + (h/2.0);
            Cells->position[1] = h*j + (h/2.0);
            
        }
        
    }
    
    //Construct cell neighbours
    for(int i = 0; i < num_cells; i++ ){
        
        for(int j = 0; j < 8; j++){
           // Cells->neighbour_index[j] = idx;
        }
      
    }
    
    
    
}


//Find neigbours (might need moving to another class)
void FindNeighbours(){
    
    int _particle_cell;
    int _num_particles_cell;
    Vec<2> _curr_position; //It should be a pointer
    //Define Cell* _curr_cell = Cell.at(_particle_cell);
    for (int i = 0 ; i <  num_particles; i++){
        
        _particle_cell = Particle->cell_ID;
        _curr_position = Particle->position;
        for(int j= 0; j < _num_cell_neighbours; j++){
            
            //_num_particles_cell = _curr_cell->num_particles;
            for(int k = 0; k < _num_particles_cell; k++){
                if(dist(_curr_cell->particle.at(k)->position, _curr_position, h )){
                    //Insert Particle in Neighbour list
                }
                
            }
            
        }
        
        
    }
    
}




#endif
