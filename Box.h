//
//  Box.h
//  SPH
//
//  Created by Nader on 06/06/2016.
//  Copyright (c) 2016 Nader. All rights reserved.
//

#ifndef SPH3_Box_h
#define SPH3_Box_h

#include "Vec.h"
#include "Particle.h"

double COLLISION_RADIUS = 0.01;
//It should be called cell

class Box{
public:
    Vec<2> max_corner;
    Vec<2> min_corner;
    
    Box(Vec<2>, Vec<2>);
    bool in_box(double, Particle*);
    Box(void){};
    ~Box(void) {};
    
    
};

Box::Box(Vec<2> _max, Vec<2> _min){
    max_corner = _max;
    min_corner = _min;
}

bool Box::in_box(double t, Particle *particle_list){
    
    bool result = true;
    double friction = 0.5;
    
    Vec<2> _pos = particle_list->position;
    Vec<2> _vel = particle_list->velocity;
    double mass = particle_list->mass;
  
    
}

#endif
