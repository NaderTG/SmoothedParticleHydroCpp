//
//  Particle.h
//  SPH
//
//  Created by Nader on 05/06/2016.
//  Copyright (c) 2016 Nader. All rights reserved.
//

#ifndef SPH_Particle_h
#define SPH_Particle_h


#include "Vec.h"
#include <vector>


//template<int N>
class Particle{
public:
    
    Vec<2> position;
    Vec<2> velocity;
    Vec<2> prev_position;
    
    double mass, density;
    int cell_ID;
    int particle_ID;
    
    std::vector<int> neighbours_ID;
    
    Particle(void);
    Particle(Vec<2>, Vec<2>, double, double);
    ~Particle(){};
    int neighbor_size();
    void setPosition(double x, double y);
    
};


Particle::Particle(){
    
    mass = 0.001; density = 1.0;
    
}

Particle::Particle(Vec<2> _pos, Vec<2> _vel, double _mass, double _density){
    mass = _mass; density = _density;
    position = _pos; velocity = _vel;
    
    
}

Particle::setPosition(double x, double y){
    Vec<2> _pos(x, y);
    position = _pos;
}

Particle::neighbor_size(){
    return neighbours_ID.size();
}

#endif
