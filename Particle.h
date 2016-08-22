//
//  Particle.h
//  SPH
//
//  Created by Nader on 05/06/2016.
//  Copyright (c) 2016 Nader. All rights reserved.
//

#ifndef SPH3_Particle_h
#define SPH3_Particle_h


#include "Vec.h"
#include <vector>


//template<int N>
class Particle{
public:
    
    Vec<2> position;
    Vec<2> velocity;
    Vec<2> velocity_h;
    Vec<2> velocity_m;
    Vec<2> force;
    Vec<2> prev_position;
    
    double mass, density, pressure;
    int cell_ID;
    int particle_ID;
    int part_cell_order;
    std::vector<int> neighbours_ID;
    
    Particle(void);
    Particle(Vec<2>, Vec<2>, double, double);
    Particle(double, double, int);
    ~Particle(){};
    int neighbour_size();
    void setPosition(double, double);
    void resetForce();
    void resetForce(double, double);
    void assignCell(int);
    void assignID(int);
    void predictor(double);
    void corrector(double);
};


Particle::Particle(){
    
    mass = 0.001; density = 1.0; pressure = 1.0;
    
}

Particle::Particle(Vec<2> _pos, Vec<2> _vel, double _mass, double _density){
    mass = _mass; density = _density; pressure = 1.0;
    position = _pos; velocity = _vel;
    velocity_h = _vel; velocity_m = _vel;
    //force[0] = 0.0; force[1] = 0.0;
    force[0] = 0.0; force[1] = 0.0;
}

Particle::Particle(double _x, double _y, int _part_ID){
    mass = 1.0; density = 1.0;
    position[0] = _x; position[1] = _y;
    force[0] = 0.0; force[1] = 0.0;
    velocity[0] = 0.0; velocity[1] = 0.0;
    velocity_h[0] = 0.0; velocity_h[1] = 0.0;
    velocity_m[0] = 0.0; velocity_m[1] = 0.0;
    particle_ID = _part_ID;
}

void Particle::setPosition(double x, double y){
    Vec<2> _pos(x, y);
    
    position = _pos;
}

int Particle::neighbour_size(){
    return (int)neighbours_ID.size();
}

void Particle::assignCell(int _cell_id){
    cell_ID = _cell_id;
}

void Particle::assignID(int _part_ID){
    particle_ID = _part_ID;
}

void Particle::resetForce(){
    force[0] = 0.0; force[1] = 0.0;
}

void Particle::resetForce(double _nu, double _lambda){
    
    force = -_nu*(velocity) - _lambda*position;
//    force[0] = -_nu*velocity[0] - _lambda*position[0];
//    force[1] = -_nu*velocity[1] - _lambda*position[1];
//    
    
}

void Particle::predictor(double _dt){
    
//    velocity_h[0] += force[0]*_dt;
//    velocity_h[1] += force[1]*_dt;
//    
//    position[0] += velocity_h[0] *_dt;
//    position[1] += velocity_h[1] *_dt;
//    
//    velocity[0] = velocity_h[0] + 0.5*force[0] * _dt;
//    velocity[1] = velocity_h[1] + 0.5*force[1] * _dt;
    
    velocity_h = velocity_m + force*_dt;
    position = position + velocity_h*_dt;
    velocity = 0.5*( velocity_h +velocity_m );
    velocity_m = velocity_h ;
    
    
}

void Particle::corrector(double _dt){
    velocity_h[0] = velocity[0] +  0.5*force[0]*_dt;
    velocity_h[1] = velocity[1] +  0.5*force[1]*_dt;
    
    position[0] += velocity_h[0] *_dt;
    position[1] += velocity_h[1] *_dt;
    
    velocity[0] += force[0] * _dt;
    velocity[1] += force[1] * _dt;
    
}

#endif
