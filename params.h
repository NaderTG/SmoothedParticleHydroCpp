#ifndef PARAMS
#define PARAMS

#include <cmath>
class Params{

    int num_particles;
    double init_mass;
    double init_density;

    int num_cells;
    double width, height, dx;
    double t_end, dt;
    int numTimeSteps;

    Params();
    Params(int, int, double, double, int); //Num_particle, num_cells, width, time, num_time

};

Params::Params(){

    int num_cell_sqt = (int) sqrt((double) num_cells);
    num_particles = 50;
    num_cells = 9;
    init_density = 1.0;
    init_mass = 1.0;
    t_end = 10.0;
    numTimeSteps = 100;
    dt = t_end/ ((double) numTimeSteps);
    width = 1.0; height = 1.0;
    dx = width / ((double) num_cell_sqt);

}

Params::Params(int _num_part, int _num_cells, double _width, double _time, int _Nt){

    int num_cell_sqt = (int) sqrt((double) _num_cells);
    num_particles = _num_part;
    num_cells = _num_cells;
    width = _width;
    height = _width;
    t_end = _time;
    numTimeSteps = _Nt;

    init_density = 1.0;
    init_mass = 1.0;
    dx = width / ((double) num_cells);

}

#endif // PARAMS

