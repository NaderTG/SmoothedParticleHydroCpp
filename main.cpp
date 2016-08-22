//
//  main.cpp
//  SPH
//  Test commit
//  Created by Nader on 05/06/2016.
//  Copyright (c) 2016 Nader. All rights reserved.
//

#include <iostream>
#include "Particle.h"
#include "Vec.h"
#include "Cell.h"
//#include "field.h"
#include "field_slow.h"
#include "params.h"
#include <vector>
#include "visualise.h"
#include <time.h>

#define PI 3.14159265

void testVec(){
    Vec<2> a(3, 2);
    
    
    Vec<2> b(3.8);
    //Vec<2> c(1);
    Vec<2> c ;
    c = b+a;
    std::cout << "Hello, World!\n";
    std::cout << "A = " << c[0] << std::endl;
    
    c = a;
    std::cout << "A = " << c[0] << std::endl;
    // a.print();
    
    
}


void testParticle(){
    Vec<2> a(1.0, 3.1);
    Vec<2> b(0.1, 0.5);
    Particle test(a, b, 2.0, 1.0);
    std::vector<Particle*> container;
    
    container.push_back(new Particle(a, b, 2.0, 1.0));
    std::cout << "Particle's density = " << test.density << std::endl;
    std::cout << "Particle's mass = " << test.mass << std::endl;
    std::cout << "Particle's position = [" << test.position[0] << ", " << test.position[1] << "]\n";
    std::cout << "Particle's velocity = [" << test.velocity[0] << ", " << test.velocity[1] << "]\n";
    
    a[0] = 4.0;
    test.position = a;
    std::cout << "Particle's position = [" << test.position[0] << ", " << test.position[1] << "]\n";
    std::cout << std::endl;
    std::cout << "Particle's density = " << container.at(0)->density<< std::endl;
    std::cout << "Particle's mass = " << container.at(0)->mass << std::endl;
    std::cout << "Particle's position = [" << container.at(0)->position[0] << ", " << container.at(0)->position[1] << "]\n";
    std::cout << "Particle's velocity = [" << container.at(0)->velocity[0] << ", " << container.at(0)->velocity[1] << "]\n";
    container.at(0)->position = a;
    std::cout << "\nParticle's position = [" << container.at(0)->position[0] << ", " << container.at(0)->position[1] << "]\n";
    
    for(int i = 0; i < 4; i++){
        container.at(0)->neighbours_ID.push_back(i);
    }
    
    std::cout << "Number of neighbours = " <<container.at(0)->neighbour_size() << std::endl;
    
}

void testParticle2(){
    
    //Initialialise particles
    std::vector<Particle*> container;
    
    Vec<2> _position(0.0);
    int num_part = 8;
    
    double radius = 0.25; double _x, _y;
    
    double theta = (2.0*PI)/ (double)  num_part ;
    
    
    std::cout << "theta = " << theta << std::endl;
    
    for(int i = 0; i <  num_part ; i++){
        _x = radius*(sin(i*theta)) + 0.5;
        _y = radius*(cos(i*theta)) + 0.5;
        //    std::cout << "x = [" << _x << " , " << _y << "]\n";
        container.push_back(new Particle(_x, _y, i));
        container.back()->mass = 1.0 + 0.1*i;
        container.back()->density = 1.0;
        
    }
    
    Vec<2> a(1.0, 3.1);
    Vec<2> b(0.1, 0.5);
    
    container[0]->position = a;
    container[2]->position = b;
    
    
    for(int i = 0; i < num_part; i++){
        std::cout << "Particle[" << i << "] = (" << container[i]->position[0] << ", " << container[i]->position[1] << ")\n";
    }
    
}
//
//bool dist(Vec<2>& _pos_a, Vec<2>& _pos_b , double radius){
//    bool result = false;
//    Vec<2> _r = _pos_a - _pos_b;
//    double _r_norm = _r.l2norm();
//    
//    result = (_r_norm < radius);
//    
//    return result;
//}

void testPartAllocation(){ //To be done
    
    
    //Initialialise particles
    std::vector<Particle*> container;
    std::vector<Cell*> cells_domain;
    Vec<2> _position(0.0);
    int num_part = 20;
    
    double radius = 0.25; double _x, _y;
    double L = 2.0;
    double theta = (2.0*PI)/ (double)  (num_part-1) ;
    
    //Cell stuff
    int num_cells_ver = 10;
    int num_cells_hor = 10;
    
    double dx = L /(double) num_cells_hor;
    double dy = L /(double) num_cells_ver;
    double dx_half = dx / 2.0;
    double dy_half = dy/ 2.0;
    int x_idx, y_idx, idx_temp, idx;
    
    Vec<2> _pos2(0.0); int _x_temp, _y_temp;
    Vec<2> _pos_other(0.0);
    int _idx_part,  _cell_idx , _cell_idx_neighbour;
    //int k1_p, k2_p; //P stands for prime k_1', k_2'
    
    for(int i = 0; i < num_cells_ver ; i++){ //Horizontal
        for(int j = 0; j < num_cells_hor ; j++){ //Vertical
            idx = num_cells_hor*i + j;
            
            _position[0] = -L/2.0 +  i*dx+ dx_half;
            _position[1] = -L/2.0  + j*dy + dy_half;
            
            cells_domain.push_back(new Cell(_position, idx,  0));
            cells_domain[idx]->cell_width = dx;
            cells_domain[idx]->cell_height = dy;
            for(int k1= -1; k1 < 2; k1++){
                for(int k2= -1; k2 < 2; k2++){
                    
                    if((k1 == 0) && (k2 == 0)){continue;}
                    
                    x_idx = (((i + k1) % num_cells_hor ) +  num_cells_hor) %  num_cells_hor;
                    y_idx = (((j + k2) % num_cells_ver) + num_cells_ver ) % num_cells_ver;
                    
                    idx_temp = num_cells_hor*x_idx + y_idx;
                    cells_domain[idx]->neighbour_index.push_back(idx_temp);
                    
                }
            }
            // cells_domain[idx]->neighbour_index.erase(cells_domain[idx]->neighbour_index.begin());
            
        }
    }
    
    
    std::cout << "theta = " << theta << std::endl;
    
    
    
    
    //std::cout << "theta = " << theta << std::endl;
    
    for(int i = 0; i <  num_part ; i++){
        radius = (double)(i + 1)/(double) (num_part);
        _x = radius*(sin((i)*theta)) ;
        _y = radius*(cos(i*theta)) ;
        //    std::cout << "x = [" << _x << " , " << _y << "]\n";
        container.push_back(new Particle(_x, _y, i));
        container.back()->mass = 1.0;
        container.back()->density = 1.0;
        
    }
    
    for(int i = 0; i < num_part; i++){
        _position = container[i]->position;
        
        x_idx = (int) floor((_position[0] + L/2.0 ) / dx );
        y_idx = (int) floor((_position[1] + L/2.0 ) / dy );
        idx_temp = num_cells_hor*x_idx + y_idx;
        
        container[i]->cell_ID = idx_temp;
        std::cout << "Particle [" << i << "] [ " << _position[0] << ", " << _position[1] << "] goes to cell[" << idx_temp << "]\n";
        container[i]->part_cell_order = (int) cells_domain[idx_temp]->parts_idx.size();
        cells_domain[idx_temp]->parts_idx.push_back(i);
    }
    
    
    std::cout << "Cell initial positions" << std::endl;
    for(int i = 0; i < num_cells_ver ; i++){ //Horizontal
        for(int j = 0; j < num_cells_hor ; j++){ //Vertical
            idx = num_cells_hor*i + j;
            std::cout << "Cell [" <<idx << "] [" <<cells_domain[idx]->position[0] << " , " << cells_domain[idx]->position[1] << "] = ";
            
            if((int) cells_domain[idx]->parts_idx.size() > 0){
                for(int kk = 0; kk < cells_domain[idx]->parts_idx.size()  ; kk++){
                std::cout << " particle [" << cells_domain[idx]->parts_idx[kk] << "]\n";
                }
            }
            else{
                std::cout << " Empty cell \n";
            }
        }
    }
    
    std:: cout << "neighbourlist \n";
    //Should add particle neighbour list test
    double _search_radius = 0.26;
    for(int i = 0; i < num_part; i++){
        
        _cell_idx = container.at(i)->cell_ID;
        _pos2 = container.at(i)->position;
        
        
        std::cout << "Part[" << i << "] in Cell[" << _cell_idx << "][" << _pos2[0] << ", " << _pos2[1] << "] = ";
        //Clear neighbour list
        container.at(i)->neighbours_ID.clear();
        
        //Loop on the particle's home cell
        //std:: cout << "neighbourlist 0.1 \n";
        for(int j = 0; j < cells_domain.at(_cell_idx)->getNumParts()  ; j++){
            //define dist function!!!!
            _idx_part = cells_domain.at(_cell_idx)->parts_idx.at(j);
            if(_idx_part == i){continue;}
            _pos_other =container.at(_idx_part)->position;
            if(dist(_pos2, _pos_other , _search_radius) && _idx_part > i ){
                container.at(i)->neighbours_ID.push_back(_idx_part);
                std::cout << "neighbour cell[" << _cell_idx << " ]= " <<  _idx_part << "\n";
           }
            
        }
        
        //std:: cout << "neighbourlist 0.2\n";
        //Loop on the neighbouring cells to find all particles in the radius;
        for(int k = 0; k < cells_domain.at(_cell_idx)->getNumNeighbours(); k++){
            _cell_idx_neighbour = cells_domain.at(_cell_idx)->neighbour_index[k];
            
            
            for(int j = 0; j < cells_domain.at(_cell_idx_neighbour)->getNumParts(); j++){
               
             //   std::cout << "There are particles\n" ;
                //                //define dist function!!!!
                _idx_part = cells_domain.at(_cell_idx_neighbour)->parts_idx.at(j);
                _pos_other =container.at(_idx_part)->position;
                if(dist(_pos2, _pos_other , _search_radius)  && _idx_part > i ){
                     std::cout << "neighbour cell[" << _cell_idx_neighbour << " ]= " <<  _idx_part << "\n";
                    container.at(i)->neighbours_ID.push_back(_idx_part);
                    //std::cout << " , " <<  _idx_part;
                }
                //
            }
            
        }
        
        
    }
    
    std:: cout << "neighbourlist2 \n";
    
    for(int i = 0; i <  num_part ; i++){
        
        
        std::cout << "Particle[" << i  << "] = ";
        for(int j = 0; j < container[i]->neighbour_size(); j++){
            
            std::cout << container[i]->neighbours_ID[j] << ", ";
            
        }
        std::cout << std::endl;
    }
    
    //Particle update positions
    radius = 0.4;
    double shift = 0.39269908169872;
    //  int _cell_idx;
    std::cout << "Update positions" << std::endl;
    for(int i = 0; i <  num_part ; i++){
        _x = container[i]->position[0]*(sin(i*theta + shift)) ;
        _y = container[i]->position[1] *(cos(i*theta + shift));
        //    std::cout << "x = [" << _x << " , " << _y << "]\n";
        container[i]->position[0] = _x;
        container[i]->position[1] = _y;
        
        
    }
    std:: cout << "Update cell \n";
    //
    //      Vec<2> _pos2(0.0); int _x_temp, _y_temp;
    //    for(int i = 0; i < num_part; i++){
    //                _pos2 =container.at(i)->position;
    //                _x_temp = (int) floor(_pos2[0] / dx);
    //                _y_temp = (int) floor(_pos2[1] / dy);
    //                _cell_idx = _x_temp*num_cells_hor + _y_temp;
    //            std::cout << "\nParticle's position = [" << _pos2[0]<< ", " << _pos2[1] << "] \t cell = " <<_cell_idx <<" \n";
    //
    //    }
    
    //update cell locations
    
    for(int i = 0; i < num_part; i++){
        
        _pos2 =container.at(i)->position;
        _x_temp = (int) floor((_pos2[0] + L/2.0 ) / dx);
        _y_temp = (int) floor((_pos2[1] + L/2.0 ) / dy);
        _cell_idx = _x_temp*num_cells_hor + _y_temp;
        
        cells_domain.at(_cell_idx)->parts_idx.push_back(i);
        container.at(i)->cell_ID = _cell_idx;
    }
    
    
    std::cout << "Cell final positions" << std::endl;
    for(int i = 0; i < num_cells_ver ; i++){ //Horizontal
        for(int j = 0; j < num_cells_hor ; j++){ //Vertical
            idx = num_cells_hor*i + j;
            std::cout << "Cell [" <<idx << "] [" <<cells_domain[idx]->position[0] << " , " << cells_domain[idx]->position[1] << "] = ";
            
            if((int) cells_domain[idx]->parts_idx.size() > 0){
                std::cout << " particle [" << cells_domain[idx]->parts_idx[0] << "]\n";
            }
            else{
                std::cout << " Empty cell \n";
            }
        }
    }
}

void testPartUpdate(); //To be done Actually, this is done previously

void testForce(); //To be done

void register_plot(Field& _a, int time){
    
    std::string _filename = "toystar3_" + std::to_string(time) + ".dat";
    
    int _num_parts = _a.num_particles;
    
    std::cout << "Number of particles = " << _num_parts << std::endl;
    
    std::ofstream output_file (_filename);
    
    if(output_file.is_open()){
        
        output_file << _a;
    }
    
    output_file.close();
    
}

void register_plot(Field_slow& _a, int time){
    
    std::string _filename = "toystar3_" + std::to_string(time) + ".dat";
    
    int _num_parts = _a.num_particles;
    
    std::cout << "Number of particles = " << _num_parts << std::endl;
    
    std::ofstream output_file (_filename);
    
    if(output_file.is_open()){
        
        output_file << _a;
    }
    
    output_file.close();
    
}

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    //  testVec();            //Done!
    //  testParticle();       //Done!
   //   testPartAllocation(); //Done!
    
    
    // testParticle2();
    
    //
        //How the code should be
    
    Params prob_params("test");
    //     //loadParams();
    //
    double dt =prob_params.dt;
    //Field_slow domain(prob_params);
      Field domain(prob_params);
    //     Kernel ker_function;
    //Main time integration
    int count = 1;
    int UPDATE_ITER = 20;
    int PLOT_ITER = 100;
    double time = 0.0;
    std::string filename_vis = "toy_star";
    visualise vis_part(filename_vis );
    int vis_count = 0;
    
    
    //vis_part.vis_time_step(domain,vis_count );
   // register_plot(domain, vis_count );
    
 
  //  Initial step
    domain.calcDensity();
    
    domain.calcPressure();
    
    domain.calcForce();
   // domain.calcInitVelocity(dt);

    
    //
    std::cout << "sim begin with end time = " <<  prob_params.t_end <<" \n";
 
//  //std::cout << "Domain = \n" << domain ;
//    for(int k = 0; k < 5; k++){
//            domain.Update_field();
//            domain.calcVelocity(dt);
//            domain.calcDensity();
//            domain.calcPressure();
//            domain.calcForce();
//        
//    }
//    std::cout << "Domain = \n" << domain ;
    clock_t start = clock();
    while (time < prob_params.t_end && count < 402){
        
    //    domain.Update_field();
        // domain.updatePositions();
       
         if( (count % UPDATE_ITER) == 0){
            //std::cout << "updating\n";
            domain.Update_field();
        }
////
         domain.calcVelocity(dt);
//        if(((count-1) % PLOT_ITER) == 0){
//            //vtk stuff
//            vis_count++;
//        //    std::cout << "plotting" << vis_count << "\n";
//          //  vis_part.vis_time_step(domain, vis_count );
//            register_plot(domain, vis_count );
//        }
       // std::cout << "Time = " << time << "\n";
        
        domain.calcDensity();
        domain.calcPressure();
        domain.calcForce();
        
        time += dt;
        count++;
    }
    
    vis_count++;
     register_plot(domain, vis_count );
    clock_t stop = clock();
    double elapsed = (double)(stop - start)/ CLOCKS_PER_SEC;
    std::cout << "sim_end in " << elapsed << " seconds \n" ;
//
//    
//    //    Vec<2> a(1.0, 3.1);
//    //    Vec<2> b(0.1, 0.5);
//    //
//    //    Vec<2> _r;
//    //    _r[0] = -1.0*a[0];
//    //    _r[1] = -20.*b[1];
//    //    
//    //    std::cout << "position = " << _r << std::endl;
//    
//
//    
    
    return 0;
}
