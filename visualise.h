//
//  visualise.h
//  SPH3
//
//  Created by Nader on 19/08/2016.
//  Copyright (c) 2016 Nader. All rights reserved.
//

#ifndef SPH3_visualise_h
#define SPH3_visualise_h

#include "Vec.h"
#include "Particle.h"
#include "Field.h"
#include <vector>
#include <string>

class visualise {
    //At the moment I can only do scatter plot
    std::string filename;
    
    visualise(std::string& _filename){
        filename = _filename;
        
    }
    visualise(){
        filename = "test";
    }
    
    void vis_time_step(Field& _field, double _time);
    
};

void visualise::vis_time_step(Field& _field, double _time){
    
    
    std::string _filename = filename + std::to_string(_time) + ".dat";
    int _num_parts = _field.num_particles;
    
    ofstream output_file (_filename);
    
    if(output_file.is_open()){
        
        output_file << _field;
    }
    
    output_file.close();
}

#endif
