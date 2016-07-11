//
//  Kernel.h
//  SPH
//
//  Created by Nader on 10/07/2016.
//  Copyright (c) 2016 Nader. All rights reserved.
//

#ifndef SPH_Kernel_h
#define SPH_Kernel_h

#include<cmath>
#include "Vec.h"
#include "Particle.h"
#include <vector>
#include <cmath>

//Here we only consider radially symmetric kernels

//Don't forget to add derivative and viscousity
class Kernel{ //It should be an abstract class
private:
    double width;
public:
    enum kernel_type{ POLY6=1, BSPLINE=2 }; //Add other types
    
    Kernel() {};
    virtual double operator() (Vec<2>& _pos_a, Vec<2>& _pos_b) = 0; //It should be operator
    virtual Vec<2> gradient() (Vec<2>& _pos_a, Vec<2>& _pos_b) = 0;
    
};

//Specific Kernel function: Poly6

class ker_poly6 :public Kernel {
    

public:
    ker_poly6() {width = 0.1;};
    double operator() (Vec<2>& _pos_a, Vec<2>& _pos_b){
        double result;
        Vec<2> _r = _pos_a - _pos_b;
        double _r_norm = _r.squaredNorm();
        
        result = (4.0/ (M_PI * POW(width, 8)))*POW((width*width - _r_norm*_r_norm),3);
        
        return result;
    }
    
    Vec<2> gradient(){
        Vec<2> _nablaW;
        Vec<2> _r = _pos_a - _pos_b;
        double _r_norm = _r.squaredNorm();
        _nablaW = _r* (24.0/ (M_PI * POW(width, 8)))*POW((width*width - _r_norm*_r_norm),2);
        return _nablaW;
    }
    
    
    
};
#endif
