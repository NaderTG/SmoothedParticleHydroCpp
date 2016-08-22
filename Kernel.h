//
//  Kernel.h
//  SPH
//
//  Created by Nader on 10/07/2016.
//  Copyright (c) 2016 Nader. All rights reserved.
//

#ifndef SPH3_Kernel_h
#define SPH3_Kernel_h


#include<cmath>
#include "Vec.h"
#include <vector>
#include <cmath>

//Here we only consider radially symmetric kernels
double POW(double _val, int _power){
    double _result = 1.0;
    for(int i = 0; i < _power; i++){
        _result *= _val;
    }
    return _result;
}
//Don't forget to add derivative and viscousity
class Kernel{ //It should be an abstract class
    
    
public:
    enum kernel_type{ POLY6=1, BSPLINE=2 }; //Add other types
    double width;
    Kernel() {};
    virtual double operator() (Vec<2>& _pos_a, Vec<2>& _pos_b) = 0; //It should be operator
    virtual Vec<2> gradient(Vec<2>& _pos_a, Vec<2>& _pos_b) = 0;
    virtual double laplacian(Vec<2>& _pos_a, Vec<2>& _pos_b) = 0;
    
    
};

//Specific Kernel function: Poly6

class ker_poly6 :public Kernel {
    
    
public:
    double C_1, C_2;
    ker_poly6() {width = 1.0;
        C_1 = 315.0/(64.0*M_PI*POW(width, 9));
        C_2 = 945.0/ (32.0*M_PI * POW(width, 9));
    }; //I should define h^9 and h^2 here
    double operator() (Vec<2>& _pos_a, Vec<2>& _pos_b){
        double result =0.0;
        Vec<2> _r = _pos_a - _pos_b;
        double _r_norm = _r.squaredNorm();
        
        //   result = (4.0/ (M_PI * POW(width, 8)))*POW((width*width - _r_norm*_r_norm),3);
        if(_r_norm <= width && _r_norm >= 0){
            result = C_1*POW((width*width - _r_norm*_r_norm),3);
        }
        return result;
    }
    
    Vec<2> gradient(Vec<2>& _pos_a, Vec<2>& _pos_b){
        Vec<2> _nablaW(0.0);
        Vec<2> _r = _pos_a - _pos_b;
        double _temp =0.0;
        double _r_norm = _r.squaredNorm();
        
        
        if(_r_norm <= width && _r_norm >= 0){
            _temp = C_2*POW((width*width - _r_norm*_r_norm),2);
            _nablaW = _temp*_r;
        }
        
        return _nablaW;
    }
    
    double laplacian(Vec<2>& _pos_a, Vec<2>& _pos_b){
        double result =0.0;
        Vec<2> _r = _pos_a - _pos_b;
        double _r_norm = _r.squaredNorm();
        
        //   result = (4.0/ (M_PI * POW(width, 8)))*POW((width*width - _r_norm*_r_norm),3);
        if(_r_norm <= width && _r_norm >= 0){
            result = C_1*((width*width - _r_norm*_r_norm))*((3.0*width*width - 7.0*_r_norm*_r_norm));
        }
        return result;
    }
    
    void changeWidth(double _width){
        width = _width;
        C_1 = 315.0/(64.0*M_PI*POW(width, 9));
        C_2 = 945.0/ (32.0*M_PI * POW(width, 9));
    }
    
};


class ker_pressure :public Kernel {
    
    
public:
    double C_1, C_2;
    ker_pressure() {width = 1.0;
        C_1 = 15.0/(M_PI*POW(width, 6));
        C_2 = -45.0/(M_PI*POW(width, 6));
    }; //I should define h^9 and h^2 here
    double operator() (Vec<2>& _pos_a, Vec<2>& _pos_b){
        double result =0.0;
        Vec<2> _r = _pos_a - _pos_b;
        double _r_norm = _r.squaredNorm();
        
        //   result = (4.0/ (M_PI * POW(width, 8)))*POW((width*width - _r_norm*_r_norm),3);
        if(_r_norm <= width && _r_norm >= 0){
            result = C_1*POW((width  -  _r_norm),3);
        }
        return result;
    }
    
    Vec<2> gradient(Vec<2>& _pos_a, Vec<2>& _pos_b){
        Vec<2> _nablaW(0.0);
        Vec<2> _r = _pos_a - _pos_b;
        double _temp =0.0;
        double _r_norm = _r.squaredNorm();
        
        
        if(_r_norm <= width && _r_norm >= 0){
            _temp = C_2*POW((width  -  _r_norm),2)*(1.0/_r_norm);;
            _nablaW = _temp*_r;
        }
        
        return _nablaW;
    }
    
    double laplacian(Vec<2>& _pos_a, Vec<2>& _pos_b){
        double result =0.0;
        Vec<2> _r = _pos_a - _pos_b;
        double _r_norm = _r.squaredNorm();
        
        //   result = (4.0/ (M_PI * POW(width, 8)))*POW((width*width - _r_norm*_r_norm),3);
        if(_r_norm <= width && _r_norm >= 0){
            result = C_1*((width*width - _r_norm*_r_norm))*((3.0*width*width - 7.0*_r_norm*_r_norm));
        }
        return result;
    }
    
    void changeWidth(double _width){
        width = _width;
        C_1 = 315.0/(64.0*M_PI*POW(width, 9));
        C_2 = 945.0/ (32.0*M_PI * POW(width, 9));
    }
    
};


class ker_viscosity :public Kernel {
    
    
public:
    double C_1, C_2;
    ker_viscosity() {width = 1.0;
        C_1 = 15.0/(2.0*M_PI*POW(width, 3));
        C_2 = 45.0/ ( M_PI * POW(width, 6));
    }; //I should define h^9 and h^2 here
    double operator() (Vec<2>& _pos_a, Vec<2>& _pos_b){
        double result =0.0;
        Vec<2> _r = _pos_a - _pos_b;
        double _r_norm = _r.squaredNorm();
        
        //   result = (4.0/ (M_PI * POW(width, 8)))*POW((width*width - _r_norm*_r_norm),3);
        if(_r_norm <= width && _r_norm >= 0){
            result = C_1*(-(POW(_r_norm, 3)/ (2.0*POW(width, 3)))  +  (POW(_r_norm/width, 2)) + (width/ (2.0*_r_norm)) - 1.0 );
        }
        return result;
    }
    
    Vec<2> gradient(Vec<2>& _pos_a, Vec<2>& _pos_b){
        Vec<2> _nablaW(0.0);
        Vec<2> _r = _pos_a - _pos_b;
        double _temp =0.0;
        double _r_norm = _r.squaredNorm();
        
        
        if(_r_norm <= width && _r_norm >= 0){
            _temp = C_2*POW((width*width - _r_norm*_r_norm),2);
            _nablaW = _temp*_r;
        }
        
        return _nablaW;
    }
    
    double laplacian(Vec<2>& _pos_a, Vec<2>& _pos_b){
        double result =0.0;
        Vec<2> _r = _pos_a - _pos_b;
        double _r_norm = _r.squaredNorm();
        
        //   result = (4.0/ (M_PI * POW(width, 8)))*POW((width*width - _r_norm*_r_norm),3);
        if(_r_norm <= width && _r_norm >= 0){
            result = C_2*((width  -  _r_norm)) ;
        }
        return result;
    }
    
    void changeWidth(double _width){
        width = _width;
        C_1 = 315.0/(64.0*M_PI*POW(width, 9));
        C_2 = 945.0/ (32.0*M_PI * POW(width, 9));
    }
    
};



class ker_spline :public Kernel {
    
    
public:
    double C_1, C_2;
    ker_spline() {width = 1.0;
        C_1 = 5.0/(14.0*M_PI*POW(width, 2));
        C_2 = 5.0/ (14.0*M_PI * POW(width, 4));
    }; //I should define h^9 and h^2 here
    double operator() (Vec<2>& _pos_a, Vec<2>& _pos_b){
        double result =0.0;
        Vec<2> _r = _pos_a - _pos_b;
        double _r_norm = _r.l2norm();
        double q = _r_norm / width;
        //   result = (4.0/ (M_PI * POW(width, 8)))*POW((width*width - _r_norm*_r_norm),3);
        if(q >= 0 && q < 1){
            result = C_1*( POW(2 - q, 3) - 4.0*POW(1 - q, 3));
        }else if (q >= 1 && q < 2){
            result = C_1*( POW(2 - q, 3) );
        }
        return result;
    }
    
    Vec<2> gradient(Vec<2>& _pos_a, Vec<2>& _pos_b){
        Vec<2> _nablaW(0.0);
        Vec<2> _r = _pos_a - _pos_b;
        double _temp =0.0;
        
        double _r_norm = _r.l2norm();
        double q = _r_norm / width;
        
      
        if(q >= 0 && q < 1){
            _temp = C_2*(-3.0*POW(2.0 - q, 2) + 12.0*POW(1.0 - q, 2))*(1.0/q);
            _nablaW = _temp*_r;
        }else if (q >= 1 && q < 2){
            _temp = C_2*(-3.0*POW(2.0 - q, 2))*(1.0/q);
            _nablaW = _temp*_r;
        }
        
        return _nablaW;
    }
    
    double laplacian(Vec<2>& _pos_a, Vec<2>& _pos_b){
        double result =0.0;
        Vec<2> _r = _pos_a - _pos_b;
        double _r_norm = _r.squaredNorm();
        
        //   result = (4.0/ (M_PI * POW(width, 8)))*POW((width*width - _r_norm*_r_norm),3);
        if(_r_norm <= width && _r_norm >= 0){
            result = C_1*((width*width - _r_norm*_r_norm))*((3.0*width*width - 7.0*_r_norm*_r_norm));
        }
        return result;
    }
    
    void changeWidth(double _width){
        width = _width;
        C_1 = 5.0/(14.0*M_PI*POW(width, 2));
        C_2 = 5.0/ (14.0*M_PI * POW(width, 4));
    }
    
};



#endif
