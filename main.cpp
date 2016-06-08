//
//  main.cpp
//  SPH
//
//  Created by Nader on 05/06/2016.
//  Copyright (c) 2016 Nader. All rights reserved.
//

#include <iostream>

#include "Vec.h"

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    
    Vec<2> a(1.0);
    a[0] = 3;
    a[1] = 2;
    
    Vec<2> b(3.8);
    //Vec<2> c(1);
    Vec<2> c = b + a;
    
    std::cout << "Hello, World!\n";
   std::cout << "A = " << c[0] << std::endl;
    
    c = a;
    std::cout << "A = " << c[0] << std::endl;
   // a.print();
    return 0;
}
