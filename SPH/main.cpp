//
//  main.cpp
//  SPH
//
//  Created by Nader on 05/06/2016.
//  Copyright (c) 2016 Nader. All rights reserved.
//

#include <iostream>
#include "Particle.h"
#include "Vec.h"
#include <vector>
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

    std::cout << "Number of neighbours = " <<container.at(0)->neighbor_size() << std::endl;

}


void testCell(){

    Vec<2> _pos(1.0, 3.0);
    std::vector<Cell*> container;
    container.push_back(new Cell(_pos, 1, 2, 0));

    container.at(0)->neighbour_index.push_back(1);
    container.at(0)->particle_idx.push_back(1);

    std::cout << "Cell Position x = " << container.at(0)->position[0] << std::endl;
    std::cout << "Cell Position y = " << container.at(0)->position[1] << std::endl;
    std::cout << "Number of particles = " << container.at(0)->getNumParts() << std::endl;
    std::cout << "Number of neighbours = " << container.at(0)->getNumNeighbours() << std::endl;

    std::cout << "First particles = " << container.at(0)->particle_idx.at(0) << std::endl;



}

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";

    //testVec();
    //testParticle();
    testCell();
    return 0;
}
