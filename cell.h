/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   cellTest.h
 * Author: Marshall
 *
 * Created on April 21, 2022, 3:57 PM
 */

#ifndef CELL_H
#define CELL_H

#include <cstdlib>
#include "structures.h"
#include "auxiliary_functions.h"

class Cell {
public:
    int index;                      // Particle index
    int numSubcells;                // Number of subcells
    int subcellBlock;
    float charge;                   // Net charge
    float unsignedcharge;           // Sum of absolute value of charge
    vec3<float> location;           // Location of minimum x, y, z
    vec3<float> centerOfCharge;     // Location of center of charge
    vec3<float> dimensions;         // Size of x, y, z dimensions
    float repLength;                // Maximum of dimensions
    float delta;                    // Distance between geometric center and center of charge
    Cell* firstChild = NULL;    // Pointer to first child
    Cell* nextSibling = NULL;   // Pointer to next sibling,
                                    // or parent's next sibling if no subsequent siblings
    Cell(vec3<float> length);
    ~Cell();
    void generateSubcells(int subcellBlock);        // Generate subcells
    void setLocationOfSubcells();   // Set subcell locations
};

#endif /* CELLTEST_H */

