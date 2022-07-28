/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   multipole_method_functions.h
 * Author: Marshall
 *
 * Created on July 4, 2022, 2:52 PM
 */

#ifndef MULTIPOLE_METHOD_FUNCTIONS_H
#define MULTIPOLE_METHOD_FUNCTIONS_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include "structures.h"
#include "constants.h"

namespace MultipoleMethodFunction {
    void morton1(vec3<int> boxIndex3D, int& boxIndex, int numLevel);
    void unmorton(int boxIndex, vec3<int>& boxIndex3D);
    void morton3(vec3<int> boxIndex3D, int& boxIndex, int numLevel);
    void unmorton3(int boxIndex, vec3<int>& boxIndex3D);
    void morton11(vec3<int> boxIndex3D, int& boxIndex, int numLevel = 1);
    void unmorton11(int boxIndex, vec3<int>& boxIndex3D);
    void mortonGeneral(vec3<int> boxIndex3D, int& boxIndex, int numLevel, int boxSize);
    void unmortonGeneral(int boxIndex, vec3<int>& boxIndex3D, int boxSize);
    void signedMorton3(vec3<int> boxIndex3D, unsigned int& boxIndex);
    void signedUnmorton3(unsigned int boxIndex, vec3<int>& boxIndex3D);

    void balancedTernary(vec3<int> boxIndex3D, int& boxIndex);
    void undoBalancedTernary(int boxIndex, vec3<int>& boxIndex3D);
}

#endif /* MULTIPOLE_METHOD_FUNCTIONS_H */

