/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   constants.h
 * Author: Marshall
 *
 * Created on June 7, 2022, 5:16 PM
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include <complex>
#include <algorithm>
#include <utility>
#include <pthread.h>

//#undef _PTHREAD_H

namespace PhysicalConstant {
    const long double PI = 3.14159265358979323846264338328L;
    const double CHARGE = 1;
    const double SOFTENING = 6.4e-19; /* softening parameter to prevent 1/r^2 singularity
                        currently using approximated size (diameter) of EMI-BF4 squared */
    const float ELEM_CHARGE = 1.6021766e-19f; // elementary charge
    const float AMU = 1.66053904e-27f; // mass / 1 amu
    const float FREE_SPACE_PERMITTIVITY = 8.85419e-12f; // permittivity of free space
    const float BOLTZMANN_CONST = 1.38065e-23f; // Boltzmann constant
    const float EMI_BF4_DISTANCE = 4.25e-10f; // distance between EMI and BF4 (for fragmentation) 
    const float EMI_ION_RADIUS = 0.5e-9f; // Radius of the EMI ion (for fragmentation)
    const float EFF_SURFACE_TENSION = 0.035; // Effective Surface tension of an EMIBF4 molecule
    const float EMI_BF4_NEUT_RADIUS = 0.8e-9f; // Radius of EMI + BF4 + Neutral
    const float EMIBF4_PERMITTIVITY = 12.8; // Relative permittivity EMIBF4
    const double EMI_MASS = 111; // Mass of EMI (in amu)
    const double BF4_MASS = 87; // Mass of BF4 (in amu)
    const double TEMPERATURE = 1100;
    
    const double K_CONST = (ELEM_CHARGE * ELEM_CHARGE) / (4 * PI * FREE_SPACE_PERMITTIVITY * AMU);
}

namespace InfoConstant {
    const int ARRAY_SIZE = 50000;
    const int BODY_STRUC_SIZE = 20;
    const int INFO_SIZE = 17;
    const int SAVE_SIZE_IONS = 9;
    const int SAVE_SIZE_NEUTRALS = 7;
    const float vz = 100;
}

namespace BHConstant {
    const float THETA = 0.1;    // Opening angle for Barnes-Hut algorithm
}

namespace FMMConstant {
    const int maxParticles = 10000000; // max of particles
    //const int numExpansions = 10;         // order of expansion in FMM
    const int maxP2PInteraction = 27; // max of P2P interacting boxes
    //const int maxM2LInteraction = 189;        // max of M2L interacting boxes
    const int maxM2LInteraction = 702; // = 27 * 27 - 27;
    //const int numRelativeBox = 512;        // max of relative box positioning
    const int targetBufferSize = 200000; // max of GPU target buffer
    const int sourceBufferSize = 100000; // max of GPU source buffer
    const int threadsPerBlockTypeA = 128; // size of GPU thread block P2P
    const int threadsPerBlockTypeB = 64; // size of GPU thread block M2L
    //const double eps = SOFTENING;       // double precision epsilon
    const int possibleRotations = 1155; /* number of possible rotations for a 11x11x11 box
                        with indices in the range [-5, 5], excluding constant factor differences */
    const int numRelativeBox = 1331; // = pow(11, 3);
    /*
     * More generally, minimum numRelativeBox = pow((4n - 1), 3), when considering nxnxn boxes;
     */
    const int minAvgPPB = 100; // Minimum average number of particles per box

    const int maxNumExpansions = 10;
    const int numExpansion2 = maxNumExpansions * maxNumExpansions;
    const int numExpansion4 = numExpansion2 * numExpansion2;
    const int numCoefficients = maxNumExpansions * (maxNumExpansions + 1) / 2;
    //const int numCoefficients = 55;
    const int DnmSize = (4 * numExpansion2 * maxNumExpansions - maxNumExpansions) / 3;


    const int windowSize = 5;
    const int crossoverSize = 10;
}

#endif /* CONSTANTS_H */

