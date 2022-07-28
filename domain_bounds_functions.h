/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   domain_bounds_functions.h
 * Author: Marshall
 *
 * Created on June 29, 2022, 2:49 PM
 */

#ifndef DOMAIN_BOUNDS_FUNCTIONS_H
#define DOMAIN_BOUNDS_FUNCTIONS_H

#include "auxiliary_functions.h"
#include "particles.h"

namespace DomainBoundsFunction {
    void isOutlier(ChargedParticles& particles, int* outlierCheck, int& numOutliers);
    void computeDomainBounds(ChargedParticles& particles, int* outlierCheck, vec3<int>& boxDimVals,
        vec3<float>& boxMin, vec3<float>& minOutliers, vec3<float>& maxOutliers,
        vec3<float>& dims, int numOutliers, int forceBoxDims, vec3<int> forcedBoxDimVals,
        DomainToggle domainToggle, int& avgPPB, vec3<int>& dimInd);
    
    void computeDomainBoundsRotations(ChargedParticles& particles, int* outlierCheck, vec3<int>& boxDimVals,
        vec3<float>& boxMin, vec3<float>& minOutliers, vec3<float>& maxOutliers,
        vec3<float>& dims, int numOutliers, int forceBoxDims, vec3<int> forcedBoxDimVals,
        DomainToggle domainToggle, int& avgPPB, vec3<int>& dimInd);
    
    void computeDomainBoundsTranslations(ChargedParticles& particles, int* outlierCheck, vec3<int>& boxDimVals,
        vec3<float>& boxMin, vec3<float>& minOutliers, vec3<float>& maxOutliers,
        vec3<float>& dims, int numOutliers, int forceBoxDims, vec3<int> forcedBoxDimVals,
        DomainToggle domainToggle, int& avgPPB, vec3<int>& dimInd, vec3<int> symmetryAxes);
    
    const float outlierFactor = 1.482602218505602; // = -1/(sqrt(2)*erfcinv(3/2));
    const float maxDeviation = 3;
    const float domainEps = 0.01; // total relative scaling factor for entire domain
}

#endif /* DOMAIN_BOUNDS_FUNCTIONS_H */

