/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "domain_bounds_functions.h"

/**
 * Determines which particles in a given system are outliers, based on the 
 *  median absolute deviation of each particle
 * @param particles - reference to ChargedParticles class being considered
 * @param outlierCheck - pointer to array indicating each particle's outlier status
 * @param numOutliers - reference to total number of outliers
 */
void DomainBoundsFunction::isOutlier(ChargedParticles& particles, 
        int* outlierCheck, int& numOutliers){
    
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    const float affectsAccelerationCheck = static_cast<float> (ParticleStatus::affectsAcceleration);

    int numParticles = particles.lastParticle - particles.firstParticle;
    numOutliers = 0;
    int i;
    float* absoluteDeviationsX = new float[numParticles];
    float* absoluteDeviationsY = new float[numParticles];
    float* absoluteDeviationsZ = new float[numParticles];
    float* tempX = new float[numParticles];
    float* tempY = new float[numParticles];
    float* tempZ = new float[numParticles];
    
    int numValidParticles = 0;

    for (i = 0; i < numParticles; i++) {
        if (particles.acceleration[i].w == inboundCheck ||
                particles.acceleration[i].w == affectsAccelerationCheck){
            tempX[numValidParticles] = particles.position[i + particles.firstParticle].x;
            tempY[numValidParticles] = particles.position[i + particles.firstParticle].y;
            tempZ[numValidParticles] = particles.position[i + particles.firstParticle].z;
            numValidParticles++;
        }
    }

    float medianX = findTrueMedian(tempX, numValidParticles);
    float medianY = findTrueMedian(tempY, numValidParticles);
    float medianZ = findTrueMedian(tempZ, numValidParticles);

    numValidParticles = 0;
    for (i = 0; i < numParticles; i++) {
        if (particles.acceleration[i].w == inboundCheck ||
                particles.acceleration[i].w == affectsAccelerationCheck){
            absoluteDeviationsX[numValidParticles] = fabsf(particles.position[i + particles.firstParticle].x - medianX);
            tempX[numValidParticles] = absoluteDeviationsX[numValidParticles];
            absoluteDeviationsY[numValidParticles] = fabsf(particles.position[i + particles.firstParticle].y - medianY);
            tempY[numValidParticles] = absoluteDeviationsY[numValidParticles];
            absoluteDeviationsZ[numValidParticles] = fabsf(particles.position[i + particles.firstParticle].z - medianZ);
            tempZ[numValidParticles] = absoluteDeviationsZ[numValidParticles];
            numValidParticles++;
        }
    }

    float madX = outlierFactor * findTrueMedian(tempX, numParticles);
    float madY = outlierFactor * findTrueMedian(tempY, numParticles);
    float madZ = outlierFactor * findTrueMedian(tempZ, numParticles);

    numValidParticles = 0;
    for (i = 0; i < numParticles; i++) {
        if (particles.acceleration[i].w == inboundCheck ||
                particles.acceleration[i].w == affectsAccelerationCheck){
            absoluteDeviationsX[numValidParticles] /= madX;
            absoluteDeviationsY[numValidParticles] /= madY;
            absoluteDeviationsZ[numValidParticles] /= madZ;
            if (absoluteDeviationsX[numValidParticles] >= maxDeviation || 
                    absoluteDeviationsY[numValidParticles] >= maxDeviation ||
                    absoluteDeviationsZ[numValidParticles] >= maxDeviation) {
                // Want to remove outliers in all dimensions, so don't reset outlierCheck
                
                // Used in acceleration calculation, but is an outlier
                outlierCheck[i] = 1;
                numOutliers++;
            } else {
                // Used in acceleration calculation and not an outlier
                outlierCheck[i] = 0;
            }
            numValidParticles++;
        }
        else {
            // Not used in acceleration calculation
            outlierCheck[i] = -1;
        }
    }

    delete[] tempX;
    delete[] tempY;
    delete[] tempZ;
    delete[] absoluteDeviationsX;
    delete[] absoluteDeviationsY;
    delete[] absoluteDeviationsZ;
}

/**
 * Computes the bounds of a domain of particles
 * @param particles - reference to chargedParticles class being considered
 * @param outlierCheck - pointer to array indicating each particle's outlier status
 * @param boxDimVals - reference to the dimensions of the base level cell structure
 * @param boxMin - reference to the minimum position values of the cell structure
 * @param minOutliers - reference to the minimum position values of all particles,
 *      including outliers
 * @param maxOutliers - reference to the maximum position values of all particles,
 *      including outliers
 * @param dims - reference to the dimensions of domain
 * @param numOutliers - number of outliers in the system
 * @param forceBoxDims - toggle for whether the dimensions of the base level cell
 *      structure will be fixed or not
 * @param forcedBoxDimVals - value to fix the dimensions of the base level cell
 *      structure, if desired
 * @param domainToggle - toggle for domain shape (cubic or rectangular prism)
 * @param avgPPB - reference to the average number of particles per box (PPB)
 * @param dimInd - reference to the indices of the directions of the domain
 */
void DomainBoundsFunction::computeDomainBounds(ChargedParticles& particles, int* outlierCheck, vec3<int>& boxDimVals,
        vec3<float>& boxMin, vec3<float>& minOutliers, vec3<float>& maxOutliers,
        vec3<float>& dims, int numOutliers, int forceBoxDims, vec3<int> forcedBoxDimVals,
        DomainToggle domainToggle, int& avgPPB, vec3<int>& dimInd){
    
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    const float affectsAccelerationCheck = static_cast<float> (ParticleStatus::affectsAcceleration);

//    float subtr = sqrt(PhysicalConstant::SOFTENING);
    float subtr = 0;
    int firstParticle = particles.firstParticle;
    int lastParticle = particles.lastParticle;
    int numParticles = 0;
    
    vec3<float> rootmax;
    
    int i = 0;
    while (outlierCheck[i] != 0){
        i++;
    }
    
    rootmax.x = particles.position[i + firstParticle].x;
    rootmax.y = particles.position[i + firstParticle].y;
    rootmax.z = particles.position[i + firstParticle].z;
    boxMin = rootmax;
    minOutliers = rootmax;
    maxOutliers = rootmax;
    
    for (i = firstParticle; i < lastParticle; i++){
        if (particles.acceleration[i].w == inboundCheck ||
                particles.acceleration[i].w == affectsAccelerationCheck){
            if (particles.position[i].x < boxMin.x && outlierCheck[i - firstParticle] == 0){
                boxMin.x = particles.position[i].x;
            }
            if (particles.position[i].x > rootmax.x && outlierCheck[i - firstParticle] == 0){
                rootmax.x = particles.position[i].x;
            }
            if (particles.position[i].y < boxMin.y && outlierCheck[i - firstParticle] == 0){
                boxMin.y = particles.position[i].y;
            }
            if (particles.position[i].y > rootmax.y && outlierCheck[i - firstParticle] == 0){
                rootmax.y = particles.position[i].y;
            }
            if (particles.position[i].z < boxMin.z && outlierCheck[i - firstParticle] == 0){
                boxMin.z = particles.position[i].z;
            }
            if (particles.position[i].z > rootmax.z && outlierCheck[i - firstParticle] == 0){
                rootmax.z = particles.position[i].z;
            }

            if (particles.position[i].x < minOutliers.x){
                minOutliers.x = particles.position[i].x;
            }
            if (particles.position[i].x > maxOutliers.x){
                maxOutliers.x = particles.position[i].x;
            }
            if (particles.position[i].y < minOutliers.y){
                minOutliers.y = particles.position[i].y;
            }
            if (particles.position[i].y > maxOutliers.y){
                maxOutliers.y = particles.position[i].y;
            }
            if (particles.position[i].z < minOutliers.z){
                minOutliers.z = particles.position[i].z;
            }
            if (particles.position[i].z > maxOutliers.z){
                maxOutliers.z = particles.position[i].z;
            }
            numParticles++;
        }
    }
    
//    minOutliers = minOutliers - subtr;
//    maxOutliers = maxOutliers + subtr;
    
    boxMin = boxMin - subtr;
    dims = rootmax - boxMin + subtr;
    
    boxDimVals[0] = 1;
    boxDimVals[1] = 1;
    boxDimVals[2] = 1;
    
    if (domainToggle == DomainToggle::cubicCells){
        float dimDiff[3] = {dims.x, dims.y, dims.z};
        dimInd = {0, 1, 2};
        int j, tempInd;
        float tempDim;

        // Sort dimDiff and dimInd, such that dimDiff is ordered from smallest to largest
        i = 0;
        while (i < 3) {
            j = i;
            while ((j > 0) && dimDiff[j - 1] > dimDiff[j]) {
                tempDim = dimDiff[j - 1];
                tempInd = dimInd[j - 1];
                dimDiff[j - 1] = dimDiff[j];
                dimInd[j - 1] = dimInd[j];
                dimDiff[j] = tempDim;
                dimInd[j] = tempInd;
                j--;
            }
            i++;
        }

        float rootBoxSize = 0;
        i = -1;
        do {
            i++;
            rootBoxSize = dimDiff[i];
            boxDimVals[dimInd[i]] = 1;
        } while (rootBoxSize == 0 && i < 2);

        rootBoxSize = std::max(rootBoxSize, (float)PhysicalConstant::SOFTENING);

        float ratio;
        int numDomains;
        int boxSizeFound = 0;
        int minPPB = FMMConstant::minAvgPPB;

        do {
            for (j = i + 1; j < 3; j++) {
                ratio = dimDiff[j] / rootBoxSize;
                if (ratio - floor(ratio) >= 0.5) {
                    boxDimVals[dimInd[j]] = ceil(ratio);
                } else {
                    if (ratio < 2) {
                        if (j == 0) boxDimVals[dimInd[j]] = 1;
                        else boxDimVals[dimInd[j]] = boxDimVals[dimInd[j - 1]];
                    } else {
                        boxDimVals[dimInd[j]] = floor(ratio);
                    }
                    if (ratio > 0.5 || j == 2) {
                        rootBoxSize = dimDiff[j] / boxDimVals[dimInd[j]];
                    }
                }
            }
            numDomains = boxDimVals[0] * boxDimVals[1] * boxDimVals[2];
            avgPPB = (numParticles - numOutliers + numDomains - 1) / numDomains;
            if (rootBoxSize == dimDiff[2] || avgPPB >= minPPB || rootBoxSize == (float)PhysicalConstant::SOFTENING) {
                // If no splitting required or average particles per box is above critical value, stop
                boxSizeFound = 1;
            } else {
                // Else, increase rootBoxSize by a factor of two
                rootBoxSize *= 2;
            }
        } while (boxSizeFound == 0);

        if (forceBoxDims == true){
            boxDimVals[0] = forcedBoxDimVals.x;
            boxDimVals[1] = forcedBoxDimVals.y;
            boxDimVals[2] = forcedBoxDimVals.z;
            numDomains = boxDimVals[0] * boxDimVals[1] * boxDimVals[2];
            rootBoxSize = std::max(dims.x / boxDimVals[0], dims.y / boxDimVals[1]);
            rootBoxSize = std::max(rootBoxSize, dims.z / boxDimVals[2]);
            rootBoxSize = std::max(rootBoxSize, (float)PhysicalConstant::SOFTENING);
            if (boxDimVals[0] >= boxDimVals[1] && boxDimVals[0] >= boxDimVals[2]){
                dimInd[2] = 0;
            }
            else if (boxDimVals[1] >= boxDimVals[0] && boxDimVals[1] >= boxDimVals[2]){
                dimInd[2] = 1;
            }
            else if (boxDimVals[2] >= boxDimVals[0] && boxDimVals[2] >= boxDimVals[1]){
                dimInd[2] = 2;
            }
        }

        rootBoxSize *= (1 + domainEps / boxDimVals[dimInd[2]]);

        boxMin.x = (rootmax.x + subtr + boxMin.x - rootBoxSize * boxDimVals[0]) / 2;
        boxMin.y = (rootmax.y + subtr + boxMin.y - rootBoxSize * boxDimVals[1]) / 2;
        boxMin.z = (rootmax.z + subtr + boxMin.z - rootBoxSize * boxDimVals[2]) / 2;
        
        dims = rootBoxSize;
    }
}

/**
 * Computes the bounds of a domain of particles for a rotationally symmetric system
 * @param particles - reference to chargedParticles class being considered
 * @param outlierCheck - pointer to array indicating each particle's outlier status
 * @param boxDimVals - reference to the dimensions of the base level cell structure
 * @param boxMin - reference to the minimum position values of the cell structure
 * @param minOutliers - reference to the minimum position values of all particles,
 *      including outliers
 * @param maxOutliers - reference to the maximum position values of all particles,
 *      including outliers
 * @param dims - reference to the dimensions of domain
 * @param numOutliers - number of outliers in the system
 * @param forceBoxDims - toggle for whether the dimensions of the base level cell
 *      structure will be fixed or not
 * @param forcedBoxDimVals - value to fix the dimensions of the base level cell
 *      structure, if desired
 * @param domainToggle - toggle for domain shape (cubic or rectangular prism)
 * @param avgPPB - reference to the average number of particles per box (PPB)
 * @param dimInd - reference to the indices of the directions of the domain
 */
void DomainBoundsFunction::computeDomainBoundsRotations(ChargedParticles& particles, int* outlierCheck, vec3<int>& boxDimVals,
        vec3<float>& boxMin, vec3<float>& minOutliers, vec3<float>& maxOutliers,
        vec3<float>& dims, int numOutliers, int forceBoxDims, vec3<int> forcedBoxDimVals,
        DomainToggle domainToggle, int& avgPPB, vec3<int>& dimInd){
    
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    const float affectsAccelerationCheck = static_cast<float> (ParticleStatus::affectsAcceleration);

//    float subtr = sqrt(PhysicalConstant::SOFTENING);
    float subtr = 0;
    int firstParticle = particles.firstParticle;
    int lastParticle = particles.lastParticle;
    int numParticles = 0;
    
    vec3<float> rootmax;
    
    int i = 0;
    while (outlierCheck[i] != 0){
        i++;
    }
    
    rootmax.x = particles.position[i + firstParticle].x;
    rootmax.y = particles.position[i + firstParticle].y;
    rootmax.z = particles.position[i + firstParticle].z;
    boxMin.x = 0.0f;
    boxMin.y = 0.0f;
    boxMin.z = rootmax.z;
    minOutliers = rootmax;
    maxOutliers = rootmax;
    
    for (i = firstParticle; i < lastParticle; i++){
        if (particles.acceleration[i].w == inboundCheck ||
                particles.acceleration[i].w == affectsAccelerationCheck){
//            if (particles.position[i].x < boxMin.x && outlierCheck[i - firstParticle] == 0){
//                boxMin.x = particles.position[i].x;
//            }
            if (particles.position[i].x > rootmax.x && outlierCheck[i - firstParticle] == 0){
                rootmax.x = particles.position[i].x;
            }
//            if (particles.position[i].y < boxMin.y && outlierCheck[i - firstParticle] == 0){
//                boxMin.y = particles.position[i].y;
//            }
            if (particles.position[i].y > rootmax.y && outlierCheck[i - firstParticle] == 0){
                rootmax.y = particles.position[i].y;
            }
            if (particles.position[i].z < boxMin.z && outlierCheck[i - firstParticle] == 0){
                boxMin.z = particles.position[i].z;
            }
            if (particles.position[i].z > rootmax.z && outlierCheck[i - firstParticle] == 0){
                rootmax.z = particles.position[i].z;
            }

            if (particles.position[i].x < minOutliers.x){
                minOutliers.x = particles.position[i].x;
            }
            if (particles.position[i].x > maxOutliers.x){
                maxOutliers.x = particles.position[i].x;
            }
            if (particles.position[i].y < minOutliers.y){
                minOutliers.y = particles.position[i].y;
            }
            if (particles.position[i].y > maxOutliers.y){
                maxOutliers.y = particles.position[i].y;
            }
            if (particles.position[i].z < minOutliers.z){
                minOutliers.z = particles.position[i].z;
            }
            if (particles.position[i].z > maxOutliers.z){
                maxOutliers.z = particles.position[i].z;
            }
            numParticles++;
        }
    }
    
//    minOutliers = minOutliers - subtr;
//    maxOutliers = maxOutliers + subtr;
    
    rootmax += subtr;
    rootmax.x = std::max(rootmax.x, rootmax.y);
    rootmax.y = rootmax.x;
    boxMin.x = -rootmax.x;
    boxMin.y = -rootmax.y;
    boxMin.z = boxMin.z - subtr;
    
    dims.x = rootmax.x * 2;
    dims.y = rootmax.y * 2;
    dims.z = rootmax.z - boxMin.z;
    
    boxDimVals[0] = 1;
    boxDimVals[1] = 1;
    boxDimVals[2] = 1;
    
    if (domainToggle == DomainToggle::cubicCells){
        float dimDiff[3] = {dims.x, dims.y, dims.z};
        dimInd = {0, 1, 2};
        int j, tempInd;
        float tempDim;

        // Sort dimDiff and dimInd, such that dimDiff is ordered from smallest to largest
        i = 0;
        while (i < 3) {
            j = i;
            while ((j > 0) && dimDiff[j - 1] > dimDiff[j]) {
                tempDim = dimDiff[j - 1];
                tempInd = dimInd[j - 1];
                dimDiff[j - 1] = dimDiff[j];
                dimInd[j - 1] = dimInd[j];
                dimDiff[j] = tempDim;
                dimInd[j] = tempInd;
                j--;
            }
            i++;
        }

        float rootBoxSize = 0;
        i = -1;
        do {
            i++;
            rootBoxSize = dimDiff[i];
            boxDimVals[dimInd[i]] = 1;
        } while (rootBoxSize == 0 && i < 2);

        rootBoxSize = std::max(rootBoxSize, (float)PhysicalConstant::SOFTENING);

        float ratio;
        int numDomains;
        int boxSizeFound = 0;
        int minPPB = FMMConstant::minAvgPPB;

        do {
            for (j = i + 1; j < 3; j++) {
                ratio = dimDiff[j] / rootBoxSize;
                if (ratio - floor(ratio) >= 0.5) {
                    boxDimVals[dimInd[j]] = ceil(ratio);
                } else {
                    if (ratio < 2) {
                        if (j == 0) boxDimVals[dimInd[j]] = 1;
                        else boxDimVals[dimInd[j]] = boxDimVals[dimInd[j - 1]];
                    } else {
                        boxDimVals[dimInd[j]] = floor(ratio);
                    }
                    if (ratio > 0.5 || j == 2) {
                        rootBoxSize = dimDiff[j] / boxDimVals[dimInd[j]];
                    }
                }
            }
            numDomains = boxDimVals[0] * boxDimVals[1] * boxDimVals[2];
            avgPPB = (numParticles - numOutliers + numDomains - 1) / numDomains;
            if (rootBoxSize == dimDiff[2] || avgPPB >= minPPB || rootBoxSize == (float)PhysicalConstant::SOFTENING) {
                // If no splitting required or average particles per box is above critical value, stop
                boxSizeFound = 1;
            } else {
                // Else, increase rootBoxSize by a factor of two
                rootBoxSize *= 2;
            }
        } while (boxSizeFound == 0);

        if (forceBoxDims == true){
            boxDimVals[0] = forcedBoxDimVals.x;
            boxDimVals[1] = forcedBoxDimVals.y;
            boxDimVals[2] = forcedBoxDimVals.z;
            numDomains = boxDimVals[0] * boxDimVals[1] * boxDimVals[2];
            rootBoxSize = std::max(dims.x / boxDimVals[0], dims.y / boxDimVals[1]);
            rootBoxSize = std::max(rootBoxSize, dims.z / boxDimVals[2]);
            rootBoxSize = std::max(rootBoxSize, (float)PhysicalConstant::SOFTENING);
            if (boxDimVals[0] >= boxDimVals[1] && boxDimVals[0] >= boxDimVals[2]){
                dimInd[2] = 0;
            }
            else if (boxDimVals[1] >= boxDimVals[0] && boxDimVals[1] >= boxDimVals[2]){
                dimInd[2] = 1;
            }
            else if (boxDimVals[2] >= boxDimVals[0] && boxDimVals[2] >= boxDimVals[1]){
                dimInd[2] = 2;
            }
        }

        rootBoxSize *= (1 + domainEps / boxDimVals[dimInd[2]]);

        boxMin.x = (rootmax.x + boxMin.x - rootBoxSize * boxDimVals[0]) / 2;
        boxMin.y = (rootmax.y + boxMin.y - rootBoxSize * boxDimVals[1]) / 2;
        boxMin.z = (rootmax.z + boxMin.z - rootBoxSize * boxDimVals[2]) / 2;
        
        dims = rootBoxSize;
    }
}

/**
 * Computes the bounds of a domain of particles for a rotationally symmetric system
 * @param particles - reference to chargedParticles class being considered
 * @param outlierCheck - pointer to array indicating each particle's outlier status
 * @param boxDimVals - reference to the dimensions of the base level cell structure
 * @param boxMin - reference to the minimum position values of the cell structure
 * @param minOutliers - reference to the minimum position values of all particles,
 *      including outliers
 * @param maxOutliers - reference to the maximum position values of all particles,
 *      including outliers
 * @param dims - reference to the dimensions of domain
 * @param numOutliers - number of outliers in the system
 * @param forceBoxDims - toggle for whether the dimensions of the base level cell
 *      structure will be fixed or not
 * @param forcedBoxDimVals - value to fix the dimensions of the base level cell
 *      structure, if desired
 * @param domainToggle - toggle for domain shape (cubic or rectangular prism)
 * @param avgPPB - reference to the average number of particles per box (PPB)
 * @param dimInd - reference to the indices of the directions of the domain
 */
void DomainBoundsFunction::computeDomainBoundsTranslations(ChargedParticles& particles, int* outlierCheck, vec3<int>& boxDimVals,
        vec3<float>& boxMin, vec3<float>& minOutliers, vec3<float>& maxOutliers,
        vec3<float>& dims, int numOutliers, int forceBoxDims, vec3<int> forcedBoxDimVals,
        DomainToggle domainToggle, int& avgPPB, vec3<int>& dimInd, vec3<int> symmetryAxes){
    
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    const float affectsAccelerationCheck = static_cast<float> (ParticleStatus::affectsAcceleration);

//    float subtr = sqrt(PhysicalConstant::SOFTENING);
    float subtr = 0;
    float rootBoxSize = 0;
    int firstParticle = particles.firstParticle;
    int lastParticle = particles.lastParticle;
    int numParticles = 0;
    vec3<float> dimDiff;
    int numDomains;
    
    vec3<float> rootmax;
    
    int i = 0;
    while (outlierCheck[i] != 0){
        i++;
    }
    
    rootmax.x = particles.position[i + firstParticle].x;
    rootmax.y = particles.position[i + firstParticle].y;
    rootmax.z = particles.position[i + firstParticle].z;
    boxMin.x = 0.0f;
    boxMin.y = 0.0f;
    boxMin.z = rootmax.z;
    minOutliers = rootmax;
    maxOutliers = rootmax;
    
    for (i = firstParticle; i < lastParticle; i++){
        if (particles.acceleration[i].w == inboundCheck ||
                particles.acceleration[i].w == affectsAccelerationCheck){
            if (particles.position[i].x < boxMin.x && outlierCheck[i - firstParticle] == 0){
                boxMin.x = particles.position[i].x;
            }
            if (particles.position[i].x > rootmax.x && outlierCheck[i - firstParticle] == 0){
                rootmax.x = particles.position[i].x;
            }
            if (particles.position[i].y < boxMin.y && outlierCheck[i - firstParticle] == 0){
                boxMin.y = particles.position[i].y;
            }
            if (particles.position[i].y > rootmax.y && outlierCheck[i - firstParticle] == 0){
                rootmax.y = particles.position[i].y;
            }
            if (particles.position[i].z < boxMin.z && outlierCheck[i - firstParticle] == 0){
                boxMin.z = particles.position[i].z;
            }
            if (particles.position[i].z > rootmax.z && outlierCheck[i - firstParticle] == 0){
                rootmax.z = particles.position[i].z;
            }

            if (particles.position[i].x < minOutliers.x){
                minOutliers.x = particles.position[i].x;
            }
            if (particles.position[i].x > maxOutliers.x){
                maxOutliers.x = particles.position[i].x;
            }
            if (particles.position[i].y < minOutliers.y){
                minOutliers.y = particles.position[i].y;
            }
            if (particles.position[i].y > maxOutliers.y){
                maxOutliers.y = particles.position[i].y;
            }
            if (particles.position[i].z < minOutliers.z){
                minOutliers.z = particles.position[i].z;
            }
            if (particles.position[i].z > maxOutliers.z){
                maxOutliers.z = particles.position[i].z;
            }
            numParticles++;
        }
    }
    
    dimDiff = rootmax - boxMin;
    
    // Ensure that the domain size is the same in all dimensions, and set rootBoxSize
    for (i = 0; i < 3; i++){
        if (rootBoxSize == 0){
            rootBoxSize = dims[i];
        }
        else {
            if (symmetryAxes[i] != 0){
                assert(rootBoxSize == dims[i]);
            }
        }
    }
    
    if (rootBoxSize == 0){
        if (dimDiff.x != 0){
            rootBoxSize = dimDiff.x;
            if (dimDiff.y != 0){
                rootBoxSize = std::min(rootBoxSize, dimDiff.y);
                if (dimDiff.z != 0){
                    rootBoxSize = std::min(rootBoxSize, dimDiff.z);
                }
            }
            else if (dimDiff.z != 0){
                rootBoxSize = std::min(rootBoxSize, dimDiff.z);
            }
        }
        else if (dimDiff.y != 0){
            rootBoxSize = std::min(rootBoxSize, dimDiff.y);
            if (dimDiff.z != 0){
                rootBoxSize = std::min(rootBoxSize, dimDiff.z);
            }
        }
        else if (dimDiff.z != 0){
            rootBoxSize = std::min(rootBoxSize, dimDiff.z);
        }
    }
    
    if (symmetryAxes.x == 1){
        boxMin.x = -rootBoxSize / 2;
        boxDimVals[0] = 1;
    }
    else {
        boxDimVals[0] = (dimDiff.x == 0) ? 1 : ceil(dimDiff.x / rootBoxSize);
        boxMin.x = (rootmax.x + boxMin.x - rootBoxSize * boxDimVals[0]) / 2;
    }
    
    if (symmetryAxes.y == 1){
        boxMin.y = -rootBoxSize / 2;
        boxDimVals[1] = 1;
    }
    else {
        boxDimVals[1] = (dimDiff.y == 0) ? 1 : ceil(dimDiff.y / rootBoxSize);
        boxMin.y = (rootmax.y + boxMin.y - rootBoxSize * boxDimVals[0]) / 2;
    }
    
    if (symmetryAxes.z == 1){
        boxMin.z = -rootBoxSize / 2;
        boxDimVals[2] = 1;
    }
    else {
        boxDimVals[2] = (dimDiff.z == 0) ? 1 : ceil(dimDiff.z / rootBoxSize);
        boxMin.z = (rootmax.z + boxMin.z - rootBoxSize * boxDimVals[0]) / 2;
    }
    
    numDomains = boxDimVals[0] * boxDimVals[1] * boxDimVals[2];
    avgPPB = (numParticles - numOutliers + numDomains - 1) / numDomains;
    
    dims = rootBoxSize;
}