/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "force_calculation_methods.h"

/**
 * Allocate memory associated with threading
 */
void ForceCalculationMethod::allocateThreadMemory_Base(){
    rangevalues = new vec4<int>[numThreads];
#ifdef _PTHREAD_H
    threads = new pthread_t[numThreads];
#endif
}

/**
 * Set the number of threads to be used and allocate necessary memory
 * @param n - number of threads
 */
void ForceCalculationMethod::setNumThreads(int n){
    numThreads = n;
    
    deallocateThreadMemory_Base();
    deallocateThreadMemory();
    allocateThreadMemory_Base();
    allocateThreadMemory();
}

/**
 * Deallocate memory associated with threading
 */
void ForceCalculationMethod::deallocateThreadMemory_Base(){
    if (rangevalues != nullptr) delete[] rangevalues;
    rangevalues = nullptr;
#ifdef _PTHREAD_H
    if (threads != nullptr) delete[] threads;
    threads = nullptr;
#endif
}

/**
 * Destructor for ForceCalculationMethod
 * Deallocate memory associated with threading
 */
ForceCalculationMethod::~ForceCalculationMethod(){
//    printf("Base class destructor called!\n");
    
    if (timeInfoFile != nullptr) fclose(timeInfoFile);
    if (L2NormFile != nullptr) fclose(L2NormFile);
    if (fullL2NormFile != nullptr) fclose(fullL2NormFile);
    timeInfoFile = nullptr;
    L2NormFile = nullptr;
    fullL2NormFile = nullptr;
    
    deallocateThreadMemory_Base();
//    deallocateThreadMemory();
    
}

/**
 * Print the number of threads used in the force calculations
 */
void ForceCalculationMethod::printNumThreads() {
    printf("Number of threads = %d.\n", numThreads);
}

/**
 * Set the index of the first particle
 * @param n - index of the first particle
 */
void ForceCalculationMethod::setFirstParticle(int n) {
    firstParticle = n;
}

/**
 * Set the index of the last particle
 * @param n - index of the last particle
 */
void ForceCalculationMethod::setLastParticle(int n) {
    lastParticle = n;
}

/**
 * Increment the index of the first particle
 * @param dN - amount to increment the index of the first particle by
 */
void ForceCalculationMethod::incrementFirstParticle(int dN) {
    firstParticle += dN;
}

/**
 * Increment the index of the first particle
 * @param dN - amount to increment the index of the last particle by
 */
void ForceCalculationMethod::incrementLastParticle(int dN) {
    lastParticle += dN;
}

/**
 * Set a pointer to a Grid class to use for symmetry
 * @param grid - pointer to Grid class to use for symmetry
 */
void ForceCalculationMethod::setGrid(Grid* grid) {
    this->grid = grid;
}

/**
 * Write the relevant information about L2-norms for a single iteration
 * @param numParticles - number of particles currently in the simulation
 * @param L2norm - L2-norm of the system of particles
 * @param maxL2norm - maximum L2-norm for a single particle for the iteration
 * @param maxL2normIndex - index of the particle with the maximum L2-norm
 */
void ForceCalculationMethod::writeL2NormInfo(int numParticles, float L2norm,
        float maxL2norm, int maxL2normIndex) {
    fwrite(&numParticles, sizeof (*&numParticles), 1, L2NormFile);
    fwrite(&L2norm, sizeof (*&L2norm), 1, L2NormFile);
    fwrite(&maxL2norm, sizeof (*&maxL2norm), 1, L2NormFile);
    fwrite(&maxL2normIndex, sizeof (*&maxL2normIndex), 1, L2NormFile);
}

/**
 * Initialize the timing file based on the current, time step, and sim time
 * @param folderName - folder in which to create timing file
 */
void ForceCalculationMethod::initializeTimingFile(const char* folderName) {
    setSaveName();
    float I0 = particleSystem->params.I0;
    float timeStep = particleSystem->params.timeStep;
    float simTime = particleSystem->params.simTime;
    char bufferTemp[512];
    int InA = roundf(I0 * 1e9);

    if (particleSystem->params.timeStep * 1e12 >= 0.999) {
        snprintf(bufferTemp, sizeof (char) * 512, "timing_%d_%.0fps_%.0fps", InA,
                timeStep * 1e12, simTime * 1e12); // creating filename depending on current
    } else {
        snprintf(bufferTemp, sizeof (char) * 512, "timing_%d_%.0ffs_%.0fps", InA,
                timeStep * 1e15, simTime * 1e12); // creating filename depending on current
    }

    char* buffer2;
    char* buffer;
    buffer2 = concat(bufferTemp, saveName);
    buffer = concat(folderName, buffer2);
    if (remove(buffer) == 0) {
        printf("\nOverwriting previous timing %s. Hope you saved the data!\n", saveName);
    } else {
        printf("\nCreating timing %s and storing your results.\n", saveName);
    }
    timeInfoFile = fopen(buffer, "ab");
    free(buffer);
    free(buffer2);
}

/**
 * Initialize the L2-norm file based on the current, time step, and sim time
 * @param folderName - folder in which to create L2-norm file
 */
void ForceCalculationMethod::initializeL2NormFile(const char* folderName) {
    setSaveName();
    float I0 = particleSystem->params.I0;
    float timeStep = particleSystem->params.timeStep;
    float simTime = particleSystem->params.simTime;
    char bufferTemp[512];
    int InA = roundf(I0 * 1e9);

    if (particleSystem->params.timeStep * 1e12 >= 0.999) {
        snprintf(bufferTemp, sizeof (char) * 512, "L2norm_%d_%.0fps_%.0fps", InA,
                timeStep * 1e12, simTime * 1e12); // creating filename depending on current
    } else {
        snprintf(bufferTemp, sizeof (char) * 512, "L2norm_%d_%.0ffs_%.0fps", InA,
                timeStep * 1e15, simTime * 1e12); // creating filename depending on current
    }

    char* buffer2;
    char* buffer;
    buffer2 = concat(bufferTemp, saveName);
    buffer = concat(folderName, buffer2);
    if (remove(buffer) == 0) {
        printf("\nOverwriting previous L2 norm %s. Hope you saved the data!\n", saveName);
    } else {
        printf("\nCreating L2 norm %s and storing your results.\n", saveName);
    }
    L2NormFile = fopen(buffer, "ab");
    free(buffer);
    free(buffer2);
}

/**
 * Initialize the full L2-norm file based on the current, time step, and sim time
 * @param folderName - folder in which to create full L2-norm file
 */
void ForceCalculationMethod::initializeFullL2NormFile(const char* folderName) {
    setSaveName();
    float I0 = particleSystem->params.I0;
    float timeStep = particleSystem->params.timeStep;
    float simTime = particleSystem->params.simTime;
    char bufferTemp[512];
    int InA = roundf(I0 * 1e9);

    if (particleSystem->params.timeStep * 1e12 >= 0.999) {
        snprintf(bufferTemp, sizeof (char) * 512, "FullL2norm_%d_%.0fps_%.0fps", InA,
                timeStep * 1e12, simTime * 1e12); // creating filename depending on current
    } else {
        snprintf(bufferTemp, sizeof (char) * 512, "FullL2norm_%d_%.0ffs_%.0fps", InA,
                timeStep * 1e15, simTime * 1e12); // creating filename depending on current
    }

    char* buffer2;
    char* buffer;
    buffer2 = concat(bufferTemp, saveName);
    buffer = concat(folderName, buffer2);
    if (remove(buffer) == 0) {
        printf("\nOverwriting previous full L2 norm %s. Hope you saved the data!\n", saveName);
    } else {
        printf("\nCreating full L2 norm %s and storing your results.\n", saveName);
    }
    fullL2NormFile = fopen(buffer, "ab");
    free(buffer);
    free(buffer2);
}

