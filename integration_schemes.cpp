/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "integration_schemes.h"

/**
 * Construct Integrator class by setting the particle system and force calculation method
 * @param particleSystem - class containing information about particles
 * @param forceMethod - class containing the force calculation method to be used
 */
Integrator::Integrator(ParticleSystem& particleSystem, ForceCalculationMethod& forceMethod) : 
particleSystem{particleSystem}, forceMethod{forceMethod} {
    if (particleSystem.params.timerName != nullptr && 
            particleSystem.params.timerName[0] != '\0'){
        fullSimTimer.setEventName(particleSystem.params.timerName);
        fullSimTimer.setDetailName("full simulation");
        intervalTimer.setEventName(particleSystem.params.timerName);
        intervalTimer.setDetailName("interval");
        iterTimer.setEventName(particleSystem.params.timerName);
        iterTimer.setDetailName("iteration");
    }
    else{
        forceMethod.setTimerName();
        fullSimTimer.setEventName(forceMethod.timerName);
        fullSimTimer.setDetailName("full simulation");
        intervalTimer.setEventName(forceMethod.timerName);
        intervalTimer.setDetailName("interval");
        iterTimer.setEventName(forceMethod.timerName);
        iterTimer.setDetailName("iteration");
    }
}

/**
 * Destructor for Integrator
 * Deallocate memory and close energy file
 */
Integrator::~Integrator(){
    if (energyFile != nullptr) fclose(energyFile);
    energyFile = nullptr;
    
    if (comparisonAcceleration != nullptr) delete[] comparisonAcceleration;
    comparisonAcceleration = nullptr;
    
    if (initialPositions != nullptr) delete[] initialPositions;
    initialPositions = nullptr;
}

/**
 * Set additional force methods to compare to the baseline
 * @param forceMethodComparisons - array of classes containing the force 
 *      calculation methods to be used as comparisons
 */
void Integrator::setComparisons(Array<ForceCalculationMethod*>& forceMethodComparisons) {
    this->forceMethodComparisons = forceMethodComparisons;
    fullSimTimerComparisons = Array<Timer>(forceMethodComparisons.numElements);
    intervalTimerComparisons = Array<Timer>(forceMethodComparisons.numElements);
    iterTimerComparisons = Array<Timer>(forceMethodComparisons.numElements);
    for (int i = 0; i < forceMethodComparisons.numElements; i++){
        forceMethodComparisons[i]->setTimerName();
        fullSimTimerComparisons[i].setEventName(forceMethodComparisons[i]->timerName);
        fullSimTimerComparisons[i].setDetailName("full simulation");
        intervalTimerComparisons[i].setEventName(forceMethodComparisons[i]->timerName);
        intervalTimerComparisons[i].setDetailName("interval");
        iterTimerComparisons[i].setEventName(forceMethodComparisons[i]->timerName);
        iterTimerComparisons[i].setDetailName("iteration");
    }
    
    if (forceMethodComparisons.numElements > 0){
        comparisonAcceleration = new vec4<float>[particleSystem.chargedParticles.maxNumParticles];
    }
}

/**
 * Toggle whether an energy conservation check will be performed 
 * @param n - true if energy conservation check will be performed, false otherwise
 * @param initSteps - number of steps in the numerical integration for the initial
 *      energy of the particles due to the external field
 * @param normalSteps - number of steps in the numerical integration for the energy
 *      of the particles due to the external field for regular time steps of the simulation
 * @param folderName - folder in which to create energy file
 */
void Integrator::toggleEnergy(int n, int initSteps, int normalSteps, const char* folderName){
    energyToggle = n;
    if (energyToggle == true){
        numEnergySteps_initialization = initSteps;
        numEnergySteps = normalSteps;
        initializeEnergyFile(folderName);
    }
    else{
        if (energyFile != nullptr) fclose(energyFile);
        energyFile = nullptr;
        energyToggle = false;
    }
}

/**
 * Initialize the energy file based on the current, time step, and sim time
 * @param folderName - folder in which to create energy file
 */
void Integrator::initializeEnergyFile(const char* folderName){
    float I0 = particleSystem.params.I0;
    float timeStep = particleSystem.params.timeStep;
    float simTime = particleSystem.params.simTime;
    char bufferTemp[512];
    int InA = roundf(I0 * 1e9);

    if (particleSystem.params.timeStep * 1e12 >= 0.999) {
        snprintf(bufferTemp, sizeof (char) * 512, "energy_%d_%.0fps_%.0fps", InA,
                timeStep * 1e12, simTime * 1e12); // creating filename depending on current
    } else {
        snprintf(bufferTemp, sizeof (char) * 512, "energy_%d_%.0ffs_%.0fps", InA,
                timeStep * 1e15, simTime * 1e12); // creating filename depending on current
    }

    char* buffer2;
    char* buffer;
    buffer2 = concat(bufferTemp, energySaveName);
    buffer = concat(folderName, buffer2);
    if (remove(buffer) == 0) {
        printf("\nOverwriting previous energy %s. Hope you saved the data!\n", energySaveName);
    } else {
        printf("\nCreating energy %s and storing your results.\n", energySaveName);
    }
    energyFile = fopen(buffer, "ab");
    free(buffer);
    free(buffer2);
}

/**
 * Set the index of the first charged particle
 * @param n - index of the first charged particle
 */
void Integrator::setFirstChargedParticle(int n){
    particleSystem.chargedParticles.setFirstParticle(n);
    forceMethod.setFirstParticle(n);
    for (int i = 0; i < forceMethodComparisons.numElements; i++){
        forceMethodComparisons[i]->setFirstParticle(n);
    }
    firstParticle = n;
}

/**
 * Set the index of the first charged particle to the smallest index containing
 *      a particle still considered in the simulation
 */
void Integrator::resetFirstChargedParticle(){
    particleSystem.resetFirstChargedParticle();
    setFirstChargedParticle(particleSystem.chargedParticles.firstParticle);
}

/**
 * Set the index of the last charged particle
 * @param n - index of the last charged particle
 */
void Integrator::setLastChargedParticle(int n){
    particleSystem.chargedParticles.setLastParticle(n);
    forceMethod.setLastParticle(n);
    for (int i = 0; i < forceMethodComparisons.numElements; i++){
        forceMethodComparisons[i]->setLastParticle(n);
    }
    lastParticle = n;
}

/**
 * Increment the index of the first charged particle
 * @param dN - amount to increment the index of the first charged particle by
 */
void Integrator::incrementFirstChargedParticle(int dN) {
    forceMethod.incrementFirstParticle(dN);
    particleSystem.chargedParticles.incrementFirstParticle(dN);
    for (int i = 0; i < forceMethodComparisons.numElements; i++){
        forceMethodComparisons[i]->incrementFirstParticle(dN);
    }
    firstParticle += dN;
}

/**
 * Increment the index of the last charged particle
 * @param dN - amount to increment the index of the last charged particle by
 */
void Integrator::incrementLastChargedParticle(int dN){
    forceMethod.incrementLastParticle(dN);
    particleSystem.chargedParticles.incrementLastParticle(dN);
    for (int i = 0; i < forceMethodComparisons.numElements; i++){
        forceMethodComparisons[i]->incrementLastParticle(dN);
    }
    lastParticle += dN;
}

/**
 * Saves the current state data of all charged particles to a given file
 * @param filename - file to write state data to
 */
void Integrator::saveFullCurrentConditions(FILE* filename){
    int gg;
    uint32_t temp;
    
    switch (particleSystem.params.saveDataToggle) {
        case (SaveDataOptions::original):
            for (gg = 0; gg < particleSystem.chargedParticles.maxNumParticles; gg++) {
//                if (particleSystem.chargedParticles.acceleration[gg].w == true) {
//                if (particleSystem.chargedParticles.acceleration[gg].w == inboundCheck ||
//                        particleSystem.chargedParticles.acceleration[gg].w == affectsAccelerationCheck){
                    fwrite(&particleSystem.chargedParticles.position[gg].x, sizeof (*&particleSystem.chargedParticles.position[gg].x), 1, filename);
                    fwrite(&particleSystem.chargedParticles.position[gg].y, sizeof (*&particleSystem.chargedParticles.position[gg].y), 1, filename);
                    fwrite(&particleSystem.chargedParticles.position[gg].z, sizeof (*&particleSystem.chargedParticles.position[gg].z), 1, filename);
                    fwrite(&particleSystem.chargedParticles.velocity[gg].x, sizeof (*&particleSystem.chargedParticles.velocity[gg].x), 1, filename);
                    fwrite(&particleSystem.chargedParticles.velocity[gg].y, sizeof (*&particleSystem.chargedParticles.velocity[gg].y), 1, filename);
                    fwrite(&particleSystem.chargedParticles.velocity[gg].z, sizeof (*&particleSystem.chargedParticles.velocity[gg].z), 1, filename);
                    
                    fwrite(&particleSystem.chargedParticles.Efield1[gg].x, sizeof (*&particleSystem.chargedParticles.Efield1[gg].x), 1, filename);
                    fwrite(&particleSystem.chargedParticles.Efield1[gg].y, sizeof (*&particleSystem.chargedParticles.Efield1[gg].y), 1, filename);
                    fwrite(&particleSystem.chargedParticles.Efield1[gg].z, sizeof (*&particleSystem.chargedParticles.Efield1[gg].z), 1, filename);
                    fwrite(&particleSystem.chargedParticles.Efield2[gg].x, sizeof (*&particleSystem.chargedParticles.Efield2[gg].x), 1, filename);
                    fwrite(&particleSystem.chargedParticles.Efield2[gg].y, sizeof (*&particleSystem.chargedParticles.Efield2[gg].y), 1, filename);
                    fwrite(&particleSystem.chargedParticles.Efield2[gg].z, sizeof (*&particleSystem.chargedParticles.Efield2[gg].z), 1, filename);
                    
                    temp = particleSystem.chargedParticles.velocity[gg].w;
                    fwrite(&temp, sizeof(*&temp), 1, filename);
//                    fwrite(&particleSystem.chargedParticles.velocity[gg].w, sizeof (*&particleSystem.chargedParticles.velocity[gg].w), 1, filename);
                    fwrite(&gg, sizeof (*&gg), 1, filename);
                    fwrite(&lastParticle, sizeof (*&lastParticle), 1, filename);
//                }
            }
            break;
        case (SaveDataOptions::newOriginal):
            temp = -1;
            fwrite(&temp, sizeof(*&temp), 1, filename);
            fwrite(&lastParticle, sizeof(*&lastParticle), 1, filename);
            fwrite(&currentTime, sizeof(*&currentTime), 1, filename);
            for (gg = 0; gg < particleSystem.chargedParticles.maxNumParticles; gg++){
//                if (particleSystem.chargedParticles.acceleration[gg].w == inboundCheck ||
//                        particleSystem.chargedParticles.acceleration[gg].w == affectsAccelerationCheck){
                    temp = particleSystem.chargedParticles.velocity[gg].w;
                    fwrite(&gg, sizeof(*&gg), 1, filename);
                    fwrite(&temp, sizeof(*&temp), 1, filename);
                    
                    fwrite(&particleSystem.chargedParticles.position[gg].x, sizeof(*&particleSystem.chargedParticles.position[gg].x), 1, filename);
                    fwrite(&particleSystem.chargedParticles.position[gg].y, sizeof(*&particleSystem.chargedParticles.position[gg].y), 1, filename);
                    fwrite(&particleSystem.chargedParticles.position[gg].z, sizeof(*&particleSystem.chargedParticles.position[gg].z), 1, filename);
                    fwrite(&particleSystem.chargedParticles.velocity[gg].x, sizeof(*&particleSystem.chargedParticles.velocity[gg].x), 1, filename);
                    fwrite(&particleSystem.chargedParticles.velocity[gg].y, sizeof(*&particleSystem.chargedParticles.velocity[gg].y), 1, filename);
                    fwrite(&particleSystem.chargedParticles.velocity[gg].z, sizeof(*&particleSystem.chargedParticles.velocity[gg].z), 1, filename);
                    
                    fwrite(&particleSystem.chargedParticles.Efield1[gg].x, sizeof(*&particleSystem.chargedParticles.Efield1[gg].x), 1, filename);
                    fwrite(&particleSystem.chargedParticles.Efield1[gg].y, sizeof(*&particleSystem.chargedParticles.Efield1[gg].y), 1, filename);
                    fwrite(&particleSystem.chargedParticles.Efield1[gg].z, sizeof(*&particleSystem.chargedParticles.Efield1[gg].z), 1, filename);
                    fwrite(&particleSystem.chargedParticles.Efield2[gg].x, sizeof(*&particleSystem.chargedParticles.Efield2[gg].x), 1, filename);
                    fwrite(&particleSystem.chargedParticles.Efield2[gg].y, sizeof(*&particleSystem.chargedParticles.Efield2[gg].y), 1, filename);
                    fwrite(&particleSystem.chargedParticles.Efield2[gg].z, sizeof(*&particleSystem.chargedParticles.Efield2[gg].z), 1, filename);
//                }
            }
            break;
    }
}

/**
 * Saves data of valid particles at a specific time step of a simulation
 */
void Integrator::saveData() {
    
    static int savedSnapNum = 0;
    float dt;
    int i = 0;
    int gg;
    uint32_t temp;
    
    const float inboundCheck = static_cast<float>(ParticleStatus::inbound);
    const float affectsAccelerationCheck = static_cast<float>(ParticleStatus::affectsAcceleration);
    
    switch (particleSystem.params.saveDataToggle) {
        case (SaveDataOptions::original):
            for (gg = firstParticle; gg < lastParticle; gg++) {
//                if (particleSystem.chargedParticles.acceleration[gg].w == true) {
                if (particleSystem.chargedParticles.acceleration[gg].w == inboundCheck ||
                        particleSystem.chargedParticles.acceleration[gg].w == affectsAccelerationCheck){
                    fwrite(&particleSystem.chargedParticles.position[gg].x, sizeof (*&particleSystem.chargedParticles.position[gg].x), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.position[gg].y, sizeof (*&particleSystem.chargedParticles.position[gg].y), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.position[gg].z, sizeof (*&particleSystem.chargedParticles.position[gg].z), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.velocity[gg].x, sizeof (*&particleSystem.chargedParticles.velocity[gg].x), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.velocity[gg].y, sizeof (*&particleSystem.chargedParticles.velocity[gg].y), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.velocity[gg].z, sizeof (*&particleSystem.chargedParticles.velocity[gg].z), 1, particleSystem.params.saveDataFiles[0]);
                    
                    fwrite(&particleSystem.chargedParticles.Efield1[gg].x, sizeof (*&particleSystem.chargedParticles.Efield1[gg].x), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.Efield1[gg].y, sizeof (*&particleSystem.chargedParticles.Efield1[gg].y), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.Efield1[gg].z, sizeof (*&particleSystem.chargedParticles.Efield1[gg].z), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.Efield2[gg].x, sizeof (*&particleSystem.chargedParticles.Efield2[gg].x), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.Efield2[gg].y, sizeof (*&particleSystem.chargedParticles.Efield2[gg].y), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.Efield2[gg].z, sizeof (*&particleSystem.chargedParticles.Efield2[gg].z), 1, particleSystem.params.saveDataFiles[0]);
                    
                    temp = particleSystem.chargedParticles.velocity[gg].w;
                    fwrite(&temp, sizeof(*&temp), 1, particleSystem.params.saveDataFiles[0]);
//                    fwrite(&particleSystem.chargedParticles.velocity[gg].w, sizeof (*&particleSystem.chargedParticles.velocity[gg].w), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&gg, sizeof (*&gg), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&lastParticle, sizeof (*&lastParticle), 1, particleSystem.params.saveDataFiles[0]);
                }
            }
            break;
        case (SaveDataOptions::newOriginal):
            temp = -1;
            fwrite(&temp, sizeof(*&temp), 1, particleSystem.params.saveDataFiles[0]);
            fwrite(&lastParticle, sizeof(*&lastParticle), 1, particleSystem.params.saveDataFiles[0]);
            fwrite(&currentTime, sizeof(*&currentTime), 1, particleSystem.params.saveDataFiles[0]);
            for (gg = firstParticle; gg < lastParticle; gg++){
                if (particleSystem.chargedParticles.acceleration[gg].w == inboundCheck ||
                        particleSystem.chargedParticles.acceleration[gg].w == affectsAccelerationCheck){
                    temp = particleSystem.chargedParticles.velocity[gg].w;
                    fwrite(&gg, sizeof(*&gg), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&temp, sizeof(*&temp), 1, particleSystem.params.saveDataFiles[0]);

                    fwrite(&particleSystem.chargedParticles.position[gg].x, sizeof(*&particleSystem.chargedParticles.position[gg].x), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.position[gg].y, sizeof(*&particleSystem.chargedParticles.position[gg].y), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.position[gg].z, sizeof(*&particleSystem.chargedParticles.position[gg].z), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.velocity[gg].x, sizeof(*&particleSystem.chargedParticles.velocity[gg].x), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.velocity[gg].y, sizeof(*&particleSystem.chargedParticles.velocity[gg].y), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.velocity[gg].z, sizeof(*&particleSystem.chargedParticles.velocity[gg].z), 1, particleSystem.params.saveDataFiles[0]);
                    
                    fwrite(&particleSystem.chargedParticles.Efield1[gg].x, sizeof(*&particleSystem.chargedParticles.Efield1[gg].x), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.Efield1[gg].y, sizeof(*&particleSystem.chargedParticles.Efield1[gg].y), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.Efield1[gg].z, sizeof(*&particleSystem.chargedParticles.Efield1[gg].z), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.Efield2[gg].x, sizeof(*&particleSystem.chargedParticles.Efield2[gg].x), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.Efield2[gg].y, sizeof(*&particleSystem.chargedParticles.Efield2[gg].y), 1, particleSystem.params.saveDataFiles[0]);
                    fwrite(&particleSystem.chargedParticles.Efield2[gg].z, sizeof(*&particleSystem.chargedParticles.Efield2[gg].z), 1, particleSystem.params.saveDataFiles[0]);
                }
            }
            break;
    }
}

/**
 * Update the acceleration of particles due to particle-particle interactions
 */
void Integrator::updateAcceleration(){
    int j;
    
    iterTimer.startTimer();
    intervalTimer.startTimer();
    fullSimTimer.startTimer();

    // Kick
    particleSystem.chargedParticles.zeroAcceleration();
    forceMethod.updateAcceleration();
    
    const double factor = PhysicalConstant::AMU / PhysicalConstant::ELEM_CHARGE;

    // Compute the electric field due to particle-particle interacitons
    if (particleSystem.chargedParticles.Efield1 != nullptr){
        for (j = firstParticle; j < lastParticle; j++){
            particleSystem.chargedParticles.Efield1[j].x = factor * particleSystem.chargedParticles.acceleration[j].x*particleSystem.chargedParticles.getMass(j)/particleSystem.chargedParticles.position[j].w;
            particleSystem.chargedParticles.Efield1[j].y = factor * particleSystem.chargedParticles.acceleration[j].y*particleSystem.chargedParticles.getMass(j)/particleSystem.chargedParticles.position[j].w;
            particleSystem.chargedParticles.Efield1[j].z = factor * particleSystem.chargedParticles.acceleration[j].z*particleSystem.chargedParticles.getMass(j)/particleSystem.chargedParticles.position[j].w;
        }
    }
    
    iterTimer.pauseTimer();
    intervalTimer.pauseTimer();
    fullSimTimer.pauseTimer();
    
    if (forceMethodComparisons.numElements > 0){
        for (j = firstParticle; j < lastParticle; j++){
            comparisonAcceleration[j] = particleSystem.chargedParticles.acceleration[j];
        }
    }
}

/**
 * Compare the particle-particle acceleration calculated with various different
 *      methods
 * @param iter - simulation iteration
 */
void Integrator::compareAcceleration(int iter){
    int j, k;
    float L2norm, maxL2norm, difference, normalizer;
    int maxL2normIndex, numParticles;
    float tempL2norm;
    int numValidParticles = 0;
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    
    float tempMedian1, tempMedian2, timesLessThan;
    
    float tempTimes[FMMConstant::windowSize];
    
    // Compare L2 norm against the baseline force calculation method
    for (j = 0; j < forceMethodComparisons.numElements; j++){
        L2norm = 0;
        maxL2norm = 0;
        maxL2normIndex = 0;
        if (iter % forceMethodComparisons[j]->compareRate == 0){
            iterTimerComparisons[j].startTimer();
            intervalTimerComparisons[j].startTimer();
            fullSimTimerComparisons[j].startTimer();

            particleSystem.chargedParticles.zeroAcceleration();

            forceMethodComparisons[j]->updateAcceleration();

            iterTimerComparisons[j].pauseTimer();
            intervalTimerComparisons[j].pauseTimer();
            fullSimTimerComparisons[j].pauseTimer();
            
            if (forceMethodComparisons[j]->fullL2NormFile != nullptr){
                numParticles = lastParticle - firstParticle;
                fwrite(&numParticles, sizeof(*&numParticles), 1, forceMethodComparisons[j]->fullL2NormFile);
            }
            
            numValidParticles = 0;

            for (k = firstParticle; k < lastParticle; k++){
                if (particleSystem.chargedParticles.acceleration[k].w == inboundCheck){
                    numValidParticles++;
                    
                    difference = (particleSystem.chargedParticles.acceleration[k].x - comparisonAcceleration[k].x) * 
                            (particleSystem.chargedParticles.acceleration[k].x - comparisonAcceleration[k].x) + 
                            (particleSystem.chargedParticles.acceleration[k].y - comparisonAcceleration[k].y) *
                            (particleSystem.chargedParticles.acceleration[k].y - comparisonAcceleration[k].y) + 
                            (particleSystem.chargedParticles.acceleration[k].z - comparisonAcceleration[k].z) *
                            (particleSystem.chargedParticles.acceleration[k].z - comparisonAcceleration[k].z);

                    normalizer = comparisonAcceleration[k].x * comparisonAcceleration[k].x + 
                            comparisonAcceleration[k].y * comparisonAcceleration[k].y + 
                            comparisonAcceleration[k].z * comparisonAcceleration[k].z;

                    if (normalizer == 0 && difference == 0){
                        tempL2norm = 0;
                        L2norm += 0;
                    }
                    else {
                        tempL2norm = difference / normalizer;
//                        L2norm += difference / normalizer / (lastParticle - firstParticle);
                        L2norm += difference / normalizer;
                    }

                    if (forceMethodComparisons[j]->fullL2NormFile != nullptr){
                        tempL2norm = sqrt(tempL2norm);
                        fwrite(&tempL2norm, sizeof(*&tempL2norm), 1, forceMethodComparisons[j]->fullL2NormFile);
                    }

                    if (difference / normalizer > maxL2norm){
                        maxL2norm = difference / normalizer;
                        maxL2normIndex = k;
                        if (maxL2norm >= 1){
                            printf("test\n");
                        }
                    }
                    if (std::isnan(difference / normalizer) || std::isnan(L2norm)){
                        printf("test\n");
                    }
                }
            }
            maxL2norm = sqrt(maxL2norm);
            L2norm /= numValidParticles;
            L2norm = sqrt(L2norm);
            if (std::isnan(L2norm)){
                printf("test\n");
            }
            printf("Error for %s: %.9g\n", forceMethodComparisons[j]->timerName, L2norm);
            printf("maxL2norm for %s: %.9g\n", forceMethodComparisons[j]->timerName, maxL2norm);
            printf("maxL2normIndex for %s: %d\n", forceMethodComparisons[j]->timerName, maxL2normIndex);
            
            numParticles = lastParticle - firstParticle;

            if (forceMethodComparisons[j]->L2NormFile != nullptr){
                forceMethodComparisons[j]->writeL2NormInfo(numParticles, L2norm, maxL2norm, maxL2normIndex);
            }
        }
    }
    
    
//    if (iter / forceMethod.compareRate >= FMMConstant::windowSize){
//        if (iter / forceMethod.compareRate == FMMConstant::windowSize){
//            for (k = 0; k < FMMConstant::crossoverSize; k++){
//                forceMethod.crossoverPoints[k] = 0;
//            }
//        }
//        if (forceMethod.olderSibling != NULL){
//            for (k = 0; k < FMMConstant::windowSize; k++){
//                tempTimes[k] = forceMethod.previousTotalTimes[k];
//            }
//            tempMedian1 = findTrueMedian(tempTimes, FMMConstant::windowSize);
//            for (k = 0; k < FMMConstant::windowSize; k++){
//                tempTimes[k] = forceMethod.olderSibling->previousTotalTimes[k];
//            }
//            tempMedian2 = findTrueMedian(tempTimes, FMMConstant::windowSize);
//            timesLessThan = 0;
//            for (k = 0; k < FMMConstant::crossoverSize - 1; k++){
//                forceMethod.crossoverPoints[k] = forceMethod.crossoverPoints[k + 1];
//                timesLessThan += forceMethod.crossoverPoints[k];
//            }
//            if (tempMedian2 < tempMedian1){
//                forceMethod.crossoverPoints[FMMConstant::crossoverSize - 1] = 1;
//                timesLessThan++;
////                        forceMethodComparisons[j]->timesLessThan++;
////                        if (forceMethodComparisons[j]->timesLessThan > 50){
////                            forceMethodComparisons[j]->timesLessThan = 0;
////                            forceMethodComparisons[j]->increaseLevel();
////                            forceMethodComparisons[j]->olderSibling->increaseLevel();
////                        }
//            }
//            else{
//                forceMethod.crossoverPoints[FMMConstant::crossoverSize - 1] = 0;
//            }
//            if (timesLessThan > FMMConstant::crossoverSize * 4 / 5){
//                forceMethod.increaseLevel();
//                forceMethod.olderSibling->increaseLevel();
//                for (k = 0; k < FMMConstant::crossoverSize; k++){
//                    forceMethod.crossoverPoints[k] = 0;
//                }
//            }
//        }
//    }
    
    // Check to see the performance between two corresponding methods
    // Mainly used for the multipole method codes to determine the maxLevel,
    // indicating the number of times a base level cell should be subdivided.
    // If one level starts outperforming the other, swap to the better method
    for (j = 0; j < forceMethodComparisons.numElements; j++){
        if (iter / forceMethodComparisons[j]->compareRate >= FMMConstant::windowSize){
            if (iter / forceMethodComparisons[j]->compareRate == FMMConstant::windowSize){
                for (k = 0; k < FMMConstant::crossoverSize; k++){
                    forceMethodComparisons[j]->crossoverPoints[k] = 0;
                }
            }
            if (forceMethodComparisons[j]->olderSibling != nullptr){
                for (k = 0; k < FMMConstant::windowSize; k++){
                    tempTimes[k] = forceMethodComparisons[j]->previousTotalTimes[k];
//                    printf("forceMethodComparisons[j]->previousTotalTimes[%d] = %.9g\n",k,forceMethodComparisons[j]->previousTotalTimes[k]);
                }
                tempMedian1 = findTrueMedian(tempTimes, FMMConstant::windowSize);
                for (k = 0; k < FMMConstant::windowSize; k++){
                    tempTimes[k] = forceMethodComparisons[j]->olderSibling->previousTotalTimes[k];
//                    printf("forceMethodComparisons[j]->olderSibling->previousTotalTimes[%d] = %.9g\n",k,forceMethodComparisons[j]->olderSibling->previousTotalTimes[k]);
//                    printf("tempTimes[%d] = %.9g\n",k,tempTimes[k]);
                }
                tempMedian2 = findTrueMedian(tempTimes, FMMConstant::windowSize);
                timesLessThan = 0;
                for (k = 0; k < FMMConstant::crossoverSize - 1; k++){
                    forceMethodComparisons[j]->crossoverPoints[k] = forceMethodComparisons[j]->crossoverPoints[k + 1];
                    timesLessThan += forceMethodComparisons[j]->crossoverPoints[k];
                }
                if (tempMedian2 < tempMedian1){
                    forceMethodComparisons[j]->crossoverPoints[FMMConstant::crossoverSize - 1] = 1;
                    timesLessThan++;
//                        forceMethodComparisons[j]->timesLessThan++;
//                        if (forceMethodComparisons[j]->timesLessThan > 50){
//                            forceMethodComparisons[j]->timesLessThan = 0;
//                            forceMethodComparisons[j]->increaseLevel();
//                            forceMethodComparisons[j]->olderSibling->increaseLevel();
//                        }
                }
                else{
                    forceMethodComparisons[j]->crossoverPoints[FMMConstant::crossoverSize - 1] = 0;
                }
                if (timesLessThan > FMMConstant::crossoverSize * 4 / 5){
//                    for (k = 0; k < FMMConstant::windowSize; k++){
//                        printf("forceMethodComparisons[j]->previousTotalTimes[%d] = %.9g\n",k,forceMethodComparisons[j]->previousTotalTimes[k]);
//                        printf("forceMethodComparisons[j]->olderSibling->previousTotalTimes[%d] = %.9g\n",k,forceMethodComparisons[j]->olderSibling->previousTotalTimes[k]);
//                    }
//                    printf("MaxLevel increased\n");
                    forceMethodComparisons[j]->increaseLevel();
                    forceMethodComparisons[j]->olderSibling->increaseLevel();
                    for (k = 0; k < FMMConstant::crossoverSize; k++){
                        forceMethodComparisons[j]->crossoverPoints[k] = 0;
                    }
                }
            }
        }
    }
    

    if (forceMethodComparisons.numElements > 0){
        for (j = firstParticle; j < lastParticle; j++){
            particleSystem.chargedParticles.acceleration[j] = comparisonAcceleration[j];
        }
    }
}

/**
 * Finish an iteration of the simulation
 * @param iter - simulation iteration
 */
void Integrator::finishIter(int iter){
    int j;
    if (particleSystem.params.timerInfo == TimerOptions::verbose){
        iterTimer.printTimeRecorded();
        for (j = 0; j < forceMethodComparisons.numElements; j++){
            iterTimerComparisons[j].printTimeRecorded();
        }
    }
    iterTimer.resetTimer();
    for (j = 0; j < forceMethodComparisons.numElements; j++){
        iterTimerComparisons[j].resetTimer();
    }
    
    
    if (iter % particleSystem.params.interval == 0){
        updateEnergyIntervalData();

        intervalCount++;
        intervalTimer.changeCount(intervalCount);
        for (j = 0; j < forceMethodComparisons.numElements; j++){
            intervalTimerComparisons[j].changeCount(intervalCount);
        }
        finishInterval();
    }
    
    int saveDataToggle = false;
    switch (saveDataRate){
        case (SaveDataRateOptions::constantRate):
            if (iter % particleSystem.params.interval == 0) saveDataToggle = true;
            break;
        case (SaveDataRateOptions::lastSnapshots):
            if (particleSystem.params.numIters - particleSystem.params.numItersStored < iter){
                saveDataToggle = true;
            }
            break;
        case (SaveDataRateOptions::constantRate_and_lastSnapshots):
            if (iter % particleSystem.params.interval == 0 || 
                    particleSystem.params.numIters - particleSystem.params.numItersStored < iter){
                saveDataToggle = true;
            }
            break;
    }
    
    if (saveDataToggle){
        saveData();
    }
}

/**
 * Finish an interval of the simulation
 */
void Integrator::finishInterval(){
    int j;
    
    if (particleSystem.params.timerInfo == TimerOptions::verbose ||
            particleSystem.params.timerInfo == TimerOptions::minimized){
        intervalTimer.printTimeRecorded();
        for (j = 0; j < forceMethodComparisons.numElements; j++){
            intervalTimerComparisons[j].printTimeRecorded();
        }
    }
    intervalTimer.resetTimer();
    for (j = 0; j < forceMethodComparisons.numElements; j++){
        intervalTimerComparisons[j].resetTimer();
    }
    printf("\n");
//    saveData();
    prevTime = currentTime;
}

/**
 * Finish the simulation and print the total time recorded
 */
void Integrator::finishFullSim(){
    int j;
    
    if (particleSystem.params.timerInfo != TimerOptions::off){
        fullSimTimer.printTimeRecorded();
        for (j = 0; j < forceMethodComparisons.numElements; j++){
            fullSimTimerComparisons[j].printTimeRecorded();
        }
    }
    
    if (forceMethodComparisons.numElements > 0){
        delete[] comparisonAcceleration;
        comparisonAcceleration = nullptr;
    }
}

/**
 * Drift step in leapfrogDKD
 *      Only drift particles which have been injected at the beginning of the time step
 * @param tempLastParticlePrev - index of the last particle in the previous iteration
 * @param tempLastParticle - index of the current last particle
 * @param dt - time step (half of the original simulation time step)
 */
void LeapfrogDKD::drift(int& tempLastParticlePrev, int& tempLastParticle, float dt){
    // Only drift particles which have been injected at the beginning of the time step
    int particleDelta;
    tempLastParticle = particleSystem.chargedParticles.findMostRecentInjected(tempPosTime);
    setLastChargedParticle(tempLastParticle + 1);
    particleDelta = tempLastParticle + 1 - tempLastParticlePrev;
    
    updateEnergyInjectionData(tempLastParticlePrev, particleDelta);

    tempLastParticlePrev = tempLastParticle + 1;

    // Drift
    particleSystem.propagatePositionNew(tempPosTime, dt);
    resetFirstChargedParticle();

    tempPosTime += dt;
}

/**
 * Drift step in leapfrogDKD
 *      Drift particles by some amount if they will be injected at some time during the time step
 * @param tempLastParticlePrev - index of the last particle in the previous iteration
 * @param tempLastParticle - index of the current last particle
 * @param dt - time step (half of the original simulation time step)
 */
void LeapfrogDKD::driftNew(int& tempLastParticlePrev, int& tempLastParticle, float dt) {
    // Drift particles by some amount if they will be injected at some time during the time step
    int particleDelta;
    tempLastParticle = particleSystem.chargedParticles.findMostRecentInjected(tempPosTime + dt);
    setLastChargedParticle(tempLastParticle + 1);
    particleDelta = tempLastParticle + 1 - tempLastParticlePrev;
    
    updateEnergyInjectionData(tempLastParticlePrev, particleDelta);

    tempLastParticlePrev = tempLastParticle + 1;

    // Drift
    particleSystem.propagatePositionNew(tempPosTime, dt);
    resetFirstChargedParticle();

    tempPosTime += dt;
}

/**
 * Kick step in leapfrogDKD
 * @param iter - simulation iteration
 * @param dt - time step (equal to the original simulation time step)
 */
void LeapfrogDKD::kick(int iter, float dt){
    updateAcceleration();
    compareAcceleration(iter);
    
//    for (int i = firstParticle; i < lastParticle; i++){
//        particleSystem.chargedParticles.printAccelInfo(i);
//    }

    particleSystem.computeBackgroundAcceleration();
    
//    for (int i = firstParticle; i < lastParticle; i++){
//        particleSystem.chargedParticles.printAccelInfo(i);
//    }
//    printf("\n");
    
//    if (lastParticle > 95){
//        particleSystem.chargedParticles.printAccelInfo(0);
//    }

    // Kick
    particleSystem.propagateVelocityNew(tempVelTime, dt);

    tempVelTime += dt;
}

/**
 * Compute the initial kinetic and field energy that each particle will have
 *      upon injection
 */
void Integrator::initializeEnergyData(){
    int i;
    if (energyToggle == 1){
        initialPositions = new vec3<float>[particleSystem.chargedParticles.maxNumParticles];
        particleSystem.chargedParticles.allocateMemory_Energy();
        
        for (i = 0; i < particleSystem.chargedParticles.maxNumParticles; i++){
            particleSystem.computeChargedInjectionFieldEnergy(i, numEnergySteps_initialization);
            particleSystem.computeChargedInjectionKineticEnergy(i);
            
            if (i%50 == 0)
                printf("%d: %.9g\n",i, particleSystem.chargedParticles.fieldEnergy[i] / PhysicalConstant::ELEM_CHARGE);
            
            initialPositions[i].x = particleSystem.chargedParticles.position[i].x;
            initialPositions[i].y = particleSystem.chargedParticles.position[i].y;
            initialPositions[i].z = particleSystem.chargedParticles.position[i].z;
        }
    }
}

/**
 * Update the total injection energy by adding the energy of newly injected particles
 * @param tempLastParticlePrev - index of the last particle in the previous iteration
 * @param particleDelta - difference between tempLastParticlePrev and the index of
 *      the current last particle
 */
void Integrator::updateEnergyInjectionData(int tempLastParticlePrev, int particleDelta){
    int i;
    if (energyToggle == 1){
        for (i = 0; i < particleDelta; i++){
            particleSystem.computeChargedInjectionElectricPotentialEnergy(tempLastParticlePrev + i);
            KE_inj += particleSystem.chargedParticles.injectionKE[tempLastParticlePrev + i];
            PE_inj += particleSystem.chargedParticles.injectionPE[tempLastParticlePrev + i];
            FE_inj += particleSystem.chargedParticles.fieldEnergy[tempLastParticlePrev + i];
        }
    }
}

/**
 * Update the current energy of all particles in the simulation, compare it to
 *      the total injected energy, and write the data to a file
 */
void Integrator::updateEnergyIntervalData(){
    int i;
    float tempTime;
    float timeStep = particleSystem.params.timeStep;
    int FE_steps;
    if (energyToggle == 1){
        KE = 0; PE = 0; FE = 0;
        for (i = firstParticle; i < lastParticle; i++){
            tempTime = std::min(timeStep * particleSystem.params.interval,
                    tempPosTime - particleSystem.chargedParticles.timeInjected[i]);
            FE_steps = std::max(tempTime * numEnergySteps / timeStep, 2.0f);
            particleSystem.updateChargedFieldEnergy(i, FE_steps, initialPositions[i]);
            
            initialPositions[i].x = particleSystem.chargedParticles.position[i].x;
            initialPositions[i].y = particleSystem.chargedParticles.position[i].y;
            initialPositions[i].z = particleSystem.chargedParticles.position[i].z;
            
            KE += particleSystem.computeChargedKineticEnergy(i);
            PE += particleSystem.computeChargedElectricPotentialEnergy(i);
            FE += particleSystem.chargedParticles.fieldEnergy[i];
        }
        PE /= 2; // Halve potential energy due to counting each interaction twice
        
        printf("Interval %d: KE = %.9g, PE = %.9g, FE = %.9g\n",intervalCount,KE / PhysicalConstant::ELEM_CHARGE, PE / PhysicalConstant::ELEM_CHARGE, FE / PhysicalConstant::ELEM_CHARGE);
        printf("Injected Energy: KE_inj = %.9g, PE_inj = %.9g, FE_inj = %.9g\n",KE_inj / PhysicalConstant::ELEM_CHARGE, PE_inj / PhysicalConstant::ELEM_CHARGE,FE_inj / PhysicalConstant::ELEM_CHARGE);
        printf("Total Energy: E = %.9g\n",(KE + PE + FE) / PhysicalConstant::ELEM_CHARGE);
        printf("Total Energy Injected: E_inj = %.9g\n",(KE_inj + PE_inj + FE_inj) / PhysicalConstant::ELEM_CHARGE);
        printf("Energy Ratio: E / E_inj = %.9g\n\n",(KE + PE + FE) / (KE_inj + PE_inj + FE_inj));

        fwrite(&KE_inj, sizeof(*&KE_inj), 1, energyFile);
        fwrite(&PE_inj, sizeof(*&PE_inj), 1, energyFile);
        fwrite(&FE_inj, sizeof(*&FE_inj), 1, energyFile);
        fwrite(&KE, sizeof(*&KE), 1, energyFile);
        fwrite(&PE, sizeof(*&PE), 1, energyFile);
        fwrite(&FE, sizeof(*&FE), 1, energyFile);
    }
}

/**
 * Initialize the simulation by zeroing the acceleration, and save the initial
 *      positions to the state file
 */
void Integrator::initializeSimulation(){
    setFirstChargedParticle(0);
    setLastChargedParticle(particleSystem.chargedParticles.maxNumParticles);
    particleSystem.chargedParticles.zeroAcceleration();
    
    saveData();
    
    initializeEnergyData();
    
    setLastChargedParticle(0);
}

/**
 * Run the leapfrog drift-kick-drift integration scheme on a system of particles
 */
void LeapfrogDKD::runSimulation() {
    int i;
    intervalCount = 0;
    prevTime = 0.0f;
    int tempLastParticle, tempLastParticlePrev;
    tempLastParticlePrev = 0;
    
    float timeStep = particleSystem.params.timeStep;
    
    initializeSimulation();
    
    for (i = 1; i <= particleSystem.params.numIters; i++){
        iterTimer.changeCount(i);
        
        tempPosTime = timeStep * (i - 1);
        tempVelTime = timeStep * (i - 1);
        
        
        driftNew(tempLastParticlePrev, tempLastParticle, timeStep / 2);
        
        kick(i, timeStep);
        
        driftNew(tempLastParticlePrev, tempLastParticle, timeStep / 2);
        
//        for (int j = 0; j < std::min(lastParticle, 10); j++){
//            particleSystem.chargedParticles.printPosInfo(j);
//            particleSystem.chargedParticles.printVelInfo(j);
//        }
//        printf("\n");
        
        currentTime = timeStep * i;
        finishIter(i);
    }
    
    finishFullSim();
}



/**
 * Update the acceleration for the particles in one of the substeps due to
 *      particle-particle interactions
 */
void LeapfrogDKDSubsteps::updateAccelerationSubsteps(){
    int j;
    
    iterTimer.startTimer();
    intervalTimer.startTimer();
    fullSimTimer.startTimer();
    
    // Kick
    particleSystem.chargedParticles.zeroAcceleration();
    substepMethod.updateAcceleration();
    
    const double factor = PhysicalConstant::AMU / PhysicalConstant::ELEM_CHARGE;

    if (particleSystem.chargedParticles.Efield1 != nullptr){
        for (j = firstParticle; j < lastParticle; j++){
            particleSystem.chargedParticles.Efield1[j].x = factor * particleSystem.chargedParticles.acceleration[j].x*particleSystem.chargedParticles.getMass(j)/particleSystem.chargedParticles.position[j].w;
            particleSystem.chargedParticles.Efield1[j].y = factor * particleSystem.chargedParticles.acceleration[j].y*particleSystem.chargedParticles.getMass(j)/particleSystem.chargedParticles.position[j].w;
            particleSystem.chargedParticles.Efield1[j].z = factor * particleSystem.chargedParticles.acceleration[j].z*particleSystem.chargedParticles.getMass(j)/particleSystem.chargedParticles.position[j].w;
        }
    }
    
    iterTimer.pauseTimer();
    intervalTimer.pauseTimer();
    fullSimTimer.pauseTimer();
}

/**
 * Check the boundary conditions for a specified particle to determine which 
 *      substep time the particle should be moved to
 * @param index - index of particle being considered
 * @param level - substep level to check
 */
void LeapfrogDKDSubsteps::checkConditions(int index, int level){
    switch (substepBoundaryConditions[level + 1].bc.boundaryConditionType){
        case (BoundaryConditionOptions::none):
            break;
        case (BoundaryConditionOptions::space):
            checkSpatialConditions(index, level);
            break;
        case (BoundaryConditionOptions::time):
            checkTemporalConditions(index, level);
            break;
        case (BoundaryConditionOptions::space_and_time):
            checkSpatialConditions(index, level);
            int temp = particleLevels[index];
            checkTemporalConditions(index, level);
            particleLevels[index] = std::max(particleLevels[index], temp);
            break;
    }
}

/**
 * Check the spatial boundary conditions for a specified particle to determine 
 *      which substep time the particle should be moved to
 * @param index - index of particle being considered
 * @param level - substep level to check
 */
void LeapfrogDKDSubsteps::checkSpatialConditions(int index, int level) {
    vec4<float> tempBC;
    vec4<float> tempPos = particleSystem.chargedParticles.position[index];
    tempPos.w = sqrt(tempPos.x * tempPos.x + tempPos.y * tempPos.y);
    
    int i = maxLevel;
    int levelFound = false;
    int possibleLevel = -1;
    float tempBoundary;
    while (i >= level && levelFound == false){
        tempBC = substepBoundaryConditions[i].bc.spatialBoundaryConditions;
        
        if (possibleLevel == -1){
            if (std::abs(tempPos.x) < tempBC.x || tempBC.x == 0){
                possibleLevel = i;
                tempBoundary = tempBC.x;
            }
        }
        else {
            if (tempBoundary == tempBC.x){
                possibleLevel--;
            }
            else {
                levelFound = true;
            }
        }
        i--;
    }
        
    particleLevels[index] = possibleLevel;
    
    i = maxLevel;
    levelFound = false;
    possibleLevel = -1;
    while (i >= level && levelFound == false){
        tempBC = substepBoundaryConditions[i].bc.spatialBoundaryConditions;
        
        if (possibleLevel == -1){
            if (std::abs(tempPos.y) < tempBC.y || tempBC.y == 0){
                possibleLevel = i;
                tempBoundary = tempBC.y;
            }
        }
        else {
            if (tempBoundary == tempBC.y){
                possibleLevel--;
            }
            else {
                levelFound = true;
            }
        }
        i--;
    }
        
    particleLevels[index] = std::max(particleLevels[index], possibleLevel);
    
    i = maxLevel;
    levelFound = false;
    possibleLevel = -1;
    while (i >= level && levelFound == false){
        tempBC = substepBoundaryConditions[i].bc.spatialBoundaryConditions;
        
        if (possibleLevel == -1){
            if (std::abs(tempPos.z) < tempBC.z || tempBC.z == 0){
                possibleLevel = i;
                tempBoundary = tempBC.z;
            }
        }
        else {
            if (tempBoundary == tempBC.z){
                possibleLevel--;
            }
            else {
                levelFound = true;
            }
        }
        i--;
    }
//    if (level == 0){
//        printf("test\n");
//    }
//    if (possibleLevel == 0){
//        printf("test\n");
//    }
        
    particleLevels[index] = std::max(particleLevels[index], possibleLevel);
    
    i = maxLevel;
    levelFound = false;
    possibleLevel = -1;
    while (i >= level && levelFound == false){
        tempBC = substepBoundaryConditions[i].bc.spatialBoundaryConditions;
        
        if (possibleLevel == -1){
            if (std::abs(tempPos.w) < tempBC.w || tempBC.w == 0){
                possibleLevel = i;
                tempBoundary = tempBC.w;
            }
        }
        else {
            if (tempBoundary == tempBC.w){
                possibleLevel--;
            }
            else {
                levelFound = true;
            }
        }
        i--;
    }
        
    particleLevels[index] = std::max(particleLevels[index], possibleLevel);
}

/**
 * Check the temporal boundary conditions for a specified particle to determine 
 *      which substep time the particle should be moved to
 * @param index - index of particle being considered
 * @param level - substep level to check
 */
void LeapfrogDKDSubsteps::checkTemporalConditions(int index, int level) {
    
    int i = maxLevel;
    int levelFound = false;
    int possibleLevel = -1;
    float tempBC;
    float tempInjectTime = particleSystem.chargedParticles.timeInjected[index];
    float tempBoundary;
    while (i >= level && levelFound == false){
        tempBC = substepBoundaryConditions[i].bc.temporalBoundaryCondition;
        
        if (possibleLevel == -1){
            if ((tempPosTime - tempInjectTime < tempBC) || tempBC == 0){
                possibleLevel = i;
                tempBoundary = tempBC;
            }
        }
        else {
            if (tempBoundary == tempBC){
                possibleLevel--;
            }
            else {
                levelFound = true;
            }
        }
        i--;
    }
        
    particleLevels[index] = possibleLevel;
}

/**
 * Change the substep level in which a particle exists
 * @param level - substep level
 */
void LeapfrogDKDSubsteps::crossoverParticles(int level){
    int i;
    for (i = firstParticle; i < lastParticle; i++){
        
        if (particleLevels[i] == level || particleLevels[i] == level + 1){
            checkConditions(i, level);
        }
    }
}

/**
 * Drift step in leapfrogDKDSubsteps
 *      Only drift particles which have been injected at the beginning of the time step
 * @param iter - iteration
 * @param tempLastParticlePrev - index of the last particle in the previous iteration
 * @param tempLastParticle - index of the current last particle
 * @param level - substep level
 * @param toggle - toggle for if second drift step in (drift-kick-drift) sequence
 */
void LeapfrogDKDSubsteps::drift(int iter, int& tempLastParticlePrev, int& tempLastParticle, int level, int toggle){
    if (level != maxLevel) {
        drift(iter, tempLastParticlePrev, tempLastParticle, level + 1, false);
        kick(iter,  tempLastParticlePrev, tempLastParticle, level + 1, level + 1);
        drift(iter, tempLastParticlePrev, tempLastParticle, level + 1, true);
        if (toggle){
            crossoverParticles(level);
        }
    }
    else {
        // Only drift particles which have been injected at the beginning of the timestep
        int particleDelta;
        float dt = substepBoundaryConditions[maxLevel].dt / 2;
        tempLastParticle = particleSystem.chargedParticles.findMostRecentInjected(tempPosTime);
        setLastChargedParticle(tempLastParticle + 1);
        particleDelta = tempLastParticle + 1 - tempLastParticlePrev;

        updateEnergyInjectionData(tempLastParticlePrev, particleDelta);

        tempLastParticlePrev = tempLastParticle + 1;

        // Drift
        particleSystem.propagatePositionNew(tempPosTime, dt);
        resetFirstChargedParticle();

        tempPosTime += dt;
    }
}

/**
 * Drift step in leapfrogDKDSubsteps
 *      Drift particles by some amount if they will be injected at some time during the time step
 * @param iter - iteration
 * @param tempLastParticlePrev - index of the last particle in the previous iteration
 * @param tempLastParticle - index of the current last particle
 * @param level - substep level
 * @param toggle - toggle for if second drift step in (drift-kick-drift) sequence
 */
void LeapfrogDKDSubsteps::driftNew(int iter, int& tempLastParticlePrev, int& tempLastParticle, int level, int toggle){
    if (level != maxLevel) {
        driftNew(iter, tempLastParticlePrev, tempLastParticle, level + 1, false);
        kick(iter, tempLastParticlePrev, tempLastParticle, level + 1, level + 1);
        driftNew(iter, tempLastParticlePrev, tempLastParticle, level + 1, true);
//        if (level == 0){
//            printf("test\n");
//        }
        if (toggle){
            crossoverParticles(level);
        }
        
//        for (int i = firstParticle; i < lastParticle; i++){
//            particleSystem.chargedParticles.printPosInfo(i);
//            particleSystem.chargedParticles.printVelInfo(i);
//        }
//        printf("\n");
    }
    else {
        // Drift particles by some amount if they will be injected at some time during the timestep
        int particleDelta;
        float dt = substepBoundaryConditions[maxLevel].dt / 2;
        tempLastParticle = particleSystem.chargedParticles.findMostRecentInjected(tempPosTime + dt);
        setLastChargedParticle(tempLastParticle + 1);
        particleDelta = tempLastParticle + 1 - tempLastParticlePrev;

        updateEnergyInjectionData(tempLastParticlePrev, particleDelta);

        tempLastParticlePrev = tempLastParticle + 1;

        // Drift
        particleSystem.propagatePositionNew(tempPosTime, dt);
        resetFirstChargedParticle();

        tempPosTime += dt;
    }
}
//
//void LeapfrogDKD::kick(int iter, float dt){
//    updateAcceleration();
//    compareAcceleration(iter);
//
//    particleSystem.computeBackgroundAcceleration();
//
//    // Kick
//    particleSystem.propagateVelocityNew(tempVelTime, dt);
//
//    tempVelTime += dt;
//}

/**
 * Set the status of the particle depending on the substep level
 * @param level - substep level
 * @param topLevel - top level that a particle can still be inbounds
 */
void LeapfrogDKDSubsteps::setParticleStatus(int level, int topLevel) {
    int i;
    const float outOfBoundsStatus = static_cast<float>(ParticleStatus::outOfBounds);
    const float crashedStatus = static_cast<float>(ParticleStatus::crashed);
    const float inboundsStatus = static_cast<float>(ParticleStatus::inbound);
    const float affectsAccelerationStatus = static_cast<float>(ParticleStatus::affectsAcceleration);
    
    for (i = firstParticle; i < lastParticle; i++){
        if (particleSystem.chargedParticles.acceleration[i].w != outOfBoundsStatus && 
                particleSystem.chargedParticles.acceleration[i].w != crashedStatus){
            if (particleLevels[i] <= level && particleLevels[i] >= topLevel){
                particleSystem.chargedParticles.acceleration[i].w = inboundsStatus;
            }
            else if (particleLevels[i] < topLevel){
                particleSystem.chargedParticles.acceleration[i].w = affectsAccelerationStatus;
            }
        }
    }
}

/**
 * Kick each particle based on the substep level that the particle exists in
 */
void LeapfrogDKDSubsteps::variableTimestepKick(){
    int i;
    float dt_shift, temp_dt;
    const float inboundCheck = static_cast<float>(ParticleStatus::inbound);
    
    float* tend = new float[maxLevel + 1];
    for (i = 0; i <= maxLevel; i++){
        tend[i] = tempVelTime + substepBoundaryConditions[i].dt;
    }
    
    for (i = firstParticle; i < lastParticle; i++){        
        if (particleSystem.chargedParticles.acceleration[i].w == inboundCheck){
            temp_dt = substepBoundaryConditions[particleLevels[i]].dt;
            dt_shift = std::min(std::max(tend[particleLevels[i]] - particleSystem.chargedParticles.timeInjected[i], 0.0f), temp_dt);
            
//            if (lastParticle > 95 && i == 0){
//                printf("dt = %.9g\n", dt_shift);
//            }
            
//            if (lastParticle > 95 && i == 0){
//                particleSystem.chargedParticles.printPosInfo(0);
//                particleSystem.chargedParticles.printVelInfo(0);
//                particleSystem.chargedParticles.printAccelInfo(0);
//            }
            
            particleSystem.chargedParticles.velocity[i].x += particleSystem.chargedParticles.acceleration[i].x * dt_shift;
            particleSystem.chargedParticles.velocity[i].y += particleSystem.chargedParticles.acceleration[i].y * dt_shift;
            particleSystem.chargedParticles.velocity[i].z += particleSystem.chargedParticles.acceleration[i].z * dt_shift;
        }
//        else{
//            printf("test\n");
//        }
    }
    
    delete[] tend;
}

/**
 * Kick step in leapfrogDKDSubsteps
 * @param iter - simulation iteration
 * @param tempLastParticlePrev - index of the last particle in the previous iteration
 * @param tempLastParticle - index of the current last particle
 * @param level - substep level
 * @param topLevel - top level that a particle can still be inbounds
 */
void LeapfrogDKDSubsteps::kick(int iter, int& tempLastParticlePrev, int& tempLastParticle, int level, int topLevel){
//    setParticleStatus(level, topLevel); 
//    // Set all particles in levels between level and toplevel to inbound (if not out-of-bounds or crashed)
//    // Set all other particles in levels above toplevel to affects acceleration
//    
//    if (level == 0){
//        updateAcceleration();
//        compareAcceleration(iter);
//    }
//    else {
//        substepMethod.updateAcceleration();
//    }
//    
//    particleSystem.computeBackgroundAcceleration();
//    // Only updating acceleration for particles that are inbounds, which agrees with the result of setParticleStatus
//    
//    // Need a way to kick particles over different timesteps based on which level they are in
//    variableTimestepKick();
    
    
    if (level != maxLevel) {
        driftNew(iter, tempLastParticlePrev, tempLastParticle, level + 1, false);
        kick(iter, tempLastParticlePrev, tempLastParticle, level + 1, topLevel);
        driftNew(iter, tempLastParticlePrev, tempLastParticle, level + 1, true);
        
//        for (int i = firstParticle; i < lastParticle; i++){
//            particleSystem.chargedParticles.printPosInfo(i);
//            particleSystem.chargedParticles.printVelInfo(i);
//        }
//        printf("\n");
    }
    else {
        setParticleStatus(level, topLevel); 
        // Set all particles in levels between level and toplevel to inbound (if not out-of-bounds or crashed)
        // Set all other particles in levels above toplevel to affects acceleration

        if (topLevel == 0){
            updateAcceleration();
            compareAcceleration(iter);
        }
        else {
            updateAccelerationSubsteps();
        }
        
//        for (int i = firstParticle; i < lastParticle; i++){
//            particleSystem.chargedParticles.printAccelInfo(i);
//        }

        particleSystem.computeBackgroundAcceleration();
        // Only updating acceleration for particles that are inbounds, which agrees with the result of setParticleStatus
//        for (int i = firstParticle; i < lastParticle; i++){
//            particleSystem.chargedParticles.printAccelInfo(i);
//        }
//        printf("\n");
        // Need a way to kick particles over different timesteps based on which level they are in
        variableTimestepKick();

        tempVelTime += substepBoundaryConditions[maxLevel].dt;
    }
}

/**
 * Deallocate the boundary conditions
 */
void LeapfrogDKDSubsteps::freeBoundaryConditions(){
    if (substepBoundaryConditions != nullptr) delete[] substepBoundaryConditions;
    substepBoundaryConditions = nullptr;
}

/**
 * Set the substep times and boundary conditions under which a particle will use
 *      a certain substep time
 * @param substepBC - array of boundary conditions for substeps
 * @param size - size of substepBC array
 */
void LeapfrogDKDSubsteps::setSubstepTimes(SubstepBoundaryConditions* substepBC, int size){
    int i;
    std::sort(substepBC, substepBC + size);
    int trueSize = 0;
    
    // Find number of 3 partitions required to satisfy smallest timestep size
    float min_dt = substepBC[0].dt;
    float temp_dt = particleSystem.params.timeStep;
    while (temp_dt > min_dt) {
        temp_dt /= 3.0;
        trueSize++;
    }
    
    maxLevel = trueSize;
    temp_dt = particleSystem.params.timeStep;
    substepBoundaryConditions = new SubstepBoundaryConditions[trueSize + 1];
    
    if (trueSize != 0){
        // If smallest substep is larger than timestep, default to using timestep
        substepBoundaryConditions[trueSize] = substepBC[0];
    }
    
    for (i = 0; i <= trueSize; i++){
        substepBoundaryConditions[i].bc.setMostStringentBoundaryConditionType();
    }
    
    int tempInd;
    
    for (i = 1; i < size; i++){
        tempInd = 0;
        temp_dt = particleSystem.params.timeStep;
        while (temp_dt > substepBC[i].dt){
            temp_dt /= 3.0;
            tempInd++;
        }
        
        if (substepBoundaryConditions[tempInd].dt == 0){
            // Set first conditions for specified level
            substepBoundaryConditions[tempInd].dt = temp_dt;
            if (tempInd != 0){
                substepBoundaryConditions[tempInd].bc = substepBC[i].bc;
            }
        }
        else {
//            // Choose the maximum set of boundary conditions for specified level
//            substepBoundaryConditions[tempInd].bc.chooseMaximumBoundaryConditionParameters(substepBC[i].bc);
            
            // Choose the minimum set of boundary conditions for specified level
            substepBoundaryConditions[tempInd].bc.chooseMinimumBoundaryConditionParameters(substepBC[i].bc);
        }
    }
    
//    for (i = trueSize - 1; i > 0; i--){
//        if (substepBoundaryConditions[i].dt == 0){
//            substepBoundaryConditions[i].bc = substepBoundaryConditions[i + 1].bc;
//        }
//    }
    
//    for (i = trueSize - 1; i > 0; i--){
//        if (substepBoundaryConditions[i].dt == 0){
//            substepBoundaryConditions[i].bc = substepBoundaryConditions[i + 1].bc;
//        }
//    }
    
    substepBoundaryConditions[0].dt = particleSystem.params.timeStep;
    for (i = 1; i <= trueSize; i++){
        substepBoundaryConditions[i].dt = substepBoundaryConditions[i - 1].dt / 3.0;
    }
    
    for (i = maxLevel - 1; i >= 0; i--){
        substepBoundaryConditions[i].bc.chooseMaximumBoundaryConditionParameters(substepBoundaryConditions[i + 1].bc);
    }
}

/**
 * Set the substep times and boundary conditions under which a particle will use
 *      a certain substep time
 * @param substepTimes - array of times for each substep
 * @param bc - array of boundary conditions for each substep
 * @param size - size of substepTimes and bc arrays
 */
void LeapfrogDKDSubsteps::setSubstepTimes(float* substepTimes, BoundaryConditionParameters* bc, int size){
    int i;
    freeBoundaryConditions();
    
    SubstepBoundaryConditions* substepBC = new SubstepBoundaryConditions[size];
    
    for (i = 0; i < size; i++){
        substepBC[i].dt = substepTimes[i];
        substepBC[i].bc = bc[i];
    }
    
    setSubstepTimes(substepBC, size);
    
    delete[] substepBC;
    
    initializeParticleLevels();
}

/**
 * Set the default substep times to have no substeps (same as LeapfrogDKD class)
 */
void LeapfrogDKDSubsteps::setDefaultSubstepTimes(){
    freeBoundaryConditions();
    
    maxLevel = 0;
    substepBoundaryConditions = new SubstepBoundaryConditions[1];
    substepBoundaryConditions[0].dt = particleSystem.params.timeStep;
    initializeParticleLevels();
}

/**
 * Initialize an array to hold the level each particle lies within
 */
void LeapfrogDKDSubsteps::initializeParticleLevels(){
    if (particleLevels != nullptr) delete[] particleLevels;
    particleLevels = new int[particleSystem.chargedParticles.maxNumParticles];
    
    int i;
    for (i = 0; i < particleSystem.chargedParticles.maxNumParticles; i++){
        particleLevels[i] = maxLevel;
    }
}

/**
 * Run the leapfrog drift-kick-drift integration scheme on a system of particles
 */
void LeapfrogDKDSubsteps::runSimulation() {
    int i;
    intervalCount = 0;
    prevTime = 0.0f;
    int tempLastParticle, tempLastParticlePrev;
    tempLastParticlePrev = 0;
    
    float timeStep = particleSystem.params.timeStep;
    
    if (substepBoundaryConditions == nullptr) setDefaultSubstepTimes();
    
    initializeSimulation();
    
    for (i = 1; i <= particleSystem.params.numIters; i++){
        iterTimer.changeCount(i);
        
        tempPosTime = timeStep * (i - 1);
        tempVelTime = timeStep * (i - 1);
        
        
        driftNew(i, tempLastParticlePrev, tempLastParticle, 0, false);
        
        kick(i, tempLastParticlePrev, tempLastParticle, 0, 0);
        
        driftNew(i, tempLastParticlePrev, tempLastParticle, 0, true);
        
//        for (int j = 0; j < std::min(lastParticle, 10); j++){
//            particleSystem.chargedParticles.printPosInfo(j);
//            particleSystem.chargedParticles.printVelInfo(j);
//        }
//        printf("\n");
        
        currentTime = timeStep * i;
        finishIter(i);
    }
    
    finishFullSim();
}