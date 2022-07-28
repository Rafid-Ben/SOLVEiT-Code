/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "particles.h"

//unsigned int randomizerSeed = std::random_device()();
unsigned int randomizerSeed = 1;
//Randomizer randomizer(randomizerSeed);  // Use for repeatable random values
////Randomizer randomizer;    // Use for random numbers
RandomizerNew randomizerNew(randomizerSeed);  // Use for repeatable random values
//RandomizerNew randomizerNew;    // Use for random numbers

//void ParticleSystem::propagatePosition(float dt){
//    chargedParticles.propagatePosition(dt);
//    neutralParticles.propagatePosition(dt);
//    
//    enforceSymmetryConditions();
//    enforceBoundaryConditions();
//}

//void ParticleSystem::propagateVelocity(float dt){
//    chargedParticles.propagateVelocity(dt);
//    neutralParticles.propagateVelocity(dt);
//}

/**
 * Propagate the position of all particles in the system using their velocities
 * @param tstart - sim time before propagation
 * @param dt - time step
 */
void ParticleSystem::propagatePositionNew(float tstart, float dt){
    chargedParticles.propagatePositionNew(tstart, dt, params.crashParameters, params.eventsFile);
    neutralParticles.propagatePositionNew(tstart, dt, params.crashParameters, params.eventsFile);
    
    enforceSymmetryConditions();
    enforceBoundaryConditions(tstart + dt);
}

/**
 * Propagate the velocity of all particles in the system using their accelerations
 * @param tstart - sim time before propagation
 * @param dt - time step
 */
void ParticleSystem::propagateVelocityNew(float tstart, float dt){
    fragment(); // Not sure where to put fragmentation? - before or after propagating velocities or elsewhere
    chargedParticles.propagateVelocityNew(tstart, dt);
    neutralParticles.propagateVelocityNew(tstart, dt);
}

/**
 * Write to an events file if a particle crashes, goes out of bounds, etc.
 * @param index - index of particle that had an event
 * @param time - time of event
 * @param particleType - type of particle (charged or neutral)
 */
void ParticleSystem::writeToEventsFile(int index, float time, ParticleType particleType){
    switch (particleType){
        case (ParticleType::charged):
            chargedParticles.writeToEventsFile(index, time, params.eventsFile);
            break;
        case (ParticleType::neutral):
            break;
    }
}

/**
 * Enforces the boundary conditions for the system of particles
 * @param currentTime - current sim time
 */
void ParticleSystem::enforceBoundaryConditions(float currentTime){
    switch (params.boundaryConditions.boundaryConditionType){
        case (BoundaryConditionOptions::none):
            break;
        case (BoundaryConditionOptions::space):
            enforceSpatialBoundaryConditions(currentTime);
            break;
        case (BoundaryConditionOptions::time):
            enforceTemporalBoundaryConditions(currentTime);
            break;
        case (BoundaryConditionOptions::space_and_time):
            enforceSpatialBoundaryConditions(currentTime);
            enforceTemporalBoundaryConditions(currentTime);
            break;
    }
}

/**
 * Enforces spatial boundary conditions for the system of particles
 * @param currentTime - current sim time
 */
void ParticleSystem::enforceSpatialBoundaryConditions(float currentTime){
    int i;
    const float inboundCheck = static_cast<float>(ParticleStatus::inbound);
    int outOfBounds;
    if (params.boundaryConditions.spatialBoundaryConditions.x == 0 &&
            params.boundaryConditions.spatialBoundaryConditions.y == 0 &&
            params.boundaryConditions.spatialBoundaryConditions.z == 0 &&
            params.boundaryConditions.spatialBoundaryConditions.w == 0){
        printf("Warning: Spatial boundary conditions not initialized!\n");
    }
    for (i = chargedParticles.firstParticle; i < chargedParticles.lastParticle; i++){
        if (chargedParticles.acceleration[i].w == inboundCheck){
            outOfBounds = checkOutOfSpatialBounds(chargedParticles.position[i].x,
                    chargedParticles.position[i].y, chargedParticles.position[i].z);
            if (outOfBounds){
                chargedParticles.acceleration[i].w = static_cast<float>(ParticleStatus::outOfBounds);
                writeToEventsFile(i, currentTime, ParticleType::charged);
            }
        }
    }
}

/**
 * Checks to see if a particle is out of bounds in the spatial domain
 * @param x - reference to x position of particle
 * @param y - reference to y position of particle
 * @param z - reference to z position of particle
 * @return - true if out of bounds, false otherwise
 */
int ParticleSystem::checkOutOfSpatialBounds(float& x, float& y, float& z){
    if (params.boundaryConditions.spatialBoundaryConditions.x != 0 && std::abs(x) > params.boundaryConditions.spatialBoundaryConditions.x){
        return true;
    }
    if (params.boundaryConditions.spatialBoundaryConditions.y != 0 && std::abs(y) > params.boundaryConditions.spatialBoundaryConditions.y){
        return true;
    }
    if (params.boundaryConditions.spatialBoundaryConditions.z != 0 && std::abs(z) > params.boundaryConditions.spatialBoundaryConditions.z){
        return true;
    }
    if (params.boundaryConditions.spatialBoundaryConditions.w != 0 && sqrt(x * x + y * y) > params.boundaryConditions.spatialBoundaryConditions.w){
        return true;
    }
//    if ((std::abs(x) > params.boundaryConditions.spatialBoundaryConditions.x) ||
//            (std::abs(y) > params.boundaryConditions.spatialBoundaryConditions.y) ||
//            (std::abs(z) > params.boundaryConditions.spatialBoundaryConditions.z)){
//        return true;
//    }
    return false;
}

/**
 * Enforces temporal boundary conditions for the system of particles
 * @param currentTime - current sim time
 */
void ParticleSystem::enforceTemporalBoundaryConditions(float currentTime){
    int i;
    const float inboundCheck = static_cast<float>(ParticleStatus::inbound);
    int outOfBounds;
    if (params.boundaryConditions.temporalBoundaryCondition == 0){
        printf("Warning: Temporal boundary condition not initialized!\n");
    }
    for (i = chargedParticles.firstParticle; i < chargedParticles.lastParticle; i++){
        if (chargedParticles.acceleration[i].w == inboundCheck){
            outOfBounds = checkOutOfTemporalBounds(chargedParticles.timeInjected[i], currentTime);
            if (outOfBounds){
                chargedParticles.acceleration[i].w = static_cast<float>(ParticleStatus::outOfBounds);
                writeToEventsFile(i, currentTime, ParticleType::charged);
            }
        }
    }
}

/**
 * Checks to see if a particle is out of bounds in the temporal domain
 * @param timeInjected - reference to the time a particle was injected
 * @param currentTime - current sim time
 * @return 
 */
int ParticleSystem::checkOutOfTemporalBounds(float& timeInjected, float currentTime){
    if (currentTime - timeInjected > params.boundaryConditions.temporalBoundaryCondition){
        return true;
    }
    return false;
}

/**
 * Set the index of the first charged particle to the smallest index containing
 *      a particle still considered in the simulation
 */
void ParticleSystem::resetFirstChargedParticle(){
    int i = chargedParticles.firstParticle;
    int firstParticleFound = false;
    while (i < chargedParticles.lastParticle && !firstParticleFound){
        if ((chargedParticles.acceleration[i].w == static_cast<float>(ParticleStatus::affectsAcceleration)) ||
                (chargedParticles.acceleration[i].w == static_cast<float>(ParticleStatus::inbound))){
            firstParticleFound = true;
            chargedParticles.firstParticle = i;
        }
        i++;
    }
}

/**
 * Enforce the symmetry conditions for the particle system
 */
void ParticleSystem::enforceSymmetryConditions(){
    switch (params.symType){
        case (SymmetryType::none):
            break;
        case (SymmetryType::rotational):
            enforceRotationalSymmetry();
            break;
        case (SymmetryType::translational):
            enforceTranslationalSymmetry();
            break;
    }
}

/**
 * Enforce rotational symmetry conditions for the particle system
 */
void ParticleSystem::enforceRotationalSymmetry(){
    int i;
    
    for (i = chargedParticles.firstParticle; i < chargedParticles.lastParticle; i++){
        rotateToFirstQuadrant(chargedParticles.position[i].x, chargedParticles.position[i].y,
                chargedParticles.velocity[i].x, chargedParticles.velocity[i].y);
    }
    
    for (i = neutralParticles.firstParticle; i < neutralParticles.lastParticle; i++){
        rotateToFirstQuadrant(neutralParticles.position[i].x, neutralParticles.position[i].y,
                neutralParticles.velocity[i].x, neutralParticles.velocity[i].y);
    }
}

/**
 * Rotates a particle to lie in the first quadrant, by rotating by a multiple of
 *      pi / 2 radians
 * @param posX - reference to x position of particle
 * @param posY - reference to y position of particle
 * @param velX - reference to x velocity of particle
 * @param velY - reference to y velocity of particle
 */
void ParticleSystem::rotateToFirstQuadrant(float& posX, float& posY, float& velX, float& velY){
    float tempX, tempY;
    if ((posX < 0) && (posY < 0)){
        posX = -posX;
        posY = -posY;
        velX = -velX;
        velY = -velY;
    }
    else if (posX < 0){
        tempX = posY;
        tempY = -posX;
        posX = tempX;
        posY = tempY;

        tempX = velY;
        tempY = -velX;
        velX = tempX;
        velY = tempY;
    }
    else if (posY < 0){
        tempX = -posY;
        tempY = posX;
        posX = tempX;
        posY = tempY;

        tempX = -velY;
        tempY = velX;
        velX = tempX;
        velY = tempY;
    }
}

/**
 * Enforce translational symmetry conditions for the particle system
 */
void ParticleSystem::enforceTranslationalSymmetry(){
    int i;
    
    for (i = chargedParticles.firstParticle; i < chargedParticles.lastParticle; i++){
        translateToSymmetryDomain(chargedParticles.position[i].x,
                chargedParticles.position[i].y, chargedParticles.position[i].z);
    }
    
    for (i = neutralParticles.firstParticle; i < neutralParticles.lastParticle; i++){
        translateToSymmetryDomain(neutralParticles.position[i].x,
                neutralParticles.position[i].y, neutralParticles.position[i].z);
    }
}

/**
 * Translates a particle to lie within the symmetry domain
 * @param posX - reference to x position of particle
 * @param posY - reference to y position of particle
 * @param posZ - reference to z position of particle
 */
void ParticleSystem::translateToSymmetryDomain(float& posX, float& posY, float& posZ){
    float temp;
    if (params.symmetryAxes.x == 1){
        temp = posX + params.domainSize.x / 2;
        if (temp >= 0){
            posX = fmod(temp, params.domainSize.x) - params.domainSize.x / 2;
        }
        else {
            temp = fmod(temp, params.domainSize.x);
            if (temp < 0){
                posX = temp + params.domainSize.x / 2;
            }
            else {
                posX = temp - params.domainSize.x / 2;
            }
        }
    }
    if (params.symmetryAxes.y == 1){
        temp = posY + params.domainSize.y / 2;
        if (temp >= 0){
            posY = fmod(temp, params.domainSize.y) - params.domainSize.y / 2;
        }
        else {
            temp = fmod(temp, params.domainSize.y);
            if (temp < 0){
                posY = temp + params.domainSize.y / 2;
            }
            else {
                posY = temp - params.domainSize.y / 2;
            }
        }
    }
    if (params.symmetryAxes.z == 1){
        temp = posZ + params.domainSize.z / 2;
        if (temp >= 0){
            posZ = fmod(temp, params.domainSize.z) - 
                    params.domainSize.z / 2;
        }
        else {
            temp = fmod(temp, params.domainSize.z);
            if (temp < 0) {
                posZ = temp + params.domainSize.z / 2;
            }
            else {
                posZ = temp - params.domainSize.z / 2;
            }
        }
    }
}

/**
 * Set the injection time for each particle
 * @param injectOptions - option for how the rate is determined
 * @param particleFractions - proportions of particles being injected for simulation
 */
void ParticleSystem::setInjectionTiming(InjectionTimingOptions injectOptions, 
        ParticleFractions particleFractions){
    
    float particleRate = params.I0 / PhysicalConstant::ELEM_CHARGE;
    int i, j, k;
    float tempCharge;
    
    float particleRateFactor = 0;
    for (i = 0; i < particleFractions.maxNumSpecies; i++){
        tempCharge = ParticleSpeciesFunction::getCharge(particleFractions.species[i]);
        tempCharge = std::fabs(tempCharge);
        switch (particleFractions.proportionType){
            case (ParticleProportionOptions::byNumber):
                particleRateFactor += particleFractions.proportion[i] * tempCharge;
                break;
            case (ParticleProportionOptions::byCharge):
                particleRateFactor += particleFractions.proportion[i] / tempCharge;
                break;
        }
    }
    
    switch (particleFractions.proportionType){
        case (ParticleProportionOptions::byNumber):
            particleRate /= particleRateFactor;
            break;
        case (ParticleProportionOptions::byCharge):
            particleRate *= particleRateFactor;
            break;
    }
    
    // Over estimate number of particles required
    int numParticlesEstimate = params.numStartingParticles + 
            1.1 * particleRate * params.totalInjectingTime;
    float* tempTimeInjected = new float[numParticlesEstimate];
    
    switch (injectOptions){
        case (InjectionTimingOptions::constantRate):
            particleRate = 1 / particleRate;
            break;
        case (InjectionTimingOptions::exponentialDistribution):
            randomizerNew.setExponentialDistributionLambda(particleRate);
            break;
    }
    
    for (i = 0; i < params.numStartingParticles; i++){
        tempTimeInjected[i] = 0;
    }
    j = params.numStartingParticles;
    for (i = 0; i < params.sizeInjectionTimes; i++){
        k = 0;
        tempTimeInjected[j] = params.startInjectionTimes[i];
        while (tempTimeInjected[j] < params.stopInjectionTimes[i]){
            k++;
            j++;
//            if (j == 109){
//                printf("test\n");
//            }
            if (j >= numParticlesEstimate){
                float* newTempTimeInjected = new float[numParticlesEstimate * 2];
                std::memcpy(newTempTimeInjected, tempTimeInjected, numParticlesEstimate * sizeof(float));
                numParticlesEstimate *= 2;
                delete[] tempTimeInjected;
                tempTimeInjected = newTempTimeInjected;
            }
            
            switch (injectOptions){
                case (InjectionTimingOptions::constantRate):
                    tempTimeInjected[j] = params.startInjectionTimes[i] + ((float)k )* particleRate;
                    break;
                case (InjectionTimingOptions::exponentialDistribution):
                    tempTimeInjected[j] = tempTimeInjected[j - 1] + randomizerNew.exponentialDistribution();
                    break;
            }
        }
    }
    
    chargedParticles.allocateMemory(j);
    
    std::memcpy(chargedParticles.timeInjected, tempTimeInjected, j * sizeof(float));
    
    delete[] tempTimeInjected;
}

/**
 * Initialize particles in particle system
 * @param initialization - type of initialization
 * @param injectOptions - option for how the rate of injection is determined
 * @param folderName - name of folder from which to read initialization info
 */
void ParticleSystem::initialize(Initialization initialization,
        InjectionTimingOptions injectOptions, const char* folderName){
    
    ParticleFractions particleFractions(initialization);
    setInjectionTiming(injectOptions, particleFractions);
    chargedParticles.initialize(initialization, particleFractions, folderName);
    neutralParticles.initialize(initialization, chargedParticles.maxNeutrals);
}

/**
 * Initialize the background electric field
 * @param bgField - background electric field type
 * @param folderName - name of folder containing necessary interpolation data
 * @param fieldsName - name of file containing interpolation scheme data
 * @param trianglesName - name of file containing triangulation scheme data
 */
void ParticleSystem::initializeBackgroundEField(BackgroundEfield bgField, 
        const char* folderName, const char* fieldsName, const char* trianglesName){
    
    interpScheme.r0 = params.r0;
    interpScheme.E0 = params.E0;
    
    interpScheme.a = params.a;
    interpScheme.nu0 = params.nu0;
    interpScheme.V = params.V;
    interpScheme.changeBackgroundField(bgField);
    interpScheme.initializeInterpolationScheme(folderName, fieldsName, trianglesName);
}

/**
 * Finish iteration
 */
void ParticleSystem::finishIteration(){
//    fragment();
    
}

/**
 * Initialize the fragmentation model used in the simulation
 * @param fragmentationModel - fragmentation model type
 */
void ParticleSystem::initializeFragmentationModel(FragmentationModel fragmentationModel){
    this->fragmentationModel = fragmentationModel;
}

/**
 * Fragment particles
 */
void ParticleSystem::fragment(){
    // If fragment, write to params.eventsFile;
    
    switch (fragmentationModel){
        case (FragmentationModel::none):
            break;
    }
}

/**
 * Compute and update the acceleration due to the background electric field
 */
void ParticleSystem::computeBackgroundAcceleration(){
    int i;
    vec3<float> ptemp;
    float inboundCheck = static_cast<float>(ParticleStatus::inbound);
    const double factor = PhysicalConstant::ELEM_CHARGE / PhysicalConstant::AMU;
    float Efactor;
    
    for (i = chargedParticles.firstParticle; i < chargedParticles.lastParticle; i++){
        if (chargedParticles.acceleration[i].w == inboundCheck){
//            if (i == 1072){
//                printf("test\n");
//            }
            Efactor = factor * chargedParticles.position[i].w / chargedParticles.getMass(i);
            ptemp = interpScheme.Efield(chargedParticles.position[i].x,
                    chargedParticles.position[i].y, chargedParticles.position[i].z);
            
            chargedParticles.Efield2[i] = ptemp;
            
            ptemp *= Efactor;
            chargedParticles.acceleration[i].x += ptemp.x;
            chargedParticles.acceleration[i].y += ptemp.y;
            chargedParticles.acceleration[i].z += ptemp.z;
        }
    }
}

/**
 * Compute the kinetic energy injected due to the n-th charged particle
 * @param n - index of the charged particle being considered
 */
void ParticleSystem::computeChargedInjectionKineticEnergy(int n){
    chargedParticles.injectionKE[n] = computeChargedKineticEnergy(n);
}

/**
 * Compute the kinetic energy of the n-th charged particle
 * @param n - index of the charged particle being considered
 * @return - kinetic energy of the particle
 */
float ParticleSystem::computeChargedKineticEnergy(int n){
    float v;
    v = chargedParticles.velocity[n].x * chargedParticles.velocity[n].x + 
            chargedParticles.velocity[n].y * chargedParticles.velocity[n].y + 
            chargedParticles.velocity[n].z * chargedParticles.velocity[n].z;
    
    return (0.5 * PhysicalConstant::AMU * chargedParticles.getMass(n) * v);
}

/**
 * Compute the potential energy injected due to the n-th charged particle
 * @param n - index of the charged particle being considered
 */
void ParticleSystem::computeChargedInjectionElectricPotentialEnergy(int n){
    chargedParticles.injectionPE[n] = computeChargedElectricPotentialEnergy(n);
}

/**
 * Compute the potential energy of the n-th charged particle
 * @param n - index of the charged particle being considered
 * @return - potential energy of the particle
 */
float ParticleSystem::computeChargedElectricPotentialEnergy(int n){
    float EPE = 0;
    float r;
//    const float ke = 8.9875517923e9;
    const float ke = 1 / (4.0 * PhysicalConstant::PI * PhysicalConstant::FREE_SPACE_PERMITTIVITY);
    const float q2 = PhysicalConstant::ELEM_CHARGE * PhysicalConstant::ELEM_CHARGE;
    
    for (int j = chargedParticles.firstParticle; j < chargedParticles.lastParticle; j++){
        if (n != j){
            r = (chargedParticles.position[n].x - chargedParticles.position[j].x) * 
                    (chargedParticles.position[n].x - chargedParticles.position[j].x) + 
                    (chargedParticles.position[n].y - chargedParticles.position[j].y) * 
                    (chargedParticles.position[n].y - chargedParticles.position[j].y) + 
                    (chargedParticles.position[n].z - chargedParticles.position[j].z) * 
                    (chargedParticles.position[n].z - chargedParticles.position[j].z) + 
                    PhysicalConstant::SOFTENING;
            r = sqrt(r);
            EPE += chargedParticles.position[n].w * chargedParticles.position[j].w / r;
        }
    }
    return (EPE * q2 * ke);
}

/**
 * Compute the field energy injected due to the n-th charged particle
 * @param index - index of the charged particle being considered
 * @param numSteps - number of integration steps to take
 */
void ParticleSystem::computeChargedInjectionFieldEnergy(int index, int numSteps){
    chargedParticles.fieldEnergy[index] = 0;
    if (interpScheme.bgField != BackgroundEfield::none){
        updateChargedFieldEnergy(index, numSteps, interpScheme.referencePosition);
    }
}

/**
 * Update the field energy of a charged particle
 * @param index - index of the charged particle being considered
 * @param numSteps - number of integration steps to take
 * @param initPos - initial position of charged particle
 */
void ParticleSystem::updateChargedFieldEnergy(int index, int numSteps, vec3<float> initPos){
    double tempFieldEnergy_z, tempFieldEnergy_r, tempFieldEnergy;
    vec3<double> dpos, posTemp;
    vec3<float> Etemp;
    int i;
    
    if (interpScheme.bgField != BackgroundEfield::none){
        dpos.x = (chargedParticles.position[index].x - initPos.x);
        dpos.y = (chargedParticles.position[index].y - initPos.y);
        dpos.z = (chargedParticles.position[index].z - initPos.z);
        
        posTemp.x = initPos.x;
        posTemp.y = initPos.y;
        posTemp.z = initPos.z;
        
        tempFieldEnergy_z = 0;
        tempFieldEnergy_r = 0;
        tempFieldEnergy = 0;
        
        for (i = 0; i < numSteps; i++){
            Etemp = interpScheme.Efield(posTemp.x, posTemp.y, posTemp.z + i * dpos.z / numSteps);
            
            tempFieldEnergy_z += Etemp.z * dpos.z / numSteps;
        }
        posTemp.z += dpos.z;
        
        for (i = 0; i < numSteps; i++){
            Etemp = interpScheme.Efield(posTemp.x + i * dpos.x / numSteps, 
                    posTemp.y + i * dpos.y / numSteps, posTemp.z);
            
            tempFieldEnergy_r += Etemp.x * dpos.x / numSteps + Etemp.y * dpos.y / numSteps;
            if (std::isnan(tempFieldEnergy_r)){
                printf("test\n");
            }
        }
        tempFieldEnergy = tempFieldEnergy_z + tempFieldEnergy_r;
        tempFieldEnergy *= -PhysicalConstant::ELEM_CHARGE * chargedParticles.position[index].w;
        
        chargedParticles.fieldEnergy[index] += tempFieldEnergy;
    }
}




/**
 * Propagate the position of particles using their velocity
 * @param dt - time step
 */
void ChargedParticles::propagatePosition(float dt) {
    int i;
    const float inboundCheck = static_cast<float>(ParticleStatus::inbound);
    
    for (i = firstParticle; i < lastParticle; i++){
        if (acceleration[i].w == inboundCheck){
            position[i].x += velocity[i].x * dt;
            position[i].y += velocity[i].y * dt;
            position[i].z += velocity[i].z * dt;
        }
    }
}

/**
 * Write information about an event to an output file
 * @param index - index of particle
 * @param time - time of event
 * @param eventsFile - pointer to event file
 */
void ChargedParticles::writeToEventsFile(int index, float time, FILE* eventsFile){
    fwrite(&time, sizeof(*&time), 1, eventsFile);
    
    fwrite(&position[index].x, sizeof(*&position[index].x), 1, eventsFile);
    fwrite(&position[index].y, sizeof(*&position[index].y), 1, eventsFile);
    fwrite(&position[index].z, sizeof(*&position[index].z), 1, eventsFile);
    
    fwrite(&velocity[index].x, sizeof(*&velocity[index].x), 1, eventsFile);
    fwrite(&velocity[index].y, sizeof(*&velocity[index].y), 1, eventsFile);
    fwrite(&velocity[index].z, sizeof(*&velocity[index].z), 1, eventsFile);
    
    fwrite(&index, sizeof(*&index), 1, eventsFile);
    
    uint32_t temp = velocity[index].w; // Species
    fwrite(&temp, sizeof(*&temp), 1, eventsFile);
    
    temp = acceleration[index].w; // Event
    fwrite(&temp, sizeof(*&temp), 1, eventsFile);
}

/**
 * Propagate the position of all particles in the system using their velocities
 * @param tstart - sim time before propagation
 * @param dt - time step
 * @param crashParameters - parameters defining conditions under which a crash occurs
 * @param eventsFile - pointer to events file
 */
void ChargedParticles::propagatePositionNew(float tstart, float dt,
        CrashParameters crashParameters, FILE* eventsFile) {
    int i;
    float dt_shift;
//    float dt_temp = tend - tstart;
    float tend = tstart + dt;
    vec3<float> initialPos, finalPos;
    float temp_dt;
    
    const float inboundCheck = static_cast<float>(ParticleStatus::inbound);
    const float affectsAccelerationCheck = static_cast<float>(ParticleStatus::affectsAcceleration);
    const float crashStatus = static_cast<float>(ParticleStatus::crashed);
    int crashCheck;
    
    for (i = firstParticle; i < lastParticle; i++){
        if (acceleration[i].w == inboundCheck || acceleration[i].w == affectsAccelerationCheck){
            initialPos.x = position[i].x;
            initialPos.y = position[i].y;
            initialPos.z = position[i].z;
//            if (i == 504){
//                printf("test\n");
//            }
//            dt_shift = std::min(std::max(tend - timeInjected[i], 0.0f), dt_temp);
//            if (tend - timeInjected[i] < 0){
//                printf("test\n");
//            }
            dt_shift = std::min(std::max(tend - timeInjected[i], 0.0f), dt);
            position[i].x += velocity[i].x * dt_shift;
            position[i].y += velocity[i].y * dt_shift;
            position[i].z += velocity[i].z * dt_shift;
            finalPos.x = position[i].x;
            finalPos.y = position[i].y;
            finalPos.z = position[i].z;
            
            crashCheck = crashParameters.checkCrash(initialPos, finalPos);
            if (crashCheck){
                acceleration[i].w = crashStatus;
                position[i].z = finalPos.z;
                temp_dt = (finalPos.z - initialPos.z) / velocity[i].z;
                position[i].x = initialPos.x + temp_dt * velocity[i].x;
                position[i].y = initialPos.y + temp_dt * velocity[i].y;
                
                writeToEventsFile(i, tstart + temp_dt, eventsFile);
            }
        }
    }
}

/**
 * Propagate the velocity of particles using their accelerations
 * @param dt - time step
 */
void ChargedParticles::propagateVelocity(float dt) {
    int i;
    const float inboundCheck = static_cast<float>(ParticleStatus::inbound);
    
    for (i = firstParticle; i < lastParticle; i++){
        if (acceleration[i].w == inboundCheck){
            velocity[i].x += acceleration[i].x * dt;
            velocity[i].y += acceleration[i].y * dt;
            velocity[i].z += acceleration[i].z * dt;
        }
    }
}

/**
 * Propagate the velocity of all particles in the system using their accelerations
 * @param tstart - sim time before propagation
 * @param dt - time step
 */
void ChargedParticles::propagateVelocityNew(float tstart, float dt) {
    int i;
    float dt_shift;
//    float dt_temp = tend - tstart;
    float tend = tstart + dt;
    
    const float inboundCheck = static_cast<float>(ParticleStatus::inbound);
    
    for (i = firstParticle; i < lastParticle; i++){
//        if (i == 1072){
//            printf("test\n");
//        }
        if (acceleration[i].w == inboundCheck){
//            dt_shift = std::min(std::max(tend - timeInjected[i], 0.0f), dt_temp);
            dt_shift = std::min(std::max(tend - timeInjected[i], 0.0f), dt);
//            if (dt_shift != dt){
//                printf("test\n");
//            }
            velocity[i].x += acceleration[i].x * dt_shift;
            velocity[i].y += acceleration[i].y * dt_shift;
            velocity[i].z += acceleration[i].z * dt_shift;
        }
    }
}

/**
 * Set the index of the last particle
 * @param n - index of the last particle
 */
void ChargedParticles::setLastParticle(int n){
    lastParticle = n;
}

/**
 * Set the index of the first particle
 * @param n - index of the first particle
 */
void ChargedParticles::setFirstParticle(int n){
    firstParticle = n;
}

/**
 * Increment the index of the first particle
 * @param dN - amount to increment the index of the first particle by
 */
void ChargedParticles::incrementFirstParticle(int dN){
    firstParticle += dN;
}

/**
 * Increment the index of the last particle
 * @param dN - amount to increment the index of the last particle by
 */
void ChargedParticles::incrementLastParticle(int dN){
    lastParticle += dN;
}

/**
 * Initialize the particles in a system
 * @param initialization - type of initialization
 * @param particleFractions - proportion of each particle type injected
 * @param folderName - name of folder from which to read initialization info
 */
void ChargedParticles::initialize(Initialization initialization,
        ParticleFractions particleFractions, const char* folderName){
    switch (initialization) {
        case (Initialization::paraxialRay):
            initParaxialRay(particleFractions);
            break;
        case (Initialization::singleEmitter):
            if (std::strcmp(folderName, "") == 0){
                printf("\nError: Initialization for single emitter requires input file name.\n");
                exit(1);
            }
            initSingleEmitter(particleFractions, folderName);
            break;
        case (Initialization::none):
            break;
    }
}

/**
 * Initialize particles for a paraxial ray simulation
 * @param particleFractions - proportion of each particle type injected
 */
void ChargedParticles::initParaxialRay(ParticleFractions particleFractions) {
    const float radius = 1e-6;
    int i;
    float r, xTemp, yTemp, charge;
    
    for (i = 0; i < maxNumParticles; i++) {
        do {
            //xTemp = 2 * (float) rand() / (float) RAND_MAX - 1;
            //yTemp = 2 * (float) rand() / (float) RAND_MAX - 1;
            xTemp = 2 * randomizerNew.uniformRealDistribution() - 1;
            yTemp = 2 * randomizerNew.uniformRealDistribution() - 1;
            r = xTemp * xTemp + yTemp*yTemp;
            //charge = 2*(float)rand()/(float)RAND_MAX - 1;
        } while (r >= 1);

        position[i].x = xTemp * radius;
        position[i].y = yTemp * radius;
        position[i].z = 0;
//        if (xTemp > 0) {
//            position[i].w = CHARGE;
//            velocity[i].w = 0;
//            position[i].x += radius;
//            velocity[i].x = -1000 * t_scale / x_scale;
//        } else {
//            position[i].w = -CHARGE;
//            velocity[i].w = 5;
//            position[i].x -= radius;
//            velocity[i].x = 1000 * t_scale / x_scale;
//        }
//        position[i].w = PhysicalConstant::CHARGE;
        velocity[i].x = 0;
        velocity[i].y = 0;
        velocity[i].z = 2000;
//        velocity[i].w = static_cast<float>(ParticleSpecies::monomerPositive); // Particle species
        
        velocity[i].w = static_cast<float>(particleFractions.selectSpecies(randomizerNew.uniformRealDistribution()));
        
        
        position[i].w = getSpeciesCharge(i);
        timesFragmented[i] = 0; //Particles start at 0 fragmentation 
        temperature[i] = getSpeciesTemperature(velocity[i].w);
        crash[i] = 0; //none are crashed at the beginning
        acceleration[i].w = static_cast<float>(ParticleStatus::inbound); //particles are inbound at the beginning
    }
};

/**
 * Initialize the particles for a single emitter simulation
 * @param particleFractions - proportion of each particle type injected
 * @param folderName - name of folder containing initial conditions for particles
 */
void ChargedParticles::initSingleEmitter(ParticleFractions particleFractions, 
        const char* folderName) {
    char* infoname = concat(folderName,"IC_test.txt");
    FILE* ICp = fopen(infoname,"r");
//    const int ARRAY_SIZE = 1201;        // Number of lines in file -- just hardcoded for now
    int size = InfoConstant::ARRAY_SIZE;
    int i, j;
    
    float** IC = new float*[size];
    for (i = 0; i < size; i++){
        IC[i] = new float[7];
    }
    
    for (i = 0; i < size; i++){
        fscanf(ICp,"%g, %g, %g, %g, %g, %g, %g",&IC[i][0],&IC[i][1],&IC[i][2],
                &IC[i][3],&IC[i][4],&IC[i][5],&IC[i][6]);
    }
    
    fclose(ICp);
    
    for (i = 0; i < maxNumParticles; i++){
//        if (i == 1072){
//            printf("test\n");
//        }
        if (i < size) j = i;
        else j = i % size;
        
        position[i].x = IC[j][0];
        position[i].y = IC[j][1];
        position[i].z = IC[j][2];
//        position[i].w = PhysicalConstant::CHARGE;
        velocity[i].x = IC[j][3];
        velocity[i].y = IC[j][4];
        velocity[i].z = IC[j][5];
//        velocity[i].w = static_cast<float>(ParticleSpecies::monomerPositive); // Particle species
        velocity[i].w = static_cast<float>(particleFractions.selectSpecies(randomizerNew.uniformRealDistribution()));
        position[i].w = getSpeciesCharge(i);
        timesFragmented[i] = 0; //Particles start at 0 fragmentation 
        temperature[i] = IC[j][6] * PhysicalConstant::TEMPERATURE;
        crash[i] = 0; //none are crashed at the beginning
        acceleration[i].w = static_cast<float>(ParticleStatus::inbound); //particles are inbound at the beginning
    }
    
    for (i = 0; i < size; i++){
        delete[] IC[i];
    }
    delete[] IC;
};

/**
 * Deconstructor for ChargedParticles
 * Deallocate all memory
 */
ChargedParticles::~ChargedParticles(){
    freeMemory();
}

/**
 * Deallocate all memory
 */
void ChargedParticles::freeMemory(){
    
    if (position != nullptr) delete[] position;
    if (velocity != nullptr) delete[] velocity;
    if (acceleration != nullptr) delete[] acceleration;
    if (temperature != nullptr) delete[] temperature;
    if (timeFlag != nullptr) delete[] timeFlag;
    if (timeInjected != nullptr) delete[] timeInjected;
    if (timesFragmented != nullptr) delete[] timesFragmented;
    if (crash != nullptr) delete[] crash;
    if (Efield1 != nullptr) delete[] Efield1;
    if (Efield2 != nullptr) delete[] Efield2;
    if (injectionKE != nullptr) delete[] injectionKE;
    if (injectionPE != nullptr) delete[] injectionPE;
    if (fieldEnergy != nullptr) delete[] fieldEnergy;
    
    position = nullptr;
    velocity = nullptr;
    acceleration = nullptr;
    temperature = nullptr;
    timeFlag = nullptr;
    timeInjected = nullptr;
    timesFragmented = nullptr;
    crash = nullptr;
    
    Efield1 = nullptr;
    Efield2 = nullptr;
    
    injectionKE = nullptr;
    injectionPE = nullptr;
    fieldEnergy = nullptr;
}

/**
 * Allocate necessary memory for the simulation
 * @param numParticles - number of particles initialized
 */
void ChargedParticles::allocateMemory(int numParticles){
    // Could add a toggle for which memory needs to be allocated
    maxNumParticles = numParticles;
    freeMemory();
    
    // Required
    position = new vec4<float>[maxNumParticles];
    velocity = new vec4<float>[maxNumParticles];
    acceleration = new vec4<float>[maxNumParticles];
    temperature = new float[maxNumParticles];
    timeInjected = new float[maxNumParticles];
    
    // Nice to have
    timesFragmented = new int[maxNumParticles];
    crash = new int[maxNumParticles];
    timeFlag = new float[maxNumParticles];
    Efield1 = new vec3<float>[maxNumParticles];
    Efield2 = new vec3<float>[maxNumParticles];
    
//    // Mainly for debugging
//    allocateMemory_Energy();
}

/**
 * Allocate memory related to energy conservation checks
 */
void ChargedParticles::allocateMemory_Energy(){
    injectionKE = new float[maxNumParticles];
    injectionPE = new float[maxNumParticles];
    fieldEnergy = new float[maxNumParticles];
}

/**
 * Set acceleration of all particles to zero
 */
void ChargedParticles::zeroAcceleration(){
    int i;
    const float inboundCheck = static_cast<float>(ParticleStatus::inbound);
    
    for (i = firstParticle; i < lastParticle; i++){
        if (acceleration[i].w == inboundCheck){
            acceleration[i].x = 0;
            acceleration[i].y = 0;
            acceleration[i].z = 0;
            Efield1[i] = 0;
            Efield2[i] = 0;
        }
    }
}

/**
 * Get the charge of a specific particle
 * @param index - index of particle being considered
 * @return - charge of particle
 */
float ChargedParticles::getCharge(int index) {
    return position[index].w;
}

/**
 * Get the temperature of a specific particle
 * @param index - index of particle being considered
 * @return - temperature of particle
 */
float ChargedParticles::getTemperature(int index){
    if (temperature != nullptr){
        return temperature[index];
    }
    else{
        ParticleSpecies spec = static_cast<ParticleSpecies>(velocity[index].w);
        if (spec != ParticleSpecies::droplet){
            return ParticleSpeciesFunction::getTemperature(spec);
        }
        else {
            // Droplet not supported
            return 0;
        }
    }
}

/**
 * Get the mass of a specific particle
 * @param index - index of particle being considered
 * @return - charge of particle
 */
long double ChargedParticles::getMass(int index) {
    ParticleSpecies spec = static_cast<ParticleSpecies>(velocity[index].w);
    if (spec != ParticleSpecies::droplet){
        return ParticleSpeciesFunction::getMass(spec);
    }
    else {
        // Droplet not supported
        return 0;
    }
}


/**
 * Get the species charge of a specific particle
 * @param index - index of particle being considered
 * @return - charge of particle
 */
float ChargedParticles::getSpeciesCharge(int index){
    ParticleSpecies spec = static_cast<ParticleSpecies>(velocity[index].w);
    if (spec != ParticleSpecies::droplet){
        return ParticleSpeciesFunction::getCharge(spec);
    }
    else {
        // Droplet not supported
        return 0;
    }
}

/**
 * Get the species temperature of a specific particle
 * @param index - index of particle being considered
 * @return - charge of particle
 */
float ChargedParticles::getSpeciesTemperature(int index){
    ParticleSpecies spec = static_cast<ParticleSpecies>(velocity[index].w);
    if (spec != ParticleSpecies::droplet){
        return ParticleSpeciesFunction::getTemperature(spec);
    }
    else {
        // Droplet not supported
        return 0;
    }
}

/**
 * Get the species mass of a specific particle
 * @param index - index of particle being considered
 * @return - charge of particle
 */
long double ChargedParticles::getSpeciesMass(int index){
    ParticleSpecies spec = static_cast<ParticleSpecies>(velocity[index].w);
    if (spec != ParticleSpecies::droplet){
        return ParticleSpeciesFunction::getMass(spec);
    }
    else {
        // Droplet not supported
        return 0;
    }
}

/**
 * Throws an error if try to access an out of bounds index
 * @param lineNumber - file line number of error
 * @param fileName - file name where error occurred
 */
void ChargedParticles::bodyErrorIndexOutOfBounds(int lineNumber, const char* fileName) {
    printf("Error: index out of bounds at line number %d of file %s.\n", lineNumber, fileName);
    exit(1);
}

/**
 * Print the total number of particles in the system
 */
void ChargedParticles::printNumParticles() {
    printf("Number of Particles: %d\n", maxNumParticles);
};

/**
 * Neatly print the position information for the n-th particle
 * @param n - index of particle being considered
 */
void ChargedParticles::printPosInfo(int n) {
    int powTen = 10;
    int i;
    int numSpaces = 0;
    if (n >= 0 && n < maxNumParticles) {
        while (powTen <= maxNumParticles - 1){
            if (n < powTen){
                numSpaces++;
            }
            powTen *= 10;
        }
//        printf("pos[%d].x = %.9e;\tpos[%d].y = %.9e;\tpos[%d].x = %.9e\n",
//                n, pos[n].x, n, pos[n].y, n, pos[n].z);
        printf("position[%d].x = ",n);
        if (std::signbit(position[n].x) == 0){
            printf(" ");
        }
        for (i = 0; i < numSpaces; i++){
            printf(" ");
        }
        printf("%.9e;\t",position[n].x);
        printf("position[%d].y = ",n);
        if (std::signbit(position[n].y) == 0){
            printf(" ");
        }
        for (i = 0; i < numSpaces; i++){
            printf(" ");
        }
        printf("%.9e;\t",position[n].y);
        printf("position[%d].z = ",n);
        if (std::signbit(position[n].z) == 0){
            printf(" ");
        }
        for (i = 0; i < numSpaces; i++){
            printf(" ");
        }
        printf("%.9e;\n",position[n].z);
    } else {
        bodyErrorIndexOutOfBounds(__LINE__, __FILE__);
    }
};

/**
 * Neatly print the velocity information for the n-th particle
 * @param n - index of particle being considered
 */
void ChargedParticles::printVelInfo(int n) {
    int powTen = 10;
    int i;
    int numSpaces = 0;
    if (n >= 0 && n < maxNumParticles) {
        while (powTen <= maxNumParticles - 1){
            if (n < powTen){
                numSpaces++;
            }
            powTen *= 10;
        }

        printf("velocity[%d].x = ",n);
        if (std::signbit(velocity[n].x) == 0){
            printf(" ");
        }
        for (i = 0; i < numSpaces; i++){
            printf(" ");
        }
        printf("%.9e;\t",velocity[n].x);
        printf("velocity[%d].y = ",n);
        if (std::signbit(velocity[n].y) == 0){
            printf(" ");
        }
        for (i = 0; i < numSpaces; i++){
            printf(" ");
        }
        printf("%.9e;\t",velocity[n].y);
        printf("velocity[%d].z = ",n);
        if (std::signbit(velocity[n].z) == 0){
            printf(" ");
        }
        for (i = 0; i < numSpaces; i++){
            printf(" ");
        }
        printf("%.9e;\n",velocity[n].z);
    } else {
        bodyErrorIndexOutOfBounds(__LINE__, __FILE__);
    }
};

/**
 * Neatly print the acceleration information for the n-th particle
 * @param n - index of particle being considered
 */
void ChargedParticles::printAccelInfo(int n) {
    int powTen = 10;
    int i;
    int numSpaces = 0;
    if (n >= 0 && n < maxNumParticles) {
        while (powTen <= maxNumParticles - 1){
            if (n < powTen){
                numSpaces++;
            }
            powTen *= 10;
        }
//        printf("accel[%d].x = %.9e;\taccel[%d].y = %.9e;\taccel[%d].z = %.9e\n",
//                n, accel[n].x, n, accel[n].y, n, accel[n].z);
        printf("acceleration[%d].x = ",n);
        if (std::signbit(acceleration[n].x) == 0){
            printf(" ");
        }
        for (i = 0; i < numSpaces; i++){
            printf(" ");
        }
        printf("%.9e;\t",acceleration[n].x);
        printf("acceleration[%d].y = ",n);
        if (std::signbit(acceleration[n].y) == 0){
            printf(" ");
        }
        for (i = 0; i < numSpaces; i++){
            printf(" ");
        }
        printf("%.9e;\t",acceleration[n].y);
        printf("acceleration[%d].z = ",n);
        if (std::signbit(acceleration[n].z) == 0){
            printf(" ");
        }
        for (i = 0; i < numSpaces; i++){
            printf(" ");
        }
        printf("%.9e;\n",acceleration[n].z);
    } else {
        bodyErrorIndexOutOfBounds(__LINE__, __FILE__);
    }
};

/**
 * Neatly print the injection time information for the n-th particle
 * @param n - index of particle being considered
 */
void ChargedParticles::printTimeInjected(int n){
    int powTen = 10;
    int i;
    int numSpaces = 0;
    if (n >= 0 && n < maxNumParticles) {
        while (powTen <= maxNumParticles - 1){
            if (n < powTen){
                numSpaces++;
            }
            powTen *= 10;
        }
        printf("timeInjected[%d] = ", n);
        if (std::signbit(timeInjected[n]) == 0){
            printf(" ");
        }
        for (i = 0; i < numSpaces; i++){
            printf(" ");
        }
        printf("%.9g\n",timeInjected[n]);
    } else {
        bodyErrorIndexOutOfBounds(__LINE__, __FILE__);
    }
}

/**
 * Finds the index of the most recently injected particle
 * @param time - current time
 * @return - index of the most recently injected particle
 */
int ChargedParticles::findMostRecentInjected(float time){
    int maxIndex = maxNumParticles - 1;
    int minIndex = 0;
    int tempIndex;
    
//    float eps = params.simTime / (maxNumParticles * 1000.0);
    
    do{
        tempIndex = (maxIndex + minIndex) / 2;
//        if (std::abs(timeInjected[tempIndex] - time) < eps){
//            return (tempIndex - 1);
//        }
        if (timeInjected[tempIndex] <= time){
            minIndex = tempIndex + 1;
        }
        else if (timeInjected[tempIndex] > time){
            maxIndex = tempIndex - 1;
        }
        
//        if (timeInjected[tempIndex] == time){
//            printf("test\n");
//        }
//        else if (timeInjected[tempIndex] == time){
//            return tempIndex;
//        }
    } while (maxIndex > minIndex);
    tempIndex = maxIndex;
    if (timeInjected[tempIndex] > time){
        tempIndex--;
    }
    return tempIndex;
}




/**
 * Propagate the position of particles using their velocity
 * @param dt - time step
 */
void NeutralParticles::propagatePosition(float dt) {
    int i;
    
    for (i = firstParticle; i < lastParticle; i++){
        position[i].x += velocity[i].x * dt;
        position[i].y += velocity[i].y * dt;
        position[i].z += velocity[i].z * dt;
    }
}

/**
 * Propagate the position of all particles in the system using their velocities
 * @param tstart - sim time before propagation
 * @param dt - time step
 * @param crashParameters - parameters defining conditions under which a crash occurs
 * @param eventsFile - pointer to events file
 */
void NeutralParticles::propagatePositionNew(float tstart, float dt, 
        CrashParameters crashParameters, FILE* eventsFile) {
    // Add crash if desired
    
    int i;
    float dt_shift;
//    float dt_temp = tend - tstart;
    float tend = tstart + dt;
            
    for (i = firstParticle; i < lastParticle; i++){
//            dt_shift = std::min(std::max(tend - timeInjected[i], 0.0f), dt_temp);
        dt_shift = std::min(std::max(tend - timeCreated[i], 0.0f), dt);
        position[i].x += velocity[i].x * dt_shift;
        position[i].y += velocity[i].y * dt_shift;
        position[i].z += velocity[i].z * dt_shift;
    }
}

/**
 * Propagate the velocity of particles
 * @param dt - time step
 */
void NeutralParticles::propagateVelocity(float dt) {
    // Neutrals currently assumed to experience no acceleration
}

/**
 * Propagate the velocity of all particles in the system using their accelerations
 * @param tstart - sim time before propagation
 * @param dt - time step
 */
void NeutralParticles::propagateVelocityNew(float tstart, float tend) {
    // Neutrals currently assumed to experience no acceleration
}

/**
 * Deallocate memory
 */
void NeutralParticles::freeMemory(){
    
    if (position != nullptr) delete[] position;
    if (velocity != nullptr) delete[] velocity;
    if (timeCreated != nullptr) delete[] timeCreated;
    if (particleIndex != nullptr) delete[] particleIndex;
    
    position = nullptr;
    velocity = nullptr;
    timeCreated = nullptr;
    particleIndex = nullptr;
}

/**
 * Allocate memory
 * @param numParticles - maximum number of particles
 */
void NeutralParticles::allocateMemory(int numParticles){
    maxNumParticles = numParticles;
    freeMemory();
    
    position = new vec3<float>[maxNumParticles];
    velocity = new vec3<float>[maxNumParticles];
    timeCreated = new float[maxNumParticles];
    particleIndex = new int[maxNumParticles];
}

/**
 * Initialize the neutrals in a system
 * @param initialization - type of initialization
 * @param numParticles - maximum number of neutrals
 */
void NeutralParticles::initialize(Initialization initialization, int numParticles){
    
    switch (initialization){
        case (Initialization::none):
            break;
        case (Initialization::paraxialRay):
            break;
        case (Initialization::singleEmitter):
            break;
    }
    
    allocateMemory(numParticles);
}

/**
 * Deconstructor for NeutralParticles
 * Deallocate memory
 */
NeutralParticles::~NeutralParticles(){
    freeMemory();
}