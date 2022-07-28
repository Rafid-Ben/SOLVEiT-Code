/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "force_calculation_methods.h"
#include "direct_force_method.h"

/**
 * Construct DirectForceCalculation class and set numThreads and naming
 * @param particleSystem - class containing information about particles
 */
DirectForceCalculation::DirectForceCalculation(ParticleSystem* particleSystem) :
ForceCalculationMethod{particleSystem}
{
    setNumThreads(particleSystem->params.numThreads);
//    switch (particleSystem->params.symType){
//        case (SymmetryType::none):
//            timerName = "Direct Force Method";
//            saveName = "_direct.bin";
//            break;
//        case (SymmetryType::rotational):
//            timerName = "Rotational Direct Force Method";
//            saveName = "_directRotational.bin";
//            break;
//        case (SymmetryType::translational):
//            timerName = "Translational Direct Force Method";
//            saveName = "_directTranslational.bin";
//            break;
//    }
}

/**
 * Destructor for DirectForceCalculation
 * Deallocate memory associated with threading
 */
DirectForceCalculation::~DirectForceCalculation() {
//    printf("Derived class destructor called!\n");
    deallocateThreadMemory();
    if (jptc != nullptr) delete[] jptc;
    jptc = nullptr;
}

/**
 * Allocate memory associated with threading
 */
void DirectForceCalculation::allocateThreadMemory() {
#ifdef _PTHREAD_H
    if (dF != nullptr) delete[] dF;
    dF = new directForceStruct[numThreads];
#endif
}

/**
 * Deallocate memory associated with threading
 */
void DirectForceCalculation::deallocateThreadMemory() {
#ifdef _PTHREAD_H
    if (dF != nullptr) delete[] dF;
    dF = nullptr;
#endif
}

/**
 * Toggle if vectorization will be used to evaluate computations and allocate 
 *  necessary memory
 * @param toggle - true if vectorized functions will be used, false otherwise
 */
 void DirectForceCalculation::setVectorizedToggle(int toggle){
    vectorizedToggle = toggle;
    if (jptc != nullptr) delete[] jptc;
    jptc = nullptr;
    if (vectorizedToggle == true){
        jptc = new Jpdata[particleSystem->chargedParticles.maxNumParticles + 2];
    }
}

/**
 * Multi-threaded function to call the direct force computation method
 * @param dF - struct that holds a pointer to the class and the thread number
 * @return nullptr
 */
void* DirectForceCalculation::directKernel(void* dF) {
    directForceStruct dfc = *((directForceStruct*) dF);
    dfc.d->directKernelCalculation(dfc.threadNum);
    return nullptr;
}

//void DirectForceFunction::direct(int i_min, int i_max, int j_min, int j_max, 
//        ChargedParticles& particles){
//    vec3<double> dist, ai;
//    double invDist, invDistCube, accel_const;
//    int i, j;
//    
//    const double eps = PhysicalConstant::SOFTENING;
//    const double K_const = PhysicalConstant::K_CONST;
//    const float inboundCheck = static_cast<float>(ParticleStatus::inbound);
//    const float affectsAccelerationCheck = static_cast<float>(ParticleStatus::affectsAcceleration);
//    
//    for (i = i_min; i < i_max; i++){
//        if (particles.acceleration[i].w == inboundCheck){
//            ai = 0;
//            for (j = j_min; j < j_max; j++){
//                if (particles.acceleration[j].w == affectsAccelerationCheck){
//                    if (j != i){
//                        dist.x = particles.position[i].x - particles.position[j].x;
//                        dist.y = particles.position[i].y - particles.position[j].y;
//                        dist.z = particles.position[i].z - particles.position[j].z;            
//                        invDist = 1.0 / sqrt(dist.x * dist.x + dist.y * dist.y + dist.z * dist.z + eps);
//                        invDistCube = particles.position[j].w * invDist * invDist * invDist;
//                        ai.x += dist.x * invDistCube;
//                        ai.y += dist.y * invDistCube;
//                        ai.z += dist.z * invDistCube;
//                    }
//                }
//            }
//            accel_const = K_const * particles.position[i].w / particles.getMass(i);
//
//            ai *= accel_const;
//
//            particles.acceleration[i].x += ai.x;
//            particles.acceleration[i].y += ai.y;
//            particles.acceleration[i].z += ai.z;
//        }
//    }
//}

void DirectForceFunction::setJpdata(int& tempNumParticles, 
        ChargedParticles& particles, Jpdata* jptc){
    int i, j;
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    const float affectsAccelerationCheck = static_cast<float> (ParticleStatus::affectsAcceleration);
    
    i = 0;
    j = 0;
    for (i = particles.firstParticle; i < particles.lastParticle; i++){
        if (particles.acceleration[i].w == inboundCheck || 
                particles.acceleration[j].w == affectsAccelerationCheck) {
            jptc[j] = {particles.position[i].x, particles.position[i].y,
                particles.position[i].z, particles.position[i].w};
            j++;
        }
    }
    tempNumParticles = j;
}

/**
 * Computes the force on all particles using the direct force method
 * @param rangevalues - bounds for indices on particles to be included in computing force
 * @param particles - reference to ChargedParticles class being considered
 */
void DirectForceFunction::direct(vec4<int> rangevalues, ChargedParticles& particles) {
    vec3<double> dist, ai;
    double invDist, invDistCube, accel_const;
    int i, j;

    const double eps = PhysicalConstant::SOFTENING;
    const double K_const = PhysicalConstant::K_CONST;
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    const float affectsAccelerationCheck = static_cast<float> (ParticleStatus::affectsAcceleration);

    for (i = rangevalues.x; i < rangevalues.y; i++) {
        if (particles.acceleration[i].w == inboundCheck) {
            ai = 0;
            for (j = rangevalues.z; j < rangevalues.w; j++) {
                if (particles.acceleration[j].w == inboundCheck || 
                        particles.acceleration[j].w == affectsAccelerationCheck) {
                    if (j != i) {
                        dist.x = particles.position[i].x - particles.position[j].x;
                        dist.y = particles.position[i].y - particles.position[j].y;
                        dist.z = particles.position[i].z - particles.position[j].z;
                        invDist = 1.0 / sqrt(dist.x * dist.x + dist.y * dist.y + dist.z * dist.z + eps);
                        invDistCube = particles.position[j].w * invDist * invDist * invDist;
                        ai.x += dist.x * invDistCube;
                        ai.y += dist.y * invDistCube;
                        ai.z += dist.z * invDistCube;
                    }
                }
            }
            accel_const = K_const * particles.position[i].w / particles.getMass(i);

            ai *= accel_const;

            particles.acceleration[i].x += ai.x;
            particles.acceleration[i].y += ai.y;
            particles.acceleration[i].z += ai.z;
        }
    }
}

/**
 * Computes the force on all particles using the direct force method using vectorized functions
 * @param rangevalues - bounds for indices on particles to be included in computing force
 * @param particles - reference to ChargedParticles class being considered
 * @param jptc - pointer to Jpdata array - an aligned array of particle position data
 */
void DirectForceFunction::direct_vectorized(vec4<int> rangevalues, ChargedParticles& particles,
        Jpdata* jptc){
    int i, j, k;//, offset;
    int numParticles = rangevalues.w - rangevalues.z;
    int numSubsetParticles = rangevalues.y - rangevalues.x;
    
    Ipdata iptc;
    Fodata fout;
    for (j = 0; j < 4; j++){
        fout.eps2[j] = PhysicalConstant::SOFTENING;
    }
    
    volatile const float K_const = PhysicalConstant::K_CONST;
    volatile const float inboundCheck = static_cast<float>(ParticleStatus::inbound);
    
//    int numFullVecs = numSubsetParticles / 4;
//    offset = rangevalues.x;
    
    int temp[4];
    
    i = rangevalues.x;
    while (i < rangevalues.y){
        j = 0;
        while (j < 4 && i < rangevalues.y){
            if (particles.acceleration[i].w == inboundCheck){
                temp[j] = i;
                iptc.x[j] = particles.position[i].x;
                iptc.y[j] = particles.position[i].y;
                iptc.z[j] = particles.position[i].z;
                iptc.q[j] = K_const * particles.position[i].w / particles.getMass(i);
                j++;
            }
            
            i++;
        }
        for (k = j; k < 4; k++){
            temp[k] = 0;
            iptc.x[k] = 0.0f;
            iptc.y[k] = 0.0f;
            iptc.z[k] = 0.0f;
            iptc.q[k] = 0.0f;
        }
        direct_vectorized_kernel(iptc, fout, jptc, numParticles, rangevalues.z);
        for (k = 0; k < j; k++){
            particles.acceleration[temp[k]].x += fout.ax[k];
            particles.acceleration[temp[k]].y += fout.ay[k];
            particles.acceleration[temp[k]].z += fout.az[k];
        }
    }
    
//    for (i = 0; i < numFullVecs; i++){
//        for (j = 0; j < 4; j++){
//            iptc.x[j] = particles.position[offset + j].x;
//            iptc.y[j] = particles.position[offset + j].y;
//            iptc.z[j] = particles.position[offset + j].z;
//            iptc.q[j] = K_const * particles.position[offset + j].w / particles.getMass(offset + j);
//        }
//        direct_vectorized_kernel(iptc, fout, jptc, numParticles, rangevalues.z);
//        for (j = 0; j < 4; j++){
//            particles.acceleration[offset + j].x = fout.ax[j];
//            particles.acceleration[offset + j].y = fout.ay[j];
//            particles.acceleration[offset + j].z = fout.az[j];
//        }
//        offset += 4;
//    }
//    
//    int numLeftover = numSubsetParticles - 4 * numFullVecs;
//    
//    if (numLeftover != 0){
//        for (j = 0; j < numLeftover; j++){
//            iptc.x[j] = particles.position[offset + j].x;
//            iptc.y[j] = particles.position[offset + j].y;
//            iptc.z[j] = particles.position[offset + j].z;
//            iptc.q[j] = K_const * particles.position[offset + j].w / particles.getMass(offset + j);
//        }
//        for (j = numLeftover; j < 4; j++){
//            iptc.x[j] = 0.0f;
//            iptc.y[j] = 0.0f;
//            iptc.z[j] = 0.0f;
//            iptc.q[j] = 0.0f;
//        }
//        direct_vectorized_kernel(iptc, fout, jptc, numParticles, rangevalues.z);
//        for (j = 0; j < numLeftover; j++){
//            particles.acceleration[offset + j].x = fout.ax[j];
//            particles.acceleration[offset + j].y = fout.ay[j];
//            particles.acceleration[offset + j].z = fout.az[j];
//        }
//    }
}

/**
 * Computes the force on all particles using the direct force method
 * @param threadNum - thread number for multi-threaded versions
 */
void DirectForceCalculation::directKernelCalculation(int threadNum) {

    if (vectorizedToggle){
        DirectForceFunction::direct_vectorized(rangevalues[threadNum],
                particleSystem->chargedParticles, jptc);
    }
    else {
        DirectForceFunction::direct(rangevalues[threadNum],
                particleSystem->chargedParticles);
    }

    //    vec3<double> dist;
    //    double invDist, invDistCube, accel_const;
    //    int i, j;
    //    
    //    const double eps = PhysicalConstant::SOFTENING;
    //    const double K_const = PhysicalConstant::K_CONST;
    //    
    //    for (i = rangevalues[threadNum].x; i < rangevalues[threadNum].y; i++) {
    //        if (particleSystem->chargedParticles.acceleration[i].w == 1){
    //            vec3<double> ai = {0.0, 0.0, 0.0};
    //            for (j = rangevalues[threadNum].z; j < rangevalues[threadNum].w; j++) {
    //                if (j != i){
    //                    dist.x = particleSystem->chargedParticles.position[i].x - particleSystem->chargedParticles.position[j].x;
    //                    dist.y = particleSystem->chargedParticles.position[i].y - particleSystem->chargedParticles.position[j].y;
    //                    dist.z = particleSystem->chargedParticles.position[i].z - particleSystem->chargedParticles.position[j].z;            
    //                    invDist = 1.0 / sqrt(dist.x * dist.x + dist.y * dist.y + dist.z * dist.z + eps);
    //                    invDistCube = particleSystem->chargedParticles.position[j].w * invDist * invDist * invDist;
    //                    ai.x += dist.x * invDistCube;
    //                    ai.y += dist.y * invDistCube;
    //                    ai.z += dist.z * invDistCube;
    //                }
    //            }
    //    //        accel_const = K_const * particleSystem->chargedParticles.position[i].w / getMass(particleSystem->chargedParticles.velocity[i].w);
    //            accel_const = K_const * particleSystem->chargedParticles.position[i].w / particleSystem->chargedParticles.getMass(i);
    //
    //            ai.x *= accel_const;
    //            ai.y *= accel_const;
    //            ai.z *= accel_const;
    //
    //            particleSystem->chargedParticles.acceleration[i].x = ai.x;
    //            particleSystem->chargedParticles.acceleration[i].y = ai.y;
    //            particleSystem->chargedParticles.acceleration[i].z = ai.z;
    //            
    ////            particleSystem->chargedParticles.printAccelInfo(i);
    ////            printf("norm: %.9g\n",sqrt(particleSystem->chargedParticles.acceleration[i].x * particleSystem->chargedParticles.acceleration[i].x +
    ////                    particleSystem->chargedParticles.acceleration[i].y * particleSystem->chargedParticles.acceleration[i].y + 
    ////                    particleSystem->chargedParticles.acceleration[i].z * particleSystem->chargedParticles.acceleration[i].z));
    //        }
    //    }
}

/**
 * Updates the acceleration of all particles due to particle-particle interactions
 */
void DirectForceCalculation::updateAcceleration() {
    int k;
    setTimerName();
    
    if (particleSystem == nullptr){
        printf("Particle system not initialized.\n");
        exit(1);
    }

    int numParticles = lastParticle - firstParticle;

    if (timeInfoFile != nullptr){
        fwrite(&numParticles, sizeof (*&numParticles), 1, timeInfoFile);
    }

    Timer t;
    t.startTimer();
    
    if (vectorizedToggle){
        int tempNumParticles;
        DirectForceFunction::setJpdata(tempNumParticles, particleSystem->chargedParticles, jptc);
        for (k = 0; k < numThreads; k++){
            rangevalues[k].z = 0;
            rangevalues[k].w = tempNumParticles;
        }
    }
    else {
        for (k = 0; k < numThreads; k++){
            rangevalues[k].z = firstParticle;
            rangevalues[k].w = lastParticle;
        }
    }

#ifdef _PTHREAD_H
    //    int numParticles = lastParticle - firstParticle;
    for (k = 0; k < numThreads; k++) {
        dF[k].d = this;
        dF[k].threadNum = k;
        rangevalues[k].x = (k * numParticles) / numThreads + firstParticle;
        rangevalues[k].y = ((k + 1) * numParticles) / numThreads + firstParticle;
//        rangevalues[k].z = firstParticle;
//        rangevalues[k].w = lastParticle;
        switch (particleSystem->params.symType) {
            case (SymmetryType::none):
                pthread_create(&threads[k], nullptr, directKernel, &dF[k]);
                break;
            case (SymmetryType::rotational):
                pthread_create(&threads[k], nullptr, directKernelRotation, &dF[k]);
                break;
            case (SymmetryType::translational):
                pthread_create(&threads[k], nullptr, directKernelTranslation, &dF[k]);
                break;
        }
    }

    for (k = 0; k < numThreads; k++) {
        pthread_join(threads[k], NULL);
    }
#else
    rangevalues[0].x = firstParticle;
    rangevalues[0].y = lastParticle;
//    rangevalues[0].z = firstParticle;
//    rangevalues[0].w = lastParticle;
    switch (particleSystem->params.symType) {
        case (SymmetryType::none):
            directKernelCalculation();
            break;
        case (SymmetryType::rotational):
            directKernelRotationCalculation();
            break;
        case (SymmetryType::translational):
            directKernelTranslationCalculation();
            break;
    }
#endif

    float accelTime = t.getTimer();
    
    if (timeInfoFile != nullptr){
        fwrite(&accelTime, sizeof (*&accelTime), 1, timeInfoFile);
    }
}

