/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "force_calculation_methods.h"
#include "direct_force_method.h"


/**
 * Multi-threaded function to call the direct force computation method for a
 *      rotationally symmetric system  
 * @param dF - struct that holds a pointer to the class and the thread number
 * @return nullptr
 */
void* DirectForceCalculation::directKernelRotation(void* dF) {
    directForceStruct dfc = *((directForceStruct*) dF);
    dfc.d->directKernelRotationCalculation(dfc.threadNum);
    return nullptr;
}

/**
 * Computes the force on all particles in a rotationally symmetric system using
 *      the direct force method
 * @param threadNum - thread number for multi-threaded versions
 */
void DirectForceCalculation::directKernelRotationCalculation(int threadNum){
    for (int k = 0; k < 4; k++){
        if (vectorizedToggle){
            DirectForceFunction::direct_vectorized_rotational(k, rangevalues[threadNum],
                    particleSystem->chargedParticles, jptc);
        }
        else {
            DirectForceFunction::direct_rotational(k, rangevalues[threadNum],
                    particleSystem->chargedParticles);
        }
    }
}

/**
 * Computes the force on all particles using the direct force method for a
 *      rotationally symmetric system
 * @param quadrant - quadrant of rotated cell relative to real cell
 * @param rangevalues - bounds for indices on particles to be included in computing force
 * @param particles - reference to ChargedParticles class being considered
 */
void DirectForceFunction::direct_rotational(int quadrant, vec4<int> rangevalues,
        ChargedParticles& particles) {
    vec3<double> dist, ai;
    double invDist, invDistCube, accel_const, distz2;
    int i, j;
    
    const double eps = PhysicalConstant::SOFTENING;
    const double K_const = PhysicalConstant::K_CONST;
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    const float affectsAccelerationCheck = static_cast<float> (ParticleStatus::affectsAcceleration);
    
    const int cosk = quadrantRotationFunction(quadrant + 1); // = cos(k * pi / 2);
    const int sink = quadrantRotationFunction(quadrant); // = sin(k * pi / 2);

    for (i = rangevalues.x; i < rangevalues.y; i++){
        if (particles.acceleration[i].w == inboundCheck) {
            ai = 0;
            for (j = rangevalues.z; j < rangevalues.w; j++){
                if (particles.acceleration[j].w == inboundCheck ||
                        particles.acceleration[j].w == affectsAccelerationCheck) {
                    dist.z = particles.position[i].z - particles.position[j].z;
                    distz2 = dist.z * dist.z;
//                    for (k = 0; k < 4; k++){
//                        cosk = quadrantRotationFunction(k + 1); // = cos(k * pi / 2);
//                        sink = quadrantRotationFunction(k); // = sin(k * pi / 2);
                        dist.x = particles.position[i].x - (cosk * particles.position[j].x - sink * particles.position[j].y);
                        dist.y = particles.position[i].y - (sink * particles.position[j].x + cosk * particles.position[j].y);
                        invDist = 1.0 / sqrt(dist.x * dist.x + dist.y * dist.y + distz2 + eps);
                        invDistCube = particles.position[j].w * invDist * invDist * invDist;
                        ai.x += dist.x * invDistCube;
                        ai.y += dist.y * invDistCube;
                        ai.z += dist.z * invDistCube;
//                    }
                }
            }
            accel_const = K_const * particles.position[i].w / particles.getMass(i);
            
            ai *= accel_const;
            
            particles.acceleration[i].x += ai.x;
            particles.acceleration[i].y += ai.y;
            particles.acceleration[i].z += ai.z;
            
//            particles.printAccelInfo(i);
        }
    }
}

/**
 * Computes the force on all particles using the direct force method for a 
 *      rotationally symmetric system using vectorized functions
 * @param quadrant - quadrant of rotated cell relative to real cell
 * @param rangevalues - bounds for indices on particles to be included in computing force
 * @param particles - reference to ChargedParticles class being considered
 * @param jptc - pointer to Jpdata array - an aligned array of particle position data
 */
void DirectForceFunction::direct_vectorized_rotational(int quadrant, 
        vec4<int> rangevalues, ChargedParticles& particles, Jpdata* jptc){
    int i, j, ii;
//    volatile int coskf, sinkf, coskb, sinkb;
    int numParticles = rangevalues.w - rangevalues.z;
    
    Ipdata iptc;
    Fodata fout;
    for (j = 0; j < 4; j++){
        fout.eps2[j] = PhysicalConstant::SOFTENING;
    }
    
    volatile const float K_const = PhysicalConstant::K_CONST;
    volatile const float inboundCheck = static_cast<float>(ParticleStatus::inbound);
    
    volatile const int coskf = quadrantRotationFunction(5 - quadrant);
    volatile const int sinkf = quadrantRotationFunction(4 - quadrant);
    volatile const int coskb = quadrantRotationFunction(quadrant + 1);
    volatile const int sinkb = quadrantRotationFunction(quadrant);
    
    int temp[4];
    
    i = rangevalues.x;
    while (i < rangevalues.y){
        j = 0;
        while (j < 4 && i < rangevalues.y){
            if (particles.acceleration[i].w == inboundCheck){
                temp[j] = i;
                j++;
            }
            i++;
        }
        
//        for (k = 0; k < 4; k++){
//            coskf = quadrantRotationFunction(5 - k);
//            sinkf = quadrantRotationFunction(4 - k);
//            coskb = quadrantRotationFunction(k + 1);
//            sinkb = quadrantRotationFunction(k);
            
            for (ii = 0; ii < 4; ii++){
                if (ii < j){
                    iptc.x[ii] = coskf * particles.position[temp[ii]].x - sinkf * particles.position[temp[ii]].y;
                    iptc.y[ii] = sinkf * particles.position[temp[ii]].x + coskf * particles.position[temp[ii]].y;
                    iptc.z[ii] = particles.position[temp[ii]].z;
                    iptc.q[ii] = K_const * particles.position[temp[ii]].w / particles.getMass(temp[ii]);
                }
                else {
                    temp[ii] = 0;
                    iptc.x[ii] = 0;
                    iptc.y[ii] = 0;
                    iptc.z[ii] = 0;
                    iptc.q[ii] = 0;
                }
            }
            
            direct_vectorized_kernel(iptc, fout, jptc, numParticles, rangevalues.z);
            for (ii = 0; ii < j; ii++){
                particles.acceleration[temp[ii]].x += (coskb * fout.ax[ii] - sinkb * fout.ay[ii]);
                particles.acceleration[temp[ii]].y += (sinkb * fout.ax[ii] + coskb * fout.ay[ii]);
                particles.acceleration[temp[ii]].z += fout.az[ii];
            }
//        }
    }
}
