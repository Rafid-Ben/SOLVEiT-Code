/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "force_calculation_methods.h"
#include "direct_force_method.h"


/**
 * Multi-threaded function to call the direct force computation method for a
 *      translationally symmetric system  
 * @param dF - struct that holds a pointer to the class and the thread number
 * @return nullptr
 */
void* DirectForceCalculation::directKernelTranslation(void* dF) {
    directForceStruct dfc = *((directForceStruct*) dF);
    dfc.d->directKernelTranslationCalculation(dfc.threadNum);
    return nullptr;
}

/**
 * Computes the force on all particles in a translationally symmetric system using
 *      the direct force method
 * @param threadNum - thread number for multi-threaded versions
 */
void DirectForceCalculation::directKernelTranslationCalculation(int threadNum){
    int i, j, k, ii, jj;
//    int numBoxes = 0;
    
    LatticePoint tempLattice;
//    volatile vec3<float> domainOffset;
    
    for (ii = 0; ii < grid->numLattice; ii++){
        tempLattice.setPoint(0, 0, 0);
        
        if (grid->symDimSum == 3){
            if (grid->lattice[ii].point.x != 0 && grid->lattice[ii].point.y != 0 && grid->lattice[ii].point.z != 0){
                for (i = 0; i < 2; i++){
                    for (j = 0; j < 2; j++){
                        for (k = 0; k < 2; k++){
                            tempLattice.point.x = ODDEVEN(i) * grid->lattice[ii].point.x;
                            tempLattice.point.y = ODDEVEN(j) * grid->lattice[ii].point.y;
                            tempLattice.point.z = ODDEVEN(k) * grid->lattice[ii].point.z;
//                            numBoxes++;
//                            tempLattice.print();
                            
//                            domainOffset.x = tempLattice.point.x * particleSystem->params.domainSize.x;
//                            domainOffset.y = tempLattice.point.y * particleSystem->params.domainSize.y;
//                            domainOffset.z = tempLattice.point.z * particleSystem->params.domainSize.z;

                            if (vectorizedToggle){
                                DirectForceFunction::direct_vectorized_translational(tempLattice, particleSystem->params.domainSize,
                                        rangevalues[threadNum], particleSystem->chargedParticles, jptc);
                            }
                            else {
                                DirectForceFunction::direct_translational(tempLattice, particleSystem->params.domainSize,
                                        rangevalues[threadNum], particleSystem->chargedParticles);
                            }
                        }
                    }
                }
            }
            else if ((grid->lattice[ii].point.x == 0 && grid->lattice[ii].point.y != 0 && grid->lattice[ii].point.z != 0) || 
                (grid->lattice[ii].point.x != 0 && grid->lattice[ii].point.y == 0 && grid->lattice[ii].point.z != 0) ||
                (grid->lattice[ii].point.x != 0 && grid->lattice[ii].point.y != 0 && grid->lattice[ii].point.z == 0)){
                for (jj = 0; jj < 3; jj++){
                    for (i = 0; i < (1 + (grid->lattice[ii].point.x != 0)); i++){
                        for (j = 0; j < (1 + (grid->lattice[ii].point.y != 0)); j++){
                            for (k = 0; k < (1 + (grid->lattice[ii].point.z != 0)); k++){
                                tempLattice.point.x = (jj != 0) * ODDEVEN(i) * grid->lattice[ii].point.x + 
                                        (grid->lattice[ii].point.x == 0) * ((jj == 1) * ODDEVEN(j) * grid->lattice[ii].point.y +
                                        (jj == 2) * ODDEVEN(k) * grid->lattice[ii].point.z);
                                tempLattice.point.y = (jj != 1) * ODDEVEN(j) * grid->lattice[ii].point.y + 
                                        (grid->lattice[ii].point.y == 0) * ((jj == 0) * ODDEVEN(i) * grid->lattice[ii].point.x + 
                                        (jj == 2) * ODDEVEN(k) * grid->lattice[ii].point.z);
                                tempLattice.point.z = (jj != 2) * ODDEVEN(k) * grid->lattice[ii].point.z + 
                                        (grid->lattice[ii].point.z == 0) * ((jj == 0) * ODDEVEN(i) * grid->lattice[ii].point.x + 
                                        (jj == 1) * ODDEVEN(j) * grid->lattice[ii].point.y);
//                                numBoxes++;
//                                tempLattice.print();
                                
//                                domainOffset.x = tempLattice.point.x * particleSystem->params.domainSize.x;
//                                domainOffset.y = tempLattice.point.y * particleSystem->params.domainSize.y;
//                                domainOffset.z = tempLattice.point.z * particleSystem->params.domainSize.z;
                                
                                if (vectorizedToggle){
                                    DirectForceFunction::direct_vectorized_translational(tempLattice, particleSystem->params.domainSize,
                                            rangevalues[threadNum], particleSystem->chargedParticles, jptc);
                                }
                                else {
                                    DirectForceFunction::direct_translational(tempLattice, particleSystem->params.domainSize,
                                            rangevalues[threadNum], particleSystem->chargedParticles);
                                }
                            }
                        }
                    }
                }
            }
            else if ((grid->lattice[ii].point.x == 0 && grid->lattice[ii].point.y == 0 && grid->lattice[ii].point.z != 0) || 
                (grid->lattice[ii].point.x == 0 && grid->lattice[ii].point.y != 0 && grid->lattice[ii].point.z == 0) ||
                (grid->lattice[ii].point.x != 0 && grid->lattice[ii].point.y == 0 && grid->lattice[ii].point.z == 0)){
                for (jj = 0; jj < 3; jj++){
                    for (i = 0; i < (1 + (grid->lattice[ii].point.x != 0)); i++){
                        for (j = 0; j < (1 + (grid->lattice[ii].point.y != 0)); j++){
                            for (k = 0; k < (1 + (grid->lattice[ii].point.z != 0)); k++){
                                tempLattice.point.x = (jj == 0) * ODDEVEN(i) * grid->lattice[ii].point.x + 
                                        (jj == 1) * ODDEVEN(j) * grid->lattice[ii].point.y +
                                        (jj == 2) * ODDEVEN(k) * grid->lattice[ii].point.z;
                                tempLattice.point.y = (jj == 0) * ODDEVEN(j) * grid->lattice[ii].point.y + 
                                        (jj == 1) * ODDEVEN(k) * grid->lattice[ii].point.z + 
                                        (jj == 2) * ODDEVEN(i) * grid->lattice[ii].point.x;
                                tempLattice.point.z = (jj == 0) * ODDEVEN(k) * grid->lattice[ii].point.z + 
                                        (jj == 1) * ODDEVEN(i) * grid->lattice[ii].point.x + 
                                        (jj == 2) * ODDEVEN(j) * grid->lattice[ii].point.y;
//                                numBoxes++;
//                                tempLattice.print();
                                
//                                domainOffset.x = tempLattice.point.x * particleSystem->params.domainSize.x;
//                                domainOffset.y = tempLattice.point.y * particleSystem->params.domainSize.y;
//                                domainOffset.z = tempLattice.point.z * particleSystem->params.domainSize.z;
                                
                                if (vectorizedToggle){
                                    DirectForceFunction::direct_vectorized_translational(tempLattice, particleSystem->params.domainSize,
                                            rangevalues[threadNum], particleSystem->chargedParticles, jptc);
                                }
                                else {
                                    DirectForceFunction::direct_translational(tempLattice, particleSystem->params.domainSize,
                                            rangevalues[threadNum], particleSystem->chargedParticles);
                                }
                            }
                        }
                    }
                }
            }
            else {
//                numBoxes++;
//                tempLattice.print();
                
//                domainOffset.x = tempLattice.point.x * particleSystem->params.domainSize.x;
//                domainOffset.y = tempLattice.point.y * particleSystem->params.domainSize.y;
//                domainOffset.z = tempLattice.point.z * particleSystem->params.domainSize.z;
                
                if (vectorizedToggle){
                    DirectForceFunction::direct_vectorized_translational(tempLattice, particleSystem->params.domainSize,
                            rangevalues[threadNum], particleSystem->chargedParticles, jptc);
                }
                else {
                    DirectForceFunction::direct_translational(tempLattice, particleSystem->params.domainSize,
                            rangevalues[threadNum], particleSystem->chargedParticles);
                }
            }
        }
        else if (grid->symDimSum == 2){
            if ((grid->lattice[ii].point.x == 0 && grid->lattice[ii].point.y != 0 && grid->lattice[ii].point.z != 0) || 
                (grid->lattice[ii].point.x != 0 && grid->lattice[ii].point.y == 0 && grid->lattice[ii].point.z != 0) ||
                (grid->lattice[ii].point.x != 0 && grid->lattice[ii].point.y != 0 && grid->lattice[ii].point.z == 0)){
                for (i = 0; i < (1 + (grid->lattice[ii].point.x != 0)); i++){
                    for (j = 0; j < (1 + (grid->lattice[ii].point.y != 0)); j++){
                        for (k = 0; k < (1 + (grid->lattice[ii].point.z != 0)); k++){
                            tempLattice.point.x = ODDEVEN(i) * grid->lattice[ii].point.x;
                            tempLattice.point.y = ODDEVEN(j) * grid->lattice[ii].point.y;
                            tempLattice.point.z = ODDEVEN(k) * grid->lattice[ii].point.z;
//                            numBoxes++;
//                            tempLattice.print();
                            
//                            domainOffset.x = tempLattice.point.x * particleSystem->params.domainSize.x;
//                            domainOffset.y = tempLattice.point.y * particleSystem->params.domainSize.y;
//                            domainOffset.z = tempLattice.point.z * particleSystem->params.domainSize.z;
                            
//                            printf("domainOffset.x = %.9g, .y = %.9g, .z = %.9g\n",domainOffset.x, domainOffset.y, domainOffset.z);
//                            printf(".x = %d, .x = %.9g, .x = %.9g\n",tempLattice.point.x, particleSystem->params.domainSize.x, ((float)tempLattice.point.x) * ((float)particleSystem->params.domainSize.x));
                
                            
                            if (vectorizedToggle){
                                DirectForceFunction::direct_vectorized_translational(tempLattice, particleSystem->params.domainSize,
                                        rangevalues[threadNum], particleSystem->chargedParticles, jptc);
                            }
                            else {
                                DirectForceFunction::direct_translational(tempLattice, particleSystem->params.domainSize,
                                        rangevalues[threadNum], particleSystem->chargedParticles);
                            }
                        }
                    }
                }
            }
            else if ((grid->lattice[ii].point.x == 0 && grid->lattice[ii].point.y == 0 && grid->lattice[ii].point.z != 0) || 
                (grid->lattice[ii].point.x == 0 && grid->lattice[ii].point.y != 0 && grid->lattice[ii].point.z == 0) ||
                (grid->lattice[ii].point.x != 0 && grid->lattice[ii].point.y == 0 && grid->lattice[ii].point.z == 0)){
                for (jj = 0; jj < 2; jj++){
                    for (i = 0; i < (1 + (grid->lattice[ii].point.x != 0)); i++){
                        for (j = 0; j < (1 + (grid->lattice[ii].point.y != 0)); j++){
                            for (k = 0; k < (1 + (grid->lattice[ii].point.z != 0)); k++){
                                tempLattice.point.x = grid->symmetryAxes.x * ((jj == 0) * ODDEVEN(i) * grid->lattice[ii].point.x + 
                                        (jj == 1) * (ODDEVEN(j) * grid->lattice[ii].point.y + ODDEVEN(k) * grid->lattice[ii].point.z));
                                tempLattice.point.y = grid->symmetryAxes.y * ((jj == 0) * ODDEVEN(j) * grid->lattice[ii].point.y + 
                                        (jj == 1) * (ODDEVEN(i) * grid->lattice[ii].point.x + ODDEVEN(k) * grid->lattice[ii].point.z));
                                tempLattice.point.z = grid->symmetryAxes.z * ((jj == 0) * ODDEVEN(k) * grid->lattice[ii].point.z + 
                                        (jj == 1) * (ODDEVEN(i) * grid->lattice[ii].point.x + ODDEVEN(j) * grid->lattice[ii].point.y));
//                                numBoxes++;
//                                tempLattice.print();
                                
//                                domainOffset.x = tempLattice.point.x * particleSystem->params.domainSize.x;
//                                domainOffset.y = tempLattice.point.y * particleSystem->params.domainSize.y;
//                                domainOffset.z = tempLattice.point.z * particleSystem->params.domainSize.z;
                                
//                                printf("domainOffset.x = %.9g, .y = %.9g, .z = %.9g\n",domainOffset.x, domainOffset.y, domainOffset.z);
//                                printf(".x = %d, .x = %.9g, .x = %.9g\n",tempLattice.point.x, particleSystem->params.domainSize.x, ((float)tempLattice.point.x) * ((float)particleSystem->params.domainSize.x));
                
                                
                                if (vectorizedToggle){
                                    DirectForceFunction::direct_vectorized_translational(tempLattice, particleSystem->params.domainSize,
                                            rangevalues[threadNum], particleSystem->chargedParticles, jptc);
                                }
                                else {
                                    DirectForceFunction::direct_translational(tempLattice, particleSystem->params.domainSize,
                                            rangevalues[threadNum], particleSystem->chargedParticles);
                                }
                            }
                        }
                    }
                }
            }
            else {
                tempLattice.point.x = 0;
                tempLattice.point.y = 0;
                tempLattice.point.z = 0;
//                numBoxes++;
//                tempLattice.print();
                
//                domainOffset.x = (float)tempLattice.point.x * particleSystem->params.domainSize.x;
//                domainOffset.y = tempLattice.point.y * particleSystem->params.domainSize.y;
//                domainOffset.z = tempLattice.point.z * particleSystem->params.domainSize.z;
                
//                printf("domainOffset.x = %.9g, .y = %.9g, .z = %.9g\n",domainOffset.x, domainOffset.y, domainOffset.z);
//                printf(".x = %d, .x = %.9g, .x = %.9g\n",tempLattice.point.x, particleSystem->params.domainSize.x, ((float)tempLattice.point.x) * ((float)particleSystem->params.domainSize.x));
                
//                float test_float = 10;
//                int test_int = 0;
//                float test_result = test_float * test_int;
//                float test_result = 0.0f * 0.0f;
//                printf("test_result = %.9g\n",test_result);
                
//                if (tempLattice.point.y * particleSystem->params.domainSize.y != 0){
//                    printf("Error\n");
//                }
                
                if (vectorizedToggle){
                    DirectForceFunction::direct_vectorized_translational(tempLattice, particleSystem->params.domainSize,
                            rangevalues[threadNum], particleSystem->chargedParticles, jptc);
                }
                else {
                    DirectForceFunction::direct_translational(tempLattice, particleSystem->params.domainSize,
                            rangevalues[threadNum], particleSystem->chargedParticles);
                }
            }
        }
        else if (grid->symDimSum == 1){
            if ((grid->lattice[ii].point.x == 0 && grid->lattice[ii].point.y == 0 && grid->lattice[ii].point.z != 0) || 
                (grid->lattice[ii].point.x == 0 && grid->lattice[ii].point.y != 0 && grid->lattice[ii].point.z == 0) ||
                (grid->lattice[ii].point.x != 0 && grid->lattice[ii].point.y == 0 && grid->lattice[ii].point.z == 0)){
                for (i = 0; i < (1 + (grid->lattice[ii].point.x != 0)); i++){
                    for (j = 0; j < (1 + (grid->lattice[ii].point.y != 0)); j++){
                        for (k = 0; k < (1 + (grid->lattice[ii].point.z != 0)); k++){
                            tempLattice.point.x = ODDEVEN(i) * grid->lattice[ii].point.x;
                            tempLattice.point.y = ODDEVEN(j) * grid->lattice[ii].point.y;
                            tempLattice.point.z = ODDEVEN(k) * grid->lattice[ii].point.z;
//                            numBoxes++;
//                            tempLattice.print();
                            
//                            domainOffset.x = tempLattice.point.x * particleSystem->params.domainSize.x;
//                            domainOffset.y = tempLattice.point.y * particleSystem->params.domainSize.y;
//                            domainOffset.z = tempLattice.point.z * particleSystem->params.domainSize.z;
                            
                            if (vectorizedToggle){
                                DirectForceFunction::direct_vectorized_translational(tempLattice, particleSystem->params.domainSize,
                                        rangevalues[threadNum], particleSystem->chargedParticles, jptc);
                            }
                            else {
                                DirectForceFunction::direct_translational(tempLattice, particleSystem->params.domainSize,
                                        rangevalues[threadNum], particleSystem->chargedParticles);
                            }
                        }
                    }
                }
            }
            else {
//                numBoxes++;
//                tempLattice.print();
                
//                domainOffset.x = tempLattice.point.x * particleSystem->params.domainSize.x;
//                domainOffset.y = tempLattice.point.y * particleSystem->params.domainSize.y;
//                domainOffset.z = tempLattice.point.z * particleSystem->params.domainSize.z;
                
                if (vectorizedToggle){
                    DirectForceFunction::direct_vectorized_translational(tempLattice, particleSystem->params.domainSize,
                            rangevalues[threadNum], particleSystem->chargedParticles, jptc);
                }
                else {
                    DirectForceFunction::direct_translational(tempLattice, particleSystem->params.domainSize,
                            rangevalues[threadNum], particleSystem->chargedParticles);
                }
            }
        }
    }
    
//    printf("\nnumBoxes = %d\n",numBoxes);
}

/**
 * Computes the force on all particles using the direct force method for a
 *      translationally symmetric system
 * @param domainOffset - offset denoting location of translated cell relative to real cell
 * @param rangevalues - bounds for indices on particles to be included in computing force
 * @param particles - reference to ChargedParticles class being considered
 */
//void DirectForceFunction::direct_translational(volatile const vec3<float> domainOffset,
//        vec4<int> rangevalues, ChargedParticles& particles) {
void DirectForceFunction::direct_translational(LatticePoint lattice, vec3<float> domainSize,
        vec4<int> rangevalues, ChargedParticles& particles){
    
    vec3<float> domainOffset;
    domainOffset.x = lattice.point.x * domainSize.x;
    domainOffset.y = lattice.point.y * domainSize.y;
    domainOffset.z = lattice.point.z * domainSize.z;
    
//    printf("domainOffset.x = %.9g, .y = %.9g, .z = %.9g\n", domainOffset.x, domainOffset.y, domainOffset.z);
    
    vec3<double> dist, ai;
    double invDist, invDistCube, accel_const;
    int i, j;
    
    const double eps = PhysicalConstant::SOFTENING;
    const double K_const = PhysicalConstant::K_CONST;
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    const float affectsAccelerationCheck = static_cast<float> (ParticleStatus::affectsAcceleration);

    for (i = rangevalues.x; i < rangevalues.y; i++){
        if (particles.acceleration[i].w == inboundCheck) {
            ai = 0;
            for (j = rangevalues.z; j < rangevalues.w; j++){
                if (particles.acceleration[j].w == inboundCheck ||
                        particles.acceleration[j].w == affectsAccelerationCheck) {
                    dist.x = (particles.position[i].x - domainOffset.x) - particles.position[j].x;
                    dist.y = (particles.position[i].y - domainOffset.y) - particles.position[j].y;
                    dist.z = (particles.position[i].z - domainOffset.z) - particles.position[j].z;
                    invDist = 1.0 / sqrt(dist.x * dist.x + dist.y * dist.y + dist.z * dist.z + eps);
                    invDistCube = particles.position[j].w * invDist * invDist * invDist;
                    ai.x += dist.x * invDistCube;
                    ai.y += dist.y * invDistCube;
                    ai.z += dist.z * invDistCube;
//                    printf("domainOffset.x = %.9g, .y = %.9g, .z = %.9g\n",domainOffset.x, domainOffset.y, domainOffset.z);
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
 *      translationally symmetric system using vectorized functions
 * @param domainOffset - offset denoting location of translated cell relative to real cell
 * @param rangevalues - bounds for indices on particles to be included in computing force
 * @param particles - reference to ChargedParticles class being considered
 * @param jptc - pointer to Jpdata array - an aligned array of particle position data
 */
//void DirectForceFunction::direct_vectorized_translational(volatile const vec3<float> domainOffset,
//        vec4<int> rangevalues, ChargedParticles& particles, Jpdata* jptc){
void DirectForceFunction::direct_vectorized_translational(LatticePoint lattice, vec3<float> domainSize,
        vec4<int> rangevalues, ChargedParticles& particles, Jpdata* jptc){
    int i, j, ii;
    int numParticles = rangevalues.w - rangevalues.z;
    
    Ipdata iptc;
    Fodata fout;
    for (j = 0; j < 4; j++){
        fout.eps2[j] = PhysicalConstant::SOFTENING;
    }
    
    volatile const float K_const = PhysicalConstant::K_CONST;
    volatile const float inboundCheck = static_cast<float>(ParticleStatus::inbound);
//    volatile const vec3<float> domainOffset = domainOffsetTemp;
    volatile const vec3<float> domainOffset = {lattice.point.x * domainSize.x,
            lattice.point.y * domainSize.y, lattice.point.z * domainSize.z};
    
//    printf("domainOffset.x = %.9g, .y = %.9g, .z = %.9g\n", domainOffset.x, domainOffset.y, domainOffset.z);
    
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
            
        for (ii = 0; ii < 4; ii++){
            if (ii < j){
                iptc.x[ii] = particles.position[temp[ii]].x - domainOffset.x;
                iptc.y[ii] = particles.position[temp[ii]].y - domainOffset.y;
                iptc.z[ii] = particles.position[temp[ii]].z - domainOffset.z;
                iptc.q[ii] = K_const * particles.position[temp[ii]].w / particles.getMass(temp[ii]);
//                printf("iptc.x[%d] = %.9g, .y = %.9g, .z = %.9g, .q = %.9g\n",ii,
//                        iptc.x[ii],iptc.y[ii],iptc.z[ii],iptc.q[ii]);
//                printf("domainOffset.x = %.9g, .y = %.9g, .z = %.9g\n",domainOffset.x, domainOffset.y, domainOffset.z);
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
            particles.acceleration[temp[ii]].x += fout.ax[ii];
            particles.acceleration[temp[ii]].y += fout.ay[ii];
            particles.acceleration[temp[ii]].z += fout.az[ii];
//            particles.printAccelInfo(temp[ii]);
        }
    }
}
