/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "force_calculation_methods.h"
#include "direct_force_method.h"


/**
 * Multi-threaded function to call the Barnes-Hut force computation method for a
 *      rotationally symmetric system
 * @param dF - struct that holds a pointer to the class and the thread number
 * @return nullptr
 */
void* BarnesHutForceCalculation::rotationalBarnesHutKernel(void* dF){
    barnesHutForceStruct dfc = *((barnesHutForceStruct*)dF);
    dfc.d->rotationalBarnesHutKernelCalculation(dfc.threadNum);
    return nullptr;
}

/**
 * Computes the force on all particles using the Barnes-Hut approximation for a
 *      rotationally symmetric system
 * @param threadNum - thread number for multi-threaded versions
 */
void BarnesHutForceCalculation::rotationalBarnesHutKernelCalculation(int threadNum){
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    vec3<float> accel, tempPos;
    double accel_const;
    int i, k, cf, sf, cb, sb;
    
    for (k = 0; k < 4; k++){
        cf = quadrantRotationFunction(5 - k);
        sf = quadrantRotationFunction(4 - k);
        cb = quadrantRotationFunction(k + 1);
        sb = quadrantRotationFunction(k);
        
        for (i = rangevalues[threadNum].x; i < rangevalues[threadNum].y; i++){
//        if (rangevalues[threadNum].y == 485 && i == 159){
//            printf("test\n");
//        }
//        if (i == 28){
//            printf("test\n");
//        }
            if (particleSystem->chargedParticles.acceleration[i].w == inboundCheck){

                tempPos.x = cf * particleSystem->chargedParticles.position[i].x -
                        sf * particleSystem->chargedParticles.position[i].y;
                tempPos.y = sf * particleSystem->chargedParticles.position[i].x + 
                        cf * particleSystem->chargedParticles.position[i].y;
                tempPos.z = particleSystem->chargedParticles.position[i].z;

                accel = 0.0f;
                computeForceFromTreeRotations(rootCell, tempPos, i, accel, cb, sb, k);

                accel_const = PhysicalConstant::K_CONST * particleSystem->chargedParticles.position[i].w / 
                        particleSystem->chargedParticles.getMass(i);
                accel *= accel_const;

                particleSystem->chargedParticles.acceleration[i].x += accel.x;
                particleSystem->chargedParticles.acceleration[i].y += accel.y;
                particleSystem->chargedParticles.acceleration[i].z += accel.z;
                
//                particleSystem->chargedParticles.printAccelInfo(i);
            }
            else {
                particleSystem->chargedParticles.acceleration[i].x += 0;
                particleSystem->chargedParticles.acceleration[i].y += 0;
                particleSystem->chargedParticles.acceleration[i].z += 0;
            }
        }
    }
}


float BarnesHutForceCalculation::computeEffectiveDistance(vec3<float> pos, Cell* cell){
    float distx = pos.x - cell->centerOfCharge.x;
    float disty = pos.y - cell->centerOfCharge.y;
    float distz = pos.z - cell->centerOfCharge.z;
    
    return sqrt(distx * distx + disty * disty + distz * distz + PhysicalConstant::SOFTENING);
}


inline void BarnesHutForceCalculation::computeForceFromCellRotations(Cell* cell,
        vec3<float> pos, vec3<float>& accel, float deff, int cb, int sb){
//    float deff = computeEffectiveDistance(index, cell);
    
    // Compute electric force according to Coulomb's law
    float f = cell->charge / (deff * deff * deff);
    
    float rx = pos.x - cell->centerOfCharge.x;
    float ry = pos.y - cell->centerOfCharge.y;
    
    accel.x += f * (cb * rx - sb * ry);
    accel.y += f * (sb * rx + cb * ry);
    accel.z += f * (pos.z - cell->centerOfCharge.z);
}

void BarnesHutForceCalculation::computeForceFromTreeRotations(Cell* cell, 
        vec3<float> pos, int index, vec3<float>& accel, int cb, int sb, int quadrant){
        
    Cell* temp = cell;
    float d;
    while (temp != NULL){
        if (temp->firstChild == NULL){
            if (temp->index != -1 && (temp->index != index || quadrant != 0)){
                d = computeEffectiveDistance(pos, temp);
//                if (index == 28){
//                    printf("Only one particle: %p, dist: %.9g, charge: %.9g\n",(void*)temp, d,temp->charge);
//                }
                computeForceFromCellRotations(temp, pos, accel, d, cb, sb);
            }
            temp = temp->nextSibling;
        }
        else {
            d = computeEffectiveDistance(pos, temp);
//            if (d > ((temp->repLength / theta) + temp->delta)){
            if ((theta * (d - temp->delta)) > temp->repLength){
//                if (index == 28){
//                    printf("Far enough away: %p, dist: %.9g, charge: %.9g\n",(void*)temp, d,temp->charge);
//                }
                computeForceFromCellRotations(temp, pos, accel, d, cb, sb);
                temp = temp->nextSibling;
            }
            else {
//                if (index == 28){
//                    printf("Too close: %p, dist: %.9g, charge: %.9g\n",(void*)temp, d,temp->charge);
//                }
                temp = temp->firstChild;
            }
        }
    }
}