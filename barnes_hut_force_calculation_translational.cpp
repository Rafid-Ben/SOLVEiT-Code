/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "force_calculation_methods.h"
#include "direct_force_method.h"

/**
 * Multi-threaded function to call the Barnes-Hut force computation method for a
 *      translationally symmetric system
 * @param dF - struct that holds a pointer to the class and the thread number
 * @return nullptr
 */
void* BarnesHutForceCalculation::translationalBarnesHutKernel(void* dF){
    barnesHutForceStruct dfc = *((barnesHutForceStruct*)dF);
    dfc.d->translationalBarnesHutKernelCalculation(dfc.threadNum);
    return nullptr;
}

/**
 * Computes the force on all particles using the Barnes-Hut approximation for a
 *      translationally symmetric system
 * @param threadNum - thread number for multi-threaded versions
 */
void BarnesHutForceCalculation::translationalBarnesHutKernelCalculation(int threadNum){
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
                            
                            translationalBarnesHut(tempLattice, threadNum);
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
                                
                                translationalBarnesHut(tempLattice, threadNum);
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
                                
                                translationalBarnesHut(tempLattice, threadNum);
                            }
                        }
                    }
                }
            }
            else {
//                numBoxes++;
//                tempLattice.print();
                
                translationalBarnesHut(tempLattice, threadNum);
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
                            
                            translationalBarnesHut(tempLattice, threadNum);
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
                                
                                translationalBarnesHut(tempLattice, threadNum);
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
                
                translationalBarnesHut(tempLattice, threadNum);
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
                            
                            translationalBarnesHut(tempLattice, threadNum);
                        }
                    }
                }
            }
            else {
//                numBoxes++;
//                tempLattice.print();
                
                translationalBarnesHut(tempLattice, threadNum);
            }
        }
    }
    
//    printf("\nnumBoxes = %d\n",numBoxes);
}

void BarnesHutForceCalculation::translationalBarnesHut(LatticePoint lattice, int threadNum){
    const float inboundCheck = static_cast<float>(ParticleStatus::inbound);
    vec3<float> accel, tempPos;
    double accel_const;
    int i;
    
    vec3<float> domainOffset;
    domainOffset.x = lattice.point.x * particleSystem->params.domainSize.x;
    domainOffset.y = lattice.point.y * particleSystem->params.domainSize.y;
    domainOffset.z = lattice.point.z * particleSystem->params.domainSize.z;
    
    int centerCheck = false;
    if (lattice.point.x == 0 && lattice.point.y == 0 && lattice.point.z == 0){
        centerCheck = true;
    }
    
    for (i = rangevalues[threadNum].x; i < rangevalues[threadNum].y; i++){
        if (particleSystem->chargedParticles.acceleration[i].w == inboundCheck){
            tempPos.x = particleSystem->chargedParticles.position[i].x - domainOffset.x;
            tempPos.y = particleSystem->chargedParticles.position[i].y - domainOffset.y;
            tempPos.z = particleSystem->chargedParticles.position[i].z - domainOffset.z;
            
            accel = 0.0f;
            
            computeForceFromTreeTranslations(rootCell, tempPos, i, accel, centerCheck);
            
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


inline void BarnesHutForceCalculation::computeForceFromCellTranslations(Cell* cell,
        vec3<float> pos, vec3<float>& accel, float deff){
//    float deff = computeEffectiveDistance(index, cell);
    
    // Compute electric force according to Coulomb's law
    float f = cell->charge / (deff * deff * deff);
    
    accel.x += f * (pos.x - cell->centerOfCharge.x);
    accel.y += f * (pos.y - cell->centerOfCharge.y);
    accel.z += f * (pos.z - cell->centerOfCharge.z);
}

void BarnesHutForceCalculation::computeForceFromTreeTranslations(Cell* cell, 
        vec3<float> pos, int index, vec3<float>& accel, int centerCheck){
        
    Cell* temp = cell;
    float d;
    while (temp != NULL){
        if (temp->firstChild == NULL){
            if (temp->index != -1 && (temp->index != index || centerCheck != true)){
                d = computeEffectiveDistance(pos, temp);
//                if (index == 28){
//                    printf("Only one particle: %p, dist: %.9g, charge: %.9g\n",(void*)temp, d,temp->charge);
//                }
                computeForceFromCellTranslations(temp, pos, accel, d);
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
                computeForceFromCellTranslations(temp, pos, accel, d);
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