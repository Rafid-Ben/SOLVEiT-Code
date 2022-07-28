/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "force_calculation_methods.h"


void* MultipoleMethodForceCalculation::direct_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->direct_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::direct_kernel_translations(int threadNum){
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

void* MultipoleMethodForceCalculation::p2p_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->p2p_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::p2p_kernel_translations(int threadNum){
    int ii, ij, jj;
    int tempDomainIndex;
    
    LatticePoint tempLattice;
    
    vec4<int> r = rangevalues[threadNum];
    
    vec4<int> particleBounds;
    
    for (ii = r.x; ii < r.y; ii++){
        particleBounds.x = variables->particleOffset[0][ii];
        particleBounds.y = variables->particleOffset[1][ii] + 1;
        for (ij = 0; ij < variables->numInteraction[ii]; ij++){
            jj = variables->interactionList[ii][ij][0];
            tempDomainIndex = variables->interactionList[ii][ij][1];
            MultipoleMethodFunction::undoBalancedTernary(tempDomainIndex, tempLattice.point);
            
            particleBounds.z = variables->particleOffset[0][jj];
            particleBounds.w = variables->particleOffset[1][jj] + 1;
            
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

void* MultipoleMethodForceCalculation::p2m_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->p2m_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::p2m_kernel_translations(int threadNum){
    int jj, j, n, m, nm, nms;
    vec3<int> boxIndex3D;
    vec3<int> domainBoxIndex3D;
    vec3<float> boxCenter, domainBoxMin;
    vec3<double> dist;
    double boxSize, rho, alpha, beta;
    double xx, s2, fact, pn, p, p1, p2;
    double rhom, rhon;
    double YnmReal[FMMConstant::numExpansion2];
    std::complex<double> MnmVector[numCoefficients], I(0.0, 1.0), eim;
        
    boxSize = rootBoxSize / pow(3.0, maxLevel);
    
    vec4<int> r = rangevalues[threadNum];
    
    int domainNum, domainMorton3;
    
    for (jj = r.x; jj < r.y; jj++){
        // Find the domain that the jj-th box lies within
        if (jj == r.x) domainNum = findDomain(jj, maxLevel);
        else{
            // If jj is in a different domain to jj - 1, find it's new domain
            if (jj > variables->domainOffset[1][domainNum][maxLevel]){
                // Increment the domain number until find a domain containing boxes with particles inside
                // Since jj will only increase by 1 between iterations, the next domain containing particles contains jj
                do{
                    domainNum++;
                }
                while (variables->domainOffset[1][domainNum][maxLevel] < variables->domainOffset[0][domainNum][maxLevel]);
            }
        }
        
        // Find the global location of the domain
        domainMorton3 = variables->domainLinkage[0][domainNum].mortonIndex;
        MultipoleMethodFunction::unmorton3(domainMorton3, domainBoxIndex3D);
        domainBoxIndex3D.x -= variables->offsetSC[0].x;
        domainBoxIndex3D.y -= variables->offsetSC[0].y;
        domainBoxIndex3D.z -= variables->offsetSC[0].z;
        
        // Find the local location of the jj-th cell within its domain
        MultipoleMethodFunction::unmorton3(variables->domainBoxIndexFull[jj], boxIndex3D);
        
        // Find the coordinates of the box center
        domainBoxMin.x = boxMin.x + domainBoxIndex3D.x * rootBoxSize;
        domainBoxMin.y = boxMin.y + domainBoxIndex3D.y * rootBoxSize;
        domainBoxMin.z = boxMin.z + domainBoxIndex3D.z * rootBoxSize;
        boxCenter.x = domainBoxMin.x + (boxIndex3D.x + 0.5) * boxSize;
        boxCenter.y = domainBoxMin.y + (boxIndex3D.y + 0.5) * boxSize;
        boxCenter.z = domainBoxMin.z + (boxIndex3D.z + 0.5) * boxSize;
        
        // Initialize MnmVector
        for (j = 0; j < numCoefficients; j++) {
            MnmVector[j] = 0;
        }
        
        // Compute MnmVector from all particles inside the jj-th cell
        for (j = variables->particleOffset[0][jj]; j <= variables->particleOffset[1][jj]; j++) {
            dist.x = particleSystem->chargedParticles.position[j].x - boxCenter.x;
            dist.y = particleSystem->chargedParticles.position[j].y - boxCenter.y;
            dist.z = particleSystem->chargedParticles.position[j].z - boxCenter.z;
            cart2sph(rho, alpha, beta, dist.x, dist.y, dist.z);
            rho /= boxSize;
            
            xx = cos(alpha);
            s2 = sqrt((1 - xx) * (1 + xx));
            fact = 1;
            pn = 1;
            rhom = 1.0;
            for (m = 0; m < numExpansions; m++) {
                p = pn;
                nm = m * m + 2 * m;
                YnmReal[nm] = rhom * constants->factorial[nm] * p;
                p1 = p;
                p = xx * (2 * m + 1) * p;
                rhom *= rho;
                rhon = rhom;
                for (n = m + 1; n < numExpansions; n++) {
                    nm = n * n + n + m;
                    YnmReal[nm] = rhon * constants->factorial[nm] * p;
                    p2 = p1;
                    p1 = p;
                    p = (xx * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
                    rhon *= rho;
                }
                pn = -pn * fact * s2;
                fact += 2;
            }
            for (n = 0; n < numExpansions; n++) {
                for (m = 0; m <= n; m++) {
                    nm = n * n + n + m;
                    nms = n * (n + 1) / 2 + m;
                    eim = exp(-m * beta * I);
                    MnmVector[nms] += ((std::complex<double>) particleSystem->chargedParticles.position[j].w) * YnmReal[nm] * eim;
                }
            }
        }
        
        // Assign Mnm
        for (j = 0; j < numCoefficients; j++) {
            variables->Mnm[jj][j] = MnmVector[j];
        }
    }
}

void* MultipoleMethodForceCalculation::m2m_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->m2m_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::m2m_kernel_translations(int threadNum){
    int ib, j, jj, nfjp, nfjc, jb, k, jk, jks, n, jnk, jnks, nm;
    int je_forward, je_backward;
    vec3<int> boxIndex3D;
    int boxIndexFullTemp;
    double rho, powThree;
    std::complex<double> cnm, MnmScalar;
    std::complex<double> MnmVectorB[numCoefficients], MnmVectorA[numCoefficients];
        
    vec4<int> r = rangevalues[threadNum];
    int numLevel = r.w;
    
    int domainNum, domainBias;
//    int domainNumBoxIndexTemp = 1 << (3 * numLevel);
    int domainNumBoxIndexTemp = pow(3, 3 * numLevel);
    
    int tempLevelOffset = (numLevel == 0) ? variables->levelOffset[maxLevel] : variables->levelOffset[numLevel - 1];
    
    int tempQuad, m;
    std::complex<double> mr;
    
    for (jj = r.x; jj < r.y; jj++) {
        jb = jj + variables->levelOffset[numLevel];
        boxIndexFullTemp = variables->domainBoxIndexFull[jb];
        
        // Find global domain that the jj-th cell lies within
        if (jj == r.x){
            domainNum = findDomain(jj, numLevel + 1);
            domainBias = domainNum * domainNumBoxIndexTemp;
        }
        else{
            // If jj is in a different domain to jj - 1, find it's new domain
            if (jj > variables->domainOffset[1][domainNum][numLevel + 1]){
                // Increment the domain number until find a domain containing boxes with particles inside
                // Since jj will only increase by 1 between iterations, the next domain containing particles contains jj
                do{
                    domainNum++;
                }
                while (variables->domainOffset[1][domainNum][numLevel + 1] < variables->domainOffset[0][domainNum][numLevel + 1]);
                domainBias = domainNum * domainNumBoxIndexTemp;
            }
        }
        
        // nfjp encodes information about the parent cell
        // nfjc encodes which of the parent's children is being considered
        nfjp = boxIndexFullTemp / 27;
        nfjc = boxIndexFullTemp % 27;
        if (numLevel == 0){
            ib = variables->domainBoxIndexMask[domainNum] + tempLevelOffset;
        }
        else{
            ib = variables->domainBoxIndexMask[domainBias + nfjp] + tempLevelOffset;
        }

        MultipoleMethodFunction::unmorton3(nfjc, boxIndex3D);
        
        rho = sqrt(abs(boxIndex3D.x - 1) + abs(boxIndex3D.y - 1) + abs(boxIndex3D.z - 1));
        
        boxIndex3D.x = 6 - boxIndex3D.x;
        boxIndex3D.y = 6 - boxIndex3D.y;
        boxIndex3D.z = 6 - boxIndex3D.z;

        rotationMorton(boxIndex3D, je_forward, je_backward);        

        // Initialize MnmVectorA with Mnm values from lower level cell
        for (j = 0; j < numCoefficients; j++) {
            MnmVectorA[j] = variables->Mnm[jb][j];
        }

        rotation(MnmVectorA, MnmVectorB, constants->Dnm[je_forward]);

        // *** Note that in the calculation of cnm, use pow(3.0, n-j) vice pow(2.0, n-j) since considering 3x3x3 boxes
        for (j = 0; j < numExpansions; j++) {
            powThree = pow(3.0, -j);
            for (k = 0; k <= j; k++) {
                jk = j * j + j + k;
                jks = j * (j + 1) / 2 + k;
                MnmScalar = 0;
                for (n = 0; n <= j - abs(k); n++) {
                    jnk = (j - n) * (j - n) + j - n + k;
                    jnks = (j - n) * (j - n + 1) / 2 + k;
                    nm = n * n + n;
//                    cnm = pow(-1.0, n) * anm[nm] * anm[jnk] / anm[jk] * pow(rho, n) * Ynm[nm] * pow(2.0, n-j);
                    cnm = pow(-1.0, n) * constants->anm[nm] * constants->anm[jnk] / constants->anm[jk] * pow(rho, n) * constants->Ynm[nm] * powThree;
                    MnmScalar += MnmVectorB[jnks] * cnm;
                }
                MnmVectorA[jks] = MnmScalar;
            }
        }

        rotation(MnmVectorA, MnmVectorB, constants->Dnm[je_backward]);

        // Ensure not writing to the same Mnm location
#ifdef _PTHREAD_H
        pthread_mutex_lock(&mutex_fmm);
#endif
        for (j = 0; j < numCoefficients; j++) {
            variables->Mnm[ib][j] += MnmVectorB[j];
            //printf("variables->Mnm[%d][%d] = %.9g + i%.9g\n",ib,j,std::real(variables->Mnm[ib][j]),std::imag(variables->Mnm[ib][j]));
        }
#ifdef _PTHREAD_H
        pthread_mutex_unlock(&mutex_fmm);
#endif
    }
}

void* MultipoleMethodForceCalculation::m2l_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->m2l_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::m2l_kernel_translations(int threadNum){
    int j, ii, ib, ix, iy, iz, ij, jj, jb, jx, jy, jz, k, jk, jks, n, nk, nks, jkn, jnk;
    vec3<int> boxIndex3D;
    vec3<double> dist;
    double boxSize, rho, rhoj, rhojk, rhojn;
    std::complex<double> LnmVectorA[numCoefficients], MnmVectorA[numCoefficients];
    std::complex<double> LnmVectorB[numCoefficients], MnmVectorB[numCoefficients];
    std::complex<double> cnm, LnmScalar;
    
    vec3<int> domainBoxIndex3D;
    
    vec4<int> r = rangevalues[threadNum];
    int numLevel = r.w;
    int jdiff = pow(3, numLevel);
    int domainNum_i, domainNum_j;
    int je_forward, je_backward;
    int gcd;
    std::complex<double> mr;
    int tempDomainIndex, m;
    vec3<int> domainSymOffset;
    
    boxSize = rootBoxSize / pow(3.0, numLevel);
    double boxSize2 = boxSize * boxSize;
    
    int tempLevelOffset = (numLevel == 0) ? variables->levelOffset[maxLevel] : variables->levelOffset[numLevel - 1];
    
    for (ii = r.x; ii < r.y; ii++) {
        ib = ii + tempLevelOffset;
        
        // Find global domain of ii-th cell
        if (ii == r.x){
            domainNum_i = findDomain(ii, numLevel);
        }
        else{
            // If ii is in a different domain to ii - 1, find it's new domain
            if (ii > variables->domainOffset[1][domainNum_i][numLevel]){
                // Increment the domain number until find a domain containing boxes with particles inside
                // Since ii will only increase by 1 between iterations, the next domain containing particles contains ii
                do{
                    domainNum_i++;
                }
                while (variables->domainOffset[1][domainNum_i][numLevel] < variables->domainOffset[0][domainNum_i][numLevel]);
            }
        }
        
        // Find global domain 3D indices
        MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][domainNum_i].mortonIndex, domainBoxIndex3D);
        domainBoxIndex3D.x -= variables->offsetSC[0].x;
        domainBoxIndex3D.y -= variables->offsetSC[0].y;
        domainBoxIndex3D.z -= variables->offsetSC[0].z;
        
        // Find local domain 3D indices
        MultipoleMethodFunction::unmorton3(variables->domainBoxIndexFull[ib], boxIndex3D);
        
        ix = boxIndex3D.x + domainBoxIndex3D.x * jdiff;
        iy = boxIndex3D.y + domainBoxIndex3D.y * jdiff;
        iz = boxIndex3D.z + domainBoxIndex3D.z * jdiff;
        
        for (ij = 0; ij < variables->numInteraction[ii]; ij++) {
            
            jj = variables->interactionList[ii][ij][0];
            tempDomainIndex = variables->interactionList[ii][ij][1];
            
            MultipoleMethodFunction::undoBalancedTernary(tempDomainIndex, domainSymOffset);
            
            jb = jj + tempLevelOffset;
            
            // Initialize MnmVectorB
            for (j = 0; j < numCoefficients; j++){
                MnmVectorB[j] = variables->Mnm[jb][j];
            }
            
            // Find global domain of jj-th cell
            domainNum_j = findDomain(jj, numLevel);
            
            // Find global domain 3D indices
            MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][domainNum_j].mortonIndex, domainBoxIndex3D);
            domainBoxIndex3D.x -= variables->offsetSC[0].x;
            domainBoxIndex3D.y -= variables->offsetSC[0].y;
            domainBoxIndex3D.z -= variables->offsetSC[0].z;
            
            domainBoxIndex3D.x += domainSymOffset.x * boxDimVals[0];
            domainBoxIndex3D.y += domainSymOffset.y * boxDimVals[1];
            domainBoxIndex3D.z += domainSymOffset.z * boxDimVals[2];
            
            // Find local domain 3D indices
            MultipoleMethodFunction::unmorton3(variables->domainBoxIndexFull[jb], boxIndex3D);
            
            jx = boxIndex3D.x + domainBoxIndex3D.x * jdiff;
            jy = boxIndex3D.y + domainBoxIndex3D.y * jdiff;            
            jz = boxIndex3D.z + domainBoxIndex3D.z * jdiff;
            
            boxIndex3D.x = ix - jx;
            boxIndex3D.y = iy - jy;
            boxIndex3D.z = iz - jz;
            
            dist.x = boxIndex3D.x;
            dist.y = boxIndex3D.y;
            dist.z = boxIndex3D.z;

            rho = sqrt(dist.x * dist.x + dist.y * dist.y + dist.z * dist.z + PhysicalConstant::SOFTENING / boxSize2);
                        
            // Remove gcd for rotation purposes
            // gcd should always be between 0 and 5 for the case of 3x3x3 boxes
            gcd = gcd3(boxIndex3D.x, boxIndex3D.y, boxIndex3D.z);
            switch (gcd){
                case 2:
                    boxIndex3D.x >>= 1;
                    boxIndex3D.y >>= 1;
                    boxIndex3D.z >>= 1;
                    break;
                case 3:
                    boxIndex3D.x = SIGNUM(boxIndex3D.x);
                    boxIndex3D.y = SIGNUM(boxIndex3D.y);
                    boxIndex3D.z = SIGNUM(boxIndex3D.z);
                    break;
                case 4:
                    boxIndex3D.x >>= 2;
                    boxIndex3D.y >>= 2;
                    boxIndex3D.z >>= 2;
                    break;
                case 5:
                    boxIndex3D.x = SIGNUM(boxIndex3D.x);
                    boxIndex3D.y = SIGNUM(boxIndex3D.y);
                    boxIndex3D.z = SIGNUM(boxIndex3D.z);
                    break;
                default:
                    if (gcd != 0){
                        // Should never happen
                        boxIndex3D.x /= gcd;
                        boxIndex3D.y /= gcd;
                        boxIndex3D.z /= gcd;
                    }
                    break;
            }
            
            boxIndex3D.x += 5;
            boxIndex3D.y += 5;
            boxIndex3D.z += 5;
            
            rotationMorton(boxIndex3D, je_forward, je_backward);
            
            rotation(MnmVectorB, MnmVectorA, constants->Dnm[je_forward]);

            rhoj = 1;
            for (j = 0; j < numExpansions; j++) {
                rhojk = rhoj;
                rhoj *= rho;
                for (k = 0; k <= j; k++) {
                    jk = j * j + j + k;
                    jks = j * (j + 1) / 2 + k;
                    LnmScalar = 0;
                    rhojn = rhojk;
                    rhojk *= rho;
                    for (n = abs(k); n < numExpansions; n++) {
                        rhojn *= rho;
                        nk = n * n + n + k;
                        nks = n * (n + 1) / 2 + k;
                        jkn = jk * FMMConstant::numExpansion2 + nk;
                        jnk = (j + n) * (j + n) + j + n;
                        cnm = constants->Anm[jkn] / rhojn * constants->Ynm[jnk];
                        LnmScalar += MnmVectorA[nks] * cnm;
                    }
                    LnmVectorA[jks] = LnmScalar;
                }
            }

            rotation(LnmVectorA, LnmVectorB, constants->Dnm[je_backward]);
            
            // Assign Lnm
            for (j = 0; j < numExpansions; j++){
                for (n = 0; n <= j; n++){
                    variables->Lnm[ii][j*(j+1)/2 + n] += LnmVectorB[j*(j+1)/2 + n];
//                    if (sqrt(std::norm(variables->Lnm[ii][j*(j+1)/2 + n])) > 1000){
//                        printf("%f + %fi\n",variables->Lnm[ii][j*(j+1)/2 + n].real(),variables->Lnm[ii][j*(j+1)/2+n].imag());
//                        printf("test\n");
//                    }
                }
            }            
        }
    }
}

void* MultipoleMethodForceCalculation::l2l_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->l2l_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::l2l_kernel_translations(int threadNum){
    int ii, ib, i, nfip, nfic, j, k, jk, jks, n, jnk, nk, nks;
    vec3<int> boxIndex3D;
    double rho;
    std::complex<double> cnm, LnmScalar;
    std::complex<double> LnmVectorA[numCoefficients], LnmVectorB[numCoefficients];
    
    vec4<int> r = rangevalues[threadNum];
    int numLevel = r.w;

    int tempLevelOffset = (numLevel == 0) ? variables->levelOffset[maxLevel] : variables->levelOffset[numLevel - 1];
    
    int domainNum, domainBias;
    int domainNumBoxIndexTemp = pow(3, 3 * numLevel);
    
    int je_forward, je_backward;
    
    for (ii = r.x; ii < r.y; ii++) {        
        ib = ii + tempLevelOffset;
        
        nfip = variables->domainBoxIndexFull[ib] / 27;
        nfic = variables->domainBoxIndexFull[ib] % 27;
        MultipoleMethodFunction::unmorton3(nfic, boxIndex3D);
        
        // If ii is in a different domain to ii - 1, find it's new domain
        if (ii == r.x){
            domainNum = findDomain(ii, numLevel);
            domainBias = domainNum * domainNumBoxIndexTemp;
        }
        else if (ii > variables->domainOffset[1][domainNum][numLevel]){
            // Increment the domain number until find a domain containing boxes with particles inside
            // Since ii will only increase by 1 between iterations, the next domain containing particles contains ii
            do{
                domainNum++;
            }
            while (variables->domainOffset[1][domainNum][numLevel] < variables->domainOffset[0][domainNum][numLevel]);
            domainBias = domainNum * domainNumBoxIndexTemp;
        }
        
        // rho is distance from center of child to center of parent cell
        // rho is normalized, such that child cell has length equal to one
        rho = sqrt(abs(boxIndex3D.x - 1) + abs(boxIndex3D.y - 1) + abs(boxIndex3D.z - 1));
        
        boxIndex3D.x = boxIndex3D.x + 4;
        boxIndex3D.y = boxIndex3D.y + 4;
        boxIndex3D.z = boxIndex3D.z + 4;
        
        rotationMorton(boxIndex3D, je_forward, je_backward);
        
        ib = variables->neo[nfip + domainBias];
        
        for (i = 0; i < numCoefficients; i++) {
            LnmVectorA[i] = variables->LnmOld[ib][i];
        }

        rotation(LnmVectorA, LnmVectorB, constants->Dnm[je_forward]);

        for (j = 0; j < numExpansions; j++){
            for (k = 0; k <= j; k++){
                jk = j * j + j + k;
                jks = j*(j+1)/2+k;
                LnmScalar = 0;
                for (n = j; n < numExpansions; n++){
                    jnk = (n-j)*(n-j)+n-j;
                    nk = n*n+n+k;
                    nks = n*(n+1)/2+k;
                    cnm = constants->anm[jnk] * constants->anm[jk] / constants->anm[nk] * pow(rho, n-j) * constants->Ynm[jnk] * pow(3.0, -n-1);
                    LnmScalar += LnmVectorB[nks] * cnm;
                }
                LnmVectorA[jks] = LnmScalar;
            }
        }
        
        rotation(LnmVectorA, LnmVectorB, constants->Dnm[je_backward]);
        for (i = 0; i < numCoefficients; i++) {
            variables->Lnm[ii][i] = LnmVectorB[i];
            //printf("variables->Lnm[%d][%d] = %.9g + i%.9g\n",ii,i,std::real(variables->Lnm[ii][i]),std::imag(variables->Lnm[ii][i]));
        }
    }
}

void* MultipoleMethodForceCalculation::l2p_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->l2p_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::l2p_kernel_translations(int threadNum){
    vec4<int> range = rangevalues[threadNum];
    
    int ii, i, n, nm, nms, m;
    vec3<int> boxIndex3D;
    vec3<int> domainBoxIndex3D;
    vec3<float> boxCenter, domainBoxMin;
    vec3<double> accel, dist;
    double boxSize, r, theta, phi, accelR, accelTheta, accelPhi, accel_const;
    double xx, yy, s2, fact, pn, p, p1, p2, rn;
    double YnmReal[FMMConstant::numExpansion2], YnmRealTheta[FMMConstant::numExpansion2];
    std::complex<double> LnmVector[numCoefficients];
    std::complex<double> rr, rtheta, rphi, I(0.0, 1.0), eim;
    int domainNum, domainMorton3;
    
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    
    boxSize = rootBoxSize / pow(3.0, maxLevel);
    for (ii = range.x; ii < range.y; ii++) {
        
        // Find global domain that ii-th cell lies within
        if (ii == range.x){
            domainNum = findDomain(ii, maxLevel);
        }
        else{
            // If ii is in a different domain to ii - 1, find it's new domain
            if (ii > variables->domainOffset[1][domainNum][maxLevel]){
                // Increment the domain number until find a domain containing boxes with particles inside
                // Since ii will only increase by 1 between iterations, the next domain containing particles contains ii
                do{
                    domainNum++;
                }
                while (variables->domainOffset[1][domainNum][maxLevel] < variables->domainOffset[0][domainNum][maxLevel]);
            }
        }
        
        // Find global and local 3D box indices
        domainMorton3 = variables->domainLinkage[0][domainNum].mortonIndex;
        MultipoleMethodFunction::unmorton3(domainMorton3, domainBoxIndex3D);
        domainBoxIndex3D.x -= variables->offsetSC[0].x;
        domainBoxIndex3D.y -= variables->offsetSC[0].y;
        domainBoxIndex3D.z -= variables->offsetSC[0].z;
        
        MultipoleMethodFunction::unmorton3(variables->domainBoxIndexFull[ii], boxIndex3D);
                
        domainBoxMin.x = boxMin.x + domainBoxIndex3D.x * rootBoxSize;
        domainBoxMin.y = boxMin.y + domainBoxIndex3D.y * rootBoxSize;
        domainBoxMin.z = boxMin.z + domainBoxIndex3D.z * rootBoxSize;
        boxCenter.x = domainBoxMin.x + (boxIndex3D.x + 0.5) * boxSize;
        boxCenter.y = domainBoxMin.y + (boxIndex3D.y + 0.5) * boxSize;
        boxCenter.z = domainBoxMin.z + (boxIndex3D.z + 0.5) * boxSize;
        
        // Initialize Lnm
        for (i = 0; i < numCoefficients; i++) LnmVector[i] = variables->Lnm[ii][i];
        
        // Compute particle acceleration
        for (i = variables->particleOffset[0][ii]; i <= variables->particleOffset[1][ii]; i++) {
            if (particleSystem->chargedParticles.acceleration[i].w == inboundCheck){
                dist.x = particleSystem->chargedParticles.position[i].x - boxCenter.x;
                dist.y = particleSystem->chargedParticles.position[i].y - boxCenter.y;
                dist.z = particleSystem->chargedParticles.position[i].z - boxCenter.z;
                cart2sph(r, theta, phi, dist.x, dist.y, dist.z);

                xx = cos(theta);
                yy = sin(theta);
                s2 = sqrt((1 - xx) * (1 + xx));
                fact = 1;
                pn = 1;
                for (m = 0; m < numExpansions; m++) {
                    p = pn;
                    nm = m * m + 2 * m;
                    YnmReal[nm] = constants->factorial[nm] * p;
                    p1 = p;
                    p = xx * (2 * m + 1) * p;
                    YnmRealTheta[nm] = constants->factorial[nm] * (p - (m + 1) * xx * p1) / yy;
                    for (n = m + 1; n < numExpansions; n++) {
                        nm = n * n + n + m;
                        YnmReal[nm] = constants->factorial[nm] * p;
                        p2 = p1;
                        p1 = p;
                        p = (xx * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
                        YnmRealTheta[nm] = constants->factorial[nm] * ((n - m + 1) * p - (n + 1) * xx * p1) / yy;
                    }
                    pn = -pn * fact * s2;
                    fact += 2;
                }

                accelR = 0;
                accelTheta = 0;
                accelPhi = 0;
                rn = 1;
                r /= boxSize;
                for (n = 0; n < numExpansions; n++) {
                    nm = n * n + n;
                    nms = n * (n + 1) / 2;
                    rr = n * rn / r * YnmReal[nm];
                    rtheta = rn * YnmRealTheta[nm];
                    accelR += real(rr * LnmVector[nms]);
                    accelTheta += real(rtheta * LnmVector[nms]);

                    for (m = 1; m <= n; m++) {
                        nm = n * n + n + m;
                        nms = n * (n + 1) / 2 + m;
                        eim = exp(m * phi * I);
                        rr = n * rn / r * YnmReal[nm] * eim;
                        rtheta = rn * YnmRealTheta[nm] * eim;
                        rphi = m * rn * YnmReal[nm] * eim * I;
                        accelR += 2 * real(rr * LnmVector[nms]);
                        accelTheta += 2 * real(rtheta * LnmVector[nms]);
                        accelPhi += 2 * real(rphi * LnmVector[nms]);
                    }
                    rn *= r;
                }
                accelR /= (boxSize * boxSize);
                accelTheta /= (boxSize * boxSize);
                accelPhi /= (boxSize * boxSize);

                accel.x = sin(theta) * cos(phi) * accelR + cos(theta) * cos(phi) / r * accelTheta - sin(phi) / r / sin(theta) * accelPhi;
                accel.y = sin(theta) * sin(phi) * accelR + cos(theta) * sin(phi) / r * accelTheta + cos(phi) / r / sin(theta) * accelPhi;
                accel.z = cos(theta) * accelR - sin(theta) / r * accelTheta;

    //            accel_const = PhysicalConstant::K_CONST * particleSystem->chargedParticles.position[i].w / getMass(particleSystem->chargedParticles.velocity[i].w);
                accel_const = PhysicalConstant::K_CONST * particleSystem->chargedParticles.position[i].w / particleSystem->chargedParticles.getMass(i);

                particleSystem->chargedParticles.acceleration[i].x -= accel_const * accel.x;
                particleSystem->chargedParticles.acceleration[i].y -= accel_const * accel.y;
                particleSystem->chargedParticles.acceleration[i].z -= accel_const * accel.z;
            }
        }
    }
}

void* MultipoleMethodForceCalculation::m2p_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->m2p_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::m2p_kernel_translations(int threadNum){
    vec4<int> range = rangevalues[threadNum];
    int numLevel = range.z;
    
    int ii, i, ij, jj, jb, j, n, nm, nms, m;
    vec3<int> boxIndex3D;
    vec3<int> domainBoxIndex3D;
    vec3<float> boxCenter, domainBoxMin;
    vec3<double> accel, dist;
    double boxSize, r, theta, phi, rn, accelR, accelTheta, accelPhi, accel_const;
    double xx, yy, s2, fact, pn, p, p1, p2;
    double YnmReal[FMMConstant::numExpansion2], YnmRealTheta[FMMConstant::numExpansion2];
    std::complex<double> MnmVector[numCoefficients];
    std::complex<double> rr, rtheta, rphi, I(0.0, 1.0), eim;
    int domainNum, domainMorton3;
    
    int tempLevelOffset = (numLevel == 0) ? variables->levelOffset[maxLevel] : variables->levelOffset[numLevel - 1];
    int zaligned;
    
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    
    int tempDomainIndex;
    vec3<int> domainSymOffset;
    
    boxSize = rootBoxSize / pow(3.0, numLevel);
    for (ii = range.x; ii < range.y; ii++) {
        for (ij = 0; ij < variables->numInteraction[ii]; ij++) {
            jj = variables->interactionList[ii][ij][0];
            tempDomainIndex = variables->interactionList[ii][ij][1];
            
            MultipoleMethodFunction::undoBalancedTernary(tempDomainIndex, domainSymOffset);
            
            jb = jj + tempLevelOffset;
            
            for (j = 0; j < numCoefficients; j++) MnmVector[j] = variables->Mnm[jb][j];
            
//            for (j = 0; j < numCoefficients; j++){
//                printf("variables->Mnm[%d][%d] = %.9g + i%.9g\n",jb,j,std::real(variables->Mnm[jb][j]),std::imag(variables->Mnm[jb][j]));
//            }

            domainNum = findDomain(jj, numLevel);
            domainMorton3 = variables->domainLinkage[0][domainNum].mortonIndex;
            MultipoleMethodFunction::unmorton3(domainMorton3, domainBoxIndex3D);
            domainBoxIndex3D.x -= variables->offsetSC[0].x;
            domainBoxIndex3D.y -= variables->offsetSC[0].y;
            domainBoxIndex3D.z -= variables->offsetSC[0].z;
            
            domainBoxIndex3D.x += domainSymOffset.x * boxDimVals[0];
            domainBoxIndex3D.y += domainSymOffset.y * boxDimVals[1];
            domainBoxIndex3D.z += domainSymOffset.z * boxDimVals[2];

            MultipoleMethodFunction::unmorton3(variables->domainBoxIndexFull[jb], boxIndex3D);

            domainBoxMin.x = boxMin.x + domainBoxIndex3D.x * rootBoxSize;
            domainBoxMin.y = boxMin.y + domainBoxIndex3D.y * rootBoxSize;
            domainBoxMin.z = boxMin.z + domainBoxIndex3D.z * rootBoxSize;
            boxCenter.x = domainBoxMin.x + (boxIndex3D.x + 0.5) * boxSize;
            boxCenter.y = domainBoxMin.y + (boxIndex3D.y + 0.5) * boxSize;
            boxCenter.z = domainBoxMin.z + (boxIndex3D.z + 0.5) * boxSize;
            for (i = variables->particleOffset[0][ii]; i <= variables->particleOffset[1][ii]; i++) {
                if (particleSystem->chargedParticles.acceleration[i].w == inboundCheck){
                    zaligned = 0;
                    dist.x = particleSystem->chargedParticles.position[i].x - boxCenter.x;
                    dist.y = particleSystem->chargedParticles.position[i].y - boxCenter.y;
                    dist.z = particleSystem->chargedParticles.position[i].z - boxCenter.z;
                    cart2sph(r, theta, phi, dist.x, dist.y, dist.z);
                    r /= boxSize;

                    xx = cos(theta);
                    yy = sin(theta);

                    if (std::abs(yy) < PhysicalConstant::SOFTENING){
                        yy = 1 / PhysicalConstant::SOFTENING;
                        zaligned = 1;
                    }

                    s2 = sqrt((1 - xx) * (1 + xx));
                    fact = 1;
                    pn = 1;
                    for (m = 0; m < numExpansions; m++) {
                        p = pn;
                        nm = m * m + 2 * m;
                        YnmReal[nm] = constants->factorial[nm] * p;
                        p1 = p;
                        p = xx * (2 * m + 1) * p;
                        YnmRealTheta[nm] = constants->factorial[nm] * (p - (m + 1) * xx * p1) / yy;
                        for (n = m + 1; n < numExpansions; n++) {
                            nm = n * n + n + m;
                            YnmReal[nm] = constants->factorial[nm] * p;
                            p2 = p1;
                            p1 = p;
                            p = (xx * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
                            YnmRealTheta[nm] = constants->factorial[nm] * ((n - m + 1) * p - (n + 1) * xx * p1) / yy;
                        }
                        pn = -pn * fact * s2;
                        fact += 2;
                    }
                    accelR = 0;
                    accelTheta = 0;
                    accelPhi = 0;
                    rn = 1 / r;
                    for (n = 0; n < numExpansions; n++) {
                        rn /= r;
                        nm = n * n + n;
                        nms = n * (n + 1) / 2;
                        rr = -(n + 1) * rn * YnmReal[nm];
                        rtheta = rn * r * YnmRealTheta[nm];
                        accelR += real(rr * MnmVector[nms]);
                        accelTheta += real(rtheta * MnmVector[nms]);
                        for (m = 1; m <= n; m++) {
                            nm = n * n + n + m;
                            nms = n * (n + 1) / 2 + m;
                            eim = exp(m * phi * I);
                            rr = -(n + 1) * rn * YnmReal[nm] * eim;
                            rtheta = rn * r * YnmRealTheta[nm] * eim;
                            rphi = m * rn * r * YnmReal[nm] * eim * I;
                            accelR += 2 * real(rr * MnmVector[nms]);
                            accelTheta += 2 * real(rtheta * MnmVector[nms]);
                            accelPhi += 2 * real(rphi * MnmVector[nms]);
                        }
                    }
                    accelR /= (boxSize * boxSize);
                    accelTheta /= (boxSize * boxSize);
                    accelPhi /= (boxSize * boxSize);

                    accel.x = sin(theta) * cos(phi) * accelR + cos(theta) * cos(phi) / r * accelTheta - sin(phi) / r / sin(theta) * accelPhi;
                    accel.y = sin(theta) * sin(phi) * accelR + cos(theta) * sin(phi) / r * accelTheta + cos(phi) / r / sin(theta) * accelPhi;
                    accel.z = cos(theta) * accelR - sin(theta) / r * accelTheta;

    //                accel_const = PhysicalConstant::K_CONST * particleSystem->chargedParticles.position[i].w / getMass(particleSystem->chargedParticles.velocity[i].w);
                    accel_const = PhysicalConstant::K_CONST * particleSystem->chargedParticles.position[i].w / particleSystem->chargedParticles.getMass(i);

                    if (zaligned == 1){
                        accel.x = 0;
                        accel.y = 0;
                    }

                    particleSystem->chargedParticles.acceleration[i].x -= accel_const * accel.x;
                    particleSystem->chargedParticles.acceleration[i].y -= accel_const * accel.y;
                    particleSystem->chargedParticles.acceleration[i].z -= accel_const * accel.z;
                }
            }
        }
    }
}

void* MultipoleMethodForceCalculation::sm2m_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sm2m_kernel_translations(dfc.threadNum);
    return nullptr;
};
void MultipoleMethodForceCalculation::sm2m_kernel_translations(int threadNum){
    int ib, j, jj, nfjc, jb, k, jk, jks, n, jnk, jnks, nm, i, ii, kk;
    int je_forward, je_backward;
    vec3<int> boxIndex3D, tempBoxIndex3D;
    int boxIndexFullTemp;
    double rho, powThree;
    std::complex<double> cnm, MnmScalar;
    std::complex<double> MnmVectorB[numCoefficients], MnmVectorA[numCoefficients];
        
    vec4<int> r = rangevalues[threadNum];
    int numLevel = r.w;
    
    int parentIndex;
    
    vec3<int> symmetryHalfWidth;
//    symmetryHalfWidth.x = particleSystem->chargedParticles.params.symmetryAxes.x * pow(3, particleSystem->chargedParticles.params.numSymmetryLevels - numLevel + 1) / 2;
//    symmetryHalfWidth.y = particleSystem->chargedParticles.params.symmetryAxes.y * pow(3, particleSystem->chargedParticles.params.numSymmetryLevels - numLevel + 1) / 2;
//    symmetryHalfWidth.z = particleSystem->chargedParticles.params.symmetryAxes.z * pow(3, particleSystem->chargedParticles.params.numSymmetryLevels - numLevel + 1) / 2;
//    
//    symmetryHalfWidth.x = (symmetryHalfWidth.x != 0);
//    symmetryHalfWidth.y = (symmetryHalfWidth.y != 0);
//    symmetryHalfWidth.z = (symmetryHalfWidth.z != 0);
//    
    symmetryHalfWidth.x = particleSystem->params.symmetryAxes.x * (particleSystem->params.numSymmetryLevels >= numLevel);
    symmetryHalfWidth.y = particleSystem->params.symmetryAxes.y * (particleSystem->params.numSymmetryLevels >= numLevel);
    symmetryHalfWidth.z = particleSystem->params.symmetryAxes.z * (particleSystem->params.numSymmetryLevels >= numLevel);
    
    for (jj = r.x; jj < r.y; jj++) {
        jb = jj + variables->levelOffset[numLevel + maxLevel - 1];
//        if (jj == 9) jb++;
//        if (jj == 10) jb--;
        
        boxIndexFullTemp = variables->domainLinkage[numLevel - 1][variables->domainBoxIndexFull[jb]].mortonIndex;
        nfjc = boxIndexFullTemp % 27;
        parentIndex = variables->domainLinkage[numLevel - 1][variables->domainBoxIndexFull[jb]].parent;

        // Assuming all supercells contain particleSystem->chargedParticles. Good assumption, but could change to boxIndexMask[parentIndex]
        ib = parentIndex + variables->levelOffset[numLevel + maxLevel];
        
        MultipoleMethodFunction::unmorton3(nfjc, tempBoxIndex3D);
        
        // If symmetry level, then have to include images of initial domain in expansion
        // Suppose symmetry in x, then have to include copies of multipole moments in
        // positive and negative x direction from real domain, which is denoted by
        // symmetryHalfWidth.x
        // Also, this assumes that when there is a symmetry level in some dimension,
        // the location of the box will always be equal to one. That is, if symmetry in x,
        // then tempBoxIndex3D.x == 1. This should always be the case unless the symmetry
        // domain size is mismatched in some directions, in which case a lot of changes
        // would be necessary for fmm to work.
        for (i = -symmetryHalfWidth.x; i <= symmetryHalfWidth.x; i++){
            for (ii = -symmetryHalfWidth.y; ii <= symmetryHalfWidth.y; ii++){
                for (kk = -symmetryHalfWidth.z; kk <= symmetryHalfWidth.z; kk++){
                    boxIndex3D.x = tempBoxIndex3D.x + i;
                    boxIndex3D.y = tempBoxIndex3D.y + ii;
                    boxIndex3D.z = tempBoxIndex3D.z + kk;
                    
                    rho = sqrt(abs(boxIndex3D.x - 1) + abs(boxIndex3D.y - 1) + abs(boxIndex3D.z - 1));
                    
                    boxIndex3D.x = 6 - boxIndex3D.x;
                    boxIndex3D.y = 6 - boxIndex3D.y;
                    boxIndex3D.z = 6 - boxIndex3D.z;
        
                    rotationMorton(boxIndex3D, je_forward, je_backward);        

                    // Initializer MnmVectorA with Mnm expansions from lower level cell
                    for (j = 0; j < numCoefficients; j++) {
                        MnmVectorA[j] = variables->Mnm[jb][j];
                    }

                    rotation(MnmVectorA, MnmVectorB, constants->Dnm[je_forward]);

                    // Scale size of parent cell to be equal to 1.
                    // Children cell will have side length of 1/3.
                    // rho is radius of smallest sphere that can cover a child cell
            //        rho = sqrt(3.0) / 6;

                    // Above is incorrect, rho is distance from center of child to center of its parent
                    // This distance is normalized such that length of parent is equal to one

                    // *** Note that in the calculation of cnm, use pow(3.0, n-j) vice pow(2.0, n-j) since considering 3x3x3 boxes
                    for (j = 0; j < numExpansions; j++) {
                        powThree = pow(3.0, -j);
                        for (k = 0; k <= j; k++) {
                            jk = j * j + j + k;
                            jks = j * (j + 1) / 2 + k;
                            MnmScalar = 0;
                            for (n = 0; n <= j - abs(k); n++) {
                                jnk = (j - n) * (j - n) + j - n + k;
                                jnks = (j - n) * (j - n + 1) / 2 + k;
                                nm = n * n + n;
            //                    cnm = pow(-1.0, n) * anm[nm] * anm[jnk] / anm[jk] * pow(rho, n) * Ynm[nm] * pow(3.0, n-j);
                                cnm = pow(-1.0, n) * constants->anm[nm] * constants->anm[jnk] / constants->anm[jk] * pow(rho, n) * constants->Ynm[nm] * powThree;
                                MnmScalar += MnmVectorB[jnks] * cnm;
                            }
                            MnmVectorA[jks] = MnmScalar;
                        }
                    }

                    rotation(MnmVectorA, MnmVectorB, constants->Dnm[je_backward]);

                    // Ensure not writing to the same Mnm location
#ifdef _PTHREAD_H
                    pthread_mutex_lock(&mutex_fmm);
#endif
                    for (j = 0; j < numCoefficients; j++) {
                        variables->Mnm[ib][j] += MnmVectorB[j];
            //            if (j < 10 && ib == 3 + variables->levelOffset[numLevel + maxLevel]){
            //                printf("variables->Mnm[0][%d] = %.9g\n",j,std::real(variables->Mnm[ib][j]));
            //            }
                    }
#ifdef _PTHREAD_H
                    pthread_mutex_unlock(&mutex_fmm);
#endif
                }
            }
        }
    }
}

void* MultipoleMethodForceCalculation::sm2l_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sm2l_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::sm2l_kernel_translations(int threadNum){
    int j, ii, ib, ix, iy, iz, ij, jj, jb, jx, jy, jz, k, jk, jks, n, nk, nks, jkn, jnk;
    vec3<int> boxIndex3D;
    vec3<double> dist;
    double boxSize, rho, rhoj, rhojk, rhojn;
    std::complex<double> LnmVectorA[numCoefficients], MnmVectorA[numCoefficients];
    std::complex<double> LnmVectorB[numCoefficients], MnmVectorB[numCoefficients];
    std::complex<double> cnm, LnmScalar;
        
    vec4<int> r = rangevalues[threadNum];
    int numLevel = r.w;
    
    boxSize = rootBoxSize * pow(3, numLevel);
    double boxSize2 = boxSize * boxSize;
    int je_forward, je_backward, gcd;
    
    int tempMortonIndex;
    
    int tempDomainIndex;
    vec3<int> domainSymOffset;
    
    for (ii = r.x; ii < r.y; ii++) {
        ib = ii + variables->levelOffset[numLevel + maxLevel];
        
        tempMortonIndex = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].mortonIndex;
        
        MultipoleMethodFunction::unmorton3(tempMortonIndex, boxIndex3D);
        ix = boxIndex3D.x;
        iy = boxIndex3D.y;
        iz = boxIndex3D.z;
        for (ij = 0; ij < variables->numInteraction[ii]; ij++) {            
            jj = variables->interactionList[ii][ij][0];
            tempDomainIndex = variables->interactionList[ii][ij][1];
            
            MultipoleMethodFunction::undoBalancedTernary(tempDomainIndex, domainSymOffset);
            
            jb = jj + variables->levelOffset[maxLevel + numLevel];
            tempMortonIndex = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[jb]].mortonIndex;
            
            for (j = 0; j < numCoefficients; j++){
                MnmVectorB[j] = variables->Mnm[jb][j];
            }
            
            MultipoleMethodFunction::unmorton3(tempMortonIndex, boxIndex3D);
            
            jx = boxIndex3D.x + domainSymOffset.x * variables->dimValsSC[numLevel].x;
            jy = boxIndex3D.y + domainSymOffset.y * variables->dimValsSC[numLevel].y;
            jz = boxIndex3D.z + domainSymOffset.z * variables->dimValsSC[numLevel].z;
            
            boxIndex3D.x = ix - jx;
            boxIndex3D.y = iy - jy;
            boxIndex3D.z = iz - jz;
            
            dist.x = boxIndex3D.x;
            dist.y = boxIndex3D.y;
            dist.z = boxIndex3D.z;

            rho = sqrt(dist.x * dist.x + dist.y * dist.y + dist.z * dist.z + PhysicalConstant::SOFTENING / boxSize2);
            
            // Remove gcd for rotation purposes
            // gcd should always be between 0 and 5 for the case of 3x3x3 boxes
            gcd = gcd3(boxIndex3D.x, boxIndex3D.y, boxIndex3D.z);
            switch (gcd){
                case 2:
                    boxIndex3D.x >>= 1;
                    boxIndex3D.y >>= 1;
                    boxIndex3D.z >>= 1;
                    break;
                case 3:
                    boxIndex3D.x = SIGNUM(boxIndex3D.x);
                    boxIndex3D.y = SIGNUM(boxIndex3D.y);
                    boxIndex3D.z = SIGNUM(boxIndex3D.z);
                    break;
                case 4:
                    boxIndex3D.x >>= 2;
                    boxIndex3D.y >>= 2;
                    boxIndex3D.z >>= 2;
                    break;
                case 5:
                    boxIndex3D.x = SIGNUM(boxIndex3D.x);
                    boxIndex3D.y = SIGNUM(boxIndex3D.y);
                    boxIndex3D.z = SIGNUM(boxIndex3D.z);
                    break;
                default:
                    if (gcd != 0){
                        // Should never happen
                        boxIndex3D.x /= gcd;
                        boxIndex3D.y /= gcd;
                        boxIndex3D.z /= gcd;
                    }
                    break;
            }
            
            boxIndex3D.x += 5;
            boxIndex3D.y += 5;
            boxIndex3D.z += 5;
            
            rotationMorton(boxIndex3D, je_forward, je_backward);

            rotation(MnmVectorB, MnmVectorA, constants->Dnm[je_forward]);

            rhoj = 1;
            for (j = 0; j < numExpansions; j++) {
                rhojk = rhoj;
                rhoj *= rho;
                for (k = 0; k <= j; k++) {
                    jk = j * j + j + k;
                    jks = j * (j + 1) / 2 + k;
                    LnmScalar = 0;
                    rhojn = rhojk;
                    rhojk *= rho;
                    for (n = abs(k); n < numExpansions; n++) {
                        rhojn *= rho;
                        nk = n * n + n + k;
                        nks = n * (n + 1) / 2 + k;
                        jkn = jk * FMMConstant::numExpansion2 + nk;
                        jnk = (j + n) * (j + n) + j + n;
                        cnm = constants->Anm[jkn] / rhojn * constants->Ynm[jnk];
                        LnmScalar += MnmVectorA[nks] * cnm;
                    }
                    LnmVectorA[jks] = LnmScalar;
                }
            }

            rotation(LnmVectorA, LnmVectorB, constants->Dnm[je_backward]);
            
            // Assign Lnm
            for (j = 0; j < numExpansions; j++){
                for (n = 0; n <= j; n++){
                    variables->Lnm[ii][j*(j+1)/2 + n] += LnmVectorB[j*(j+1)/2 + n];
                }
            }            
        }
    }
}

//void* MultipoleMethodForceCalculation::sl2l_threaded_translations(void* dF){
//    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
//    dfc.d->sl2l_kernel_translations(dfc.threadNum);
//    return nullptr;
//}
//
//void MultipoleMethodForceCalculation::sl2l_kernel_translations(int threadNum){
//    int ii, ib, i, nfip, nfic, j, k, jk, jks, n, jnk, nk, nks;
//    vec3<int> boxIndex3D;
//    double rho;
//    std::complex<double> cnm, LnmScalar;
//    std::complex<double> LnmVectorA[numCoefficients], LnmVectorB[numCoefficients];
//    
//    vec4<int> r = rangevalues[threadNum];
//    int numLevel = r.w;
//    
//    int je_forward, je_backward;
//    
//    for (ii = r.x; ii < r.y; ii++) {
//        ib = ii + variables->levelOffset[numLevel + maxLevel];
//        
//        nfip = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].parent;
//        nfic = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].mortonIndex % 27;
//        MultipoleMethodFunction::unmorton3(nfic, boxIndex3D);
//        
//        // rho is distance from center of child to center of parent cell
//        // rho is normalized, such that child cell has length equal to one
//        rho = sqrt(abs(boxIndex3D.x - 1) + abs(boxIndex3D.y - 1) + abs(boxIndex3D.z - 1));
//
//        boxIndex3D.x = boxIndex3D.x + 4;
//        boxIndex3D.y = boxIndex3D.y + 4;
//        boxIndex3D.z = boxIndex3D.z + 4;
//        
//        rotationMorton(boxIndex3D, je_forward, je_backward);
//        
//        ib = variables->neo[nfip];
//        
//        for (i = 0; i < numCoefficients; i++) {
//            LnmVectorA[i] = variables->LnmOld[ib][i];
//        }
//
//        rotation(LnmVectorA, LnmVectorB, constants->Dnm[je_forward]);
//
//        //rho = sqrt(3.0) / 2;
//        
//        // *** NOTE - in the calculation for cnm, use pow(3.0, -n-1) vice pow(2.0, -n-1) since considering 3x3x3 cells
//        for (j = 0; j < numExpansions; j++){
//            for (k = 0; k <= j; k++){
//                jk = j * j + j + k;
//                jks = j*(j+1)/2+k;
//                LnmScalar = 0;
//                for (n = j; n < numExpansions; n++){
//                    jnk = (n-j)*(n-j)+n-j;
//                    nk = n*n+n+k;
//                    nks = n*(n+1)/2+k;
//                    cnm = constants->anm[jnk] * constants->anm[jk] / constants->anm[nk] * pow(rho, n-j) * constants->Ynm[jnk] * pow(3.0, -n-1);
//                    LnmScalar += LnmVectorB[nks] * cnm;
//                }
//                LnmVectorA[jks] = LnmScalar;
//            }
//        }
//        
//        rotation(LnmVectorA, LnmVectorB, constants->Dnm[je_backward]);
//        for (i = 0; i < numCoefficients; i++) {
//            variables->Lnm[ii][i] = LnmVectorB[i];
////            if (i < 10) {
////                printf("variables->Lnm[%d][%d] = %.9g + i %.9g\n", ii, i, std::real(variables->Lnm[ii][i]), std::imag(variables->Lnm[ii][i]));
////            }
//        }
//    }
//}

void* MultipoleMethodForceCalculation::sm2p_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sm2p_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::sm2p_kernel_translations(int threadNum){
    vec4<int> range = rangevalues[threadNum];
    int numLevel = range.z;
    
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    
    // Could store temporary variables for boxSize * boxSize, and the trig functions to reduce repeat calculations
    
    int ii, i, ib, ij, jj, jb, j, n, nm, nms, m;
    vec3<int> boxIndex3D;
    vec3<float> boxCenter;
    vec3<double> accel, dist;
    double boxSize, r, theta, phi, rn, accelR, accelTheta, accelPhi, accel_const;
    double xx, yy, s2, fact, pn, p, p1, p2;
    double YnmReal[FMMConstant::numExpansion2], YnmRealTheta[FMMConstant::numExpansion2];
    std::complex<double> MnmVector[numCoefficients], MnmVectorB[numCoefficients];
    std::complex<double> rr, rtheta, rphi, I(0.0, 1.0), eim;
    int domainMorton3;

    boxSize = rootBoxSize * pow(3, numLevel);
    int zaligned, je_forward, je_backward;
    
    int tempDomainIndex;
    vec3<int> domainSymOffset;
    
    for (ii = range.x; ii < range.y; ii++) {
        ib = ii + variables->levelOffset[numLevel + maxLevel];
//        printf("first particle : %d, last particle : %d\n",variables->domainLinkage[numLevel][ii].firstParticle,variables->domainLinkage[numLevel][ii].lastParticle);
        for (ij = 0; ij < variables->numInteraction[ii]; ij++) {
            jj = variables->interactionList[ii][ij][0];
            tempDomainIndex = variables->interactionList[ii][ij][1];
            
            MultipoleMethodFunction::undoBalancedTernary(tempDomainIndex, domainSymOffset);
            
            jb = jj + variables->levelOffset[numLevel + maxLevel];

            for (j = 0; j < numCoefficients; j++) MnmVector[j] = variables->Mnm[jb][j];

            domainMorton3 = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[jb]].mortonIndex;

            MultipoleMethodFunction::unmorton3(domainMorton3, boxIndex3D);
            boxIndex3D.x -= variables->offsetSC[numLevel].x;
            boxIndex3D.y -= variables->offsetSC[numLevel].y;
            boxIndex3D.z -= variables->offsetSC[numLevel].z;
            
            boxIndex3D.x += domainSymOffset.x * variables->dimValsSC[numLevel].x;
            boxIndex3D.y += domainSymOffset.y * variables->dimValsSC[numLevel].y;
            boxIndex3D.z += domainSymOffset.z * variables->dimValsSC[numLevel].z;

            boxCenter.x = variables->boxMinSC[numLevel].x + (boxIndex3D.x + 0.5) * boxSize;
            boxCenter.y = variables->boxMinSC[numLevel].y + (boxIndex3D.y + 0.5) * boxSize;
            boxCenter.z = variables->boxMinSC[numLevel].z + (boxIndex3D.z + 0.5) * boxSize;
            for (i = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].firstParticle; i <= variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].lastParticle; i++) {
                if (particleSystem->chargedParticles.acceleration[i].w == inboundCheck){
                    zaligned = 0;

                    dist.x = particleSystem->chargedParticles.position[i].x - boxCenter.x;
                    dist.y = particleSystem->chargedParticles.position[i].y - boxCenter.y;
                    dist.z = particleSystem->chargedParticles.position[i].z - boxCenter.z;
                    cart2sph(r, theta, phi, dist.x, dist.y, dist.z);
                    r /= boxSize;
                    xx = cos(theta);
                    yy = sin(theta);
                    if (std::abs(yy) < PhysicalConstant::SOFTENING){
                        yy = 1 / PhysicalConstant::SOFTENING;
                        zaligned = 1;
                    }
                    s2 = sqrt((1 - xx) * (1 + xx));
                    fact = 1;
                    pn = 1;
                    for (m = 0; m < numExpansions; m++) {
                        p = pn;
                        nm = m * m + 2 * m;
                        YnmReal[nm] = constants->factorial[nm] * p;
                        p1 = p;
                        p = xx * (2 * m + 1) * p;
                        YnmRealTheta[nm] = constants->factorial[nm] * (p - (m + 1) * xx * p1) / yy;
                        for (n = m + 1; n < numExpansions; n++) {
                            nm = n * n + n + m;
                            YnmReal[nm] = constants->factorial[nm] * p;
                            p2 = p1;
                            p1 = p;
                            p = (xx * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
                            YnmRealTheta[nm] = constants->factorial[nm] * ((n - m + 1) * p - (n + 1) * xx * p1) / yy;
                        }
                        pn = -pn * fact * s2;
                        fact += 2;
                    }
                    accelR = 0;
                    accelTheta = 0;
                    accelPhi = 0;
                    rn = 1 / r;
                    for (n = 0; n < numExpansions; n++) {
                        rn /= r;
                        nm = n * n + n;
                        nms = n * (n + 1) / 2;
                        rr = -(n + 1) * rn * YnmReal[nm];
                        rtheta = rn * r * YnmRealTheta[nm];
                        accelR += real(rr * MnmVector[nms]);
                        accelTheta += real(rtheta * MnmVector[nms]);
                        for (m = 1; m <= n; m++) {
                            nm = n * n + n + m;
                            nms = n * (n + 1) / 2 + m;
                            eim = exp(m * phi * I);
                            rr = -(n + 1) * rn * YnmReal[nm] * eim;
                            rtheta = rn * r * YnmRealTheta[nm] * eim;
                            rphi = m * rn * r * YnmReal[nm] * eim * I;
                            accelR += 2 * real(rr * MnmVector[nms]);
                            accelTheta += 2 * real(rtheta * MnmVector[nms]);
                            accelPhi += 2 * real(rphi * MnmVector[nms]);
                        }
                    }
                    accelR /= (boxSize * boxSize);
                    accelTheta /= (boxSize * boxSize);
                    accelPhi /= (boxSize * boxSize);

                    accel.x = sin(theta) * cos(phi) * accelR + cos(theta) * cos(phi) / r * accelTheta - sin(phi) / r / sin(theta) * accelPhi;
                    accel.y = sin(theta) * sin(phi) * accelR + cos(theta) * sin(phi) / r * accelTheta + cos(phi) / r / sin(theta) * accelPhi;
                    accel.z = cos(theta) * accelR - sin(theta) / r * accelTheta;
    //                accel_const = PhysicalConstant::K_CONST * particleSystem->chargedParticles.position[i].w / getMass(particleSystem->chargedParticles.velocity[i].w);
                    accel_const = PhysicalConstant::K_CONST * particleSystem->chargedParticles.position[i].w / particleSystem->chargedParticles.getMass(i);

                    if (zaligned == 1){

    //                    for (j = 0; j < numCoefficients; j++){
    //                        MnmVector[j] = MnmVectorB[j];
    //                    }
    //                    
    //                    xx = accel.x;
    //                    accel.x = accel.z;
    //                    accel.z = xx;

                        accel.x = 0;
                        accel.y = 0;
                    }

                    particleSystem->chargedParticles.acceleration[i].x -= accel_const * accel.x;
                    particleSystem->chargedParticles.acceleration[i].y -= accel_const * accel.y;
                    particleSystem->chargedParticles.acceleration[i].z -= accel_const * accel.z;

                }
            }
        }
    }
}

void* MultipoleMethodForceCalculation::sp2p_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sp2p_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::sp2p_kernel_translations(int threadNum){
    int ii, ij, jj, i, j;
    vec3<double> dist;
    double accel_const, invDist, invDistCube, s;
    vec3<double> ai;
    
    int tempDomainIndex;
    LatticePoint tempLattice;
    
    vec4<int> r = rangevalues[threadNum];

    vec4<int> particleBounds;
    
    for (ii = r.x; ii < r.y; ii++){
        particleBounds.x = variables->domainLinkage[0][variables->domainBoxIndexFull[ii]].firstParticle;
        particleBounds.y = variables->domainLinkage[0][variables->domainBoxIndexFull[ii]].lastParticle + 1;
        for (ij = 0; ij < variables->numInteraction[ii]; ij++){
            jj = variables->interactionList[ii][ij][0];
            tempDomainIndex = variables->interactionList[ii][ij][1];
            MultipoleMethodFunction::undoBalancedTernary(tempDomainIndex, tempLattice.point);
            
            particleBounds.z = variables->domainLinkage[0][variables->domainBoxIndexFull[jj]].firstParticle;
            particleBounds.w = variables->domainLinkage[0][variables->domainBoxIndexFull[jj]].lastParticle + 1;
            
            if (vectorizedToggle){
                DirectForceFunction::direct_vectorized_translational(tempLattice, particleSystem->params.domainSize,
                        particleBounds, particleSystem->chargedParticles, jptc);
            }
            else {
                DirectForceFunction::direct_translational(tempLattice, particleSystem->params.domainSize,
                        particleBounds, particleSystem->chargedParticles);
            }
        }
    }
}

void* MultipoleMethodForceCalculation::p2outlier_p_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->p2outlier_p_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::p2outlier_p_kernel_translations(int threadNum){
    int ii, ij, jj, i, j;
    vec3<double> dist;
    double accel_const, invDist, invDistCube, s;
    vec3<double> ai;
    
    int tempDomainIndex;
    LatticePoint tempLattice;
    
    vec4<int> r = rangevalues[threadNum];
    
    vec4<int> particleBounds;
    
    for (ii = r.x; ii < r.y; ii++){
        particleBounds.x = variables->outlierParticleOffset[0][ii];
        particleBounds.y = variables->outlierParticleOffset[1][ii] + 1;
        for (ij = 0; ij < variables->numInteraction[ii]; ij++){
            jj = variables->interactionList[ii][ij][0];
            tempDomainIndex = variables->interactionList[ii][ij][1];
            MultipoleMethodFunction::undoBalancedTernary(tempDomainIndex, tempLattice.point);
            
            particleBounds.z = variables->particleOffset[0][jj];
            particleBounds.w = variables->particleOffset[1][jj] + 1;
            
            if (vectorizedToggle){
                DirectForceFunction::direct_vectorized_translational(tempLattice, particleSystem->params.domainSize,
                        particleBounds, particleSystem->chargedParticles, jptc);
            }
            else {
                DirectForceFunction::direct_translational(tempLattice, particleSystem->params.domainSize,
                        particleBounds, particleSystem->chargedParticles);
            }
        }
    }
}

void* MultipoleMethodForceCalculation::sp2outlier_p_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sp2outlier_p_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::sp2outlier_p_kernel_translations(int threadNum){
    int ii, ij, jj, i, j;
    vec3<double> dist;
    double accel_const, invDist, invDistCube, s;
    vec3<double> ai;
    
    int tempDomainIndex;
    LatticePoint tempLattice;
    
    vec4<int> r = rangevalues[threadNum];

    vec4<int> particleBounds;
    
    for (ii = r.x; ii < r.y; ii++){
        particleBounds.x = variables->outlierDomainLinkage[ii].firstParticle;
        particleBounds.y = variables->outlierDomainLinkage[ii].lastParticle + 1;
        for (ij = 0; ij < variables->numInteraction[ii]; ij++){
            jj = variables->interactionList[ii][ij][0];
            tempDomainIndex = variables->interactionList[ii][ij][1];
            MultipoleMethodFunction::undoBalancedTernary(tempDomainIndex, tempLattice.point);
            
            particleBounds.z = variables->domainLinkage[0][variables->domainBoxIndexFull[jj]].firstParticle;
            particleBounds.w = variables->domainLinkage[0][variables->domainBoxIndexFull[jj]].lastParticle + 1;
            
            if (vectorizedToggle){
                DirectForceFunction::direct_vectorized_translational(tempLattice, particleSystem->params.domainSize,
                        particleBounds, particleSystem->chargedParticles, jptc);
            }
            else {
                DirectForceFunction::direct_translational(tempLattice, particleSystem->params.domainSize,
                        particleBounds, particleSystem->chargedParticles);
            }
        }
    }
}

void* MultipoleMethodForceCalculation::m2outlier_p_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->m2outlier_p_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::m2outlier_p_kernel_translations(int threadNum){
    vec4<int> range = rangevalues[threadNum];
    int numLevel = range.z;
    
    int ii, i, ij, jj, jb, j, n, nm, nms, m;
    vec3<int> boxIndex3D;
    vec3<int> domainBoxIndex3D;
    vec3<float> boxCenter, domainBoxMin;
    vec3<double> accel, dist;
    double boxSize, r, theta, phi, rn, accelR, accelTheta, accelPhi, accel_const;
    double xx, yy, s2, fact, pn, p, p1, p2;
    double YnmReal[FMMConstant::numExpansion2], YnmRealTheta[FMMConstant::numExpansion2];
    std::complex<double> MnmVector[numCoefficients];
    std::complex<double> rr, rtheta, rphi, I(0.0, 1.0), eim;
    int domainNum, domainMorton3;
    
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    
    int tempLevelOffset = (numLevel == 0) ? variables->levelOffset[maxLevel] : variables->levelOffset[numLevel - 1];
    int zaligned;

    int tempDomainIndex;
    vec3<int> domainSymOffset;
    
    boxSize = rootBoxSize / pow(3.0, numLevel);
    for (ii = range.x; ii < range.y; ii++) {
        for (ij = 0; ij < variables->numInteraction[ii]; ij++) {
            jj = variables->interactionList[ii][ij][0];
            tempDomainIndex = variables->interactionList[ii][ij][1];
            MultipoleMethodFunction::undoBalancedTernary(tempDomainIndex, domainSymOffset);
            jb = jj + tempLevelOffset;

            for (j = 0; j < numCoefficients; j++) MnmVector[j] = variables->Mnm[jb][j];

            domainNum = findDomain(jj, numLevel);
            domainMorton3 = variables->domainLinkage[0][domainNum].mortonIndex;
            MultipoleMethodFunction::unmorton3(domainMorton3, domainBoxIndex3D);
            domainBoxIndex3D.x -= variables->offsetSC[0].x;
            domainBoxIndex3D.y -= variables->offsetSC[0].y;
            domainBoxIndex3D.z -= variables->offsetSC[0].z;
            
            domainBoxIndex3D.x += domainSymOffset.x * boxDimVals[0];
            domainBoxIndex3D.y += domainSymOffset.y * boxDimVals[1];
            domainBoxIndex3D.z += domainSymOffset.z * boxDimVals[2];

            MultipoleMethodFunction::unmorton3(variables->domainBoxIndexFull[jb], boxIndex3D);

            domainBoxMin.x = boxMin.x + domainBoxIndex3D.x * rootBoxSize;
            domainBoxMin.y = boxMin.y + domainBoxIndex3D.y * rootBoxSize;
            domainBoxMin.z = boxMin.z + domainBoxIndex3D.z * rootBoxSize;
            boxCenter.x = domainBoxMin.x + (boxIndex3D.x + 0.5) * boxSize;
            boxCenter.y = domainBoxMin.y + (boxIndex3D.y + 0.5) * boxSize;
            boxCenter.z = domainBoxMin.z + (boxIndex3D.z + 0.5) * boxSize;
            for (i = variables->outlierParticleOffset[0][ii]; i <= variables->outlierParticleOffset[1][ii]; i++) {
                if (particleSystem->chargedParticles.acceleration[i].w == inboundCheck){
//                
//                if (i == 370){
//                    printf("test\n");
//                }
                
                    zaligned = 0;
                    dist.x = particleSystem->chargedParticles.position[i].x - boxCenter.x;
                    dist.y = particleSystem->chargedParticles.position[i].y - boxCenter.y;
                    dist.z = particleSystem->chargedParticles.position[i].z - boxCenter.z;
                    cart2sph(r, theta, phi, dist.x, dist.y, dist.z);
                    r /= boxSize;

                    xx = cos(theta);
                    yy = sin(theta);

                    if (std::abs(yy) < PhysicalConstant::SOFTENING){
                        yy = 1 / PhysicalConstant::SOFTENING;
                        zaligned = 1;
                    }

                    s2 = sqrt((1 - xx) * (1 + xx));
                    fact = 1;
                    pn = 1;
                    for (m = 0; m < numExpansions; m++) {
                        p = pn;
                        nm = m * m + 2 * m;
                        YnmReal[nm] = constants->factorial[nm] * p;
                        p1 = p;
                        p = xx * (2 * m + 1) * p;
                        YnmRealTheta[nm] = constants->factorial[nm] * (p - (m + 1) * xx * p1) / yy;
                        for (n = m + 1; n < numExpansions; n++) {
                            nm = n * n + n + m;
                            YnmReal[nm] = constants->factorial[nm] * p;
                            p2 = p1;
                            p1 = p;
                            p = (xx * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
                            YnmRealTheta[nm] = constants->factorial[nm] * ((n - m + 1) * p - (n + 1) * xx * p1) / yy;
                        }
                        pn = -pn * fact * s2;
                        fact += 2;
                    }
                    accelR = 0;
                    accelTheta = 0;
                    accelPhi = 0;
                    rn = 1 / r;
                    for (n = 0; n < numExpansions; n++) {
                        rn /= r;
                        nm = n * n + n;
                        nms = n * (n + 1) / 2;
                        rr = -(n + 1) * rn * YnmReal[nm];
                        rtheta = rn * r * YnmRealTheta[nm];
                        accelR += real(rr * MnmVector[nms]);
                        accelTheta += real(rtheta * MnmVector[nms]);
                        for (m = 1; m <= n; m++) {
                            nm = n * n + n + m;
                            nms = n * (n + 1) / 2 + m;
                            eim = exp(m * phi * I);
                            rr = -(n + 1) * rn * YnmReal[nm] * eim;
                            rtheta = rn * r * YnmRealTheta[nm] * eim;
                            rphi = m * rn * r * YnmReal[nm] * eim * I;
                            accelR += 2 * real(rr * MnmVector[nms]);
                            accelTheta += 2 * real(rtheta * MnmVector[nms]);
                            accelPhi += 2 * real(rphi * MnmVector[nms]);
                        }
                    }
                    accelR /= (boxSize * boxSize);
                    accelTheta /= (boxSize * boxSize);
                    accelPhi /= (boxSize * boxSize);

                    accel.x = sin(theta) * cos(phi) * accelR + cos(theta) * cos(phi) / r * accelTheta - sin(phi) / r / sin(theta) * accelPhi;
                    accel.y = sin(theta) * sin(phi) * accelR + cos(theta) * sin(phi) / r * accelTheta + cos(phi) / r / sin(theta) * accelPhi;
                    accel.z = cos(theta) * accelR - sin(theta) / r * accelTheta;
    //                accel_const = PhysicalConstant::K_CONST * particleSystem->chargedParticles.position[i].w / getMass(particleSystem->chargedParticles.velocity[i].w);
                    accel_const = PhysicalConstant::K_CONST * particleSystem->chargedParticles.position[i].w / particleSystem->chargedParticles.getMass(i);

                    if (zaligned == 1){
                        accel.x = 0;
                        accel.y = 0;
                    }

                    particleSystem->chargedParticles.acceleration[i].x -= accel_const * accel.x;
                    particleSystem->chargedParticles.acceleration[i].y -= accel_const * accel.y;
                    particleSystem->chargedParticles.acceleration[i].z -= accel_const * accel.z;
                }
            }
        }
    }
}

void* MultipoleMethodForceCalculation::sm2outlier_p_threaded_translations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sm2outlier_p_kernel_translations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::sm2outlier_p_kernel_translations(int threadNum){
    vec4<int> range = rangevalues[threadNum];
    int numLevel = range.z;
    
    // Could store temporary variables for boxSize * boxSize, and the trig functions to reduce repeat calculations
    
    int ii, i, ib, ij, jj, jb, j, n, nm, nms, m;
    vec3<int> boxIndex3D;
    vec3<float> boxCenter;
    vec3<double> accel, dist;
    double boxSize, r, theta, phi, rn, accelR, accelTheta, accelPhi, accel_const;
    double xx, yy, s2, fact, pn, p, p1, p2;
    double YnmReal[FMMConstant::numExpansion2], YnmRealTheta[FMMConstant::numExpansion2];
    std::complex<double> MnmVector[numCoefficients], MnmVectorB[numCoefficients];
    std::complex<double> rr, rtheta, rphi, I(0.0, 1.0), eim;
    int domainMorton3;
    
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);

    boxSize = rootBoxSize * pow(3, numLevel);
    int zaligned, je_forward, je_backward;
    
    int tempDomainIndex;
    vec3<int> domainSymOffset;
    
    for (ii = range.x; ii < range.y; ii++) {
//        printf("first particle : %d, last particle : %d\n",variables->domainLinkage[numLevel][ii].firstParticle,variables->domainLinkage[numLevel][ii].lastParticle);
        for (ij = 0; ij < variables->numInteraction[ii]; ij++) {
            jj = variables->interactionList[ii][ij][0];
            tempDomainIndex = variables->interactionList[ii][ij][1];
            MultipoleMethodFunction::undoBalancedTernary(tempDomainIndex, domainSymOffset);
            jb = jj + variables->levelOffset[numLevel + maxLevel];

            for (j = 0; j < numCoefficients; j++) MnmVector[j] = variables->Mnm[jb][j];

            domainMorton3 = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[jb]].mortonIndex;

            MultipoleMethodFunction::unmorton3(domainMorton3, boxIndex3D);
            boxIndex3D.x -= variables->offsetSC[numLevel].x;
            boxIndex3D.y -= variables->offsetSC[numLevel].y;
            boxIndex3D.z -= variables->offsetSC[numLevel].z;
            
            boxIndex3D.x += domainSymOffset.x * variables->dimValsSC[numLevel].x;
            boxIndex3D.y += domainSymOffset.y * variables->dimValsSC[numLevel].y;
            boxIndex3D.z += domainSymOffset.z * variables->dimValsSC[numLevel].z;

            boxCenter.x = variables->boxMinSC[numLevel].x + (boxIndex3D.x + 0.5) * boxSize;
            boxCenter.y = variables->boxMinSC[numLevel].y + (boxIndex3D.y + 0.5) * boxSize;
            boxCenter.z = variables->boxMinSC[numLevel].z + (boxIndex3D.z + 0.5) * boxSize;
            for (i = variables->outlierDomainLinkage[ii].firstParticle; i <= variables->outlierDomainLinkage[ii].lastParticle; i++) {
                if (particleSystem->chargedParticles.acceleration[i].w == inboundCheck){
                    zaligned = 0;

                    dist.x = particleSystem->chargedParticles.position[i].x - boxCenter.x;
                    dist.y = particleSystem->chargedParticles.position[i].y - boxCenter.y;
                    dist.z = particleSystem->chargedParticles.position[i].z - boxCenter.z;
                    cart2sph(r, theta, phi, dist.x, dist.y, dist.z);
                    r /= boxSize;
                    xx = cos(theta);
                    yy = sin(theta);
                    if (std::abs(yy) < PhysicalConstant::SOFTENING){
                        printf("test\n");
                        yy = 1 / PhysicalConstant::SOFTENING;
                        zaligned = 1;
                    }
                    s2 = sqrt((1 - xx) * (1 + xx));
                    fact = 1;
                    pn = 1;
                    for (m = 0; m < numExpansions; m++) {
                        p = pn;
                        nm = m * m + 2 * m;
                        YnmReal[nm] = constants->factorial[nm] * p;
                        p1 = p;
                        p = xx * (2 * m + 1) * p;
                        YnmRealTheta[nm] = constants->factorial[nm] * (p - (m + 1) * xx * p1) / yy;
                        for (n = m + 1; n < numExpansions; n++) {
                            nm = n * n + n + m;
                            YnmReal[nm] = constants->factorial[nm] * p;
                            p2 = p1;
                            p1 = p;
                            p = (xx * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
                            YnmRealTheta[nm] = constants->factorial[nm] * ((n - m + 1) * p - (n + 1) * xx * p1) / yy;
                        }
                        pn = -pn * fact * s2;
                        fact += 2;
                    }
                    accelR = 0;
                    accelTheta = 0;
                    accelPhi = 0;
                    rn = 1 / r;
                    for (n = 0; n < numExpansions; n++) {
                        rn /= r;
                        nm = n * n + n;
                        nms = n * (n + 1) / 2;
                        rr = -(n + 1) * rn * YnmReal[nm];
                        rtheta = rn * r * YnmRealTheta[nm];
                        accelR += real(rr * MnmVector[nms]);
                        accelTheta += real(rtheta * MnmVector[nms]);
                        for (m = 1; m <= n; m++) {
                            nm = n * n + n + m;
                            nms = n * (n + 1) / 2 + m;
                            eim = exp(m * phi * I);
                            rr = -(n + 1) * rn * YnmReal[nm] * eim;
                            rtheta = rn * r * YnmRealTheta[nm] * eim;
                            rphi = m * rn * r * YnmReal[nm] * eim * I;
                            accelR += 2 * real(rr * MnmVector[nms]);
                            accelTheta += 2 * real(rtheta * MnmVector[nms]);
                            accelPhi += 2 * real(rphi * MnmVector[nms]);
                        }
                    }
                    accelR /= (boxSize * boxSize);
                    accelTheta /= (boxSize * boxSize);
                    accelPhi /= (boxSize * boxSize);

                    accel.x = sin(theta) * cos(phi) * accelR + cos(theta) * cos(phi) / r * accelTheta - sin(phi) / r / sin(theta) * accelPhi;
                    accel.y = sin(theta) * sin(phi) * accelR + cos(theta) * sin(phi) / r * accelTheta + cos(phi) / r / sin(theta) * accelPhi;
                    accel.z = cos(theta) * accelR - sin(theta) / r * accelTheta;
    //                accel_const = PhysicalConstant::K_CONST * particleSystem->chargedParticles.position[i].w / getMass(particleSystem->chargedParticles.velocity[i].w);
                    accel_const = PhysicalConstant::K_CONST * particleSystem->chargedParticles.position[i].w / particleSystem->chargedParticles.getMass(i);

                    if (zaligned == 1){
                        accel.x = 0;
                        accel.y = 0;
                    }


                    particleSystem->chargedParticles.acceleration[i].x -= accel_const * accel.x;
                    particleSystem->chargedParticles.acceleration[i].y -= accel_const * accel.y;
                    particleSystem->chargedParticles.acceleration[i].z -= accel_const * accel.z;
                }
            }
        }
    }
}
