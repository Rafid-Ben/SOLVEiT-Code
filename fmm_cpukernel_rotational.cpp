/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "force_calculation_methods.h"

void* MultipoleMethodForceCalculation::direct_threaded_rotations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->direct_kernel_rotations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::direct_kernel_rotations(int threadNum){
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

void* MultipoleMethodForceCalculation::p2p_threaded_rotations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->p2p_kernel_rotations(dfc.threadNum);
    return NULL;
}

void MultipoleMethodForceCalculation::p2p_kernel_rotations(int threadNum){
    int ii, ij, jj;
    int tempQuad;
    
    vec4<int> r = rangevalues[threadNum];
    
    vec4<int> particleBounds;
    
    for (ii = r.x; ii < r.y; ii++){
        particleBounds.x = variables->particleOffset[0][ii];
        particleBounds.y = variables->particleOffset[1][ii] + 1;
        for (ij = 0; ij < variables->numInteraction[ii]; ij++){
            jj = variables->interactionList[ii][ij][0];
            tempQuad = variables->interactionList[ii][ij][1];
            particleBounds.z = variables->particleOffset[0][jj];
            particleBounds.w = variables->particleOffset[1][jj] + 1;
            
            if (vectorizedToggle){
                DirectForceFunction::direct_vectorized_rotational(tempQuad, particleBounds,
                        particleSystem->chargedParticles, jptc);
            }
            else {
                DirectForceFunction::direct_rotational(tempQuad, particleBounds,
                        particleSystem->chargedParticles);
            }
        }
    }
}

void* MultipoleMethodForceCalculation::sp2p_threaded_rotations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sp2p_kernel_rotations(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::sp2p_kernel_rotations(int threadNum){
    int ii, ij, jj, i, j;
    vec3<double> dist;
    double accel_const, invDist, invDistCube, s;
    vec3<double> ai;
    
    vec4<int> r = rangevalues[threadNum];

    vec4<int> particleBounds;
    
    for (ii = r.x; ii < r.y; ii++){
        particleBounds.x = variables->domainLinkage[0][variables->domainBoxIndexFull[ii]].firstParticle;
        particleBounds.y = variables->domainLinkage[0][variables->domainBoxIndexFull[ii]].lastParticle + 1;
        for (ij = 0; ij < variables->numInteraction[ii]; ij++){
            jj = variables->interactionList[ii][ij][0];
            particleBounds.z = variables->domainLinkage[0][variables->domainBoxIndexFull[jj]].firstParticle;
            particleBounds.w = variables->domainLinkage[0][variables->domainBoxIndexFull[jj]].lastParticle + 1;
            
            for (i = 0; i < 4; i++){
                if (vectorizedToggle){
                    DirectForceFunction::direct_vectorized_rotational(i, particleBounds,
                            particleSystem->chargedParticles, jptc);
                }
                else {
                    DirectForceFunction::direct_rotational(i, particleBounds,
                            particleSystem->chargedParticles);
                }
            }
        }
    }
}

void* MultipoleMethodForceCalculation::m2m_threaded_rotations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->m2m_kernel_rotations(dfc.threadNum);
    return NULL;
}

void MultipoleMethodForceCalculation::m2m_kernel_rotations(int threadNum){
    int ib, j, jj, nfjp, nfjc, jb, k, jk, jks, n, jnk, jnks, nm;
    int je_forward, je_backward;
    vec3<int> boxIndex3D;
    int boxIndexFullTemp;
    double rho, powTwo;
    std::complex<double> cnm, MnmScalar;
    std::complex<double> MnmVectorB[numCoefficients], MnmVectorA[numCoefficients];
        
    vec4<int> r = rangevalues[threadNum];
    int numLevel = r.w;
    
    int domainNum, domainBias;
    int domainNumBoxIndexTemp = 1 << (3 * numLevel);
    
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
        nfjp = boxIndexFullTemp / 8;
        nfjc = boxIndexFullTemp % 8;
        if (numLevel == 0){
            ib = variables->domainBoxIndexMask[domainNum] + tempLevelOffset;
        }
        else{
            ib = variables->domainBoxIndexMask[domainBias + nfjp] + tempLevelOffset;
        }
        
        if (numLevel == 0){
            ib = variables->domainBoxIndexMask[domainNum] + tempLevelOffset;
            MultipoleMethodFunction::unmorton(nfjc, boxIndex3D);
            boxIndex3D.z = 6 - boxIndex3D.z * 2;
            for (tempQuad = 0; tempQuad < 4; tempQuad++){
                switch (tempQuad){
                    case 0:
                        boxIndex3D.x = 1;
                        boxIndex3D.y = 1;
                        break;
                    case 1:
                        boxIndex3D.x = 0;
                        boxIndex3D.y = 1;
                        break;
                    case 2:
                        boxIndex3D.x = 0;
                        boxIndex3D.y = 0;
                        break;
                    case 3:
                        boxIndex3D.x = 1;
                        boxIndex3D.y = 0;
                        break;
                }
                boxIndex3D.x = 6 - boxIndex3D.x * 2;
                boxIndex3D.y = 6 - boxIndex3D.y * 2;
                
                rotationMorton(boxIndex3D, je_forward, je_backward);
                
                for (j = 0; j < numExpansions; j++){
                    for (m = 0; m <= j; m++){
                        mr = qrf_complex(-m * tempQuad);
                        MnmVectorA[j * (j + 1) / 2 + m] = variables->Mnm[jb][j * (j + 1) / 2 + m] * mr;
                    }
                }
                
                rotation(MnmVectorA, MnmVectorB, constants->Dnm[je_forward]);
                
                rho = sqrt(3.0)/2;
        
                for (j = 0; j < numExpansions; j++) {
                    powTwo = pow(2.0, -j);
                    for (k = 0; k <= j; k++) {
                        jk = j * j + j + k;
                        jks = j * (j + 1) / 2 + k;
                        MnmScalar = 0;
                        for (n = 0; n <= j - abs(k); n++) {
                            jnk = (j - n) * (j - n) + j - n + k;
                            jnks = (j - n) * (j - n + 1) / 2 + k;
                            nm = n * n + n;
        //                    cnm = pow(-1.0, n) * anm[nm] * anm[jnk] / anm[jk] * pow(rho, n) * Ynm[nm] * pow(2.0, n-j);
                            cnm = pow(-1.0, n) * constants->anm[nm] * constants->anm[jnk] / constants->anm[jk] * pow(rho, n) * constants->Ynm[nm] * powTwo;
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
        else {
            MultipoleMethodFunction::unmorton(nfjc, boxIndex3D);
            boxIndex3D.x = 6 - boxIndex3D.x * 2;
            boxIndex3D.y = 6 - boxIndex3D.y * 2;
            boxIndex3D.z = 6 - boxIndex3D.z * 2;

            rotationMorton(boxIndex3D, je_forward, je_backward);        

            // Initialize MnmVectorA with Mnm values from lower level cell
            for (j = 0; j < numCoefficients; j++) {
                MnmVectorA[j] = variables->Mnm[jb][j];
            }

            rotation(MnmVectorA, MnmVectorB, constants->Dnm[je_forward]);

            // Scale size of parent cell to be equal to 1.
            // Children cell will have side length of 1/2.
            // rho is radius of smallest sphere that can cover a child cell


            // Above is incorrect, rho is distance from center of child to center of its parent
            // This distance is normalized such that length of parent is equal to one
    //        rho = sqrt(3.0)/4;
            rho = sqrt(3.0)/2;

            for (j = 0; j < numExpansions; j++) {
                powTwo = pow(2.0, -j);
                for (k = 0; k <= j; k++) {
                    jk = j * j + j + k;
                    jks = j * (j + 1) / 2 + k;
                    MnmScalar = 0;
                    for (n = 0; n <= j - abs(k); n++) {
                        jnk = (j - n) * (j - n) + j - n + k;
                        jnks = (j - n) * (j - n + 1) / 2 + k;
                        nm = n * n + n;
    //                    cnm = pow(-1.0, n) * anm[nm] * anm[jnk] / anm[jk] * pow(rho, n) * Ynm[nm] * pow(2.0, n-j);
                        cnm = pow(-1.0, n) * constants->anm[nm] * constants->anm[jnk] / constants->anm[jk] * pow(rho, n) * constants->Ynm[nm] * powTwo;
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
}

void* MultipoleMethodForceCalculation::m2l_threaded_rotations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->m2l_kernel_rotations(dfc.threadNum);
    return NULL;
}

void MultipoleMethodForceCalculation::m2l_kernel_rotations(int threadNum){
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
    int jdiff = (1 << numLevel);
    int domainNum_i, domainNum_j;
    int je_forward, je_backward;
    int gcd;
    std::complex<double> mr;
    int tempQuad, m;
    
    boxSize = rootBoxSize / (1 << numLevel);
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
        MultipoleMethodFunction::unmorton(variables->domainBoxIndexFull[ib], boxIndex3D);
        
        ix = boxIndex3D.x + domainBoxIndex3D.x * jdiff;
        iy = boxIndex3D.y + domainBoxIndex3D.y * jdiff;
        iz = boxIndex3D.z + domainBoxIndex3D.z * jdiff;
        
        for (ij = 0; ij < variables->numInteraction[ii]; ij++) {
            
            jj = variables->interactionList[ii][ij][0];
            tempQuad = variables->interactionList[ii][ij][1];
            jb = jj + tempLevelOffset;
            
//            // Initialize MnmVectorB
//            for (j = 0; j < numCoefficients; j++){
//                MnmVectorB[j] = variables->Mnm[jb][j];
//            }
            for (j = 0; j < numExpansions; j++){
                for (m = 0; m <= j; m++){
                    mr = qrf_complex(-m * tempQuad);
                    MnmVectorB[j * (j + 1) / 2 + m] = variables->Mnm[jb][j * (j + 1) / 2 + m] * mr;
                }
            }
            
            // Find global domain of jj-th cell
            domainNum_j = findDomain(jj, numLevel);
            
            // Find global domain 3D indices
            MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][domainNum_j].mortonIndex, domainBoxIndex3D);
            domainBoxIndex3D.x -= variables->offsetSC[0].x;
            domainBoxIndex3D.y -= variables->offsetSC[0].y;
            domainBoxIndex3D.z -= variables->offsetSC[0].z;
            
            // Find local domain 3D indices
            MultipoleMethodFunction::unmorton(variables->domainBoxIndexFull[jb], boxIndex3D);
            
//            jx = boxIndex3D.x + domainBoxIndex3D.x * jdiff;
//            jy = boxIndex3D.y + domainBoxIndex3D.y * jdiff;
            switch (tempQuad){
                case 0:
                    jx = boxIndex3D.x;
                    jy = boxIndex3D.y;
                    break;
                case 1:
                    jx = jdiff - boxIndex3D.y - 1; //-(boxIndex3D.y - jdiff + 1);
                    jy = boxIndex3D.x;
                    break;
                case 2:
                    jx = jdiff - boxIndex3D.x - 1;
                    jy = jdiff - boxIndex3D.y - 1;
                    break;
                case 3:
                    jx = boxIndex3D.y;
                    jy = jdiff - boxIndex3D.x - 1;
                    break;
            }
            
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

void* MultipoleMethodForceCalculation::m2p_threaded_rotations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->m2p_kernel_rotations(dfc.threadNum);
    return NULL;
}

void MultipoleMethodForceCalculation::m2p_kernel_rotations(int threadNum){
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
    
    int tempQuad, cosk, sink;
    std::complex<double> mr;
    
    boxSize = rootBoxSize / (1 << numLevel);
    for (ii = range.x; ii < range.y; ii++) {
        for (ij = 0; ij < variables->numInteraction[ii]; ij++) {
            jj = variables->interactionList[ii][ij][0];
            tempQuad = variables->interactionList[ii][ij][1];
            jb = jj + tempLevelOffset;
            
            cosk = quadrantRotationFunction(tempQuad + 1);
            sink = quadrantRotationFunction(tempQuad);

            for (j = 0; j < numExpansions; j++){
                for (m = 0; m <= j; m++){
                    mr = qrf_complex(-m * tempQuad);
                    MnmVector[j * (j + 1) / 2 + m] = variables->Mnm[jb][j * (j + 1) / 2 + m] * mr;
                }
            }
//            for (j = 0; j < numCoefficients; j++) MnmVector[j] = variables->Mnm[jb][j];
            
//            for (j = 0; j < numCoefficients; j++){
//                printf("variables->Mnm[%d][%d] = %.9g + i%.9g\n",jb,j,std::real(variables->Mnm[jb][j]),std::imag(variables->Mnm[jb][j]));
//            }

            domainNum = findDomain(jj, numLevel);
            domainMorton3 = variables->domainLinkage[0][domainNum].mortonIndex;
            MultipoleMethodFunction::unmorton3(domainMorton3, domainBoxIndex3D);
            domainBoxIndex3D.x -= variables->offsetSC[0].x;
            domainBoxIndex3D.y -= variables->offsetSC[0].y;
            domainBoxIndex3D.z -= variables->offsetSC[0].z;

            MultipoleMethodFunction::unmorton(variables->domainBoxIndexFull[jb], boxIndex3D);

            domainBoxMin.x = boxMin.x + domainBoxIndex3D.x * rootBoxSize;
            domainBoxMin.y = boxMin.y + domainBoxIndex3D.y * rootBoxSize;
            domainBoxMin.z = boxMin.z + domainBoxIndex3D.z * rootBoxSize;
            boxCenter.x = domainBoxMin.x + (boxIndex3D.x + 0.5) * boxSize;
            boxCenter.y = domainBoxMin.y + (boxIndex3D.y + 0.5) * boxSize;
            boxCenter.z = domainBoxMin.z + (boxIndex3D.z + 0.5) * boxSize;
            for (i = variables->particleOffset[0][ii]; i <= variables->particleOffset[1][ii]; i++) {
                if (particleSystem->chargedParticles.acceleration[i].w == inboundCheck){
                    zaligned = 0;
                    dist.x = particleSystem->chargedParticles.position[i].x - (cosk * boxCenter.x - sink * boxCenter.y);
                    dist.y = particleSystem->chargedParticles.position[i].y - (sink * boxCenter.x + cosk * boxCenter.y);
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

void* MultipoleMethodForceCalculation::p2outlier_p_threaded_rotations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->p2outlier_p_kernel_rotations(dfc.threadNum);
    return NULL;
}

void MultipoleMethodForceCalculation::p2outlier_p_kernel_rotations(int threadNum){
    int ii, ij, jj, tempQuad;
        
    vec4<int> r = rangevalues[threadNum];
    
    vec4<int> particleBounds;
    
    for (ii = r.x; ii < r.y; ii++){
        particleBounds.x = variables->outlierParticleOffset[0][ii];
        particleBounds.y = variables->outlierParticleOffset[1][ii] + 1;
        for (ij = 0; ij < variables->numInteraction[ii]; ij++){
            jj = variables->interactionList[ii][ij][0];
            tempQuad = variables->interactionList[ii][ij][1];
            particleBounds.z = variables->particleOffset[0][jj];
            particleBounds.w = variables->particleOffset[1][jj] + 1;
            
            if (vectorizedToggle){
                DirectForceFunction::direct_vectorized_rotational(tempQuad, particleBounds,
                        particleSystem->chargedParticles, jptc);
            }
            else {
                DirectForceFunction::direct_rotational(tempQuad, particleBounds,
                        particleSystem->chargedParticles);
            }
        }
    }
}

void* MultipoleMethodForceCalculation::m2outlier_p_threaded_rotations(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->m2outlier_p_kernel_rotations(dfc.threadNum);
    return NULL;
}

void MultipoleMethodForceCalculation::m2outlier_p_kernel_rotations(int threadNum){
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
    int tempQuad, cosk, sink;
    std::complex<double> mr;
    
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    
    boxSize = rootBoxSize / (1 << numLevel);
    for (ii = range.x; ii < range.y; ii++) {
        for (ij = 0; ij < variables->numInteraction[ii]; ij++) {
            jj = variables->interactionList[ii][ij][0];
            tempQuad = variables->interactionList[ii][ij][1];
            jb = jj + tempLevelOffset;

            cosk = quadrantRotationFunction(tempQuad + 1);
            sink = quadrantRotationFunction(tempQuad);

            for (j = 0; j < numExpansions; j++){
                for (m = 0; m <= j; m++){
                    mr = qrf_complex(-m * tempQuad);
                    MnmVector[j * (j + 1) / 2 + m] = variables->Mnm[jb][j * (j + 1) / 2 + m] * mr;
                }
            }
//            for (j = 0; j < numCoefficients; j++) MnmVector[j] = variables->Mnm[jb][j];

            domainNum = findDomain(jj, numLevel);
            domainMorton3 = variables->domainLinkage[0][domainNum].mortonIndex;
            MultipoleMethodFunction::unmorton3(domainMorton3, domainBoxIndex3D);
            domainBoxIndex3D.x -= variables->offsetSC[0].x;
            domainBoxIndex3D.y -= variables->offsetSC[0].y;
            domainBoxIndex3D.z -= variables->offsetSC[0].z;

            MultipoleMethodFunction::unmorton(variables->domainBoxIndexFull[jb], boxIndex3D);

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
                    dist.x = particleSystem->chargedParticles.position[i].x - (cosk * boxCenter.x - sink * boxCenter.y);
                    dist.y = particleSystem->chargedParticles.position[i].y - (sink * boxCenter.x + cosk * boxCenter.y);
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