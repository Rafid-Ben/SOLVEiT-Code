/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "force_calculation_methods.h"
#include "direct_force_method.h"

void MultipoleMethodForceCalculation::rotation(std::complex<double>* Cnm,
        std::complex<double>* CnmOut, std::complex<double>** Dnm){
    int n, m, nms, k, nk, nks;
    std::complex<double> CnmScalar;

    for (n = 0; n < numExpansions; n++) {
        for (m = 0; m <= n; m++) {
            nms = n * (n + 1) / 2 + m;
            CnmScalar = 0;
            for (k = -n; k <= -1; k++) {
                nk = n * (n + 1) + k;
                nks = n * (n + 1) / 2 - k;
                CnmScalar += Dnm[m][nk] * conj(Cnm[nks]);
            }
            for (k = 0; k <= n; k++) {
                nk = n * (n + 1) + k;
                nks = n * (n + 1) / 2 + k;
                CnmScalar += Dnm[m][nk] * Cnm[nks];
            }
            CnmOut[nms] = CnmScalar;
        }
    }
}

void MultipoleMethodForceCalculation::direct(){
    timing.timer.startTimer();
    
    int numParticles = lastParticle - firstParticle;
    int k;
#ifdef _PTHREAD_H
    for (k = 0; k < numThreads; k++){
        dF[k].d = this;
        dF[k].threadNum = k;
        rangevalues[k].x = (k * numParticles) / numThreads + firstParticle;
        rangevalues[k].y = ((k + 1) * numParticles) / numThreads + firstParticle;
//        rangevalues[k].z = firstParticle;
//        rangevalues[k].w = lastParticle;
        switch (particleSystem->params.symType){
            case (SymmetryType::none):
                pthread_create(&threads[k], nullptr, direct_threaded, &dF[k]);
                break;
            case (SymmetryType::rotational):
                pthread_create(&threads[k], nullptr, direct_threaded_rotations, &dF[k]);
                break;
            case (SymmetryType::translational):
                pthread_create(&threads[k], nullptr, direct_threaded_translations, &dF[k]);
                break;
        }
    }
    
    for (k = 0; k < numThreads; k++){
        pthread_join(threads[k], nullptr);
    }
#else
    rangevalues[0].x = firstParticle;
    rangevalues[0].y = lastParticle;
    rangevalues[0].z = firstParticle;
    rangevalues[0].w = lastParticle;
    switch (particleSystem->params.symType){
        case (SymmetryType::none):
            direct_kernel();
            break;
        case (SymmetryType::rotational):
            direct_kernel_rotations();
            break;
        case (SymmetryType::translational):
            direct_kernel_translations();
            break;
    }
#endif    
    
    timing.direct_time += timing.timer.getTimer();
    timing.timer.resetTimer();
}

void* MultipoleMethodForceCalculation::direct_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->direct_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::direct_kernel(int threadNum){
    
    if (vectorizedToggle){
        DirectForceFunction::direct_vectorized(rangevalues[threadNum],
                particleSystem->chargedParticles, jptc);
    }
    else {
        DirectForceFunction::direct(rangevalues[threadNum],
                particleSystem->chargedParticles);
    }
}

void MultipoleMethodForceCalculation::p2p(int numBoxIndex, int toggle){
    timing.timer.startTimer();
    int k;
    
#ifdef _PTHREAD_H
    for (k = 0; k < numThreads; k++){
        dF[k].d = this;
        dF[k].threadNum = k;
        rangevalues[k].x = (k * numBoxIndex) / numThreads;
        rangevalues[k].y = ((k + 1) * numBoxIndex) / numThreads;
        rangevalues[k].z = numBoxIndex;
        rangevalues[k].w = 0;
        if (toggle == 0){
            switch (particleSystem->params.symType){
                case (SymmetryType::none):
                    pthread_create(&threads[k], nullptr, p2p_threaded, &dF[k]);
                    break;
                case (SymmetryType::rotational):
                    pthread_create(&threads[k], nullptr, p2p_threaded_rotations, &dF[k]);
                    break;
                case (SymmetryType::translational):
                    pthread_create(&threads[k], nullptr, p2p_threaded_translations, &dF[k]);
                    break;
            }
        }
        else if (toggle == 1){
            switch (particleSystem->params.symType){
                case (SymmetryType::none):
                    pthread_create(&threads[k], nullptr, sp2p_threaded, &dF[k]);
                    break;
                case (SymmetryType::rotational):
                    pthread_create(&threads[k], nullptr, sp2p_threaded_rotations, &dF[k]);
                    break;
                case (SymmetryType::translational):
                    pthread_create(&threads[k], nullptr, sp2p_threaded_translations, &dF[k]);
                    break;
            }
        }
    }
    
    for (int k = 0; k < numThreads; k++){
        pthread_join(threads[k], nullptr);
    }
#else
    rangevalues[0].x = 0;
    rangevalues[0].y = numBoxIndex;
    rangevalues[0].z = numBoxIndex;
    rangevalues[0].w = 0;
    if (toggle == 0){
        switch (particleSystem->params.symType){
            case (SymmetryType::none):
                p2p_kernel();
                break;
            case (SymmetryType::rotational):
                p2p_kernel_rotations();
                break;
            case (SymmetryType::translational):
                p2p_kernel_translations();
                break;
        }
    }
    else if (toggle == 1){
        switch (particleSystem->params.symType){
            case (SymmetryType::none):
                sp2p_kernel();
                break;
            case (SymmetryType::rotational):
                sp2p_kernel_rotations();
                break;
            case (SymmetryType::translational):
                sp2p_kernel_translations();
                break;
        }
    }
#endif   
    
    if (toggle == 0){
        timing.p2p_time += timing.timer.getTimer();
    }
    else if (toggle == 1){
        timing.sp2p_time += timing.timer.getTimer();
    }
    timing.timer.resetTimer();
}

void* MultipoleMethodForceCalculation::p2p_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->p2p_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::p2p_kernel(int threadNum){
    int ii, ij, jj;
    
    vec4<int> r = rangevalues[threadNum];
    
    vec4<int> particleBounds;
    
    for (ii = r.x; ii < r.y; ii++){
        particleBounds.x = variables->particleOffset[0][ii];
        particleBounds.y = variables->particleOffset[1][ii] + 1;
        for (ij = 0; ij < variables->numInteraction[ii]; ij++){
            jj = variables->interactionList[ii][ij][0];
            particleBounds.z = variables->particleOffset[0][jj];
            particleBounds.w = variables->particleOffset[1][jj] + 1;
            
            if (vectorizedToggle){
                DirectForceFunction::direct_vectorized(particleBounds,
                        particleSystem->chargedParticles, jptc);
            }
            else {
                DirectForceFunction::direct(particleBounds,
                        particleSystem->chargedParticles);
            }
        }
    }
}

void MultipoleMethodForceCalculation::p2m(int numBoxIndex, int toggle){
    timing.timer.startTimer();
    int k;
#ifdef _PTHREAD_H
    for (k = 0; k < numThreads; k++){
        dF[k].d = this;
        dF[k].threadNum = k;
        rangevalues[k].x = (k * numBoxIndex) / numThreads;
        rangevalues[k].y = ((k + 1) * numBoxIndex) / numThreads;
        rangevalues[k].z = numBoxIndex;
        rangevalues[k].w = 0;
        if (toggle == 0){
            if (particleSystem->params.symType == SymmetryType::translational){
                pthread_create(&threads[k], nullptr, p2m_threaded_translations, &dF[k]);
            }
            else {
                pthread_create(&threads[k], nullptr, p2m_threaded, &dF[k]);
            }
        }
        else if (toggle == 1){
            pthread_create(&threads[k], nullptr, sp2m_threaded, &dF[k]);
        }
    }
    
    for (k = 0; k < numThreads; k++){
        pthread_join(threads[k], nullptr);
    }
#else
    rangevalues[0].x = 0;
    rangevalues[0].y = numBoxIndex;
    rangevalues[0].z = numBoxIndex;
    rangevalues[0].w = 0;
    if (toggle == 0){
        if (particleSystem->params.symType == SymmetryType::translational){
            p2m_kernel_translations();
        }
        else {
            p2m_kernel();
        }
    }
    else if (toggle == 1){
        sp2m_kernel();
    }
#endif
    
    if (toggle == 0){
        timing.p2m_time += timing.timer.getTimer();
    }
    else if (toggle == 1){
        timing.sp2m_time += timing.timer.getTimer();
    }
    timing.timer.resetTimer();
}

void* MultipoleMethodForceCalculation::p2m_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->p2m_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::p2m_kernel(int threadNum){
    int jj, j, n, m, nm, nms;
    vec3<int> boxIndex3D;
    vec3<int> domainBoxIndex3D;
    vec3<float> boxCenter, domainBoxMin;
    vec3<double> dist;
    double boxSize, rho, alpha, beta;
    double xx, s2, fact, pn, p, p1, p2;
    double rhom, rhon;
    double YnmReal[FMMConstant::numExpansion2];
    std::complex<double> MnmVector[FMMConstant::numCoefficients], I(0.0, 1.0), eim;
        
    boxSize = rootBoxSize / (1 << maxLevel);
    
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
        MultipoleMethodFunction::unmorton(variables->domainBoxIndexFull[jj], boxIndex3D);
        
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
            cart2sph_exact(rho, alpha, beta, dist.x, dist.y, dist.z);
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

void MultipoleMethodForceCalculation::m2m(int numBoxIndex, int numBoxIndexOld, int numLevel, int toggle){
    timing.timer.startTimer();
    
    int ii, ib, j;
    
    if (toggle == 0){
        // Use when considering subcell starting domain
        
        if (numLevel > 0){
            for (ii = 0; ii < numBoxIndex; ii++) {
                ib = ii + variables->levelOffset[numLevel - 1];
                for (j = 0; j < numCoefficients; j++) {
                    variables->Mnm[ib][j] = 0;
                }
            }
        }
        else{
            for (ii = 0; ii < numBoxIndex; ii++){
                ib = ii + variables->levelOffset[maxLevel];
                for (j = 0; j < numCoefficients; j++){
                    variables->Mnm[ib][j] = 0;
                }
            }
        }
        
        int k;
#ifdef _PTHREAD_H
        for (k = 0; k < numThreads; k++){
            dF[k].d = this;
            dF[k].threadNum = k;
            rangevalues[k].x = (k * numBoxIndexOld) / numThreads;
            rangevalues[k].y = ((k + 1) * numBoxIndexOld) / numThreads;
            rangevalues[k].z = numBoxIndexOld;
            rangevalues[k].w = numLevel;
            switch (particleSystem->params.symType){
                case (SymmetryType::none):
                    pthread_create(&threads[k], nullptr, m2m_threaded, &dF[k]);
                    break;
                case (SymmetryType::rotational):
                    pthread_create(&threads[k], nullptr, m2m_threaded_rotations, &dF[k]);
                    break;
                case (SymmetryType::translational):
                    pthread_create(&threads[k], nullptr, m2m_threaded_translations, &dF[k]);
                    break;
            }
        }

        for (int k = 0; k < numThreads; k++){
            pthread_join(threads[k], nullptr);
        }
#else
        rangevalues[0].x = 0;
        rangevalues[0].y = numBoxIndexOld;
        rangevalues[0].z = numBoxIndexOld;
        rangevalues[0].w = numLevel;
        switch (particleSystem->params.symType){
            case (SymmetryType::none):
                m2m_kernel();
                break;
            case (SymmetryType::rotational):
                m2m_kernel_rotations();
                break;
            case (SymmetryType::translational):
                m2m_kernel_translations();
                break;
        }
#endif
    }
    else if (toggle == 1){
        // Used when considering supercell starting domain
        
        for (ii = 0; ii < numBoxIndex; ii++){
            ib = ii + variables->levelOffset[numLevel + maxLevel];
            for (j = 0; j < numCoefficients; j++){
                variables->Mnm[ib][j] = 0;
            }
        }
        
#ifdef _PTHREAD_H
        int k;
        for (k = 0; k < numThreads; k++){
            dF[k].d = this;
            dF[k].threadNum = k;
            rangevalues[k].x = (k * numBoxIndexOld) / numThreads;
            rangevalues[k].y = ((k + 1) * numBoxIndexOld) / numThreads;
            rangevalues[k].z = numBoxIndexOld;
            rangevalues[k].w = numLevel;
            if (particleSystem->params.symType == SymmetryType::translational){
                pthread_create(&threads[k], nullptr, sm2m_threaded_translations, &dF[k]);
            }
            else {
                pthread_create(&threads[k], nullptr, sm2m_threaded, &dF[k]);
            }
        }

        for (int k = 0; k < numThreads; k++){
            pthread_join(threads[k], nullptr);
        }
#else
        rangevalues[0].x = 0;
        rangevalues[0].y = numBoxIndexOld;
        rangevalues[0].z = numBoxIndexOld;
        rangevalues[0].w = numLevel;
        if (particleSystem->params.symType == SymmetryType::translational){
            sm2m_kernel_translations();
        }
        else {
            sm2m_kernel();
        }
#endif
    }
    
    if (toggle == 0){
        timing.m2m_time += timing.timer.getTimer();
    }
    else if (toggle == 1){
        timing.sm2m_time += timing.timer.getTimer();
    }
    timing.timer.resetTimer();
}

void* MultipoleMethodForceCalculation::m2m_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->m2m_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::m2m_kernel(int threadNum){
    int ib, j, jj, nfjp, nfjc, jb, k, jk, jks, n, jnk, jnks, nm;
    int je_forward, je_backward;
    vec3<int> boxIndex3D;
    int boxIndexFullTemp;
    double rho, powTwo;
    std::complex<double> cnm, MnmScalar;
    std::complex<double> MnmVectorB[FMMConstant::numCoefficients], MnmVectorA[FMMConstant::numCoefficients];
        
    vec4<int> r = rangevalues[threadNum];
    int numLevel = r.w;
    
    int domainNum, domainBias;
    int domainNumBoxIndexTemp = 1 << (3 * numLevel);
    
    int tempLevelOffset = (numLevel == 0) ? variables->levelOffset[maxLevel] : variables->levelOffset[numLevel - 1];
    
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

void MultipoleMethodForceCalculation::m2l(int numBoxIndex, int numLevel, int toggle, int topLevel){
    timing.timer.startTimer();
    
    int i, j, jj, jb;
    
    if (topLevel == 1){
        for (i = 0; i < numBoxIndex; i++){
            for (j = 0; j < numCoefficients; j++){
                variables->Lnm[i][j] = 0;
            }
        }
    }
    
    if (toggle == 0){
        int tempLevelOffset = (numLevel == 0) ? variables->levelOffset[maxLevel] : variables->levelOffset[numLevel - 1];
        int k;
#ifdef _PTHREAD_H
        for (k = 0; k < numThreads; k++){
            dF[k].d = this;
            dF[k].threadNum = k;
            rangevalues[k].x = (k * numBoxIndex) / numThreads;
            rangevalues[k].y = ((k + 1) * numBoxIndex) / numThreads;
            rangevalues[k].z = numBoxIndex;
            rangevalues[k].w = numLevel;
            
            switch (particleSystem->params.symType){
                case (SymmetryType::none):
                    pthread_create(&threads[k], nullptr, m2l_threaded, &dF[k]);
                    break;
                case (SymmetryType::rotational):
                    pthread_create(&threads[k], nullptr, m2l_threaded_rotations, &dF[k]);
                    break;
                case (SymmetryType::translational):
                    pthread_create(&threads[k], nullptr, m2l_threaded_translations, &dF[k]);
                    break;
            }
        }

        for (k = 0; k < numThreads; k++){
            pthread_join(threads[k], nullptr);
        }
#else
        rangevalues[0].x = 0;
        rangevalues[0].y = numBoxIndex;
        rangevalues[0].z = numBoxIndex;
        rangevalues[0].w = numLevel;
        switch (particleSystem->params.symType){
            case (SymmetryType::none):
                m2l_kernel();
                break;
            case (SymmetryType::rotational):
                m2l_kernel_rotations();
                break;
            case (SymmetryType::translational):
                m2l_kernel_translations();
                break;
        }
#endif
        for (jj = 0; jj < numBoxIndex; jj++){
            jb = jj + tempLevelOffset;
            for (j = 0; j < numCoefficients; j++) {
                variables->Mnm[jb][j] = 0;
            }
        }
    }
    else if (toggle == 1){
#ifdef _PTHREAD_H
        int k;
        for (k = 0; k < numThreads; k++){
            dF[k].d = this;
            dF[k].threadNum = k;
            rangevalues[k].x = (k * numBoxIndex) / numThreads;
            rangevalues[k].y = ((k + 1) * numBoxIndex) / numThreads;
            rangevalues[k].z = numBoxIndex;
            rangevalues[k].w = numLevel;
            if (particleSystem->params.symType == SymmetryType::translational){
                pthread_create(&threads[k], nullptr, sm2l_threaded_translations, &dF[k]);
            }
            else {
                pthread_create(&threads[k], nullptr, sm2l_threaded, &dF[k]);
            }
        }

        for (k = 0; k < numThreads; k++){
            pthread_join(threads[k], nullptr);
        }
#else
        rangevalues[0].x = 0;
        rangevalues[0].y = numBoxIndex;
        rangevalues[0].z = numBoxIndex;
        rangevalues[0].w = numLevel;
        if (particleSystem->params.symType == SymmetryType::translational){
            sm2l_kernel_translations();
        }
        else {
            sm2l_kernel();
        }
#endif
        for (jj = 0; jj < numBoxIndex; jj++) {
            jb = jj + variables->levelOffset[numLevel + maxLevel];
            for (j = 0; j < numCoefficients; j++) {
                variables->Mnm[jb][j] = 0;
            }
        }
    }
    
    if (toggle == 0){
        timing.m2l_time += timing.timer.getTimer();
    }
    else if (toggle == 1){
        timing.sm2l_time += timing.timer.getTimer();
    }
    timing.timer.resetTimer();
}

void* MultipoleMethodForceCalculation::m2l_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->m2l_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::m2l_kernel(int threadNum){
    int j, ii, ib, ix, iy, iz, ij, jj, jb, jx, jy, jz, k, jk, jks, n, nk, nks, jkn, jnk;
    vec3<int> boxIndex3D;
    vec3<double> dist;
    double boxSize, rho, rhoj, rhojk, rhojn;
    std::complex<double> LnmVectorA[FMMConstant::numCoefficients], MnmVectorA[FMMConstant::numCoefficients];
    std::complex<double> LnmVectorB[FMMConstant::numCoefficients], MnmVectorB[FMMConstant::numCoefficients];
    std::complex<double> cnm, LnmScalar;
    
    vec3<int> domainBoxIndex3D;
    
    vec4<int> r = rangevalues[threadNum];
    int numLevel = r.w;
    int jdiff = (1 << numLevel);
    int domainNum_i, domainNum_j;
    int je_forward, je_backward;
    int gcd;
    
    boxSize = rootBoxSize / (1 << numLevel);
    double boxSize2 = boxSize * boxSize;
    
    int tempLevelOffset = (numLevel == 0) ? variables->levelOffset[maxLevel] : variables->levelOffset[numLevel - 1];
//    int numExpansion2 = numExpansions * numExpansions;
    
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
            
            // Find local domain 3D indices
            MultipoleMethodFunction::unmorton(variables->domainBoxIndexFull[jb], boxIndex3D);
            
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
                }
            }            
        }
    }
}

void MultipoleMethodForceCalculation::l2l(int numBoxIndex, int numLevel, int toggle){
    timing.timer.startTimer();
    
    int numBoxIndexOld, ii, i, ib, domainNum, domainBias, nbc;
    int numBoxesSubcell;
    int domainNumBoxIndexTemp;
    if (particleSystem->params.symType == SymmetryType::translational){
        numBoxesSubcell = 27;
        domainNumBoxIndexTemp = pow(3, 3 * numLevel);
    }
    else {
        numBoxesSubcell = 8;
        domainNumBoxIndexTemp = 1 << (3 * numLevel);
    }

    numBoxIndexOld = numBoxIndex;
    for (ii = 0; ii < numBoxIndexOld; ii++) {
        for (i = 0; i < numCoefficients; i++) {
            variables->LnmOld[ii][i] = variables->Lnm[ii][i];
        }
    }
    numBoxIndexOld = 0;
    
    if (toggle == 0){
        int tempLevelOffset = (numLevel == 0) ? variables->levelOffset[maxLevel] : variables->levelOffset[numLevel - 1];
//        int domainNumBoxIndexTemp = 1 << (3 * numLevel);
        nbc = -1;
//        domainBias = 0;
        domainNum = findDomain(0, numLevel);
        domainBias = domainNum * domainNumBoxIndexTemp;
        
        for (i = 0; i < numBoxIndexFull; i++) variables->neo[i] = -1;
        
        for (ii = 0; ii < numBoxIndex; ii++){
            ib = ii + tempLevelOffset;

            // If ii is in a different domain to ii - 1, find it's new domain
            if (ii > variables->domainOffset[1][domainNum][numLevel]){
                // Increment the domain number until find a domain containing boxes with particles inside
                // Since ii will only increase by 1 between iterations, the next domain containing particles contains ii
                do{
                    domainNum++;
                }
                while (variables->domainOffset[1][domainNum][numLevel] < variables->domainOffset[0][domainNum][numLevel]);
                domainBias = domainNum * domainNumBoxIndexTemp;
                nbc = -1;
            }

            if (nbc != variables->domainBoxIndexFull[ib] / numBoxesSubcell){
                nbc = variables->domainBoxIndexFull[ib] / numBoxesSubcell;
                variables->neo[nbc + domainBias] = numBoxIndexOld;
                numBoxIndexOld++;
            }
        }
        int k;
#ifdef _PTHREAD_H
        for (k = 0; k < numThreads; k++){
            dF[k].d = this;
            dF[k].threadNum = k;
            rangevalues[k].x = (k * numBoxIndex) / numThreads;
            rangevalues[k].y = ((k + 1) * numBoxIndex) / numThreads;
            rangevalues[k].z = numBoxIndexOld;
            rangevalues[k].w = numLevel;
            if (particleSystem->params.symType == SymmetryType::translational){
                pthread_create(&threads[k], nullptr, l2l_threaded_translations, &dF[k]);
            }
            else{
                pthread_create(&threads[k], nullptr, l2l_threaded, &dF[k]);
            }
        }

        for (k = 0; k < numThreads; k++){
            pthread_join(threads[k], nullptr);
        }
#else
        rangevalues[0].x = 0;
        rangevalues[0].y = numBoxIndex;
        rangevalues[0].z = numBoxIndexOld;
        rangevalues[0].w = numLevel;
        if (particleSystem->params.symType == SymmetryType::translational){
            l2l_kernel_translations();
        }
        else{
            l2l_kernel();
        }
#endif
    }
    else if (toggle == 1){
        nbc = -1;
        for (ii = 0; ii < numBoxIndex; ii++){
            ib = ii + variables->levelOffset[numLevel + maxLevel];

//            if (nbc != variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].parent){
//                nbc = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].parent;
            if (nbc != variables->domainBoxIndexMask[variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].parent]){
                nbc = variables->domainBoxIndexMask[variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].parent];
                variables->neo[nbc] = numBoxIndexOld;
                numBoxIndexOld++;
            }
        }
        int k;
#ifdef _PTHREAD_H
        for (k = 0; k < numThreads; k++){
            dF[k].d = this;
            dF[k].threadNum = k;
            rangevalues[k].x = (k * numBoxIndex) / numThreads;
            rangevalues[k].y = ((k + 1) * numBoxIndex) / numThreads;
            rangevalues[k].z = numBoxIndexOld;
            rangevalues[k].w = numLevel;
            pthread_create(&threads[k], nullptr, sl2l_threaded, &dF[k]);
        }

        for (k = 0; k < numThreads; k++){
            pthread_join(threads[k], nullptr);
        }
#else
        rangevalues[0].x = 0;
        rangevalues[0].y = numBoxIndex;
        rangevalues[0].z = numBoxIndexOld;
        rangevalues[0].w = numLevel;
        sl2l_kernel();
#endif
    }
    
    if (toggle == 0){
        timing.l2l_time += timing.timer.getTimer();
    }
    else if (toggle == 1){
        timing.sl2l_time += timing.timer.getTimer();
    }
    timing.timer.resetTimer();
}

void* MultipoleMethodForceCalculation::l2l_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->l2l_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::l2l_kernel(int threadNum){
    int ii, ib, i, nfip, nfic, j, k, jk, jks, n, jnk, nk, nks;
    vec3<int> boxIndex3D;
    double rho;
    std::complex<double> cnm, LnmScalar;
    std::complex<double> LnmVectorA[FMMConstant::numCoefficients], LnmVectorB[FMMConstant::numCoefficients];
    
    vec4<int> r = rangevalues[threadNum];
    int numLevel = r.w;

    int tempLevelOffset = (numLevel == 0) ? variables->levelOffset[maxLevel] : variables->levelOffset[numLevel - 1];
    
    int domainNum, domainBias;
    int domainNumBoxIndexTemp = 1 << (3 * numLevel);
    
    int je_forward, je_backward;
    
    for (ii = r.x; ii < r.y; ii++) {        
        ib = ii + tempLevelOffset;
        
        nfip = variables->domainBoxIndexFull[ib] / 8;
        nfic = variables->domainBoxIndexFull[ib] % 8;
        MultipoleMethodFunction::unmorton(nfic, boxIndex3D);
        
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
        
        boxIndex3D.x = boxIndex3D.x * 2 + 4;
        boxIndex3D.y = boxIndex3D.y * 2 + 4;
        boxIndex3D.z = boxIndex3D.z * 2 + 4;
        
        rotationMorton(boxIndex3D, je_forward, je_backward);
        
        ib = variables->neo[nfip + domainBias];
        
        for (i = 0; i < numCoefficients; i++) {
            LnmVectorA[i] = variables->LnmOld[ib][i];
        }

        rotation(LnmVectorA, LnmVectorB, constants->Dnm[je_forward]);

        rho = sqrt(3.0) / 2;
        for (j = 0; j < numExpansions; j++){
            for (k = 0; k <= j; k++){
                jk = j * j + j + k;
                jks = j*(j+1)/2+k;
                LnmScalar = 0;
                for (n = j; n < numExpansions; n++){
                    jnk = (n-j)*(n-j)+n-j;
                    nk = n*n+n+k;
                    nks = n*(n+1)/2+k;
                    cnm = constants->anm[jnk] * constants->anm[jk] / constants->anm[nk] * pow(rho, n-j) * constants->Ynm[jnk] * pow(2.0, -n-1);
                    LnmScalar += LnmVectorB[nks] * cnm;
                }
                LnmVectorA[jks] = LnmScalar;
            }
        }
        
        rotation(LnmVectorA, LnmVectorB, constants->Dnm[je_backward]);
        for (i = 0; i < numCoefficients; i++) {
            variables->Lnm[ii][i] = LnmVectorB[i];
//            printf("variables->Lnm[%d][%d] = %.9g + i%.9g\n",ii,i,std::real(variables->Lnm[ii][i]),std::imag(variables->Lnm[ii][i]));
        }
    }
}

void MultipoleMethodForceCalculation::l2p(int numBoxIndex, int toggle){
    timing.timer.startTimer();
    int k;
#ifdef _PTHREAD_H
        for (k = 0; k < numThreads; k++){
            dF[k].d = this;
            dF[k].threadNum = k;
            rangevalues[k].x = (k * numBoxIndex) / numThreads;
            rangevalues[k].y = ((k + 1) * numBoxIndex) / numThreads;
            rangevalues[k].z = numBoxIndex;
            rangevalues[k].w = 0;
            if (toggle == 0){
                if (particleSystem->params.symType == SymmetryType::translational){
                    pthread_create(&threads[k], nullptr, l2p_threaded_translations, &dF[k]);
                }
                else {
                    pthread_create(&threads[k], nullptr, l2p_threaded, &dF[k]);
                }
            }
            else if (toggle == 1){
                pthread_create(&threads[k], nullptr, sl2p_threaded, &dF[k]);
            }
        }

        for (int k = 0; k < numThreads; k++){
            pthread_join(threads[k], nullptr);
        }
#else
        rangevalues[0].x = 0;
        rangevalues[0].y = numBoxIndex;
        rangevalues[0].z = numBoxIndex;
        rangevalues[0].w = 0;
        if (toggle == 0){
            if (particleSystem->params.symType == SymmetryType::translational){
                l2p_kernel_translations();
            }
            else {
                l2p_kernel();
            }
        }
        else if (toggle == 1){
            sl2p_kernel();
        }
#endif
        
    if (toggle == 0){
        timing.l2p_time += timing.timer.getTimer();
    }
    else if (toggle == 1){
        timing.sl2p_time += timing.timer.getTimer();
    }
    timing.timer.resetTimer();
}

void* MultipoleMethodForceCalculation::l2p_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->l2p_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::l2p_kernel(int threadNum){
    vec4<int> range = rangevalues[threadNum];
    
    int ii, i, n, nm, nms, m;
    vec3<int> boxIndex3D;
    vec3<int> domainBoxIndex3D;
    vec3<float> boxCenter, domainBoxMin;
    vec3<double> accel, dist;
    double boxSize, r, theta, phi, accelR, accelTheta, accelPhi, accel_const;
    double xx, yy, s2, fact, pn, p, p1, p2, rn;
    double YnmReal[FMMConstant::numExpansion2], YnmRealTheta[FMMConstant::numExpansion2];
    std::complex<double> LnmVector[FMMConstant::numCoefficients];
    std::complex<double> rr, rtheta, rphi, I(0.0, 1.0), eim;
    int domainNum, domainMorton3;
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    
    boxSize = rootBoxSize / (1 << maxLevel);
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
        
        MultipoleMethodFunction::unmorton(variables->domainBoxIndexFull[ii], boxIndex3D);
                
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
                cart2sph_exact(r, theta, phi, dist.x, dist.y, dist.z);

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

void MultipoleMethodForceCalculation::m2p(int numBoxIndex, int numLevel, int toggle){
    timing.timer.startTimer();
    int k;
#ifdef _PTHREAD_H
        for (k = 0; k < numThreads; k++){
            dF[k].d = this;
            dF[k].threadNum = k;
            rangevalues[k].x = (k * numBoxIndex) / numThreads;
            rangevalues[k].y = ((k + 1) * numBoxIndex) / numThreads;
            rangevalues[k].z = numLevel;
            rangevalues[k].w = 0;
            if (toggle == 0){
                switch (particleSystem->params.symType){
                    case (SymmetryType::none):
                        pthread_create(&threads[k], nullptr, m2p_threaded, &dF[k]);
                        break;
                    case (SymmetryType::rotational):
                        pthread_create(&threads[k], nullptr, m2p_threaded_rotations, &dF[k]);
                        break;
                    case (SymmetryType::translational):
                        pthread_create(&threads[k], nullptr, m2p_threaded_translations, &dF[k]);
                        break;
                }
            }
            else if (toggle == 1){
                if (particleSystem->params.symType == SymmetryType::translational){
                    pthread_create(&threads[k], nullptr, sm2p_threaded_translations, &dF[k]);
                }
                else {
                    pthread_create(&threads[k], nullptr, sm2p_threaded, &dF[k]);
                }
            }
        }

        for (int k = 0; k < numThreads; k++){
            pthread_join(threads[k], nullptr);
        }
#else
        rangevalues[0].x = 0;
        rangevalues[0].y = numBoxIndex;
        rangevalues[0].z = numLevel;
        rangevalues[0].w = 0;
        if (toggle == 0){
            switch (particleSystem->params.symType){
                case (SymmetryType::none):
                    m2p_kernel();
                    break;
                case (SymmetryType::rotational):
                    m2p_kernel_rotations();
                    break;
                case (SymmetryType::translational):
                    m2p_kernel_translations();
                    break;
            }
        }
        else if (toggle == 1){
            if (particleSystem->params.symType == SymmetryType::translational){
                sm2p_kernel_translations();
            }
            else {
                sm2p_kernel();
            }
        }
#endif
        
    if (toggle == 0){
        timing.m2p_time += timing.timer.getTimer();
    }
    else if (toggle == 1){
        timing.sm2p_time += timing.timer.getTimer();
    }
    timing.timer.resetTimer();
}

void* MultipoleMethodForceCalculation::m2p_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->m2p_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::m2p_kernel(int threadNum){
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
    std::complex<double> MnmVector[FMMConstant::numCoefficients];
    std::complex<double> rr, rtheta, rphi, I(0.0, 1.0), eim;
    int domainNum, domainMorton3;
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    
    int tempLevelOffset = (numLevel == 0) ? variables->levelOffset[maxLevel] : variables->levelOffset[numLevel - 1];
    int zaligned;
    
    boxSize = rootBoxSize / (1 << numLevel);
    for (ii = range.x; ii < range.y; ii++) {
        for (ij = 0; ij < variables->numInteraction[ii]; ij++) {
            jj = variables->interactionList[ii][ij][0];
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
                    dist.x = particleSystem->chargedParticles.position[i].x - boxCenter.x;
                    dist.y = particleSystem->chargedParticles.position[i].y - boxCenter.y;
                    dist.z = particleSystem->chargedParticles.position[i].z - boxCenter.z;
                    cart2sph_exact(r, theta, phi, dist.x, dist.y, dist.z);
                    r = sqrt(r * r + PhysicalConstant::SOFTENING);
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

void* MultipoleMethodForceCalculation::sm2m_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sm2m_kernel(dfc.threadNum);
    return nullptr;
};
void MultipoleMethodForceCalculation::sm2m_kernel(int threadNum){
    int ib, j, jj, nfjc, jb, k, jk, jks, n, jnk, jnks, nm;
    int je_forward, je_backward;
    vec3<int> boxIndex3D;
    int boxIndexFullTemp;
    double rho, powThree;
    std::complex<double> cnm, MnmScalar;
    std::complex<double> MnmVectorB[FMMConstant::numCoefficients], MnmVectorA[FMMConstant::numCoefficients];
        
    vec4<int> r = rangevalues[threadNum];
    int numLevel = r.w;
    
    int parentIndex;
        
    for (jj = r.x; jj < r.y; jj++) {
        jb = jj + variables->levelOffset[numLevel + maxLevel - 1];
//        if (jj == 9) jb++;
//        if (jj == 10) jb--;
        
        boxIndexFullTemp = variables->domainLinkage[numLevel - 1][variables->domainBoxIndexFull[jb]].mortonIndex;
        nfjc = boxIndexFullTemp % 27;
        parentIndex = variables->domainLinkage[numLevel - 1][variables->domainBoxIndexFull[jb]].parent;

        // Assuming all supercells contain particleSystem->chargedParticles. Good assumption, but could change to boxIndexMask[parentIndex]
//        ib = parentIndex + variables->levelOffset[numLevel + maxLevel];
        ib = variables->domainBoxIndexMask[parentIndex] + variables->levelOffset[numLevel + maxLevel];
        
        MultipoleMethodFunction::unmorton3(nfjc, boxIndex3D);
        
//        rho = sqrt(abs(boxIndex3D.x - 1) + abs(boxIndex3D.y - 1) + abs(boxIndex3D.z - 1)) / 3;
        rho = sqrt(abs(boxIndex3D.x - 1) + abs(boxIndex3D.y - 1) + abs(boxIndex3D.z - 1));
        // More generally,
        // rho = sqrt(pow((boxIndex3D.x-center.x),2) + pow((boxIndex3D.y-center.y),2) + pow((boxIndex3D.z-center.z),2)) / n
        // where considering nxnxn boxes
        
        
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

void* MultipoleMethodForceCalculation::sm2l_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sm2l_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::sm2l_kernel(int threadNum){
    int j, ii, ib, ix, iy, iz, ij, jj, jb, jx, jy, jz, k, jk, jks, n, nk, nks, jkn, jnk;
    vec3<int> boxIndex3D;
    vec3<double> dist;
    double boxSize, rho, rhoj, rhojk, rhojn;
    std::complex<double> LnmVectorA[FMMConstant::numCoefficients], MnmVectorA[FMMConstant::numCoefficients];
    std::complex<double> LnmVectorB[FMMConstant::numCoefficients], MnmVectorB[FMMConstant::numCoefficients];
    std::complex<double> cnm, LnmScalar;
        
    vec4<int> r = rangevalues[threadNum];
    int numLevel = r.w;
    
    boxSize = rootBoxSize * pow(3, numLevel);
    double boxSize2 = boxSize * boxSize;
    int je_forward, je_backward, gcd;
    
//    int numExpansion2 = numExpansions * numExpansions;
    
    int tempMortonIndex;
    
    for (ii = r.x; ii < r.y; ii++) {
        ib = ii + variables->levelOffset[numLevel + maxLevel];
        
        tempMortonIndex = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].mortonIndex;
        
        MultipoleMethodFunction::unmorton3(tempMortonIndex, boxIndex3D);
        ix = boxIndex3D.x;
        iy = boxIndex3D.y;
        iz = boxIndex3D.z;
        for (ij = 0; ij < variables->numInteraction[ii]; ij++) {            
            jj = variables->interactionList[ii][ij][0];
            jb = jj + variables->levelOffset[maxLevel + numLevel];
            tempMortonIndex = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[jb]].mortonIndex;
            
            for (j = 0; j < numCoefficients; j++){
                MnmVectorB[j] = variables->Mnm[jb][j];
            }
            
            MultipoleMethodFunction::unmorton3(tempMortonIndex, boxIndex3D);
            
            jx = boxIndex3D.x;
            jy = boxIndex3D.y;
            jz = boxIndex3D.z;
            
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

void* MultipoleMethodForceCalculation::sl2l_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sl2l_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::sl2l_kernel(int threadNum){
    int ii, ib, i, nfip, nfic, j, k, jk, jks, n, jnk, nk, nks;
    vec3<int> boxIndex3D;
    double rho;
    std::complex<double> cnm, LnmScalar;
    std::complex<double> LnmVectorA[FMMConstant::numCoefficients], LnmVectorB[FMMConstant::numCoefficients];
    
    vec4<int> r = rangevalues[threadNum];
    int numLevel = r.w;
    
    int je_forward, je_backward;
    
    for (ii = r.x; ii < r.y; ii++) {
        ib = ii + variables->levelOffset[numLevel + maxLevel];
        
//        nfip = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].parent;
        nfip = variables->domainBoxIndexMask[variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].parent];
        nfic = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].mortonIndex % 27;
        MultipoleMethodFunction::unmorton3(nfic, boxIndex3D);
        
        // rho is distance from center of child to center of parent cell
        // rho is normalized, such that child cell has length equal to one
        rho = sqrt(abs(boxIndex3D.x - 1) + abs(boxIndex3D.y - 1) + abs(boxIndex3D.z - 1));

        boxIndex3D.x = boxIndex3D.x + 4;
        boxIndex3D.y = boxIndex3D.y + 4;
        boxIndex3D.z = boxIndex3D.z + 4;
        
        rotationMorton(boxIndex3D, je_forward, je_backward);
        
        ib = variables->neo[nfip];
        
        for (i = 0; i < numCoefficients; i++) {
            LnmVectorA[i] = variables->LnmOld[ib][i];
        }

        rotation(LnmVectorA, LnmVectorB, constants->Dnm[je_forward]);

        //rho = sqrt(3.0) / 2;
        
        // *** NOTE - in the calculation for cnm, use pow(3.0, -n-1) vice pow(2.0, -n-1) since considering 3x3x3 cells
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
//            if (i < 10) {
//                printf("variables->Lnm[%d][%d] = %.9g + i %.9g\n", ii, i, std::real(variables->Lnm[ii][i]), std::imag(variables->Lnm[ii][i]));
//            }
        }
    }
}

void* MultipoleMethodForceCalculation::sm2p_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sm2p_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::sm2p_kernel(int threadNum){
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
    std::complex<double> MnmVector[FMMConstant::numCoefficients], MnmVectorB[FMMConstant::numCoefficients];
    std::complex<double> rr, rtheta, rphi, I(0.0, 1.0), eim;
    int domainMorton3;
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);

    boxSize = rootBoxSize * pow(3, numLevel);
    int zaligned, je_forward, je_backward;
    
    for (ii = range.x; ii < range.y; ii++) {
        ib = ii + variables->levelOffset[numLevel + maxLevel];
//        printf("first particle : %d, last particle : %d\n",variables->domainLinkage[numLevel][ii].firstParticle,variables->domainLinkage[numLevel][ii].lastParticle);
        for (ij = 0; ij < variables->numInteraction[ii]; ij++) {
            jj = variables->interactionList[ii][ij][0];
            jb = jj + variables->levelOffset[numLevel + maxLevel];

            for (j = 0; j < numCoefficients; j++) MnmVector[j] = variables->Mnm[jb][j];

            domainMorton3 = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[jb]].mortonIndex;

            MultipoleMethodFunction::unmorton3(domainMorton3, boxIndex3D);
            boxIndex3D.x -= variables->offsetSC[numLevel].x;
            boxIndex3D.y -= variables->offsetSC[numLevel].y;
            boxIndex3D.z -= variables->offsetSC[numLevel].z;

            boxCenter.x = variables->boxMinSC[numLevel].x + (boxIndex3D.x + 0.5) * boxSize;
            boxCenter.y = variables->boxMinSC[numLevel].y + (boxIndex3D.y + 0.5) * boxSize;
            boxCenter.z = variables->boxMinSC[numLevel].z + (boxIndex3D.z + 0.5) * boxSize;
            for (i = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].firstParticle; i <= variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].lastParticle; i++) {
                if (particleSystem->chargedParticles.acceleration[i].w == inboundCheck){
                    zaligned = 0;

                    dist.x = particleSystem->chargedParticles.position[i].x - boxCenter.x;
                    dist.y = particleSystem->chargedParticles.position[i].y - boxCenter.y;
                    dist.z = particleSystem->chargedParticles.position[i].z - boxCenter.z;
                    cart2sph_exact(r, theta, phi, dist.x, dist.y, dist.z);
                    r = sqrt(r * r + PhysicalConstant::SOFTENING);
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

void* MultipoleMethodForceCalculation::sp2p_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sp2p_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::sp2p_kernel(int threadNum){
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
            
            if (vectorizedToggle){
                DirectForceFunction::direct_vectorized(particleBounds,
                        particleSystem->chargedParticles, jptc);
            }
            else {
                DirectForceFunction::direct(particleBounds,
                        particleSystem->chargedParticles);
            }
        }
    }
}

void* MultipoleMethodForceCalculation::sp2m_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sp2m_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::sp2m_kernel(int threadNum){
 int jj, j, jb, n, m, nm, nms;
    vec3<int> boxIndex3D;
    vec3<int> domainBoxIndex3D;
    vec3<float> boxCenter, domainBoxMin;
    vec3<double> dist;
    double boxSize, rho, alpha, beta;
    double xx, s2, fact, pn, p, p1, p2;
    double rhom, rhon;
    double YnmReal[FMMConstant::numExpansion2];
    std::complex<double> MnmVector[FMMConstant::numCoefficients], I(0.0, 1.0), eim;
        
    boxSize = rootBoxSize / (1 << maxLevel);
    
    vec4<int> r = rangevalues[threadNum];
    
    int domainNum, domainMorton3;
    
    for (jj = r.x; jj < r.y; jj++){
        jb = jj + variables->levelOffset[maxLevel];
        // Find the global location of the domain
        //domainMorton3 = variables->domainLinkage[0][jj].mortonIndex;
        domainMorton3 = variables->domainLinkage[0][variables->domainBoxIndexFull[jb]].mortonIndex;
        MultipoleMethodFunction::unmorton3(domainMorton3, domainBoxIndex3D);
        domainBoxIndex3D.x -= variables->offsetSC[0].x;
        domainBoxIndex3D.y -= variables->offsetSC[0].y;
        domainBoxIndex3D.z -= variables->offsetSC[0].z;
        
        // Find the coordinates of the box center
        boxCenter.x = boxMin.x + (domainBoxIndex3D.x + 0.5) * boxSize;
        boxCenter.y = boxMin.y + (domainBoxIndex3D.y + 0.5) * boxSize;
        boxCenter.z = boxMin.z + (domainBoxIndex3D.z + 0.5) * boxSize;
        
        // Initialize MnmVector
        for (j = 0; j < numCoefficients; j++) {
            MnmVector[j] = 0;
        }
        
        // Compute MnmVector from all particles inside the jj-th cell
        for (j = variables->domainLinkage[0][variables->domainBoxIndexFull[jb]].firstParticle; j <= variables->domainLinkage[0][variables->domainBoxIndexFull[jb]].lastParticle; j++) {
            dist.x = particleSystem->chargedParticles.position[j].x - boxCenter.x;
            dist.y = particleSystem->chargedParticles.position[j].y - boxCenter.y;
            dist.z = particleSystem->chargedParticles.position[j].z - boxCenter.z;
            cart2sph_exact(rho, alpha, beta, dist.x, dist.y, dist.z);
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

void* MultipoleMethodForceCalculation::sl2p_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sl2p_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::sl2p_kernel(int threadNum){
    vec4<int> range = rangevalues[threadNum];
    
    int ii, i, ib, n, nm, nms, m;
    vec3<int> boxIndex3D;
    vec3<int> domainBoxIndex3D;
    vec3<float> boxCenter, domainBoxMin;
    vec3<double> accel, dist;
    double boxSize, r, theta, phi, accelR, accelTheta, accelPhi, accel_const;
    double xx, yy, s2, fact, pn, p, p1, p2, rn;
    double YnmReal[FMMConstant::numExpansion2], YnmRealTheta[FMMConstant::numExpansion2];
    std::complex<double> LnmVector[FMMConstant::numCoefficients];
    std::complex<double> rr, rtheta, rphi, I(0.0, 1.0), eim;
    int domainNum, domainMorton3;
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    
    boxSize = rootBoxSize / (1 << maxLevel);
    for (ii = range.x; ii < range.y; ii++) {
        ib = ii + variables->levelOffset[maxLevel];
        // Find global and local 3D box indices
        //domainMorton3 = variables->domainLinkage[0][ii].mortonIndex;
        domainMorton3 = variables->domainLinkage[0][variables->domainBoxIndexFull[ib]].mortonIndex;
        MultipoleMethodFunction::unmorton3(domainMorton3, domainBoxIndex3D);
        domainBoxIndex3D.x -= variables->offsetSC[0].x;
        domainBoxIndex3D.y -= variables->offsetSC[0].y;
        domainBoxIndex3D.z -= variables->offsetSC[0].z;
                
        boxCenter.x = boxMin.x + (domainBoxIndex3D.x + 0.5) * boxSize;
        boxCenter.y = boxMin.y + (domainBoxIndex3D.y + 0.5) * boxSize;
        boxCenter.z = boxMin.z + (domainBoxIndex3D.z + 0.5) * boxSize;
        
        // Initialize Lnm
        for (i = 0; i < numCoefficients; i++) LnmVector[i] = variables->Lnm[ii][i];
        
        // Compute particle acceleration
        for (i = variables->domainLinkage[0][variables->domainBoxIndexFull[ib]].firstParticle; i <= variables->domainLinkage[0][variables->domainBoxIndexFull[ib]].lastParticle; i++) {
            if (particleSystem->chargedParticles.acceleration[i].w == inboundCheck){
                dist.x = particleSystem->chargedParticles.position[i].x - boxCenter.x;
                dist.y = particleSystem->chargedParticles.position[i].y - boxCenter.y;
                dist.z = particleSystem->chargedParticles.position[i].z - boxCenter.z;
                cart2sph_exact(r, theta, phi, dist.x, dist.y, dist.z);

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

void MultipoleMethodForceCalculation::outlier_p2p(){
    timing.timer.startTimer();
    
#ifdef _PTHREAD_H
    int numParticles = lastParticle - firstParticle;
    int k;
    for (k = 0; k < numThreads; k++){
        dF[k].d = this;
        dF[k].threadNum = k;
        rangevalues[k].x = (k * numParticles) / numThreads + firstParticle;
        rangevalues[k].y = ((k + 1) * numParticles) / numThreads + firstParticle;
        rangevalues[k].z = lastParticle - numOutliers;
        rangevalues[k].w = lastParticle;
        
        switch (particleSystem->params.symType){
            case (SymmetryType::none):
                pthread_create(&threads[k], nullptr, direct_threaded, &dF[k]);
                break;
            case (SymmetryType::rotational):
                pthread_create(&threads[k], nullptr, direct_threaded_rotations, &dF[k]);
                break;
            case (SymmetryType::translational):
                pthread_create(&threads[k], nullptr, direct_threaded_translations, &dF[k]);
                break;
        }
    }
    
    for (int k = 0; k < numThreads; k++){
        pthread_join(threads[k], nullptr);
    }
#else
    rangevalues[0].x = firstParticle;
    rangevalues[0].y = lastParticle;
    rangevalues[0].z = lastParticle - numOutliers;
    rangevalues[0].w = lastParticle;
    switch (particleSystem->params.symType){
        case (SymmetryType::none):
            direct_kernel();
            break;
        case (SymmetryType::rotational):
            direct_kernel_rotations();
            break;
        case (SymmetryType::translational):
            direct_kernel_translations();
            break;
    }
#endif
    
    timing.outlier_p2p_time += timing.timer.getTimer();
    timing.timer.resetTimer();
}

void MultipoleMethodForceCalculation::p2outlier_p(int numBoxIndex, int toggle){
    timing.timer.startTimer();
    
#ifdef _PTHREAD_H
    int k;
    for (k = 0; k < numThreads; k++){
        dF[k].d = this;
        dF[k].threadNum = k;
        rangevalues[k].x = (k * numBoxIndex) / numThreads;
        rangevalues[k].y = ((k + 1) * numBoxIndex) / numThreads;
        rangevalues[k].z = numBoxIndex;
        rangevalues[k].w = 0;
        if (toggle == 0){
            switch (particleSystem->params.symType){
                case (SymmetryType::none):
                    pthread_create(&threads[k], nullptr, p2outlier_p_threaded, &dF[k]);
                    break;
                case (SymmetryType::rotational):
                    pthread_create(&threads[k], nullptr, p2outlier_p_threaded_rotations, &dF[k]);
                    break;
                case (SymmetryType::translational):
                    pthread_create(&threads[k], nullptr, p2outlier_p_threaded_translations, &dF[k]);
                    break;
            }
        }
        else if (toggle == 1){
            if (particleSystem->params.symType == SymmetryType::translational){
                pthread_create(&threads[k], nullptr, sp2outlier_p_threaded_translations, &dF[k]);
            }
            else {
                pthread_create(&threads[k], nullptr, sp2outlier_p_threaded, &dF[k]);
            }
        }
    }
    
    for (int k = 0; k < numThreads; k++){
        pthread_join(threads[k], nullptr);
    }
#else
    rangevalues[0].x = 0;
    rangevalues[0].y = numBoxIndex;
    rangevalues[0].z = numBoxIndex;
    rangevalues[0].w = 0;
    if (toggle == 0){
        switch (particleSystem->params.symType){
            case (SymmetryType::none):
                p2outlier_p_kernel();
                break;
            case (SymmetryType::rotational):
                p2outlier_p_kernel_rotations();
                break;
            case (SymmetryType::translational):
                p2outlier_p_kernel_translations();
                break;
        }
    }
    else if (toggle == 1){
        if (particleSystem->params.symType == SymmetryType::translational){
            sp2outlier_p_kernel_translations();
        }
        else {
            sp2outlier_p_kernel();
        }
    }
#endif
    
    if (toggle == 0){
        timing.p2outlier_p_time += timing.timer.getTimer();
    }
    else if (toggle == 1){
        timing.sp2outlier_p_time += timing.timer.getTimer();
    }
    timing.timer.resetTimer();
}

void* MultipoleMethodForceCalculation::p2outlier_p_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->p2outlier_p_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::p2outlier_p_kernel(int threadNum){
    int ii, ij, jj, i, j;
    vec3<double> dist;
    double accel_const, invDist, invDistCube, s;
    vec3<double> ai;
    
    vec4<int> r = rangevalues[threadNum];
    
    vec4<int> particleBounds;
    
    for (ii = r.x; ii < r.y; ii++){
        particleBounds.x = variables->outlierParticleOffset[0][ii];
        particleBounds.y = variables->outlierParticleOffset[1][ii] + 1;
        for (ij = 0; ij < variables->numInteraction[ii]; ij++){
            jj = variables->interactionList[ii][ij][0];
            particleBounds.z = variables->particleOffset[0][jj];
            particleBounds.w = variables->particleOffset[1][jj] + 1;
            
            if (vectorizedToggle){
                DirectForceFunction::direct_vectorized(particleBounds,
                        particleSystem->chargedParticles, jptc);
            }
            else {
                DirectForceFunction::direct(particleBounds,
                        particleSystem->chargedParticles);
            }
        }
    }
}

void* MultipoleMethodForceCalculation::sp2outlier_p_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sp2outlier_p_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::sp2outlier_p_kernel(int threadNum){
    int ii, ij, jj, i, j;
    vec3<double> dist;
    double accel_const, invDist, invDistCube, s;
    vec3<double> ai;
    
    vec4<int> r = rangevalues[threadNum];

    vec4<int> particleBounds;
    
    for (ii = r.x; ii < r.y; ii++){
        particleBounds.x = variables->outlierDomainLinkage[ii].firstParticle;
        particleBounds.y = variables->outlierDomainLinkage[ii].lastParticle + 1;
        for (ij = 0; ij < variables->numInteraction[ii]; ij++){
            jj = variables->interactionList[ii][ij][0];
            particleBounds.z = variables->domainLinkage[0][variables->domainBoxIndexFull[jj]].firstParticle;
            particleBounds.w = variables->domainLinkage[0][variables->domainBoxIndexFull[jj]].lastParticle + 1;
            
            if (vectorizedToggle){
                DirectForceFunction::direct_vectorized(particleBounds,
                        particleSystem->chargedParticles, jptc);
            }
            else {
                DirectForceFunction::direct(particleBounds,
                        particleSystem->chargedParticles);
            }
        }
    }
}

void MultipoleMethodForceCalculation::m2outlier_p(int numBoxIndex, int numLevel, int toggle){
    timing.timer.startTimer();
    
#ifdef _PTHREAD_H
    int k;
    for (k = 0; k < numThreads; k++){
        dF[k].d = this;
        dF[k].threadNum = k;
        rangevalues[k].x = (k * numBoxIndex) / numThreads;
        rangevalues[k].y = ((k + 1) * numBoxIndex) / numThreads;
        rangevalues[k].z = numLevel;
        rangevalues[k].w = 0;
        if (toggle == 0){
            switch (particleSystem->params.symType){
                case (SymmetryType::none):
                    pthread_create(&threads[k], nullptr, m2outlier_p_threaded, &dF[k]);
                    break;
                case (SymmetryType::rotational):
                    pthread_create(&threads[k], nullptr, m2outlier_p_threaded_rotations, &dF[k]);
                    break;
                case (SymmetryType::translational):
                    pthread_create(&threads[k], nullptr, m2outlier_p_threaded_translations, &dF[k]);
                    break;
            }
        }
        else if (toggle == 1){
            if (particleSystem->params.symType == SymmetryType::translational){
                pthread_create(&threads[k], nullptr, sm2outlier_p_threaded_translations, &dF[k]);
            }
            else {
                pthread_create(&threads[k], nullptr, sm2outlier_p_threaded, &dF[k]);
            }
        }
    }

    for (int k = 0; k < numThreads; k++){
        pthread_join(threads[k], nullptr);
    }
#else
    rangevalues[0].x = 0;
    rangevalues[0].y = numBoxIndex;
    rangevalues[0].z = numLevel;
    rangevalues[0].w = 0;
    if (toggle == 0){
        switch (particleSystem->params.symType){
                case (SymmetryType::none):
                    m2outlier_p_kernel();
                    break;
                case (SymmetryType::rotational):
                    m2outlier_p_kernel_rotations();
                    break;
                case (SymmetryType::translational):
                    m2outlier_p_kernel_translations();
                    break;
            }
    }
    else if (toggle == 1){
        if (particleSystem->params.symType == SymmetryType::translational){
            sm2outlier_p_kernel_translations();
        }
        else {
            sm2outlier_p_kernel();
        }
    }
#endif
    
    if (toggle == 0){
        timing.m2outlier_p_time += timing.timer.getTimer();
    }
    else if (toggle == 1){
        timing.sm2outlier_p_time += timing.timer.getTimer();
    }
    timing.timer.resetTimer();
}

void* MultipoleMethodForceCalculation::m2outlier_p_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->m2outlier_p_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::m2outlier_p_kernel(int threadNum){
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
    std::complex<double> MnmVector[FMMConstant::numCoefficients];
    std::complex<double> rr, rtheta, rphi, I(0.0, 1.0), eim;
    int domainNum, domainMorton3;
    
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    
    int tempLevelOffset = (numLevel == 0) ? variables->levelOffset[maxLevel] : variables->levelOffset[numLevel - 1];
    int zaligned;
    
    boxSize = rootBoxSize / (1 << numLevel);
    for (ii = range.x; ii < range.y; ii++) {
        for (ij = 0; ij < variables->numInteraction[ii]; ij++) {
            jj = variables->interactionList[ii][ij][0];
            jb = jj + tempLevelOffset;

            for (j = 0; j < numCoefficients; j++) MnmVector[j] = variables->Mnm[jb][j];

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
                    dist.x = particleSystem->chargedParticles.position[i].x - boxCenter.x;
                    dist.y = particleSystem->chargedParticles.position[i].y - boxCenter.y;
                    dist.z = particleSystem->chargedParticles.position[i].z - boxCenter.z;
                    cart2sph_exact(r, theta, phi, dist.x, dist.y, dist.z);
                    r = sqrt(r * r + PhysicalConstant::SOFTENING);
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

void* MultipoleMethodForceCalculation::sm2outlier_p_threaded(void* dF){
    multipoleMethodForceStruct dfc = *((multipoleMethodForceStruct*)dF);
    dfc.d->sm2outlier_p_kernel(dfc.threadNum);
    return nullptr;
}

void MultipoleMethodForceCalculation::sm2outlier_p_kernel(int threadNum){
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
    std::complex<double> MnmVector[FMMConstant::numCoefficients], MnmVectorB[FMMConstant::numCoefficients];
    std::complex<double> rr, rtheta, rphi, I(0.0, 1.0), eim;
    int domainMorton3;
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);

    boxSize = rootBoxSize * pow(3, numLevel);
    int zaligned, je_forward, je_backward;
    
    for (ii = range.x; ii < range.y; ii++) {
//        printf("first particle : %d, last particle : %d\n",variables->domainLinkage[numLevel][ii].firstParticle,variables->domainLinkage[numLevel][ii].lastParticle);
        for (ij = 0; ij < variables->numInteraction[ii]; ij++) {
            jj = variables->interactionList[ii][ij][0];
            jb = jj + variables->levelOffset[numLevel + maxLevel];

            for (j = 0; j < numCoefficients; j++) MnmVector[j] = variables->Mnm[jb][j];

            domainMorton3 = variables->domainLinkage[numLevel][variables->domainBoxIndexFull[jb]].mortonIndex;

            MultipoleMethodFunction::unmorton3(domainMorton3, boxIndex3D);
            boxIndex3D.x -= variables->offsetSC[numLevel].x;
            boxIndex3D.y -= variables->offsetSC[numLevel].y;
            boxIndex3D.z -= variables->offsetSC[numLevel].z;

            boxCenter.x = variables->boxMinSC[numLevel].x + (boxIndex3D.x + 0.5) * boxSize;
            boxCenter.y = variables->boxMinSC[numLevel].y + (boxIndex3D.y + 0.5) * boxSize;
            boxCenter.z = variables->boxMinSC[numLevel].z + (boxIndex3D.z + 0.5) * boxSize;
            for (i = variables->outlierDomainLinkage[ii].firstParticle; i <= variables->outlierDomainLinkage[ii].lastParticle; i++) {
                if (particleSystem->chargedParticles.acceleration[i].w == inboundCheck){
                    zaligned = 0;

                    dist.x = particleSystem->chargedParticles.position[i].x - boxCenter.x;
                    dist.y = particleSystem->chargedParticles.position[i].y - boxCenter.y;
                    dist.z = particleSystem->chargedParticles.position[i].z - boxCenter.z;
                    cart2sph_exact(r, theta, phi, dist.x, dist.y, dist.z);
                    r = sqrt(r * r + PhysicalConstant::SOFTENING);
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
