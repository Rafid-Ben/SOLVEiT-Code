#include "force_calculation_methods.h"

MultipoleMethodConstants::MultipoleMethodConstants(){
    int i, j;
    
    // Storage for variables related to FMM and rotations
    Anm = new float[FMMConstant::numExpansion4];
    anm = new float[4 * FMMConstant::numExpansion2];
    
    factorial = new float[4 * FMMConstant::numExpansion2];
    Ynm = new std::complex<double>[4 * FMMConstant::numExpansion2];
    Dnm = new std::complex<double>**[FMMConstant::possibleRotations];
    for (i = 0; i < FMMConstant::possibleRotations; i++) {
        Dnm[i] = new std::complex<double>*[FMMConstant::maxNumExpansions];
        for (j = 0; j < FMMConstant::maxNumExpansions; j++) Dnm[i][j] = new std::complex<double>[FMMConstant::numExpansion2];
    }

    rotationLink = new int[FMMConstant::numRelativeBox];
    
    precalc();
}

MultipoleMethodConstants::~MultipoleMethodConstants(){
    int i, j;
    
    delete[] Anm;
    delete[] anm;
    
    delete[] factorial;
    delete[] Ynm;
    for (i = 0; i < FMMConstant::possibleRotations; i++) {
        for (j = 0; j < FMMConstant::maxNumExpansions; j++) delete[] Dnm[i][j];
        delete[] Dnm[i];
    }
    delete[] Dnm;

    delete[] rotationLink;
}

void MultipoleMethodConstants::precalc(){
    int n, m, nm, nabsm, j, k, nk, npn, nmn, npm, nmm, nmk, i, nmk1, nm1k, nmk2, gcd, numRotations;
    vec3<int> boxIndex3D;
    vec3<double> dist;
    auto anmk = new double[FMMConstant::numExpansion4][2];
    double *Dnmd = new double[FMMConstant::numExpansion4];
    double fnma, fnpa, pn, p, p1, p2, anmd, anmkd, rho, alpha, beta, sc, ank, ek;
    std::complex<double> expBeta[FMMConstant::numExpansion2], I(0.0, 1.0);
    
    int count;

    int jk, jkn, jnk;
    double fnmm, fnpm, fad;

    for (n = 0; n < 2 * FMMConstant::maxNumExpansions; n++) {
        for (m = -n; m <= n; m++) {
            nm = n * n + n + m;
            nabsm = abs(m);
            fnmm = 1.0;
            for (i = 1; i <= n - m; i++) fnmm *= i;
            fnpm = 1.0;
            for (i = 1; i <= n + m; i++) fnpm *= i;
            fnma = 1.0;
            for (i = 1; i <= n - nabsm; i++) fnma *= i;
            fnpa = 1.0;
            for (i = 1; i <= n + nabsm; i++) fnpa *= i;
            factorial[nm] = sqrt(fnma / fnpa);
            fad = sqrt(fnmm * fnpm);
            anm[nm] = pow(-1.0, n) / fad;
        }
    }

    for (j = 0; j < FMMConstant::maxNumExpansions; j++) {
        for (k = -j; k <= j; k++) {
            jk = j * j + j + k;
            for (n = abs(k); n < FMMConstant::maxNumExpansions; n++) {
                nk = n * n + n + k;
                jkn = jk * FMMConstant::numExpansion2 + nk;
                jnk = (j + n) * (j + n) + j + n;
                Anm[jkn] = pow(-1.0, j + k) * anm[nk] * anm[jk] / anm[jnk];
            }
        }
    }

    pn = 1;
    for (m = 0; m < 2 * FMMConstant::maxNumExpansions; m++) {
        p = pn;
        npn = m * m + 2 * m;
        nmn = m * m;
        Ynm[npn] = factorial[npn] * p;
        Ynm[nmn] = conj(Ynm[npn]);
        p1 = p;
        p = (2 * m + 1) * p;
        for (n = m + 1; n < 2 * FMMConstant::maxNumExpansions; n++) {
            npm = n * n + n + m;
            nmm = n * n + n - m;
            Ynm[npm] = factorial[npm] * p;
            Ynm[nmm] = conj(Ynm[npm]);
            p2 = p1;
            p1 = p;
            p = ((2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
        }
        pn = 0;
    }

    for (n = 0; n < FMMConstant::maxNumExpansions; n++) {
        for (m = 1; m <= n; m++) {
            anmd = n * (n + 1) - m * (m - 1);
            for (k = 1 - m; k < m; k++) {
                nmk = (4 * n * n * n + 6 * n * n + 5 * n) / 3 + m * (2 * n + 1) + k;
                anmkd = ((double)(n * (n + 1) - k * (k + 1))) / (n * (n + 1) - m * (m - 1));
                anmk[nmk][0] = -(m + k) / sqrt(anmd);
                anmk[nmk][1] = sqrt(anmkd);
            }
        }
    }
    
    numRotations = 0;
    
    for (i = 0; i < FMMConstant::numRelativeBox; i++) {
        // Use unmorton7 if only considering 2x2x2 boxes
        // Use unmorton11 if considering 3x3x3 boxes
        // Use unmorton(4n - 1) if considering nxnxn boxes (3*n + (n - 1) = 4n - 1)
        MultipoleMethodFunction::unmorton11(i, boxIndex3D);
        
        // Use for 2x2x2 boxes
//        dist.x = boxIndex3D.x - 3;
//        dist.y = boxIndex3D.y - 3;
//        dist.z = boxIndex3D.z - 3;
        // Use for 3x3x3 boxes
        dist.x = boxIndex3D.x - 5;
        dist.y = boxIndex3D.y - 5;
        dist.z = boxIndex3D.z - 5;
        // In general, dist.i = boxIndex3D.i - (2n - 1), if considering nxnxn boxes
        
        gcd = gcd3(dist.x, dist.y, dist.z);
        
        if (gcd > 1){
            // If gcd > 1, there is an equivalent relation for the rotation already stored in Dnm
            rotationLink[i] = -1;
            continue;
        }
        rotationLink[i] = numRotations;

        cart2sph_exact(rho, alpha, beta, dist.x, dist.y, dist.z);
        if (rho == 0) alpha = 0;
        sc = sin(alpha) / (1 + cos(alpha));

        for (n = 0; n < 4 * FMMConstant::maxNumExpansions - 3; n++) {
            expBeta[n] = exp((n - 2 * FMMConstant::maxNumExpansions + 2) * beta * I);
        }

        for (n = 0; n < FMMConstant::maxNumExpansions; n++) {
            nmk = (4 * n * n * n + 6 * n * n + 5 * n) / 3 + n * (2 * n + 1) + n;
            Dnmd[nmk] = pow(cos(alpha * 0.5), 2 * n);
            if ((rho == -dist.z) && (rho > 0) && (n > 0)){
                Dnmd[nmk] = 0;
            }
            for (k = n; k >= 1 - n; k--) {
                nmk = (4 * n * n * n + 6 * n * n + 5 * n) / 3 + n * (2 * n + 1) + k;
                nmk1 = (4 * n * n * n + 6 * n * n + 5 * n) / 3 + n * (2 * n + 1) + k - 1;
                ank = ((double)n + k) / (n - k + 1);
                Dnmd[nmk1] = sqrt(ank) * tan(alpha * 0.5) * Dnmd[nmk];
                if ((rho == -dist.z) && (rho > 0)){
                    Dnmd[nmk1] = 0;
                    if (k == 1 - n){
                        Dnmd[nmk1] = 1;
                    }
                }
            }
            count = 0;
            for (m = n; m >= 1; m--) {
                for (k = m - 1; k >= 1 - m; k--) {
                    nmk = (4 * n * n * n + 6 * n * n + 5 * n) / 3 + m * (2 * n + 1) + k;
                    nmk1 = (4 * n * n * n + 6 * n * n + 5 * n) / 3 + m * (2 * n + 1) + k + 1;
                    nm1k = (4 * n * n * n + 6 * n * n + 5 * n) / 3 + (m - 1) * (2 * n + 1) + k;
                    Dnmd[nm1k] = anmk[nmk][1] * Dnmd[nmk1] + anmk[nmk][0] * sc * Dnmd[nmk];
                    if ((rho == -dist.z) && (rho > 0)){
                        Dnmd[nm1k] = 0;
                        count++;
                        if (k == 1 - m){
                            Dnmd[nm1k] = pow(-1.0,count);
                        }
                    }
                }
            }
        }

        for (n = 1; n < FMMConstant::maxNumExpansions; n++) {
            for (m = 0; m <= n; m++) {
                for (k = -m; k <= -1; k++) {
                    ek = pow(-1.0, k);
                    nmk = (4 * n * n * n + 6 * n * n + 5 * n) / 3 + m * (2 * n + 1) + k;
                    nmk1 = (4 * n * n * n + 6 * n * n + 5 * n) / 3 - k * (2 * n + 1) - m;
                    Dnmd[nmk] = ek * Dnmd[nmk];
                    Dnmd[nmk1] = pow(-1.0, m + k) * Dnmd[nmk];
                }
                for (k = 0; k <= m; k++) {
                    nmk = (4 * n * n * n + 6 * n * n + 5 * n) / 3 + m * (2 * n + 1) + k;
                    nmk1 = (4 * n * n * n + 6 * n * n + 5 * n) / 3 + k * (2 * n + 1) + m;
                    nmk2 = (4 * n * n * n + 6 * n * n + 5 * n) / 3 - k * (2 * n + 1) - m;
                    Dnmd[nmk1] = pow(-1.0, m + k) * Dnmd[nmk];
                    Dnmd[nmk2] = Dnmd[nmk1];
                }
            }
        }

        for (n = 0; n < FMMConstant::maxNumExpansions; n++) {
            for (m = 0; m <= n; m++) {
                for (k = -n; k <= n; k++) {
                    nmk = (4 * n * n * n + 6 * n * n + 5 * n) / 3 + m * (2 * n + 1) + k;
                    nk = n * (n + 1) + k;
                    //Dnm[i][m][nk] = Dnmd[nmk] * expBeta[k + m + 2 * maxNumExpansions - 2];
                    Dnm[numRotations][m][nk] = Dnmd[nmk] * expBeta[k + m + 2 * FMMConstant::maxNumExpansions - 2];
                }
            }
        }
        
        numRotations++;
        
    }
    // Flipping the signs of alpha and beta simultaneously is equivalent to flipping the sign of the x-coordinate
    // This is already stored in Dnm, so there is no need to double its size, as in previous versions

    delete[] Dnmd;
    delete[] anmk;
}



void MultipoleMethodVariables::allocateMMvariables(int maxNumParticles,
        int numBoxIndexLeaf, int numBoxIndexTotal, int numBoxIndexFull,
        int numExpansions, int maxLevel, int numLevelSC, int numDomains,
        int numOutlierDomains){
    
    if (maxNumParticles > maxMaxNumParticles){
        if (sortValue != nullptr) delete[] sortValue;
        if (sortIndex != nullptr) delete[] sortIndex;
        if (sortValueBuffer != nullptr) delete[] sortValueBuffer;
        if (permutation != nullptr) delete[] permutation;
        if (outlierCheck != nullptr) delete[] outlierCheck;
        if (domainMortonIndex != nullptr) delete[] domainMortonIndex;
        
        sortValue = new int[maxNumParticles];
        sortIndex = new int[maxNumParticles];
        sortValueBuffer = new int[maxNumParticles];
        permutation = new int[maxNumParticles];
        outlierCheck = new int[maxNumParticles];
        domainMortonIndex = new int[maxNumParticles];
    }
    
    
    int i, j;

    // particleOffset stores indices of first and last particle in some subset of the domain
    // particleOffset[0][i] gives the index of the first particle in domain i
    // particleOffset[1][i] gives the index of the last particle in domain i
    if (numBoxIndexLeaf > maxNumBoxIndexLeaf){
        if (particleOffset != nullptr){
            for (i = 0; i < 2; i++) delete[] particleOffset[i];
            delete[] particleOffset;
        }
        particleOffset = new int* [2];
        for (i = 0; i < 2; i++) particleOffset[i] = new int[numBoxIndexLeaf];
    }

    // levelOffset stores offset values for the number of cells in each level of the domain
    if (maxLevel > maxMaxLevel || numLevelSC > maxNumLevelSC){
        if (levelOffset != nullptr) delete[] levelOffset;
        if (std::max(maxLevel,maxMaxLevel) + std::max(numLevelSC,maxNumLevelSC) > 0){
            levelOffset = new int[std::max(maxLevel,maxMaxLevel) + std::max(numLevelSC,maxNumLevelSC)];
        }
    }
    
    if (numLevelSC > maxNumLevelSC){
        if (dimValsSC != nullptr) delete[] dimValsSC;
        if (boxMinSC != nullptr) delete[] boxMinSC;
        if (offsetSC != nullptr) delete[] offsetSC;
        if (numDomainsLevel != nullptr) delete[] numDomainsLevel;
        if (boxMinIndexSC != nullptr) delete[] boxMinIndexSC;
        
        dimValsSC = new vec3<int>[numLevelSC];
        boxMinSC = new vec3<float>[numLevelSC];
        offsetSC = new vec3<int>[numLevelSC];
        numDomainsLevel = new int[numLevelSC];
        boxMinIndexSC = new vec3<int>[numLevelSC];
        
        if (domainLinkage != nullptr){
            for (i = 0; i < maxNumLevelSC; i++) delete[] domainLinkage[i];
            delete[] domainLinkage;
        }
        domainLinkage = new domainInfo*[numLevelSC];
    }

    // numInteraction stores the number of interactions each cell has
    if (numBoxIndexLeaf > maxNumBoxIndexLeaf){
        if (numInteraction != nullptr) delete[] numInteraction;
        numInteraction = new int[numBoxIndexLeaf];
    }


    // interactionList[i][j] stores the index of the j-th cell interacting with cell i;
    if (numBoxIndexLeaf > maxNumBoxIndexLeaf){
        if (interactionList != nullptr){
            for (i = 0; i < maxNumBoxIndexLeaf; i++){
                for (j = 0; j < FMMConstant::maxM2LInteraction; j++){
                    delete[] interactionList[i][j];
                }
                delete[] interactionList[i];
            }
            delete[] interactionList;
        }
        
        interactionList = new int**[numBoxIndexLeaf];
        for (i = 0; i < numBoxIndexLeaf; i++) {
            interactionList[i] = new int*[FMMConstant::maxM2LInteraction];
            for (j = 0; j < FMMConstant::maxM2LInteraction; j++){
                interactionList[i][j] = new int[2];
            }
        }
    }
    
    if (numBoxIndexFull > maxNumBoxIndexFull){
        if (sortIndexBuffer != nullptr) delete[] sortIndexBuffer;
        sortIndexBuffer = new int[numBoxIndexFull];
    }
    
    int tempNumExpansions = std::max(numExpansions, maxNumExpansions);
    int tempNumCoefficents = tempNumExpansions * (tempNumExpansions + 1) / 2;

    if (numBoxIndexLeaf > maxNumBoxIndexLeaf || numExpansions > maxNumExpansions){
        if (Lnm != nullptr){
            if (maxNumBoxIndexLeaf > 0){
                for (i = 0; i < maxNumBoxIndexLeaf; i++){
                    delete[] Lnm[i];
                }
                delete[] Lnm;
            }
        }
        if (LnmOld != nullptr){
            if (maxNumBoxIndexLeaf > 0){
                for (i = 0; i < maxNumBoxIndexLeaf; i++){
                    delete[] LnmOld[i];
                }
                delete[] LnmOld;
            }
        }
        
        if (std::max(numBoxIndexLeaf, maxNumBoxIndexLeaf) > 0){
            Lnm = new std::complex<double>*[std::max(numBoxIndexLeaf,maxNumBoxIndexLeaf)];
            LnmOld = new std::complex<double>*[std::max(numBoxIndexLeaf,maxNumBoxIndexLeaf)];
            for (i = 0; i < std::max(numBoxIndexLeaf,maxNumBoxIndexLeaf); i++){
                if (tempNumCoefficents > 0){
                    Lnm[i] = new std::complex<double>[tempNumCoefficents];
                    LnmOld[i] = new std::complex<double>[tempNumCoefficents];
                }
            }
        }
    }
    
    if (numBoxIndexTotal > maxNumBoxIndexTotal || numExpansions > maxNumExpansions){
        if (Mnm != nullptr){
            if (maxNumBoxIndexTotal > 0){
                for (i = 0; i < maxNumBoxIndexTotal; i++){
                    delete[] Mnm[i];
                }
                delete[] Mnm;
            }
        }
        if (std::max(numBoxIndexTotal, maxNumBoxIndexTotal) > 0){
            Mnm = new std::complex<double>*[std::max(numBoxIndexTotal, maxNumBoxIndexTotal)];
            for (i = 0; i < std::max(numBoxIndexTotal, maxNumBoxIndexTotal); i++){
                if (tempNumCoefficents > 0){
                    Mnm[i] = new std::complex<double>[tempNumCoefficents];
                    for (j = 0; j < tempNumCoefficents; j++){
                        Mnm[i][j] = 0;
                    }
                }
            }
        }
    }
//    Lnm = new std::complex<double>[numBoxIndexLeaf][numCoefficients];
//    LnmOld = new std::complex<double>[numBoxIndexLeaf][numCoefficients];
//    Mnm = new std::complex<double>[numBoxIndexTotal][numCoefficients];
////    Lnm = new std::complex<double>*[numBoxIndexLeaf];
////    LnmOld = new std::complex<double>*[numBoxIndexLeaf];
////    for (i = 0; i < numBoxIndexLeaf; i++){
////        Lnm[i] = new std::complex<double>[numCoefficients];
////        LnmOld[i] = new std::complex<double>[numCoefficients];
////    }
////    Mnm = new std::complex<double>*[numBoxIndexTotal];
////    for (i = 0; i < numBoxIndexTotal; i++){
////        Mnm[i] = new std::complex<double>[numCoefficients];
////    }
//
//    for (j = 0; j < numBoxIndexTotal; j++) {
//        for (i = 0; i < numCoefficients; i++) {
//            Mnm[j][i] = 0;
//        }
//    }

    // Links local morton index to global box index
    if (numBoxIndexFull > maxNumBoxIndexFull){
        if (domainBoxIndexMask != nullptr) delete[] domainBoxIndexMask;
        domainBoxIndexMask = new int[numBoxIndexFull];
    }

    // Links global box index to local morton index
    if (numBoxIndexTotal > maxNumBoxIndexTotal){
        if (domainBoxIndexFull != nullptr) delete[] domainBoxIndexFull;
        domainBoxIndexFull = new int[numBoxIndexTotal];
    }

    // domainOffset stores indices of the first and last next lowest level cell in some subset of the domain
    // domainOffset[0][i] gives the index of the first cell in domain i
    // domainOffset[1][i] gives the index of the last cell in domain i
    if (numDomains > maxNumDomains || maxLevel > maxMaxLevel){
        if (domainOffset != nullptr){
            for (i = 0; i < 2; i++){
                if (maxNumDomains > 0){
                    for (j = 0; j < maxNumDomains; j++){
                        if (maxMaxLevel + 1 > 0){
                            delete[] domainOffset[i][j];
                        }
                    }
                    delete[] domainOffset[i];
                }
            }
            delete[] domainOffset;
        }
        
        domainOffset = new int** [2];
        for (i = 0; i < 2; i++) {
            if (std::max(numDomains, maxNumDomains) > 0){
                domainOffset[i] = new int* [std::max(numDomains, maxNumDomains)];
                for (j = 0; j < std::max(numDomains, maxNumDomains); j++) {
                    if (std::max(maxLevel,maxMaxLevel) + 1 > 0){
                        domainOffset[i][j] = new int[std::max(maxLevel,maxMaxLevel) + 1];
                    }
                }
            }
        }
    }

    if (numOutlierDomains > maxNumOutlierDomains || maxLevel > maxMaxLevel){
        if (outlierDomainOffset != nullptr){
            for (i = 0; i < 2; i++){
                if (maxNumOutlierDomains > 0){
                    for (j = 0; j < maxNumOutlierDomains; j++){
                        if (maxMaxLevel + 1 > 0){
                            delete[] outlierDomainOffset[i][j];
                        }
                    }
                    delete[] outlierDomainOffset[i];
                }
            }
            delete[] outlierDomainOffset;
        }
        
        outlierDomainOffset = new int** [2];
        for (i = 0; i < 2; i++) {
            if (std::max(numOutlierDomains, maxNumOutlierDomains) > 0){
                outlierDomainOffset[i] = new int* [std::max(numOutlierDomains, maxNumOutlierDomains)];
                for (j = 0; j < std::max(numOutlierDomains, maxNumOutlierDomains); j++) {
                    if (std::max(maxLevel,maxMaxLevel) + 1 > 0){
                        outlierDomainOffset[i][j] = new int[std::max(maxLevel,maxMaxLevel) + 1];
                    }
                }
            }
        }
    }

    if (numBoxIndexLeaf > maxNumBoxIndexLeaf){
        if (outlierParticleOffset != nullptr){
            for (i = 0; i < 2; i++) delete[] outlierParticleOffset[i];
            delete[] outlierParticleOffset;
        }
        
        outlierParticleOffset = new int* [2];
        for (i = 0; i < 2; i++) outlierParticleOffset[i] = new int[numBoxIndexLeaf];
    }

    if (numBoxIndexFull > maxNumBoxIndexFull){
        if (neo != nullptr) delete[] neo;
        neo = new int[numBoxIndexFull];
    }
        
    maxMaxNumParticles = std::max(maxMaxNumParticles, maxNumParticles);
    maxNumBoxIndexLeaf = std::max(maxNumBoxIndexLeaf, numBoxIndexLeaf);
    maxNumBoxIndexTotal = std::max(maxNumBoxIndexTotal, numBoxIndexTotal);
    maxNumBoxIndexFull = std::max(maxNumBoxIndexFull, numBoxIndexFull);
    maxNumExpansions = std::max(maxNumExpansions, numExpansions);
    maxMaxLevel = std::max(maxMaxLevel, maxLevel);
    maxNumLevelSC = std::max(maxNumLevelSC, numLevelSC);
    maxNumDomains = std::max(maxNumDomains, numDomains);
    maxNumOutlierDomains = std::max(maxNumOutlierDomains, numOutlierDomains);
}

void MultipoleMethodVariables::allocateDomainVariables(int numLevelSC, 
        int numDomains, int avgPPB){
    
    int i;
    if (numLevelSC > maxNumLevelSC){
        if (dimValsSC != nullptr) delete[] dimValsSC;
        if (boxMinSC != nullptr) delete[] boxMinSC;
        if (offsetSC != nullptr) delete[] offsetSC;
        if (numDomainsLevel != nullptr) delete[] numDomainsLevel;
        if (boxMinIndexSC != nullptr) delete[] boxMinIndexSC;
        
        dimValsSC = new vec3<int>[numLevelSC];
        boxMinSC = new vec3<float>[numLevelSC];
        offsetSC = new vec3<int>[numLevelSC];
        numDomainsLevel = new int[numLevelSC];
        boxMinIndexSC = new vec3<int>[numLevelSC];
        
        if (domainLinkage != nullptr){
            for (i = 0; i < maxNumLevelSC; i++) delete[] domainLinkage[i];
            delete[] domainLinkage;
        }
        domainLinkage = new domainInfo*[numLevelSC];
        
        if (levelOffset != nullptr) delete[] levelOffset;
        levelOffset = new int[maxMaxLevel + numLevelSC];
    }
    
    if (numDomains > maxNumDomains){
        if (tempParticles != nullptr){
            for (i = 0; i < maxNumDomains + 1; i++) delete[] tempParticles[i];
            delete[] tempParticles;
        }
        if (tempParticlesPerDomain != nullptr) delete[] tempParticlesPerDomain;
        
        tempParticles = new int*[numDomains + 1];
        tempParticlesPerDomain = new int*[numDomains + 1];
        for (i = 0; i < numDomains + 1; i++){
            tempParticles[i] = new int[avgPPB];
            tempParticlesPerDomain[i] = new int[2];
            tempParticlesPerDomain[i][0] = avgPPB;
            tempParticlesPerDomain[i][1] = 0;
        }
    }
    
    for (i = 0; i < numDomains + 1; i++){
        tempParticlesPerDomain[i][1] = 0;
    }
    
    int j;
    if (numDomains > maxNumDomains){
        if (domainOffset != nullptr){
            for (i = 0; i < 2; i++){
                if (maxNumDomains > 0){
                    for (j = 0; j < maxNumDomains; j++){
                        if (maxMaxLevel + 1 > 0){
                            delete[] domainOffset[i][j];
                        }
                    }
                    delete[] domainOffset[i];
                }
            }
            delete[] domainOffset;
        }
        
        domainOffset = new int** [2];
        for (i = 0; i < 2; i++) {
            if (numDomains > 0){
                domainOffset[i] = new int* [numDomains];
                for (j = 0; j < numDomains; j++) {
                    if (maxMaxLevel + 1 > 0){
                        domainOffset[i][j] = new int[maxMaxLevel + 1];
                    }
                }
            }
        }
    }
    
    maxNumLevelSC = std::max(maxNumLevelSC, numLevelSC);
    maxNumDomains = std::max(maxNumDomains, numDomains);
}

void MultipoleMethodVariables::reallocateTempParticles(int tempIndex, int additionalParticles){
    tempParticlesPerDomain[tempIndex][0] += std::min(additionalParticles, tempParticlesPerDomain[tempIndex][0]);
    tempParticles[tempIndex] = (int*) realloc(tempParticles[tempIndex], tempParticlesPerDomain[tempIndex][0] * sizeof (int));
}

void MultipoleMethodVariables::allocateDomainLinkage(){
    int i;
    
    for (i = 0; i < maxNumLevelSC; i++){
        domainLinkage[i] = new domainInfo[numDomainsLevel[i]];
    }
}

void MultipoleMethodVariables::allocateOutlierVariables(int numOutlierDomains, 
        int numOutliers, int outlierNumBoxIndexLeaf){
    
    if (numOutlierDomains > maxNumOutlierDomains){
        if (outlierDomainLinkage != nullptr) delete[] outlierDomainLinkage;
        outlierDomainLinkage = new outlierDomainInfo[numOutlierDomains];
    }
    if (numOutliers > maxNumOutliers){
        if (outlierList != nullptr) delete[] outlierList;
        outlierList = new outlierInfo[numOutliers];
    }
    if (outlierNumBoxIndexLeaf > maxOutlierNumBoxIndexLeaf){
        if (outlierBoxIndexFull != nullptr) delete[] outlierBoxIndexFull;
        outlierBoxIndexFull = new int[outlierNumBoxIndexLeaf];
    }
}

void MultipoleMethodVariables::deallocateAll(){
    if (sortValue != nullptr) delete[] sortValue;
    if (sortIndex != nullptr) delete[] sortIndex;
    if (sortValueBuffer != nullptr) delete[] sortValueBuffer;
    if (permutation != nullptr) delete[] permutation;
    if (outlierCheck != nullptr) delete[] outlierCheck;
    if (domainMortonIndex != nullptr) delete[] domainMortonIndex;

    int i, j;
    if (particleOffset != nullptr){
        for (i = 0; i < 2; i++) delete[] particleOffset[i];
        delete[] particleOffset;
    }

    if (levelOffset != nullptr) delete[] levelOffset;

    if (numInteraction != nullptr) delete[] numInteraction;

    if (interactionList != nullptr){
        for (i = 0; i < maxNumBoxIndexLeaf; i++){
            for (j = 0; j < FMMConstant::maxM2LInteraction; j++){
                delete[] interactionList[i][j];
            }
            delete[] interactionList[i];
        }
        delete[] interactionList;
    }

    if (sortIndexBuffer != nullptr) delete[] sortIndexBuffer;

    if (Lnm != nullptr){
        if (maxNumBoxIndexLeaf > 0){
            for (i = 0; i < maxNumBoxIndexLeaf; i++){
                delete[] Lnm[i];
            }
            delete[] Lnm;
        }
    }
    if (LnmOld != nullptr){
        if (maxNumBoxIndexLeaf > 0){
            for (i = 0; i < maxNumBoxIndexLeaf; i++){
                delete[] LnmOld[i];
            }
            delete[] LnmOld;
        }
    }
    if (Mnm != nullptr){
        if (maxNumBoxIndexTotal > 0){
            for (i = 0; i < maxNumBoxIndexTotal; i++){
                delete[] Mnm[i];
            }
            delete[] Mnm;
        }
    }

    if (domainBoxIndexMask != nullptr) delete[] domainBoxIndexMask;

    if (domainBoxIndexFull != nullptr) delete[] domainBoxIndexFull;

    if (domainOffset != nullptr){
        for (i = 0; i < 2; i++){
            if (maxNumDomains > 0){
                for (j = 0; j < maxNumDomains; j++){
                    if (maxMaxLevel + 1 > 0){
                        delete[] domainOffset[i][j];
                    }
                }
                delete[] domainOffset[i];
            }
        }
        delete[] domainOffset;
    }

    if (outlierDomainOffset != nullptr){
        for (i = 0; i < 2; i++){
            if (maxNumOutlierDomains > 0){
                for (j = 0; j < maxNumOutlierDomains; j++){
                    if (maxMaxLevel + 1 > 0){
                        delete[] outlierDomainOffset[i][j];
                    }
                }
                delete[] outlierDomainOffset[i];
            }
        }
        delete[] outlierDomainOffset;
    }

    if (outlierParticleOffset != nullptr){
        for (i = 0; i < 2; i++) delete[] outlierParticleOffset[i];
        delete[] outlierParticleOffset;
    }

    if (neo != nullptr) delete[] neo;
    
    if (dimValsSC != nullptr) delete[] dimValsSC;
    if (boxMinSC != nullptr) delete[] boxMinSC;
    if (offsetSC != nullptr) delete[] offsetSC;
    if (numDomainsLevel != nullptr) delete[] numDomainsLevel;
    if (boxMinIndexSC != nullptr) delete[] boxMinIndexSC;

    if (domainLinkage != nullptr){
        for (i = 0; i < maxNumLevelSC; i++) delete[] domainLinkage[i];
        delete[] domainLinkage;
    }

    if (tempParticles != nullptr){
        for (i = 0; i < maxNumDomains + 1; i++) delete[] tempParticles[i];
        delete[] tempParticles;
    }
    if (tempParticlesPerDomain != nullptr) delete[] tempParticlesPerDomain;
    
    if (outlierDomainLinkage != nullptr) delete[] outlierDomainLinkage;

    if (outlierList != nullptr) delete[] outlierList;

    if (outlierBoxIndexFull != nullptr) delete[] outlierBoxIndexFull;
}

MultipoleMethodVariables::~MultipoleMethodVariables(){
    deallocateAll();
}
