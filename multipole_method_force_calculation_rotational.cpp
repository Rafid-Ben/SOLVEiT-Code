/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "force_calculation_methods.h"

void MultipoleMethodForceCalculation::updateAccelerationRotational(){
    int numParticles = lastParticle - firstParticle;
    int i, numLevel, numBoxIndex, numBoxIndexOld;
    
    if (timeInfoFile != nullptr){
        fwrite(&numParticles, sizeof(*&numParticles), 1, timeInfoFile);
    }
    
    timing.resetTiming();
    
    timing.total_timer.startTimer();
    timing.timer.startTimer();
    
    numBoxIndexLeaf = 0;
    numBoxIndexTotal = 0;
    numBoxIndexFull = 0;
    maxLevel = 0;
    numLevelSC = 0;
    numDomains = 0;
    numOutlierDomains = 0;
    numOutliers = 0;
    outlierNumBoxIndexLeaf = 0;
    
    variables->allocateMMvariables(particleSystem->chargedParticles.maxNumParticles, 
            numBoxIndexLeaf, numBoxIndexTotal, numBoxIndexFull, numExpansions,
            maxLevel, numLevelSC, numDomains, numOutlierDomains);
    
    if (outlierToggle == 1){
        DomainBoundsFunction::isOutlier(particleSystem->chargedParticles,
                variables->outlierCheck, numOutliers);
    }
    else if (outlierToggle == 0){
        numOutliers = 0;
        for (i = 0; i < numParticles; i++){
            variables->outlierCheck[i] = 0;
        }
    }

    variables->allocateOutlierVariables(numOutlierDomains, numOutliers, outlierNumBoxIndexLeaf);
    
    setDomainSizeRotations();
    
    std::sort(variables->outlierList, variables->outlierList + numOutliers, compareOutlierInfo);
    
    numOutlierDomains = 0;
    for (i = 0; i < numOutliers; i++) {
        if (i == 0 || variables->outlierList[i].signedMorton != variables->outlierList[i - 1].signedMorton) {
            numOutlierDomains++;
        }
    }
    variables->allocateOutlierVariables(numOutlierDomains, numOutliers, outlierNumBoxIndexLeaf);

    setOptimumLevelRotations();
    
    if (timeInfoFile != nullptr){
        fwrite(&numOutliers, sizeof(*&numOutliers), 1, timeInfoFile);
        fwrite(&boxDimVals[0], sizeof(*&boxDimVals[0]), 1, timeInfoFile);
        fwrite(&boxDimVals[1], sizeof(*&boxDimVals[1]), 1, timeInfoFile);
        fwrite(&boxDimVals[2], sizeof(*&boxDimVals[2]), 1, timeInfoFile);
        fwrite(&maxLevel, sizeof(*&maxLevel), 1, timeInfoFile);
    }

    if (boxInfo != nullptr){
        fwrite(&boxMin.x, sizeof (*&boxMin.x), 1, boxInfo);
        fwrite(&boxMin.y, sizeof (*&boxMin.y), 1, boxInfo);
        fwrite(&boxMin.z, sizeof (*&boxMin.z), 1, boxInfo);
        fwrite(&rootBoxSize, sizeof (*&rootBoxSize), 1, boxInfo);
        fwrite(&boxDimVals[0], sizeof (*&boxDimVals[0]), 1, boxInfo);
        fwrite(&boxDimVals[1], sizeof (*&boxDimVals[1]), 1, boxInfo);
        fwrite(&boxDimVals[2], sizeof (*&boxDimVals[2]), 1, boxInfo);
        fwrite(&maxLevel, sizeof (*&maxLevel), 1, boxInfo);
    }
    
    if (maxLevel == -1) {
        
        timing.setup_time = timing.timer.getTimer();
        timing.timer.resetTimer();
        
        int k;
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
        
        // Small number of particles -- use direct method
        direct();
    } else {
        variables->allocateMMvariables(particleSystem->chargedParticles.maxNumParticles, 
                numBoxIndexLeaf, numBoxIndexTotal, numBoxIndexFull, numExpansions,
                maxLevel, numLevelSC, numDomains, numOutlierDomains);

        sortParticles();
        
        if (vectorizedToggle){
            int tempNumParticles;
            DirectForceFunction::setJpdata(tempNumParticles, particleSystem->chargedParticles, jptc);
        }
        
        countNonEmptyBoxes();

        variables->allocateMMvariables(particleSystem->chargedParticles.maxNumParticles, 
                numBoxIndexLeaf, numBoxIndexTotal, numBoxIndexFull, numExpansions,
                maxLevel, numLevelSC, numDomains, numOutlierDomains);

        numLevel = maxLevel;

        for (i = firstParticle; i < lastParticle; i++) {
            particleSystem->chargedParticles.acceleration[i].x = 0;
            particleSystem->chargedParticles.acceleration[i].y = 0;
            particleSystem->chargedParticles.acceleration[i].z = 0;
        }
        
        timing.setup_time = timing.timer.getTimer();
        timing.timer.resetTimer();
        
        if (numLevel > 0) {
            variables->levelOffset[numLevel - 1] = 0;

            getBoxData(numBoxIndex, 0);

            // P2P
            getInteractionListRotations(numBoxIndex, numLevel, FmmInteractionType::P2P);

            if (numOutliers > 0) {
                outlier_p2p();
            }
            
            p2p(numBoxIndex, 0);
            
            // Interaction list for interactions between outliers and particles close by
            if (numOutliers > 0) {
                getInteractionListRotations(outlierNumBoxIndex, numLevel, FmmInteractionType::P2P_Outliers);
                p2outlier_p(outlierNumBoxIndex, 0);
            }

            // P2M
            p2m(numBoxIndex, 0);
        } else {
            variables->levelOffset[0] = 0;
            getBoxData(numBoxIndex, 1);
            
//            printf("variables->domainBoxIndexFull[0] = %d\n",variables->domainBoxIndexFull[0]);
            
            // P2P for supercells
            getInteractionListRotations(numBoxIndex, numLevel, FmmInteractionType::P2P_Supercell);
            if (numOutliers > 0) {
                outlier_p2p();
            }

            p2p(numBoxIndex, 1);
            
//            printf("variables->domainBoxIndexFull[0] = %d\n",variables->domainBoxIndexFull[0]);
            

            // Interaction list for interactions between outliers and particles close by
            if (numOutliers > 0) {
                getInteractionListRotations(outlierNumBoxIndex, numLevel, FmmInteractionType::P2P_SupercellOutliers);
                p2outlier_p(outlierNumBoxIndex, 1);
            }

            // P2M for supercells
            p2m(numBoxIndex, 1);
        }
        
        numLevel = maxLevel;

        int tempLevel = maxLevel;
        int tempToggle = (maxLevel == 0) ? 1 : 0;

        //if (maxLevel + numLevelSC > 3){
        if (maxLevel + numLevelSC <= 0) {
            getBoxIndexMask(numBoxIndex, numLevel, 0);
        }
        
        for (numLevel = maxLevel - 1; numLevel >= std::max(0, 2 - numLevelSC); numLevel--){
            if (treeOrFMM == 0) {

                // M2P at lower levels

                getInteractionListRotations(numBoxIndex, numLevel + 1, FmmInteractionType::M2L_Subcell);

                m2p(numBoxIndex, numLevel + 1, 0);

            }
            if (numOutliers > 0) {
                getInteractionListRotations(outlierNumBoxIndex, numLevel + 1, FmmInteractionType::M2L_SubcellOutliers);
                m2outlier_p(outlierNumBoxIndex, numLevel + 1, 0);
            }

            numBoxIndexOld = numBoxIndex;

            if (numLevel == 0) {
                getBoxDataOfParent(numBoxIndex, numLevel, 1);
            } else {
                getBoxDataOfParent(numBoxIndex, numLevel, 0);
            }

            m2m(numBoxIndex, numBoxIndexOld, numLevel, 0);

            tempLevel = numLevel;
            tempToggle = 0;
        }
        
        for (numLevel = 1; numLevel <= numLevelSC - 1; numLevel++){
            if (treeOrFMM == 0) {
                // M2P at supercell levels

                getInteractionListRotations(numBoxIndex, numLevel - 1, FmmInteractionType::M2L_Supercell);

                m2p(numBoxIndex, numLevel - 1, 1);
            }
            if (numOutliers > 0) {
                getInteractionListRotations(outlierNumBoxIndex, numLevel - 1, FmmInteractionType::M2L_SupercellOutliers);
                m2outlier_p(outlierNumBoxIndex, numLevel - 1, 1);
            }

            numBoxIndexOld = numBoxIndex;

            getBoxDataOfParent(numBoxIndex, numLevel, 1);

            m2m(numBoxIndex, numBoxIndexOld, numLevel, 1);

            tempLevel = numLevel;
            tempToggle = 1;
        }
        
        numLevel = tempLevel;
        FmmInteractionType tempInteractionType;

        if (numOutliers > 0) {
            if (numLevel == 0) tempToggle = 1;
            if (tempToggle == 0){
                tempInteractionType = FmmInteractionType::M2L_SubcellTopLevelOutliers;
            }
            else{
                tempInteractionType = FmmInteractionType::M2L_SupercellTopLevelOutliers;
            }
            getInteractionListRotations(numBoxIndex, tempLevel, tempInteractionType); // Using numBoxIndex instead of outlierNumBoxIndex is correct
            m2outlier_p(outlierNumBoxIndex, tempLevel, tempToggle);
        }

        if (treeOrFMM == 0) {
            if (numLevel == 0) tempToggle = 1;
            if (tempToggle == 0){
                tempInteractionType = FmmInteractionType::M2L_SubcellTopLevel;
            }
            else{
                tempInteractionType = FmmInteractionType::M2L_SupercellTopLevel;
            }
            getInteractionListRotations(numBoxIndex, tempLevel, tempInteractionType);

            m2p(numBoxIndex, tempLevel, tempToggle);
            
        } else if (treeOrFMM == 1) {
            if (numLevel == 0) tempToggle = 1;
            if (tempToggle == 0){
                tempInteractionType = FmmInteractionType::M2L_SubcellTopLevel;
            }
            else{
                tempInteractionType = FmmInteractionType::M2L_SupercellTopLevel;
            }
            getInteractionListRotations(numBoxIndex, numLevel, tempInteractionType);

            m2l(numBoxIndex, numLevel, tempToggle, 1);

            // L2L
            if (maxLevel + numLevelSC >= 2){

                for (numLevel = numLevelSC - 2; numLevel >= 0; numLevel--){
                    numBoxIndex = variables->levelOffset[maxLevel + numLevel + 1] - variables->levelOffset[maxLevel + numLevel];

                    l2l(numBoxIndex, numLevel, 1);

                    getBoxIndexMask(numBoxIndex, numLevel, 1);

                    // M2L at lower levels
                    getInteractionListRotations(numBoxIndex, numLevel, FmmInteractionType::M2L_Supercell);

                    m2l(numBoxIndex, numLevel, 1);
                }

                for (numLevel = std::max(1, 3 - numLevelSC); numLevel <= maxLevel; numLevel++){

                    if (numLevel == 1) numBoxIndex = variables->levelOffset[maxLevel] - variables->levelOffset[0];
                    else numBoxIndex = variables->levelOffset[numLevel - 2] - variables->levelOffset[numLevel - 1];

                    l2l(numBoxIndex, numLevel, 0);

                    getBoxIndexMask(numBoxIndex, numLevel, 0);

                    // M2L at lower levels

                    getInteractionListRotations(numBoxIndex, numLevel, FmmInteractionType::M2L_Subcell);

                    m2l(numBoxIndex, numLevel, 0);
                }

            }

            // L2P
            if (maxLevel == 0) {
                l2p(numBoxIndex, 1);
            } else {
                l2p(numBoxIndex, 0);
            }
        }

        unsortParticles();
    }
    
    timing.total_time = timing.total_timer.getTimer();
    timing.writeTiming(timeInfoFile);
}

void MultipoleMethodForceCalculation::setDomainSizeRotations(){
    int i, j, nx, ny, nz;
    vec3<float> minOutliers, maxOutliers, dims;
    vec3<int> dimInd;
    int avgPPB;
    DomainToggle domainToggle = DomainToggle::cubicCells;
    vec3<int> tempBoxDimVals;
    
    const float inboundCheck = static_cast<float> (ParticleStatus::inbound);
    const float affectsAccelerationCheck = static_cast<float> (ParticleStatus::affectsAcceleration);
    
    DomainBoundsFunction::computeDomainBoundsRotations(particleSystem->chargedParticles, 
            variables->outlierCheck, tempBoxDimVals, boxMin, minOutliers, maxOutliers,
            dims, numOutliers, forceBoxDims, forcedBoxDimVals, domainToggle,
            avgPPB, dimInd);
    
    rootBoxSize = dims.x;
    
    boxDimVals[0] = tempBoxDimVals[0];
    boxDimVals[1] = tempBoxDimVals[1];
    boxDimVals[2] = tempBoxDimVals[2];
    numDomains = boxDimVals[0] * boxDimVals[1] * boxDimVals[2];
    
    // Not sure if this would ever lead to rounding errors, so below is safer
    //numLevelSC = ceil(log(boxDimVals[dimInd[2]]) / log(3)) + 1;

    int powThree = 1;
    numLevelSC = 1;
    while (boxDimVals[dimInd[2]] > powThree) {
        numLevelSC++;
        powThree *= 3;
    }

    variables->allocateDomainVariables(numLevelSC, numDomains, avgPPB);
    variables->numDomainsLevel[0] = numDomains;
    float boxSizeSC = rootBoxSize;
    int temp = 3;
    variables->dimValsSC[0].x = boxDimVals[0];
    variables->dimValsSC[0].y = boxDimVals[1];
    variables->dimValsSC[0].z = boxDimVals[2];
    variables->boxMinSC[0].x = boxMin.x;
    variables->boxMinSC[0].y = boxMin.y;
    variables->boxMinSC[0].z = boxMin.z;

    variables->boxMinIndexSC[0].x = 0;
    variables->boxMinIndexSC[0].y = 0;
    variables->boxMinIndexSC[0].z = 0;
    variables->offsetSC[0].x = (variables->dimValsSC[0].x % 3) % 2;
    variables->offsetSC[0].y = (variables->dimValsSC[0].y % 3) % 2;
    variables->offsetSC[0].z = (variables->dimValsSC[0].z % 3) % 2;

    // Calculate domain offsets and dimensions for higher level supercells
    for (i = 1; i < numLevelSC; i++) {
        variables->dimValsSC[i].x = (variables->dimValsSC[i - 1].x + 2) / 3;
        variables->dimValsSC[i].y = (variables->dimValsSC[i - 1].y + 2) / 3;
        variables->dimValsSC[i].z = (variables->dimValsSC[i - 1].z + 2) / 3;
        variables->numDomainsLevel[i] = variables->dimValsSC[i].x * variables->dimValsSC[i].y * variables->dimValsSC[i].z;
        variables->offsetSC[i].x = (variables->dimValsSC[i].x % 3) % 2;
        variables->offsetSC[i].y = (variables->dimValsSC[i].y % 3) % 2;
        variables->offsetSC[i].z = (variables->dimValsSC[i].z % 3) % 2;
        variables->boxMinSC[i].x = variables->boxMinSC[i - 1].x - variables->offsetSC[i - 1].x * boxSizeSC;
        variables->boxMinSC[i].y = variables->boxMinSC[i - 1].y - variables->offsetSC[i - 1].y * boxSizeSC;
        variables->boxMinSC[i].z = variables->boxMinSC[i - 1].z - variables->offsetSC[i - 1].z * boxSizeSC;
        boxSizeSC *= 3;

        variables->boxMinIndexSC[i].x = variables->boxMinIndexSC[i - 1].x - temp * variables->offsetSC[i].x;
        variables->boxMinIndexSC[i].y = variables->boxMinIndexSC[i - 1].y - temp * variables->offsetSC[i].y;
        variables->boxMinIndexSC[i].z = variables->boxMinIndexSC[i - 1].z - temp * variables->offsetSC[i].z;
        temp *= 3;
    }
    
    createDomainMortonLink();
    
    int tempIndex;

    vec3<int> tempDomainIndex = {-1, -1, -1};
    int tempMorton3;

    numOutliers = 0;
    unsigned int tempSignedMorton3;
    vec3<int> tempOutlierDomainIndex = variables->offsetSC[0];
    float tempPos;
    
    int numParticles = lastParticle - firstParticle;

    // Assign each particle to a certain domain
    for (i = 0; i < numParticles; i++) {
//        if (i == 113){
//            printf("test\n");
//        }
        if (particleSystem->chargedParticles.acceleration[i + firstParticle].w == inboundCheck ||
                particleSystem->chargedParticles.acceleration[i + firstParticle].w == affectsAccelerationCheck){
        
            tempPos = (particleSystem->chargedParticles.position[i + firstParticle].x - boxMin.x) / rootBoxSize;
            nx = int(tempPos);
            if (tempPos < 0) {
                nx--;
            }
            tempPos = (particleSystem->chargedParticles.position[i + firstParticle].y - boxMin.y) / rootBoxSize;
            ny = int(tempPos);
            if (tempPos < 0) {
                ny--;
            }
            tempPos = (particleSystem->chargedParticles.position[i + firstParticle].z - boxMin.z) / rootBoxSize;
            nz = int(tempPos);
            if (tempPos < 0) {
                nz--;
            }
            if (rootBoxSize == 0){
                nx = 0;
                ny = 0;
                nz = 0;
            }

            //        nx = int((particles.pos[i + initial].x - boxMin.x) / rootBoxSize);
            //        ny = int((particles.pos[i + initial].y - boxMin.y) / rootBoxSize);
            //        nz = int((particles.pos[i + initial].z - boxMin.z) / rootBoxSize);

            if (nx >= boxDimVals[0] || nx < 0 || ny >= boxDimVals[1] || ny < 0 || nz >= boxDimVals[2] || nz < 0) {
                nx += variables->offsetSC[0].x;
                ny += variables->offsetSC[0].y;
                nz += variables->offsetSC[0].z;
                if (nx != tempOutlierDomainIndex.x || ny != tempOutlierDomainIndex.y || nz != tempOutlierDomainIndex.z) {
                    tempOutlierDomainIndex.x = nx;
                    tempOutlierDomainIndex.y = ny;
                    tempOutlierDomainIndex.z = nz;
                    MultipoleMethodFunction::signedMorton3(tempOutlierDomainIndex, tempSignedMorton3);
                }
                variables->outlierList[numOutliers].index = i + firstParticle;
                variables->outlierList[numOutliers].signedMorton = tempSignedMorton3;
                numOutliers++;
            } else {
                if (nx >= boxDimVals[0]) nx--;
                if (ny >= boxDimVals[1]) ny--;
                if (nz >= boxDimVals[2]) nz--;
                nx += variables->offsetSC[0].x;
                ny += variables->offsetSC[0].y;
                nz += variables->offsetSC[0].z;
                if (nx != tempDomainIndex.x || ny != tempDomainIndex.y || nz != tempDomainIndex.z) {
                    tempDomainIndex.x = nx;
                    tempDomainIndex.y = ny;
                    tempDomainIndex.z = nz;
                    MultipoleMethodFunction::morton3(tempDomainIndex, tempMorton3, numLevelSC);
                    tempIndex = findMortonIndex(tempMorton3, 0);
                }

                if (variables->tempParticlesPerDomain[tempIndex][1] >= variables->tempParticlesPerDomain[tempIndex][0]) {
                    variables->reallocateTempParticles(tempIndex, numParticles - i);
                }
                variables->tempParticles[tempIndex][variables->tempParticlesPerDomain[tempIndex][1]] = i + firstParticle;
                variables->tempParticlesPerDomain[tempIndex][1]++;
            }
        }
        else{
            // If not in bounds, place in dummy domain
            if (variables->tempParticlesPerDomain[numDomains][1] >= variables->tempParticlesPerDomain[numDomains][0]) {
                variables->reallocateTempParticles(numDomains, numParticles - i);
            }
            variables->tempParticles[numDomains][variables->tempParticlesPerDomain[numDomains][1]] = i + firstParticle;
            variables->tempParticlesPerDomain[numDomains][1]++;
        }
    }

    tempIndex = firstParticle;
    for (i = 0; i < numDomains; i++) {
        variables->domainLinkage[0][i].firstParticle = tempIndex;
        tempIndex += variables->tempParticlesPerDomain[i][1];
        variables->domainLinkage[0][i].lastParticle = tempIndex - 1;
    }

//    printf("xdims = %d, ydims = %d, zdims = %d\n", boxDimVals[0], boxDimVals[1], boxDimVals[2]);
//    printf("xdiff = %.9g, ydiff = %.9g, zdiff = %.9g\n", xdiff, ydiff, zdiff);
//    if (numLevelSC > 2) {
//        printf("test\n");
//    }
//    if (numDomains == 3) {
//        printf("test\n");
//    }
//        if (numDomains > 1){
//            printf("test\n");
//        }
    //    if (numDomains > 4){
    //        printf("test\n");
    //    }
}

// This function should definitely have some parameter tuning to find optimum values for level_switch
void MultipoleMethodForceCalculation::setOptimumLevelRotations(){
    //    float level_switch[7]={100,2e4,1.7e5,1.3e6,1e7,7e7,5e8}; // cpu-tree
    //    float level_switch[8]={10,100,1.3e4,1e5,7e5,5e6,3e7,1.5e8}; // cpu-fmm
    //    float level_switch[7]={minAvgPPB,1.3e4,1e5,7e5,5e6,3e7,1.5e8}; // cpu-fmm
    //    float level_switch[7]={1000,1e5,5e5,5e6,3e7,2e8,1.5e9}; // gpu-tree
    //    float level_switch[7] = {1000,1e5,7e5,7e6,5e7,3e8,2e9}; // gpu-fmm

    float level_switch[8] = {static_cast<float>(2 * FMMConstant::minAvgPPB), 1000, 1.3e4, 1e5, 7e5, 5e6, 3e7, 1.5e8};

    int domainMaxLevel, domainNumParticles;
    numBoxIndexFull = 0;
    maxLevel = 0;
    for (int i = 0; i < numDomains; i++) {
        domainMaxLevel = 0;
        domainNumParticles = variables->domainLinkage[0][i].lastParticle - variables->domainLinkage[0][i].firstParticle + 1;
        if (domainNumParticles < level_switch[0]) {
            domainMaxLevel = 0;
        } else if (domainNumParticles < level_switch[1]) {
            domainMaxLevel += 1;
        } else if (domainNumParticles < level_switch[2]) {
            domainMaxLevel += 2;
        } else if (domainNumParticles < level_switch[3]) {
            domainMaxLevel += 3;
        } else if (domainNumParticles < level_switch[4]) {
            domainMaxLevel += 4;
        } else if (domainNumParticles < level_switch[5]) {
            domainMaxLevel += 5;
        } else if (domainNumParticles < level_switch[6]) {
            domainMaxLevel += 6;
        } else if (domainNumParticles < level_switch[7]) {
            domainMaxLevel += 7;
        } else {
            domainMaxLevel += 8;
        }
        maxLevel = std::max(maxLevel, domainMaxLevel);
    }
    
    
    if (forcedMaxLevel != -2){
        maxLevel = forcedMaxLevel;
    }
    

    domainNumBoxIndexFull = 1 << (3 * maxLevel);
    numBoxIndexFull = numDomains * domainNumBoxIndexFull;
//    printf("level   : %d\n", maxLevel);
}

void MultipoleMethodForceCalculation::getInteractionListRotations(int numBoxIndex, int numLevel, FmmInteractionType interactionType){
    int jxmin, jxmax, jymin, jymax, jzmin, jzmax, ii, ib, jj, jb, ix, iy, iz, jx, jy, jz, boxIndex;
    int ixp, iyp, izp, jxp, jyp, jzp;
    vec3<int> boxIndex3D, domainIndex3D;
    int i, j, jdiff, tempQuad, halfDiff, jtemp;
    
    if (interactionType == FmmInteractionType::M2L_Subcell || 
            interactionType == FmmInteractionType::M2L_SubcellOutliers ||
            interactionType == FmmInteractionType::M2L_SubcellTopLevel ||
            interactionType == FmmInteractionType::M2L_SubcellTopLevelOutliers ||
            interactionType == FmmInteractionType::P2P ||
            interactionType == FmmInteractionType::P2P_Outliers){
        jdiff = (1 << numLevel);
    }
    else if (interactionType == FmmInteractionType::M2L_Supercell ||
            interactionType == FmmInteractionType::M2L_SupercellOutliers ||
            interactionType == FmmInteractionType::M2L_SupercellTopLevel ||
            interactionType == FmmInteractionType::M2L_SupercellTopLevelOutliers ||
            interactionType == FmmInteractionType::P2P_Supercell ||
            interactionType == FmmInteractionType::P2P_SupercellOutliers){
        jdiff = pow(3, numLevel);
    }
    halfDiff = (jdiff >> 1);

    int tempDomainNumBoxIndexFull = 1 << (3 * numLevel);

    int domainBias, domainBoxIndex, guessIndex;
    vec3<int> tempDomainIndex3D;

    // P2P
    if (interactionType == FmmInteractionType::P2P) {
        // Loop through every domain
        for (i = 0; i < numDomains; i++) {
            MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][i].mortonIndex, domainIndex3D);
            // Loop through all cells within domain i
            for (ii = variables->domainOffset[0][i][numLevel]; ii <= variables->domainOffset[1][i][numLevel]; ii++) {
                ib = ii + variables->levelOffset[numLevel - 1];
                variables->numInteraction[ii] = 0;

                // Find local box index
                MultipoleMethodFunction::unmorton(variables->domainBoxIndexFull[ib], boxIndex3D);
                ix = boxIndex3D.x;
                iy = boxIndex3D.y;
                iz = boxIndex3D.z;
                
                for (jx = ix - 1; jx <= ix + 1; jx++){
                    for (jy = iy - 1; jy <= iy + 1; jy++){
                        for (jz = iz - 1; jz <= iz + 1; jz++){
                            boxIndex3D.x = jx;
                            boxIndex3D.y = jy;
                            boxIndex3D.z = jz;
                            tempDomainIndex3D = domainIndex3D;
                            
                            if (jx >= halfDiff && jy >= halfDiff){
                                tempQuad = 0;
                            }
                            else if (jx < halfDiff && jy >= halfDiff){
                                tempQuad = 1;
                                boxIndex3D.x = jy;
                                boxIndex3D.y = jdiff - jx - 1;
                            }
                            else if (jx < halfDiff && jy < halfDiff){
                                tempQuad = 2;
                                boxIndex3D.x = jdiff - jx - 1;
                                boxIndex3D.y = jdiff - jy - 1;
                            }
                            else if (jx >= halfDiff && jy < halfDiff){
                                tempQuad = 3;
                                boxIndex3D.x = jdiff - jy - 1;
                                boxIndex3D.y = jx;
                            }
                            
                            if (jz < 0) {
                                do {
                                    boxIndex3D.z += jdiff;
                                    tempDomainIndex3D.z--;
                                } while (boxIndex3D.z < 0);
                            } else if (jz >= jdiff) {
                                do {
                                    boxIndex3D.z -= jdiff;
                                    tempDomainIndex3D.z++;
                                } while (boxIndex3D.z >= jdiff);
                            }
                            
                            // Check to see if domain is valid and cell contains any particles
                            if (boxIndex3D.x >= 0 && boxIndex3D.x < jdiff && boxIndex3D.y >= 0 && boxIndex3D.y < jdiff
                                    && tempDomainIndex3D.z >= variables->offsetSC[0].z && tempDomainIndex3D.z < boxDimVals[2] + variables->offsetSC[0].z) {
                                MultipoleMethodFunction::morton1(boxIndex3D, boxIndex, numLevel);
                                if (tempDomainIndex3D.x == domainIndex3D.x && tempDomainIndex3D.y == domainIndex3D.y && tempDomainIndex3D.z == domainIndex3D.z) {
                                    domainBias = i * tempDomainNumBoxIndexFull;
                                } else {
                                    MultipoleMethodFunction::morton3(tempDomainIndex3D, domainBoxIndex, numLevelSC);
                                    guessIndex = findMortonIndex(domainBoxIndex, 0);
                                    domainBias = guessIndex * tempDomainNumBoxIndexFull;
                                }
                                jj = variables->domainBoxIndexMask[boxIndex + domainBias];
                                if (jj != -1) {
                                    variables->interactionList[ii][variables->numInteraction[ii]][0] = jj;
                                    variables->interactionList[ii][variables->numInteraction[ii]][1] = tempQuad;
                                    variables->numInteraction[ii]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }// P2P at domain level
    else if (interactionType == FmmInteractionType::P2P_Supercell) {
        // Loop through every domain
        //for (ii = 0; ii < numDomains; ii++){
        for (ii = 0; ii < numBoxIndex; ii++) {
            ib = ii + variables->levelOffset[numLevel + maxLevel];
            variables->numInteraction[ii] = 0;
            //MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][ii].mortonIndex, boxIndex3D);
            MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][variables->domainBoxIndexFull[ib]].mortonIndex, boxIndex3D);
            ix = boxIndex3D.x;
            iy = boxIndex3D.y;
            iz = boxIndex3D.z;
            // Check cells within a 3x3x3 box centered around box being considered
            for (jx = std::max(ix - 1, variables->offsetSC[0].x); jx <= std::min(ix + 1, variables->offsetSC[0].x + boxDimVals[0] - 1); jx++) {
                for (jy = std::max(iy - 1, variables->offsetSC[0].y); jy <= std::min(iy + 1, variables->offsetSC[0].y + boxDimVals[1] - 1); jy++) {
                    for (jz = std::max(iz - 1, variables->offsetSC[0].z); jz <= std::min(iz + 1, variables->offsetSC[0].z + boxDimVals[2] - 1); jz++) {
                        boxIndex3D.x = jx;
                        boxIndex3D.y = jy;
                        boxIndex3D.z = jz;

                        MultipoleMethodFunction::morton3(boxIndex3D, domainBoxIndex, numLevelSC - numLevel);
                        guessIndex = findMortonIndex(domainBoxIndex, numLevel);
                        guessIndex = variables->domainBoxIndexMask[guessIndex];
                        if (guessIndex != -1) {
                            variables->interactionList[ii][variables->numInteraction[ii]][0] = guessIndex;
                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 0; // Supercells don't have rotational symmetry
                            variables->numInteraction[ii]++;
                        }
                    }
                }
            }
        }
    }// P2P for outliers
    else if (interactionType == FmmInteractionType::P2P_Outliers) {
        // Loop through every domain
        for (i = 0; i < numOutlierDomains; i++) {
            MultipoleMethodFunction::signedUnmorton3(variables->outlierDomainLinkage[i].mortonIndex, domainIndex3D);
            // Loop through all cells within domain i
            for (ii = variables->outlierDomainOffset[0][i][numLevel]; ii <= variables->outlierDomainOffset[1][i][numLevel]; ii++) {
                variables->numInteraction[ii] = 0;

                // Find local box index
                MultipoleMethodFunction::unmorton(variables->outlierBoxIndexFull[ii], boxIndex3D);
                ix = boxIndex3D.x;
                iy = boxIndex3D.y;
                iz = boxIndex3D.z;

                // Check cells within a 3x3x3 box centered around box being considered
                for (jx = ix - 1; jx <= ix + 1; jx++) {
                    for (jy = iy - 1; jy <= iy + 1; jy++) {
                        for (jz = iz - 1; jz <= iz + 1; jz++) {
                            boxIndex3D.x = jx;
                            boxIndex3D.y = jy;
                            boxIndex3D.z = jz;
                            tempDomainIndex3D = domainIndex3D;
                            
                            // Shift cell to lie within some domain
                            if (jx < 0) {
                                do {
                                    boxIndex3D.x += jdiff;
                                    tempDomainIndex3D.x--;
                                } while (boxIndex3D.x < 0);
                            } else if (jx >= jdiff) {
                                do {
                                    boxIndex3D.x -= jdiff;
                                    tempDomainIndex3D.x++;
                                } while (boxIndex3D.x >= jdiff);
                            }
                            if (jy < 0) {
                                do {
                                    boxIndex3D.y += jdiff;
                                    tempDomainIndex3D.y--;
                                } while (boxIndex3D.y < 0);
                            } else if (jy >= jdiff) {
                                do {
                                    boxIndex3D.y -= jdiff;
                                    tempDomainIndex3D.y++;
                                } while (boxIndex3D.y >= jdiff);
                            }
                            
                            if (boxIndex3D.x >= halfDiff && boxIndex3D.y >= halfDiff){
                                tempQuad = 0;
                            }
                            else if (boxIndex3D.x < halfDiff && boxIndex3D.y >= halfDiff){
                                tempQuad = 1;
                                jtemp = boxIndex3D.x;
                                boxIndex3D.x = boxIndex3D.y;
                                boxIndex3D.y = jdiff - jtemp - 1;
                            }
                            else if (boxIndex3D.x < halfDiff && boxIndex3D.y < halfDiff){
                                tempQuad = 2;
                                boxIndex3D.x = jdiff - boxIndex3D.x - 1;
                                boxIndex3D.y = jdiff - boxIndex3D.y - 1;
                            }
                            else if (boxIndex3D.x >= halfDiff && boxIndex3D.y < halfDiff){
                                tempQuad = 3;
                                jtemp = boxIndex3D.x;
                                boxIndex3D.x = jdiff - boxIndex3D.y - 1;
                                boxIndex3D.y = jtemp;
                            }
                            
                            if (jz < 0) {
                                do {
                                    boxIndex3D.z += jdiff;
                                    tempDomainIndex3D.z--;
                                } while (boxIndex3D.z < 0);
                            } else if (jz >= jdiff) {
                                do {
                                    boxIndex3D.z -= jdiff;
                                    tempDomainIndex3D.z++;
                                } while (boxIndex3D.z >= jdiff);
                            }

                            // Check to see if domain is valid and cell contains any particles
                            if (boxIndex3D.x >= 0 && boxIndex3D.x < jdiff && boxIndex3D.y >= 0 && boxIndex3D.y < jdiff
                                    && tempDomainIndex3D.x >= variables->offsetSC[0].x && tempDomainIndex3D.x < boxDimVals[0] + variables->offsetSC[0].x
                                    && tempDomainIndex3D.y >= variables->offsetSC[0].y && tempDomainIndex3D.y < boxDimVals[1] + variables->offsetSC[0].y
                                    && tempDomainIndex3D.z >= variables->offsetSC[0].z && tempDomainIndex3D.z < boxDimVals[2] + variables->offsetSC[0].z) {
                                MultipoleMethodFunction::morton1(boxIndex3D, boxIndex, numLevel);
                                if (tempDomainIndex3D.x == domainIndex3D.x && tempDomainIndex3D.y == domainIndex3D.y && tempDomainIndex3D.z == domainIndex3D.z) {
                                    domainBias = i * tempDomainNumBoxIndexFull;
                                } else {
                                    MultipoleMethodFunction::morton3(tempDomainIndex3D, domainBoxIndex, numLevelSC);
                                    guessIndex = findMortonIndex(domainBoxIndex, 0);
                                    domainBias = guessIndex * tempDomainNumBoxIndexFull;
                                }
                                jj = variables->domainBoxIndexMask[boxIndex + domainBias];
                                if (jj != -1) {
                                    variables->interactionList[ii][variables->numInteraction[ii]][0] = jj;
                                    variables->interactionList[ii][variables->numInteraction[ii]][1] = tempQuad;
                                    variables->numInteraction[ii]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }// P2P at domain level for outliers
    else if (interactionType == FmmInteractionType::P2P_SupercellOutliers) {
        // Loop through every domain
        //for (ii = 0; ii < numDomains; ii++){
        for (ii = 0; ii < numBoxIndex; ii++) {
            variables->numInteraction[ii] = 0;
            //MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][ii].mortonIndex, boxIndex3D);
            MultipoleMethodFunction::signedUnmorton3(variables->outlierDomainLinkage[ii].mortonIndex, boxIndex3D);
            ix = boxIndex3D.x;
            iy = boxIndex3D.y;
            iz = boxIndex3D.z;
            // Check cells within a 3x3x3 box centered around box being considered
            for (jx = std::max(ix - 1, variables->offsetSC[0].x); jx <= std::min(ix + 1, variables->offsetSC[0].x + boxDimVals[0] - 1); jx++) {
                for (jy = std::max(iy - 1, variables->offsetSC[0].y); jy <= std::min(iy + 1, variables->offsetSC[0].y + boxDimVals[1] - 1); jy++) {
                    for (jz = std::max(iz - 1, variables->offsetSC[0].z); jz <= std::min(iz + 1, variables->offsetSC[0].z + boxDimVals[2] - 1); jz++) {
                        boxIndex3D.x = jx;
                        boxIndex3D.y = jy;
                        boxIndex3D.z = jz;

                        MultipoleMethodFunction::morton3(boxIndex3D, domainBoxIndex, numLevelSC - numLevel);
                        guessIndex = findMortonIndex(domainBoxIndex, numLevel);
                        guessIndex = variables->domainBoxIndexMask[guessIndex];
                        if (guessIndex != -1) {
                            variables->interactionList[ii][variables->numInteraction[ii]][0] = guessIndex;
                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 0;
                            variables->numInteraction[ii]++;
                        }
                    }
                }
            }
        }
    }// M2L at level 1
    else if (interactionType == FmmInteractionType::M2L_SubcellTopLevel) {
        // Loop through every domain
        for (i = 0; i < numDomains; i++) {
            MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][i].mortonIndex, domainIndex3D);
            // Loop through all cells within domain i
            for (ii = variables->domainOffset[0][i][numLevel]; ii <= variables->domainOffset[1][i][numLevel]; ii++) {
                ib = ii + variables->levelOffset[numLevel - 1];
                variables->numInteraction[ii] = 0;

                // Find local box index
                MultipoleMethodFunction::unmorton(variables->domainBoxIndexFull[ib], boxIndex3D);
                ix = boxIndex3D.x;
                iy = boxIndex3D.y;
                iz = boxIndex3D.z;

                // Loop through every domain
                for (j = 0; j < numDomains; j++) {
                    MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][j].mortonIndex, tempDomainIndex3D);
                    domainBias = j * (1 << (3 * numLevel));

                    // Loop through all cells within domain j
                    for (jj = variables->domainOffset[0][j][numLevel]; jj <= variables->domainOffset[1][j][numLevel]; jj++) {
                        jb = jj + variables->levelOffset[numLevel - 1];
                        MultipoleMethodFunction::unmorton(variables->domainBoxIndexFull[jb], boxIndex3D);
                        jx = boxIndex3D.x;
                        jy = boxIndex3D.y;
                        jz = boxIndex3D.z;
                        jx += jdiff * (tempDomainIndex3D.x - domainIndex3D.x);  // Should add zero
                        jy += jdiff * (tempDomainIndex3D.y - domainIndex3D.y);  // Should add zero
                        jz += jdiff * (tempDomainIndex3D.z - domainIndex3D.z);

                        // Check if jj-th cell lies within 3x3x3 box of ii-th cell
                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                            variables->interactionList[ii][variables->numInteraction[ii]][0] = variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[jb]];
                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 0;
                            variables->numInteraction[ii]++;
                        }
                        
                        // Second quadrant
                        jy = boxIndex3D.x;
                        jx = jdiff - boxIndex3D.y - 1;
                        jx += jdiff * (tempDomainIndex3D.x - domainIndex3D.x);
                        jy += jdiff * (tempDomainIndex3D.y - domainIndex3D.y);
                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                            variables->interactionList[ii][variables->numInteraction[ii]][0] = variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[jb]];
                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 1;
                            variables->numInteraction[ii]++;
                        }
                        
                        // Third quadrant
                        jx = jdiff - boxIndex3D.x - 1;
                        jy = jdiff - boxIndex3D.y - 1;
                        jx += jdiff * (tempDomainIndex3D.x - domainIndex3D.x);
                        jy += jdiff * (tempDomainIndex3D.y - domainIndex3D.y);
                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                            variables->interactionList[ii][variables->numInteraction[ii]][0] = variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[jb]];
                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 2;
                            variables->numInteraction[ii]++;
                        }
                        
                        // Fourth quadrant
                        jx = boxIndex3D.y;
                        jy = jdiff - boxIndex3D.x - 1;
                        jx += jdiff * (tempDomainIndex3D.x - domainIndex3D.x);
                        jy += jdiff * (tempDomainIndex3D.y - domainIndex3D.y);
                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                            variables->interactionList[ii][variables->numInteraction[ii]][0] = variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[jb]];
                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 3;
                            variables->numInteraction[ii]++;
                        }
                    }
                }
            }
        }
    }// M2L at level 1 for outliers
    else if (interactionType == FmmInteractionType::M2L_SubcellTopLevelOutliers) {
        // Loop through every domain
        for (i = 0; i < numOutlierDomains; i++) {
            MultipoleMethodFunction::signedUnmorton3(variables->outlierDomainLinkage[i].mortonIndex, domainIndex3D);
            // Loop through all cells within domain i
            for (ii = variables->outlierDomainOffset[0][i][numLevel]; ii <= variables->outlierDomainOffset[1][i][numLevel]; ii++) {
                variables->numInteraction[ii] = 0;

                // Find local box index
                MultipoleMethodFunction::unmorton(variables->outlierBoxIndexFull[ii], boxIndex3D);
                ix = boxIndex3D.x;
                iy = boxIndex3D.y;
                iz = boxIndex3D.z;

                // Loop through every domain
                for (j = 0; j < numDomains; j++) {
                    MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][j].mortonIndex, tempDomainIndex3D);
                    domainBias = j * (1 << (3 * numLevel));

                    // Loop through all cells within domain j
                    for (jj = variables->domainOffset[0][j][numLevel]; jj <= variables->domainOffset[1][j][numLevel]; jj++) {
                        jb = jj + variables->levelOffset[numLevel - 1];
                        MultipoleMethodFunction::unmorton(variables->domainBoxIndexFull[jb], boxIndex3D);
                        jx = boxIndex3D.x;
                        jy = boxIndex3D.y;
                        jz = boxIndex3D.z;
                        jx += jdiff * (tempDomainIndex3D.x - domainIndex3D.x);
                        jy += jdiff * (tempDomainIndex3D.y - domainIndex3D.y);
                        jz += jdiff * (tempDomainIndex3D.z - domainIndex3D.z);

                        // Check if jj-th cell lies within 3x3x3 box of ii-th cell
                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                            variables->interactionList[ii][variables->numInteraction[ii]][0] = variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[jb]];
                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 0;
                            variables->numInteraction[ii]++;
                        }
                        
                        // Second quadrant
                        jy = boxIndex3D.x;
                        jx = jdiff - boxIndex3D.y - 1;
                        jx += jdiff * (tempDomainIndex3D.x - domainIndex3D.x);
                        jy += jdiff * (tempDomainIndex3D.y - domainIndex3D.y);
                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                            variables->interactionList[ii][variables->numInteraction[ii]][0] = variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[jb]];
                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 1;
                            variables->numInteraction[ii]++;
                        }
                        
                        // Third quadrant
                        jx = jdiff - boxIndex3D.x - 1;
                        jy = jdiff - boxIndex3D.y - 1;
                        jx += jdiff * (tempDomainIndex3D.x - domainIndex3D.x);
                        jy += jdiff * (tempDomainIndex3D.y - domainIndex3D.y);
                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                            variables->interactionList[ii][variables->numInteraction[ii]][0] = variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[jb]];
                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 2;
                            variables->numInteraction[ii]++;
                        }
                        
                        // Fourth quadrant
                        jx = boxIndex3D.y;
                        jy = jdiff - boxIndex3D.x - 1;
                        jx += jdiff * (tempDomainIndex3D.x - domainIndex3D.x);
                        jy += jdiff * (tempDomainIndex3D.y - domainIndex3D.y);
                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                            variables->interactionList[ii][variables->numInteraction[ii]][0] = variables->domainBoxIndexMask[domainBias + variables->domainBoxIndexFull[jb]];
                            variables->interactionList[ii][variables->numInteraction[ii]][1] = 3;
                            variables->numInteraction[ii]++;
                        }
                    }
                }
            }
        }
    }// M2L at subcell level
    else if (interactionType == FmmInteractionType::M2L_Subcell) {
        // Loop through every domain
        for (i = 0; i < numDomains; i++) {
            MultipoleMethodFunction::unmorton3(variables->domainLinkage[0][i].mortonIndex, domainIndex3D);
            // Loop through all cells within domain i
            for (ii = variables->domainOffset[0][i][numLevel]; ii <= variables->domainOffset[1][i][numLevel]; ii++) {
                ib = ii + variables->levelOffset[numLevel - 1];
                variables->numInteraction[ii] = 0;

                // Find local box index
                MultipoleMethodFunction::unmorton(variables->domainBoxIndexFull[ib], boxIndex3D);
                ix = boxIndex3D.x;
                iy = boxIndex3D.y;
                iz = boxIndex3D.z;
                ixp = (ix + 2) / 2;
                iyp = (iy + 2) / 2;
                izp = (iz + 2) / 2;

                // Loop through all cells within a 3x3x3 box centered around the parent cell
                for (jxp = ixp - 1; jxp <= ixp + 1; jxp++) {
                    for (jyp = iyp - 1; jyp <= iyp + 1; jyp++) {
                        for (jzp = izp - 1; jzp <= izp + 1; jzp++) {
                            for (jx = 2 * jxp - 2; jx <= 2 * jxp - 1; jx++) {
                                for (jy = 2 * jyp - 2; jy <= 2 * jyp - 1; jy++) {
                                    for (jz = 2 * jzp - 2; jz <= 2 * jzp - 1; jz++) {
                                        // Cells must be outside of a 3x3x3 box centered around ii-th cell
                                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                                            boxIndex3D.x = jx;
                                            boxIndex3D.y = jy;
                                            boxIndex3D.z = jz;
                                            tempDomainIndex3D = domainIndex3D;

                                            if (jx >= halfDiff && jy >= halfDiff){
                                                tempQuad = 0;
                                            }
                                            else if (jx < halfDiff && jy >= halfDiff){
                                                tempQuad = 1;
                                                boxIndex3D.x = jy;
                                                boxIndex3D.y = jdiff - jx - 1;
                                            }
                                            else if (jx < halfDiff && jy < halfDiff){
                                                tempQuad = 2;
                                                boxIndex3D.x = jdiff - jx - 1;
                                                boxIndex3D.y = jdiff - jy - 1;
                                            }
                                            else if (jx >= halfDiff && jy < halfDiff){
                                                tempQuad = 3;
                                                boxIndex3D.x = jdiff - jy - 1;
                                                boxIndex3D.y = jx;
                                            }
                                                                                        
                                            if (jz < 0) {
                                                do {
                                                    boxIndex3D.z += jdiff;
                                                    tempDomainIndex3D.z--;
                                                } while (boxIndex3D.z < 0);
                                            } else if (jz >= jdiff) {
                                                do {
                                                    boxIndex3D.z -= jdiff;
                                                    tempDomainIndex3D.z++;
                                                } while (boxIndex3D.z >= jdiff);
                                            }

                                            // Check to see if domain is valid and cell contains any particles
                                            if (boxIndex3D.x >= 0 && boxIndex3D.x < jdiff && boxIndex3D.y >= 0 && boxIndex3D.y < jdiff
                                                    && tempDomainIndex3D.z >= variables->offsetSC[0].z && tempDomainIndex3D.z < boxDimVals[2] + variables->offsetSC[0].z) {
                                                MultipoleMethodFunction::morton1(boxIndex3D, boxIndex, numLevel);
                                                if (tempDomainIndex3D.x == domainIndex3D.x && tempDomainIndex3D.y == domainIndex3D.y && tempDomainIndex3D.z == domainIndex3D.z) {
                                                    domainBias = i * tempDomainNumBoxIndexFull;
                                                } else {
                                                    MultipoleMethodFunction::morton3(tempDomainIndex3D, domainBoxIndex, numLevelSC);
                                                    guessIndex = findMortonIndex(domainBoxIndex, 0);
                                                    domainBias = guessIndex * tempDomainNumBoxIndexFull;
                                                }
                                                jj = variables->domainBoxIndexMask[boxIndex + domainBias];
                                                if (jj != -1) {
                                                    variables->interactionList[ii][variables->numInteraction[ii]][0] = jj;
                                                    variables->interactionList[ii][variables->numInteraction[ii]][1] = tempQuad;
                                                    variables->numInteraction[ii]++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }// M2L at subcell level for outliers
    else if (interactionType == FmmInteractionType::M2L_SubcellOutliers) {
        // Loop through every domain
        for (i = 0; i < numOutlierDomains; i++) {
            MultipoleMethodFunction::signedUnmorton3(variables->outlierDomainLinkage[i].mortonIndex, domainIndex3D);
            // Loop through all cells within domain i
            for (ii = variables->outlierDomainOffset[0][i][numLevel]; ii <= variables->outlierDomainOffset[1][i][numLevel]; ii++) {
                variables->numInteraction[ii] = 0;

                // Find local box index
                MultipoleMethodFunction::unmorton(variables->outlierBoxIndexFull[ii], boxIndex3D);
                ix = boxIndex3D.x;
                iy = boxIndex3D.y;
                iz = boxIndex3D.z;
                ixp = (ix + 2) / 2;
                iyp = (iy + 2) / 2;
                izp = (iz + 2) / 2;

                // Loop through all cells within a 3x3x3 box centered around the parent cell
                for (jxp = ixp - 1; jxp <= ixp + 1; jxp++) {
                    for (jyp = iyp - 1; jyp <= iyp + 1; jyp++) {
                        for (jzp = izp - 1; jzp <= izp + 1; jzp++) {
                            for (jx = 2 * jxp - 2; jx <= 2 * jxp - 1; jx++) {
                                for (jy = 2 * jyp - 2; jy <= 2 * jyp - 1; jy++) {
                                    for (jz = 2 * jzp - 2; jz <= 2 * jzp - 1; jz++) {
                                        // Cells must be outside of a 3x3x3 box centered around ii-th cell
                                        if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                                            boxIndex3D.x = jx;
                                            boxIndex3D.y = jy;
                                            boxIndex3D.z = jz;
                                            tempDomainIndex3D = domainIndex3D;
                                            
                                            // Shift cell to lie within some domain
                                            if (jx < 0) {
                                                do {
                                                    boxIndex3D.x += jdiff;
                                                    tempDomainIndex3D.x--;
                                                } while (boxIndex3D.x < 0);
                                            } else if (jx >= jdiff) {
                                                do {
                                                    boxIndex3D.x -= jdiff;
                                                    tempDomainIndex3D.x++;
                                                } while (boxIndex3D.x >= jdiff);
                                            }
                                            if (jy < 0) {
                                                do {
                                                    boxIndex3D.y += jdiff;
                                                    tempDomainIndex3D.y--;
                                                } while (boxIndex3D.y < 0);
                                            } else if (jy >= jdiff) {
                                                do {
                                                    boxIndex3D.y -= jdiff;
                                                    tempDomainIndex3D.y++;
                                                } while (boxIndex3D.y >= jdiff);
                                            }

                                            if (boxIndex3D.x >= halfDiff && boxIndex3D.y >= halfDiff){
                                                tempQuad = 0;
                                            }
                                            else if (boxIndex3D.x < halfDiff && boxIndex3D.y >= halfDiff){
                                                tempQuad = 1;
                                                jtemp = boxIndex3D.x;
                                                boxIndex3D.x = boxIndex3D.y;
                                                boxIndex3D.y = jdiff - jtemp - 1;
                                            }
                                            else if (boxIndex3D.x < halfDiff && boxIndex3D.y < halfDiff){
                                                tempQuad = 2;
                                                boxIndex3D.x = jdiff - boxIndex3D.x - 1;
                                                boxIndex3D.y = jdiff - boxIndex3D.y - 1;
                                            }
                                            else if (boxIndex3D.x >= halfDiff && boxIndex3D.y < halfDiff){
                                                tempQuad = 3;
                                                jtemp = boxIndex3D.x;
                                                boxIndex3D.x = jdiff - boxIndex3D.y - 1;
                                                boxIndex3D.y = jtemp;
                                            }
                                            
                                            if (jz < 0) {
                                                do {
                                                    boxIndex3D.z += jdiff;
                                                    tempDomainIndex3D.z--;
                                                } while (boxIndex3D.z < 0);
                                            } else if (jz >= jdiff) {
                                                do {
                                                    boxIndex3D.z -= jdiff;
                                                    tempDomainIndex3D.z++;
                                                } while (boxIndex3D.z >= jdiff);
                                            }

                                            // Check to see if domain is valid and cell contains any particles
                                            if (boxIndex3D.x >= 0 && boxIndex3D.x < jdiff && boxIndex3D.y >= 0 && boxIndex3D.y < jdiff
                                                    && tempDomainIndex3D.x >= variables->offsetSC[0].x && tempDomainIndex3D.x < boxDimVals[0] + variables->offsetSC[0].x
                                                    && tempDomainIndex3D.y >= variables->offsetSC[0].y && tempDomainIndex3D.y < boxDimVals[1] + variables->offsetSC[0].y
                                                    && tempDomainIndex3D.z >= variables->offsetSC[0].z && tempDomainIndex3D.z < boxDimVals[2] + variables->offsetSC[0].z) {
                                                MultipoleMethodFunction::morton1(boxIndex3D, boxIndex, numLevel);
                                                if (tempDomainIndex3D.x == domainIndex3D.x && tempDomainIndex3D.y == domainIndex3D.y && tempDomainIndex3D.z == domainIndex3D.z) {
                                                    domainBias = i * tempDomainNumBoxIndexFull;
                                                } else {
                                                    MultipoleMethodFunction::morton3(tempDomainIndex3D, domainBoxIndex, numLevelSC);
                                                    guessIndex = findMortonIndex(domainBoxIndex, 0);
                                                    domainBias = guessIndex * tempDomainNumBoxIndexFull;
                                                }
                                                jj = variables->domainBoxIndexMask[boxIndex + domainBias];
                                                if (jj != -1) {
                                                    variables->interactionList[ii][variables->numInteraction[ii]][0] = jj;
                                                    variables->interactionList[ii][variables->numInteraction[ii]][1] = tempQuad;
                                                    variables->numInteraction[ii]++;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }// M2L for supercells
    else if (interactionType == FmmInteractionType::M2L_Supercell) {
        // Set bounds for supercell indices
        jxmin = variables->offsetSC[numLevel].x;
        jxmax = variables->offsetSC[numLevel].x + variables->dimValsSC[numLevel].x - 1;
        jymin = variables->offsetSC[numLevel].y;
        jymax = variables->offsetSC[numLevel].y + variables->dimValsSC[numLevel].y - 1;
        jzmin = variables->offsetSC[numLevel].z;
        jzmax = variables->offsetSC[numLevel].z + variables->dimValsSC[numLevel].z - 1;

        // Loop through every domain at level numLevel
        //for (ii = 0; ii < numDomainsLevel[numLevel]; ii++){
        for (ii = 0; ii < numBoxIndex; ii++) {
            ib = ii + variables->levelOffset[maxLevel + numLevel];
            variables->numInteraction[ii] = 0;

            // Find box index of domain
            //MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][ii].mortonIndex, boxIndex3D);
            MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].mortonIndex, boxIndex3D);
            ix = boxIndex3D.x;
            iy = boxIndex3D.y;
            iz = boxIndex3D.z;
            ixp = (ix + 3) / 3;
            iyp = (iy + 3) / 3;
            izp = (iz + 3) / 3;

            // Loop through all cells within a 3x3x3 box centered around the parent cell
            for (jxp = ixp - 1; jxp <= ixp + 1; jxp++) {
                for (jyp = iyp - 1; jyp <= iyp + 1; jyp++) {
                    for (jzp = izp - 1; jzp <= izp + 1; jzp++) {
                        for (jx = std::max(3 * jxp - 3, jxmin); jx <= std::min(3 * jxp - 1, jxmax); jx++) {
                            for (jy = std::max(3 * jyp - 3, jymin); jy <= std::min(3 * jyp - 1, jymax); jy++) {
                                for (jz = std::max(3 * jzp - 3, jzmin); jz <= std::min(3 * jzp - 1, jzmax); jz++) {
                                    // Cells must be outside of a 3x3x3 box centered around ii-th cell
                                    if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                                        boxIndex3D.x = jx;
                                        boxIndex3D.y = jy;
                                        boxIndex3D.z = jz;
                                        MultipoleMethodFunction::morton3(boxIndex3D, domainBoxIndex, numLevelSC - numLevel);
                                        guessIndex = findMortonIndex(domainBoxIndex, numLevel);
                                        guessIndex = variables->domainBoxIndexMask[guessIndex];
                                        if (guessIndex != -1) {
                                            variables->interactionList[ii][variables->numInteraction[ii]][0] = guessIndex;
                                            variables->numInteraction[ii]++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }// M2L at supercell level for outliers
    else if (interactionType == FmmInteractionType::M2L_SupercellOutliers) {
        // Set bounds for supercell indices
        jxmin = variables->offsetSC[numLevel].x;
        jxmax = variables->offsetSC[numLevel].x + variables->dimValsSC[numLevel].x - 1;
        jymin = variables->offsetSC[numLevel].y;
        jymax = variables->offsetSC[numLevel].y + variables->dimValsSC[numLevel].y - 1;
        jzmin = variables->offsetSC[numLevel].z;
        jzmax = variables->offsetSC[numLevel].z + variables->dimValsSC[numLevel].z - 1;

        // Loop through every domain at level numLevel
        //for (ii = 0; ii < numDomainsLevel[numLevel]; ii++){
        for (ii = 0; ii < numBoxIndex; ii++) {
            variables->numInteraction[ii] = 0;

            // Find box index of domain
            //MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][ii].mortonIndex, boxIndex3D);
            //            MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].mortonIndex, boxIndex3D);
            //            ix = boxIndex3D.x;
            //            iy = boxIndex3D.y;
            //            iz = boxIndex3D.z;
            MultipoleMethodFunction::signedUnmorton3(variables->outlierDomainLinkage[ii].mortonIndex, boxIndex3D);
            ix = boxIndex3D.x;
            iy = boxIndex3D.y;
            iz = boxIndex3D.z;

//            if (ix < variables->boxMinIndexSC[numLevel].x) {
            if (ix < variables->boxMinIndexSC[numLevel].x && ((ix - variables->boxMinIndexSC[numLevel].x) % jdiff != 0)) {
                ix = (ix - variables->boxMinIndexSC[numLevel].x) / jdiff - 1;
            } else {
                ix = (ix - variables->boxMinIndexSC[numLevel].x) / jdiff;
            }

//            if (iy < variables->boxMinIndexSC[numLevel].y) {
            if (iy < variables->boxMinIndexSC[numLevel].y && ((iy - variables->boxMinIndexSC[numLevel].y) % jdiff != 0)) {
                iy = (iy - variables->boxMinIndexSC[numLevel].y) / jdiff - 1;
            } else {
                iy = (iy - variables->boxMinIndexSC[numLevel].y) / jdiff;
            }

//            if (iz < variables->boxMinIndexSC[numLevel].z) {
            if (iz < variables->boxMinIndexSC[numLevel].z && ((iz - variables->boxMinIndexSC[numLevel].z) % jdiff != 0)) {
                iz = (iz - variables->boxMinIndexSC[numLevel].z) / jdiff - 1;
            } else {
                iz = (iz - variables->boxMinIndexSC[numLevel].z) / jdiff;
            }

            ixp = (ix + 3) / 3;
            iyp = (iy + 3) / 3;
            izp = (iz + 3) / 3;

            // Loop through all cells within a 3x3x3 box centered around the parent cell
            for (jxp = ixp - 1; jxp <= ixp + 1; jxp++) {
                for (jyp = iyp - 1; jyp <= iyp + 1; jyp++) {
                    for (jzp = izp - 1; jzp <= izp + 1; jzp++) {
                        for (jx = std::max(3 * jxp - 3, jxmin); jx <= std::min(3 * jxp - 1, jxmax); jx++) {
                            for (jy = std::max(3 * jyp - 3, jymin); jy <= std::min(3 * jyp - 1, jymax); jy++) {
                                for (jz = std::max(3 * jzp - 3, jzmin); jz <= std::min(3 * jzp - 1, jzmax); jz++) {
                                    // Cells must be outside of a 3x3x3 box centered around ii-th cell
                                    if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                                        boxIndex3D.x = jx;
                                        boxIndex3D.y = jy;
                                        boxIndex3D.z = jz;
                                        MultipoleMethodFunction::morton3(boxIndex3D, domainBoxIndex, numLevelSC - numLevel);
                                        guessIndex = findMortonIndex(domainBoxIndex, numLevel);
                                        guessIndex = variables->domainBoxIndexMask[guessIndex];
                                        if (guessIndex != -1) {
                                            variables->interactionList[ii][variables->numInteraction[ii]][0] = guessIndex;
                                            variables->numInteraction[ii]++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }// M2L for supercells
    else if (interactionType == FmmInteractionType::M2L_SupercellTopLevel) {
        // Loop through every domain at level numLevel
        //for (ii = 0; ii < numDomainsLevel[numLevel]; ii++){
        for (ii = 0; ii < numBoxIndex; ii++) {
            ib = ii + variables->levelOffset[maxLevel + numLevel];
            variables->numInteraction[ii] = 0;

            // Find box index of domain
            //MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][ii].mortonIndex, boxIndex3D);
            MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][variables->domainBoxIndexFull[ib]].mortonIndex, boxIndex3D);
            ix = boxIndex3D.x;
            iy = boxIndex3D.y;
            iz = boxIndex3D.z;

            // Loop through every domain at level numLevel
            //for (jj = 0; jj < numDomainsLevel[numLevel]; jj++){
            for (jj = 0; jj < numBoxIndex; jj++) {
                jb = jj + variables->levelOffset[maxLevel + numLevel];
                //MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][jj].mortonIndex, boxIndex3D);
                MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][variables->domainBoxIndexFull[jb]].mortonIndex, boxIndex3D);
                jx = boxIndex3D.x;
                jy = boxIndex3D.y;
                jz = boxIndex3D.z;
                if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                    variables->interactionList[ii][variables->numInteraction[ii]][0] = jj;
                    variables->numInteraction[ii]++;
                }
            }
        }
    }// M2L for supercells for outliers
    else if (interactionType == FmmInteractionType::M2L_SupercellTopLevelOutliers) {
        // Loop through every domain at level numLevel

        //        vec3<int> netOffset = {0, 0, 0};
        //        int powThree = pow(3, numLevel);
        //        for (ii = 0; ii < numLevel; ii++){
        //            netOffset.x += powThree * variables->offsetSC[ii].x;
        //            netOffset.y += powThree * variables->offsetSC[ii].y;
        //            netOffset.z += powThree * variables->offsetSC[ii].z;
        //            powThree *= 3;
        //        }
        //int jOffset = jdiff * 3 / 2;

        //for (ii = 0; ii < numDomainsLevel[numLevel]; ii++){
        for (ii = 0; ii < outlierNumBoxIndex; ii++) {
            variables->numInteraction[ii] = 0;

            // Find box index of domain
            //MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][ii].mortonIndex, boxIndex3D);
            MultipoleMethodFunction::signedUnmorton3(variables->outlierDomainLinkage[ii].mortonIndex, boxIndex3D);
            ix = boxIndex3D.x;
            iy = boxIndex3D.y;
            iz = boxIndex3D.z;

//            if (ix < variables->boxMinIndexSC[numLevel].x) {
            if (ix < variables->boxMinIndexSC[numLevel].x && ((ix - variables->boxMinIndexSC[numLevel].x) % jdiff != 0)) {
                ix = (ix - variables->boxMinIndexSC[numLevel].x) / jdiff - 1;
            } else {
                ix = (ix - variables->boxMinIndexSC[numLevel].x) / jdiff;
            }

//            if (iy < variables->boxMinIndexSC[numLevel].y) {
            if (iy < variables->boxMinIndexSC[numLevel].y && ((iy - variables->boxMinIndexSC[numLevel].y) % jdiff != 0)) {
                iy = (iy - variables->boxMinIndexSC[numLevel].y) / jdiff - 1;
            } else {
                iy = (iy - variables->boxMinIndexSC[numLevel].y) / jdiff;
            }

//            if (iz < variables->boxMinIndexSC[numLevel].z) {
            if (iz < variables->boxMinIndexSC[numLevel].z && ((iz - variables->boxMinIndexSC[numLevel].z) % jdiff != 0)) {
                iz = (iz - variables->boxMinIndexSC[numLevel].z) / jdiff - 1;
            } else {
                iz = (iz - variables->boxMinIndexSC[numLevel].z) / jdiff;
            }


            // Loop through every domain at level numLevel
            //for (jj = 0; jj < numDomainsLevel[numLevel]; jj++){
            for (jj = 0; jj < numBoxIndex; jj++) {
                //                jb = jj + variables->levelOffset[maxLevel + numLevel];
                //                //MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][jj].mortonIndex, boxIndex3D);
                //                MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][variables->domainBoxIndexFull[jb]].mortonIndex, boxIndex3D);
                //                jx = variables->boxMinIndexSC[i].x + jdiff * boxIndex3D.x + jdiff / 2;
                //                jy = variables->boxMinIndexSC[i].y + jdiff * boxIndex3D.y + jdiff / 2;
                //                jz = variables->boxMinIndexSC[i].z + jdiff * boxIndex3D.z + jdiff / 2;
                //                if (jx < ix - jOffset || ix + jOffset < jx || jy < iy - jOffset || iy + jOffset < jy || jz < iz - jOffset || iz + jOffset < jz){
                //                    variables->interactionList[ii][variables->numInteraction[ii]] = jj;
                //                    variables->numInteraction[ii]++;
                //                }
                jb = jj + variables->levelOffset[maxLevel + numLevel];
                //MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][jj].mortonIndex, boxIndex3D);
                MultipoleMethodFunction::unmorton3(variables->domainLinkage[numLevel][variables->domainBoxIndexFull[jb]].mortonIndex, boxIndex3D);
                jx = boxIndex3D.x;
                jy = boxIndex3D.y;
                jz = boxIndex3D.z;
                if (jx < ix - 1 || ix + 1 < jx || jy < iy - 1 || iy + 1 < jy || jz < iz - 1 || iz + 1 < jz) {
                    variables->interactionList[ii][variables->numInteraction[ii]][0] = jj;
                    variables->numInteraction[ii]++;
                }
            }
            //            if (numDomains == 3){
            //                printf("M2P variables->numInteraction[%d] = %d\n",ii,variables->numInteraction[ii]);
            //            }
        }
    }
}